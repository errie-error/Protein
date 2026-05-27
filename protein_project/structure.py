import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from sklearn.cluster import KMeans

from protein_project.constants import AA3_TO_AA1


LOW_CONFIDENCE_TOKEN = "#"


def _get_parser(structure_path: str | Path):
    path = str(structure_path)
    if path.endswith(".cif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def extract_residue_table(structure_path: str | Path, chain_id: str) -> pd.DataFrame:
    parser = _get_parser(structure_path)
    structure = parser.get_structure("protein", structure_path)
    model = structure[0]
    if chain_id not in model:
        available = [chain.id for chain in model]
        raise ValueError(f"Chain {chain_id} not found. Available chains: {available}")
    chain = model[chain_id]
    rows: list[dict[str, Any]] = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        if "CA" not in residue:
            continue
        residue_name = residue.get_resname().title()
        residue_aa = AA3_TO_AA1.get(residue_name)
        if residue_aa is None:
            continue
        coord = residue["CA"].coord
        plddt = float(np.mean([atom.get_bfactor() for atom in residue.get_atoms()]))
        rows.append(
            {
                "position": int(residue.id[1]),
                "structure_residue": residue_aa,
                "x": float(coord[0]),
                "y": float(coord[1]),
                "z": float(coord[2]),
                "plddt": plddt,
            }
        )
    return pd.DataFrame(rows)


def build_saprot_sequences(
    foldseek_path: str | Path,
    structure_path: str | Path,
    chain_id: str,
    plddt_threshold: float,
    random_seed: int = 0,
) -> dict[str, str]:
    structure_path = Path(structure_path)
    residue_table = extract_residue_table(structure_path, chain_id)
    with tempfile.TemporaryDirectory() as temp_dir:
        descriptor_path = Path(temp_dir) / "foldseek_descriptor.tsv"
        command = [
            str(foldseek_path),
            "structureto3didescriptor",
            "-v",
            "0",
            "--threads",
            "1",
            "--chain-name-mode",
            "1",
            str(structure_path),
            str(descriptor_path),
        ]
        subprocess.run(command, check=True)
        lines = descriptor_path.read_text().strip().splitlines()
    protein_name = structure_path.name
    selected = None
    for line in lines:
        desc, sequence, structure_sequence, *_ = line.split("\t")
        parsed_chain = desc.split(" ")[0].replace(protein_name, "").split("_")[-1]
        if parsed_chain == chain_id:
            selected = (sequence, structure_sequence)
            break
    if selected is None:
        raise ValueError(f"No Foldseek descriptor found for chain {chain_id}")
    sequence, structure_sequence = selected
    if len(sequence) != len(structure_sequence):
        raise ValueError("Foldseek output length mismatch")
    if len(residue_table) != len(structure_sequence):
        raise ValueError("Residue table length does not match Foldseek output")
    full_combined = "".join(a + b.lower() for a, b in zip(sequence, structure_sequence))
    masked_chars = list(structure_sequence)
    low_confidence_positions = residue_table.index[residue_table["plddt"] < plddt_threshold].tolist()
    for index in low_confidence_positions:
        masked_chars[index] = LOW_CONFIDENCE_TOKEN
    masked_structure = "".join(masked_chars)
    masked_combined = "".join(a + b.lower() for a, b in zip(sequence, masked_structure))

    rng = np.random.default_rng(random_seed)
    random_masked_chars = list(structure_sequence)
    if low_confidence_positions:
        random_positions = sorted(
            rng.choice(len(structure_sequence), size=len(low_confidence_positions), replace=False).tolist()
        )
    else:
        random_positions = []
    for index in random_positions:
        random_masked_chars[index] = LOW_CONFIDENCE_TOKEN
    random_masked_structure = "".join(random_masked_chars)
    random_masked_combined = "".join(a + b.lower() for a, b in zip(sequence, random_masked_structure))

    high_masked_chars = list(structure_sequence)
    high_confidence_positions = residue_table.index[residue_table["plddt"] > plddt_threshold].tolist()
    for index in high_confidence_positions:
        high_masked_chars[index] = LOW_CONFIDENCE_TOKEN
    high_masked_structure = "".join(high_masked_chars)
    high_masked_combined = "".join(a + b.lower() for a, b in zip(sequence, high_masked_structure))
    return {
        "sequence": sequence,
        "structure_sequence": structure_sequence.lower(),
        "full_combined_seq": full_combined,
        "masked_structure_sequence": masked_structure.lower(),
        "masked_combined_seq": masked_combined,
        "random_masked_structure_sequence": random_masked_structure.lower(),
        "random_masked_combined_seq": random_masked_combined,
        "high_masked_structure_sequence": high_masked_structure.lower(),
        "high_masked_combined_seq": high_masked_combined,
        "mask_metadata": {
            "plddt_threshold": float(plddt_threshold),
            "random_seed": int(random_seed),
            "low_confidence_count": len(low_confidence_positions),
            "random_mask_count": len(random_positions),
            "high_confidence_count": len(high_confidence_positions),
        },
    }


def assign_spatial_clusters(
    mutation_df: pd.DataFrame,
    residue_df: pd.DataFrame,
    n_clusters: int = 2,
    random_state: int = 0,
) -> pd.DataFrame:
    if mutation_df.empty:
        return mutation_df.copy()
    residue_subset = residue_df.loc[
        residue_df["position"].isin(mutation_df["position"].unique()),
        ["position", "x", "y", "z"],
    ].drop_duplicates(subset=["position"])
    if len(residue_subset) < n_clusters:
        mapping = {position: 0 for position in residue_subset["position"]}
    else:
        model = KMeans(n_clusters=n_clusters, random_state=random_state, n_init="auto")
        labels = model.fit_predict(residue_subset[["x", "y", "z"]])
        mapping = dict(zip(residue_subset["position"], labels))
    result = mutation_df.copy()
    result["spatial_cluster"] = result["position"].map(mapping).fillna(-1).astype(int)
    return result


def save_saprot_sequences(payload: dict[str, Any], save_path: str | Path) -> Path:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path
