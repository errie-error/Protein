import re
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from protein_project.constants import AA3_TO_AA1


AA3_CHANGE_PATTERN = re.compile(r"p\.\[?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)\]?")
AA1_CHANGE_PATTERN = re.compile(r"p\.\[?([A-Z])(\d+)([A-Z]|\*|=)\]?")
PROTEIN_CHANGE_FIELDS = [
    "Protein_change",
    "ProteinChange",
    "protein_change",
    "AminoAcidChange",
    "Amino_acid_change",
    "Name",
]
VARIATION_ID_FIELDS = ["VariationID", "AlleleID", "#AlleleID"]
TYPE_FIELDS = ["Type", "VariantType"]


def download_file(url: str, destination: str | Path, force: bool = False) -> Path:
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not force:
        return path
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return path


def read_fasta_sequence(fasta_path: str | Path) -> str:
    lines = []
    with Path(fasta_path).open() as handle:
        for line in handle:
            if not line.startswith(">"):
                lines.append(line.strip())
    return "".join(lines)


def load_variant_summary(summary_path: str | Path) -> pd.DataFrame:
    return pd.read_csv(summary_path, sep="\t", compression="infer", dtype=str, low_memory=False)


def _first_present(row: dict[str, Any], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key)
        if value is None:
            continue
        value = str(value).strip()
        if value and value.lower() != "nan":
            return value
    return ""


def parse_protein_change(raw_value: str) -> tuple[str, int, str] | None:
    if not raw_value:
        return None
    aa3_match = AA3_CHANGE_PATTERN.search(raw_value)
    if aa3_match:
        wt3, position, mut3 = aa3_match.groups()
        wt = AA3_TO_AA1.get(wt3)
        mut = AA3_TO_AA1.get(mut3)
        if wt and mut and mut not in {"*"}:
            return wt, int(position), mut
        return None
    aa1_match = AA1_CHANGE_PATTERN.search(raw_value)
    if aa1_match:
        wt, position, mut = aa1_match.groups()
        if mut not in {"*", "="}:
            return wt, int(position), mut
    return None


def clinical_significance_to_label(raw_value: str) -> int | None:
    value = str(raw_value or "").lower()
    blocked_terms = [
        "conflict",
        "uncertain",
        "not provided",
        "association",
        "drug response",
        "risk factor",
        "protective",
        "other",
        "affects",
    ]
    if any(term in value for term in blocked_terms):
        return None
    if "pathogenic" in value and "benign" in value:
        return None
    if "pathogenic" in value:
        return 1
    if "benign" in value:
        return 0
    return None


def clean_clinvar_table(
    clinvar_df: pd.DataFrame,
    gene: str,
    reference_sequence: str,
    domain_start: int,
    domain_end: int,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for row in clinvar_df.to_dict(orient="records"):
        if str(row.get("GeneSymbol", "")).strip() != gene:
            continue
        label = clinical_significance_to_label(row.get("ClinicalSignificance", ""))
        if label is None:
            continue
        raw_change = _first_present(row, PROTEIN_CHANGE_FIELDS)
        parsed = parse_protein_change(raw_change)
        if parsed is None:
            continue
        wild_type, position, mutant = parsed
        if position < domain_start or position > domain_end:
            continue
        if position > len(reference_sequence):
            continue
        if reference_sequence[position - 1] != wild_type:
            continue
        mutation = f"{wild_type}{position}{mutant}"
        records.append(
            {
                "gene": gene,
                "mutation": mutation,
                "wild_type": wild_type,
                "position": position,
                "mutant": mutant,
                "label": label,
                "clinical_significance": row.get("ClinicalSignificance", ""),
                "review_status": row.get("ReviewStatus", ""),
                "variant_type": _first_present(row, TYPE_FIELDS),
                "protein_change": raw_change,
                "variation_id": _first_present(row, VARIATION_ID_FIELDS),
            }
        )
    cleaned = pd.DataFrame(records)
    if cleaned.empty:
        return cleaned
    inconsistent = cleaned.groupby("mutation")["label"].nunique()
    inconsistent = inconsistent[inconsistent > 1].index
    cleaned = cleaned.loc[~cleaned["mutation"].isin(inconsistent)].copy()
    cleaned = cleaned.sort_values(["position", "mutation", "variation_id"])
    cleaned = cleaned.drop_duplicates(subset=["mutation"], keep="first")
    cleaned = cleaned.reset_index(drop=True)
    return cleaned
