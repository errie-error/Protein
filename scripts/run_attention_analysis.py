import argparse
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from transformers import EsmForMaskedLM, EsmTokenizer

from protein_project.config import ensure_project_dirs, get_processed_path, get_raw_path, get_results_path, load_config
from protein_project.structure import extract_residue_table
from protein_project.zero_shot import choose_device, load_saprot_sequences, parse_mutation


def _use_domain_local_scoring(config: dict[str, Any]) -> bool:
    return bool(config.get("scoring", {}).get("domain_local", False))


def _slice_combined_sequence(combined_sequence: str, domain_start: int, domain_end: int) -> str:
    return combined_sequence[(domain_start - 1) * 2 : domain_end * 2]


def _localize_dataset(dataset: pd.DataFrame, domain_start: int) -> pd.DataFrame:
    localized = dataset.copy()
    localized["full_position"] = localized["position"]
    localized["full_mutation"] = localized["mutation"]
    localized["local_position"] = localized["position"] - domain_start + 1
    localized["local_mutation"] = localized["wild_type"] + localized["local_position"].astype(str) + localized["mutant"]
    localized["mutation"] = localized["local_mutation"]
    return localized


def _mask_saprot_position(tokenizer: EsmTokenizer, combined_sequence: str, position: int) -> str:
    tokens = tokenizer.tokenize(combined_sequence)
    tokens[position - 1] = "#" + tokens[position - 1][-1]
    return " ".join(tokens)


def select_variants(dataframe: pd.DataFrame, max_variants: int) -> pd.DataFrame:
    df = dataframe.copy()
    if {"saprot_masked_score", "saprot_full_score", "label"}.issubset(df.columns):
        score_delta = df["saprot_masked_score"] - df["saprot_full_score"]
        df["correct_direction_delta"] = np.where(df["label"] == 1, score_delta, -score_delta)
        df = df.sort_values("correct_direction_delta", ascending=False)
    if len(df) <= max_variants:
        return df.reset_index(drop=True)
    per_class = max_variants // 2
    selected = []
    for label in [0, 1]:
        part = df.loc[df["label"] == label].head(per_class)
        selected.append(part)
    result = pd.concat(selected, ignore_index=True)
    if len(result) < max_variants:
        remaining = df.drop(index=result.index, errors="ignore").head(max_variants - len(result))
        result = pd.concat([result, remaining], ignore_index=True)
    return result.reset_index(drop=True)


def get_attention_profile(
    model: EsmForMaskedLM,
    tokenizer: EsmTokenizer,
    device: str,
    combined_sequence: str,
    mutation: str,
) -> tuple[np.ndarray, np.ndarray]:
    _, position, _ = parse_mutation(mutation)
    masked_sequence = _mask_saprot_position(tokenizer, combined_sequence, position)
    encoded = tokenizer(masked_sequence, return_tensors="pt")
    encoded = {key: value.to(device) for key, value in encoded.items()}
    with torch.no_grad():
        output = model(**encoded, output_attentions=True)
    # layers x heads x query_len x key_len, excluding batch dimension.
    attentions = torch.stack([layer_attention[0].detach().cpu() for layer_attention in output.attentions])
    # Token index equals residue position because token 0 is BOS.
    query_attention = attentions[:, :, position, 1 : len(tokenizer.tokenize(combined_sequence)) + 1]
    mean_profile = query_attention.mean(dim=(0, 1)).numpy()
    layer_profile = query_attention.mean(dim=1).numpy()
    return mean_profile, layer_profile


def summarize_profile(profile: np.ndarray, plddt: np.ndarray, position: int, threshold: float, local_window: int) -> dict[str, float]:
    low_mask = plddt < threshold
    high_mask = ~low_mask
    positions = np.arange(1, len(profile) + 1)
    local_mask = np.abs(positions - position) <= local_window
    profile_sum = float(profile.sum())
    if profile_sum > 0:
        normalized = profile / profile_sum
    else:
        normalized = profile
    entropy = float(-(normalized * np.log(normalized + 1e-12)).sum())
    return {
        "low_plddt_attention_mass": float(profile[low_mask].sum()),
        "high_plddt_attention_mass": float(profile[high_mask].sum()),
        "local_attention_mass": float(profile[local_mask].sum()),
        "nonlocal_attention_mass": float(profile[~local_mask].sum()),
        "attention_entropy": entropy,
        "top_attention_position": int(np.argmax(profile) + 1),
        "top_attention_plddt": float(plddt[int(np.argmax(profile))]),
    }


def plot_attention_case(
    full_profile: np.ndarray,
    masked_profile: np.ndarray,
    plddt: np.ndarray,
    mutation: str,
    label: int,
    save_path: Path,
    window: int,
) -> None:
    _, position, _ = parse_mutation(mutation)
    start = max(1, position - window)
    end = min(len(full_profile), position + window)
    indices = np.arange(start, end + 1)
    full_slice = full_profile[start - 1 : end]
    masked_slice = masked_profile[start - 1 : end]
    plddt_slice = plddt[start - 1 : end]

    fig, axes = plt.subplots(3, 1, figsize=(12, 5.8), sharex=True, gridspec_kw={"height_ratios": [1, 1, 0.8]})
    vmax = max(float(full_slice.max()), float(masked_slice.max()), 1e-8)
    axes[0].imshow(full_slice[np.newaxis, :], aspect="auto", cmap="magma", vmin=0, vmax=vmax)
    axes[0].set_ylabel("Full")
    axes[1].imshow(masked_slice[np.newaxis, :], aspect="auto", cmap="magma", vmin=0, vmax=vmax)
    axes[1].set_ylabel("Masked")
    axes[2].plot(indices, plddt_slice, color="#2563eb", linewidth=1.8)
    axes[2].axhline(90, color="#ef4444", linestyle="--", linewidth=1)
    axes[2].set_ylabel("pLDDT")
    axes[2].set_xlabel("Residue position")
    for axis in axes[:2]:
        axis.axvline(position - start, color="#22c55e", linestyle="--", linewidth=1.5)
        axis.set_yticks([])
    axes[2].axvline(position, color="#22c55e", linestyle="--", linewidth=1.5)
    tick_step = max(1, len(indices) // 12)
    axes[2].set_xticks(indices[::tick_step])
    label_text = "pathogenic" if label == 1 else "benign"
    fig.suptitle(f"Attention from mutation site: {mutation} ({label_text})")
    fig.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=220)
    plt.close(fig)


def save_json(payload: dict[str, Any], path: Path) -> None:
    import json

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/tp53_plddt90.yaml")
    parser.add_argument("--device", default=None)
    parser.add_argument("--max-variants", type=int, default=80)
    parser.add_argument("--local-window", type=int, default=10)
    parser.add_argument("--plot-window", type=int, default=40)
    parser.add_argument("--n-cases", type=int, default=2)
    parser.add_argument("--output-prefix", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_project_dirs(config)
    output_prefix = args.output_prefix or Path(args.config).stem
    threshold = float(config["structure"].get("plddt_mask_threshold", 90.0))

    scores = pd.read_csv(get_results_path(config, "zero_shot_scores"))
    selected = select_variants(scores, args.max_variants)
    domain_start = int(config["structure"].get("domain_start", 1))
    domain_end = int(config["structure"].get("domain_end", int(scores["position"].max())))
    if _use_domain_local_scoring(config):
        selected = _localize_dataset(selected, domain_start)

    saprot_payload = load_saprot_sequences(get_processed_path(config, "saprot_sequences"))
    full_combined_seq = saprot_payload["full_combined_seq"]
    masked_combined_seq = saprot_payload["masked_combined_seq"]
    if _use_domain_local_scoring(config):
        full_combined_seq = _slice_combined_sequence(full_combined_seq, domain_start, domain_end)
        masked_combined_seq = _slice_combined_sequence(masked_combined_seq, domain_start, domain_end)

    residue_table = extract_residue_table(get_raw_path(config, "alphafold_pdb"), config["structure"]["chain_id"])
    residue_table = residue_table.sort_values("position")
    plddt_by_position = residue_table["plddt"].to_numpy(dtype=float)
    if _use_domain_local_scoring(config):
        plddt_by_position = plddt_by_position[domain_start - 1 : domain_end]

    device = choose_device(args.device)
    tokenizer = EsmTokenizer.from_pretrained(config["models"]["saprot"])
    model = EsmForMaskedLM.from_pretrained(config["models"]["saprot"]).to(device)
    model.eval()

    rows: list[dict[str, Any]] = []
    case_profiles: list[tuple[pd.Series, np.ndarray, np.ndarray]] = []
    for _, row in tqdm(selected.iterrows(), total=len(selected), desc="Attention analysis"):
        full_profile, _ = get_attention_profile(model, tokenizer, device, full_combined_seq, row["mutation"])
        masked_profile, _ = get_attention_profile(model, tokenizer, device, masked_combined_seq, row["mutation"])
        _, local_position, _ = parse_mutation(row["mutation"])
        full_summary = summarize_profile(full_profile, plddt_by_position, local_position, threshold, args.local_window)
        masked_summary = summarize_profile(masked_profile, plddt_by_position, local_position, threshold, args.local_window)
        record = row.to_dict()
        for key, value in full_summary.items():
            record[f"full_{key}"] = value
        for key, value in masked_summary.items():
            record[f"masked_{key}"] = value
            record[f"delta_{key}"] = value - full_summary[key]
        rows.append(record)
        if len(case_profiles) < args.n_cases:
            case_profiles.append((row.copy(), full_profile, masked_profile))

    result = pd.DataFrame(rows)
    metric_columns = [column for column in result.columns if column.startswith("delta_")]
    aggregate = {
        "config": args.config,
        "n_selected_variants": int(len(result)),
        "selection": "balanced by label, prioritized by correct-direction SaProt masked-minus-full score change",
        "plddt_threshold": threshold,
        "local_window": args.local_window,
        "mean_full_low_plddt_attention_mass": float(result["full_low_plddt_attention_mass"].mean()),
        "mean_masked_low_plddt_attention_mass": float(result["masked_low_plddt_attention_mass"].mean()),
        "mean_delta_low_plddt_attention_mass": float(result["delta_low_plddt_attention_mass"].mean()),
        "mean_full_high_plddt_attention_mass": float(result["full_high_plddt_attention_mass"].mean()),
        "mean_masked_high_plddt_attention_mass": float(result["masked_high_plddt_attention_mass"].mean()),
        "mean_delta_high_plddt_attention_mass": float(result["delta_high_plddt_attention_mass"].mean()),
        "mean_full_local_attention_mass": float(result["full_local_attention_mass"].mean()),
        "mean_masked_local_attention_mass": float(result["masked_local_attention_mass"].mean()),
        "mean_delta_local_attention_mass": float(result["delta_local_attention_mass"].mean()),
        "mean_full_attention_entropy": float(result["full_attention_entropy"].mean()),
        "mean_masked_attention_entropy": float(result["masked_attention_entropy"].mean()),
        "mean_delta_attention_entropy": float(result["delta_attention_entropy"].mean()),
        "deltas": {column: float(result[column].mean()) for column in metric_columns},
    }

    results_dir = Path(config["paths"]["results_dir"])
    csv_path = results_dir / f"{output_prefix}_attention_analysis.csv"
    json_path = results_dir / f"{output_prefix}_attention_analysis_metrics.json"
    result.to_csv(csv_path, index=False)
    save_json(aggregate, json_path)

    for case_index, (case_row, full_profile, masked_profile) in enumerate(case_profiles, start=1):
        mutation_name = str(case_row.get("full_mutation", case_row["mutation"])).replace("*", "X")
        plot_path = results_dir / f"{output_prefix}_attention_case{case_index}_{mutation_name}.png"
        plot_attention_case(
            full_profile,
            masked_profile,
            plddt_by_position,
            case_row["mutation"],
            int(case_row["label"]),
            plot_path,
            args.plot_window,
        )

    print(f"Saved attention table to {csv_path}")
    print(f"Saved attention metrics to {json_path}")
    print(aggregate)


if __name__ == "__main__":
    main()
