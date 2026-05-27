import argparse
from pathlib import Path

import pandas as pd

from protein_project.benchmarks import add_simple_baselines, plot_roc_curves, save_json, summarize_zero_shot
from protein_project.config import ensure_project_dirs, get_processed_path, get_raw_path, get_results_path, load_config
from protein_project.data import read_fasta_sequence
from protein_project.zero_shot import Esm2ZeroShotScorer, SaProtZeroShotScorer, load_saprot_sequences, score_dataframe


def _use_domain_local_scoring(config: dict) -> bool:
    return bool(config.get("scoring", {}).get("domain_local", False))


def _slice_combined_sequence(combined_sequence: str, domain_start: int, domain_end: int) -> str:
    return combined_sequence[(domain_start - 1) * 2 : domain_end * 2]


def _localize_dataset(dataset: pd.DataFrame, domain_start: int) -> pd.DataFrame:
    localized = dataset.copy()
    localized["full_position"] = localized["position"]
    localized["full_mutation"] = localized["mutation"]
    localized["local_position"] = localized["position"] - domain_start + 1
    localized["local_mutation"] = (
        localized["wild_type"] + localized["local_position"].astype(str) + localized["mutant"]
    )
    localized["mutation"] = localized["local_mutation"]
    return localized


def _restore_full_coordinates(dataframe: pd.DataFrame) -> pd.DataFrame:
    restored = dataframe.copy()
    if {"full_position", "full_mutation"}.issubset(restored.columns):
        restored["position"] = restored["full_position"]
        restored["mutation"] = restored["full_mutation"]
    return restored


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/tp53.yaml")
    parser.add_argument("--device", default=None)
    parser.add_argument("--esm-batch-size", type=int, default=16)
    parser.add_argument("--saprot-batch-size", type=int, default=8)
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_project_dirs(config)

    dataset_path = get_processed_path(config, "cleaned_dataset")
    dataset = pd.read_csv(dataset_path)

    fasta_path = get_raw_path(config, "reference_fasta")
    sequence = read_fasta_sequence(fasta_path)
    domain_start = int(config["structure"].get("domain_start", 1))
    domain_end = int(config["structure"].get("domain_end", len(sequence)))
    if _use_domain_local_scoring(config):
        sequence = sequence[domain_start - 1 : domain_end]
        dataset = _localize_dataset(dataset, domain_start)

    esm2 = Esm2ZeroShotScorer(config["models"]["esm2"], device=args.device)
    scored = score_dataframe(dataset, esm2, sequence, score_column="esm2_score", batch_size=args.esm_batch_size)

    score_columns = ["esm2_score"]
    saprot_path = get_processed_path(config, "saprot_sequences")
    if Path(saprot_path).exists():
        saprot_payload = load_saprot_sequences(saprot_path)
        full_combined_seq = saprot_payload["full_combined_seq"]
        masked_combined_seq = saprot_payload["masked_combined_seq"]
        if _use_domain_local_scoring(config):
            full_combined_seq = _slice_combined_sequence(full_combined_seq, domain_start, domain_end)
            masked_combined_seq = _slice_combined_sequence(masked_combined_seq, domain_start, domain_end)
        saprot = SaProtZeroShotScorer(config["models"]["saprot"], device=args.device)
        scored = score_dataframe(scored, saprot, full_combined_seq, score_column="saprot_full_score", batch_size=args.saprot_batch_size)
        scored = score_dataframe(scored, saprot, masked_combined_seq, score_column="saprot_masked_score", batch_size=args.saprot_batch_size)
        score_columns.extend(["saprot_full_score", "saprot_masked_score"])

    scored = _restore_full_coordinates(scored)
    scored = add_simple_baselines(scored)
    for baseline_column in ["plddt_score", "region_score"]:
        if baseline_column in scored.columns:
            score_columns.append(baseline_column)

    results_path = get_results_path(config, "zero_shot_scores")
    metrics_path = get_results_path(config, "zero_shot_metrics")
    roc_path = results_path.with_name(f"{results_path.stem}_roc.png")

    scored.to_csv(results_path, index=False)
    metrics = summarize_zero_shot(scored, score_columns)
    save_json(metrics, metrics_path)
    plot_roc_curves(scored, score_columns, roc_path)

    print(f"Saved zero-shot scores to {results_path}")
    print(f"Saved metrics to {metrics_path}")
    print(metrics)


if __name__ == "__main__":
    main()
