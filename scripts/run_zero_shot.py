import argparse
from pathlib import Path

import pandas as pd

from protein_project.benchmarks import add_simple_baselines, plot_roc_curves, save_json, summarize_zero_shot
from protein_project.config import ensure_project_dirs, get_processed_path, get_raw_path, get_results_path, load_config
from protein_project.data import read_fasta_sequence
from protein_project.zero_shot import Esm2ZeroShotScorer, SaProtZeroShotScorer, load_saprot_sequences, score_dataframe


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

    esm2 = Esm2ZeroShotScorer(config["models"]["esm2"], device=args.device)
    scored = score_dataframe(dataset, esm2, sequence, score_column="esm2_score", batch_size=args.esm_batch_size)

    score_columns = ["esm2_score"]
    saprot_path = get_processed_path(config, "saprot_sequences")
    if Path(saprot_path).exists():
        saprot_payload = load_saprot_sequences(saprot_path)
        saprot = SaProtZeroShotScorer(config["models"]["saprot"], device=args.device)
        scored = score_dataframe(scored, saprot, saprot_payload["full_combined_seq"], score_column="saprot_full_score", batch_size=args.saprot_batch_size)
        scored = score_dataframe(scored, saprot, saprot_payload["masked_combined_seq"], score_column="saprot_masked_score", batch_size=args.saprot_batch_size)
        score_columns.extend(["saprot_full_score", "saprot_masked_score"])

    scored = add_simple_baselines(scored)
    for baseline_column in ["plddt_score", "region_score"]:
        if baseline_column in scored.columns:
            score_columns.append(baseline_column)

    results_path = get_results_path(config, "zero_shot_scores")
    metrics_path = get_results_path(config, "zero_shot_metrics")
    roc_path = results_path.with_name("tp53_zero_shot_roc.png")

    scored.to_csv(results_path, index=False)
    metrics = summarize_zero_shot(scored, score_columns)
    save_json(metrics, metrics_path)
    plot_roc_curves(scored, score_columns, roc_path)

    print(f"Saved zero-shot scores to {results_path}")
    print(f"Saved metrics to {metrics_path}")
    print(metrics)


if __name__ == "__main__":
    main()
