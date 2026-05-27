import argparse

import pandas as pd

from protein_project.benchmarks import run_spatial_linear_probe, save_json
from protein_project.config import ensure_project_dirs, get_results_path, load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/tp53.yaml")
    parser.add_argument("--scores", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_project_dirs(config)

    score_path = args.scores or get_results_path(config, "zero_shot_scores")
    dataframe = pd.read_csv(score_path)

    preferred = ["esm2_score", "saprot_full_score", "saprot_masked_score", "plddt"]
    metrics, prediction_df = run_spatial_linear_probe(dataframe, preferred)

    metrics_path = get_results_path(config, "linear_probe_metrics")
    predictions_path = metrics_path.with_name(f"{metrics_path.stem.replace('_metrics', '')}_predictions.csv")
    save_json(metrics, metrics_path)
    prediction_df.to_csv(predictions_path, index=False)

    print(f"Saved linear probe metrics to {metrics_path}")
    print(f"Saved predictions to {predictions_path}")
    print(metrics)


if __name__ == "__main__":
    main()
