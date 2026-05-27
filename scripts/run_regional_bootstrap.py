import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from protein_project.benchmarks import compute_binary_metrics, save_json


DEFAULT_REGIONS = [
    ("n_terminal_1_93", 1, 93),
    ("structured_core_94_312", 94, 312),
    ("c_terminal_313_393", 313, 393),
]


def parse_region(raw_region: str) -> tuple[str, int, int]:
    name, bounds = raw_region.split(":", maxsplit=1)
    start, end = bounds.split("-", maxsplit=1)
    return name, int(start), int(end)


def bootstrap_difference_ci(
    labels: pd.Series,
    left_scores: pd.Series,
    right_scores: pd.Series,
    metric: str,
    n_bootstrap: int,
    seed: int,
) -> list[float]:
    labels_array = np.asarray(labels)
    left_array = np.asarray(left_scores)
    right_array = np.asarray(right_scores)
    if len(labels_array) == 0 or len(np.unique(labels_array)) < 2:
        return [float("nan"), float("nan")]

    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(labels_array), size=len(labels_array))
        sampled_labels = labels_array[indices]
        if len(np.unique(sampled_labels)) < 2:
            continue
        sampled_left = left_array[indices]
        sampled_right = right_array[indices]
        if metric == "roc_auc":
            left_value = roc_auc_score(sampled_labels, sampled_left)
            right_value = roc_auc_score(sampled_labels, sampled_right)
        elif metric == "average_precision":
            left_value = average_precision_score(sampled_labels, sampled_left)
            right_value = average_precision_score(sampled_labels, sampled_right)
        else:
            raise ValueError(f"Unsupported metric: {metric}")
        values.append(float(left_value - right_value))

    if not values:
        return [float("nan"), float("nan")]
    return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]


def safe_metrics(dataframe: pd.DataFrame, score_column: str) -> dict[str, float]:
    if score_column not in dataframe.columns or dataframe["label"].nunique() < 2:
        return {
            "roc_auc": float("nan"),
            "average_precision": float("nan"),
            "accuracy": float("nan"),
            "precision": float("nan"),
            "recall": float("nan"),
            "f1": float("nan"),
        }
    return compute_binary_metrics(dataframe["label"], dataframe[score_column], threshold=0.0)


def summarize_region(region_df: pd.DataFrame, n_bootstrap: int, seed: int) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "n": int(len(region_df)),
        "pathogenic": int((region_df["label"] == 1).sum()),
        "benign": int((region_df["label"] == 0).sum()),
        "mean_plddt": float(region_df["plddt"].mean()) if "plddt" in region_df.columns and len(region_df) else float("nan"),
        "plddt_lt_70": int((region_df["plddt"] < 70).sum()) if "plddt" in region_df.columns else 0,
        "plddt_lt_90": int((region_df["plddt"] < 90).sum()) if "plddt" in region_df.columns else 0,
        "metrics": {},
        "masked_minus_full": {},
    }
    for score_column in ["esm2_score", "saprot_full_score", "saprot_masked_score", "plddt_score", "region_score"]:
        if score_column in region_df.columns:
            summary["metrics"][score_column] = safe_metrics(region_df, score_column)

    if {"saprot_masked_score", "saprot_full_score"}.issubset(region_df.columns) and region_df["label"].nunique() >= 2:
        masked = summary["metrics"]["saprot_masked_score"]
        full = summary["metrics"]["saprot_full_score"]
        summary["masked_minus_full"] = {
            "roc_auc_diff": float(masked["roc_auc"] - full["roc_auc"]),
            "roc_auc_diff_ci": bootstrap_difference_ci(
                region_df["label"],
                region_df["saprot_masked_score"],
                region_df["saprot_full_score"],
                "roc_auc",
                n_bootstrap,
                seed,
            ),
            "average_precision_diff": float(masked["average_precision"] - full["average_precision"]),
            "average_precision_diff_ci": bootstrap_difference_ci(
                region_df["label"],
                region_df["saprot_masked_score"],
                region_df["saprot_full_score"],
                "average_precision",
                n_bootstrap,
                seed,
            ),
        }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--region",
        action="append",
        default=None,
        help="Region as name:start-end. May be passed multiple times.",
    )
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    dataframe = pd.read_csv(args.scores)
    results: dict[str, Any] = {
        "scores": str(Path(args.scores).resolve()),
        "n_bootstrap": args.n_bootstrap,
        "regions": {},
    }
    rows: list[dict[str, Any]] = []
    regions = [parse_region(region) for region in args.region] if args.region else DEFAULT_REGIONS
    for region_name, start, end in regions:
        region_df = dataframe.loc[dataframe["position"].between(start, end)].copy()
        region_summary = summarize_region(region_df, args.n_bootstrap, args.seed)
        results["regions"][region_name] = region_summary
        metrics = region_summary["metrics"]
        diff = region_summary["masked_minus_full"]
        rows.append(
            {
                "region": region_name,
                "start": start,
                "end": end,
                "n": region_summary["n"],
                "pathogenic": region_summary["pathogenic"],
                "benign": region_summary["benign"],
                "mean_plddt": region_summary["mean_plddt"],
                "plddt_lt_70": region_summary["plddt_lt_70"],
                "plddt_lt_90": region_summary["plddt_lt_90"],
                "esm2_auc": metrics.get("esm2_score", {}).get("roc_auc", float("nan")),
                "full_auc": metrics.get("saprot_full_score", {}).get("roc_auc", float("nan")),
                "masked_auc": metrics.get("saprot_masked_score", {}).get("roc_auc", float("nan")),
                "masked_minus_full_auc": diff.get("roc_auc_diff", float("nan")),
                "masked_minus_full_auc_ci_low": diff.get("roc_auc_diff_ci", [float("nan"), float("nan")])[0],
                "masked_minus_full_auc_ci_high": diff.get("roc_auc_diff_ci", [float("nan"), float("nan")])[1],
                "full_ap": metrics.get("saprot_full_score", {}).get("average_precision", float("nan")),
                "masked_ap": metrics.get("saprot_masked_score", {}).get("average_precision", float("nan")),
                "masked_minus_full_ap": diff.get("average_precision_diff", float("nan")),
                "masked_minus_full_ap_ci_low": diff.get("average_precision_diff_ci", [float("nan"), float("nan")])[0],
                "masked_minus_full_ap_ci_high": diff.get("average_precision_diff_ci", [float("nan"), float("nan")])[1],
            }
        )

    output_path = Path(args.output)
    save_json(results, output_path)
    table_path = output_path.with_suffix(".csv")
    pd.DataFrame(rows).to_csv(table_path, index=False)
    print(f"Saved regional bootstrap JSON to {output_path}")
    print(f"Saved regional bootstrap table to {table_path}")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
