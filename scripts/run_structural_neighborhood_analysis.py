import argparse
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from protein_project.config import ensure_project_dirs, get_raw_path, get_results_path, load_config
from protein_project.structure import extract_residue_table


def compute_nearest_low_plddt_distance(residue_table: pd.DataFrame, threshold: float) -> pd.DataFrame:
    residues = residue_table.sort_values("position").copy()
    coordinates = residues[["x", "y", "z"]].to_numpy(dtype=float)
    low_coordinates = residues.loc[residues["plddt"] < threshold, ["x", "y", "z"]].to_numpy(dtype=float)
    if len(low_coordinates) == 0:
        distances = np.full(len(residues), np.nan)
    else:
        distances = cdist(coordinates, low_coordinates).min(axis=1)
    residues["nearest_low_plddt_distance"] = distances
    return residues[["position", "nearest_low_plddt_distance"]]


def score_gain(dataframe: pd.DataFrame) -> pd.Series:
    raw_delta = dataframe["saprot_masked_score"] - dataframe["saprot_full_score"]
    # SaProt scores are oriented so larger values imply stronger deleterious evidence.
    return pd.Series(np.where(dataframe["label"] == 1, raw_delta, -raw_delta), index=dataframe.index)


def compute_bin_metrics(dataframe: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    result = dataframe.copy()
    result = result.dropna(subset=["nearest_low_plddt_distance"])
    result["distance_bin"] = pd.qcut(
        result["nearest_low_plddt_distance"],
        q=min(n_bins, result["nearest_low_plddt_distance"].nunique()),
        duplicates="drop",
    )
    rows: list[dict[str, Any]] = []
    for interval, group in result.groupby("distance_bin", observed=True):
        row: dict[str, Any] = {
            "distance_bin": str(interval),
            "n": int(len(group)),
            "n_pathogenic": int(group["label"].sum()),
            "n_benign": int((group["label"] == 0).sum()),
            "mean_distance": float(group["nearest_low_plddt_distance"].mean()),
            "median_distance": float(group["nearest_low_plddt_distance"].median()),
            "mean_correct_direction_gain": float(group["correct_direction_score_gain"].mean()),
        }
        if group["label"].nunique() >= 2:
            full_auc = float(roc_auc_score(group["label"], group["saprot_full_score"]))
            masked_auc = float(roc_auc_score(group["label"], group["saprot_masked_score"]))
            full_ap = float(average_precision_score(group["label"], group["saprot_full_score"]))
            masked_ap = float(average_precision_score(group["label"], group["saprot_masked_score"]))
            row.update(
                {
                    "full_auc": full_auc,
                    "masked_auc": masked_auc,
                    "masked_minus_full_auc": masked_auc - full_auc,
                    "full_ap": full_ap,
                    "masked_ap": masked_ap,
                    "masked_minus_full_ap": masked_ap - full_ap,
                }
            )
        else:
            row.update(
                {
                    "full_auc": float("nan"),
                    "masked_auc": float("nan"),
                    "masked_minus_full_auc": float("nan"),
                    "full_ap": float("nan"),
                    "masked_ap": float("nan"),
                    "masked_minus_full_ap": float("nan"),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def compute_summary(dataframe: pd.DataFrame, bins: pd.DataFrame, config_path: str, threshold: float) -> dict[str, Any]:
    valid = dataframe.dropna(subset=["nearest_low_plddt_distance", "correct_direction_score_gain"])
    rho_gain, p_gain = spearmanr(valid["nearest_low_plddt_distance"], valid["correct_direction_score_gain"])
    bin_valid = bins.dropna(subset=["mean_distance", "masked_minus_full_auc"])
    if len(bin_valid) >= 2:
        rho_auc, p_auc = spearmanr(bin_valid["mean_distance"], bin_valid["masked_minus_full_auc"])
    else:
        rho_auc, p_auc = float("nan"), float("nan")
    return {
        "config": config_path,
        "plddt_threshold": threshold,
        "n_variants": int(len(dataframe)),
        "n_pathogenic": int(dataframe["label"].sum()),
        "n_benign": int((dataframe["label"] == 0).sum()),
        "mean_nearest_low_plddt_distance": float(valid["nearest_low_plddt_distance"].mean()),
        "median_nearest_low_plddt_distance": float(valid["nearest_low_plddt_distance"].median()),
        "mean_correct_direction_score_gain": float(valid["correct_direction_score_gain"].mean()),
        "spearman_distance_vs_correct_direction_gain": float(rho_gain),
        "spearman_distance_vs_correct_direction_gain_pvalue": float(p_gain),
        "spearman_bin_distance_vs_auc_gain": float(rho_auc),
        "spearman_bin_distance_vs_auc_gain_pvalue": float(p_auc),
    }


def plot_scatter(dataframe: pd.DataFrame, save_path: Path, title: str) -> None:
    colors = np.where(dataframe["label"] == 1, "#ef4444", "#3b82f6")
    fig, axis = plt.subplots(figsize=(7, 5))
    axis.scatter(
        dataframe["nearest_low_plddt_distance"],
        dataframe["correct_direction_score_gain"],
        c=colors,
        alpha=0.7,
        s=22,
        edgecolors="none",
    )
    if dataframe["nearest_low_plddt_distance"].nunique() > 1:
        x = dataframe["nearest_low_plddt_distance"].to_numpy(dtype=float)
        y = dataframe["correct_direction_score_gain"].to_numpy(dtype=float)
        coefficients = np.polyfit(x, y, deg=1)
        xs = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), 100)
        axis.plot(xs, coefficients[0] * xs + coefficients[1], color="#111827", linewidth=1.5, label="linear trend")
    axis.axhline(0, color="#6b7280", linewidth=1, linestyle="--")
    axis.set_xlabel("Distance to nearest low-pLDDT residue (Angstrom)")
    axis.set_ylabel("Correct-direction score gain (masked - full)")
    axis.set_title(title)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#3b82f6", label="benign", markersize=7),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#ef4444", label="pathogenic", markersize=7),
    ]
    axis.legend(handles=handles, loc="best")
    fig.tight_layout()
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=220)
    plt.close(fig)


def plot_binned_auc(bins: pd.DataFrame, save_path: Path, title: str) -> None:
    fig, axis = plt.subplots(figsize=(7, 5))
    valid = bins.dropna(subset=["masked_minus_full_auc"])
    axis.plot(valid["mean_distance"], valid["masked_minus_full_auc"], marker="o", linewidth=2, color="#7c3aed")
    axis.axhline(0, color="#6b7280", linewidth=1, linestyle="--")
    for _, row in valid.iterrows():
        axis.annotate(f"n={int(row['n'])}", (row["mean_distance"], row["masked_minus_full_auc"]), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
    axis.set_xlabel("Mean distance to nearest low-pLDDT residue (Angstrom)")
    axis.set_ylabel("Binned AUC gain (masked - full)")
    axis.set_title(title)
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
    parser.add_argument("--output-prefix", default=None)
    parser.add_argument("--n-bins", type=int, default=4)
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_project_dirs(config)
    output_prefix = args.output_prefix or Path(args.config).stem
    threshold = float(config["structure"].get("plddt_mask_threshold", 90.0))

    scores = pd.read_csv(get_results_path(config, "zero_shot_scores"))
    residue_table = extract_residue_table(get_raw_path(config, "alphafold_pdb"), config["structure"]["chain_id"])
    distances = compute_nearest_low_plddt_distance(residue_table, threshold)
    merged = scores.merge(distances, on="position", how="left")
    merged["score_delta_masked_minus_full"] = merged["saprot_masked_score"] - merged["saprot_full_score"]
    merged["correct_direction_score_gain"] = score_gain(merged)

    bins = compute_bin_metrics(merged, args.n_bins)
    summary = compute_summary(merged, bins, args.config, threshold)

    results_dir = Path(config["paths"]["results_dir"])
    detail_path = results_dir / f"{output_prefix}_structural_neighborhood.csv"
    bin_path = results_dir / f"{output_prefix}_structural_neighborhood_bins.csv"
    metrics_path = results_dir / f"{output_prefix}_structural_neighborhood_metrics.json"
    scatter_path = results_dir / f"{output_prefix}_structural_neighborhood_scatter.png"
    bin_plot_path = results_dir / f"{output_prefix}_structural_neighborhood_binned_auc.png"

    merged.to_csv(detail_path, index=False)
    bins.to_csv(bin_path, index=False)
    save_json(summary, metrics_path)
    plot_scatter(merged, scatter_path, f"{config['gene']}: distance to low-pLDDT region vs masking gain")
    plot_binned_auc(bins, bin_plot_path, f"{config['gene']}: binned AUC gain by structural distance")

    print(f"Saved structural neighborhood table to {detail_path}")
    print(f"Saved binned metrics to {bin_path}")
    print(f"Saved summary metrics to {metrics_path}")
    print(f"Saved scatter plot to {scatter_path}")
    print(f"Saved binned AUC plot to {bin_plot_path}")
    print(summary)


if __name__ == "__main__":
    main()
