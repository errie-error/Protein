import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, precision_recall_fscore_support, roc_auc_score, roc_curve


def compute_binary_metrics(labels: pd.Series, scores: pd.Series, threshold: float = 0.0) -> dict[str, float]:
    if labels.nunique() < 2:
        return {
            "roc_auc": float("nan"),
            "average_precision": float("nan"),
            "accuracy": float("nan"),
            "precision": float("nan"),
            "recall": float("nan"),
            "f1": float("nan"),
        }
    predictions = (scores >= threshold).astype(int)
    precision, recall, f1, _ = precision_recall_fscore_support(labels, predictions, average="binary", zero_division=0)
    return {
        "roc_auc": float(roc_auc_score(labels, scores)),
        "average_precision": float(average_precision_score(labels, scores)),
        "accuracy": float(accuracy_score(labels, predictions)),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
    }


def _bootstrap_metric_interval(
    labels: pd.Series,
    scores: pd.Series,
    metric_name: str,
    n_bootstrap: int = 2000,
    seed: int = 0,
) -> list[float]:
    labels_array = np.asarray(labels)
    scores_array = np.asarray(scores)
    if len(labels_array) == 0 or len(np.unique(labels_array)) < 2:
        return [float("nan"), float("nan")]

    rng = np.random.default_rng(seed)
    values: list[float] = []
    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(labels_array), size=len(labels_array))
        sampled_labels = labels_array[indices]
        if len(np.unique(sampled_labels)) < 2:
            continue
        sampled_scores = scores_array[indices]
        if metric_name == "roc_auc":
            values.append(float(roc_auc_score(sampled_labels, sampled_scores)))
        elif metric_name == "average_precision":
            values.append(float(average_precision_score(sampled_labels, sampled_scores)))
        else:
            raise ValueError(f"Unsupported metric for bootstrap interval: {metric_name}")

    if not values:
        return [float("nan"), float("nan")]
    return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]


def _bootstrap_metric_difference_interval(
    labels: pd.Series,
    left_scores: pd.Series,
    right_scores: pd.Series,
    metric_name: str,
    n_bootstrap: int = 2000,
    seed: int = 0,
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
        if metric_name == "roc_auc":
            left_value = roc_auc_score(sampled_labels, sampled_left)
            right_value = roc_auc_score(sampled_labels, sampled_right)
        elif metric_name == "average_precision":
            left_value = average_precision_score(sampled_labels, sampled_left)
            right_value = average_precision_score(sampled_labels, sampled_right)
        else:
            raise ValueError(f"Unsupported metric for bootstrap interval: {metric_name}")
        values.append(float(left_value - right_value))

    if not values:
        return [float("nan"), float("nan")]
    return [float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))]


def add_simple_baselines(dataframe: pd.DataFrame) -> pd.DataFrame:
    result = dataframe.copy()
    if "plddt" in result.columns and "plddt_score" not in result.columns:
        result["plddt_score"] = result["plddt"].astype(float)
    if "position" in result.columns and "region_score" not in result.columns:
        # Region-only baseline: the structured TP53 core (94-312) is scored above the termini.
        result["region_score"] = result["position"].between(94, 312).astype(float)
    return result


def summarize_zero_shot(dataframe: pd.DataFrame, score_columns: list[str]) -> dict[str, dict[str, float]]:
    summary: dict[str, Any] = {}
    labels = dataframe["label"]
    for score_column in score_columns:
        if score_column not in dataframe.columns:
            continue
        threshold = 0.0 if score_column != "plddt_score" else 70.0
        metrics = compute_binary_metrics(labels, dataframe[score_column], threshold=threshold)
        metrics["roc_auc_ci"] = _bootstrap_metric_interval(labels, dataframe[score_column], "roc_auc")
        metrics["average_precision_ci"] = _bootstrap_metric_interval(labels, dataframe[score_column], "average_precision")
        summary[score_column] = metrics

    pairwise_differences: dict[str, Any] = {}
    pairs = [
        ("saprot_full_score", "esm2_score"),
        ("saprot_masked_score", "saprot_full_score"),
        ("saprot_masked_score", "esm2_score"),
    ]
    for left_column, right_column in pairs:
        if left_column not in dataframe.columns or right_column not in dataframe.columns:
            continue
        pairwise_differences[f"{left_column}_minus_{right_column}"] = {
            "roc_auc_diff": float(summary[left_column]["roc_auc"] - summary[right_column]["roc_auc"]),
            "roc_auc_diff_ci": _bootstrap_metric_difference_interval(
                labels,
                dataframe[left_column],
                dataframe[right_column],
                "roc_auc",
            ),
            "average_precision_diff": float(
                summary[left_column]["average_precision"] - summary[right_column]["average_precision"]
            ),
            "average_precision_diff_ci": _bootstrap_metric_difference_interval(
                labels,
                dataframe[left_column],
                dataframe[right_column],
                "average_precision",
            ),
        }

    if pairwise_differences:
        summary["pairwise_differences"] = pairwise_differences
    return summary


def save_json(payload: dict[str, Any], path: str | Path) -> Path:
    save_path = Path(path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    save_path.write_text(json.dumps(payload, indent=2))
    return save_path


def plot_roc_curves(dataframe: pd.DataFrame, score_columns: list[str], path: str | Path) -> Path:
    plt.figure(figsize=(7, 6))
    labels = dataframe["label"]
    for score_column in score_columns:
        if score_column not in dataframe.columns:
            continue
        fpr, tpr, _ = roc_curve(labels, dataframe[score_column])
        auc = roc_auc_score(labels, dataframe[score_column])
        plt.plot(fpr, tpr, label=f"{score_column} (AUC={auc:.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("TP53 ClinVar zero-shot ROC")
    plt.legend()
    save_path = Path(path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()
    return save_path


def run_spatial_linear_probe(dataframe: pd.DataFrame, feature_columns: list[str]) -> tuple[dict[str, Any], pd.DataFrame]:
    available = [column for column in feature_columns if column in dataframe.columns]
    if not available:
        raise ValueError("No feature columns available for linear probing")
    clusters = sorted(cluster for cluster in dataframe["spatial_cluster"].dropna().unique() if cluster >= 0)
    if len(clusters) < 2:
        raise ValueError("Need at least two spatial clusters")
    test_cluster = clusters[-1]
    train_df = dataframe.loc[dataframe["spatial_cluster"] != test_cluster].copy()
    test_df = dataframe.loc[dataframe["spatial_cluster"] == test_cluster].copy()
    if train_df["label"].nunique() < 2 or test_df["label"].nunique() < 2:
        raise ValueError("Train/test clusters must both contain both classes")
    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    model.fit(train_df[available], train_df["label"])
    train_scores = pd.Series(model.predict_proba(train_df[available])[:, 1], index=train_df.index)
    test_scores = pd.Series(model.predict_proba(test_df[available])[:, 1], index=test_df.index)
    metrics = {
        "feature_columns": available,
        "test_cluster": int(test_cluster),
        "train_metrics": compute_binary_metrics(train_df["label"], train_scores, threshold=0.5),
        "test_metrics": compute_binary_metrics(test_df["label"], test_scores, threshold=0.5),
        "coefficients": {column: float(value) for column, value in zip(available, model.coef_[0])},
        "intercept": float(model.intercept_[0]),
    }
    prediction_df = pd.concat(
        [
            train_df.assign(split="train", linear_probe_score=train_scores.values),
            test_df.assign(split="test", linear_probe_score=test_scores.values),
        ],
        ignore_index=True,
    )
    return metrics, prediction_df
