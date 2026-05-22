import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
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


def summarize_zero_shot(dataframe: pd.DataFrame, score_columns: list[str]) -> dict[str, dict[str, float]]:
    return {
        score_column: compute_binary_metrics(dataframe["label"], dataframe[score_column], threshold=0.0)
        for score_column in score_columns
        if score_column in dataframe.columns
    }


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
