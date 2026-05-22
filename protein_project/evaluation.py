import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, precision_recall_fscore_support, roc_auc_score, roc_curve


def compute_binary_metrics(labels: pd.Series, scores: pd.Series) -> dict[str, float]:
    if labels.nunique() < 2:
        return {
            "roc_auc": float("nan"),
            "average_precision": float("nan"),
            "accuracy": float("nan"),
            "precision": float("nan"),
            "recall": float("nan"),
            "f1": float("nan"),
        }
    predictions = (scores >= 0).astype(int)
    precision, recall, f1, _ = precision_recall_fscore_support(labels, predictions, average="binary", zero_division=0)
    return {
        "roc_auc": float(roc_auc_score(labels, scores)),
        "average_precision": float(average_precision_score(labels, scores)),
        "accuracy": float(accuracy_score(labels, predictions)),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
    }


def summarize_zero_shot(dataframe: pd.DataFrame, label_column: str, score_columns: list[str]) -> dict[str, dict[str, float]]:
    return {
        score_column: compute_binary_metrics(dataframe[label_column], dataframe[score_column])
        for score_column in score_columns
        if score_column in dataframe.columns
    }


def save_json(payload: dict[str, Any], save_path: str | Path) -> Path:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    return path


def plot_roc_curves(dataframe: pd.DataFrame, label_column: str, score_columns: list[str], save_path: str | Path) -> Path:
    plt.figure(figsize=(7, 6))
    labels = dataframe[label_column]
    for score_column in score_columns:
        if score_column not in dataframe.columns:
            continue
        if labels.nunique() < 2:
            continue
        fpr, tpr, _ = roc_curve(labels, dataframe[score_column])
        auc = roc_auc_score(labels, dataframe[score_column])
        plt.plot(fpr, tpr, label=f"{score_column} (AUC={auc:.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("TP53 ClinVar zero-shot ROC")
    plt.legend()
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    return path


def run_spatial_linear_probe(
    dataframe: pd.DataFrame,
    feature_columns: list[str],
    label_column: str = "label",
    cluster_column: str = "spatial_cluster",
) -> tuple[dict[str, Any], pd.DataFrame]:
    available = [column for column in feature_columns if column in dataframe.columns]
    if len(available) == 0:
        raise ValueError("No feature columns are available for linear probing")
    clusters = sorted(cluster for cluster in dataframe[cluster_column].dropna().unique() if cluster >= 0)
    if len(clusters) < 2:
        raise ValueError("At least two spatial clusters are required for linear probing")
    test_cluster = clusters[-1]
    train_df = dataframe.loc[dataframe[cluster_column] != test_cluster].copy()
    test_df = dataframe.loc[dataframe[cluster_column] == test_cluster].copy()
    if train_df[label_column].nunique() < 2 or test_df[label_column].nunique() < 2:
        raise ValueError("Train and test splits must both contain positive and negative samples")
    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    model.fit(train_df[available], train_df[label_column])
    train_scores = model.predict_proba(train_df[available])[:, 1]
    test_scores = model.predict_proba(test_df[available])[:, 1]
    metrics = {
        "feature_columns": available,
        "train_cluster_count": int(train_df[cluster_column].nunique()),
        "test_cluster": int(test_cluster),
        "train_metrics": compute_binary_metrics(train_df[label_column], pd.Series(train_scores)),
        "test_metrics": compute_binary_metrics(test_df[label_column], pd.Series(test_scores)),
        "coefficients": {column: float(value) for column, value in zip(available, model.coef_[0])},
        "intercept": float(model.intercept_[0]),
    }
    prediction_df = pd.concat(
        [
            train_df.assign(split="train", linear_probe_score=train_scores),
            test_df.assign(split="test", linear_probe_score=test_scores),
        ],
        ignore_index=True,
    )
    return metrics, prediction_df
