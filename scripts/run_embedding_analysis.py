import argparse
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import average_precision_score, roc_auc_score, silhouette_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from transformers import EsmForMaskedLM, EsmTokenizer

from protein_project.config import ensure_project_dirs, get_processed_path, get_results_path, load_config
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


def extract_mutation_site_embeddings(
    model: EsmForMaskedLM,
    tokenizer: EsmTokenizer,
    device: str,
    combined_sequence: str,
    mutations: list[str],
    batch_size: int,
    desc: str,
) -> np.ndarray:
    combined_tokens = tokenizer.tokenize(combined_sequence)
    embeddings: list[np.ndarray] = []
    for start in tqdm(range(0, len(mutations), batch_size), desc=desc):
        batch = mutations[start : start + batch_size]
        masked_sequences = []
        positions = []
        for mutation in batch:
            wild_type, position, _ = parse_mutation(mutation)
            if combined_tokens[position - 1][0] != wild_type:
                raise ValueError(f"Mutation {mutation} does not match the structure-aware sequence")
            masked_sequences.append(_mask_saprot_position(tokenizer, combined_sequence, position))
            positions.append(position)
        encoded = tokenizer(masked_sequences, return_tensors="pt", padding=True)
        encoded = {key: value.to(device) for key, value in encoded.items()}
        with torch.no_grad():
            output = model(**encoded, output_hidden_states=True)
            hidden = output.hidden_states[-1].detach().cpu()
        for row_index, position in enumerate(positions):
            embeddings.append(hidden[row_index, position].numpy())
    return np.vstack(embeddings)


def compute_embedding_metrics(embeddings: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    benign = embeddings[labels == 0]
    pathogenic = embeddings[labels == 1]
    benign_center = benign.mean(axis=0)
    pathogenic_center = pathogenic.mean(axis=0)
    centroid_distance = float(np.linalg.norm(pathogenic_center - benign_center))
    benign_within = float(np.mean(np.sum((benign - benign_center) ** 2, axis=1)))
    pathogenic_within = float(np.mean(np.sum((pathogenic - pathogenic_center) ** 2, axis=1)))
    within_variance = float((benign_within + pathogenic_within) / 2.0)
    fisher_ratio = float((centroid_distance**2) / (within_variance + 1e-12))
    normalized_centroid_distance = float(centroid_distance / (np.sqrt(within_variance) + 1e-12))
    if len(np.unique(labels)) < 2 or min(np.bincount(labels.astype(int))) < 2:
        linear_auc = float("nan")
        linear_ap = float("nan")
        silhouette = float("nan")
    else:
        n_splits = min(5, int(np.bincount(labels.astype(int)).min()))
        classifier = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight="balanced"),
        )
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=0)
        probabilities = cross_val_predict(classifier, embeddings, labels, cv=cv, method="predict_proba")[:, 1]
        linear_auc = float(roc_auc_score(labels, probabilities))
        linear_ap = float(average_precision_score(labels, probabilities))
        silhouette = float(silhouette_score(embeddings, labels, metric="euclidean"))
    return {
        "centroid_distance": centroid_distance,
        "benign_within_variance": benign_within,
        "pathogenic_within_variance": pathogenic_within,
        "mean_within_variance": within_variance,
        "normalized_centroid_distance": normalized_centroid_distance,
        "fisher_ratio": fisher_ratio,
        "silhouette_score": silhouette,
        "linear_probe_cv_roc_auc": linear_auc,
        "linear_probe_cv_average_precision": linear_ap,
    }


def plot_tsne(embeddings_by_strategy: dict[str, np.ndarray], labels: np.ndarray, save_path: Path) -> None:
    strategy_names = list(embeddings_by_strategy)
    combined = np.vstack([embeddings_by_strategy[name] for name in strategy_names])
    n_components = min(50, combined.shape[1], combined.shape[0] - 1)
    reduced = PCA(n_components=n_components, random_state=0).fit_transform(combined)
    perplexity = max(5, min(30, (combined.shape[0] - 1) // 3))
    coords = TSNE(n_components=2, perplexity=perplexity, init="pca", learning_rate="auto", random_state=0).fit_transform(reduced)

    fig, axes = plt.subplots(1, len(strategy_names), figsize=(6 * len(strategy_names), 5), sharex=True, sharey=True)
    if len(strategy_names) == 1:
        axes = [axes]
    offset = 0
    colors = np.array(["#3b82f6" if label == 0 else "#ef4444" for label in labels])
    for axis, name in zip(axes, strategy_names):
        current = coords[offset : offset + len(labels)]
        offset += len(labels)
        axis.scatter(current[:, 0], current[:, 1], c=colors, s=18, alpha=0.78, edgecolors="none")
        axis.set_title(name)
        axis.set_xlabel("t-SNE 1")
        axis.set_ylabel("t-SNE 2")
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#3b82f6", label="benign", markersize=7),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#ef4444", label="pathogenic", markersize=7),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
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
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--output-prefix", default=None)
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_project_dirs(config)
    output_prefix = args.output_prefix or Path(args.config).stem

    dataset = pd.read_csv(get_processed_path(config, "cleaned_dataset"))
    domain_start = int(config["structure"].get("domain_start", 1))
    domain_end = int(config["structure"].get("domain_end", int(dataset["position"].max())))
    if _use_domain_local_scoring(config):
        dataset = _localize_dataset(dataset, domain_start)

    saprot_payload = load_saprot_sequences(get_processed_path(config, "saprot_sequences"))
    full_combined_seq = saprot_payload["full_combined_seq"]
    masked_combined_seq = saprot_payload["masked_combined_seq"]
    if _use_domain_local_scoring(config):
        full_combined_seq = _slice_combined_sequence(full_combined_seq, domain_start, domain_end)
        masked_combined_seq = _slice_combined_sequence(masked_combined_seq, domain_start, domain_end)

    device = choose_device(args.device)
    tokenizer = EsmTokenizer.from_pretrained(config["models"]["saprot"])
    model = EsmForMaskedLM.from_pretrained(config["models"]["saprot"]).to(device)
    model.eval()

    mutations = dataset["mutation"].tolist()
    labels = dataset["label"].astype(int).to_numpy()
    embeddings_by_strategy = {
        "SaProt full": extract_mutation_site_embeddings(
            model, tokenizer, device, full_combined_seq, mutations, args.batch_size, "SaProt full embeddings"
        ),
        "SaProt pLDDT-mask": extract_mutation_site_embeddings(
            model, tokenizer, device, masked_combined_seq, mutations, args.batch_size, "SaProt masked embeddings"
        ),
    }

    results_dir = Path(config["paths"]["results_dir"])
    metrics = {
        "config": args.config,
        "n_variants": int(len(dataset)),
        "n_pathogenic": int(labels.sum()),
        "n_benign": int((labels == 0).sum()),
        "embedding_definition": "last hidden state at the masked mutation-site token",
        "strategies": {
            name: compute_embedding_metrics(embeddings, labels) for name, embeddings in embeddings_by_strategy.items()
        },
    }
    full_metrics = metrics["strategies"]["SaProt full"]
    masked_metrics = metrics["strategies"]["SaProt pLDDT-mask"]
    metrics["masked_minus_full"] = {
        key: float(masked_metrics[key] - full_metrics[key]) for key in full_metrics if np.isfinite(masked_metrics[key])
    }

    metrics_path = results_dir / f"{output_prefix}_embedding_analysis_metrics.json"
    plot_path = results_dir / f"{output_prefix}_embedding_tsne.png"
    save_json(metrics, metrics_path)
    plot_tsne(embeddings_by_strategy, labels, plot_path)

    embedding_table = dataset.copy()
    for name, embeddings in embeddings_by_strategy.items():
        key = name.lower().replace(" ", "_").replace("-", "_")
        embedding_table[f"{key}_pc1"] = PCA(n_components=1, random_state=0).fit_transform(embeddings)[:, 0]
    scores_path = results_dir / f"{output_prefix}_embedding_analysis_summary.csv"
    embedding_table.to_csv(scores_path, index=False)

    print(f"Saved embedding metrics to {metrics_path}")
    print(f"Saved t-SNE plot to {plot_path}")
    print(f"Saved embedding summary table to {scores_path}")
    print(metrics)


if __name__ == "__main__":
    main()
