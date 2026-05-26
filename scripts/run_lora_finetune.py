import argparse
from pathlib import Path

import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm
from transformers import EsmForSequenceClassification, EsmTokenizer

from protein_project.benchmarks import compute_binary_metrics, save_json
from protein_project.config import ensure_project_dirs, get_processed_path, get_raw_path, load_config
from protein_project.data import read_fasta_sequence
from protein_project.zero_shot import load_saprot_sequences, parse_mutation


try:
    from peft import LoraConfig, TaskType, get_peft_model
except ModuleNotFoundError as exc:
    raise SystemExit("Missing dependency: install LoRA support with `pip install peft>=0.13`.") from exc


class MutationClassificationDataset(Dataset):
    def __init__(self, dataframe: pd.DataFrame, tokenizer: EsmTokenizer, sequence_column: str):
        self.dataframe = dataframe.reset_index(drop=True)
        self.tokenizer = tokenizer
        self.sequence_column = sequence_column

    def __len__(self) -> int:
        return len(self.dataframe)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        row = self.dataframe.iloc[index]
        encoded = self.tokenizer(
            row[self.sequence_column],
            return_tensors="pt",
            truncation=True,
            max_length=1024,
        )
        return {
            "input_ids": encoded["input_ids"].squeeze(0),
            "attention_mask": encoded["attention_mask"].squeeze(0),
            "labels": torch.tensor(int(row["label"]), dtype=torch.long),
        }


def collate_batch(batch: list[dict[str, torch.Tensor]], tokenizer: EsmTokenizer) -> dict[str, torch.Tensor]:
    input_ids = [item["input_ids"] for item in batch]
    attention_mask = [item["attention_mask"] for item in batch]
    labels = torch.stack([item["labels"] for item in batch])
    padded = tokenizer.pad(
        {"input_ids": input_ids, "attention_mask": attention_mask},
        return_tensors="pt",
    )
    padded["labels"] = labels
    return padded


def build_mutant_sequence(sequence: str, mutation: str) -> str:
    wild_type, position, mutant = parse_mutation(mutation)
    tokens = list(sequence)
    if tokens[position - 1] != wild_type:
        raise ValueError(f"Mutation {mutation} does not match the reference sequence")
    tokens[position - 1] = mutant
    return " ".join(tokens)


def build_mutant_saprot_sequence(combined_sequence: str, mutation: str, tokenizer: EsmTokenizer) -> str:
    wild_type, position, mutant = parse_mutation(mutation)
    tokens = tokenizer.tokenize(combined_sequence)
    if tokens[position - 1][0] != wild_type:
        raise ValueError(f"Mutation {mutation} does not match the structure-aware sequence")
    tokens[position - 1] = mutant + tokens[position - 1][-1]
    return " ".join(tokens)


def add_model_inputs(
    dataframe: pd.DataFrame,
    mode: str,
    tokenizer: EsmTokenizer,
    reference_sequence: str,
    saprot_payload: dict | None,
) -> tuple[pd.DataFrame, str]:
    result = dataframe.copy()
    sequence_column = f"{mode}_input"
    if mode == "esm2":
        result[sequence_column] = [build_mutant_sequence(reference_sequence, mutation) for mutation in result["mutation"]]
        return result, sequence_column

    if saprot_payload is None:
        raise ValueError("SaProt payload is required for saprot_full and saprot_masked modes")
    payload_key = "masked_combined_seq" if mode == "saprot_masked" else "full_combined_seq"
    combined_sequence = saprot_payload[payload_key]
    result[sequence_column] = [
        build_mutant_saprot_sequence(combined_sequence, mutation, tokenizer) for mutation in result["mutation"]
    ]
    return result, sequence_column


def split_spatial(dataframe: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    clusters = sorted(cluster for cluster in dataframe["spatial_cluster"].dropna().unique() if cluster >= 0)
    if len(clusters) < 2:
        raise ValueError("At least two spatial clusters are required for LoRA fine-tuning")
    test_cluster = int(clusters[-1])
    train_df = dataframe.loc[dataframe["spatial_cluster"] != test_cluster].copy()
    test_df = dataframe.loc[dataframe["spatial_cluster"] == test_cluster].copy()
    if train_df["label"].nunique() < 2 or test_df["label"].nunique() < 2:
        raise ValueError("Train and test splits must both contain positive and negative samples")
    return train_df, test_df, test_cluster


def evaluate_model(model, dataloader: DataLoader, device: str) -> tuple[dict[str, float], pd.Series]:
    model.eval()
    probabilities: list[float] = []
    labels: list[int] = []
    with torch.no_grad():
        for batch in dataloader:
            batch = {key: value.to(device) for key, value in batch.items()}
            outputs = model(**batch)
            probs = torch.softmax(outputs.logits, dim=-1)[:, 1]
            probabilities.extend(probs.detach().cpu().tolist())
            labels.extend(batch["labels"].detach().cpu().tolist())
    score_series = pd.Series(probabilities)
    label_series = pd.Series(labels)
    return compute_binary_metrics(label_series, score_series, threshold=0.5), score_series


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/tp53.yaml")
    parser.add_argument("--mode", choices=["saprot_masked", "saprot_full", "esm2"], default="saprot_masked")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=0.01)
    parser.add_argument("--lora-r", type=int, default=8)
    parser.add_argument("--lora-alpha", type=int, default=16)
    parser.add_argument("--lora-dropout", type=float, default=0.05)
    parser.add_argument("--max-train-steps", type=int, default=None)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    config = load_config(args.config)
    ensure_project_dirs(config)

    model_path = config["models"]["esm2"] if args.mode == "esm2" else config["models"]["saprot"]
    tokenizer = EsmTokenizer.from_pretrained(model_path)
    dataset = pd.read_csv(get_processed_path(config, "cleaned_dataset"))
    reference_sequence = read_fasta_sequence(get_raw_path(config, "reference_fasta"))
    saprot_payload = None
    if args.mode != "esm2":
        saprot_payload = load_saprot_sequences(get_processed_path(config, "saprot_sequences"))
    dataset, sequence_column = add_model_inputs(dataset, args.mode, tokenizer, reference_sequence, saprot_payload)
    train_df, test_df, test_cluster = split_spatial(dataset)

    train_dataset = MutationClassificationDataset(train_df, tokenizer, sequence_column)
    test_dataset = MutationClassificationDataset(test_df, tokenizer, sequence_column)
    train_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=True,
        collate_fn=lambda batch: collate_batch(batch, tokenizer),
    )
    test_loader = DataLoader(
        test_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        collate_fn=lambda batch: collate_batch(batch, tokenizer),
    )

    model = EsmForSequenceClassification.from_pretrained(model_path, num_labels=2).to(args.device)
    lora_config = LoraConfig(
        task_type=TaskType.SEQ_CLS,
        r=args.lora_r,
        lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout,
        target_modules=["query", "value"],
    )
    model = get_peft_model(model, lora_config)
    model.print_trainable_parameters()

    class_counts = train_df["label"].value_counts().sort_index()
    class_weights = torch.tensor(
        [len(train_df) / (2 * class_counts.get(label, 1)) for label in [0, 1]],
        dtype=torch.float32,
        device=args.device,
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay)
    loss_fn = torch.nn.CrossEntropyLoss(weight=class_weights)

    train_losses: list[float] = []
    global_step = 0
    for epoch in range(args.epochs):
        model.train()
        epoch_losses: list[float] = []
        for batch in tqdm(train_loader, desc=f"LoRA fine-tuning epoch {epoch + 1}/{args.epochs}"):
            batch = {key: value.to(args.device) for key, value in batch.items()}
            labels = batch.pop("labels")
            optimizer.zero_grad(set_to_none=True)
            outputs = model(**batch)
            loss = loss_fn(outputs.logits, labels)
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.detach().cpu()))
            global_step += 1
            if args.max_train_steps is not None and global_step >= args.max_train_steps:
                break
        train_losses.append(float(sum(epoch_losses) / max(len(epoch_losses), 1)))
        if args.max_train_steps is not None and global_step >= args.max_train_steps:
            break

    test_metrics, test_scores = evaluate_model(model, test_loader, args.device)
    train_eval_loader = DataLoader(
        train_dataset,
        batch_size=args.batch_size,
        shuffle=False,
        collate_fn=lambda batch: collate_batch(batch, tokenizer),
    )
    train_metrics, train_scores = evaluate_model(model, train_eval_loader, args.device)

    metrics = {
        "mode": args.mode,
        "model_path": model_path,
        "test_cluster": test_cluster,
        "epochs_requested": args.epochs,
        "train_steps": global_step,
        "train_size": int(len(train_df)),
        "test_size": int(len(test_df)),
        "train_losses": train_losses,
        "train_metrics": train_metrics,
        "test_metrics": test_metrics,
        "lora": {
            "r": args.lora_r,
            "alpha": args.lora_alpha,
            "dropout": args.lora_dropout,
            "target_modules": ["query", "value"],
        },
    }

    results_dir = Path(config["paths"]["results_dir"])
    output_prefix = f"{config['gene'].lower()}_lora_{args.mode}"
    metrics_path = results_dir / f"{output_prefix}_metrics.json"
    predictions_path = results_dir / f"{output_prefix}_predictions.csv"
    save_json(metrics, metrics_path)

    prediction_df = pd.concat(
        [
            train_df.assign(split="train", lora_score=train_scores.values),
            test_df.assign(split="test", lora_score=test_scores.values),
        ],
        ignore_index=True,
    )
    prediction_df.drop(columns=[sequence_column]).to_csv(predictions_path, index=False)

    print(f"Saved LoRA metrics to {metrics_path}")
    print(f"Saved LoRA predictions to {predictions_path}")
    print(metrics)


if __name__ == "__main__":
    main()
