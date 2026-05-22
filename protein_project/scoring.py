import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import torch
from tqdm import tqdm
from transformers import EsmForMaskedLM, EsmTokenizer

from protein_project.constants import AA_ALPHABET


FOLDSEEK_STRUC_VOCAB = "pynwrqhgdlvtmfsaeikc#"


def choose_device(explicit_device: str | None = None) -> str:
    if explicit_device:
        return explicit_device
    return "cuda" if torch.cuda.is_available() else "cpu"


def parse_mutation(mutation: str) -> tuple[str, int, str]:
    wild_type = mutation[0]
    mutant = mutation[-1]
    position = int(mutation[1:-1])
    return wild_type, position, mutant


class Esm2ZeroShotScorer:
    def __init__(self, model_name: str, device: str | None = None):
        self.device = choose_device(device)
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmForMaskedLM.from_pretrained(model_name).to(self.device)
        self.model.eval()
        self.aa_token_ids = {aa: self.tokenizer.convert_tokens_to_ids(aa) for aa in AA_ALPHABET}

    def _build_masked_sequence(self, sequence: str, position: int) -> str:
        tokens = list(sequence)
        tokens[position - 1] = self.tokenizer.mask_token
        return " ".join(tokens)

    def score_mutations(
        self,
        sequence: str,
        mutations: Iterable[str],
        batch_size: int = 16,
        show_progress: bool = True,
    ) -> list[float]:
        mutation_list = list(mutations)
        scores: list[float] = []
        iterator = range(0, len(mutation_list), batch_size)
        if show_progress:
            iterator = tqdm(iterator, desc="ESM-2 scoring")
        for start in iterator:
            batch = mutation_list[start : start + batch_size]
            masked_sequences = [self._build_masked_sequence(sequence, parse_mutation(mutation)[1]) for mutation in batch]
            encoded = self.tokenizer.batch_encode_plus(masked_sequences, return_tensors="pt", padding=True)
            encoded = {key: value.to(self.device) for key, value in encoded.items()}
            with torch.no_grad():
                log_probs = torch.log_softmax(self.model(**encoded).logits, dim=-1)
            for row_index, mutation in enumerate(batch):
                wild_type, position, mutant = parse_mutation(mutation)
                if sequence[position - 1] != wild_type:
                    raise ValueError(f"Mutation {mutation} does not match the reference sequence")
                wild_type_id = self.aa_token_ids[wild_type]
                mutant_id = self.aa_token_ids[mutant]
                score = (log_probs[row_index, position, mutant_id] - log_probs[row_index, position, wild_type_id]).item()
                scores.append(score)
        return scores


class SaProtZeroShotScorer:
    def __init__(self, model_name: str, device: str | None = None):
        self.device = choose_device(device)
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmForMaskedLM.from_pretrained(model_name).to(self.device)
        self.model.eval()
        vocab = self.tokenizer.get_vocab()
        self.aa_struct_token_ids = {
            aa: [vocab[token] for token in [f"{aa}{struct_char}" for struct_char in FOLDSEEK_STRUC_VOCAB] if token in vocab]
            for aa in AA_ALPHABET
        }

    def _build_masked_sequence(self, combined_sequence: str, position: int) -> str:
        tokens = self.tokenizer.tokenize(combined_sequence)
        tokens[position - 1] = "#" + tokens[position - 1][-1]
        return " ".join(tokens)

    def score_mutations(
        self,
        combined_sequence: str,
        mutations: Iterable[str],
        batch_size: int = 8,
        show_progress: bool = True,
    ) -> list[float]:
        mutation_list = list(mutations)
        tokens = self.tokenizer.tokenize(combined_sequence)
        scores: list[float] = []
        iterator = range(0, len(mutation_list), batch_size)
        if show_progress:
            iterator = tqdm(iterator, desc="SaProt scoring")
        for start in iterator:
            batch = mutation_list[start : start + batch_size]
            masked_sequences = []
            parsed_batch = []
            for mutation in batch:
                wild_type, position, mutant = parse_mutation(mutation)
                if tokens[position - 1][0] != wild_type:
                    raise ValueError(f"Mutation {mutation} does not match the combined sequence")
                masked_sequences.append(self._build_masked_sequence(sequence, position))
                parsed_batch.append((wild_type, position, mutant))
            encoded = self.tokenizer.batch_encode_plus(masked_sequences, return_tensors="pt", padding=True)
            encoded = {key: value.to(self.device) for key, value in encoded.items()}
            with torch.no_grad():
                log_probs = torch.log_softmax(self.model(**encoded).logits, dim=-1)
            for row_index, (wild_type, position, mutant) in enumerate(parsed_batch):
                wild_ids = self.aa_struct_token_ids[wild_type]
                mutant_ids = self.aa_struct_token_ids[mutant]
                wild_log_prob = torch.logsumexp(log_probs[row_index, position, wild_ids], dim=0)
                mutant_log_prob = torch.logsumexp(log_probs[row_index, position, mutant_ids], dim=0)
                scores.append((mutant_log_prob - wild_log_prob).item())
        return scores


def load_saprot_sequences(sequence_path: str | Path) -> dict:
    return json.loads(Path(sequence_path).read_text())


def score_dataframe(
    dataframe: pd.DataFrame,
    scorer,
    sequence: str,
    mutation_column: str = "mutation",
    score_column: str = "score",
    batch_size: int = 16,
) -> pd.DataFrame:
    result = dataframe.copy()
    result[score_column] = scorer.score_mutations(
        sequence=sequence,
        mutations=result[mutation_column].tolist(),
        batch_size=batch_size,
    )
    return result
