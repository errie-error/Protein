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
    return mutation[0], int(mutation[1:-1]), mutation[-1]


class Esm2ZeroShotScorer:
    def __init__(self, model_name: str, device: str | None = None):
        self.device = choose_device(device)
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmForMaskedLM.from_pretrained(model_name).to(self.device)
        self.model.eval()
        self.aa_token_ids = {aa: self.tokenizer.convert_tokens_to_ids(aa) for aa in AA_ALPHABET}

    def _mask_sequence(self, sequence: str, position: int) -> str:
        tokens = list(sequence)
        tokens[position - 1] = self.tokenizer.mask_token
        return " ".join(tokens)

    def score_mutations(self, sequence: str, mutations: Iterable[str], batch_size: int = 16) -> list[float]:
        mutation_list = list(mutations)
        scores: list[float] = []
        for start in tqdm(range(0, len(mutation_list), batch_size), desc="ESM-2 scoring"):
            batch = mutation_list[start : start + batch_size]
            masked = [self._mask_sequence(sequence, parse_mutation(mutation)[1]) for mutation in batch]
            encoded = self.tokenizer(masked, return_tensors="pt", padding=True)
            encoded = {key: value.to(self.device) for key, value in encoded.items()}
            with torch.no_grad():
                log_probs = torch.log_softmax(self.model(**encoded).logits, dim=-1)
            for row_index, mutation in enumerate(batch):
                wild_type, position, mutant = parse_mutation(mutation)
                if sequence[position - 1] != wild_type:
                    raise ValueError(f"Mutation {mutation} does not match the sequence")
                wt_id = self.aa_token_ids[wild_type]
                mut_id = self.aa_token_ids[mutant]
                # Higher score means the wild-type residue is preferred over the mutant,
                # which aligns larger values with more deleterious mutations.
                scores.append((log_probs[row_index, position, wt_id] - log_probs[row_index, position, mut_id]).item())
        return scores


class SaProtZeroShotScorer:
    def __init__(self, model_name: str, device: str | None = None):
        self.device = choose_device(device)
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmForMaskedLM.from_pretrained(model_name).to(self.device)
        self.model.eval()
        vocab = self.tokenizer.get_vocab()
        self.aa_struct_ids = {
            aa: [vocab[token] for token in [f"{aa}{struct_char}" for struct_char in FOLDSEEK_STRUC_VOCAB] if token in vocab]
            for aa in AA_ALPHABET
        }

    def _mask_sequence(self, combined_sequence: str, position: int) -> str:
        tokens = self.tokenizer.tokenize(combined_sequence)
        tokens[position - 1] = "#" + tokens[position - 1][-1]
        return " ".join(tokens)

    def score_mutations(self, combined_sequence: str, mutations: Iterable[str], batch_size: int = 8) -> list[float]:
        mutation_list = list(mutations)
        combined_tokens = self.tokenizer.tokenize(combined_sequence)
        scores: list[float] = []
        for start in tqdm(range(0, len(mutation_list), batch_size), desc="SaProt scoring"):
            batch = mutation_list[start : start + batch_size]
            masked = []
            parsed = []
            for mutation in batch:
                wild_type, position, mutant = parse_mutation(mutation)
                if combined_tokens[position - 1][0] != wild_type:
                    raise ValueError(f"Mutation {mutation} does not match the structure-aware sequence")
                masked.append(self._mask_sequence(combined_sequence, position))
                parsed.append((wild_type, position, mutant))
            encoded = self.tokenizer(masked, return_tensors="pt", padding=True)
            encoded = {key: value.to(self.device) for key, value in encoded.items()}
            with torch.no_grad():
                log_probs = torch.log_softmax(self.model(**encoded).logits, dim=-1)
            for row_index, (wild_type, position, mutant) in enumerate(parsed):
                wt_log_prob = torch.logsumexp(log_probs[row_index, position, self.aa_struct_ids[wild_type]], dim=0)
                mut_log_prob = torch.logsumexp(log_probs[row_index, position, self.aa_struct_ids[mutant]], dim=0)
                # Higher score means the wild-type residue is preferred over the mutant,
                # which aligns larger values with more deleterious mutations.
                scores.append((wt_log_prob - mut_log_prob).item())
        return scores


def score_dataframe(dataframe: pd.DataFrame, scorer, sequence: str, score_column: str, batch_size: int) -> pd.DataFrame:
    result = dataframe.copy()
    result[score_column] = scorer.score_mutations(sequence, result["mutation"].tolist(), batch_size=batch_size)
    return result


def load_saprot_sequences(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())
