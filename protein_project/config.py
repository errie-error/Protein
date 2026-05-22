from pathlib import Path
from typing import Any

import yaml


RESULT_IMAGE_NAME = "tp53_zero_shot_roc.png"
RESULT_PREDICTIONS_NAME = "tp53_linear_probe_predictions.csv"
RESIDUE_TABLE_NAME = "tp53_residue_table.csv"


def load_config(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path).resolve()
    config = yaml.safe_load(path.read_text())
    repo_root = path.parent.parent
    config["repo_root"] = str(repo_root)
    for key, value in config.get("paths", {}).items():
        config["paths"][key] = str((repo_root / value).resolve())
    return config


def ensure_project_dirs(config: dict[str, Any]) -> dict[str, Path]:
    resolved = {}
    for key, value in config.get("paths", {}).items():
        directory = Path(value)
        directory.mkdir(parents=True, exist_ok=True)
        resolved[key] = directory
    return resolved


def get_raw_path(config: dict[str, Any], filename_key: str) -> Path:
    return Path(config["paths"]["raw_dir"]) / config["filenames"][filename_key]


def get_processed_path(config: dict[str, Any], filename_key: str) -> Path:
    return Path(config["paths"]["processed_dir"]) / config["filenames"][filename_key]


def get_results_path(config: dict[str, Any], filename_key: str) -> Path:
    return Path(config["paths"]["results_dir"]) / config["filenames"][filename_key]
