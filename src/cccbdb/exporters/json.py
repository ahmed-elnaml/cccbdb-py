"""Structured JSON exporter for cccbdb-py models and datasets."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

from pydantic import BaseModel

if TYPE_CHECKING:
    from cccbdb.models.molecule import BenchmarkDataset


def export_model_to_json(model: BaseModel, path: Union[str, Path], indent: int = 2) -> Path:
    """Serialize any Pydantic model to a formatted JSON file."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    json_str = model.model_dump_json(indent=indent) if hasattr(model, "model_dump_json") else model.json(indent=indent)
    out_path.write_text(json_str, encoding="utf-8")
    return out_path


def export_dataset_to_json(dataset: "BenchmarkDataset", path: Union[str, Path]) -> Path:
    """Export a :class:`BenchmarkDataset` with full provenance and failure report to JSON."""
    return export_model_to_json(dataset, path, indent=2)
