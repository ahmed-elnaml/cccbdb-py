"""CSV and reproducibility manifest exporter for cccbdb-py."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Union

from cccbdb.config import __version__

if TYPE_CHECKING:
    from cccbdb.models.molecule import BenchmarkDataset


def write_dicts_to_csv(rows: List[Dict[str, Any]], path: Union[str, Path]) -> Path:
    """Write a list of dictionaries to a UTF-8 CSV file, preserving missing values as empty cells."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        out_path.write_text("", encoding="utf-8")
        return out_path

    fieldnames: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in fieldnames})
    return out_path


def export_dataset_to_csv(
    dataset: "BenchmarkDataset",
    path: Union[str, Path],
    write_manifest: bool = True,
) -> Path:
    """Export a :class:`BenchmarkDataset` to CSV and optionally write a companion `.manifest.json`."""
    csv_path = write_dicts_to_csv(dataset.to_dicts(), path)
    if write_manifest:
        manifest_path = csv_path.with_suffix(csv_path.suffix + ".manifest.json")
        manifest = {
            "package": "cccbdb-py",
            "version": __version__,
            "created_at": dataset.created_at,
            "method": dataset.method,
            "basis_set": dataset.basis_set,
            "requested_properties": dataset.requested_properties,
            "molecule_count": len(dataset.records),
            "partial_failures": dataset.partial_failures,
            "attribution": "NIST Computational Chemistry Comparison and Benchmark Database (SRD 101), https://cccbdb.nist.gov/",
        }
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return csv_path
