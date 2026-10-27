"""Normalized SQLite storage backend for molecules, energies, geometries, and properties."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from cccbdb.models.geometry import MolecularGeometry
from cccbdb.models.molecule import SpeciesSummary
from cccbdb.models.property import DipoleMoment, EnergyRecord


class SQLiteStorage:
    """Normalized relational storage in SQLite for reproducible local benchmarking."""

    def __init__(self, db_path: Union[str, Path]) -> None:
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.db_path))
        self._conn.row_factory = sqlite3.Row
        self._init_schema()

    def _init_schema(self) -> None:
        with self._conn:
            self._conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS species (
                    id TEXT NOT NULL,
                    charge INTEGER NOT NULL DEFAULT 0,
                    cas_rn TEXT NOT NULL,
                    name TEXT NOT NULL,
                    formula TEXT NOT NULL,
                    hill_formula TEXT NOT NULL,
                    point_group TEXT,
                    electronic_state TEXT,
                    inchi TEXT,
                    inchikey TEXT,
                    smiles TEXT,
                    molecular_weight REAL NOT NULL,
                    PRIMARY KEY (id, charge)
                );

                CREATE TABLE IF NOT EXISTS energies (
                    species_id TEXT NOT NULL,
                    formula TEXT NOT NULL,
                    method TEXT NOT NULL,
                    basis_set TEXT NOT NULL,
                    energy_type TEXT NOT NULL,
                    conformer_index INTEGER NOT NULL DEFAULT 1,
                    value REAL NOT NULL,
                    unit TEXT NOT NULL,
                    method_category TEXT,
                    core_correlation TEXT,
                    point_group TEXT,
                    electronic_state TEXT,
                    source_url TEXT,
                    retrieved_at TEXT,
                    PRIMARY KEY (species_id, method, basis_set, energy_type, conformer_index)
                );

                CREATE TABLE IF NOT EXISTS geometries (
                    species_id TEXT NOT NULL,
                    formula TEXT NOT NULL,
                    method TEXT NOT NULL,
                    basis_set TEXT NOT NULL,
                    conformer_index INTEGER NOT NULL DEFAULT 1,
                    point_group TEXT,
                    unit TEXT NOT NULL,
                    xyz_json TEXT NOT NULL,
                    source_url TEXT,
                    PRIMARY KEY (species_id, method, basis_set, conformer_index)
                );

                CREATE TABLE IF NOT EXISTS dipoles (
                    species_id TEXT NOT NULL,
                    formula TEXT NOT NULL,
                    method TEXT NOT NULL,
                    basis_set TEXT NOT NULL,
                    x REAL,
                    y REAL,
                    z REAL,
                    total REAL NOT NULL,
                    unit TEXT NOT NULL,
                    point_group TEXT,
                    source_url TEXT,
                    PRIMARY KEY (species_id, method, basis_set)
                );
                """
            )

    def save_species(self, species: SpeciesSummary) -> None:
        """Insert or update a :class:`SpeciesSummary` in the SQLite database."""
        with self._conn:
            self._conn.execute(
                """
                INSERT OR REPLACE INTO species
                (id, charge, cas_rn, name, formula, hill_formula, point_group,
                 electronic_state, inchi, inchikey, smiles, molecular_weight)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    species.id,
                    species.charge,
                    species.cas_rn,
                    species.name,
                    species.formula,
                    species.hill_formula,
                    species.point_group,
                    species.electronic_state,
                    species.inchi,
                    species.inchikey,
                    species.smiles,
                    species.molecular_weight,
                ),
            )

    def get_species(self, species_id: str, charge: int = 0) -> Optional[SpeciesSummary]:
        """Retrieve a :class:`SpeciesSummary` by ``species_id`` and ``charge``."""
        row = self._conn.execute(
            "SELECT * FROM species WHERE id = ? AND charge = ?", (species_id, charge)
        ).fetchone()
        if not row:
            return None
        return SpeciesSummary(**dict(row))

    def save_energies(self, records: List[EnergyRecord]) -> int:
        """Persist a batch of :class:`EnergyRecord` items to SQLite."""
        with self._conn:
            for r in records:
                self._conn.execute(
                    """
                    INSERT OR REPLACE INTO energies
                    (species_id, formula, method, basis_set, energy_type, conformer_index,
                     value, unit, method_category, core_correlation, point_group,
                     electronic_state, source_url, retrieved_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        r.species_id,
                        r.formula,
                        r.method,
                        r.basis_set or "Predefined",
                        r.energy_type,
                        r.conformer_index,
                        r.value,
                        r.unit,
                        r.method_category,
                        r.core_correlation,
                        r.point_group,
                        r.electronic_state,
                        r.provenance.source_url if r.provenance else r.detail_url,
                        r.provenance.retrieved_at if r.provenance else None,
                    ),
                )
        return len(records)

    def query_energies(
        self,
        species_id: Optional[str] = None,
        method: Optional[str] = None,
        basis_set: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """Query stored energy records with optional filters."""
        clauses: List[str] = []
        params: List[Any] = []
        if species_id:
            clauses.append("species_id = ?")
            params.append(species_id)
        if method:
            clauses.append("LOWER(method) = LOWER(?)")
            params.append(method)
        if basis_set:
            clauses.append("LOWER(basis_set) = LOWER(?)")
            params.append(basis_set)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        rows = self._conn.execute(f"SELECT * FROM energies {where}", params).fetchall()
        return [dict(r) for r in rows]

    def save_geometry(self, geom: MolecularGeometry) -> None:
        """Persist a :class:`MolecularGeometry` record to SQLite."""
        xyz_json = json.dumps([s.model_dump() for s in geom.sites])
        with self._conn:
            self._conn.execute(
                """
                INSERT OR REPLACE INTO geometries
                (species_id, formula, method, basis_set, conformer_index, point_group, unit, xyz_json, source_url)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    geom.species_id,
                    geom.formula,
                    geom.method,
                    geom.basis_set or "Predefined",
                    geom.conformer_index,
                    geom.point_group,
                    geom.unit,
                    xyz_json,
                    geom.provenance.source_url if geom.provenance else None,
                ),
            )

    def save_dipole(self, dipole: DipoleMoment) -> None:
        """Persist a :class:`DipoleMoment` record to SQLite."""
        with self._conn:
            self._conn.execute(
                """
                INSERT OR REPLACE INTO dipoles
                (species_id, formula, method, basis_set, x, y, z, total, unit, point_group, source_url)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    dipole.species_id,
                    dipole.formula,
                    dipole.method,
                    dipole.basis_set or "Predefined",
                    dipole.x,
                    dipole.y,
                    dipole.z,
                    dipole.total,
                    dipole.unit,
                    dipole.point_group,
                    dipole.provenance.source_url if dipole.provenance else None,
                ),
            )

    def close(self) -> None:
        self._conn.close()
