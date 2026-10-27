"""Molecule species models and interactive Molecule handle for NIST CCCBDB."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Union

from pydantic import BaseModel as _PydanticBaseModel, Field
try:
    from pydantic import ConfigDict
except ImportError:
    def ConfigDict(**kwargs): return kwargs

class BaseModel(_PydanticBaseModel):
    if not hasattr(_PydanticBaseModel, "model_dump"):
        def model_dump(self, **kwargs):
            return self.dict(**kwargs)
    if not hasattr(_PydanticBaseModel, "model_dump_json"):
        def model_dump_json(self, **kwargs):
            return self.json(**kwargs)
    if not hasattr(_PydanticBaseModel, "model_copy"):
        def model_copy(self, **kwargs):
            return self.copy(**kwargs)


from cccbdb.chemistry.formula import calculate_molecular_weight, format_cas_hyphenated, to_hill_formula
from cccbdb.chemistry.symmetry import PointGroupInfo, parse_point_group
from cccbdb.models.calculation import CalculationResult
from cccbdb.models.geometry import MolecularGeometry
from cccbdb.models.property import (
    DipoleMoment,
    EnergyCollection,
    PolarizabilityTensor,
    QuadrupoleMoment,
    RotationalConstants,
    ThermochemicalProperty,
    VibrationalSpectrum,
)
from cccbdb.models.provenance import DataProvenance, ReferenceCitation

if TYPE_CHECKING:
    from cccbdb.client import CCCBDB


class SpeciesSummary(BaseModel):
    """Lightweight discovery record for a chemical species in NIST CCCBDB."""

    model_config = ConfigDict(frozen=True)

    id: str = Field(..., description="NIST CCCBDB casno identifier (digits only, e.g. '7732185')")
    cas_rn: str = Field(..., description="Formatted CAS Registry Number (e.g. '7732-18-5')")
    name: str = Field(..., description="Chemical species name (e.g. 'Water', 'Carbon monoxide')")
    formula: str = Field(..., description="Chemical formula (e.g. 'H2O', 'CO')")
    hill_formula: str = Field(..., description="Standardized Hill system formula")
    charge: int = Field(default=0, description="Net molecular charge")
    spin_multiplicity: Optional[int] = Field(default=None, description="Spin multiplicity (2S+1)")
    electronic_state: Optional[str] = Field(default=None, description="Ground electronic state term symbol")
    point_group: Optional[str] = Field(default=None, description="Schoenflies symmetry point group")
    inchi: Optional[str] = Field(default=None, description="IUPAC InChI string")
    inchikey: Optional[str] = Field(default=None, description="IUPAC InChIKey")
    smiles: Optional[str] = Field(default=None, description="SMILES representation")
    molecular_weight: float = Field(..., description="Molecular weight in g/mol")


class ExperimentalDataSummary(BaseModel):
    """Complete experimental reference dataset for a species (parsed from exp2x.asp)."""

    species_id: str
    cas_rn: str
    name: str
    formula: str
    charge: int = 0
    synonyms: List[str] = Field(default_factory=list)
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    smiles: Optional[str] = None
    iupac_name: Optional[str] = None
    electronic_state: Optional[str] = None
    point_group: Optional[str] = None
    thermochemistry: List[ThermochemicalProperty] = Field(default_factory=list)
    geometry: Optional[MolecularGeometry] = None
    vibrations: Optional[VibrationalSpectrum] = None
    rotational_constants: Optional[RotationalConstants] = None
    dipole: Optional[DipoleMoment] = None
    quadrupoles: List[QuadrupoleMoment] = Field(default_factory=list)
    polarizability: Optional[PolarizabilityTensor] = None
    ionization_energy_ev: Optional[float] = None
    proton_affinities_kj_mol: List[Dict[str, Any]] = Field(default_factory=list)
    references: Dict[str, ReferenceCitation] = Field(default_factory=dict)
    provenance: Optional[DataProvenance] = None


class Molecule:
    """Rich interactive handle representing a specific chemical species in NIST CCCBDB.

    Provides bound convenience methods (`get_energies()`, `get_geometry()`,
    `get_vibrations()`, `get_dipole()`, `get_calculation()`, `get_experimental()`)
    matching the proposed Python API.
    """

    def __init__(
        self,
        summary: SpeciesSummary,
        client: Optional["CCCBDB"] = None,
        experimental: Optional[ExperimentalDataSummary] = None,
    ) -> None:
        self._summary = summary
        self._client = client
        self._experimental = experimental

    @property
    def id(self) -> str:
        """Return the NIST CCCBDB species identifier (`casno`)."""
        return self._summary.id

    @property
    def cas_rn(self) -> str:
        """Return the hyphenated CAS Registry Number."""
        return self._summary.cas_rn

    @property
    def name(self) -> str:
        """Return the primary chemical name."""
        return self._summary.name

    @property
    def formula(self) -> str:
        """Return the chemical formula."""
        return self._summary.formula

    @property
    def hill_formula(self) -> str:
        """Return the Hill-ordered chemical formula."""
        return self._summary.hill_formula

    @property
    def charge(self) -> int:
        """Return the net molecular electric charge."""
        return self._summary.charge

    @property
    def point_group(self) -> Optional[str]:
        """Return the molecular Schoenflies point group."""
        if self._experimental and self._experimental.point_group:
            return self._experimental.point_group
        return self._summary.point_group

    @property
    def symmetry(self) -> Optional[PointGroupInfo]:
        """Return detailed symmetry point group metadata."""
        pg = self.point_group
        if not pg:
            return None
        try:
            return parse_point_group(pg)
        except ValueError:
            return None

    @property
    def inchi(self) -> Optional[str]:
        return (self._experimental.inchi if self._experimental else None) or self._summary.inchi

    @property
    def inchikey(self) -> Optional[str]:
        return (self._experimental.inchikey if self._experimental else None) or self._summary.inchikey

    @property
    def smiles(self) -> Optional[str]:
        return (self._experimental.smiles if self._experimental else None) or self._summary.smiles

    @property
    def molecular_weight(self) -> float:
        return self._summary.molecular_weight

    @property
    def summary(self) -> SpeciesSummary:
        return self._summary

    def _require_client(self) -> "CCCBDB":
        if self._client is None:
            raise RuntimeError("Molecule instance is not bound to an active CCCBDB client.")
        return self._client

    def get_experimental(self) -> ExperimentalDataSummary:
        """Retrieve all experimental reference data (`exp2x.asp`) for this molecule."""
        if self._experimental is None:
            self._experimental = self._require_client().get_experimental_data(
                self.id, charge=self.charge
            )
        return self._experimental

    def get_energies(
        self,
        method: Optional[str] = None,
        basis_set: Optional[str] = None,
        energy_type: Optional[str] = None,
        unit: str = "Hartree",
    ) -> EnergyCollection:
        """Retrieve electronic and thermochemical energies for this molecule."""
        return self._require_client().get_energies(
            species_id=self.id,
            charge=self.charge,
            method=method,
            basis_set=basis_set,
            energy_type=energy_type,
            unit=unit,
        )

    def get_geometry(
        self,
        method: str = "Experimental",
        basis_set: Optional[str] = None,
        unit: str = "Angstrom",
    ) -> MolecularGeometry:
        """Retrieve the experimental or calculated optimized molecular geometry."""
        return self._require_client().get_geometry(
            species_id=self.id,
            charge=self.charge,
            method=method,
            basis_set=basis_set,
            unit=unit,
        )

    def get_vibrations(
        self,
        method: str = "Experimental",
        basis_set: Optional[str] = None,
    ) -> VibrationalSpectrum:
        """Retrieve the experimental or calculated vibrational frequency spectrum."""
        return self._require_client().get_vibrations(
            species_id=self.id,
            charge=self.charge,
            method=method,
            basis_set=basis_set,
        )

    def get_dipole(
        self,
        method: str = "Experimental",
        basis_set: Optional[str] = None,
        unit: str = "Debye",
    ) -> Optional[DipoleMoment]:
        """Retrieve the experimental or calculated electric dipole moment."""
        return self._require_client().get_dipole(
            species_id=self.id,
            charge=self.charge,
            method=method,
            basis_set=basis_set,
            unit=unit,
        )

    def get_calculation(
        self,
        method: str,
        basis_set: Optional[str] = None,
    ) -> CalculationResult:
        """Retrieve the full calculation dossier (`energy3x.asp`) for a specific method and basis set."""
        return self._require_client().get_calculation(
            species_id=self.id,
            method=method,
            basis_set=basis_set,
            charge=self.charge,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert molecule metadata to a dictionary."""
        data = self._summary.model_dump()
        if self._experimental:
            data["experimental"] = self._experimental.model_dump()
        return data

    def __repr__(self) -> str:
        return (
            f"Molecule(id='{self.id}', formula='{self.formula}', "
            f"name='{self.name}', point_group='{self.point_group}')"
        )


class DatasetRecord(BaseModel):
    """Row record in a collected multi-molecule benchmark dataset."""

    species_id: str
    cas_rn: str
    name: str
    formula: str
    charge: int = 0
    point_group: Optional[str] = None
    molecular_weight: float
    method: str
    basis_set: Optional[str] = None
    properties: Dict[str, Any] = Field(default_factory=dict)
    source_urls: List[str] = Field(default_factory=list)
    errors: List[str] = Field(default_factory=list)


class BenchmarkDataset(BaseModel):
    """Multi-molecule dataset returned by :meth:`CCCBDB.collect` with manifest & export support."""

    records: List[DatasetRecord] = Field(default_factory=list)
    requested_properties: List[str] = Field(default_factory=list)
    method: str = "B3LYP"
    basis_set: Optional[str] = "6-31G*"
    partial_failures: Dict[str, str] = Field(default_factory=dict)
    created_at: str = Field(default="")

    def to_dicts(self) -> List[Dict[str, Any]]:
        """Flatten dataset records into dictionaries for tabular analysis."""
        rows: List[Dict[str, Any]] = []
        for rec in self.records:
            row: Dict[str, Any] = {
                "species_id": rec.species_id,
                "cas_rn": rec.cas_rn,
                "name": rec.name,
                "formula": rec.formula,
                "charge": rec.charge,
                "point_group": rec.point_group,
                "molecular_weight": rec.molecular_weight,
                "method": rec.method,
                "basis_set": rec.basis_set,
            }
            row.update(rec.properties)
            row["source_urls"] = ";".join(rec.source_urls)
            row["errors"] = ";".join(rec.errors) if rec.errors else None
            rows.append(row)
        return rows

    def to_dataframe(self) -> Any:
        """Convert dataset into a :class:`pandas.DataFrame`."""
        import pandas as pd

        return pd.DataFrame(self.to_dicts())

    def to_csv(self, path: Union[str, Path], write_manifest: bool = True) -> Path:
        """Export dataset to CSV and optionally write a reproducibility JSON manifest alongside it."""
        from cccbdb.exporters.csv import export_dataset_to_csv

        return export_dataset_to_csv(self, path, write_manifest=write_manifest)

    def to_json(self, path: Union[str, Path]) -> Path:
        """Export dataset and provenance manifest to a structured JSON file."""
        from cccbdb.exporters.json import export_dataset_to_json

        return export_dataset_to_json(self, path)
