"""Property models for energies, vibrational frequencies, thermochemistry, and electrostatics."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

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


from cccbdb.chemistry.units import (
    DipoleUnit,
    EnergyUnit,
    PolarizabilityUnit,
    RotationalConstantUnit,
    convert_dipole,
    convert_energy,
    convert_polarizability,
    convert_rotational_constant,
)
from cccbdb.models.provenance import DataProvenance


class EnergyRecord(BaseModel):
    """Single electronic or thermochemical energy record for a species, method, and basis set."""

    model_config = ConfigDict(frozen=True)

    species_id: str = Field(..., description="CCCBDB species identifier (CAS RN)")
    formula: str = Field(..., description="Chemical formula")
    method: str = Field(..., description="Computational method (e.g. 'B3LYP', 'CCSD(T)', 'G4')")
    basis_set: Optional[str] = Field(default=None, description="Basis set name (e.g. '6-31G*', 'cc-pVTZ')")
    energy_type: str = Field(
        default="electronic_0K",
        description="Type of energy ('electronic_0K', 'energy_298K', 'hf_energy', 'nuclear_repulsion', 'zpe')",
    )
    value: float = Field(..., description="Numerical energy value")
    unit: str = Field(default=EnergyUnit.HARTREE.value, description="Energy unit (default: 'Hartree')")
    method_category: Optional[str] = Field(default=None, description="Category of quantum chemical method")
    core_correlation: Optional[str] = Field(default=None, description="'frozen', 'full', or 'n/a'")
    conformer_index: int = Field(default=1, description="Conformer index (1-based)")
    point_group: Optional[str] = Field(default=None, description="Molecular symmetry point group")
    electronic_state: Optional[str] = Field(default=None, description="Electronic state term symbol")
    detail_url: Optional[str] = Field(default=None, description="Link to calculation detail page (energy3x.asp)")
    provenance: Optional[DataProvenance] = Field(default=None, description="Data provenance metadata")

    def to_unit(self, target_unit: Union[str, EnergyUnit]) -> "EnergyRecord":
        """Return a new :class:`EnergyRecord` with ``value`` converted to ``target_unit``."""
        dst = target_unit.value if isinstance(target_unit, EnergyUnit) else target_unit
        converted = convert_energy(self.value, self.unit, dst)
        return self.model_copy(update={"value": round(converted, 8), "unit": dst})


class EnergyCollection(BaseModel):
    """Collection of :class:`EnergyRecord` items with filtering, unit conversion, and DataFrame export."""

    records: List[EnergyRecord] = Field(default_factory=list)

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self) -> Any:
        return iter(self.records)

    def __getitem__(self, idx: int) -> EnergyRecord:
        return self.records[idx]

    def filter(
        self,
        method: Optional[str] = None,
        basis_set: Optional[str] = None,
        energy_type: Optional[str] = None,
    ) -> "EnergyCollection":
        """Filter energy records by method, basis set, and/or energy_type."""
        from cccbdb.chemistry.basis_sets import resolve_basis_set
        from cccbdb.chemistry.methods import resolve_method

        target_method = None
        if method:
            try:
                target_method = resolve_method(method).name.lower()
            except ValueError:
                target_method = method.strip().lower()

        target_basis = None
        if basis_set:
            try:
                target_basis = resolve_basis_set(basis_set).name.lower()
            except ValueError:
                target_basis = basis_set.strip().lower()

        filtered: List[EnergyRecord] = []
        for rec in self.records:
            if target_method and rec.method.lower() != target_method:
                continue
            if target_basis:
                rec_basis = (rec.basis_set or "Predefined").lower()
                if rec_basis != target_basis:
                    continue
            if energy_type and rec.energy_type.lower() != energy_type.lower():
                continue
            filtered.append(rec)
        return EnergyCollection(records=filtered)

    def convert_units(self, target_unit: Union[str, EnergyUnit]) -> "EnergyCollection":
        """Return a new :class:`EnergyCollection` with all records converted to ``target_unit``."""
        return EnergyCollection(records=[r.to_unit(target_unit) for r in self.records])

    def to_dicts(self) -> List[Dict[str, Any]]:
        """Return a list of flat dictionaries suitable for JSON/CSV/DataFrame export."""
        rows: List[Dict[str, Any]] = []
        for r in self.records:
            rows.append(
                {
                    "species_id": r.species_id,
                    "formula": r.formula,
                    "method": r.method,
                    "basis_set": r.basis_set,
                    "method_category": r.method_category,
                    "core_correlation": r.core_correlation,
                    "energy_type": r.energy_type,
                    "value": r.value,
                    "unit": r.unit,
                    "point_group": r.point_group,
                    "electronic_state": r.electronic_state,
                    "conformer_index": r.conformer_index,
                    "source_url": r.provenance.source_url if r.provenance else r.detail_url,
                    "original_unit": r.provenance.original_unit if r.provenance else r.unit,
                }
            )
        return rows

    def to_dataframe(self) -> Any:
        """Convert energy records into a :class:`pandas.DataFrame`."""
        import pandas as pd

        return pd.DataFrame(self.to_dicts())


class VibrationalMode(BaseModel):
    """Single harmonic or fundamental normal vibrational mode."""

    model_config = ConfigDict(frozen=True)

    mode_number: int = Field(..., description="1-based vibrational mode index")
    symmetry: str = Field(..., description="Irreducible representation / Mulliken symbol (e.g. 'A1', 'B2', 'Σg', 'Πu')")
    frequency_cm_inv: float = Field(
        ...,
        description="Harmonic/calculated frequency in cm^-1 (negative value indicates imaginary transition-state mode)",
    )
    scaled_frequency_cm_inv: Optional[float] = Field(
        default=None, description="Empirically scaled harmonic frequency in cm^-1"
    )
    fundamental_frequency_cm_inv: Optional[float] = Field(
        default=None, description="Experimental fundamental frequency in cm^-1 if experimental"
    )
    ir_intensity_km_mol: Optional[float] = Field(
        default=None, description="Infrared absorption intensity in km/mol"
    )
    raman_activity: Optional[float] = Field(
        default=None, description="Raman scattering activity in Angstrom^4/u"
    )
    depolarization_polar: Optional[float] = Field(default=None, description="Depolarization ratio (polar)")
    depolarization_unpolar: Optional[float] = Field(default=None, description="Depolarization ratio (unpolar)")
    is_imaginary: bool = Field(
        default=False, description="True if mode is an imaginary frequency (written with 'i' in CCCBDB, e.g. '872i')"
    )
    reference: Optional[str] = Field(default=None, description="Literature reference squib if experimental")
    description: Optional[str] = Field(default=None, description="Mode description if provided")


class VibrationalSpectrum(BaseModel):
    """Complete vibrational spectrum and zero-point vibrational energy (ZPE) for a species."""

    species_id: str = Field(..., description="CCCBDB species identifier (CAS RN)")
    formula: str = Field(..., description="Chemical formula")
    method: str = Field(default="Experimental", description="Computational method or 'Experimental'")
    basis_set: Optional[str] = Field(default=None, description="Basis set name if calculated")
    point_group: Optional[str] = Field(default=None, description="Molecular point group")
    conformer_index: int = Field(default=1, description="Conformer index")
    modes: List[VibrationalMode] = Field(default_factory=list, description="List of normal vibrational modes")
    zpe_unscaled_cm_inv: Optional[float] = Field(
        default=None, description="Unscaled zero-point vibrational energy in cm^-1"
    )
    zpe_scaled_cm_inv: Optional[float] = Field(
        default=None, description="Scaled zero-point vibrational energy in cm^-1"
    )
    scaling_factor: Optional[float] = Field(
        default=None, description="Empirical vibrational frequency scaling factor applied"
    )
    provenance: Optional[DataProvenance] = Field(default=None, description="Data provenance metadata")

    @property
    def imaginary_modes(self) -> List[VibrationalMode]:
        """Return all imaginary vibrational modes (saddle points / transition states)."""
        return [m for m in self.modes if m.is_imaginary or m.frequency_cm_inv < 0]

    @property
    def is_transition_state(self) -> bool:
        """Return True if the spectrum has at least one imaginary frequency."""
        return len(self.imaginary_modes) > 0

    def zpe_in(self, unit: Union[str, EnergyUnit] = EnergyUnit.HARTREE, scaled: bool = True) -> Optional[float]:
        """Return the zero-point vibrational energy converted to ``unit`` (e.g. 'Hartree', 'kJ/mol')."""
        val = self.zpe_scaled_cm_inv if (scaled and self.zpe_scaled_cm_inv is not None) else self.zpe_unscaled_cm_inv
        if val is None:
            return None
        return round(convert_energy(val, EnergyUnit.CM_INV, unit), 8)

    def to_dataframe(self) -> Any:
        """Convert vibrational modes into a :class:`pandas.DataFrame`."""
        import pandas as pd

        rows = [
            {
                "mode_number": m.mode_number,
                "symmetry": m.symmetry,
                "frequency_cm_inv": m.frequency_cm_inv,
                "scaled_frequency_cm_inv": m.scaled_frequency_cm_inv,
                "fundamental_frequency_cm_inv": m.fundamental_frequency_cm_inv,
                "ir_intensity_km_mol": m.ir_intensity_km_mol,
                "raman_activity": m.raman_activity,
                "is_imaginary": m.is_imaginary,
                "method": self.method,
                "basis_set": self.basis_set,
            }
            for m in self.modes
        ]
        return pd.DataFrame(rows)


class DipoleMoment(BaseModel):
    """Electric dipole moment vector and magnitude."""

    model_config = ConfigDict(frozen=True)

    species_id: str
    formula: str
    method: str = "Experimental"
    basis_set: Optional[str] = None
    x: Optional[float] = None
    y: Optional[float] = None
    z: Optional[float] = None
    total: float = Field(..., description="Total dipole moment magnitude")
    unit: str = Field(default=DipoleUnit.DEBYE.value, description="Dipole unit (default: 'Debye')")
    point_group: Optional[str] = None
    electronic_state: Optional[str] = None
    reference: Optional[str] = None
    comment: Optional[str] = None
    provenance: Optional[DataProvenance] = None

    def to_unit(self, target_unit: Union[str, DipoleUnit]) -> "DipoleMoment":
        """Convert dipole components and magnitude to ``target_unit`` ('Debye', 'C*m', or 'au')."""
        dst = target_unit.value if isinstance(target_unit, DipoleUnit) else target_unit
        return self.model_copy(
            update={
                "x": round(convert_dipole(self.x, self.unit, dst), 8) if self.x is not None else None,
                "y": round(convert_dipole(self.y, self.unit, dst), 8) if self.y is not None else None,
                "z": round(convert_dipole(self.z, self.unit, dst), 8) if self.z is not None else None,
                "total": round(convert_dipole(self.total, self.unit, dst), 8),
                "unit": dst,
            }
        )


class QuadrupoleMoment(BaseModel):
    """Electric quadrupole moment tensor (primitive or traceless) in Debye*Angstrom."""

    model_config = ConfigDict(frozen=True)

    tensor_type: str = Field(default="Traceless", description="'Primitive', 'Traceless', or 'Experimental'")
    xx: Optional[float] = None
    yy: Optional[float] = None
    zz: Optional[float] = None
    xy: Optional[float] = 0.0
    xz: Optional[float] = 0.0
    yz: Optional[float] = 0.0
    unit: str = "Debye*Angstrom"
    reference: Optional[str] = None


class PolarizabilityTensor(BaseModel):
    """Static electric dipole polarizability tensor and isotropic mean in Angstrom^3."""

    model_config = ConfigDict(frozen=True)

    xx: Optional[float] = None
    yy: Optional[float] = None
    zz: Optional[float] = None
    xy: Optional[float] = 0.0
    xz: Optional[float] = 0.0
    yz: Optional[float] = 0.0
    isotropic: float = Field(..., description="Mean isotropic polarizability alpha = (xx + yy + zz)/3")
    uncertainty: Optional[float] = None
    unit: str = Field(default=PolarizabilityUnit.ANGSTROM3.value)
    reference: Optional[str] = None

    def to_unit(self, target_unit: Union[str, PolarizabilityUnit]) -> "PolarizabilityTensor":
        """Convert polarizability tensor to ``target_unit`` ('Angstrom^3' or 'Bohr^3')."""
        dst = target_unit.value if isinstance(target_unit, PolarizabilityUnit) else target_unit
        conv = lambda v: round(convert_polarizability(v, self.unit, dst), 6) if v is not None else None
        return self.model_copy(
            update={
                "xx": conv(self.xx),
                "yy": conv(self.yy),
                "zz": conv(self.zz),
                "xy": conv(self.xy),
                "xz": conv(self.xz),
                "yz": conv(self.yz),
                "isotropic": round(convert_polarizability(self.isotropic, self.unit, dst), 6),
                "unit": dst,
            }
        )


class AtomicChargeRecord(BaseModel):
    """Partial atomic charges (Mulliken, CHELPG, AIM, ESP) for an atom in a molecule."""

    model_config = ConfigDict(frozen=True)

    atom_index: int
    element: str
    mulliken: Optional[float] = None
    chelpg: Optional[float] = None
    aim: Optional[float] = None
    esp: Optional[float] = None
    unit: str = "e"


class RotationalConstants(BaseModel):
    """Rotational constants (A, B, C) for a symmetric, asymmetric, spherical, or linear rotor."""

    model_config = ConfigDict(frozen=True)

    A: Optional[float] = None
    B: Optional[float] = None
    C: Optional[float] = None
    unit: str = Field(default=RotationalConstantUnit.CM_INV.value)
    reference: Optional[str] = None
    comment: Optional[str] = None

    def to_unit(self, target_unit: Union[str, RotationalConstantUnit]) -> "RotationalConstants":
        """Convert rotational constants between cm^-1, GHz, and MHz."""
        dst = target_unit.value if isinstance(target_unit, RotationalConstantUnit) else target_unit
        conv = lambda v: round(convert_rotational_constant(v, self.unit, dst), 6) if v is not None else None
        return self.model_copy(
            update={"A": conv(self.A), "B": conv(self.B), "C": conv(self.C), "unit": dst}
        )


class ThermochemicalProperty(BaseModel):
    """Experimental or calculated thermochemical datum (enthalpy of formation, entropy, heat capacity)."""

    model_config = ConfigDict(frozen=True)

    property_name: str = Field(..., description="e.g., 'Hfg(298.15K)', 'Hfg(0K)', 'Entropy (298.15K)', 'Heat Capacity (298.15K)'")
    category: str = Field(..., description="'enthalpy_of_formation', 'entropy', 'integrated_heat_capacity', 'heat_capacity'")
    temperature_k: Optional[float] = Field(default=298.15, description="Reference temperature in Kelvin")
    value: float = Field(..., description="Numerical property value")
    uncertainty: Optional[float] = Field(default=None, description="Experimental uncertainty")
    unit: str = Field(..., description="Physical unit (e.g., 'kJ/mol', 'J/(K*mol)')")
    reference: Optional[str] = Field(default=None, description="Reference squib (e.g., 'CODATA', 'Gurvich')")
    comment: Optional[str] = Field(default=None)
    provenance: Optional[DataProvenance] = None
