"""Complete single-calculation summary model (representing NIST CCCBDB energy3x.asp)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

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


from cccbdb.models.geometry import MolecularGeometry
from cccbdb.models.property import (
    AtomicChargeRecord,
    DipoleMoment,
    EnergyRecord,
    PolarizabilityTensor,
    QuadrupoleMoment,
    RotationalConstants,
    VibrationalSpectrum,
)
from cccbdb.models.provenance import DataProvenance


class ConformerCalculation(BaseModel):
    """All computed quantum chemical properties for a single conformer/electronic state."""

    conformer_index: int = 1
    state_index: int = 1
    is_minimum: bool = True
    point_group: Optional[str] = None
    electronic_state: Optional[str] = None
    energies: Dict[str, EnergyRecord] = Field(default_factory=dict)
    geometry: Optional[MolecularGeometry] = None
    vibrations: Optional[VibrationalSpectrum] = None
    rotational_constants: Optional[RotationalConstants] = None
    atomic_charges: List[AtomicChargeRecord] = Field(default_factory=list)
    dipole: Optional[DipoleMoment] = None
    quadrupoles: Dict[str, QuadrupoleMoment] = Field(default_factory=dict)
    polarizability: Optional[PolarizabilityTensor] = None
    r2_angstrom2: Optional[float] = Field(default=None, description="<r^2> spatial extent in Angstrom^2")


class CalculationResult(BaseModel):
    """Complete result from a given model chemistry calculation (method + basis set) in NIST CCCBDB."""

    species_id: str
    formula: str
    name: Optional[str] = None
    method: str
    basis_set: Optional[str] = None
    method_id: Optional[int] = None
    basis_id: Optional[int] = None
    conformers: List[ConformerCalculation] = Field(default_factory=list)
    provenance: Optional[DataProvenance] = None

    @property
    def primary(self) -> Optional[ConformerCalculation]:
        """Return the minimum-energy or primary conformer calculation."""
        for c in self.conformers:
            if c.is_minimum:
                return c
        return self.conformers[0] if self.conformers else None

    @property
    def electronic_energy(self) -> Optional[float]:
        """Return the 0K electronic energy (in Hartree) of the primary conformer."""
        p = self.primary
        if not p:
            return None
        rec = p.energies.get("Energy at 0K") or p.energies.get("electronic_0K") or p.energies.get("HF Energy")
        return rec.value if rec else None

    @property
    def geometry(self) -> Optional[MolecularGeometry]:
        """Return the optimized geometry of the primary conformer."""
        p = self.primary
        return p.geometry if p else None

    @property
    def vibrations(self) -> Optional[VibrationalSpectrum]:
        """Return the vibrational spectrum of the primary conformer."""
        p = self.primary
        return p.vibrations if p else None

    @property
    def dipole(self) -> Optional[DipoleMoment]:
        """Return the electric dipole moment of the primary conformer."""
        p = self.primary
        return p.dipole if p else None
