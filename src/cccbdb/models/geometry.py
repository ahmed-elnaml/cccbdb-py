"""Molecular geometry models (Cartesian coordinates, internal coordinates, distance matrices)."""

from __future__ import annotations

import math
from pathlib import Path
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


from cccbdb.chemistry.formula import ATOMIC_WEIGHTS
from cccbdb.chemistry.symmetry import PointGroupInfo, normalize_point_group_symbol, parse_point_group
from cccbdb.chemistry.units import LengthUnit, convert_length
from cccbdb.exceptions import IncompleteGeometryError
from cccbdb.models.provenance import DataProvenance


class AtomCoordinate(BaseModel):
    """Single atomic site in 3D Cartesian space."""

    model_config = ConfigDict(frozen=True)

    index: int = Field(..., description="1-based atom index in the molecule (e.g., 1 for C1)")
    label: str = Field(..., description="Atom label from CCCBDB (e.g., 'C1', 'H2', 'X14')")
    symbol: str = Field(..., description="Chemical element symbol (e.g., 'C', 'H', 'D', or 'X' for dummy)")
    x: float = Field(..., description="Cartesian x coordinate")
    y: float = Field(..., description="Cartesian y coordinate")
    z: float = Field(..., description="Cartesian z coordinate")
    is_dummy: bool = Field(default=False, description="True if site is a dummy atom ('X')")


class InternalCoordinate(BaseModel):
    """Internal coordinate (bond length, bond angle, or dihedral angle)."""

    model_config = ConfigDict(frozen=True)

    name: str = Field(..., description="Internal coordinate label (e.g., 'rCO', 'rOH', 'aHOH')")
    coord_type: str = Field(..., description="'distance', 'angle', or 'dihedral'")
    value: float = Field(..., description="Numerical value (Angstrom for distance, degrees for angle/dihedral)")
    uncertainty: Optional[float] = Field(default=None, description="Experimental uncertainty if reported")
    unit: str = Field(..., description="Unit ('Angstrom' or 'degrees')")
    atoms: Tuple[Union[int, str], ...] = Field(
        default=(), description="Connected atom indices or labels defining the coordinate"
    )
    reference: Optional[str] = Field(default=None, description="Reference squib")
    comment: Optional[str] = Field(default=None, description="Additional comment (e.g., 're', 'r0', 'rs')")


class MolecularGeometry(BaseModel):
    """Complete representation of an experimental or calculated molecular geometry."""

    species_id: str = Field(..., description="CCCBDB species identifier (CAS RN)")
    formula: str = Field(..., description="Chemical formula")
    name: Optional[str] = Field(default=None, description="Chemical species name")
    method: str = Field(default="Experimental", description="Computational method or 'Experimental'")
    basis_set: Optional[str] = Field(default=None, description="Basis set name if calculated")
    point_group: Optional[str] = Field(default=None, description="Normalized Schoenflies point group")
    electronic_state: Optional[str] = Field(default=None, description="Electronic state symmetry (e.g. '1A1', '1Σ')")
    conformer_index: int = Field(default=1, description="Conformer number (1-indexed)")
    is_minimum_conformation: bool = Field(default=True, description="Whether this is the global minimum conformer")
    unit: str = Field(default=LengthUnit.ANGSTROM.value, description="Coordinate length unit (default: Angstrom)")
    sites: List[AtomCoordinate] = Field(default_factory=list, description="All atomic sites including any dummy atoms")
    internal_coordinates: List[InternalCoordinate] = Field(
        default_factory=list, description="Internal coordinates (bonds, angles, dihedrals)"
    )
    rotational_constants_cm_inv: Dict[str, float] = Field(
        default_factory=dict, description="Rotational constants {'A': ..., 'B': ..., 'C': ...} in cm^-1"
    )
    provenance: Optional[DataProvenance] = Field(default=None, description="Data provenance metadata")

    @property
    def physical_sites(self) -> List[AtomCoordinate]:
        """Return all real physical atoms, excluding dummy ('X') sites."""
        return [s for s in self.sites if not s.is_dummy and s.symbol.upper() != "X"]

    @property
    def atoms(self) -> List[str]:
        """Return list of chemical element symbols for all physical atoms."""
        return [s.symbol for s in self.physical_sites]

    @property
    def coordinates(self) -> List[Tuple[float, float, float]]:
        """Return list of (x, y, z) Cartesian coordinate tuples for all physical atoms."""
        return [(s.x, s.y, s.z) for s in self.physical_sites]

    @property
    def has_complete_cartesian(self) -> bool:
        """Return True if unambiguous 3D Cartesian coordinates exist for at least one physical atom."""
        return len(self.physical_sites) > 0

    @property
    def symmetry_info(self) -> Optional[PointGroupInfo]:
        """Return structured Schoenflies point group details if point_group is known."""
        if not self.point_group:
            return None
        try:
            return parse_point_group(self.point_group)
        except ValueError:
            return None

    def convert_units(self, target_unit: Union[str, LengthUnit]) -> "MolecularGeometry":
        """Return a new :class:`MolecularGeometry` with Cartesian coordinates converted to ``target_unit``."""
        dst = LengthUnit(target_unit).value if isinstance(target_unit, LengthUnit) else target_unit
        new_sites = [
            AtomCoordinate(
                index=s.index,
                label=s.label,
                symbol=s.symbol,
                x=round(convert_length(s.x, self.unit, dst), 8),
                y=round(convert_length(s.y, self.unit, dst), 8),
                z=round(convert_length(s.z, self.unit, dst), 8),
                is_dummy=s.is_dummy,
            )
            for s in self.sites
        ]
        return self.model_copy(update={"unit": dst, "sites": new_sites})

    def to_numpy(self) -> Any:
        """Return an (N, 3) NumPy float64 array of physical Cartesian coordinates."""
        import numpy as np

        if not self.has_complete_cartesian:
            raise IncompleteGeometryError(
                f"Geometry for '{self.formula}' ({self.method}/{self.basis_set}) "
                "does not contain complete 3D Cartesian coordinates."
            )
        return np.array(self.coordinates, dtype=np.float64)

    def distance(self, i: int, j: int) -> float:
        """Compute Euclidean distance between physical atom indices ``i`` and ``j`` (0-indexed)."""
        phys = self.physical_sites
        a, b = phys[i], phys[j]
        return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)

    def distance_matrix(self) -> List[List[float]]:
        """Compute full N x N interatomic distance matrix for all physical atoms."""
        phys = self.physical_sites
        n = len(phys)
        mat: List[List[float]] = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d = round(self.distance(i, j), 6)
                mat[i][j] = d
                mat[j][i] = d
        return mat

    def center_of_mass(self) -> Tuple[float, float, float]:
        """Compute the mass-weighted center of mass (x, y, z) of the physical atoms."""
        phys = self.physical_sites
        if not phys:
            raise IncompleteGeometryError("Cannot compute center of mass without Cartesian coordinates.")
        total_mass = 0.0
        cx, cy, cz = 0.0, 0.0, 0.0
        for site in phys:
            mass = ATOMIC_WEIGHTS.get(site.symbol, 12.0)
            total_mass += mass
            cx += mass * site.x
            cy += mass * site.y
            cz += mass * site.z
        return (round(cx / total_mass, 6), round(cy / total_mass, 6), round(cz / total_mass, 6))

    def to_xyz_string(self, comment: Optional[str] = None) -> str:
        """Format geometry as a standard XMol `.xyz` string in Angstroms.

        Raises :class:`IncompleteGeometryError` if full 3D Cartesian coordinates are unavailable.
        """
        from cccbdb.exporters.xyz import format_xyz

        return format_xyz(self, comment=comment)

    def to_xyz(self, path: Union[str, Path], comment: Optional[str] = None) -> Path:
        """Write geometry to a standard `.xyz` file on disk."""
        from cccbdb.exporters.xyz import write_xyz

        return write_xyz(self, path, comment=comment)

    def to_ase(self) -> Any:
        """Convert geometry to an Atomic Simulation Environment (:class:`ase.Atoms`) object."""
        from cccbdb.exporters.ase_export import to_ase_atoms

        return to_ase_atoms(self)

    def to_dataframe(self) -> Any:
        """Convert physical atomic coordinates into a :class:`pandas.DataFrame`."""
        import pandas as pd

        rows = [
            {
                "atom_index": s.index,
                "label": s.label,
                "symbol": s.symbol,
                "x": s.x,
                "y": s.y,
                "z": s.z,
                "unit": self.unit,
                "method": self.method,
                "basis_set": self.basis_set,
                "point_group": self.point_group,
            }
            for s in self.physical_sites
        ]
        return pd.DataFrame(rows)
