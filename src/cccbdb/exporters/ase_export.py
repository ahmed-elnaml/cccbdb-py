"""Optional Atomic Simulation Environment (ASE) and RDKit structure adapters."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict

from cccbdb.chemistry.units import LengthUnit
from cccbdb.exceptions import IncompleteGeometryError

if TYPE_CHECKING:
    from cccbdb.models.geometry import MolecularGeometry


def to_ase_atoms(geometry: "MolecularGeometry") -> Any:
    """Convert a :class:`MolecularGeometry` into an :class:`ase.Atoms` object (if ASE is installed).

    If `ase` is not installed in the active Python environment, raises an informative
    `ImportError` instructing the user to install `cccbdb-py[chemistry]`.
    """
    if not geometry.has_complete_cartesian:
        raise IncompleteGeometryError(
            f"Cannot convert '{geometry.formula}' to ASE Atoms without 3D Cartesian coordinates."
        )
    geom_ang = (
        geometry.convert_units(LengthUnit.ANGSTROM)
        if geometry.unit != LengthUnit.ANGSTROM.value
        else geometry
    )
    symbols = [ "H" if s == "D" else s for s in geom_ang.atoms ]
    positions = geom_ang.coordinates
    info: Dict[str, Any] = {
        "species_id": geom_ang.species_id,
        "formula": geom_ang.formula,
        "method": geom_ang.method,
        "basis_set": geom_ang.basis_set,
        "point_group": geom_ang.point_group,
    }
    try:
        from ase import Atoms  # type: ignore[import-not-found]

        atoms = Atoms(symbols=symbols, positions=positions, info=info)
        # Set deuterium masses explicitly if present
        for idx, raw_sym in enumerate(geom_ang.atoms):
            if raw_sym == "D":
                atoms[idx].mass = 2.014101778
        return atoms
    except ImportError as exc:
        raise ImportError(
            "ASE (Atomic Simulation Environment) is required for to_ase(). "
            "Install it via `pip install cccbdb-py[chemistry]` or `pip install ase`."
        ) from exc
