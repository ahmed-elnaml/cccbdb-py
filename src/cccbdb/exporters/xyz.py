"""Standard XMol (.xyz) molecular structure exporter with Cartesian completeness validation."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from cccbdb.chemistry.units import LengthUnit
from cccbdb.exceptions import IncompleteGeometryError

if TYPE_CHECKING:
    from cccbdb.models.geometry import MolecularGeometry


def format_xyz(geometry: "MolecularGeometry", comment: Optional[str] = None) -> str:
    """Format a :class:`MolecularGeometry` into an XMol `.xyz` string in Angstroms.

    Strictly validates that unambiguous 3D Cartesian coordinates exist and excludes
    any non-physical dummy atoms (`X`).
    """
    if not geometry.has_complete_cartesian:
        raise IncompleteGeometryError(
            f"Cannot export XYZ for '{geometry.formula}' ({geometry.method}/{geometry.basis_set}): "
            "complete 3D Cartesian coordinates are not available (internal coordinates alone do not "
            "uniquely orient the structure without a Z-matrix conversion)."
        )

    geom_ang = (
        geometry.convert_units(LengthUnit.ANGSTROM)
        if geometry.unit != LengthUnit.ANGSTROM.value
        else geometry
    )
    phys = geom_ang.physical_sites
    if comment is None:
        model_str = (
            f"{geom_ang.method}/{geom_ang.basis_set}"
            if geom_ang.basis_set
            else geom_ang.method
        )
        pg_str = f" | point_group={geom_ang.point_group}" if geom_ang.point_group else ""
        src_str = (
            f" | source={geom_ang.provenance.source_url}"
            if geom_ang.provenance
            else ""
        )
        comment = f"{geom_ang.formula} ({geom_ang.name or geom_ang.species_id}) | {model_str}{pg_str} | unit=Angstrom{src_str}"

    lines = [str(len(phys)), comment]
    for site in phys:
        lines.append(f"{site.symbol:<3s} {site.x:14.8f} {site.y:14.8f} {site.z:14.8f}")
    return "\n".join(lines) + "\n"


def write_xyz(
    geometry: "MolecularGeometry",
    path: Union[str, Path],
    comment: Optional[str] = None,
) -> Path:
    """Write a :class:`MolecularGeometry` to an `.xyz` file on disk."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    xyz_text = format_xyz(geometry, comment=comment)
    out_path.write_text(xyz_text, encoding="utf-8")
    return out_path
