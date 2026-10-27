"""Parsers for experimental and calculated molecular geometries (Cartesians and internal coordinates)."""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple, Union

from bs4 import BeautifulSoup, Tag

from cccbdb.chemistry.symmetry import normalize_point_group_symbol
from cccbdb.chemistry.units import LengthUnit
from cccbdb.models.geometry import AtomCoordinate, InternalCoordinate, MolecularGeometry
from cccbdb.models.provenance import DataProvenance
from cccbdb.parsers._utils import clean_text, make_soup, parse_float_or_none


def _parse_atom_label(raw_label: str, fallback_index: int) -> Tuple[str, int, bool]:
    """Split a CCCBDB atom label like ``C1``, ``H13``, or ``X14`` into ``(symbol, index, is_dummy)``."""
    cleaned = raw_label.strip()
    m = re.match(r"^([A-Za-z]+)(\d+)?$", cleaned)
    if m:
        sym = m.group(1).capitalize()
        idx = int(m.group(2)) if m.group(2) else fallback_index
        is_dummy = sym.upper() == "X"
        return sym, idx, is_dummy
    return cleaned, fallback_index, (cleaned.upper() == "X")


def _extract_cartesian_tables(soup: BeautifulSoup) -> List[List[AtomCoordinate]]:
    """Find all Cartesian coordinate tables (``Atom | x (Å) | y (Å) | z (Å)``) in the document."""
    all_tables: List[List[AtomCoordinate]] = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        # Match either ['atom', 'x (å)', 'y (å)', 'z (å)'] or ['atom', 'internal', 'x (å)', 'y (å)', 'z (å)']
        if "atom" in headers and any("x" in h for h in headers) and any("z" in h for h in headers):
            atom_col = headers.index("atom")
            x_col = next(i for i, h in enumerate(headers) if "x" in h)
            y_col = next(i for i, h in enumerate(headers) if "y" in h)
            z_col = next(i for i, h in enumerate(headers) if "z" in h)

            sites: List[AtomCoordinate] = []
            for r_idx, r in enumerate(rows[1:], start=1):
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) <= max(atom_col, x_col, y_col, z_col):
                    continue
                lbl = cols[atom_col]
                x = parse_float_or_none(cols[x_col])
                y = parse_float_or_none(cols[y_col])
                z = parse_float_or_none(cols[z_col])
                if lbl and x is not None and y is not None and z is not None:
                    sym, idx, is_dummy = _parse_atom_label(lbl, r_idx)
                    sites.append(
                        AtomCoordinate(
                            index=idx,
                            label=lbl,
                            symbol=sym,
                            x=x,
                            y=y,
                            z=z,
                            is_dummy=is_dummy,
                        )
                    )
            if sites:
                all_tables.append(sites)
    return all_tables


def _extract_internal_coordinates(soup: BeautifulSoup) -> List[InternalCoordinate]:
    """Extract experimental internal coordinates (``rCO``, ``rOH``, ``aHOH``) and calculated bond angles."""
    internals: List[InternalCoordinate] = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]

        # Format 1: Experimental internal coordinates table (Description | Value | unc. | Connectivity... | Reference | Comment)
        if "description" in headers and "value" in headers:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) < 2:
                    continue
                desc = cols[0]
                if desc.lower() == "description":
                    continue
                val = parse_float_or_none(cols[1])
                if val is None:
                    continue
                unc = parse_float_or_none(cols[2]) if len(cols) > 2 else None
                atoms_list: List[Union[int, str]] = []
                for a_col in cols[3:7]:
                    if a_col:
                        atoms_list.append(int(a_col) if a_col.isdigit() else a_col)
                ref = cols[7] if len(cols) > 7 and cols[7] else None
                comment = cols[8] if len(cols) > 8 and cols[8] else None
                if desc.startswith("r"):
                    c_type, unit = "distance", LengthUnit.ANGSTROM.value
                elif desc.startswith("a"):
                    c_type, unit = "angle", "degrees"
                else:
                    c_type, unit = "dihedral", "degrees"
                internals.append(
                    InternalCoordinate(
                        name=desc,
                        coord_type=c_type,
                        value=val,
                        uncertainty=unc,
                        unit=unit,
                        atoms=tuple(atoms_list),
                        reference=ref,
                        comment=comment,
                    )
                )

        # Format 2: Calculated Bond Angles table (atom1 | atom2 | atom3 | angle | [blank] | atom1 | atom2 | atom3 | angle)
        if "atom1" in headers and "atom2" in headers and "angle" in headers:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                # Check left group (cols 0..3) and optional right group (cols 4..7 or 5..8)
                for offset in (0, 4, 5):
                    if offset + 3 < len(cols):
                        a1, a2, a3, ang_str = cols[offset : offset + 4]
                        ang = parse_float_or_none(ang_str)
                        if a1 and a2 and a3 and ang is not None:
                            internals.append(
                                InternalCoordinate(
                                    name=f"a({a1}-{a2}-{a3})",
                                    coord_type="angle",
                                    value=ang,
                                    unit="degrees",
                                    atoms=(a1, a2, a3),
                                )
                            )
    return internals


def parse_experimental_geometry(
    soup: BeautifulSoup,
    species_id: str,
    formula: str,
    name: Optional[str] = None,
    point_group: Optional[str] = None,
    electronic_state: Optional[str] = None,
    provenance: Optional[DataProvenance] = None,
) -> Optional[MolecularGeometry]:
    """Extract the experimental molecular geometry from an ``exp2x.asp`` or ``expgeom2x.asp`` BeautifulSoup tree."""
    cart_tables = _extract_cartesian_tables(soup)
    internals = _extract_internal_coordinates(soup)
    sites = cart_tables[0] if cart_tables else []
    if not sites and not internals:
        return None
    return MolecularGeometry(
        species_id=species_id,
        formula=formula,
        name=name,
        method="Experimental",
        basis_set=None,
        point_group=point_group,
        electronic_state=electronic_state,
        conformer_index=1,
        is_minimum_conformation=True,
        unit=LengthUnit.ANGSTROM.value,
        sites=sites,
        internal_coordinates=internals,
        provenance=provenance,
    )


def parse_calculated_geometry_section(
    soup: BeautifulSoup,
    species_id: str,
    formula: str,
    name: Optional[str],
    method: str,
    basis_set: Optional[str],
    point_group: Optional[str],
    electronic_state: Optional[str],
    conformer_index: int = 1,
    is_minimum: bool = True,
    occurrence_index: int = 0,
    provenance: Optional[DataProvenance] = None,
) -> Optional[MolecularGeometry]:
    """Extract the ``occurrence_index``-th calculated geometry section from ``energy3x.asp`` or ``geom3x.asp``."""
    cart_tables = _extract_cartesian_tables(soup)
    internals = _extract_internal_coordinates(soup)

    sites = cart_tables[occurrence_index] if occurrence_index < len(cart_tables) else []

    if not point_group:
        pg_matches = re.findall(r"Point Group is\s+([A-Za-z0-9\*∞]+)", soup.get_text(" ", strip=True))
        if occurrence_index < len(pg_matches):
            try:
                point_group = normalize_point_group_symbol(pg_matches[occurrence_index])
            except ValueError:
                point_group = pg_matches[occurrence_index]

    if not sites and not internals:
        return None

    return MolecularGeometry(
        species_id=species_id,
        formula=formula,
        name=name,
        method=method,
        basis_set=basis_set,
        point_group=point_group,
        electronic_state=electronic_state,
        conformer_index=conformer_index,
        is_minimum_conformation=is_minimum,
        unit=LengthUnit.ANGSTROM.value,
        sites=sites,
        internal_coordinates=internals if occurrence_index == 0 else [],
        provenance=provenance,
    )
