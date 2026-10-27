"""Parsers for thermochemistry, dipole moments, quadrupole tensors, polarizability, charges, and rotational constants."""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from bs4 import BeautifulSoup

from cccbdb.chemistry.symmetry import normalize_point_group_symbol
from cccbdb.chemistry.units import DipoleUnit, PolarizabilityUnit, RotationalConstantUnit
from cccbdb.models.property import (
    AtomicChargeRecord,
    DipoleMoment,
    PolarizabilityTensor,
    QuadrupoleMoment,
    RotationalConstants,
    ThermochemicalProperty,
)
from cccbdb.models.provenance import DataProvenance
from cccbdb.parsers._utils import clean_text, parse_float_or_none


def parse_experimental_thermochemistry(
    soup: BeautifulSoup,
    provenance: Optional[DataProvenance] = None,
) -> List[ThermochemicalProperty]:
    """Parse the experimental thermochemistry table (``Hfg(298.15K)``, ``Hfg(0K)``, ``Entropy``, ``Heat Capacity``)."""
    props: List[ThermochemicalProperty] = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if "property" in headers and "value" in headers and "units" in headers:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) < 4:
                    continue
                prop_name = cols[0]
                val = parse_float_or_none(cols[1])
                if val is None:
                    continue
                unc = parse_float_or_none(cols[2])
                unit_raw = cols[3].replace("mol -1", "mol^-1").replace("K -1", "K^-1")
                ref = cols[4] if len(cols) > 4 and cols[4] else None
                comment = cols[5] if len(cols) > 5 and cols[5] else None

                lower_p = prop_name.lower()
                if "hfg" in lower_p or "enthalpy" in lower_p:
                    cat = "enthalpy_of_formation"
                elif "entropy" in lower_p:
                    cat = "entropy"
                elif "integrated" in lower_p:
                    cat = "integrated_heat_capacity"
                else:
                    cat = "heat_capacity"

                temp = 0.0 if "0k" in lower_p.replace(" ", "") and "298" not in lower_p else 298.15
                props.append(
                    ThermochemicalProperty(
                        property_name=prop_name,
                        category=cat,
                        temperature_k=temp,
                        value=val,
                        uncertainty=unc,
                        unit=unit_raw,
                        reference=ref,
                        comment=comment,
                        provenance=provenance,
                    )
                )
    return props


def parse_experimental_rotational_constants(soup: BeautifulSoup) -> Optional[RotationalConstants]:
    """Parse experimental rotational constants (A, B, C in cm^-1) from ``exp2x.asp``."""
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c) for c in rows[0].find_all(["th", "td"])]
        if headers[:3] == ["A", "B", "C"]:
            cols = [clean_text(c) for c in rows[1].find_all(["td", "th"])]
            a_val = parse_float_or_none(cols[0]) if len(cols) > 0 else None
            b_val = parse_float_or_none(cols[1]) if len(cols) > 1 else None
            c_val = parse_float_or_none(cols[2]) if len(cols) > 2 else None
            ref = cols[3] if len(cols) > 3 and cols[3] else None
            comment = cols[4] if len(cols) > 4 and cols[4] else None
            if any(v is not None for v in (a_val, b_val, c_val)):
                return RotationalConstants(
                    A=a_val,
                    B=b_val,
                    C=c_val,
                    unit=RotationalConstantUnit.CM_INV.value,
                    reference=ref,
                    comment=comment,
                )
    return None


def parse_calculated_rotational_constants(
    soup: BeautifulSoup,
    occurrence_index: int = 0,
) -> Optional[RotationalConstants]:
    """Parse calculated rotational constants (A, B, C or linear B in cm^-1) from ``energy3x.asp``."""
    rot_tables = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) != 2:
            continue
        headers = [clean_text(c) for c in rows[0].find_all(["th", "td"])]
        if headers in (["A", "B", "C"], ["B"]):
            rot_tables.append((headers, rows[1]))

    if occurrence_index >= len(rot_tables):
        return None
    headers, val_row = rot_tables[occurrence_index]
    cols = [clean_text(c) for c in val_row.find_all(["td", "th"])]
    if headers == ["B"] and cols:
        b_val = parse_float_or_none(cols[0])
        return RotationalConstants(A=None, B=b_val, C=None, unit=RotationalConstantUnit.CM_INV.value)
    if headers == ["A", "B", "C"] and len(cols) >= 3:
        return RotationalConstants(
            A=parse_float_or_none(cols[0]),
            B=parse_float_or_none(cols[1]),
            C=parse_float_or_none(cols[2]),
            unit=RotationalConstantUnit.CM_INV.value,
        )
    return None


def parse_experimental_dipole(
    soup: BeautifulSoup,
    species_id: str,
    formula: str,
    point_group: Optional[str] = None,
    provenance: Optional[DataProvenance] = None,
) -> Optional[DipoleMoment]:
    """Parse the experimental electric dipole moment table from ``exp2x.asp``."""
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        header_text = " ".join(clean_text(r).lower() for r in rows[:2])
        if ("dipole" in header_text or "state description" in header_text) and "total" in header_text:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) < 9 or not cols[0].isdigit():
                    continue
                el_state = cols[2] or None
                pg_raw = cols[3] or point_group
                try:
                    pg = normalize_point_group_symbol(pg_raw) if pg_raw else point_group
                except ValueError:
                    pg = pg_raw
                x = parse_float_or_none(cols[5])
                y = parse_float_or_none(cols[6])
                z = parse_float_or_none(cols[7])
                tot = parse_float_or_none(cols[8])
                ref = cols[9] if len(cols) > 9 and cols[9] else None
                comment = cols[10] if len(cols) > 10 and cols[10] else None
                if tot is not None:
                    return DipoleMoment(
                        species_id=species_id,
                        formula=formula,
                        method="Experimental",
                        basis_set=None,
                        x=x,
                        y=y,
                        z=z,
                        total=tot,
                        unit=DipoleUnit.DEBYE.value,
                        point_group=pg,
                        electronic_state=el_state,
                        reference=ref,
                        comment=comment,
                        provenance=provenance,
                    )
    return None


def parse_calculated_dipole_table(
    soup: BeautifulSoup,
    species_id: str,
    formula: str,
    method: str,
    basis_set: Optional[str],
    point_group: Optional[str] = None,
    electronic_state: Optional[str] = None,
    occurrence_index: int = 0,
    provenance: Optional[DataProvenance] = None,
) -> Optional[DipoleMoment]:
    """Parse the ``occurrence_index``-th calculated electric dipole moment table from ``energy3x.asp``."""
    dip_tables = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if "x" in headers and "y" in headers and "z" in headers and "total" in headers:
            dip_tables.append(table)

    if occurrence_index >= len(dip_tables):
        return None

    dtable = dip_tables[occurrence_index]
    rows = dtable.find_all("tr")
    for r in rows[1:]:
        cols = [clean_text(c) for c in r.find_all(["td", "th"])]
        if len(cols) >= 4:
            # Could be [x, y, z, Total] or [label, x, y, z, Total]
            offset = 1 if len(cols) >= 5 and not parse_float_or_none(cols[0]) and cols[0] != "0.000" else 0
            if len(cols) >= offset + 4:
                x = parse_float_or_none(cols[offset])
                y = parse_float_or_none(cols[offset + 1])
                z = parse_float_or_none(cols[offset + 2])
                tot = parse_float_or_none(cols[offset + 3])
                if tot is not None:
                    return DipoleMoment(
                        species_id=species_id,
                        formula=formula,
                        method=method,
                        basis_set=basis_set,
                        x=x,
                        y=y,
                        z=z,
                        total=tot,
                        unit=DipoleUnit.DEBYE.value,
                        point_group=point_group,
                        electronic_state=electronic_state,
                        provenance=provenance,
                    )
    return None


def parse_atomic_charges_table(
    soup: BeautifulSoup,
    occurrence_index: int = 0,
) -> List[AtomicChargeRecord]:
    """Parse Mulliken, CHELPG, AIM, and ESP partial atomic charges from ``energy3x.asp``."""
    charge_tables = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if "element" in headers and "mulliken" in headers:
            charge_tables.append(table)

    if occurrence_index >= len(charge_tables):
        return []

    ctable = charge_tables[occurrence_index]
    records: List[AtomicChargeRecord] = []
    for r in ctable.find_all("tr")[1:]:
        cols = [clean_text(c) for c in r.find_all(["td", "th"])]
        if len(cols) >= 3 and cols[0].isdigit():
            records.append(
                AtomicChargeRecord(
                    atom_index=int(cols[0]),
                    element=cols[1],
                    mulliken=parse_float_or_none(cols[2]) if len(cols) > 2 else None,
                    chelpg=parse_float_or_none(cols[3]) if len(cols) > 3 else None,
                    aim=parse_float_or_none(cols[4]) if len(cols) > 4 else None,
                    esp=parse_float_or_none(cols[5]) if len(cols) > 5 else None,
                )
            )
    return records


def parse_experimental_quadrupoles(soup: BeautifulSoup) -> List[QuadrupoleMoment]:
    """Parse experimental electric quadrupole moment rows from ``exp2x.asp``."""
    quads: List[QuadrupoleMoment] = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        header_text = " ".join(clean_text(r).lower() for r in rows[:2])
        if "quadrupole" in header_text and "xx" in header_text and "zz" in header_text:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) >= 8 and cols[0].isdigit():
                    xx = parse_float_or_none(cols[5])
                    yy = parse_float_or_none(cols[6])
                    zz = parse_float_or_none(cols[7])
                    ref = cols[8] if len(cols) > 8 and cols[8] else None
                    if any(v is not None for v in (xx, yy, zz)):
                        quads.append(
                            QuadrupoleMoment(
                                tensor_type="Experimental",
                                xx=xx,
                                yy=yy,
                                zz=zz,
                                reference=ref,
                            )
                        )
    return quads


def parse_calculated_quadrupole_tables(
    soup: BeautifulSoup,
    occurrence_index: int = 0,
) -> Dict[str, QuadrupoleMoment]:
    """Parse Primitive and Traceless 3x3 quadrupole moment tensors from ``energy3x.asp``."""
    primitive_tables = []
    traceless_tables = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) != 4:
            continue
        h_text = clean_text(rows[0]).lower()
        r1 = [clean_text(c).lower() for c in rows[1].find_all(["td", "th"])]
        if "x" in h_text and "y" in h_text and "z" in h_text and r1 and r1[0] == "x":
            # Distinguish primitive vs traceless vs polarizability by preceding text or header
            prev_text = clean_text(table.find_previous(string=re.compile(r"Primitive|Traceless|Units are|Angstrom", re.I)))
            if "primitive" in h_text or "primitive" in prev_text.lower():
                primitive_tables.append(table)
            elif "traceless" in h_text or "traceless" in prev_text.lower():
                traceless_tables.append(table)

    result: Dict[str, QuadrupoleMoment] = {}
    for label, t_list in [("Primitive", primitive_tables), ("Traceless", traceless_tables)]:
        if occurrence_index < len(t_list):
            tbl = t_list[occurrence_index]
            rows = tbl.find_all("tr")
            r_x = [clean_text(c) for c in rows[1].find_all(["td", "th"])]
            r_y = [clean_text(c) for c in rows[2].find_all(["td", "th"])]
            r_z = [clean_text(c) for c in rows[3].find_all(["td", "th"])]
            if len(r_x) >= 4 and len(r_y) >= 4 and len(r_z) >= 4:
                result[label] = QuadrupoleMoment(
                    tensor_type=label,
                    xx=parse_float_or_none(r_x[1]),
                    xy=parse_float_or_none(r_x[2]),
                    xz=parse_float_or_none(r_x[3]),
                    yy=parse_float_or_none(r_y[2]),
                    yz=parse_float_or_none(r_y[3]),
                    zz=parse_float_or_none(r_z[3]),
                )
    return result


def parse_experimental_polarizability(soup: BeautifulSoup) -> Optional[PolarizabilityTensor]:
    """Parse experimental isotropic dipole polarizability (in Angstrom^3) from ``exp2x.asp``."""
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if "alpha" in headers and "reference" in headers:
            cols = [clean_text(c) for c in rows[1].find_all(["td", "th"])]
            val = parse_float_or_none(cols[0]) if cols else None
            unc = parse_float_or_none(cols[1]) if len(cols) > 1 else None
            ref = cols[2] if len(cols) > 2 and cols[2] else None
            if val is not None:
                return PolarizabilityTensor(
                    isotropic=val,
                    uncertainty=unc,
                    unit=PolarizabilityUnit.ANGSTROM3.value,
                    reference=ref,
                )
    return None


def parse_calculated_polarizability_table(
    soup: BeautifulSoup,
    occurrence_index: int = 0,
) -> Tuple[Optional[PolarizabilityTensor], Optional[float]]:
    """Parse calculated 3x3 polarizability tensor (Angstrom^3) and ``<r^2>`` (Angstrom^2) from ``energy3x.asp``."""
    pol_tables = []
    r2_values: List[float] = []

    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) == 4:
            r0 = [clean_text(c).lower() for c in rows[0].find_all(["td", "th"])]
            r1 = [clean_text(c).lower() for c in rows[1].find_all(["td", "th"])]
            if r0 in (["x", "y", "z"], ["", "x", "y", "z"]) and r1 and r1[0] == "x":
                prev_text = clean_text(
                    table.find_previous(string=re.compile(r"Primitive|Traceless|Polarizabilit|Angstrom", re.I))
                ).lower()
                if "polarizabilit" in prev_text or "angstrom" in prev_text or "change units" in prev_text:
                    pol_tables.append(table)
        elif len(rows) == 2:
            r0 = [clean_text(c).lower() for c in rows[0].find_all(["td", "th"])]
            if r0 and "<r" in r0[0]:
                val = parse_float_or_none(r0[1]) if len(r0) > 1 else None
                if val is not None:
                    r2_values.append(val)

    pol_tensor: Optional[PolarizabilityTensor] = None
    if occurrence_index < len(pol_tables):
        tbl = pol_tables[occurrence_index]
        rows = tbl.find_all("tr")
        r_x = [clean_text(c) for c in rows[1].find_all(["td", "th"])]
        r_y = [clean_text(c) for c in rows[2].find_all(["td", "th"])]
        r_z = [clean_text(c) for c in rows[3].find_all(["td", "th"])]
        if len(r_x) >= 4 and len(r_y) >= 4 and len(r_z) >= 4:
            xx = parse_float_or_none(r_x[1])
            xy = parse_float_or_none(r_x[2])
            xz = parse_float_or_none(r_x[3])
            yy = parse_float_or_none(r_y[2])
            yz = parse_float_or_none(r_y[3])
            zz = parse_float_or_none(r_z[3])
            if xx is not None and yy is not None and zz is not None:
                iso = round((xx + yy + zz) / 3.0, 6)
                pol_tensor = PolarizabilityTensor(
                    xx=xx,
                    yy=yy,
                    zz=zz,
                    xy=xy,
                    xz=xz,
                    yz=yz,
                    isotropic=iso,
                    unit=PolarizabilityUnit.ANGSTROM3.value,
                )

    r2_val = r2_values[occurrence_index] if occurrence_index < len(r2_values) else None
    return pol_tensor, r2_val
