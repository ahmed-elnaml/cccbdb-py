"""Parsers for NIST CCCBDB electronic and thermochemical energy pages (energy2x.asp and energy3x.asp)."""

from __future__ import annotations

import re
from typing import Dict, List, Optional
from urllib.parse import parse_qs, urljoin, urlparse

from bs4 import BeautifulSoup, Tag

from cccbdb.chemistry.basis_sets import get_basis_set_by_nist_id, resolve_basis_set
from cccbdb.chemistry.methods import get_method_by_nist_id, resolve_method
from cccbdb.chemistry.symmetry import normalize_point_group_symbol
from cccbdb.chemistry.units import EnergyUnit
from cccbdb.exceptions import StructureChangedError
from cccbdb.models.calculation import CalculationResult, ConformerCalculation
from cccbdb.models.property import EnergyCollection, EnergyRecord
from cccbdb.models.provenance import DataProvenance
from cccbdb.parsers._utils import (
    clean_text,
    extract_casno_from_url,
    extract_species_header,
    make_soup,
    parse_float_or_none,
)
from cccbdb.parsers.frequencies import parse_calculated_vibrations_table
from cccbdb.parsers.geometries import parse_calculated_geometry_section
from cccbdb.parsers.properties import (
    parse_atomic_charges_table,
    parse_calculated_dipole_table,
    parse_calculated_polarizability_table,
    parse_calculated_quadrupole_tables,
    parse_calculated_rotational_constants,
)


def parse_energies_matrix_page(
    html: str,
    species_id: str,
    formula: Optional[str] = None,
    source_url: str = "https://cccbdb.nist.gov/energy2x.asp",
) -> EnergyCollection:
    """Parse a NIST CCCBDB ``energy2x.asp`` 2D method-by-basis-set matrix page into an :class:`EnergyCollection`."""
    soup = make_soup(html)
    hdr_formula, _ = extract_species_header(soup)
    mol_formula = formula or hdr_formula or "Unknown"
    prov = DataProvenance.from_html(
        source_url=source_url,
        html=html,
        original_unit=EnergyUnit.HARTREE.value,
    )

    records: List[EnergyRecord] = []

    for table in soup.find_all("table"):
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue

        # Detect if this table has basis set headers across columns (e.g., STO-3G, 3-21G, 6-31G*, cc-pVDZ...)
        basis_headers: Dict[int, str] = {}
        for row in rows:
            cells = row.find_all(["th", "td"])
            texts = [clean_text(c) for c in cells]
            # Check if at least 2 cells match known basis sets
            matched_cols: Dict[int, str] = {}
            for idx, txt in enumerate(texts):
                if not txt:
                    continue
                try:
                    b_spec = resolve_basis_set(txt)
                    if b_spec.nist_id != 0:
                        matched_cols[idx] = b_spec.name
                except ValueError:
                    continue

            if len(matched_cols) >= 2:
                basis_headers = matched_cols
                continue

            # Parse data rows
            # Case 1: Standard method x basis-set grid row
            if basis_headers:
                # Find the method cell (usually column 0 or column 1 if column 0 is category rowspan)
                method_name: Optional[str] = None
                method_col_idx = -1
                for idx in range(min(3, len(texts))):
                    cand = texts[idx]
                    if not cand:
                        continue
                    try:
                        m_spec = resolve_method(cand)
                        method_name = m_spec.name
                        method_col_idx = idx
                        break
                    except ValueError:
                        continue

                if method_name is not None:
                    m_spec = resolve_method(method_name)
                    # Align cells with basis_headers (handling optional category column offset)
                    header_min_col = min(basis_headers.keys())
                    offset = (method_col_idx + 1) - header_min_col
                    for h_col, b_name in basis_headers.items():
                        c_idx = h_col + offset
                        if 0 <= c_idx < len(cells):
                            cell = cells[c_idx]
                            val = parse_float_or_none(clean_text(cell))
                            if val is None:
                                # Also check if anchor text inside cell has the float
                                a_tag = cell.find("a")
                                if a_tag:
                                    val = parse_float_or_none(clean_text(a_tag))
                            if val is not None:
                                a_tag = cell.find("a", href=True)
                                detail_url = (
                                    urljoin("https://cccbdb.nist.gov/", a_tag["href"])
                                    if a_tag
                                    else None
                                )
                                records.append(
                                    EnergyRecord(
                                        species_id=species_id,
                                        formula=mol_formula,
                                        method=m_spec.name,
                                        basis_set=b_name,
                                        energy_type="electronic_0K",
                                        value=val,
                                        unit=EnergyUnit.HARTREE.value,
                                        method_category=m_spec.category.value,
                                        core_correlation=m_spec.core_correlation,
                                        detail_url=detail_url,
                                        provenance=prov,
                                    )
                                )
                continue

            # Case 2: Predefined basis set table (e.g., AM1, PM3, G1, G2, G3, G4, CBS-Q)
            for idx, cell in enumerate(cells):
                txt = clean_text(cell)
                if not txt:
                    continue
                try:
                    m_spec = resolve_method(txt)
                except ValueError:
                    continue
                if m_spec.predefined_basis and idx + 1 < len(cells):
                    val_cell = cells[idx + 1]
                    val = parse_float_or_none(clean_text(val_cell))
                    if val is not None:
                        a_tag = val_cell.find("a", href=True)
                        detail_url = (
                            urljoin("https://cccbdb.nist.gov/", a_tag["href"])
                            if a_tag
                            else None
                        )
                        records.append(
                            EnergyRecord(
                                species_id=species_id,
                                formula=mol_formula,
                                method=m_spec.name,
                                basis_set=None,
                                energy_type="electronic_0K",
                                value=val,
                                unit=EnergyUnit.HARTREE.value,
                                method_category=m_spec.category.value,
                                core_correlation=m_spec.core_correlation,
                                detail_url=detail_url,
                                provenance=prov,
                            )
                        )

    if not records:
        raise StructureChangedError(
            f"No electronic energy records could be extracted from '{source_url}'."
        )

    return EnergyCollection(records=records)


def parse_calculation_result_page(
    html: str,
    source_url: str = "https://cccbdb.nist.gov/energy3x.asp",
    species_id: Optional[str] = None,
) -> CalculationResult:
    """Parse a NIST CCCBDB single-calculation dossier page (``energy3x.asp``) into a :class:`CalculationResult`.

    Handles both single-conformer pages (e.g., H2O, CO, N(CH3)3) and multi-conformer/multi-state
    pages (e.g., CH2 triplet C2v and Dinfh conformers).
    """
    soup = make_soup(html)
    full_text = soup.get_text("\n", strip=True)

    formula, name = extract_species_header(soup)
    casno = species_id or extract_casno_from_url(source_url) or "00000"

    # Extract "using model chemistry: METHOD/BASIS" or "using model chemistry: AM1"
    model_match = re.search(r"using model chemistry:\s*([^\n\r]+)", full_text, re.IGNORECASE)
    if not model_match:
        raise StructureChangedError(
            f"Could not locate 'using model chemistry:' banner on '{source_url}'."
        )
    model_str = model_match.group(1).strip()
    if "/" in model_str:
        raw_method, raw_basis = model_str.split("/", 1)
        method_name = raw_method.strip()
        basis_name: Optional[str] = raw_basis.strip()
    else:
        method_name = model_str.strip()
        basis_name = None

    try:
        m_spec = resolve_method(method_name)
        method_name = m_spec.name
        m_cat = m_spec.category.value
        core_corr = m_spec.core_correlation
        m_id = m_spec.nist_id
    except ValueError:
        m_cat = None
        core_corr = None
        m_id = None

    b_id: Optional[int] = None
    if basis_name:
        try:
            b_spec = resolve_basis_set(basis_name)
            basis_name = b_spec.name
            b_id = b_spec.nist_id
        except ValueError:
            pass

    prov = DataProvenance.from_html(
        source_url=source_url,
        html=html,
        original_unit=EnergyUnit.HARTREE.value,
    )

    # Parse "States and conformations" summary table if present
    conformer_metas: List[Dict[str, object]] = []
    for table in soup.find_all("table"):
        rows = table.find_all("tr")
        if not rows:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if "state" in headers and "conformation" in headers and "conformer description" in headers:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) >= 5 and cols[0].isdigit() and cols[1].isdigit():
                    pg_raw = cols[3]
                    try:
                        pg_norm = normalize_point_group_symbol(pg_raw)
                    except ValueError:
                        pg_norm = pg_raw
                    conformer_metas.append(
                        {
                            "state_index": int(cols[0]),
                            "conformer_index": int(cols[1]),
                            "is_minimum": cols[2].lower() in ("yes", "true", "1"),
                            "point_group": pg_norm,
                            "electronic_state": cols[4],
                        }
                    )
            break

    if not conformer_metas:
        conformer_metas = [
            {
                "state_index": 1,
                "conformer_index": 1,
                "is_minimum": True,
                "point_group": None,
                "electronic_state": None,
            }
        ]

    # Extract energy tables (each conformer has an energy table with "Energy at 0K", "Energy at 298.15K", "HF Energy", "Nuclear repulsion energy")
    energy_tables: List[Tag] = []
    for table in soup.find_all("table"):
        txt = table.get_text(" ", strip=True)
        if "Energy at 0K" in txt or "Nuclear repulsion energy" in txt:
            # Ensure it's the innermost table, not a layout wrapper table
            if not table.find("table"):
                energy_tables.append(table)

    energy_key_map = {
        "energy at 0k": "electronic_0K",
        "energy at 298.15k": "energy_298K",
        "hf energy": "hf_energy",
        "nuclear repulsion energy": "nuclear_repulsion",
    }

    conformers: List[ConformerCalculation] = []
    for idx, meta in enumerate(conformer_metas):
        c_idx = int(meta["conformer_index"])  # type: ignore[arg-type]
        s_idx = int(meta["state_index"])  # type: ignore[arg-type]
        is_min = bool(meta["is_minimum"])
        pg = str(meta["point_group"]) if meta.get("point_group") else None
        el_state = str(meta["electronic_state"]) if meta.get("electronic_state") else None

        energies_dict: Dict[str, EnergyRecord] = {}
        if idx < len(energy_tables):
            etable = energy_tables[idx]
            for r in etable.find_all("tr"):
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) >= 2:
                    label = cols[0].strip()
                    val = parse_float_or_none(cols[1])
                    if val is not None and label.lower() in energy_key_map:
                        e_type = energy_key_map[label.lower()]
                        rec = EnergyRecord(
                            species_id=casno,
                            formula=formula or "Unknown",
                            method=method_name,
                            basis_set=basis_name,
                            energy_type=e_type,
                            value=val,
                            unit=EnergyUnit.HARTREE.value,
                            method_category=m_cat,
                            core_correlation=core_corr,
                            conformer_index=c_idx,
                            point_group=pg,
                            electronic_state=el_state,
                            detail_url=source_url,
                            provenance=prov,
                        )
                        energies_dict[e_type] = rec
                        energies_dict[label] = rec

        geom = parse_calculated_geometry_section(
            soup=soup,
            species_id=casno,
            formula=formula or "Unknown",
            name=name,
            method=method_name,
            basis_set=basis_name,
            point_group=pg,
            electronic_state=el_state,
            conformer_index=c_idx,
            is_minimum=is_min,
            occurrence_index=idx,
            provenance=prov,
        )
        if geom and geom.point_group and not pg:
            pg = geom.point_group

        vibs = parse_calculated_vibrations_table(
            soup=soup,
            species_id=casno,
            formula=formula or "Unknown",
            method=method_name,
            basis_set=basis_name,
            point_group=pg,
            conformer_index=c_idx,
            occurrence_index=idx,
            provenance=prov,
        )
        rot = parse_calculated_rotational_constants(soup, occurrence_index=idx)
        charges = parse_atomic_charges_table(soup, occurrence_index=idx)
        dipole = parse_calculated_dipole_table(
            soup=soup,
            species_id=casno,
            formula=formula or "Unknown",
            method=method_name,
            basis_set=basis_name,
            point_group=pg,
            electronic_state=el_state,
            occurrence_index=idx,
            provenance=prov,
        )
        quads = parse_calculated_quadrupole_tables(soup, occurrence_index=idx)
        polar, r2_val = parse_calculated_polarizability_table(soup, occurrence_index=idx)

        conformers.append(
            ConformerCalculation(
                conformer_index=c_idx,
                state_index=s_idx,
                is_minimum=is_min,
                point_group=pg,
                electronic_state=el_state,
                energies=energies_dict,
                geometry=geom,
                vibrations=vibs,
                rotational_constants=rot,
                atomic_charges=charges,
                dipole=dipole,
                quadrupoles=quads,
                polarizability=polar,
                r2_angstrom2=r2_val,
            )
        )

    return CalculationResult(
        species_id=casno,
        formula=formula or "Unknown",
        name=name,
        method=method_name,
        basis_set=basis_name,
        method_id=m_id,
        basis_id=b_id,
        conformers=conformers,
        provenance=prov,
    )
