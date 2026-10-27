"""Parsers for experimental and calculated vibrational frequency spectra and ZPE."""

from __future__ import annotations

import re
from typing import List, Optional

from bs4 import BeautifulSoup

from cccbdb.models.property import VibrationalMode, VibrationalSpectrum
from cccbdb.models.provenance import DataProvenance
from cccbdb.parsers._utils import clean_text, parse_float_or_none, parse_frequency_value


def parse_experimental_vibrations(
    soup: BeautifulSoup,
    species_id: str,
    formula: str,
    point_group: Optional[str] = None,
    provenance: Optional[DataProvenance] = None,
) -> Optional[VibrationalSpectrum]:
    """Parse experimental vibrational frequencies and ZPE from ``exp2x.asp``."""
    modes: List[VibrationalMode] = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if len(rows) < 2:
            continue
        header_text = " ".join(clean_text(r).lower() for r in rows[:2])
        if "mode number" in header_text and "symmetry" in header_text and "fundamental" in header_text:
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) < 3 or not cols[0].isdigit():
                    continue
                mode_num = int(cols[0])
                sym = cols[1]
                fund_val, is_imag = parse_frequency_value(cols[2])
                harm_val, _ = parse_frequency_value(cols[3]) if len(cols) > 3 else (None, False)
                ref = cols[4] if len(cols) > 4 and cols[4] else None
                intensity = parse_float_or_none(cols[5]) if len(cols) > 5 else None
                desc = cols[-1] if len(cols) > 8 and cols[-1] else None

                freq_primary = harm_val if harm_val is not None else fund_val
                if freq_primary is not None:
                    modes.append(
                        VibrationalMode(
                            mode_number=mode_num,
                            symmetry=sym,
                            frequency_cm_inv=freq_primary,
                            fundamental_frequency_cm_inv=fund_val,
                            ir_intensity_km_mol=intensity,
                            is_imaginary=is_imag,
                            reference=ref,
                            description=desc,
                        )
                    )

    full_text = soup.get_text(" ", strip=True)
    zpe_match = re.search(
        r"vibrational zero-point energy:\s*([0-9]+\.?[0-9]*)\s*cm-1", full_text, re.IGNORECASE
    )
    zpe_val = float(zpe_match.group(1)) if zpe_match else None

    if not modes and zpe_val is None:
        return None

    return VibrationalSpectrum(
        species_id=species_id,
        formula=formula,
        method="Experimental",
        basis_set=None,
        point_group=point_group,
        conformer_index=1,
        modes=modes,
        zpe_unscaled_cm_inv=zpe_val,
        zpe_scaled_cm_inv=zpe_val,
        scaling_factor=1.0,
        provenance=provenance,
    )


def parse_calculated_vibrations_table(
    soup: BeautifulSoup,
    species_id: str,
    formula: str,
    method: str,
    basis_set: Optional[str],
    point_group: Optional[str] = None,
    conformer_index: int = 1,
    occurrence_index: int = 0,
    provenance: Optional[DataProvenance] = None,
) -> Optional[VibrationalSpectrum]:
    """Parse the ``occurrence_index``-th calculated vibrational frequency table from ``energy3x.asp`` or ``vibs2x.asp``."""
    vib_tables = []
    for table in soup.find_all("table"):
        if table.find("table"):
            continue
        rows = table.find_all("tr")
        if not rows:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if any("mode" in h for h in headers) and any("frequency" in h for h in headers):
            if "fundamental" not in " ".join(headers):
                vib_tables.append(table)

    modes: List[VibrationalMode] = []
    if occurrence_index < len(vib_tables):
        vtable = vib_tables[occurrence_index]
        for r in vtable.find_all("tr")[1:]:
            cols = [clean_text(c) for c in r.find_all(["td", "th"])]
            if len(cols) < 3 or not cols[0].isdigit():
                continue
            mode_num = int(cols[0])
            sym = cols[1]
            freq_val, is_imag = parse_frequency_value(cols[2])
            scaled_val, _ = parse_frequency_value(cols[3]) if len(cols) > 3 else (None, False)
            ir_val = parse_float_or_none(cols[4]) if len(cols) > 4 else None
            raman_val = parse_float_or_none(cols[5]) if len(cols) > 5 else None
            dep_p = parse_float_or_none(cols[6]) if len(cols) > 6 else None
            dep_u = parse_float_or_none(cols[7]) if len(cols) > 7 else None
            if freq_val is not None:
                modes.append(
                    VibrationalMode(
                        mode_number=mode_num,
                        symmetry=sym,
                        frequency_cm_inv=freq_val,
                        scaled_frequency_cm_inv=scaled_val,
                        ir_intensity_km_mol=ir_val,
                        raman_activity=raman_val,
                        depolarization_polar=dep_p,
                        depolarization_unpolar=dep_u,
                        is_imaginary=is_imag,
                    )
                )

    full_text = soup.get_text("\n", strip=True)
    unscaled_matches = re.findall(
        r"Unscaled Zero Point Vibrational Energy\s*\(zpe\)\s*([0-9]+\.?[0-9]*)\s*cm-1",
        full_text,
        re.IGNORECASE,
    )
    scaled_matches = re.findall(
        r"Scaled\s*\(by\s*([0-9]+\.?[0-9]*)\)\s*Zero Point Vibrational Energy\s*\(zpe\)\s*([0-9]+\.?[0-9]*)\s*cm-1",
        full_text,
        re.IGNORECASE,
    )

    zpe_unscaled = (
        float(unscaled_matches[occurrence_index])
        if occurrence_index < len(unscaled_matches)
        else None
    )
    scale_factor: Optional[float] = None
    zpe_scaled: Optional[float] = None
    if occurrence_index < len(scaled_matches):
        scale_factor = float(scaled_matches[occurrence_index][0])
        zpe_scaled = float(scaled_matches[occurrence_index][1])

    if not modes and zpe_unscaled is None:
        return None

    return VibrationalSpectrum(
        species_id=species_id,
        formula=formula,
        method=method,
        basis_set=basis_set,
        point_group=point_group,
        conformer_index=conformer_index,
        modes=modes,
        zpe_unscaled_cm_inv=zpe_unscaled,
        zpe_scaled_cm_inv=zpe_scaled,
        scaling_factor=scale_factor,
        provenance=provenance,
    )
