"""Species discovery, chemical identifier, and experimental reference parser."""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from cccbdb.chemistry.formula import (
    calculate_molecular_weight,
    format_cas_hyphenated,
    normalize_cas,
    parse_formula,
    to_hill_formula,
)
from cccbdb.chemistry.symmetry import normalize_point_group_symbol
from cccbdb.exceptions import StructureChangedError
from cccbdb.models.molecule import ExperimentalDataSummary, SpeciesSummary
from cccbdb.models.provenance import DataProvenance, ReferenceCitation
from cccbdb.parsers._utils import clean_text, extract_casno_from_url, extract_species_header, make_soup
from cccbdb.parsers.frequencies import parse_experimental_vibrations
from cccbdb.parsers.geometries import parse_experimental_geometry
from cccbdb.parsers.properties import (
    parse_experimental_dipole,
    parse_experimental_polarizability,
    parse_experimental_quadrupoles,
    parse_experimental_rotational_constants,
    parse_experimental_thermochemistry,
)


# Curated built-in benchmark catalog of common NIST CCCBDB molecules
# Enables instant zero-network resolution for standard benchmark species and isomer searches
BUILTIN_SPECIES_CATALOG: List[SpeciesSummary] = [
    SpeciesSummary(
        id="7732185",
        cas_rn="7732-18-5",
        name="Water",
        formula="H2O",
        hill_formula="H2O",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C2v",
        inchi="InChI=1S/H2O/h1H2",
        inchikey="XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        smiles="O",
        molecular_weight=calculate_molecular_weight("H2O"),
    ),
    SpeciesSummary(
        id="7789200",
        cas_rn="7789-20-0",
        name="Water-d2",
        formula="D2O",
        hill_formula="D2O",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C2v",
        inchi="InChI=1S/H2O/h1H2/i1+1D2",
        inchikey="XLYOFNOQVPJJNP-ZSJDYOACSA-N",
        smiles="[2H]O[2H]",
        molecular_weight=calculate_molecular_weight("D2O"),
    ),
    SpeciesSummary(
        id="630080",
        cas_rn="630-08-0",
        name="Carbon monoxide",
        formula="CO",
        hill_formula="CO",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σ+",
        point_group="Cinfv",
        inchi="InChI=1S/CO/c1-2",
        inchikey="UGFAIRIUMAVXCW-UHFFFAOYSA-N",
        smiles="[C-]#[O+]",
        molecular_weight=calculate_molecular_weight("CO"),
    ),
    SpeciesSummary(
        id="124389",
        cas_rn="124-38-9",
        name="Carbon dioxide",
        formula="CO2",
        hill_formula="CO2",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σg+",
        point_group="Dinfh",
        inchi="InChI=1S/CO2/c2-1-3",
        inchikey="CURLTUGMZLYLDI-UHFFFAOYSA-N",
        smiles="O=C=O",
        molecular_weight=calculate_molecular_weight("CO2"),
    ),
    SpeciesSummary(
        id="1333740",
        cas_rn="1333-74-0",
        name="Hydrogen",
        formula="H2",
        hill_formula="H2",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σg+",
        point_group="Dinfh",
        inchi="InChI=1S/H2/h1H",
        inchikey="UFHFLCQGNIYNRP-UHFFFAOYSA-N",
        smiles="[H][H]",
        molecular_weight=calculate_molecular_weight("H2"),
    ),
    SpeciesSummary(
        id="7727379",
        cas_rn="7727-37-9",
        name="Nitrogen",
        formula="N2",
        hill_formula="N2",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σg+",
        point_group="Dinfh",
        inchi="InChI=1S/N2/c1-2",
        inchikey="IJGRMHOSHXDMSA-UHFFFAOYSA-N",
        smiles="N#N",
        molecular_weight=calculate_molecular_weight("N2"),
    ),
    SpeciesSummary(
        id="7782447",
        cas_rn="7782-44-7",
        name="Oxygen",
        formula="O2",
        hill_formula="O2",
        charge=0,
        spin_multiplicity=3,
        electronic_state="3Σg-",
        point_group="Dinfh",
        inchi="InChI=1S/O2/c1-2",
        inchikey="MYMOFIZGZYHOMD-UHFFFAOYSA-N",
        smiles="O=O",
        molecular_weight=calculate_molecular_weight("O2"),
    ),
    SpeciesSummary(
        id="7664417",
        cas_rn="7664-41-7",
        name="Ammonia",
        formula="NH3",
        hill_formula="H3N",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C3v",
        inchi="InChI=1S/H3N/h1H3",
        inchikey="QGZKDVFQNNGYKY-UHFFFAOYSA-N",
        smiles="N",
        molecular_weight=calculate_molecular_weight("NH3"),
    ),
    SpeciesSummary(
        id="74828",
        cas_rn="74-82-8",
        name="Methane",
        formula="CH4",
        hill_formula="CH4",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="Td",
        inchi="InChI=1S/CH4/h1H4",
        inchikey="VNWKTOKETHGBQD-UHFFFAOYSA-N",
        smiles="C",
        molecular_weight=calculate_molecular_weight("CH4"),
    ),
    SpeciesSummary(
        id="676551",
        cas_rn="676-55-1",
        name="Methane-d2",
        formula="CH2D2",
        hill_formula="CH2D2",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C2v",
        inchi="InChI=1S/CH4/h1H4/i1D2",
        inchikey="VNWKTOKETHGBQD-PICCSMPSSA-N",
        smiles="[2H]C([2H])",
        molecular_weight=calculate_molecular_weight("CH2D2"),
    ),
    SpeciesSummary(
        id="2465567",
        cas_rn="2465-56-7",
        name="Methylene",
        formula="CH2",
        hill_formula="CH2",
        charge=0,
        spin_multiplicity=3,
        electronic_state="3B1",
        point_group="C2v",
        inchi="InChI=1S/CH2/h1H2",
        inchikey="HZVOZRGWRWCICA-UHFFFAOYSA-N",
        smiles="[CH2]",
        molecular_weight=calculate_molecular_weight("CH2"),
    ),
    SpeciesSummary(
        id="50000",
        cas_rn="50-00-0",
        name="Formaldehyde",
        formula="H2CO",
        hill_formula="CH2O",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C2v",
        inchi="InChI=1S/CH2O/c1-2/h1H2",
        inchikey="WSFSSNUMVMOOMR-UHFFFAOYSA-N",
        smiles="C=O",
        molecular_weight=calculate_molecular_weight("H2CO"),
    ),
    SpeciesSummary(
        id="74862",
        cas_rn="74-86-2",
        name="Acetylene",
        formula="C2H2",
        hill_formula="C2H2",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σg+",
        point_group="Dinfh",
        inchi="InChI=1S/C2H2/c1-2/h1-2H",
        inchikey="HSFWRNGVRCDJHI-UHFFFAOYSA-N",
        smiles="C#C",
        molecular_weight=calculate_molecular_weight("C2H2"),
    ),
    SpeciesSummary(
        id="74851",
        cas_rn="74-85-1",
        name="Ethylene",
        formula="C2H4",
        hill_formula="C2H4",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Ag",
        point_group="D2h",
        inchi="InChI=1S/C2H4/c1-2/h1-2H2",
        inchikey="VGGSQFUCUMXWEO-UHFFFAOYSA-N",
        smiles="C=C",
        molecular_weight=calculate_molecular_weight("C2H4"),
    ),
    SpeciesSummary(
        id="74840",
        cas_rn="74-84-0",
        name="Ethane",
        formula="C2H6",
        hill_formula="C2H6",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Ag",
        point_group="D3d",
        inchi="InChI=1S/C2H6/c1-2/h1-2H3",
        inchikey="OTMSDBZUPAUEDD-UHFFFAOYSA-N",
        smiles="CC",
        molecular_weight=calculate_molecular_weight("C2H6"),
    ),
    # C2H6O Isomers (Ethanol vs. Dimethyl ether) for explicit isomer disambiguation
    SpeciesSummary(
        id="64175",
        cas_rn="64-17-5",
        name="Ethanol",
        formula="CH3CH2OH",
        hill_formula="C2H6O",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A'",
        point_group="Cs",
        inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        smiles="CCO",
        molecular_weight=calculate_molecular_weight("C2H6O"),
    ),
    SpeciesSummary(
        id="115106",
        cas_rn="115-10-6",
        name="Dimethyl ether",
        formula="CH3OCH3",
        hill_formula="C2H6O",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C2v",
        inchi="InChI=1S/C2H6O/c1-3-2/h1-2H3",
        inchikey="LCGLNKUTAGEVQW-UHFFFAOYSA-N",
        smiles="COC",
        molecular_weight=calculate_molecular_weight("C2H6O"),
    ),
    # C6H6 Isomers (Benzene vs. Dewar benzene)
    SpeciesSummary(
        id="71432",
        cas_rn="71-43-2",
        name="Benzene",
        formula="C6H6",
        hill_formula="C6H6",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1g",
        point_group="D6h",
        inchi="InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
        inchikey="UHOVQNZJYSORNB-UHFFFAOYSA-N",
        smiles="c1ccccc1",
        molecular_weight=calculate_molecular_weight("C6H6"),
    ),
    SpeciesSummary(
        id="5649956",
        cas_rn="5649-95-6",
        name="Dewar benzene",
        formula="C6H6",
        hill_formula="C6H6",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C2v",
        inchi="InChI=1S/C6H6/c1-2-5-3-4-6(1)5/h1-6H",
        inchikey="QDWQMPDJWGRJAA-UHFFFAOYSA-N",
        smiles="C1=CC2C1C=C2",
        molecular_weight=calculate_molecular_weight("C6H6"),
    ),
    SpeciesSummary(
        id="75503",
        cas_rn="75-50-3",
        name="Trimethylamine",
        formula="N(CH3)3",
        hill_formula="C3H9N",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1A1",
        point_group="C3v",
        inchi="InChI=1S/C3H9N/c1-4(2)3/h1-3H3",
        inchikey="GETQZCLCWQTVFV-UHFFFAOYSA-N",
        smiles="CN(C)C",
        molecular_weight=calculate_molecular_weight("C3H9N"),
    ),
    SpeciesSummary(
        id="7447418",
        cas_rn="7447-41-8",
        name="Lithium chloride",
        formula="LiCl",
        hill_formula="ClLi",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σ+",
        point_group="Cinfv",
        inchi="InChI=1S/ClH.Li/h1H;/q;+1/p-1",
        inchikey="GWKQBVSWPJMOQI-UHFFFAOYSA-M",
        smiles="[Li]Cl",
        molecular_weight=calculate_molecular_weight("LiCl"),
    ),
    SpeciesSummary(
        id="7647010",
        cas_rn="7647-01-0",
        name="Hydrogen chloride",
        formula="HCl",
        hill_formula="ClH",
        charge=0,
        spin_multiplicity=1,
        electronic_state="1Σ+",
        point_group="Cinfv",
        inchi="InChI=1S/ClH/h1H",
        inchikey="VEXZGXHMUGYJMC-UHFFFAOYSA-N",
        smiles="Cl",
        molecular_weight=calculate_molecular_weight("HCl"),
    ),
]


def parse_experimental_page(html: str, source_url: str = "https://cccbdb.nist.gov/exp2x.asp") -> ExperimentalDataSummary:
    """Parse a NIST CCCBDB ``exp2x.asp`` HTML page into a structured :class:`ExperimentalDataSummary`."""
    soup = make_soup(html)
    formula, name = extract_species_header(soup)
    if not formula:
        raise StructureChangedError(
            f"Could not locate 'Experimental data for <Formula> (<Name>)' header in {source_url}"
        )

    casno = extract_casno_from_url(source_url) or "00000"
    # Also check links on page for casno
    if casno == "00000":
        for a in soup.find_all("a", href=True):
            c = extract_casno_from_url(a["href"])
            if c:
                casno = c
                break

    cas_rn = format_cas_hyphenated(casno) if len(casno) >= 5 else casno
    prov = DataProvenance.from_html(source_url=source_url, html=html)

    # Extract synonyms, InChI, InChIKey, SMILES, State, Conformation
    synonyms: List[str] = []
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    smiles: Optional[str] = None
    iupac_name: Optional[str] = None
    electronic_state: Optional[str] = None
    point_group: Optional[str] = None

    for table in soup.find_all("table"):
        rows = table.find_all("tr")
        for r_idx, row in enumerate(rows):
            cells = [clean_text(c) for c in row.find_all(["td", "th"])]
            if not cells:
                continue
            if cells[0].lower().startswith("other names") and len(cells) >= 2:
                synonyms = [s.strip() for s in cells[1].split(";") if s.strip()]
            # Check header row for INChI / SMILES
            lower_cells = [c.lower() for c in cells]
            if "inchi" in lower_cells and r_idx + 1 < len(rows):
                val_cells = [clean_text(c) for c in rows[r_idx + 1].find_all(["td", "th"])]
                for col_i, h_name in enumerate(lower_cells):
                    if col_i >= len(val_cells):
                        break
                    val = val_cells[col_i]
                    if h_name == "inchi":
                        inchi = val
                    elif h_name == "inchikey":
                        inchikey = val
                    elif h_name == "smiles":
                        smiles = val
                    elif "iupac" in h_name:
                        iupac_name = val
            if "state" in lower_cells and "conformation" in lower_cells and r_idx + 1 < len(rows):
                val_cells = [clean_text(c) for c in rows[r_idx + 1].find_all(["td", "th"])]
                if len(val_cells) >= 2:
                    electronic_state = val_cells[0] or None
                    try:
                        point_group = normalize_point_group_symbol(val_cells[1])
                    except ValueError:
                        point_group = val_cells[1] or None

    # Also check explicit "Point Group C∞v" text on page
    if not point_group:
        pg_match = re.search(r"Point Group(?:\s+is)?\s+([A-Za-z0-9\*∞]+)", soup.get_text(" ", strip=True))
        if pg_match:
            try:
                point_group = normalize_point_group_symbol(pg_match.group(1))
            except ValueError:
                point_group = pg_match.group(1)

    thermo = parse_experimental_thermochemistry(soup, prov)
    geom = parse_experimental_geometry(
        soup,
        species_id=casno,
        formula=formula,
        name=name,
        point_group=point_group,
        electronic_state=electronic_state,
        provenance=prov,
    )
    vibs = parse_experimental_vibrations(
        soup,
        species_id=casno,
        formula=formula,
        point_group=point_group,
        provenance=prov,
    )
    rot = parse_experimental_rotational_constants(soup)
    dipole = parse_experimental_dipole(
        soup,
        species_id=casno,
        formula=formula,
        point_group=point_group,
        provenance=prov,
    )
    quads = parse_experimental_quadrupoles(soup)
    polar = parse_experimental_polarizability(soup)

    # Extract bibliographic references table ("squib", "reference", "DOI")
    references: Dict[str, ReferenceCitation] = {}
    for table in soup.find_all("table"):
        rows = table.find_all("tr")
        if not rows:
            continue
        headers = [clean_text(c).lower() for c in rows[0].find_all(["th", "td"])]
        if "squib" in headers and "reference" in headers:
            sq_idx = headers.index("squib")
            ref_idx = headers.index("reference")
            doi_idx = headers.index("doi") if "doi" in headers else None
            for r in rows[1:]:
                cols = [clean_text(c) for c in r.find_all(["td", "th"])]
                if len(cols) > max(sq_idx, ref_idx) and cols[sq_idx]:
                    squib = cols[sq_idx]
                    citation = cols[ref_idx] if len(cols) > ref_idx else None
                    doi = cols[doi_idx] if (doi_idx is not None and len(cols) > doi_idx and cols[doi_idx]) else None
                    references[squib] = ReferenceCitation(squib=squib, citation=citation, doi=doi)

    _, charge = parse_formula(formula)
    return ExperimentalDataSummary(
        species_id=casno,
        cas_rn=cas_rn,
        name=name or formula,
        formula=formula,
        charge=charge,
        synonyms=synonyms,
        inchi=inchi,
        inchikey=inchikey,
        smiles=smiles,
        iupac_name=iupac_name,
        electronic_state=electronic_state,
        point_group=point_group,
        thermochemistry=thermo,
        geometry=geom,
        vibrations=vibs,
        rotational_constants=rot,
        dipole=dipole,
        quadrupoles=quads,
        polarizability=polar,
        references=references,
        provenance=prov,
    )


def parse_species_disambiguation_page(html: str, source_url: str = "https://cccbdb.nist.gov/getformx.asp") -> List[SpeciesSummary]:
    """Parse a NIST CCCBDB species disambiguation / index HTML table into :class:`SpeciesSummary` records."""
    soup = make_soup(html)
    results: List[SpeciesSummary] = []
    seen_keys: set[tuple[str, int]] = set()

    for a in soup.find_all("a", href=True):
        href = a["href"]
        casno = extract_casno_from_url(href)
        if not casno:
            continue
        charge_match = re.search(r"charge=([+-]?\d+)", href)
        charge = int(charge_match.group(1)) if charge_match else 0
        key = (casno, charge)
        if key in seen_keys:
            continue

        row = a.find_parent("tr")
        formula_str = ""
        name_str = clean_text(a)
        pg_str: Optional[str] = None
        inchi_str: Optional[str] = None

        if row:
            cells = [clean_text(c) for c in row.find_all(["td", "th"])]
            for c_text in cells:
                if c_text.startswith("InChI="):
                    inchi_str = c_text
                elif re.match(r"^[CDTIOK][0-9\*∞sihvd]{0,4}$", c_text, re.IGNORECASE):
                    try:
                        pg_str = normalize_point_group_symbol(c_text)
                    except ValueError:
                        pass
            if len(cells) >= 2 and not formula_str:
                # Often [formula, name, casno, ...] or [name, formula, ...]
                for cand in cells[:2]:
                    try:
                        parse_formula(cand)
                        formula_str = cand
                        break
                    except ValueError:
                        continue

        if not formula_str:
            try:
                parse_formula(name_str)
                formula_str = name_str
            except ValueError:
                formula_str = "H2O"

        try:
            hill = to_hill_formula(formula_str)
            mw = calculate_molecular_weight(formula_str)
        except ValueError:
            hill = formula_str
            mw = 0.0

        seen_keys.add(key)
        results.append(
            SpeciesSummary(
                id=casno,
                cas_rn=format_cas_hyphenated(casno) if len(casno) >= 5 else casno,
                name=name_str or formula_str,
                formula=formula_str,
                hill_formula=hill,
                charge=charge,
                point_group=pg_str,
                inchi=inchi_str,
                molecular_weight=mw,
            )
        )
    return results
