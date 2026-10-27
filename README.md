# cccbdb-py

[![Tests](https://github.com/your-username/cccbdb-py/actions/workflows/tests.yml/badge.svg)](https://github.com/your-username/cccbdb-py/actions/workflows/tests.yml)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

**An open-source Python toolkit for accessing, parsing, validating, and analyzing computational chemistry data from the NIST CCCBDB.**

`cccbdb-py` is a community-maintained Python library designed to simplify programmatic access to the National Institute of Standards and Technology (NIST) [Computational Chemistry Comparison and Benchmark Database (CCCBDB, SRD 101)](https://cccbdb.nist.gov/). It provides a clean, typed, reproducible, and researcher-friendly interface for discovering molecular species, retrieving quantum chemical and experimental reference data, converting physical units, analyzing molecular point group symmetry, caching results locally, and exporting machine-learning-ready datasets.

> **Important Notice**: `cccbdb-py` is an independent, community-maintained client and is **not** an official NIST API, nor is it affiliated with, endorsed by, or maintained by NIST.

---

## Architecture & Data-Access Pipeline

```text
+-------------------------------------------------------------------+
|                        NIST CCCBDB (SRD 101)                      |
|      exp2x.asp | energy2x.asp | energy3x.asp | geom3x.asp         |
+-------------------------------------------------------------------+
                                  |
                                  v
+-------------------------------------------------------------------+
|             Transport Layer (cccbdb.transport.http)               |
|  - Polite RateLimiter (configurable delay >= 1.5s)                |
|  - ASP Session Cookie Manager & Auto-Recovery ("lost molecule")   |
+-------------------------------------------------------------------+
                                  |
                    +-------------+-------------+
                    |                           |
                    v                           v
+-----------------------------------+ +-----------------------------+
|   Storage Layer (cccbdb.storage)  | | Isolated HTML Parsers       |
| - SHA-256 Content-Addressable     | | (cccbdb.parsers.*)          |
|   DiskCache + Offline Mode        | | - Species & Identifiers     |
| - Normalized SQLite Database      | | - 2D Energy Matrices        |
+-----------------------------------+ | - Cartesians & Z-Matrices   |
                                      | - Vibrational Spectra & ZPE |
                                      | - Dipoles, Charges & Thermo |
                                      +-----------------------------+
                                                    |
                                                    v
+-------------------------------------------------------------------+
|         Computational Chemistry Ontology (cccbdb.chemistry)       |
| - Units: CODATA 2018 (Hartree, kJ/mol, kcal/mol, eV, cm^-1, Bohr) |
| - Symmetry: Schoenflies Point Groups, Rotational Symmetry (sigma) |
| - Methods (45+) & Basis Sets (45+) + Vibrational Scaling Factors  |
+-------------------------------------------------------------------+
                                  |
                                  v
+-------------------------------------------------------------------+
|               Public Typed API & CLI (cccbdb.CCCBDB)              |
|     pandas DataFrame | JSON | CSV + Manifest | XMol .xyz | ASE    |
+-------------------------------------------------------------------+
```

---

## Core Features

1. **Molecular Search & Isomer Disambiguation (`cccbdb.parsers.species`)**
   - Search by chemical formula (including Hill ordering and condensed structural formulas like `CH3(CH2)4CH3`), chemical name, CAS Registry Number (with check-digit verification), InChI, SMILES, Schoenflies point group, or constituent elements.
   - Explicit isomer disambiguation: querying `db.get_molecule("C2H6O")` raises `AmbiguousSpeciesError` listing `Ethanol` (`64-17-5`) and `Dimethyl ether` (`115-10-6`), preventing silent scientific misattribution.

2. **Quantum Chemical & Experimental Property Retrieval**
   - **Energies**: 0 K electronic energies, 298.15 K enthalpies/energies, Hartree-Fock reference energies, nuclear repulsion energies, and zero-point vibrational energies (ZPE) across 45+ methods (`HF`, `MP2`, `MP2=FULL`, `MP3`, `MP4`, `CCSD`, `CCSD(T)`, `QCISD(T)`, `B3LYP`, `M06-2X`, `wB97X-D`, `PBE1PBE`, `B2PLYP`, `G4`, `CBS-QB3`, etc.) and 45+ basis sets (`STO-3G`, `6-31G*`, `6-311+G(3df,2p)`, `cc-pVDZ`, `cc-pVTZ`, `cc-pVQZ`, `aug-cc-pVTZ`, `cc-pCVTZ`, `Def2TZVPP`, `SDD`, `LANL2DZ`, etc.).
   - **Geometries**: 3D Cartesian coordinates ($\text{\AA}$ or $\text{Bohr}$), internal coordinates ($r_e, r_0$, bond angles, dihedrals), automatic dummy-atom (`X`) filtering, interatomic distance matrices, and center-of-mass calculation.
   - **Vibrational Frequencies**: Normal-mode symmetries ($\Sigma, \Pi, A_1, B_2, E$), harmonic and fundamental wavenumbers ($\text{cm}^{-1}$), empirical Scott-Radom / NIST frequency scaling factors, IR intensities ($\text{km/mol}$), Raman activities ($\text{\AA}^4/\text{u}$), and transition-state imaginary mode (`i`) detection.
   - **Electrostatics, Rotation & Thermochemistry**: Mulliken/CHELPG/AIM/ESP partial atomic charges, dipole moment vectors ($\text{Debye}, \text{C}\cdot\text{m}, \text{a.u.}$), primitive and traceless quadrupole tensors ($\text{Debye}\cdot\text{\AA}$), anisotropic polarizability tensors ($\text{\AA}^3, \text{Bohr}^3$), rotational constants ($A, B, C$ in $\text{cm}^{-1}, \text{GHz}, \text{MHz}$), and experimental thermochemistry ($\Delta_f H^\circ_{298}, \Delta_f H^\circ_{0}, S^\circ_{298}, C_p$).

3. **Scientific Correctness & Data Provenance**
   - Missing or uncomputed entries (`-` or empty cells) are never coerced to `0.0`.
   - Every extracted datum embeds a `DataProvenance` object recording `source_url`, `retrieved_at` (ISO-8601 UTC), `parser_version`, `raw_content_sha256`, `original_unit`, and `reference_squib`.
   - `.xyz` structure export is strictly gated on complete 3D Cartesian coordinates (`IncompleteGeometryError` is raised if only internal bond lengths/angles exist).

---

## Installation

```bash
# Minimal core installation (httpx, beautifulsoup4, lxml, pydantic)
pip install cccbdb-py

# With scientific data analysis stack (pandas, numpy)
pip install "cccbdb-py[data]"

# With computational chemistry structure adapters (ASE, RDKit)
pip install "cccbdb-py[chemistry]"

# Development installation from source
git clone https://github.com/your-username/cccbdb-py.git
cd cccbdb-py
pip install -e ".[dev]"
```

---

## Quickstart: Python API

### Example A — Search for Molecules & Handle Isomers

```python
from cccbdb import CCCBDB

with CCCBDB() as db:
    # Search by chemical formula
    molecules = db.search(formula="C2H6O")
    for molecule in molecules:
        print(molecule.id, molecule.formula, molecule.name, molecule.point_group)
    # Output:
    # 64175 CH3CH2OH Ethanol Cs
    # 115106 CH3OCH3 Dimethyl ether C2v
```

### Example B — Retrieve Electronic Energies & Convert Units

```python
from cccbdb import CCCBDB, EnergyUnit

with CCCBDB() as db:
    water = db.get_molecule("H2O")
    energies = water.get_energies(
        method="B3LYP",
        basis_set="6-31G*"
    )
    print(energies.to_dataframe())

    # Convert from Hartree to kJ/mol or kcal/mol using CODATA 2018 constants
    energies_kj = energies.convert_units(EnergyUnit.KJ_MOL)
    print(energies_kj.to_dicts()[0])
```

### Example C — Retrieve Optimized Molecular Geometry & Export `.xyz`

```python
from cccbdb import CCCBDB

with CCCBDB() as db:
    water = db.get_molecule("H2O")
    geometry = water.get_geometry(
        method="B3LYP",
        basis_set="6-31G*"
    )
    print("Atoms:", geometry.atoms)
    print("Coordinates (Angstrom):", geometry.coordinates)
    print("Symmetry Point Group:", geometry.point_group)
    print("Rotational Symmetry Number (sigma):", geometry.symmetry_info.rotational_symmetry_number)
    geometry.to_xyz("water.xyz")
```

### Example D — Collect & Export a Reproducible Benchmark Dataset

```python
from cccbdb import CCCBDB

with CCCBDB(cache=True, use_sqlite=True) as db:
    molecules = db.search(elements=["C", "H", "O"])
    dataset = db.collect(
        molecules=molecules[:20],
        properties=["energy", "dipole", "zpe"],
        method="B3LYP",
        basis_set="6-31G*",
        delay=2.0,
    )
    # Writes molecular_data.csv and molecular_data.csv.manifest.json
    dataset.to_csv("molecular_data.csv")
```

---

## Command-Line Interface (`cccbdb`)

`cccbdb-py` installs a rich terminal CLI for shell scripts and rapid inspection:

```bash
# Search molecules by formula, point group, or constituent elements
cccbdb search --formula C6H6
cccbdb search --point-group D6h

# Inspect molecular symmetry properties and rotational symmetry numbers
cccbdb symmetry C3v

# List supported quantum chemical methods and Gaussian basis sets
cccbdb methods --category "Coupled Cluster"
cccbdb basis-sets --family Dunning

# Convert between quantum chemical units
cccbdb convert-units 1.0 Hartree kcal/mol

# Retrieve energies, geometries, or vibrations from terminal
cccbdb get H2O --property geometry --method B3LYP --basis "6-31G*" --xyz water.xyz
```

---

## Data Source, Attribution & License

- **Data Source**: Data retrieved by this package originates from the **NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB)**, NIST Standard Reference Database Number 101, Release 22 (May 2022), Editor: Russell D. Johnson III, https://cccbdb.nist.gov/, DOI: [10.18434/T47C7Z](https://doi.org/10.18434/T47C7Z).
- **Disclaimer**: This project is independent of NIST and is not affiliated with, endorsed by, or maintained by the National Institute of Standards and Technology. Users are responsible for complying with applicable NIST data access, attribution, and redistribution policies.
- **License**: The original Python source code in `cccbdb-py` is released under the [MIT License](LICENSE).
