# cccbdb-py

[Tests](https://www.google.com/search?q=https://github.com/your-username/cccbdb-py/actions/workflows/tests.yml)
[](https://www.python.org/downloads/)
[](https://opensource.org/licenses/MIT)
[Typed: Pydantic](https://www.google.com/search?q=https://docs.pydantic.dev/)
[](https://github.com/astral-sh/ruff)

**An open-source Python toolkit for programmatic access, parsing, validation, unit conversion, and dataset generation from the NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB, SRD 101).**

`cccbdb-py` is a community-maintained scientific software library that bridges the [NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB, Release 22)](https://cccbdb.nist.gov/) and modern Python workflows in quantum chemistry, thermochemistry, spectroscopy, cheminformatics, and molecular machine learning. Rather than acting as a fragile script scraper, `cccbdb-py` provides a modular, strictly typed data-access architecture with CODATA 2018 unit conversions, Schoenflies point group symmetry analysis, explicit isomer disambiguation, stateful ASP session recovery, SHA-256 content-addressable caching, normalized SQLite storage, and scientific structure exporters.

> **Important Attribution & Disclaimer**: `cccbdb-py` is an independent, community-maintained client and is **not** an official NIST API, nor is it affiliated with, endorsed by, or maintained by the National Institute of Standards and Technology (NIST). All retrieved reference data originates from [NIST SRD 101](https://cccbdb.nist.gov/) (DOI: [10.18434/T47C7Z](https://www.google.com/search?q=https://doi.org/10.18434/T47C7Z)).

---

## Why `cccbdb-py`?

Researchers routinely consult [NIST CCCBDB](https://cccbdb.nist.gov/) to benchmark Density Functional Theory (DFT) functionals, *ab initio* wavefunction hierarchy ($E_{\text{HF}} \rightarrow \text{MP2} \rightarrow \text{CCSD(T)}$), basis set convergence, harmonic vibrational scaling factors, and experimental gas-phase thermochemistry ($\Delta_f H^\circ$, $S^\circ$, $C_p$). However, assembling multi-molecule benchmark sets manually through the browser is repetitive and error-prone because:

1. **Stateful Classic ASP Sessions**: Pages such as `energy2x.asp` and `vibs2x.asp` rely on server-side `ASPSESSIONID` cookies initialized on earlier search forms (`exp2x.asp` / `getformx.asp`) and return `"Sorry I lost your molecule"` when queried directly without priming.
2. **Constitutional & Geometric Isomers**: A chemical formula such as `C2H6O` (ethanol vs. dimethyl ether) or `C6H6` (benzene vs. Dewar benzene) maps to multiple distinct species that must never be conflated.
3. **Multi-Conformer & Transition-State Records**: Open-shell and flexible species (e.g., $\text{CH}_2$ triplet $C_{2v}$ minimum vs. $D_{\infty h}$ saddle point) contain multiple conformers, auxiliary Z-matrix dummy atoms (`X`), and imaginary vibrational modes (`872i`) on a single calculation page.
4. **Scientific Context & Provenance**: A raw electronic energy number without its exact method, frozen/full-core setting, basis set, point group, unit, and source URL is scientifically ambiguous.

`cccbdb-py` solves each of these challenges with a clean, researcher-friendly Python API and CLI.

---

## System Architecture & Data-Access Pipeline

```text
+-------------------------------------------------------------------------+
|                    NIST CCCBDB Web Server (SRD 101)                     |
|   exp2x.asp | energy2x.asp | energy3x.asp | geom3x.asp | vibs2x.asp     |
+-------------------------------------------------------------------------+
                                    |
                                    v
+-------------------------------------------------------------------------+
|              Transport Layer (cccbdb.transport.http)                    |
|  - Polite RateLimiter (configurable delay >= 1.5s between requests)     |
|  - Automatic ASP Session Cookie Priming & "Lost Molecule" Auto-Recovery |
|  - Exponential Backoff Retries on Transient Network Errors              |
+-------------------------------------------------------------------------+
                                    |
                     +--------------+--------------+
                     |                             |
                     v                             v
+-------------------------------------+ +---------------------------------+
|   Storage Layer (cccbdb.storage)    | |   Isolated HTML DOM Parsers     |
| - Content-Addressable SHA-256       | |   (cccbdb.parsers.*)            |
|   DiskCache (~/.cache/cccbdb-py)    | | - Species & Identifiers         |
| - Strict Offline Mode (offline=True)| | - 2D Method x Basis Matrices    |
| - Normalized SQLite Database        | | - Single Calculation Dossiers   |
|   (species, energies, geometries,   | | - Cartesians & Z-Matrices       |
|    dipoles tables)                  | | - Vibrational Spectra & ZPE     |
+-------------------------------------+ | - Charges, Dipoles, Quadrupoles |
                                        +---------------------------------+
                                                          |
                                                          v
+-------------------------------------------------------------------------+
|          Computational Chemistry Ontology (cccbdb.chemistry)            |
| - Units: CODATA 2018 (Hartree, eV, kJ/mol, kcal/mol, cm^-1, K, Bohr, D) |
| - Symmetry: Schoenflies Point Groups, Rotational Symmetry Number (sigma)|
| - Methods (45+) & Basis Sets (45+) + Empirical Vib. Scaling Factors     |
| - Hill System Formula Parser, Molecular Weights & CAS RN Checksums      |
+-------------------------------------------------------------------------+
                                    |
                                    v
+-------------------------------------------------------------------------+
|                 Public Typed API & CLI (cccbdb.CCCBDB)                  |
|   pandas DataFrame | NumPy | JSON | CSV + Manifest | XMol .xyz | ASE    |
+-------------------------------------------------------------------------+
```

---

## Comprehensive Feature Matrix

### 1. Supported Data Queries & NIST CCCBDB Endpoints

| Query Domain | NIST Endpoint(s) | Python API Method | Extracted Scientific Quantities & Metadata |
| :--- | :--- | :--- | :--- |
| **Species Discovery & Identifiers** | `getformx.asp`, `exp2x.asp`, `inchix.asp` | `db.search(...)`, `db.get_molecule(...)` | CAS RN (with check-digit verification), Hill formula, IUPAC name, synonyms, InChI, InChIKey, SMILES, net charge, spin state, ground Schoenflies point group, exact molecular weight ($\text{g/mol}$). |
| **2D Electronic Energy Matrix** | `energy2x.asp` | `mol.get_energies(...)` | Full method $\times$ basis-set matrix of 0 K total electronic energies in $\text{Hartree}$ ($E_{\text{h}}$), plus composite and semi-empirical predefined-basis energies. |
| **Single-Calculation Dossier** | `energy3x.asp` | `mol.get_calculation(method, basis)` | Multi-conformer state table, $E(0\text{ K})$, $E(298.15\text{ K})$, $E_{\text{HF}}$, $E_{\text{nuc}}$, Cartesian geometry, bond angles, harmonic frequencies, ZPE, rotational constants ($A, B, C$), Mulliken/CHELPG/AIM/ESP charges, dipole vector, primitive/traceless quadrupole tensors, polarizability tensor, and $\langle r^2 \rangle$. |
| **Molecular Geometries** | `exp2x.asp`, `energy3x.asp`, `geom3x.asp` | `mol.get_geometry(method, basis)` | 3D Cartesian coordinates ($\AA$, $\text{Bohr}$, $\text{pm}$, $\text{nm}$), automatic filtering of Z-matrix dummy atoms (`X`), internal coordinates ($r_e, r_0$, angles, dihedrals), $N \times N$ distance matrix, and mass-weighted center of mass. |
| **Vibrational Spectroscopy & ZPE** | `exp2x.asp`, `energy3x.asp`, `vibs2x.asp` | `mol.get_vibrations(method, basis)` | Normal mode indices, Mulliken symmetry labels ($\Sigma_g, \Pi_u, A_1, B_2, E$), harmonic and fundamental wavenumbers ($\text{cm}^{-1}$), empirical Scott-Radom scaled frequencies, IR intensities ($\text{km/mol}$), Raman activities ($\AA^4/\text{u}$), depolarization ratios, transition-state imaginary mode (`872i` $\rightarrow -872.0\text{ cm}^{-1}$) detection, and unscaled/scaled ZPE. |
| **Electrostatics & Multipoles** | `exp2x.asp`, `energy3x.asp` | `mol.get_dipole(...)`, `calc.primary` | Dipole components ($\mu_x, \mu_y, \mu_z, \mu_{\text{total}}$ in $\text{Debye}$, $\text{C}\cdot\text{m}$, or $\text{a.u.}$), $3 \times 3$ primitive and traceless quadrupole tensors ($\text{Debye}\cdot\AA$), $3 \times 3$ polarizability tensor and isotropic mean $\alpha_{\text{iso}}$ ($\AA^3$ or $\text{Bohr}^3$), and atomic partial charges ($e$). |
| **Experimental Thermochemistry** | `exp2x.asp` | `mol.get_experimental()` | Standard enthalpy of formation $\Delta_f H^\circ(298.15\text{ K})$ and $\Delta_f H^\circ(0\text{ K})$ ($\text{kJ/mol}$), standard molar entropy $S^\circ(298.15\text{ K})$ ($\text{J}\cdot\text{K}^{-1}\cdot\text{mol}^{-1}$), integrated heat capacity $H(298.15\text{ K}) - H(0\text{ K})$, heat capacity $C_p(298.15\text{ K})$, experimental uncertainties, and bibliographic DOIs/squibs. |

---

### 2. Computational Chemistry Domain Ontology (`cccbdb.chemistry`)

#### A. Physical Units & Exact CODATA 2018 Conversions (`cccbdb.chemistry.units`)

| Physical Quantity | Canonical Base Unit | Supported Target Units & Aliases | Exact CODATA 2018 Conversion Factor |
| :--- | :--- | :--- | :--- |
| **Energy / Thermochemistry** | `Hartree` ($E_{\text{h}}$) | `kJ/mol`, `kcal/mol`, `eV`, `cm^-1`, `K`, `J`, `MHz` | $1\ E_{\text{h}} = 2625.49963948\ \text{kJ/mol} = 627.50947406\ \text{kcal/mol} = 27.21138625\ \text{eV} = 219474.63136\ \text{cm}^{-1}$ |
| **Length / Coordinates** | `Angstrom` ($\AA$) | `Bohr` ($a_0$), `pm`, `nm` | $1\ a_0 = 0.529177210903\ \AA = 52.91772109\ \text{pm}$ |
| **Electric Dipole Moment** | `Debye` ($\text{D}$) | `au` ($e a_0$), `C*m` | $1\ e a_0 = 2.541746473\ \text{D}$, $1\ \text{D} = 3.33564095 \times 10^{-30}\ \text{C}\cdot\text{m}$ |
| **Rotational Constants** | `cm^-1` | `GHz`, `MHz` | $1\ \text{cm}^{-1} = 29.9792458\ \text{GHz} = 29979.2458\ \text{MHz}$ |
| **Moment of Inertia** | `amu*Angstrom^2` | `g*cm^2`, `au` ($m_e a_0^2$) | $1\ \text{amu}\cdot\AA^2 = 1.6605390666 \times 10^{-40}\ \text{g}\cdot\text{cm}^2$ |
| **Polarizability Volume** | `Angstrom^3` ($\AA^3$) | `Bohr^3` ($a_0^3$) | $1\ \AA^3 = 6.74833449\ a_0^3$ |

#### B. Molecular Point Group Symmetry & Statistical Thermodynamics (`cccbdb.chemistry.symmetry`)

Normalizes raw CCCBDB ASCII point group labels (`C*V` $\rightarrow$ `Cinfv`, `D*H` $\rightarrow$ `Dinfh`, `C2V` $\rightarrow$ `C2v`) and provides symmetry metadata required for rotational partition functions ($q_{\text{rot}} \propto 1/\sigma$) and spectroscopic selection rules:

| Point Group Class | Supported Schoenflies Symbols | Rotational Symmetry Number ($\sigma$) | Inversion Center ($i$) | Permanent Dipole Allowed | Chirality |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Non-Axial** | $C_1, C_s, C_i$ | $\sigma = 1$ | $C_i$ only | $C_1, C_s$ | $C_1$ only |
| **Linear** | $C_{\infty v}, D_{\infty h}$ | $\sigma = 1$ ($C_{\infty v}$), $\sigma = 2$ ($D_{\infty h}$) | $D_{\infty h}$ only | $C_{\infty v}$ only | Achiral |
| **Cyclic / Pyramidal** | $C_n, C_{nv}, C_{nh}$ | $\sigma = n$ | Even $n$ in $C_{nh}$ | $C_n, C_{nv}$ | $C_n$ only |
| **Dihedral** | $D_n, D_{nh}, D_{nd}$ | $\sigma = 2n$ | Even $n$ ($D_{nh}$), Odd $n$ ($D_{nd}$) | Forbidden ($\mu = 0$) | $D_n$ only |
| **Improper Rotation** | $S_n$ (even $n$) | $\sigma = n/2$ | $n \equiv 2 \pmod 4$ | Forbidden ($\mu = 0$) | Achiral |
| **Cubic & Icosahedral** | $T, T_d, T_h, O, O_h, I, I_h$ | $\sigma = 12$ ($T/T_d/T_h$), $24$ ($O/O_h$), $60$ ($I/I_h$) | $T_h, O_h, I_h$ | Forbidden ($\mu = 0$) | $T, O, I$ |
| **Spherical (Atoms)** | $K_h$ | $\sigma = 1$ | Yes | Forbidden ($\mu = 0$) | Achiral |

#### C. Quantum Chemical Methods Registry (`cccbdb.chemistry.methods`)

Includes 45+ computational methods mapped to their numeric NIST CCCBDB `method` IDs, core-correlation settings (`frozen` vs. `full`), and grid specifications:

- **Hartree-Fock**: `HF` (RHF/UHF, ID `1`), `ROHF` (ID `9`).
- **Density Functional Theory (Local, GGA, Meta-GGA, Hybrid, Range-Separated)**: `LSDA` (`51`), `BLYP` (`13`), `B1B95` (`24`), `B3LYP` (`8`), `B3LYPultrafine` (`68`), `B3PW91` (`14`), `mPW1PW91` (`25`), `M06-2X` (`55`), `PBEPBE` (`52`), `PBEPBEultrafine` (`69`), `PBE1PBE` / `PBE0` (`54`), `HSEh1PBE` / `HSE06` (`56`), `TPSSh` (`58`), `wB97X-D` (`60`), `B97D3` (`66`), `B97D3ultrafine` (`74`).
- **Double-Hybrid DFT**: `B2PLYP` (`59`), `B2PLYP=FULL` (`67`), `B2PLYP=FULLultrafine` (`75`).
- **Møller-Plesset Perturbation Theory**: `MP2` (`3`), `MP2=FULL` (`2`), `PMP2` (`10`), `ROMP2` (`11`), `MP3` (`36`), `MP3=FULL` (`37`), `MP4` (`4`), `MP4=FULL` (`44`).
- **Configuration Interaction & Quadratic CI**: `CID` (`16`), `CID=FULL` (`45`), `CISD` (`15`), `QCISD` (`5`), `QCISD(T)` (`22`), `QCISD(T)=FULL` (`61`), `QCISD(TQ)` (`64`), `QCISD(TQ)=FULL` (`65`).
- **Coupled Cluster Theory**: `CCD` (`6`), `CCSD` (`7`), `CCSD=FULL` (`62`), `CCSD(T)` (`12`), `CCSD(T)=FULL` (`63`).
- **Composite Thermochemical Protocols** *(predefined basis sets)*: `G1` (`19`), `G2` (`20`), `G2MP2` (`21`), `G3` (`23`), `G3B3` (`28`), `G3MP2` (`29`), `G3MP2B3` (`30`), `G4` (`53`), `G4MP2` (`57`), `CBS-Q` (`26`), `CBS-QB3` (`27`), `W1U` (`70`).
- **Semi-Empirical MO Methods**: `AM1` (`17`), `PM3` (`18`), `PM6` (`50`), `MNDOd` (`40`).

#### D. Gaussian Basis Sets & ECP Registry (`cccbdb.chemistry.basis_sets`)

Maps all 45 numeric basis set indices from `selectbasis1.asp` along with empirical harmonic vibrational scaling factors (`vsfx.asp`):

- **Pople Split-Valence**: `STO-3G` (`20`), `3-21G` (`2`), `3-21G*` (`8`), `6-31G` (`9`), `6-31G*` (`1`), `6-31G**` (`3`), `6-31+G**` (`4`), `6-311G*` (`7`), `6-311G**` (`21`), `6-31G(2df,p)` (`19`), `6-311+G(3df,2p)` (`10`), `6-311+G(3df,2pd)` (`22`).
- **Dunning Correlation-Consistent**: `cc-pVDZ` (`5`), `cc-pVTZ` (`6`), `cc-pVQZ` (`25`), `cc-pV5Z` (`39`), `aug-cc-pVDZ` (`17`), `aug-cc-pVTZ` (`18`), `aug-cc-pVQZ` (`26`), `aug-cc-pV5Z` (`40`).
- **Tight-$d$ Second-Row Correlation-Consistent**: `cc-pV(D+d)Z` (`31`), `cc-pV(T+d)Z` (`27`), `cc-pV(Q+d)Z` (`34`), `aug-cc-pV(T+d)Z` (`28`), `aug-cc-p(Q+d)Z` (`45`).
- **Core-Valence & Doubly-Augmented Dunning**: `cc-pCVDZ` (`23`), `cc-pCVTZ` (`24`), `cc-pCVQZ` (`30`), `aug-cc-pCVTZ` (`33`), `aug-cc-pCVQZ` (`44`), `aug-cc-pCV5Z` (`43`), `daug-cc-pVDZ` (`36`), `daug-cc-pVTZ` (`35`), `daug-cc-pVQZ` (`37`).
- **Ahlrichs & Sadlej**: `TZVP` (`32`), `Def2TZVPP` (`38`), `Sadlej_pVTZ` (`29`).
- **Effective Core Potentials (ECPs)**: `CEP-31G` (`11`), `CEP-31G*` (`12`), `CEP-121G` (`13`), `CEP-121G*` (`14`), `LANL2DZ` (`15`), `SDD` (`16`), `cc-pVTZ-PP` (`41`), `aug-cc-pVTZ-PP` (`42`).

---

## Installation

```bash
# 1. Core installation (httpx, beautifulsoup4, lxml, pydantic)
pip install cccbdb-py

# 2. With scientific tabular & array stack (pandas, numpy)
pip install "cccbdb-py[data]"

# 3. With computational chemistry structure adapters (ASE, RDKit)
pip install "cccbdb-py[chemistry]"

# 4. Full contributor & development installation from local clone
git clone https://github.com/your-username/cccbdb-py.git
cd cccbdb-py
pip install -e ".[dev]"
```

---

## Detailed Python API Guide

### 1. Molecular Discovery & Explicit Isomer Disambiguation

```python
from cccbdb import CCCBDB, AmbiguousSpeciesError

with CCCBDB() as db:
    # Search by Hill or structural formula (returns all matching isomers)
    isomers = db.search(formula="C2H6O")
    for mol in isomers:
        print(f"{mol.id:<8} | {mol.cas_rn:<10} | {mol.formula:<10} | {mol.name:<15} | {mol.point_group}")
    # 64175    | 64-17-5    | CH3CH2OH   | Ethanol         | Cs
    # 115106   | 115-10-6   | CH3OCH3    | Dimethyl ether  | C2v

    # Calling get_molecule() on an ambiguous formula raises AmbiguousSpeciesError
    try:
        db.get_molecule("C2H6O")
    except AmbiguousSpeciesError as err:
        print("Ambiguous formula:", err.candidates)

    # Resolve unambiguously via CAS RN, species ID, or exact name
    ethanol = db.get_molecule("64-17-5")
    benzene = db.get_molecule("Benzene")
```

### 2. Electronic & Thermochemical Energies (`Hartree`, `kJ/mol`, `kcal/mol`, `eV`)

```python
from cccbdb import CCCBDB, EnergyUnit

with CCCBDB() as db:
    water = db.get_molecule("H2O")

    # Retrieve B3LYP/6-31G* energies (0K, 298.15K, HF, Nuclear Repulsion)
    energies = water.get_energies(method="B3LYP", basis_set="6-31G*")
    df = energies.to_dataframe()
    print(df[["method", "basis_set", "energy_type", "value", "unit"]])

    # Convert entire collection to kJ/mol or kcal/mol
    energies_kcal = energies.convert_units(EnergyUnit.KCAL_MOL)
    for rec in energies_kcal:
        print(f"{rec.energy_type:<18}: {rec.value:14.4f} {rec.unit} (source: {rec.provenance.source_url})")
```

### 3. Optimized Molecular Geometries, `.xyz` Export & ASE Integration

`cccbdb-py` automatically strips non-physical Z-matrix dummy atoms (`X`) when accessing `.atoms`, `.coordinates`, `.physical_sites`, `.to_xyz()`, or `.to_ase()`, and raises `IncompleteGeometryError` if a record only contains internal bond lengths/angles without 3D Cartesian coordinates:

```python
from cccbdb import CCCBDB, LengthUnit

with CCCBDB() as db:
    water = db.get_molecule("H2O")
    geom = water.get_geometry(method="B3LYP", basis_set="6-31G*")

    print("Physical atoms:", geom.atoms)               # ['O', 'H', 'H']
    print("Coordinates (A):", geom.coordinates)        # [(0.0, 0.0, 0.1197), ...]
    print("O-H1 bond length (A):", geom.distance(0, 1))
    print("Center of mass (A):", geom.center_of_mass())
    print("Distance matrix (A):", geom.distance_matrix())

    # Convert Cartesian coordinates to atomic units (Bohr)
    geom_bohr = geom.convert_units(LengthUnit.BOHR)
    print("O-H1 bond length (Bohr):", geom_bohr.distance(0, 1))

    # Export validated XMol .xyz file (always standardized in Angstroms)
    geom.to_xyz("water_b3lyp_631gs.xyz")

    # Optional: Convert to an ASE Atoms object (requires `pip install ase`)
    # atoms = geom.to_ase()
```

### 4. Vibrational Frequencies, Transition States & Zero-Point Energy (ZPE)

```python
from cccbdb import CCCBDB, EnergyUnit

with CCCBDB() as db:
    water = db.get_molecule("H2O")
    vibs = water.get_vibrations(method="B3LYP", basis_set="6-31G*")

    print(vibs.to_dataframe())
    print("Is transition state (imaginary modes > 0)?", vibs.is_transition_state)
    print("Scaling factor:", vibs.scaling_factor)      # 0.9603
    print("Scaled ZPE (cm^-1):", vibs.zpe_scaled_cm_inv)
    print("Scaled ZPE (kcal/mol):", vibs.zpe_in(EnergyUnit.KCAL_MOL, scaled=True))
    print("Scaled ZPE (Hartree):", vibs.zpe_in(EnergyUnit.HARTREE, scaled=True))
```

### 5. Single-Calculation Full Dossier (Charges, Dipoles, Quadrupoles, Polarizability)

```python
from cccbdb import CCCBDB, DipoleUnit

with CCCBDB() as db:
    water = db.get_molecule("H2O")
    calc = water.get_calculation(method="B3LYP", basis_set="6-31G*")
    conf = calc.primary

    # Partial atomic charges (Mulliken, CHELPG, AIM, ESP)
    for chg in conf.atomic_charges:
        print(f"Atom {chg.atom_index} ({chg.element}): Mulliken={chg.mulliken}, CHELPG={chg.chelpg}, AIM={chg.aim}")

    # Electric dipole moment in Debye and atomic units
    print("Dipole (Debye):", conf.dipole.total)
    print("Dipole (a.u.):", conf.dipole.to_unit(DipoleUnit.ATOMIC_UNIT).total)

    # Traceless quadrupole tensor (Debye*Angstrom) & polarizability tensor (Angstrom^3)
    print("Traceless Quadrupole zz:", conf.quadrupoles["Traceless"].zz)
    print("Isotropic Polarizability (A^3):", conf.polarizability.isotropic)
    print("Electronic spatial extent <r^2> (A^2):", conf.r2_angstrom2)
```

### 6. Batch Dataset Collection, Local Caching & SQLite Persistence

```python
from cccbdb import CCCBDB

with CCCBDB(cache=True, use_sqlite=True, delay=1.5) as db:
    molecules = db.search(elements=["C", "H", "O"])
    dataset = db.collect(
        molecules=molecules[:10],
        properties=["energy", "dipole", "zpe"],
        method="B3LYP",
        basis_set="6-31G*",
    )
    # Exports CSV + companion JSON reproducibility manifest (.csv.manifest.json)
    dataset.to_csv("cho_b3lyp_benchmark.csv")

    # Query normalized local SQLite database directly
    rows = db.sqlite.query_energies(method="B3LYP", basis_set="6-31G*")
    print(f"Stored {len(rows)} energy records in local SQLite database.")
```

---

## Command-Line Interface (`cccbdb`)

The `cccbdb` terminal tool provides immediate access to the ontology, unit converter, search engine, and data exporters without writing Python code:

```bash
# 1. Search species by formula, name, CAS RN, point group, or element set
cccbdb search --formula C2H6O
cccbdb search --point-group D6h
cccbdb search --elements C,H,O --json

# 2. Inspect molecular point group symmetry and rotational symmetry numbers (sigma)
cccbdb symmetry C3v
cccbdb symmetry Dinfh

# 3. Browse supported NIST CCCBDB quantum chemistry methods and basis sets
cccbdb methods --category "Coupled Cluster"
cccbdb methods --category "Density Functional"
cccbdb basis-sets --family Dunning
cccbdb basis-sets --family Pople

# 4. Convert physical chemistry quantities on the command line
cccbdb convert-units 1.0 Hartree kJ/mol
cccbdb convert-units -76.408953 Hartree kcal/mol
cccbdb convert-units 1.0 Bohr Angstrom
cccbdb convert-units 1.855 Debye au

# 5. Retrieve molecular properties or export .xyz structures
cccbdb get H2O --property summary
cccbdb get H2O --property energies --method B3LYP --basis "6-31G*" --unit kJ/mol
cccbdb get H2O --property geometry --method B3LYP --basis "6-31G*" --xyz water.xyz
cccbdb get H2O --property vibrations --method B3LYP --basis "6-31G*"

# 6. Export a multi-molecule dataset to CSV or JSON
cccbdb export --elements C,H,O --method Experimental --properties dipole,zpe,enthalpy -o exp_benchmark.csv
```

---

## Configuration & Environment Variables

`CCCBDBConfig` can be configured either via `CCCBDB(...)` constructor arguments or environment variables:

| Environment Variable | Constructor Parameter | Default Value | Description |
| :--- | :--- | :--- | :--- |
| `CCCBDB_BASE_URL` | `base_url` | `https://cccbdb.nist.gov` | Base URL for the NIST CCCBDB web server. |
| `CCCBDB_CACHE` | `cache` | `1` (`True`) | Enable or disable content-addressable local HTML caching. |
| `CCCBDB_CACHE_DIR` | `cache_dir` | `~/.cache/cccbdb-py` | Directory for cached HTML pages (`pages/`) and SQLite DB. |
| `CCCBDB_OFFLINE` | `offline` | `0` (`False`) | When `1`/`True`, never make network calls; raise `CacheMissError` on uncached URLs. |
| `CCCBDB_DELAY` | `delay` | `1.5` (seconds) | Polite minimum interval between consecutive HTTP requests. |
| `CCCBDB_TIMEOUT` | `timeout` | `30.0` (seconds) | HTTP socket timeout per request. |

---

## Repository Structure

```text
cccbdb-py/
├── src/
│   └── cccbdb/
│       ├── __init__.py              # Public package exports
│       ├── client.py                # High-level CCCBDB context-managed client
│       ├── cli.py                   # Command-line interface (`cccbdb`)
│       ├── config.py                # CCCBDBConfig & environment variable loader
│       ├── exceptions.py            # Custom scientific & transport exception hierarchy
│       ├── chemistry/
│       │   ├── units.py             # CODATA 2018 unit conversion engine
│       │   ├── symmetry.py          # Schoenflies point group & rotational symmetry (sigma)
│       │   ├── methods.py           # 45+ quantum chemical methods registry & NIST IDs
│       │   ├── basis_sets.py        # 45+ Gaussian basis sets, ECPs & vib. scaling factors
│       │   └── formula.py           # Hill formula parser, atomic weights & CAS checksums
│       ├── transport/
│       │   ├── http.py              # Stateful httpx transport & ASP session recovery
│       │   └── rate_limit.py        # Thread-safe polite RateLimiter
│       ├── parsers/
│       │   ├── _utils.py            # BeautifulSoup/lxml table helpers & imaginary mode parser
│       │   ├── species.py           # Species discovery & exp2x.asp parser
│       │   ├── energies.py          # energy2x.asp 2D matrix & energy3x.asp dossier parsers
│       │   ├── geometries.py        # Cartesian, Z-matrix & internal coordinate parsers
│       │   ├── frequencies.py       # Harmonic/fundamental vibrational spectra & ZPE parsers
│       │   └── properties.py        # Charges, dipoles, quadrupoles, polarizability & thermo
│       ├── models/
│       │   ├── molecule.py          # SpeciesSummary, ExperimentalDataSummary, Molecule, BenchmarkDataset
│       │   ├── calculation.py       # CalculationResult & ConformerCalculation
│       │   ├── geometry.py          # AtomCoordinate, InternalCoordinate & MolecularGeometry
│       │   ├── property.py          # EnergyRecord, VibrationalSpectrum, DipoleMoment, etc.
│       │   └── provenance.py        # DataProvenance (SHA-256, timestamp, URL) & ReferenceCitation
│       ├── storage/
│       │   ├── cache.py             # SHA-256 content-addressable local DiskCache
│       │   └── sqlite.py            # Normalized SQLite relational database backend
│       └── exporters/
│           ├── csv.py               # CSV + .manifest.json reproducibility exporter
│           ├── json.py              # Structured JSON exporter
│           ├── xyz.py               # Validated XMol .xyz Cartesian structure exporter
│           └── ase_export.py        # Optional ASE Atoms structure adapter
├── tests/
│   ├── fixtures/                    # Saved NIST CCCBDB HTML fixtures for 100% offline testing
│   ├── test_chemistry_ontology.py   # Unit tests for units, point groups, methods, basis sets
│   ├── test_species.py              # Species discovery & isomer disambiguation tests
│   ├── test_energies.py             # 2D matrix & single-calculation energy tests
│   ├── test_geometries_and_vibrations.py # XYZ validation, dummy atoms & imaginary mode tests
│   └── test_client_and_storage.py   # End-to-end client, cache, SQLite & CLI integration tests
├── examples/
│   ├── search_molecules.py
│   ├── retrieve_energies.py
│   └── export_dataset.py
├── docs/
│   ├── getting-started.md
│   ├── api-reference.md
│   └── data-provenance.md
├── .github/workflows/
│   ├── tests.yml                    # Cross-platform CI (Ubuntu, macOS, Windows; Python 3.9–3.12)
│   └── publish.yml                  # Trusted Publisher PyPI release workflow
├── pyproject.toml
├── README.md
├── LICENSE
├── CONTRIBUTING.md
└── CHANGELOG.md
```

---

## Testing & Verification

All unit and integration tests run offline against deterministic HTML fixtures in `tests/fixtures/` so that continuous integration never burdens NIST servers:

```bash
pytest
# Or with coverage report:
pytest --cov=cccbdb --cov-report=term-missing
```

---

## How to Cite

If you use `cccbdb-py` in published computational chemistry or machine learning research, please cite **both** the underlying NIST database (as required by NIST SRD 101) and this software package:

```bibtex
@misc{nist_cccbdb_srd101,
  author       = {Johnson III, Russell D.},
  title        = {{NIST Computational Chemistry Comparison and Benchmark Database, NIST Standard Reference Database Number 101}},
  howpublished = {Release 22, National Institute of Standards and Technology},
  year         = {2022},
  month        = may,
  url          = {https://cccbdb.nist.gov/},
  doi          = {10.18434/T47C7Z}
}
```

---

## License

The original Python source code of `cccbdb-py` is released under the [MIT License](https://www.google.com/search?q=LICENSE). Data retrieved from the NIST Computational Chemistry Comparison and Benchmark Database (SRD 101) is subject to [NIST copyright, fair use, and data attribution policies](https://www.nist.gov/open/copyright-fair-use-and-licensing-statements-srd-data-software-and-technical-series-publications).
