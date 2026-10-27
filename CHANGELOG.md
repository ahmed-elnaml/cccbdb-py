# Changelog

All notable changes to `cccbdb-py` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-10-27

### Added
- **Core Client (`CCCBDB`)**: Context-managed client supporting species discovery (`search`, `get_molecule`), single-calculation and matrix-level property retrieval, and polite batch dataset collection (`collect`).
- **Computational Chemistry Ontology (`cccbdb.chemistry`)**:
  - Exact CODATA 2018 / IUPAC physical unit conversion engine (`Hartree`, `eV`, `kJ/mol`, `kcal/mol`, `cm^-1`, `K`, `Angstrom`, `Bohr`, `pm`, `Debye`, `C*m`, `atomic_units`, `GHz`, `MHz`, `amu*A^2`, `g*cm^2`).
  - Molecular point group parser & symmetry analyzer (`C1`, `Cs`, `Ci`, `Cn`, `Cnv`, `Cnh`, `Dn`, `Dnh`, `Dnd`, `Sn`, `Td`, `Oh`, `Ih`, `Cinfv`, `Dinfh`, `Kh`) with rotational symmetry numbers ($\sigma$) and chirality/inversion flags.
  - NIST CCCBDB computational methods registry (45+ ab initio, DFT, composite, semi-empirical, and force-field methods with numeric NIST IDs, categories, and frozen/full-core metadata).
  - NIST CCCBDB basis sets registry (45+ Pople, Dunning correlation-consistent, Ahlrichs def2, Sadlej, and ECP basis sets with numeric NIST IDs, zeta quality, diffuse/polarization flags, and Scott-Radom / NIST vibrational scaling factors).
  - Hill system chemical formula parser, exact molecular weight calculator, stoichiometry analyzer, and CAS RN checksum validator.
- **Isolated HTML Parsers (`cccbdb.parsers`)**:
  - Species & experimental summary parser (`exp2x.asp`, `alldata2x.asp`, `inchix.asp`).
  - Electronic & thermochemical energy parser (`energy2x.asp` 2D method/basis matrices and `energy3x.asp` single-calculation pages).
  - Optimized geometry parser (`geom2x.asp`, `geom3x.asp`, `expgeom2x.asp`, `energy3x.asp`) with Cartesian coordinates, Z-matrix / internal coordinates, interatomic distance matrices, and dummy atom (`X`) filtering.
  - Vibrational frequency parser (`vibs2x.asp`, `vibs3x.asp`, `energy3x.asp`) with mode symmetries, imaginary frequency (`i`) transition-state detection, IR intensities, Raman activities, and ZPE.
  - Electrostatic & rotational property parser (Mulliken/CHELPG/AIM/ESP partial charges, dipole vectors, primitive & traceless quadrupole tensors, anisotropic polarizability tensors, rotational constants $A, B, C$, and experimental thermochemistry $\Delta_f H^\circ$, $S^\circ$, $C_p$).
- **Storage & Caching (`cccbdb.storage`)**:
  - Content-addressable SHA-256 local disk cache (`DiskCache`) with TTL, offline mode, and request provenance.
  - Normalized SQLite database backend (`SQLiteStorage`) for querying and archiving benchmark datasets locally.
- **Scientific Exporters (`cccbdb.exporters`)**:
  - Validated `.xyz` Cartesian structure exporter (strictly rejecting incomplete internal-coordinate-only geometries and filtering dummy atoms).
  - CSV, JSON, pandas DataFrame, NumPy array, and optional ASE `Atoms` / RDKit `Mol` adapters.
- **Command-Line Interface (`cccbdb.cli`)**:
  - `cccbdb search`, `cccbdb get`, `cccbdb export`, `cccbdb methods`, `cccbdb basis-sets`, and `cccbdb convert-units`.
