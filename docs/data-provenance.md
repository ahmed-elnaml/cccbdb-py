# Data Provenance, Units & Scientific Validation

## 1. Provenance Tracking (`DataProvenance`)

Every scientific datum extracted by `cccbdb-py` carries a `DataProvenance` object with:

- `source_url`: Exact NIST CCCBDB URL (`exp2x.asp`, `energy2x.asp`, `energy3x.asp`, etc.).
- `retrieved_at`: UTC ISO-8601 timestamp when the HTML page was parsed.
- `parser_version`: Version of `cccbdb-py` that performed the extraction.
- `raw_content_sha256`: SHA-256 cryptographic hash of the raw HTML source document.
- `original_unit`: Physical unit reported on the original NIST page (`Hartree`, `Angstrom`, `cm^-1`, `Debye`, `kJ/mol`).
- `reference_squib`: Literature citation squib (`1979HUB/HER`, `CODATA`, `Gurvich`, etc.) linked to full citations and DOIs in `ExperimentalDataSummary.references`.

## 2. Strict Handling of Missing Values & Dummy Atoms

- **No Silent Zero Coercion**: Uncomputed cells (`-`, blank cells, `N/A`) are omitted or represented as `None` (`NaN` in pandas), never `0.0`.
- **Dummy Atom (`X`) Filtering**: NIST CCCBDB Z-matrix/Cartesian listings frequently include auxiliary dummy atoms (such as `X14` in $\text{N(CH}_3)_3$ or `X4` in $\text{CH}_2$) used to define symmetry axes. `MolecularGeometry.sites` preserves the full raw table, while `MolecularGeometry.physical_sites`, `.atoms`, `.coordinates`, and `.to_xyz()` automatically filter out dummy atoms so downstream quantum chemistry codes and ASE/RDKit never encounter non-physical element symbols.
- **Imaginary Vibrational Frequencies**: Saddle points and transition states in CCCBDB report imaginary modes with an `i` suffix (e.g., `872i` in linear $\text{D}_{\infty\text{h}}$ methylene). `cccbdb-py` sets `is_imaginary=True` and stores the wavenumber as a negative float (`-872.0`), following standard quantum chemistry conventions (`VibrationalSpectrum.is_transition_state`).
