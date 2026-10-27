# Contributing to `cccbdb-py`

Thank you for your interest in contributing to `cccbdb-py`! We welcome contributions from computational chemists, quantum chemistry software developers, cheminformaticians, and open-source contributors worldwide.

## Guiding Principles

1. **Scientific Correctness First**: Every parsed numerical value must remain strictly bound to its computational method, basis set, spin state, conformer/symmetry point group, physical unit, and NIST reference squib/DOI. Never silently coerce `None` or missing entries into `0.0`.
2. **Offline-First Testing**: All HTML parser tests must run against saved, deterministic HTML fixtures in `tests/fixtures/` so continuous integration and local development do not send unnecessary traffic to NIST servers.
3. **Polite Network Access**: Live queries must respect rate limits (`delay >= 1.5s` recommended), include a transparent `User-Agent`, and utilize local caching.
4. **Parser Isolation**: Keep HTML DOM extraction inside `src/cccbdb/parsers/`. Public models (`src/cccbdb/models/`) and client methods (`src/cccbdb/client.py`) must not depend on raw BeautifulSoup tags.

## Development Setup

```bash
git clone https://github.com/your-username/cccbdb-py.git
cd cccbdb-py
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Running Tests

```bash
# Run with pytest
pytest

# Or run with standard library unittest
python3 -m unittest discover -s tests
```

## Adding a New Property Parser

1. Save a representative NIST CCCBDB HTML snippet or page in `tests/fixtures/`.
2. Define or extend the typed Pydantic model in `src/cccbdb/models/` with explicit physical units and `DataProvenance`.
3. Implement the extraction function in `src/cccbdb/parsers/` and raise `StructureChangedError` if required table headers or landmarks are missing.
4. Add unit tests in `tests/` verifying both valid parsing and unit conversions.
