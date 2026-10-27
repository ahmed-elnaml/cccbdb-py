# Getting Started with `cccbdb-py`

`cccbdb-py` provides a clean, typed Python interface between the [NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB, SRD 101)](https://cccbdb.nist.gov/) and computational chemistry workflows.

## 1. Installation

```bash
pip install cccbdb-py
# Optional scientific stack (pandas, numpy, ASE)
pip install "cccbdb-py[all]"
```

## 2. Context-Managed Client

Always use `CCCBDB` as a context manager (`with CCCBDB() as db:`) so HTTP sessions and optional SQLite database connections are closed cleanly:

```python
from cccbdb import CCCBDB

with CCCBDB(cache=True, delay=1.5) as db:
    water = db.get_molecule("H2O")
    print(water.id, water.cas_rn, water.name, water.point_group)
```

## 3. Caching & Offline Reproducibility

By default, `cccbdb-py` caches downloaded HTML pages under `~/.cache/cccbdb-py/pages` using deterministic SHA-256 keys of the canonical URL and molecule session context.

- Pass `offline=True` to `CCCBDB(offline=True)` to guarantee zero network requests during analysis or automated testing.
- Pass `use_sqlite=True` to automatically normalize and store retrieved species, energies, geometries, and dipoles in `~/.cache/cccbdb-py/cccbdb_normalized.sqlite`.

## 4. Physical Units & Symmetry Ontology

`cccbdb-py` includes a dedicated `cccbdb.chemistry` subpackage with CODATA 2018 unit conversions, Schoenflies point group analysis, and registries of NIST CCCBDB computational methods and basis sets:

```python
from cccbdb import (
    EnergyUnit,
    convert_energy,
    parse_point_group,
    resolve_method,
    resolve_basis_set,
    get_vibrational_scaling_factor,
)

# Convert -76.408953 Hartree to kJ/mol
print(convert_energy(-76.408953, EnergyUnit.HARTREE, EnergyUnit.KJ_MOL))

# Inspect point group symmetry and rotational symmetry number (sigma)
pg = parse_point_group("D3h")
print(pg.unicode_symbol, pg.rotational_symmetry_number, pg.allows_permanent_dipole)

# Look up method, basis set, and harmonic frequency scaling factor
print(resolve_method("CCSD(T)"))
print(resolve_basis_set("aug-cc-pVTZ"))
print(get_vibrational_scaling_factor("B3LYP", "6-31G*"))  # 0.9603
```
