# API Reference

## Client (`cccbdb.CCCBDB`)

### `CCCBDB(cache=True, cache_dir=None, offline=False, delay=1.5, timeout=30.0, use_sqlite=False)`
Main client entry point.

- `db.search(formula=None, name=None, cas=None, inchi=None, smiles=None, point_group=None, elements=None, charge=None) -> list[Molecule]`
- `db.get_molecule(identifier: str, charge: int = 0) -> Molecule`
- `db.get_energies(species_id, charge=0, method=None, basis_set=None, energy_type=None, unit="Hartree") -> EnergyCollection`
- `db.get_geometry(species_id, charge=0, method="Experimental", basis_set=None, unit="Angstrom") -> MolecularGeometry`
- `db.get_vibrations(species_id, charge=0, method="Experimental", basis_set=None) -> VibrationalSpectrum`
- `db.get_dipole(species_id, charge=0, method="Experimental", basis_set=None, unit="Debye") -> Optional[DipoleMoment]`
- `db.get_calculation(species_id, method, basis_set=None, charge=0) -> CalculationResult`
- `db.get_experimental_data(species_id, charge=0) -> ExperimentalDataSummary`
- `db.collect(molecules, properties=("energy", "dipole", "zpe"), method="B3LYP", basis_set="6-31G*", delay=None) -> BenchmarkDataset`

---

## Models (`cccbdb.models`)

- **`Molecule`**: Interactive species handle with `.get_energies()`, `.get_geometry()`, `.get_vibrations()`, `.get_dipole()`, `.get_calculation()`, and `.get_experimental()`.
- **`EnergyCollection` / `EnergyRecord`**: Supports `.filter(method=..., basis_set=..., energy_type=...)`, `.convert_units(target_unit)`, `.to_dataframe()`, and `.to_dicts()`.
- **`MolecularGeometry` / `AtomCoordinate` / `InternalCoordinate`**: Supports `.atoms`, `.coordinates`, `.physical_sites` (automatically excluding dummy `X` atoms), `.distance(i, j)`, `.distance_matrix()`, `.center_of_mass()`, `.convert_units(target_unit)`, `.to_xyz(path)`, `.to_numpy()`, `.to_ase()`, and `.to_dataframe()`.
- **`VibrationalSpectrum` / `VibrationalMode`**: Supports `.modes`, `.imaginary_modes`, `.is_transition_state`, `.zpe_in(unit, scaled=True)`, and `.to_dataframe()`.
- **`DipoleMoment`**, **`QuadrupoleMoment`**, **`PolarizabilityTensor`**, **`AtomicChargeRecord`**, **`RotationalConstants`**, **`ThermochemicalProperty`**.

---

## Exceptions (`cccbdb.exceptions`)

- `CCCBDBError`: Base exception class.
- `SpeciesNotFoundError`: Raised when no species matches a query.
- `AmbiguousSpeciesError`: Raised when `get_molecule()` is called with an ambiguous formula matching multiple isomers.
- `CalculationNotFoundError`: Raised when a requested `method`/`basis_set` combination is not present for a molecule.
- `IncompleteGeometryError`: Raised when attempting `.to_xyz()` or `.to_ase()` on a record lacking 3D Cartesian coordinates.
- `StructureChangedError`: Raised if NIST CCCBDB HTML layout changes break parser landmarks.
- `SessionLostError`: Raised if NIST CCCBDB loses server-side ASP session state after retry.
- `CacheMissError`: Raised in `offline=True` mode when a page is not cached locally.
