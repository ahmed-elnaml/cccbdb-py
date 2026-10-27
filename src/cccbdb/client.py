"""High-level CCCBDB Python client for computational chemistry workflows."""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import httpx

from cccbdb.chemistry.basis_sets import BasisSetSpec, list_basis_sets, resolve_basis_set
from cccbdb.chemistry.formula import (
    calculate_molecular_weight,
    format_cas_hyphenated,
    normalize_cas,
    parse_formula,
    to_hill_formula,
)
from cccbdb.chemistry.methods import ComputationalMethodSpec, list_methods, resolve_method
from cccbdb.chemistry.symmetry import normalize_point_group_symbol
from cccbdb.chemistry.units import DipoleUnit, EnergyUnit, LengthUnit
from cccbdb.config import CCCBDBConfig
from cccbdb.exceptions import (
    AmbiguousSpeciesError,
    CalculationNotFoundError,
    CCCBDBError,
    SpeciesNotFoundError,
)
from cccbdb.models.calculation import CalculationResult
from cccbdb.models.geometry import MolecularGeometry
from cccbdb.models.molecule import (
    BenchmarkDataset,
    DatasetRecord,
    ExperimentalDataSummary,
    Molecule,
    SpeciesSummary,
)
from cccbdb.models.property import (
    DipoleMoment,
    EnergyCollection,
    EnergyRecord,
    VibrationalSpectrum,
)
from cccbdb.parsers.energies import parse_calculation_result_page, parse_energies_matrix_page
from cccbdb.parsers.species import (
    BUILTIN_SPECIES_CATALOG,
    parse_experimental_page,
    parse_species_disambiguation_page,
)
from cccbdb.storage.cache import DiskCache
from cccbdb.storage.sqlite import SQLiteStorage
from cccbdb.transport.http import CCCBDBTransport


class CCCBDB:
    """Primary interface for discovering molecules and retrieving NIST CCCBDB data.

    Examples
    --------
    >>> from cccbdb import CCCBDB
    >>> with CCCBDB() as db:
    ...     molecules = db.search(formula="H2O")
    ...     for mol in molecules:
    ...         print(mol.id, mol.name)
    """

    def __init__(
        self,
        cache: bool = True,
        cache_dir: Optional[Union[str, Path]] = None,
        offline: bool = False,
        delay: float = 1.5,
        timeout: float = 30.0,
        base_url: str = "https://cccbdb.nist.gov",
        use_sqlite: bool = False,
        http_client: Optional[httpx.Client] = None,
    ) -> None:
        self.config = CCCBDBConfig(
            base_url=base_url,
            cache_enabled=cache,
            cache_dir=Path(cache_dir) if cache_dir else (Path.home() / ".cache" / "cccbdb-py"),
            offline=offline,
            rate_limit_delay=delay,
            timeout_seconds=timeout,
        )
        self.disk_cache = (
            DiskCache(self.config.cache_dir, ttl_seconds=self.config.cache_ttl_seconds)
            if self.config.cache_enabled
            else None
        )
        self.transport = CCCBDBTransport(
            config=self.config,
            cache=self.disk_cache,
            http_client=http_client,
        )
        self.sqlite = SQLiteStorage(self.config.sqlite_path) if use_sqlite else None
        self._species_index: List[SpeciesSummary] = list(BUILTIN_SPECIES_CATALOG)

    def register_species(self, species: SpeciesSummary) -> None:
        """Register or update a :class:`SpeciesSummary` in the client's local catalog."""
        self._species_index = [
            s for s in self._species_index if not (s.id == species.id and s.charge == species.charge)
        ]
        self._species_index.append(species)
        if self.sqlite:
            self.sqlite.save_species(species)

    def search(
        self,
        formula: Optional[str] = None,
        name: Optional[str] = None,
        cas: Optional[str] = None,
        inchi: Optional[str] = None,
        smiles: Optional[str] = None,
        point_group: Optional[str] = None,
        elements: Optional[Sequence[str]] = None,
        charge: Optional[int] = None,
    ) -> List[Molecule]:
        """Search for chemical species by formula, name, CAS RN, InChI, SMILES, point group, or element set.

        Always returns a list of :class:`Molecule` objects because a single chemical formula
        (such as ``C2H6O`` or ``C6H6``) can represent multiple constitutional or geometric isomers.
        """
        target_hill: Optional[str] = None
        if formula:
            target_hill = to_hill_formula(formula, include_charge=False)
            _, parsed_charge = parse_formula(formula)
            if charge is None and ("+" in formula or "-" in formula):
                charge = parsed_charge

        target_cas = normalize_cas(cas) if cas else None
        target_pg = normalize_point_group_symbol(point_group) if point_group else None
        target_elements = {el.strip() for el in elements} if elements else None

        matched: List[SpeciesSummary] = []
        for sp in self._species_index:
            if target_cas and sp.id != target_cas:
                continue
            if charge is not None and sp.charge != charge:
                continue
            if target_hill and sp.hill_formula.lower() != target_hill.lower():
                continue
            if name and name.strip().lower() not in sp.name.lower():
                continue
            if inchi and (not sp.inchi or inchi.strip().lower() not in sp.inchi.lower()):
                continue
            if smiles and (not sp.smiles or smiles.strip() != sp.smiles):
                continue
            if target_pg and (not sp.point_group or normalize_point_group_symbol(sp.point_group) != target_pg):
                continue
            if target_elements:
                comp, _ = parse_formula(sp.formula)
                if not set(comp.keys()).issubset(target_elements):
                    continue
            matched.append(sp)

        # If no local match and online mode is active, query NIST CCCBDB getformx.asp
        if not matched and not self.config.offline and (formula or cas or name):
            q_str = formula or cas or name or ""
            try:
                url, html = self.transport.get_page("getformx.asp", {"formula": q_str})
                if "Experimental data for" in html:
                    exp = parse_experimental_page(html, source_url=url)
                    sp_new = SpeciesSummary(
                        id=exp.species_id,
                        cas_rn=exp.cas_rn,
                        name=exp.name,
                        formula=exp.formula,
                        hill_formula=to_hill_formula(exp.formula, include_charge=False),
                        charge=exp.charge,
                        electronic_state=exp.electronic_state,
                        point_group=exp.point_group,
                        inchi=exp.inchi,
                        inchikey=exp.inchikey,
                        smiles=exp.smiles,
                        molecular_weight=calculate_molecular_weight(exp.formula),
                    )
                    self.register_species(sp_new)
                    matched.append(sp_new)
                else:
                    discovered = parse_species_disambiguation_page(html, source_url=url)
                    for sp_new in discovered:
                        self.register_species(sp_new)
                        matched.append(sp_new)
            except CCCBDBError:
                pass

        return [Molecule(summary=sp, client=self) for sp in matched]

    def get_molecule(
        self,
        identifier: str,
        charge: int = 0,
    ) -> Molecule:
        """Retrieve a unique :class:`Molecule` by CAS RN, unambiguous formula, or exact species name.

        Raises
        ------
        AmbiguousSpeciesError
            If ``identifier`` is a formula matching multiple isomers (e.g., ``C2H6O`` or ``C6H6``).
        SpeciesNotFoundError
            If no matching species exists.
        """
        cleaned = identifier.strip()

        # 1. Check if identifier is a numeric CAS RN (e.g., "7732185" or "7732-18-5")
        if re.match(r"^\d{2,7}-\d{2}-\d$", cleaned) or cleaned.isdigit():
            cas_digits = normalize_cas(cleaned)
            matches = [
                sp for sp in self._species_index if sp.id == cas_digits and sp.charge == charge
            ]
            if len(matches) == 1:
                return Molecule(summary=matches[0], client=self)
            if not self.config.offline:
                exp = self.get_experimental_data(cas_digits, charge=charge)
                sp_new = SpeciesSummary(
                    id=exp.species_id,
                    cas_rn=exp.cas_rn,
                    name=exp.name,
                    formula=exp.formula,
                    hill_formula=to_hill_formula(exp.formula, include_charge=False),
                    charge=exp.charge,
                    electronic_state=exp.electronic_state,
                    point_group=exp.point_group,
                    inchi=exp.inchi,
                    inchikey=exp.inchikey,
                    smiles=exp.smiles,
                    molecular_weight=calculate_molecular_weight(exp.formula),
                )
                self.register_species(sp_new)
                return Molecule(summary=sp_new, client=self, experimental=exp)

        # 2. Check exact name match (case-insensitive)
        name_matches = [
            sp
            for sp in self._species_index
            if sp.name.lower() == cleaned.lower() and sp.charge == charge
        ]
        if len(name_matches) == 1:
            return Molecule(summary=name_matches[0], client=self)

        # 3. Check formula match
        try:
            hill = to_hill_formula(cleaned, include_charge=False)
            formula_matches = [
                sp
                for sp in self._species_index
                if (sp.formula.lower() == cleaned.lower() or sp.hill_formula.lower() == hill.lower())
                and sp.charge == charge
            ]
            # If the user gave a structural formula like CH3CH2OH that matches exactly 1 entry, return it
            exact_structural = [
                sp for sp in formula_matches if sp.formula.lower() == cleaned.lower()
            ]
            if len(exact_structural) == 1 and len(formula_matches) > 1 and cleaned.lower() != hill.lower():
                return Molecule(summary=exact_structural[0], client=self)
            if len(formula_matches) == 1:
                return Molecule(summary=formula_matches[0], client=self)
            if len(formula_matches) > 1:
                candidates = [f"{sp.name} (id='{sp.id}', formula='{sp.formula}')" for sp in formula_matches]
                raise AmbiguousSpeciesError(cleaned, candidates=candidates)
        except ValueError:
            pass

        raise SpeciesNotFoundError(cleaned)

    def get_experimental_data(
        self,
        species_id: str,
        charge: int = 0,
        force_refresh: bool = False,
    ) -> ExperimentalDataSummary:
        """Retrieve and parse ``exp2x.asp`` experimental reference data for ``species_id``."""
        casno = normalize_cas(species_id)
        params: Dict[str, object] = {"casno": casno}
        if charge != 0:
            params["charge"] = charge
        url, html = self.transport.get_page(
            "exp2x.asp",
            params=params,
            force_refresh=force_refresh,
        )
        return parse_experimental_page(html, source_url=url)

    def get_calculation(
        self,
        species_id: str,
        method: str,
        basis_set: Optional[str] = None,
        charge: int = 0,
        force_refresh: bool = False,
    ) -> CalculationResult:
        """Retrieve a complete single-calculation dossier (`energy3x.asp`) for a species, method, and basis set."""
        casno = normalize_cas(species_id)
        m_spec = resolve_method(method)
        b_spec = (
            resolve_basis_set(basis_set)
            if not m_spec.predefined_basis and basis_set is not None
            else resolve_basis_set(None)
        )
        if not m_spec.predefined_basis and basis_set is None:
            raise ValueError(
                f"Method '{m_spec.name}' requires an explicit basis_set (e.g., '6-31G*' or 'cc-pVTZ')."
            )

        params: Dict[str, object] = {
            "casno": casno,
            "method": m_spec.nist_id,
            "basis": b_spec.nist_id,
        }
        if charge != 0:
            params["charge"] = charge

        url, html = self.transport.get_page(
            "energy3x.asp",
            params=params,
            force_refresh=force_refresh,
        )
        calc = parse_calculation_result_page(html, source_url=url, species_id=casno)
        if self.sqlite and calc.primary:
            if calc.primary.energies:
                self.sqlite.save_energies(list(calc.primary.energies.values()))
            if calc.primary.geometry:
                self.sqlite.save_geometry(calc.primary.geometry)
            if calc.primary.dipole:
                self.sqlite.save_dipole(calc.primary.dipole)
        return calc

    def get_energies(
        self,
        species_id: str,
        charge: int = 0,
        method: Optional[str] = None,
        basis_set: Optional[str] = None,
        energy_type: Optional[str] = None,
        unit: Union[str, EnergyUnit] = EnergyUnit.HARTREE,
        force_refresh: bool = False,
    ) -> EnergyCollection:
        """Retrieve electronic and thermochemical energies for a molecule.

        If both ``method`` and ``basis_set`` (or a predefined-basis composite method) are provided,
        attempts to fetch the detailed ``energy3x.asp`` page first (which includes 0K, 298.15K, HF,
        and nuclear repulsion energies), falling back to the ``energy2x.asp`` matrix page.
        """
        casno = normalize_cas(species_id)
        sp_match = next((s for s in self._species_index if s.id == casno), None)
        mol_formula = sp_match.formula if sp_match else None

        collection: Optional[EnergyCollection] = None
        if method is not None:
            m_spec = resolve_method(method)
            if m_spec.predefined_basis or basis_set is not None:
                try:
                    calc = self.get_calculation(
                        species_id=casno,
                        method=m_spec.name,
                        basis_set=basis_set,
                        charge=charge,
                        force_refresh=force_refresh,
                    )
                    recs: List[EnergyRecord] = []
                    seen_types: set[tuple[int, str]] = set()
                    for conf in calc.conformers:
                        for e_rec in conf.energies.values():
                            key = (conf.conformer_index, e_rec.energy_type)
                            if key not in seen_types:
                                seen_types.add(key)
                                recs.append(e_rec)
                    if recs:
                        collection = EnergyCollection(records=recs)
                except CCCBDBError:
                    collection = None

        if collection is None:
            url, html = self.transport.get_page(
                "energy2x.asp",
                species_context=casno,
                charge=charge,
                force_refresh=force_refresh,
            )
            collection = parse_energies_matrix_page(
                html,
                species_id=casno,
                formula=mol_formula,
                source_url=url,
            )

        filtered = collection.filter(method=method, basis_set=basis_set, energy_type=energy_type)
        if not filtered.records:
            raise CalculationNotFoundError(species_id=casno, method=method, basis_set=basis_set)

        if self.sqlite:
            self.sqlite.save_energies(filtered.records)

        target_unit_str = unit.value if isinstance(unit, EnergyUnit) else unit
        if target_unit_str.lower() not in ("hartree", "hartrees"):
            return filtered.convert_units(unit)
        return filtered

    def get_geometry(
        self,
        species_id: str,
        charge: int = 0,
        method: str = "Experimental",
        basis_set: Optional[str] = None,
        unit: Union[str, LengthUnit] = LengthUnit.ANGSTROM,
        force_refresh: bool = False,
    ) -> MolecularGeometry:
        """Retrieve an experimental or calculated optimized molecular geometry."""
        casno = normalize_cas(species_id)
        if method.strip().lower() in ("experimental", "exp"):
            exp = self.get_experimental_data(casno, charge=charge, force_refresh=force_refresh)
            if not exp.geometry:
                raise CalculationNotFoundError(species_id=casno, method="Experimental")
            geom = exp.geometry
        else:
            calc = self.get_calculation(
                species_id=casno,
                method=method,
                basis_set=basis_set,
                charge=charge,
                force_refresh=force_refresh,
            )
            if not calc.geometry:
                raise CalculationNotFoundError(species_id=casno, method=method, basis_set=basis_set)
            geom = calc.geometry

        target_unit_str = unit.value if isinstance(unit, LengthUnit) else unit
        if target_unit_str.lower() not in ("angstrom", "angstroms", "a", "å"):
            return geom.convert_units(unit)
        return geom

    def get_vibrations(
        self,
        species_id: str,
        charge: int = 0,
        method: str = "Experimental",
        basis_set: Optional[str] = None,
        force_refresh: bool = False,
    ) -> VibrationalSpectrum:
        """Retrieve experimental or calculated vibrational frequencies and ZPE."""
        casno = normalize_cas(species_id)
        if method.strip().lower() in ("experimental", "exp"):
            exp = self.get_experimental_data(casno, charge=charge, force_refresh=force_refresh)
            if not exp.vibrations:
                raise CalculationNotFoundError(species_id=casno, method="Experimental")
            return exp.vibrations

        calc = self.get_calculation(
            species_id=casno,
            method=method,
            basis_set=basis_set,
            charge=charge,
            force_refresh=force_refresh,
        )
        if not calc.vibrations:
            raise CalculationNotFoundError(species_id=casno, method=method, basis_set=basis_set)
        return calc.vibrations

    def get_dipole(
        self,
        species_id: str,
        charge: int = 0,
        method: str = "Experimental",
        basis_set: Optional[str] = None,
        unit: Union[str, DipoleUnit] = DipoleUnit.DEBYE,
        force_refresh: bool = False,
    ) -> Optional[DipoleMoment]:
        """Retrieve the experimental or calculated electric dipole moment."""
        casno = normalize_cas(species_id)
        if method.strip().lower() in ("experimental", "exp"):
            exp = self.get_experimental_data(casno, charge=charge, force_refresh=force_refresh)
            dip = exp.dipole
        else:
            calc = self.get_calculation(
                species_id=casno,
                method=method,
                basis_set=basis_set,
                charge=charge,
                force_refresh=force_refresh,
            )
            dip = calc.dipole

        if dip is None:
            return None
        target_unit_str = unit.value if isinstance(unit, DipoleUnit) else unit
        if target_unit_str.lower() not in ("debye", "d"):
            return dip.to_unit(unit)
        return dip

    def collect(
        self,
        molecules: Sequence[Union[Molecule, SpeciesSummary, str]],
        properties: Sequence[str] = ("energy", "dipole", "zpe"),
        method: str = "B3LYP",
        basis_set: Optional[str] = "6-31G*",
        delay: Optional[float] = None,
    ) -> BenchmarkDataset:
        """Collect a reproducible multi-molecule benchmark dataset across requested properties.

        Respects polite rate limits, caches downloaded pages, and records partial failures
        explicitly rather than silently skipping failed queries or substituting 0.0.
        """
        old_delay = self.transport.rate_limiter.delay_seconds
        if delay is not None:
            self.transport.rate_limiter.delay_seconds = max(0.0, float(delay))

        dataset_records: List[DatasetRecord] = []
        partial_failures: Dict[str, str] = {}
        req_props = [p.strip().lower() for p in properties]

        try:
            for item in molecules:
                if isinstance(item, Molecule):
                    mol = item
                elif isinstance(item, SpeciesSummary):
                    mol = Molecule(summary=item, client=self)
                else:
                    mol = self.get_molecule(str(item))

                props_out: Dict[str, object] = {}
                urls: List[str] = []
                errs: List[str] = []

                try:
                    if method.strip().lower() in ("experimental", "exp"):
                        exp = mol.get_experimental()
                        if exp.provenance:
                            urls.append(exp.provenance.source_url)
                        if "dipole" in req_props:
                            props_out["dipole_debye"] = exp.dipole.total if exp.dipole else None
                        if "zpe" in req_props:
                            props_out["zpe_cm_inv"] = (
                                exp.vibrations.zpe_unscaled_cm_inv if exp.vibrations else None
                            )
                        if "enthalpy" in req_props or "hfg" in req_props:
                            h298 = next(
                                (t.value for t in exp.thermochemistry if "298" in t.property_name and "hfg" in t.property_name.lower()),
                                None,
                            )
                            props_out["hfg_298k_kj_mol"] = h298
                    else:
                        # Try single calculation dossier first
                        try:
                            calc = mol.get_calculation(method=method, basis_set=basis_set)
                            if calc.provenance:
                                urls.append(calc.provenance.source_url)
                            if "energy" in req_props:
                                props_out["energy_hartree"] = calc.electronic_energy
                            if "dipole" in req_props:
                                props_out["dipole_debye"] = calc.dipole.total if calc.dipole else None
                            if "zpe" in req_props:
                                props_out["zpe_cm_inv"] = (
                                    calc.vibrations.zpe_unscaled_cm_inv if calc.vibrations else None
                                )
                                props_out["zpe_scaled_cm_inv"] = (
                                    calc.vibrations.zpe_scaled_cm_inv if calc.vibrations else None
                                )
                        except CCCBDBError:
                            if "energy" in req_props:
                                e_col = mol.get_energies(method=method, basis_set=basis_set)
                                props_out["energy_hartree"] = e_col[0].value if len(e_col) > 0 else None
                                if len(e_col) > 0 and e_col[0].provenance:
                                    urls.append(e_col[0].provenance.source_url)
                except Exception as exc:
                    err_msg = str(exc)
                    errs.append(err_msg)
                    partial_failures[mol.id] = err_msg

                dataset_records.append(
                    DatasetRecord(
                        species_id=mol.id,
                        cas_rn=mol.cas_rn,
                        name=mol.name,
                        formula=mol.formula,
                        charge=mol.charge,
                        point_group=mol.point_group,
                        molecular_weight=mol.molecular_weight,
                        method=method,
                        basis_set=basis_set,
                        properties=props_out,
                        source_urls=urls,
                        errors=errs,
                    )
                )
        finally:
            self.transport.rate_limiter.delay_seconds = old_delay

        return BenchmarkDataset(
            records=dataset_records,
            requested_properties=list(properties),
            method=method,
            basis_set=basis_set,
            partial_failures=partial_failures,
            created_at=datetime.now(timezone.utc).isoformat(),
        )

    @staticmethod
    def available_methods() -> List[ComputationalMethodSpec]:
        """List all quantum chemical and empirical methods in the NIST CCCBDB ontology."""
        return list_methods()

    @staticmethod
    def available_basis_sets() -> List[BasisSetSpec]:
        """List all Gaussian basis sets and ECPs in the NIST CCCBDB ontology."""
        return list_basis_sets()

    def close(self) -> None:
        """Close network sessions and SQLite connections."""
        self.transport.close()
        if self.sqlite:
            self.sqlite.close()

    def __enter__(self) -> "CCCBDB":
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        self.close()
