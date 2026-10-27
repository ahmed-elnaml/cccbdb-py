"""Command-line interface (`cccbdb`) for researchers and terminal workflows."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

from cccbdb.chemistry.basis_sets import list_basis_sets
from cccbdb.chemistry.methods import list_methods
from cccbdb.chemistry.symmetry import parse_point_group
from cccbdb.chemistry.units import convert_value
from cccbdb.client import CCCBDB
from cccbdb.config import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cccbdb",
        description="cccbdb-py: Open-source Python toolkit for programmatic access to NIST CCCBDB.",
    )
    parser.add_argument("--version", action="version", version=f"cccbdb-py {__version__}")
    parser.add_argument("--cache-dir", type=str, default=None, help="Custom local cache directory")
    parser.add_argument("--offline", action="store_true", help="Operate strictly from local cache")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. cccbdb search
    p_search = subparsers.add_parser("search", help="Search for chemical species in CCCBDB")
    p_search.add_argument("--formula", type=str, help="Chemical formula (e.g. H2O, C2H6O, C6H6)")
    p_search.add_argument("--name", type=str, help="Species name substring (e.g. Water, Benzene)")
    p_search.add_argument("--cas", type=str, help="CAS Registry Number (e.g. 7732-18-5)")
    p_search.add_argument("--point-group", type=str, help="Schoenflies point group (e.g. C2v, D6h, Td)")
    p_search.add_argument("--elements", type=str, help="Comma-separated allowed elements (e.g. C,H,O)")
    p_search.add_argument("--json", action="store_true", help="Output results as JSON")

    # 2. cccbdb get
    p_get = subparsers.add_parser("get", help="Retrieve molecular properties for a species")
    p_get.add_argument("species", type=str, help="CAS RN, unambiguous formula, or species name")
    p_get.add_argument(
        "--property",
        choices=["summary", "energies", "geometry", "vibrations", "dipole", "experimental"],
        default="summary",
        help="Property type to retrieve",
    )
    p_get.add_argument("--method", type=str, default=None, help="Computational method (e.g. B3LYP, CCSD(T))")
    p_get.add_argument("--basis", type=str, default=None, help="Basis set (e.g. 6-31G*, cc-pVTZ)")
    p_get.add_argument("--unit", type=str, default=None, help="Target physical unit for conversion")
    p_get.add_argument("--xyz", type=str, default=None, help="Export geometry to .xyz file path")

    # 3. cccbdb export
    p_export = subparsers.add_parser("export", help="Collect and export a benchmark dataset")
    p_export.add_argument("--elements", type=str, default="C,H,O", help="Comma-separated element filter")
    p_export.add_argument("--method", type=str, default="B3LYP", help="Computational method")
    p_export.add_argument("--basis", type=str, default="6-31G*", help="Basis set")
    p_export.add_argument("--properties", type=str, default="energy,dipole", help="Comma-separated properties")
    p_export.add_argument("--output", "-o", type=str, required=True, help="Output file path (.csv or .json)")

    # 4. cccbdb methods
    p_methods = subparsers.add_parser("methods", help="List supported NIST CCCBDB quantum chemistry methods")
    p_methods.add_argument("--category", type=str, default=None, help="Filter by category substring")

    # 5. cccbdb basis-sets
    p_basis = subparsers.add_parser("basis-sets", help="List supported NIST CCCBDB Gaussian basis sets")
    p_basis.add_argument("--family", type=str, default=None, help="Filter by family (Pople, Dunning, Ahlrichs, ECP)")

    # 6. cccbdb symmetry
    p_sym = subparsers.add_parser("symmetry", help="Inspect a Schoenflies molecular point group")
    p_sym.add_argument("symbol", type=str, help="Point group symbol (e.g. C2v, D3h, Td, Dinfh)")

    # 7. cccbdb convert-units
    p_conv = subparsers.add_parser("convert-units", help="Convert physical chemistry quantities between units")
    p_conv.add_argument("value", type=float, help="Numerical value")
    p_conv.add_argument("from_unit", type=str, help="Source unit (e.g. Hartree, kJ/mol, Bohr, Debye)")
    p_conv.add_argument("to_unit", type=str, help="Target unit (e.g. kcal/mol, eV, Angstrom, au)")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "methods":
        methods = list_methods()
        if args.category:
            methods = [m for m in methods if args.category.lower() in m.category.value.lower()]
        print(f"{'Method':<22} {'NIST ID':<8} {'Category':<32} {'Core':<8}")
        print("-" * 74)
        for m in methods:
            print(f"{m.name:<22} {str(m.nist_id):<8} {m.category.value:<32} {m.core_correlation:<8}")
        return 0

    if args.command == "basis-sets":
        bases = list_basis_sets()
        if args.family:
            bases = [b for b in bases if args.family.lower() in b.family.value.lower()]
        print(f"{'Basis Set':<20} {'NIST ID':<8} {'Family':<12} {'Zeta':<8} {'Size':<8}")
        print("-" * 60)
        for b in bases:
            print(f"{b.name:<20} {b.nist_id:<8} {b.family.value:<12} {b.zeta_quality:<8} {b.size_class:<8}")
        return 0

    if args.command == "symmetry":
        info = parse_point_group(args.symbol)
        print(json.dumps(info.__dict__, indent=2, ensure_ascii=False))
        return 0

    if args.command == "convert-units":
        res = convert_value(args.value, args.from_unit, args.to_unit)
        print(f"{args.value} {args.from_unit} = {res:.8g} {args.to_unit}")
        return 0

    with CCCBDB(cache_dir=args.cache_dir, offline=args.offline) as db:
        if args.command == "search":
            elems = [e.strip() for e in args.elements.split(",")] if args.elements else None
            mols = db.search(
                formula=args.formula,
                name=args.name,
                cas=args.cas,
                point_group=args.point_group,
                elements=elems,
            )
            if args.json:
                print(json.dumps([m.to_dict() for m in mols], indent=2))
            else:
                print(f"{'ID (CAS)':<12} {'Formula':<12} {'Name':<24} {'Point Group':<12} {'MW (g/mol)':<10}")
                print("-" * 74)
                for m in mols:
                    print(
                        f"{m.id:<12} {m.formula:<12} {m.name:<24} "
                        f"{(m.point_group or '-'):<12} {m.molecular_weight:<10.4f}"
                    )
            return 0

        if args.command == "get":
            mol = db.get_molecule(args.species)
            if args.property == "summary":
                print(json.dumps(mol.to_dict(), indent=2))
            elif args.property == "experimental":
                exp = mol.get_experimental()
                print(exp.model_dump_json(indent=2))
            elif args.property == "energies":
                unit = args.unit or "Hartree"
                energies = mol.get_energies(method=args.method, basis_set=args.basis, unit=unit)
                print(json.dumps(energies.to_dicts(), indent=2))
            elif args.property == "geometry":
                method = args.method or "Experimental"
                unit = args.unit or "Angstrom"
                geom = mol.get_geometry(method=method, basis_set=args.basis, unit=unit)
                if args.xyz:
                    out = geom.to_xyz(args.xyz)
                    print(f"Wrote XYZ coordinates to {out}")
                else:
                    print(geom.to_xyz_string(), end="")
            elif args.property == "vibrations":
                method = args.method or "Experimental"
                vibs = mol.get_vibrations(method=method, basis_set=args.basis)
                print(vibs.model_dump_json(indent=2))
            elif args.property == "dipole":
                method = args.method or "Experimental"
                unit = args.unit or "Debye"
                dip = mol.get_dipole(method=method, basis_set=args.basis, unit=unit)
                print(dip.model_dump_json(indent=2) if dip else "null")
            return 0

        if args.command == "export":
            elems = [e.strip() for e in args.elements.split(",")]
            props = [p.strip() for p in args.properties.split(",")]
            mols = db.search(elements=elems)
            dataset = db.collect(
                molecules=mols,
                properties=props,
                method=args.method,
                basis_set=args.basis,
            )
            out_path = Path(args.output)
            if out_path.suffix.lower() == ".json":
                dataset.to_json(out_path)
            else:
                dataset.to_csv(out_path)
            print(f"Exported {len(dataset.records)} molecules to {out_path}")
            return 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
