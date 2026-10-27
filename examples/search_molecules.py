"""Example 1: Searching molecules, inspecting symmetry point groups, and resolving isomers."""

from cccbdb import CCCBDB, AmbiguousSpeciesError, parse_point_group


def main() -> None:
    with CCCBDB(offline=True) as db:
        print("=== 1. Search by Formula (H2O) ===")
        for mol in db.search(formula="H2O"):
            sym = mol.symmetry
            print(
                f"ID: {mol.id} | CAS: {mol.cas_rn} | Name: {mol.name} | "
                f"Formula: {mol.formula} | Point Group: {mol.point_group} "
                f"(sigma={sym.rotational_symmetry_number if sym else 'N/A'})"
            )

        print("\n=== 2. Search Constitutional Isomers (C2H6O) ===")
        for mol in db.search(formula="C2H6O"):
            print(f"ID: {mol.id} | Formula: {mol.formula:<10} | Name: {mol.name:<16} | SMILES: {mol.smiles}")

        print("\n=== 3. Explicit Isomer Disambiguation Safeguard ===")
        try:
            db.get_molecule("C2H6O")
        except AmbiguousSpeciesError as exc:
            print(f"Caught expected AmbiguousSpeciesError:\n  {exc}")

        # Resolve specific isomer by CAS RN
        ethanol = db.get_molecule("64-17-5")
        print(f"Resolved unambiguously by CAS RN: {ethanol}")

        print("\n=== 4. Point Group Symmetry Analysis ===")
        for pg in ["C2v", "C3v", "D6h", "Td", "Dinfh"]:
            info = parse_point_group(pg)
            print(
                f"Point Group {info.unicode_symbol:<5} ({info.family:<28}): "
                f"sigma={info.rotational_symmetry_number:<2} | "
                f"dipole_allowed={str(info.allows_permanent_dipole):<5} | "
                f"inversion={info.has_inversion_center}"
            )


if __name__ == "__main__":
    main()
