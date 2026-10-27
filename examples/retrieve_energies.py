"""Example 2: Retrieving electronic energies, geometries, and vibrations using local fixtures."""

from pathlib import Path
import tempfile

from cccbdb import CCCBDB, EnergyUnit, LengthUnit

FIXTURES = Path(__file__).resolve().parent.parent / "tests" / "fixtures"


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        with CCCBDB(cache=True, cache_dir=tmpdir, offline=True) as db:
            # Seed cache from repository test fixtures so example runs 100% offline
            db.disk_cache.set(
                db.transport.build_url("energy2x.asp", {"_ctx_casno": "7732185", "_ctx_charge": 0}),
                (FIXTURES / "energy2x_h2o.html").read_text(encoding="utf-8"),
            )
            db.disk_cache.set(
                db.transport.build_url("energy3x.asp", {"basis": 1, "casno": "7732185", "method": 8}),
                (FIXTURES / "energy3x_h2o_b3lyp_631gs.html").read_text(encoding="utf-8"),
            )

            water = db.get_molecule("H2O")

            print("=== 1. B3LYP/6-31G* Energies (Hartree & kJ/mol) ===")
            energies = water.get_energies(method="B3LYP", basis_set="6-31G*")
            for rec in energies:
                rec_kj = rec.to_unit(EnergyUnit.KJ_MOL)
                print(f"  {rec.energy_type:<18}: {rec.value:14.6f} {rec.unit}  ({rec_kj.value:14.2f} {rec_kj.unit})")

            print("\n=== 2. Optimized Molecular Geometry (Angstrom & Bohr) ===")
            geom = water.get_geometry(method="B3LYP", basis_set="6-31G*")
            print(geom.to_xyz_string())
            geom_bohr = geom.convert_units(LengthUnit.BOHR)
            print("Center of mass (Angstrom):", geom.center_of_mass())
            print("O-H1 bond distance (Angstrom):", round(geom.distance(0, 1), 5))
            print("O-H1 bond distance (Bohr):", round(geom_bohr.distance(0, 1), 5))

            print("\n=== 3. Harmonic Vibrational Frequencies & ZPE ===")
            vibs = water.get_vibrations(method="B3LYP", basis_set="6-31G*")
            for m in vibs.modes:
                print(
                    f"  Mode {m.mode_number} ({m.symmetry}): "
                    f"unscaled={m.frequency_cm_inv:.1f} cm^-1 | "
                    f"scaled={m.scaled_frequency_cm_inv:.1f} cm^-1 | "
                    f"IR={m.ir_intensity_km_mol} km/mol"
                )
            print(
                f"  Scaled ZPE: {vibs.zpe_scaled_cm_inv} cm^-1 "
                f"({vibs.zpe_in(EnergyUnit.KCAL_MOL):.3f} kcal/mol)"
            )


if __name__ == "__main__":
    main()
