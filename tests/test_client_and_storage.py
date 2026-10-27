"""End-to-end integration tests for CCCBDB client, DiskCache, SQLiteStorage, and Dataset export."""

import tempfile
import unittest
from pathlib import Path

from cccbdb import CCCBDB
from cccbdb.cli import main as cli_main

FIXTURES = Path(__file__).parent / "fixtures"


class TestClientStorageAndCLI(unittest.TestCase):
    def test_client_with_populated_cache_and_sqlite(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir) / "cache"
            with CCCBDB(cache=True, cache_dir=cache_dir, offline=True, use_sqlite=True) as db:
                # Pre-populate cache with fixtures for H2O and CO
                url_exp_h2o = db.transport.build_url("exp2x.asp", {"casno": "7732185"})
                url_e2_h2o = db.transport.build_url(
                    "energy2x.asp", {"_ctx_casno": "7732185", "_ctx_charge": 0}
                )
                url_e3_h2o = db.transport.build_url(
                    "energy3x.asp", {"basis": 1, "casno": "7732185", "method": 8}
                )
                url_exp_co = db.transport.build_url("exp2x.asp", {"casno": "630080"})

                db.disk_cache.set(url_exp_h2o, (FIXTURES / "exp2x_h2o.html").read_text(encoding="utf-8"))
                db.disk_cache.set(url_e2_h2o, (FIXTURES / "energy2x_h2o.html").read_text(encoding="utf-8"))
                db.disk_cache.set(
                    url_e3_h2o, (FIXTURES / "energy3x_h2o_b3lyp_631gs.html").read_text(encoding="utf-8")
                )
                db.disk_cache.set(url_exp_co, (FIXTURES / "exp2x_co.html").read_text(encoding="utf-8"))

                # Example A: Search for molecules
                molecules = db.search(formula="H2O")
                self.assertEqual(len(molecules), 1)
                self.assertEqual(molecules[0].id, "7732185")

                # Example B: Retrieve electronic energies
                water = db.get_molecule("H2O")
                energies = water.get_energies(method="B3LYP", basis_set="6-31G*")
                df = energies.to_dataframe()
                self.assertGreaterEqual(len(df), 1)
                self.assertAlmostEqual(energies.filter(energy_type="electronic_0K")[0].value, -76.408953)

                # Example C: Retrieve optimized molecular geometry
                geometry = water.get_geometry(method="B3LYP", basis_set="6-31G*")
                self.assertEqual(geometry.atoms, ["O", "H", "H"])
                xyz_out = Path(tmpdir) / "water.xyz"
                geometry.to_xyz(xyz_out)
                self.assertTrue(xyz_out.exists())

                # Verify SQLite persistence
                stored = db.sqlite.query_energies(species_id="7732185", method="B3LYP", basis_set="6-31G*")
                self.assertGreaterEqual(len(stored), 1)

                # Example D: Collect and export dataset
                dataset = db.collect(
                    molecules=[water, "CO"],
                    properties=["dipole", "zpe", "enthalpy"],
                    method="Experimental",
                    basis_set=None,
                    delay=0.0,
                )
                csv_path = Path(tmpdir) / "benchmark.csv"
                dataset.to_csv(csv_path)
                self.assertTrue(csv_path.exists())
                self.assertTrue(Path(str(csv_path) + ".manifest.json").exists())

    def test_cli_commands(self) -> None:
        self.assertEqual(cli_main(["search", "--formula", "C2H6O"]), 0)
        self.assertEqual(cli_main(["methods", "--category", "Coupled"]), 0)
        self.assertEqual(cli_main(["basis-sets", "--family", "Dunning"]), 0)
        self.assertEqual(cli_main(["symmetry", "C3v"]), 0)
        self.assertEqual(cli_main(["convert-units", "1.0", "Hartree", "kcal/mol"]), 0)


if __name__ == "__main__":
    unittest.main()
