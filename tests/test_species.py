"""Unit tests for species discovery, isomer handling, and experimental reference parsing."""

import unittest
from pathlib import Path

from cccbdb import CCCBDB, AmbiguousSpeciesError, SpeciesNotFoundError
from cccbdb.parsers.species import parse_experimental_page

FIXTURES = Path(__file__).parent / "fixtures"


class TestSpeciesDiscoveryAndParsing(unittest.TestCase):
    def test_parse_experimental_co_fixture(self) -> None:
        html = (FIXTURES / "exp2x_co.html").read_text(encoding="utf-8")
        exp = parse_experimental_page(html, source_url="https://cccbdb.nist.gov/exp2x.asp?casno=630080&charge=0")
        self.assertEqual(exp.species_id, "630080")
        self.assertEqual(exp.cas_rn, "630-08-0")
        self.assertEqual(exp.formula, "CO")
        self.assertEqual(exp.point_group, "Cinfv")
        self.assertEqual(exp.inchikey, "UGFAIRIUMAVXCW-UHFFFAOYSA-N")
        self.assertEqual(len(exp.thermochemistry), 5)
        self.assertAlmostEqual(exp.thermochemistry[0].value, -110.53)
        self.assertIsNotNone(exp.geometry)
        self.assertEqual(exp.geometry.atoms, ["C", "O"])
        self.assertAlmostEqual(exp.geometry.distance(0, 1), 1.1282, places=4)
        self.assertIsNotNone(exp.vibrations)
        self.assertEqual(len(exp.vibrations.modes), 1)
        self.assertAlmostEqual(exp.vibrations.modes[0].fundamental_frequency_cm_inv, 2143.0)
        self.assertAlmostEqual(exp.vibrations.zpe_unscaled_cm_inv, 1071.6)
        self.assertIsNotNone(exp.dipole)
        self.assertAlmostEqual(exp.dipole.total, 0.112)
        self.assertIsNotNone(exp.polarizability)
        self.assertAlmostEqual(exp.polarizability.isotropic, 1.953)
        self.assertIn("CODATA", exp.references)

    def test_search_and_isomer_disambiguation(self) -> None:
        with CCCBDB(offline=True) as db:
            # Search H2O
            water_list = db.search(formula="H2O")
            self.assertEqual(len(water_list), 1)
            self.assertEqual(water_list[0].id, "7732185")
            self.assertEqual(water_list[0].name, "Water")

            # Search C2H6O isomers (Ethanol and Dimethyl ether)
            isomers = db.search(formula="C2H6O")
            self.assertEqual(len(isomers), 2)
            names = {m.name for m in isomers}
            self.assertEqual(names, {"Ethanol", "Dimethyl ether"})

            # Ambiguous get_molecule("C2H6O") must raise AmbiguousSpeciesError
            with self.assertRaises(AmbiguousSpeciesError):
                db.get_molecule("C2H6O")

            # Specific CAS RN or structural formula succeeds unambiguously
            ethanol = db.get_molecule("64-17-5")
            self.assertEqual(ethanol.name, "Ethanol")

            # Non-existent species raises SpeciesNotFoundError
            with self.assertRaises(SpeciesNotFoundError):
                db.get_molecule("C99H198O99")


if __name__ == "__main__":
    unittest.main()
