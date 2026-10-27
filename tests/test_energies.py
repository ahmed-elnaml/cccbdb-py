"""Unit tests for electronic and thermochemical energy extraction, filtering, and DataFrame conversion."""

import unittest
from pathlib import Path

from cccbdb.chemistry.units import EnergyUnit
from cccbdb.parsers.energies import parse_calculation_result_page, parse_energies_matrix_page

FIXTURES = Path(__file__).parent / "fixtures"


class TestEnergiesParsing(unittest.TestCase):
    def test_parse_energies_matrix_page(self) -> None:
        html = (FIXTURES / "energy2x_h2o.html").read_text(encoding="utf-8")
        collection = parse_energies_matrix_page(
            html,
            species_id="7732185",
            formula="H2O",
            source_url="https://cccbdb.nist.gov/energy2x.asp",
        )
        # 4 predefined + 7 HF + 7 B3LYP + 6 MP2 (STO-3G is '-') + 5 CCSD(T) = 29 records
        self.assertEqual(len(collection), 29)

        # Verify missing '-' cell for MP2/STO-3G was NOT coerced to 0.0
        mp2_sto3g = collection.filter(method="MP2", basis_set="STO-3G")
        self.assertEqual(len(mp2_sto3g), 0)

        # Filter B3LYP / 6-31G*
        b3lyp_631gs = collection.filter(method="B3LYP", basis_set="6-31G*")
        self.assertEqual(len(b3lyp_631gs), 1)
        rec = b3lyp_631gs[0]
        self.assertAlmostEqual(rec.value, -76.408953, places=6)
        self.assertEqual(rec.unit, "Hartree")
        self.assertIsNotNone(rec.provenance)
        self.assertEqual(rec.provenance.original_unit, "Hartree")

        # Convert to kJ/mol and DataFrame
        in_kj = b3lyp_631gs.convert_units(EnergyUnit.KJ_MOL)
        self.assertEqual(in_kj[0].unit, "kJ/mol")
        self.assertAlmostEqual(in_kj[0].value, -76.408953 * 2625.49963948, places=2)

        df = collection.to_dataframe()
        self.assertEqual(len(df), 29)
        self.assertIn("method", df.columns)
        self.assertIn("basis_set", df.columns)

    def test_parse_single_calculation_page(self) -> None:
        html = (FIXTURES / "energy3x_h2o_b3lyp_631gs.html").read_text(encoding="utf-8")
        calc = parse_calculation_result_page(
            html,
            source_url="https://cccbdb.nist.gov/energy3x.asp?method=8&basis=1&casno=7732185",
            species_id="7732185",
        )
        self.assertEqual(calc.method, "B3LYP")
        self.assertEqual(calc.basis_set, "6-31G*")
        self.assertAlmostEqual(calc.electronic_energy, -76.408953, places=6)
        self.assertIsNotNone(calc.primary)
        self.assertAlmostEqual(calc.primary.energies["energy_298K"].value, -76.406118, places=6)
        self.assertAlmostEqual(calc.primary.energies["nuclear_repulsion"].value, 9.157175, places=6)
        self.assertEqual(len(calc.primary.atomic_charges), 3)
        self.assertAlmostEqual(calc.primary.atomic_charges[0].mulliken, -0.788)
        self.assertIsNotNone(calc.dipole)
        self.assertAlmostEqual(calc.dipole.total, 2.099)
        self.assertIsNotNone(calc.primary.polarizability)
        self.assertAlmostEqual(calc.primary.polarizability.isotropic, (0.812 + 1.415 + 1.083) / 3.0, places=5)
        self.assertAlmostEqual(calc.primary.r2_angstrom2, 19.124)


if __name__ == "__main__":
    unittest.main()
