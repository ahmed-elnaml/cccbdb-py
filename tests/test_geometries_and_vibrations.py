"""Unit tests for geometries, XYZ export, dummy-atom filtering, and imaginary vibrational modes."""

import tempfile
import unittest
from pathlib import Path

from cccbdb.chemistry.units import EnergyUnit, LengthUnit
from cccbdb.exceptions import IncompleteGeometryError
from cccbdb.models.geometry import InternalCoordinate, MolecularGeometry
from cccbdb.parsers.energies import parse_calculation_result_page

FIXTURES = Path(__file__).parent / "fixtures"


class TestGeometriesAndVibrations(unittest.TestCase):
    def test_geometry_xyz_export_and_unit_conversion(self) -> None:
        html = (FIXTURES / "energy3x_h2o_b3lyp_631gs.html").read_text(encoding="utf-8")
        calc = parse_calculation_result_page(html, species_id="7732185")
        geom = calc.geometry
        self.assertIsNotNone(geom)
        self.assertEqual(geom.atoms, ["O", "H", "H"])
        self.assertEqual(geom.point_group, "C2v")
        self.assertEqual(geom.symmetry_info.rotational_symmetry_number, 2)

        # Convert Angstrom -> Bohr
        geom_bohr = geom.convert_units(LengthUnit.BOHR)
        self.assertEqual(geom_bohr.unit, "Bohr")
        self.assertAlmostEqual(geom_bohr.sites[0].z, 0.1197 / 0.529177210903, places=5)

        # Export to XYZ (always standardizes to Angstroms)
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "water.xyz"
            geom_bohr.to_xyz(xyz_path)
            lines = xyz_path.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(lines[0], "3")
            self.assertIn("B3LYP/6-31G*", lines[1])
            self.assertTrue(lines[2].startswith("O "))

    def test_incomplete_geometry_rejects_xyz_export(self) -> None:
        # Internal coordinates only (no 3D Cartesians) must raise IncompleteGeometryError
        geom_internal_only = MolecularGeometry(
            species_id="7732185",
            formula="H2O",
            method="Experimental",
            sites=[],
            internal_coordinates=[
                InternalCoordinate(name="rOH", coord_type="distance", value=0.9578, unit="Angstrom"),
                InternalCoordinate(name="aHOH", coord_type="angle", value=104.48, unit="degrees"),
            ],
        )
        self.assertFalse(geom_internal_only.has_complete_cartesian)
        with self.assertRaises(IncompleteGeometryError):
            geom_internal_only.to_xyz_string()

    def test_multi_conformer_dummy_atoms_and_imaginary_frequencies(self) -> None:
        html = (FIXTURES / "energy3x_ch2_multi_conformer.html").read_text(encoding="utf-8")
        calc = parse_calculation_result_page(html, species_id="2465567")
        self.assertEqual(len(calc.conformers), 2)

        # Conformer 1: C2v minimum, has dummy atom X4 in raw sites which must be excluded from physical_sites
        c1 = calc.conformers[0]
        self.assertTrue(c1.is_minimum)
        self.assertEqual(c1.point_group, "C2v")
        self.assertEqual(len(c1.geometry.sites), 4)
        self.assertEqual(len(c1.geometry.physical_sites), 3)
        self.assertEqual(c1.geometry.atoms, ["C", "H", "H"])
        self.assertFalse(c1.vibrations.is_transition_state)

        # Conformer 2: Dinfh saddle point with two degenerate 872i imaginary modes
        c2 = calc.conformers[1]
        self.assertFalse(c2.is_minimum)
        self.assertEqual(c2.point_group, "Dinfh")
        self.assertTrue(c2.vibrations.is_transition_state)
        self.assertEqual(len(c2.vibrations.imaginary_modes), 2)
        self.assertAlmostEqual(c2.vibrations.imaginary_modes[0].frequency_cm_inv, -872.0)
        self.assertAlmostEqual(c2.vibrations.zpe_in(EnergyUnit.KJ_MOL, scaled=True), 28.978, places=2)


if __name__ == "__main__":
    unittest.main()
