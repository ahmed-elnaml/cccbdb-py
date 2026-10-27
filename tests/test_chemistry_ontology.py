"""Unit tests for physical units, point group symmetry, methods, basis sets, and formulas."""

import unittest

from cccbdb.chemistry import (
    BasisSetFamily,
    DipoleUnit,
    EnergyUnit,
    LengthUnit,
    MethodCategory,
    RotationalConstantUnit,
    calculate_molecular_weight,
    convert_dipole,
    convert_energy,
    convert_length,
    convert_rotational_constant,
    convert_value,
    format_cas_hyphenated,
    get_vibrational_scaling_factor,
    list_basis_sets,
    list_methods,
    normalize_cas,
    normalize_point_group_symbol,
    parse_formula,
    parse_point_group,
    resolve_basis_set,
    resolve_method,
    to_hill_formula,
    validate_cas_checksum,
)


class TestChemistryOntology(unittest.TestCase):
    def test_energy_unit_conversions(self) -> None:
        # 1 Hartree = 2625.4996394799 kJ/mol = 627.5094740631 kcal/mol = 27.211386245988 eV
        self.assertAlmostEqual(convert_energy(1.0, EnergyUnit.HARTREE, EnergyUnit.KJ_MOL), 2625.499639, places=5)
        self.assertAlmostEqual(convert_energy(1.0, "Hartree", "kcal/mol"), 627.509474, places=5)
        self.assertAlmostEqual(convert_energy(1.0, "eh", "eV"), 27.211386, places=5)
        self.assertAlmostEqual(convert_energy(1.0, "Hartree", "cm^-1"), 219474.63136, places=4)
        # Round-trip conversion
        val = -76.408953
        in_kj = convert_energy(val, "Hartree", "kJ/mol")
        back = convert_energy(in_kj, "kJ/mol", "Hartree")
        self.assertAlmostEqual(val, back, places=8)

    def test_length_dipole_rotational_conversions(self) -> None:
        self.assertAlmostEqual(convert_length(1.0, LengthUnit.BOHR, LengthUnit.ANGSTROM), 0.52917721, places=7)
        self.assertAlmostEqual(convert_length(1.0, "Angstrom", "pm"), 100.0, places=6)
        self.assertAlmostEqual(convert_dipole(1.0, DipoleUnit.ATOMIC_UNIT, DipoleUnit.DEBYE), 2.54174647, places=7)
        self.assertAlmostEqual(convert_rotational_constant(1.0, RotationalConstantUnit.CM_INV, "GHz"), 29.9792458, places=6)
        self.assertAlmostEqual(convert_value(1.0, "Bohr", "Angstrom"), 0.52917721, places=7)

    def test_point_group_symmetry_normalization_and_properties(self) -> None:
        self.assertEqual(normalize_point_group_symbol("C*V"), "Cinfv")
        self.assertEqual(normalize_point_group_symbol("D*H"), "Dinfh")
        self.assertEqual(normalize_point_group_symbol("D∞h"), "Dinfh")
        self.assertEqual(normalize_point_group_symbol("C2V"), "C2v")

        c2v = parse_point_group("C2V")
        self.assertEqual(c2v.symbol, "C2v")
        self.assertEqual(c2v.rotational_symmetry_number, 2)
        self.assertTrue(c2v.allows_permanent_dipole)
        self.assertFalse(c2v.has_inversion_center)
        self.assertFalse(c2v.is_chiral)

        d6h = parse_point_group("D6h")
        self.assertEqual(d6h.rotational_symmetry_number, 12)
        self.assertTrue(d6h.has_inversion_center)
        self.assertFalse(d6h.allows_permanent_dipole)

        td = parse_point_group("Td")
        self.assertEqual(td.rotational_symmetry_number, 12)
        self.assertFalse(td.allows_permanent_dipole)

        dinfh = parse_point_group("D*H")
        self.assertTrue(dinfh.is_linear)
        self.assertEqual(dinfh.rotational_symmetry_number, 2)

    def test_methods_and_basis_sets_registries(self) -> None:
        b3lyp = resolve_method("b3lyp")
        self.assertEqual(b3lyp.name, "B3LYP")
        self.assertEqual(b3lyp.nist_id, 8)
        self.assertEqual(b3lyp.category, MethodCategory.DFT)

        mp2_full = resolve_method("MP2(full)")
        self.assertEqual(mp2_full.name, "MP2=FULL")
        self.assertEqual(mp2_full.core_correlation, "full")

        g4 = resolve_method("G4")
        self.assertTrue(g4.predefined_basis)

        basis = resolve_basis_set("6-31G(d)")
        self.assertEqual(basis.name, "6-31G*")
        self.assertEqual(basis.nist_id, 1)
        self.assertEqual(basis.family, BasisSetFamily.POPLE)

        sf = get_vibrational_scaling_factor("B3LYP", "6-31G*")
        self.assertAlmostEqual(sf, 0.9603, places=4)

        self.assertGreaterEqual(len(list_methods()), 40)
        self.assertGreaterEqual(len(list_basis_sets()), 40)

    def test_formula_parser_and_cas_checksum(self) -> None:
        comp, charge = parse_formula("CH3(CH2)4CH3")
        self.assertEqual(comp, {"C": 6, "H": 14})
        self.assertEqual(charge, 0)
        self.assertEqual(to_hill_formula("CH3(CH2)4CH3"), "C6H14")

        comp_ion, charge_ion = parse_formula("CO3--")
        self.assertEqual(comp_ion, {"C": 1, "O": 3})
        self.assertEqual(charge_ion, -2)

        self.assertAlmostEqual(calculate_molecular_weight("H2O"), 18.01528, places=4)
        self.assertEqual(normalize_cas("7732-18-5"), "7732185")
        self.assertEqual(format_cas_hyphenated("7732185"), "7732-18-5")
        self.assertTrue(validate_cas_checksum("7732-18-5"))
        self.assertTrue(validate_cas_checksum("630-08-0"))
        self.assertFalse(validate_cas_checksum("7732-18-9"))


if __name__ == "__main__":
    unittest.main()
