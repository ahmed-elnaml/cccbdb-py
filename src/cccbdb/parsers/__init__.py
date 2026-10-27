"""Isolated HTML parsers for NIST CCCBDB pages."""

from cccbdb.parsers.energies import parse_calculation_result_page, parse_energies_matrix_page
from cccbdb.parsers.frequencies import parse_calculated_vibrations_table, parse_experimental_vibrations
from cccbdb.parsers.geometries import parse_calculated_geometry_section, parse_experimental_geometry
from cccbdb.parsers.properties import (
    parse_atomic_charges_table,
    parse_calculated_dipole_table,
    parse_calculated_polarizability_table,
    parse_calculated_quadrupole_tables,
    parse_calculated_rotational_constants,
    parse_experimental_dipole,
    parse_experimental_polarizability,
    parse_experimental_quadrupoles,
    parse_experimental_rotational_constants,
    parse_experimental_thermochemistry,
)
from cccbdb.parsers.species import (
    BUILTIN_SPECIES_CATALOG,
    parse_experimental_page,
    parse_species_disambiguation_page,
)

__all__ = [
    "BUILTIN_SPECIES_CATALOG",
    "parse_calculation_result_page",
    "parse_energies_matrix_page",
    "parse_calculated_vibrations_table",
    "parse_experimental_vibrations",
    "parse_calculated_geometry_section",
    "parse_experimental_geometry",
    "parse_atomic_charges_table",
    "parse_calculated_dipole_table",
    "parse_calculated_polarizability_table",
    "parse_calculated_quadrupole_tables",
    "parse_calculated_rotational_constants",
    "parse_experimental_dipole",
    "parse_experimental_polarizability",
    "parse_experimental_quadrupoles",
    "parse_experimental_rotational_constants",
    "parse_experimental_thermochemistry",
    "parse_experimental_page",
    "parse_species_disambiguation_page",
]
