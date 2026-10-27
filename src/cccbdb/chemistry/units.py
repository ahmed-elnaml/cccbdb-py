"""Physical constants and unit conversion engine for computational chemistry.

All conversion factors are grounded in CODATA 2018 recommended values and
NIST CCCBDB standard unit definitions (Section VII.A of NIST SRD 101).
"""

from __future__ import annotations

from enum import Enum
from typing import Dict, Tuple, Union


class EnergyUnit(str, Enum):
    """Supported energy and thermochemical units."""

    HARTREE = "Hartree"
    EV = "eV"
    KJ_MOL = "kJ/mol"
    KCAL_MOL = "kcal/mol"
    CM_INV = "cm^-1"
    KELVIN = "K"
    JOULE = "J"
    MHZ = "MHz"


class LengthUnit(str, Enum):
    """Supported atomic and molecular length units."""

    ANGSTROM = "Angstrom"
    BOHR = "Bohr"
    PICOMETER = "pm"
    NANOMETER = "nm"


class DipoleUnit(str, Enum):
    """Supported electric dipole moment units."""

    DEBYE = "Debye"
    COULOMB_METER = "C*m"
    ATOMIC_UNIT = "au"  # e * a0


class RotationalConstantUnit(str, Enum):
    """Supported rotational constant units."""

    CM_INV = "cm^-1"
    GHZ = "GHz"
    MHZ = "MHz"


class MomentOfInertiaUnit(str, Enum):
    """Supported moment of inertia units."""

    AMU_ANGSTROM2 = "amu*Angstrom^2"
    G_CM2 = "g*cm^2"
    ATOMIC_UNIT = "au"  # me * a0^2


class PolarizabilityUnit(str, Enum):
    """Supported electric dipole polarizability volume units."""

    ANGSTROM3 = "Angstrom^3"
    BOHR3 = "Bohr^3"


# CODATA 2018 Exact / Recommended Constants
# 1 Hartree in various units
HARTREE_TO_EV: float = 27.211386245988
HARTREE_TO_KJ_MOL: float = 2625.4996394799
HARTREE_TO_KCAL_MOL: float = 627.5094740631
HARTREE_TO_CM_INV: float = 219474.63136320
HARTREE_TO_KELVIN: float = 315775.02480407
HARTREE_TO_JOULE: float = 4.3597447222071e-18
HARTREE_TO_MHZ: float = 6579683920.502

# Length conversions
BOHR_TO_ANGSTROM: float = 0.529177210903
ANGSTROM_TO_BOHR: float = 1.0 / BOHR_TO_ANGSTROM
ANGSTROM_TO_PM: float = 100.0
ANGSTROM_TO_NM: float = 0.1

# Dipole conversions (1 a.u. = e * a0 = 2.541746473 Debye)
AU_TO_DEBYE: float = 2.541746473
DEBYE_TO_AU: float = 1.0 / AU_TO_DEBYE
DEBYE_TO_COULOMB_METER: float = 3.33564095e-30

# Rotational constant conversions (1 cm^-1 = 29.9792458 GHz = 29979.2458 MHz)
CM_INV_TO_GHZ: float = 29.9792458
CM_INV_TO_MHZ: float = 29979.2458

# Moment of inertia conversions (1 amu * A^2 = 1.66053906660e-40 g * cm^2)
# Note: NIST CCCBDB reports gm cm^2 (where 1 amu = 1.66053906660e-24 g and 1 A^2 = 1e-16 cm^2)
AMU_ANGSTROM2_TO_G_CM2: float = 1.66053906660e-40
AMU_ANGSTROM2_TO_AU: float = 1822.888486209 * (ANGSTROM_TO_BOHR ** 2)

# Polarizability volume conversions
ANGSTROM3_TO_BOHR3: float = ANGSTROM_TO_BOHR ** 3

# Normalization aliases for flexible string input
_UNIT_ALIASES: Dict[str, str] = {
    # Energy
    "hartree": EnergyUnit.HARTREE.value,
    "hartrees": EnergyUnit.HARTREE.value,
    "eh": EnergyUnit.HARTREE.value,
    "ha": EnergyUnit.HARTREE.value,
    "a.u.": EnergyUnit.HARTREE.value,
    "ev": EnergyUnit.EV.value,
    "electronvolt": EnergyUnit.EV.value,
    "kj/mol": EnergyUnit.KJ_MOL.value,
    "kj mol-1": EnergyUnit.KJ_MOL.value,
    "kj mol^-1": EnergyUnit.KJ_MOL.value,
    "kcal/mol": EnergyUnit.KCAL_MOL.value,
    "kcal mol-1": EnergyUnit.KCAL_MOL.value,
    "kcal mol^-1": EnergyUnit.KCAL_MOL.value,
    "cm-1": EnergyUnit.CM_INV.value,
    "cm^-1": EnergyUnit.CM_INV.value,
    "wavenumber": EnergyUnit.CM_INV.value,
    "k": EnergyUnit.KELVIN.value,
    "kelvin": EnergyUnit.KELVIN.value,
    "j": EnergyUnit.JOULE.value,
    "joule": EnergyUnit.JOULE.value,
    "mhz": EnergyUnit.MHZ.value,
    # Length
    "angstrom": LengthUnit.ANGSTROM.value,
    "angstroms": LengthUnit.ANGSTROM.value,
    "a": LengthUnit.ANGSTROM.value,
    "å": LengthUnit.ANGSTROM.value,
    "bohr": LengthUnit.BOHR.value,
    "bohrs": LengthUnit.BOHR.value,
    "a0": LengthUnit.BOHR.value,
    "pm": LengthUnit.PICOMETER.value,
    "picometer": LengthUnit.PICOMETER.value,
    "nm": LengthUnit.NANOMETER.value,
    "nanometer": LengthUnit.NANOMETER.value,
    # Dipole
    "debye": DipoleUnit.DEBYE.value,
    "d": DipoleUnit.DEBYE.value,
    "c*m": DipoleUnit.COULOMB_METER.value,
    "c m": DipoleUnit.COULOMB_METER.value,
    "au": DipoleUnit.ATOMIC_UNIT.value,
    # Rotational
    "ghz": RotationalConstantUnit.GHZ.value,
    # Moment of Inertia
    "amu*angstrom^2": MomentOfInertiaUnit.AMU_ANGSTROM2.value,
    "amu a2": MomentOfInertiaUnit.AMU_ANGSTROM2.value,
    "amu å2": MomentOfInertiaUnit.AMU_ANGSTROM2.value,
    "amu*a^2": MomentOfInertiaUnit.AMU_ANGSTROM2.value,
    "g*cm^2": MomentOfInertiaUnit.G_CM2.value,
    "gm cm2": MomentOfInertiaUnit.G_CM2.value,
    "g cm2": MomentOfInertiaUnit.G_CM2.value,
    # Polarizability
    "angstrom^3": PolarizabilityUnit.ANGSTROM3.value,
    "a3": PolarizabilityUnit.ANGSTROM3.value,
    "å3": PolarizabilityUnit.ANGSTROM3.value,
    "bohr^3": PolarizabilityUnit.BOHR3.value,
}

_ENERGY_TO_HARTREE: Dict[str, float] = {
    EnergyUnit.HARTREE.value: 1.0,
    EnergyUnit.EV.value: 1.0 / HARTREE_TO_EV,
    EnergyUnit.KJ_MOL.value: 1.0 / HARTREE_TO_KJ_MOL,
    EnergyUnit.KCAL_MOL.value: 1.0 / HARTREE_TO_KCAL_MOL,
    EnergyUnit.CM_INV.value: 1.0 / HARTREE_TO_CM_INV,
    EnergyUnit.KELVIN.value: 1.0 / HARTREE_TO_KELVIN,
    EnergyUnit.JOULE.value: 1.0 / HARTREE_TO_JOULE,
    EnergyUnit.MHZ.value: 1.0 / HARTREE_TO_MHZ,
}

_LENGTH_TO_ANGSTROM: Dict[str, float] = {
    LengthUnit.ANGSTROM.value: 1.0,
    LengthUnit.BOHR.value: BOHR_TO_ANGSTROM,
    LengthUnit.PICOMETER.value: 1.0 / ANGSTROM_TO_PM,
    LengthUnit.NANOMETER.value: 1.0 / ANGSTROM_TO_NM,
}

_DIPOLE_TO_DEBYE: Dict[str, float] = {
    DipoleUnit.DEBYE.value: 1.0,
    DipoleUnit.COULOMB_METER.value: 1.0 / DEBYE_TO_COULOMB_METER,
    DipoleUnit.ATOMIC_UNIT.value: AU_TO_DEBYE,
}

_ROTATIONAL_TO_CM_INV: Dict[str, float] = {
    RotationalConstantUnit.CM_INV.value: 1.0,
    RotationalConstantUnit.GHZ.value: 1.0 / CM_INV_TO_GHZ,
    RotationalConstantUnit.MHZ.value: 1.0 / CM_INV_TO_MHZ,
}

_INERTIA_TO_AMU_A2: Dict[str, float] = {
    MomentOfInertiaUnit.AMU_ANGSTROM2.value: 1.0,
    MomentOfInertiaUnit.G_CM2.value: 1.0 / AMU_ANGSTROM2_TO_G_CM2,
    MomentOfInertiaUnit.ATOMIC_UNIT.value: 1.0 / AMU_ANGSTROM2_TO_AU,
}

_POLARIZABILITY_TO_A3: Dict[str, float] = {
    PolarizabilityUnit.ANGSTROM3.value: 1.0,
    PolarizabilityUnit.BOHR3.value: 1.0 / ANGSTROM3_TO_BOHR3,
}


def normalize_unit_name(unit: Union[str, Enum]) -> str:
    """Normalize a user- or HTML-supplied unit string to its canonical representation."""
    if isinstance(unit, Enum):
        return str(unit.value)
    cleaned = unit.strip()
    if cleaned in _UNIT_ALIASES:
        return _UNIT_ALIASES[cleaned]
    lowered = cleaned.lower()
    if lowered in _UNIT_ALIASES:
        return _UNIT_ALIASES[lowered]
    raise ValueError(f"Unsupported or unrecognized physical unit: '{unit}'")


def convert_energy(value: float, from_unit: Union[str, EnergyUnit], to_unit: Union[str, EnergyUnit]) -> float:
    """Convert an energy value between any supported quantum/thermochemical units."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src not in _ENERGY_TO_HARTREE or dst not in _ENERGY_TO_HARTREE:
        raise ValueError(f"Cannot convert energy between '{from_unit}' and '{to_unit}'")
    in_hartree = value * _ENERGY_TO_HARTREE[src]
    return in_hartree / _ENERGY_TO_HARTREE[dst]


def convert_length(value: float, from_unit: Union[str, LengthUnit], to_unit: Union[str, LengthUnit]) -> float:
    """Convert a length/distance value between Angstrom, Bohr, pm, or nm."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src not in _LENGTH_TO_ANGSTROM or dst not in _LENGTH_TO_ANGSTROM:
        raise ValueError(f"Cannot convert length between '{from_unit}' and '{to_unit}'")
    in_angstrom = value * _LENGTH_TO_ANGSTROM[src]
    return in_angstrom / _LENGTH_TO_ANGSTROM[dst]


def convert_dipole(value: float, from_unit: Union[str, DipoleUnit], to_unit: Union[str, DipoleUnit]) -> float:
    """Convert an electric dipole moment between Debye, C*m, and atomic units (e*a0)."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src not in _DIPOLE_TO_DEBYE or dst not in _DIPOLE_TO_DEBYE:
        raise ValueError(f"Cannot convert dipole between '{from_unit}' and '{to_unit}'")
    in_debye = value * _DIPOLE_TO_DEBYE[src]
    return in_debye / _DIPOLE_TO_DEBYE[dst]


def convert_rotational_constant(
    value: float,
    from_unit: Union[str, RotationalConstantUnit],
    to_unit: Union[str, RotationalConstantUnit],
) -> float:
    """Convert a rotational constant between cm^-1, GHz, and MHz."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src not in _ROTATIONAL_TO_CM_INV or dst not in _ROTATIONAL_TO_CM_INV:
        raise ValueError(f"Cannot convert rotational constant between '{from_unit}' and '{to_unit}'")
    in_cm = value * _ROTATIONAL_TO_CM_INV[src]
    return in_cm / _ROTATIONAL_TO_CM_INV[dst]


def convert_moment_of_inertia(
    value: float,
    from_unit: Union[str, MomentOfInertiaUnit],
    to_unit: Union[str, MomentOfInertiaUnit],
) -> float:
    """Convert a moment of inertia between amu*Angstrom^2, g*cm^2, and atomic units."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src not in _INERTIA_TO_AMU_A2 or dst not in _INERTIA_TO_AMU_A2:
        raise ValueError(f"Cannot convert moment of inertia between '{from_unit}' and '{to_unit}'")
    in_amu_a2 = value * _INERTIA_TO_AMU_A2[src]
    return in_amu_a2 / _INERTIA_TO_AMU_A2[dst]


def convert_polarizability(
    value: float,
    from_unit: Union[str, PolarizabilityUnit],
    to_unit: Union[str, PolarizabilityUnit],
) -> float:
    """Convert electric dipole polarizability volume between Angstrom^3 and Bohr^3."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src not in _POLARIZABILITY_TO_A3 or dst not in _POLARIZABILITY_TO_A3:
        raise ValueError(f"Cannot convert polarizability between '{from_unit}' and '{to_unit}'")
    in_a3 = value * _POLARIZABILITY_TO_A3[src]
    return in_a3 / _POLARIZABILITY_TO_A3[dst]


def convert_value(value: float, from_unit: str, to_unit: str) -> float:
    """Universal dispatcher that converts a numerical quantity between compatible units."""
    src = normalize_unit_name(from_unit)
    dst = normalize_unit_name(to_unit)
    if src == dst:
        return value
    for table, fn in [
        (_ENERGY_TO_HARTREE, convert_energy),
        (_LENGTH_TO_ANGSTROM, convert_length),
        (_DIPOLE_TO_DEBYE, convert_dipole),
        (_ROTATIONAL_TO_CM_INV, convert_rotational_constant),
        (_INERTIA_TO_AMU_A2, convert_moment_of_inertia),
        (_POLARIZABILITY_TO_A3, convert_polarizability),
    ]:
        if src in table and dst in table:
            return fn(value, src, dst)  # type: ignore[arg-type]
    raise ValueError(f"Incompatible or unsupported unit conversion from '{from_unit}' to '{to_unit}'")
