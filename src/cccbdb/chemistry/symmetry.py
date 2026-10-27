"""Molecular symmetry and Schoenflies point group analysis for NIST CCCBDB.

Handles normalization of NIST CCCBDB ASCII/HTML point group notations
(e.g., ``C*V`` -> ``Cinfv``, ``D*H`` -> ``Dinfh``, ``C2V`` -> ``C2v``),
rotational symmetry numbers (sigma) for statistical thermodynamics / rotational
entropy, inversion center presence, polarity rules, and chirality classification.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass(frozen=True)
class PointGroupInfo:
    """Structured representation of a molecular Schoenflies point group."""

    symbol: str  # Canonical Schoenflies symbol (e.g., "C2v", "Dinfh", "Td", "Oh")
    unicode_symbol: str  # Pretty Unicode symbol (e.g., "C₂ᵥ", "D∞h")
    family: str  # Non-axial, Cyclic, Dihedral, Linear, Cubic, Icosahedral, Spherical
    principal_axis_order: Optional[int]  # n in Cn/Dn, None for linear/spherical/non-axial
    rotational_symmetry_number: int  # sigma for rotational partition function q_rot
    is_linear: bool
    has_inversion_center: bool
    allows_permanent_dipole: bool
    is_chiral: bool
    hermann_mauguin: Optional[str] = None


_SUBSCRIPT_MAP = str.maketrans("0123456789vhsd", "₀₁₂₃₄₅₆₇₈₉ᵥₕₛ𝒹")


def normalize_point_group_symbol(raw: str) -> str:
    """Normalize NIST CCCBDB point group strings to canonical Schoenflies notation.

    Examples
    --------
    >>> normalize_point_group_symbol("C*V")
    'Cinfv'
    >>> normalize_point_group_symbol("D∞h")
    'Dinfh'
    >>> normalize_point_group_symbol("C2V")
    'C2v'
    >>> normalize_point_group_symbol("D3H")
    'D3h'
    """
    cleaned = raw.strip()
    if not cleaned:
        raise ValueError("Empty point group string")

    # Replace infinity representations (* or ∞ or inf)
    cleaned = cleaned.replace("∞", "inf").replace("*", "inf")
    upper = cleaned.upper()

    if upper in ("CINFV", "C_INF_V"):
        return "Cinfv"
    if upper in ("DINFH", "D_INF_H"):
        return "Dinfh"
    if upper in ("KH", "K_H", "ATOM", "R3", "SO(3)"):
        return "Kh"
    if upper == "CS":
        return "Cs"
    if upper == "CI":
        return "Ci"
    if upper == "C1":
        return "C1"
    if upper in ("TD", "T_D"):
        return "Td"
    if upper in ("TH", "T_H"):
        return "Th"
    if upper == "T":
        return "T"
    if upper in ("OH", "O_H"):
        return "Oh"
    if upper == "O":
        return "O"
    if upper in ("IH", "I_H"):
        return "Ih"
    if upper == "I":
        return "I"

    # Match axial groups: C2, C2v, C2h, D3, D3h, D3d, S4, etc.
    match = re.match(r"^([CDS])(\d+)([VHD]?)$", upper)
    if match:
        prefix, order_str, suffix = match.groups()
        return f"{prefix}{int(order_str)}{suffix.lower()}"

    raise ValueError(f"Unrecognized Schoenflies point group symbol: '{raw}'")


def parse_point_group(raw: str) -> PointGroupInfo:
    """Parse any point group string into a rich :class:`PointGroupInfo` descriptor."""
    symbol = normalize_point_group_symbol(raw)

    if symbol == "Kh":
        return PointGroupInfo(
            symbol="Kh",
            unicode_symbol="Kₕ",
            family="Spherical (Single Atom)",
            principal_axis_order=None,
            rotational_symmetry_number=1,
            is_linear=False,
            has_inversion_center=True,
            allows_permanent_dipole=False,
            is_chiral=False,
            hermann_mauguin="∞/∞m",
        )
    if symbol == "Cinfv":
        return PointGroupInfo(
            symbol="Cinfv",
            unicode_symbol="C∞v",
            family="Linear Asymmetric",
            principal_axis_order=None,
            rotational_symmetry_number=1,
            is_linear=True,
            has_inversion_center=False,
            allows_permanent_dipole=True,
            is_chiral=False,
            hermann_mauguin="∞m",
        )
    if symbol == "Dinfh":
        return PointGroupInfo(
            symbol="Dinfh",
            unicode_symbol="D∞h",
            family="Linear Centrosymmetric",
            principal_axis_order=None,
            rotational_symmetry_number=2,
            is_linear=True,
            has_inversion_center=True,
            allows_permanent_dipole=False,
            is_chiral=False,
            hermann_mauguin="∞/mm",
        )
    if symbol == "C1":
        return PointGroupInfo(
            symbol="C1",
            unicode_symbol="C₁",
            family="Non-axial",
            principal_axis_order=1,
            rotational_symmetry_number=1,
            is_linear=False,
            has_inversion_center=False,
            allows_permanent_dipole=True,
            is_chiral=True,
            hermann_mauguin="1",
        )
    if symbol == "Cs":
        return PointGroupInfo(
            symbol="Cs",
            unicode_symbol="Cₛ",
            family="Non-axial (Mirror)",
            principal_axis_order=1,
            rotational_symmetry_number=1,
            is_linear=False,
            has_inversion_center=False,
            allows_permanent_dipole=True,
            is_chiral=False,
            hermann_mauguin="m",
        )
    if symbol == "Ci":
        return PointGroupInfo(
            symbol="Ci",
            unicode_symbol="Cᵢ",
            family="Non-axial (Inversion)",
            principal_axis_order=1,
            rotational_symmetry_number=1,
            is_linear=False,
            has_inversion_center=True,
            allows_permanent_dipole=False,
            is_chiral=False,
            hermann_mauguin="-1",
        )

    # Cubic & Icosahedral groups
    cubic_specs = {
        "T": ("Tetrahedral", 12, False, True, "23"),
        "Td": ("Tetrahedral", 12, False, False, "-43m"),
        "Th": ("Tetrahedral", 12, True, False, "m-3"),
        "O": ("Octahedral", 24, False, True, "432"),
        "Oh": ("Octahedral", 24, True, False, "m-3m"),
        "I": ("Icosahedral", 60, False, True, "532"),
        "Ih": ("Icosahedral", 60, True, False, "-5-3m"),
    }
    if symbol in cubic_specs:
        fam, sigma, inv, chiral, hm = cubic_specs[symbol]
        u_sym = symbol[0] + (symbol[1:].translate(_SUBSCRIPT_MAP) if len(symbol) > 1 else "")
        return PointGroupInfo(
            symbol=symbol,
            unicode_symbol=u_sym,
            family=fam,
            principal_axis_order=None,
            rotational_symmetry_number=sigma,
            is_linear=False,
            has_inversion_center=inv,
            allows_permanent_dipole=False,
            is_chiral=chiral,
            hermann_mauguin=hm,
        )

    # Axial groups Cn, Cnv, Cnh, Dn, Dnh, Dnd, Sn
    match = re.match(r"^([CDS])(\d+)([vhd]?)$", symbol)
    if not match:
        raise ValueError(f"Unable to parse Schoenflies symbol: '{symbol}'")

    prefix, n_str, suffix = match.groups()
    n = int(n_str)
    u_sym = f"{prefix}{n_str.translate(_SUBSCRIPT_MAP)}{suffix.translate(_SUBSCRIPT_MAP)}"

    if prefix == "C":
        sigma = n
        if suffix == "":
            return PointGroupInfo(
                symbol=symbol,
                unicode_symbol=u_sym,
                family="Cyclic",
                principal_axis_order=n,
                rotational_symmetry_number=sigma,
                is_linear=False,
                has_inversion_center=False,
                allows_permanent_dipole=True,
                is_chiral=True,
                hermann_mauguin=str(n),
            )
        if suffix == "v":
            return PointGroupInfo(
                symbol=symbol,
                unicode_symbol=u_sym,
                family="Pyramidal / Cyclic-Vertical",
                principal_axis_order=n,
                rotational_symmetry_number=sigma,
                is_linear=False,
                has_inversion_center=False,
                allows_permanent_dipole=True,
                is_chiral=False,
                hermann_mauguin=f"{n}mm" if n % 2 == 0 else f"{n}m",
            )
        if suffix == "h":
            return PointGroupInfo(
                symbol=symbol,
                unicode_symbol=u_sym,
                family="Cyclic-Horizontal",
                principal_axis_order=n,
                rotational_symmetry_number=sigma,
                is_linear=False,
                has_inversion_center=(n % 2 == 0),
                allows_permanent_dipole=False,
                is_chiral=False,
                hermann_mauguin=f"{n}/m",
            )

    if prefix == "D":
        sigma = 2 * n
        if suffix == "":
            return PointGroupInfo(
                symbol=symbol,
                unicode_symbol=u_sym,
                family="Dihedral",
                principal_axis_order=n,
                rotational_symmetry_number=sigma,
                is_linear=False,
                has_inversion_center=False,
                allows_permanent_dipole=False,
                is_chiral=True,
                hermann_mauguin=f"{n}22" if n % 2 == 0 else f"{n}2",
            )
        if suffix == "h":
            return PointGroupInfo(
                symbol=symbol,
                unicode_symbol=u_sym,
                family="Dihedral-Horizontal",
                principal_axis_order=n,
                rotational_symmetry_number=sigma,
                is_linear=False,
                has_inversion_center=(n % 2 == 0),
                allows_permanent_dipole=False,
                is_chiral=False,
                hermann_mauguin=f"{n}/mmm" if n % 2 == 0 else f"-{2*n}m2",
            )
        if suffix == "d":
            return PointGroupInfo(
                symbol=symbol,
                unicode_symbol=u_sym,
                family="Dihedral-Dihedral",
                principal_axis_order=n,
                rotational_symmetry_number=sigma,
                is_linear=False,
                has_inversion_center=(n % 2 == 1),
                allows_permanent_dipole=False,
                is_chiral=False,
                hermann_mauguin=f"-{2*n}2m" if n % 2 == 0 else f"-{n}m",
            )

    if prefix == "S":
        sigma = n // 2 if n % 2 == 0 else n
        return PointGroupInfo(
            symbol=symbol,
            unicode_symbol=u_sym,
            family="Improper Rotation",
            principal_axis_order=n,
            rotational_symmetry_number=sigma,
            is_linear=False,
            has_inversion_center=(n % 4 == 2),
            allows_permanent_dipole=False,
            is_chiral=False,
            hermann_mauguin=f"-{n}",
        )

    raise ValueError(f"Unhandled point group: '{symbol}'")
