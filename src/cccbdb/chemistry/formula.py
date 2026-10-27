"""Chemical formula parser, Hill system formatter, atomic weights, and CAS RN validator."""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

# IUPAC / NIST Standard Atomic Weights (u / Da) for elements in CCCBDB (H to Xe + D/T)
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.00794,
    "D": 2.014101778,
    "T": 3.01604927,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3050,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973762,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955912,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938045,
    "Fe": 55.845,
    "Co": 58.933195,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.92160,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.96,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
}

ATOMIC_NUMBERS: Dict[str, int] = {
    "H": 1, "D": 1, "T": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35,
    "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Ru": 44, "Rh": 45,
    "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
}


def parse_formula(formula: str) -> Tuple[Dict[str, int], int]:
    """Parse a chemical formula (including parentheses and terminal ion charges) into composition and charge.

    Examples
    --------
    >>> parse_formula("CH3(CH2)4CH3")
    ({'C': 6, 'H': 14}, 0)
    >>> parse_formula("CO3--")
    ({'C': 1, 'O': 3}, -2)
    >>> parse_formula("NH4+")
    ({'N': 1, 'H': 4}, 1)
    """
    raw = formula.strip().replace(" ", "")
    if not raw:
        raise ValueError("Formula cannot be empty")

    # Extract trailing ionic charge (+, -, ++, --, +2, -2, 2+, 2-)
    charge = 0
    charge_match = re.search(r"([+-]+|\d+[+-]|[+-]\d+)$", raw)
    if charge_match:
        c_str = charge_match.group(1)
        raw = raw[: charge_match.start()]
        if c_str.startswith("+") or c_str.endswith("+"):
            digits = re.sub(r"\D", "", c_str)
            charge = int(digits) if digits else len(c_str)
        else:
            digits = re.sub(r"\D", "", c_str)
            charge = -(int(digits) if digits else len(c_str))

    stack: List[Counter[str]] = [Counter()]
    i = 0
    n = len(raw)
    while i < n:
        if raw[i] == "(":
            stack.append(Counter())
            i += 1
        elif raw[i] == ")":
            if len(stack) < 2:
                raise ValueError(f"Unmatched closing parenthesis in formula '{formula}'")
            group = stack.pop()
            i += 1
            m_match = re.match(r"(\d+)", raw[i:])
            mult = 1
            if m_match:
                mult = int(m_match.group(1))
                i += len(m_match.group(1))
            for elem, count in group.items():
                stack[-1][elem] += count * mult
        else:
            e_match = re.match(r"([A-Z][a-z]?)(\d*)", raw[i:])
            if not e_match:
                raise ValueError(f"Invalid token at '{raw[i:]}' in formula '{formula}'")
            elem, count_str = e_match.groups()
            if elem not in ATOMIC_WEIGHTS:
                raise ValueError(f"Unknown chemical element '{elem}' in formula '{formula}'")
            count = int(count_str) if count_str else 1
            stack[-1][elem] += count
            i += len(e_match.group(0))

    if len(stack) != 1:
        raise ValueError(f"Unclosed parenthesis in formula '{formula}'")

    return dict(stack[0]), charge


def to_hill_formula(formula: str, include_charge: bool = True) -> str:
    """Convert any chemical formula into standard Hill system order (C first, H/D second, then alphabetical)."""
    comp, charge = parse_formula(formula)
    elements: List[str] = []
    if "C" in comp:
        elements.append("C")
        for h_iso in ("H", "D", "T"):
            if h_iso in comp:
                elements.append(h_iso)
        for elem in sorted(k for k in comp if k not in ("C", "H", "D", "T")):
            elements.append(elem)
    else:
        elements = sorted(comp.keys())

    parts: List[str] = []
    for elem in elements:
        cnt = comp[elem]
        parts.append(f"{elem}{cnt if cnt > 1 else ''}")

    result = "".join(parts)
    if include_charge and charge != 0:
        sign = "+" if charge > 0 else "-"
        result += sign * abs(charge) if abs(charge) <= 2 else f"{sign}{abs(charge)}"
    return result


def calculate_molecular_weight(formula: str) -> float:
    """Calculate the molecular weight (in g/mol or u) for a chemical formula."""
    comp, _ = parse_formula(formula)
    return round(sum(ATOMIC_WEIGHTS[el] * cnt for el, cnt in comp.items()), 5)


def count_heavy_atoms(formula: str) -> int:
    """Return the number of non-hydrogen (Z > 2) heavy atoms in the formula."""
    comp, _ = parse_formula(formula)
    return sum(cnt for el, cnt in comp.items() if ATOMIC_NUMBERS.get(el, 0) > 2)


def normalize_cas(cas: str) -> str:
    """Normalize a CAS Registry Number to NIST CCCBDB integer string (no hyphens) and dashed string.

    Returns the integer string representation used in CCCBDB URLs (e.g., ``7732185`` for ``7732-18-5``).
    """
    digits = re.sub(r"[^0-9]", "", cas.strip())
    if len(digits) < 5 or len(digits) > 10:
        raise ValueError(f"Invalid CAS Registry Number length: '{cas}'")
    return digits


def format_cas_hyphenated(cas: str) -> str:
    """Format a CAS Registry Number with standard hyphens (``XXXXXXX-YY-Z``)."""
    digits = normalize_cas(cas)
    return f"{digits[:-3]}-{digits[-3:-1]}-{digits[-1]}"


def validate_cas_checksum(cas: str) -> bool:
    """Verify the CAS Registry Number check digit according to ACS CAS rules."""
    try:
        digits = normalize_cas(cas)
    except ValueError:
        return False
    check_digit = int(digits[-1])
    body = digits[:-1][::-1]
    total = sum((idx + 1) * int(ch) for idx, ch in enumerate(body))
    return (total % 10) == check_digit
