"""Shared BeautifulSoup parsing helpers for NIST CCCBDB HTML tables."""

from __future__ import annotations

import re
from typing import List, Optional, Tuple
from urllib.parse import parse_qs, urlparse

from bs4 import BeautifulSoup, Tag


def make_soup(html: str) -> BeautifulSoup:
    """Create a BeautifulSoup DOM tree using lxml if available, falling back to html.parser."""
    try:
        return BeautifulSoup(html, "lxml")
    except Exception:
        return BeautifulSoup(html, "html.parser")


def clean_text(node: Optional[Tag | str]) -> str:
    """Extract normalized, single-spaced text from a BeautifulSoup tag or string."""
    if node is None:
        return ""
    raw = node.get_text(" ", strip=True) if isinstance(node, Tag) else str(node)
    return re.sub(r"\s+", " ", raw).strip()


def parse_float_or_none(text: Optional[str]) -> Optional[float]:
    """Parse a floating-point number from text, returning ``None`` for blanks, dashes, or non-numerics.

    Never coerces missing values to ``0.0``.
    """
    if text is None:
        return None
    cleaned = text.strip().replace(",", "")
    if not cleaned or cleaned in ("-", "--", "N/A", "n/a", "?", "None"):
        return None
    # Strip trailing footnote markers or units if numeric prefix is unambiguous
    match = re.match(r"^([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)$", cleaned)
    if match:
        return float(match.group(1))
    return None


def parse_frequency_value(text: Optional[str]) -> Tuple[Optional[float], bool]:
    """Parse a vibrational frequency string, detecting imaginary modes suffixed with 'i'.

    Examples
    --------
    >>> parse_frequency_value("3052")
    (3052.0, False)
    >>> parse_frequency_value("872i")
    (-872.0, True)
    >>> parse_frequency_value("-872.4")
    (-872.4, True)
    """
    if text is None:
        return None, False
    cleaned = text.strip().replace(",", "")
    if not cleaned:
        return None, False
    imag_match = re.match(r"^([+-]?(?:\d+\.?\d*|\.\d+))\s*[iI]$", cleaned)
    if imag_match:
        val = abs(float(imag_match.group(1)))
        return -val, True
    val = parse_float_or_none(cleaned)
    if val is not None:
        return val, (val < 0.0)
    return None, False


def extract_species_header(soup: BeautifulSoup) -> Tuple[Optional[str], Optional[str]]:
    """Extract ``(formula, name)`` from CCCBDB page headings such as:
    - ``Experimental data for CO (Carbon monoxide)``
    - ``All results from a given calculation for CH2 (Methylene)``
    - ``Energies for H2O (Water)``
    """
    full_text = soup.get_text("\n", strip=True)
    patterns = [
        r"(?:Experimental data|All results from a given calculation|Energies|Geometry|Vibrational Frequencies|Dipole moments?)\s+for\s+([^\(\n]+?)\s*\(([^\)\n]*)\)",
        r"for\s+([A-Za-z0-9\+\-\(\)]+)\s*\(([^\)\n]+)\)",
    ]
    for pat in patterns:
        m = re.search(pat, full_text, re.IGNORECASE)
        if m:
            raw_formula = re.sub(r"\s+", "", m.group(1))
            raw_name = m.group(2).strip()
            if raw_formula:
                return raw_formula, (raw_name or raw_formula)
    return None, None


def extract_casno_from_url(url: str) -> Optional[str]:
    """Extract ``casno`` query parameter from a CCCBDB URL if present."""
    parsed = urlparse(url)
    qs = parse_qs(parsed.query)
    if "casno" in qs and qs["casno"]:
        return re.sub(r"\D", "", qs["casno"][0])
    if "_ctx_casno" in qs and qs["_ctx_casno"]:
        return re.sub(r"\D", "", qs["_ctx_casno"][0])
    return None
