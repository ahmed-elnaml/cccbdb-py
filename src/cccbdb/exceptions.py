"""Custom exception hierarchy for cccbdb-py."""

from __future__ import annotations

from typing import List, Optional


class CCCBDBError(Exception):
    """Base exception for all cccbdb-py errors."""


class SpeciesNotFoundError(CCCBDBError):
    """Raised when a molecular formula, name, or identifier is not found in NIST CCCBDB."""

    def __init__(self, query: str, message: Optional[str] = None) -> None:
        self.query = query
        super().__init__(
            message or f"No chemical species matching '{query}' was found in NIST CCCBDB."
        )


class AmbiguousSpeciesError(CCCBDBError):
    """Raised when a formula or name matches multiple isomers, conformers, or charge states."""

    def __init__(self, query: str, candidates: List[str], message: Optional[str] = None) -> None:
        self.query = query
        self.candidates = candidates
        cand_str = ", ".join(candidates[:10])
        super().__init__(
            message
            or (
                f"Ambiguous query '{query}' matched {len(candidates)} chemical species "
                f"({cand_str}). Use db.search() or pass an explicit CAS RN / species ID and charge."
            )
        )


class CalculationNotFoundError(CCCBDBError):
    """Raised when a requested method/basis-set combination is unavailable for a species."""

    def __init__(
        self,
        species_id: str,
        method: Optional[str] = None,
        basis_set: Optional[str] = None,
        message: Optional[str] = None,
    ) -> None:
        self.species_id = species_id
        self.method = method
        self.basis_set = basis_set
        model_chem = f"{method}/{basis_set}" if basis_set else (method or "requested model")
        super().__init__(
            message
            or f"No calculation data found for species '{species_id}' at {model_chem}."
        )


class IncompleteGeometryError(CCCBDBError):
    """Raised when exporting to XYZ or ASE from a record lacking unambiguous 3D Cartesian coordinates."""


class StructureChangedError(CCCBDBError):
    """Raised when NIST CCCBDB HTML structure differs from expected parser landmarks."""


class SessionLostError(CCCBDBError):
    """Raised when NIST CCCBDB returns 'Sorry I lost your molecule' due to expired ASP session."""


class TransportError(CCCBDBError):
    """Raised when HTTP requests to NIST CCCBDB fail or time out."""


class CacheMissError(CCCBDBError):
    """Raised in offline mode when a requested page or record is not present in the local cache."""


class ValidationError(CCCBDBError):
    """Raised when scientific data validation (units, stoichiometry, or bounds) fails."""
