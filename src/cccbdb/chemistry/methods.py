"""Registry of computational quantum chemistry methods in NIST CCCBDB.

Maps canonical method names, aliases, method families (Hartree-Fock, DFT,
Moller-Plesset Perturbation Theory, Configuration Interaction, Coupled Cluster,
Composite Thermochemical, Semi-Empirical, Molecular Mechanics), core-correlation
settings (frozen-core vs. full-core), and verified NIST CCCBDB numeric method IDs.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional


class MethodCategory(str, Enum):
    """Quantum chemical or empirical method category in NIST CCCBDB."""

    HARTREE_FOCK = "Hartree-Fock"
    DFT = "Density Functional Theory"
    MP_PERTURBATION = "Moller-Plesset Perturbation"
    DOUBLE_HYBRID_DFT = "Double-Hybrid DFT"
    CONFIGURATION_INTERACTION = "Configuration Interaction"
    QUADRATIC_CI = "Quadratic Configuration Interaction"
    COUPLED_CLUSTER = "Coupled Cluster"
    COMPOSITE = "Composite Thermochemical"
    SEMI_EMPIRICAL = "Semi-Empirical"
    GROUP_ADDITIVITY = "Group Additivity"
    MOLECULAR_MECHANICS = "Molecular Mechanics"


@dataclass(frozen=True)
class ComputationalMethodSpec:
    """Specification of a computational method supported by NIST CCCBDB."""

    name: str  # Canonical CCCBDB method name (e.g., "B3LYP", "CCSD(T)", "MP2=FULL")
    category: MethodCategory
    nist_id: Optional[int]  # Numeric method parameter used in energy3x.asp / geom3x.asp
    predefined_basis: bool  # True for G1, G2, G3, G4, CBS-Q, AM1, PM3, PM6, UFF, etc.
    core_correlation: str  # "frozen", "full", or "n/a"
    ultrafine_grid: bool = False
    description: str = ""


# Verified NIST CCCBDB Method Catalog (Release 22)
METHODS_CATALOG: Dict[str, ComputationalMethodSpec] = {
    # Hartree-Fock
    "HF": ComputationalMethodSpec(
        name="HF",
        category=MethodCategory.HARTREE_FOCK,
        nist_id=1,
        predefined_basis=False,
        core_correlation="n/a",
        description="Restricted (closed-shell) or Unrestricted (open-shell) Hartree-Fock SCF",
    ),
    "ROHF": ComputationalMethodSpec(
        name="ROHF",
        category=MethodCategory.HARTREE_FOCK,
        nist_id=9,
        predefined_basis=False,
        core_correlation="n/a",
        description="Restricted Open-Shell Hartree-Fock",
    ),
    # Density Functional Theory (DFT)
    "LSDA": ComputationalMethodSpec(
        name="LSDA",
        category=MethodCategory.DFT,
        nist_id=51,
        predefined_basis=False,
        core_correlation="n/a",
        description="Local Spin Density Approximation (SVWN5)",
    ),
    "BLYP": ComputationalMethodSpec(
        name="BLYP",
        category=MethodCategory.DFT,
        nist_id=13,
        predefined_basis=False,
        core_correlation="n/a",
        description="Becke 88 exchange + Lee-Yang-Parr correlation GGA functional",
    ),
    "B1B95": ComputationalMethodSpec(
        name="B1B95",
        category=MethodCategory.DFT,
        nist_id=24,
        predefined_basis=False,
        core_correlation="n/a",
        description="One-parameter hybrid meta-GGA functional (Becke 88 + B95)",
    ),
    "B3LYP": ComputationalMethodSpec(
        name="B3LYP",
        category=MethodCategory.DFT,
        nist_id=8,
        predefined_basis=False,
        core_correlation="n/a",
        description="Becke 3-parameter hybrid exchange + Lee-Yang-Parr correlation functional",
    ),
    "B3LYPultrafine": ComputationalMethodSpec(
        name="B3LYPultrafine",
        category=MethodCategory.DFT,
        nist_id=68,
        predefined_basis=False,
        core_correlation="n/a",
        ultrafine_grid=True,
        description="B3LYP evaluated with UltraFine integration grid (99,590)",
    ),
    "B3PW91": ComputationalMethodSpec(
        name="B3PW91",
        category=MethodCategory.DFT,
        nist_id=14,
        predefined_basis=False,
        core_correlation="n/a",
        description="Becke 3-parameter hybrid exchange + Perdew-Wang 91 correlation",
    ),
    "mPW1PW91": ComputationalMethodSpec(
        name="mPW1PW91",
        category=MethodCategory.DFT,
        nist_id=25,
        predefined_basis=False,
        core_correlation="n/a",
        description="Modified Perdew-Wang 1-parameter hybrid functional",
    ),
    "M06-2X": ComputationalMethodSpec(
        name="M06-2X",
        category=MethodCategory.DFT,
        nist_id=55,
        predefined_basis=False,
        core_correlation="n/a",
        description="Minnesota 2006 global hybrid meta-GGA with 54% HF exchange",
    ),
    "PBEPBE": ComputationalMethodSpec(
        name="PBEPBE",
        category=MethodCategory.DFT,
        nist_id=52,
        predefined_basis=False,
        core_correlation="n/a",
        description="Perdew-Burke-Ernzerhof (PBE) non-empirical GGA functional",
    ),
    "PBEPBEultrafine": ComputationalMethodSpec(
        name="PBEPBEultrafine",
        category=MethodCategory.DFT,
        nist_id=69,
        predefined_basis=False,
        core_correlation="n/a",
        ultrafine_grid=True,
        description="PBE GGA evaluated with UltraFine integration grid",
    ),
    "PBE1PBE": ComputationalMethodSpec(
        name="PBE1PBE",
        category=MethodCategory.DFT,
        nist_id=54,
        predefined_basis=False,
        core_correlation="n/a",
        description="PBE0 hybrid functional (25% exact HF exchange + PBE)",
    ),
    "HSEh1PBE": ComputationalMethodSpec(
        name="HSEh1PBE",
        category=MethodCategory.DFT,
        nist_id=56,
        predefined_basis=False,
        core_correlation="n/a",
        description="Heyd-Scuseria-Ernzerhof range-separated screened hybrid functional (HSE06)",
    ),
    "TPSSh": ComputationalMethodSpec(
        name="TPSSh",
        category=MethodCategory.DFT,
        nist_id=58,
        predefined_basis=False,
        core_correlation="n/a",
        description="Tao-Perdew-Staroverov-Scuseria hybrid meta-GGA functional (10% HF exchange)",
    ),
    "wB97X-D": ComputationalMethodSpec(
        name="wB97X-D",
        category=MethodCategory.DFT,
        nist_id=60,
        predefined_basis=False,
        core_correlation="n/a",
        description="Chai-Head-Gordon long-range corrected hybrid DFT with empirical dispersion",
    ),
    "B97D3": ComputationalMethodSpec(
        name="B97D3",
        category=MethodCategory.DFT,
        nist_id=66,
        predefined_basis=False,
        core_correlation="n/a",
        description="Grimme B97 GGA with D3(BJ) dispersion correction",
    ),
    "B97D3ultrafine": ComputationalMethodSpec(
        name="B97D3ultrafine",
        category=MethodCategory.DFT,
        nist_id=74,
        predefined_basis=False,
        core_correlation="n/a",
        ultrafine_grid=True,
        description="B97D3 evaluated with UltraFine integration grid",
    ),
    # Double-Hybrid DFT
    "B2PLYP": ComputationalMethodSpec(
        name="B2PLYP",
        category=MethodCategory.DOUBLE_HYBRID_DFT,
        nist_id=59,
        predefined_basis=False,
        core_correlation="frozen",
        description="Grimme double-hybrid functional with perturbative MP2 correlation (frozen core)",
    ),
    "B2PLYP=FULL": ComputationalMethodSpec(
        name="B2PLYP=FULL",
        category=MethodCategory.DOUBLE_HYBRID_DFT,
        nist_id=67,
        predefined_basis=False,
        core_correlation="full",
        description="Grimme B2PLYP double-hybrid functional including all core electrons",
    ),
    "B2PLYP=FULLultrafine": ComputationalMethodSpec(
        name="B2PLYP=FULLultrafine",
        category=MethodCategory.DOUBLE_HYBRID_DFT,
        nist_id=75,
        predefined_basis=False,
        core_correlation="full",
        ultrafine_grid=True,
        description="B2PLYP=FULL evaluated with UltraFine integration grid",
    ),
    # Moller-Plesset Perturbation Theory
    "MP2": ComputationalMethodSpec(
        name="MP2",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=3,
        predefined_basis=False,
        core_correlation="frozen",
        description="Second-order Moller-Plesset perturbation theory (frozen core)",
    ),
    "MP2=FULL": ComputationalMethodSpec(
        name="MP2=FULL",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=2,
        predefined_basis=False,
        core_correlation="full",
        description="Second-order Moller-Plesset perturbation theory (all electrons correlated)",
    ),
    "PMP2": ComputationalMethodSpec(
        name="PMP2",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=10,
        predefined_basis=False,
        core_correlation="frozen",
        description="Spin-projected second-order Moller-Plesset perturbation theory",
    ),
    "ROMP2": ComputationalMethodSpec(
        name="ROMP2",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=11,
        predefined_basis=False,
        core_correlation="frozen",
        description="Restricted Open-Shell MP2",
    ),
    "MP3": ComputationalMethodSpec(
        name="MP3",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=36,
        predefined_basis=False,
        core_correlation="frozen",
        description="Third-order Moller-Plesset perturbation theory (frozen core)",
    ),
    "MP3=FULL": ComputationalMethodSpec(
        name="MP3=FULL",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=37,
        predefined_basis=False,
        core_correlation="full",
        description="Third-order Moller-Plesset perturbation theory (full core)",
    ),
    "MP4": ComputationalMethodSpec(
        name="MP4",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=4,
        predefined_basis=False,
        core_correlation="frozen",
        description="Fourth-order Moller-Plesset perturbation theory SDTQ (frozen core)",
    ),
    "MP4=FULL": ComputationalMethodSpec(
        name="MP4=FULL",
        category=MethodCategory.MP_PERTURBATION,
        nist_id=44,
        predefined_basis=False,
        core_correlation="full",
        description="Fourth-order Moller-Plesset perturbation theory SDTQ (full core)",
    ),
    # Configuration Interaction
    "CID": ComputationalMethodSpec(
        name="CID",
        category=MethodCategory.CONFIGURATION_INTERACTION,
        nist_id=16,
        predefined_basis=False,
        core_correlation="frozen",
        description="Configuration Interaction with Double excitations",
    ),
    "CID=FULL": ComputationalMethodSpec(
        name="CID=FULL",
        category=MethodCategory.CONFIGURATION_INTERACTION,
        nist_id=45,
        predefined_basis=False,
        core_correlation="full",
        description="Configuration Interaction with Double excitations (full core)",
    ),
    "CISD": ComputationalMethodSpec(
        name="CISD",
        category=MethodCategory.CONFIGURATION_INTERACTION,
        nist_id=15,
        predefined_basis=False,
        core_correlation="frozen",
        description="Configuration Interaction with Single and Double excitations",
    ),
    # Quadratic Configuration Interaction
    "QCISD": ComputationalMethodSpec(
        name="QCISD",
        category=MethodCategory.QUADRATIC_CI,
        nist_id=5,
        predefined_basis=False,
        core_correlation="frozen",
        description="Quadratic Configuration Interaction with Singles and Doubles",
    ),
    "QCISD(T)": ComputationalMethodSpec(
        name="QCISD(T)",
        category=MethodCategory.QUADRATIC_CI,
        nist_id=22,
        predefined_basis=False,
        core_correlation="frozen",
        description="QCISD with perturbative Triples correction (frozen core)",
    ),
    "QCISD(T)=FULL": ComputationalMethodSpec(
        name="QCISD(T)=FULL",
        category=MethodCategory.QUADRATIC_CI,
        nist_id=61,
        predefined_basis=False,
        core_correlation="full",
        description="QCISD(T) with all electrons correlated",
    ),
    "QCISD(TQ)": ComputationalMethodSpec(
        name="QCISD(TQ)",
        category=MethodCategory.QUADRATIC_CI,
        nist_id=64,
        predefined_basis=False,
        core_correlation="frozen",
        description="QCISD with perturbative Triples and Quadruples",
    ),
    "QCISD(TQ)=FULL": ComputationalMethodSpec(
        name="QCISD(TQ)=FULL",
        category=MethodCategory.QUADRATIC_CI,
        nist_id=65,
        predefined_basis=False,
        core_correlation="full",
        description="QCISD(TQ) with all electrons correlated",
    ),
    # Coupled Cluster
    "CCD": ComputationalMethodSpec(
        name="CCD",
        category=MethodCategory.COUPLED_CLUSTER,
        nist_id=6,
        predefined_basis=False,
        core_correlation="frozen",
        description="Coupled Cluster Doubles",
    ),
    "CCSD": ComputationalMethodSpec(
        name="CCSD",
        category=MethodCategory.COUPLED_CLUSTER,
        nist_id=7,
        predefined_basis=False,
        core_correlation="frozen",
        description="Coupled Cluster Singles and Doubles (frozen core)",
    ),
    "CCSD=FULL": ComputationalMethodSpec(
        name="CCSD=FULL",
        category=MethodCategory.COUPLED_CLUSTER,
        nist_id=62,
        predefined_basis=False,
        core_correlation="full",
        description="Coupled Cluster Singles and Doubles (all electrons correlated)",
    ),
    "CCSD(T)": ComputationalMethodSpec(
        name="CCSD(T)",
        category=MethodCategory.COUPLED_CLUSTER,
        nist_id=12,
        predefined_basis=False,
        core_correlation="frozen",
        description="Coupled Cluster Singles, Doubles, and perturbative Triples (gold standard, frozen core)",
    ),
    "CCSD(T)=FULL": ComputationalMethodSpec(
        name="CCSD(T)=FULL",
        category=MethodCategory.COUPLED_CLUSTER,
        nist_id=63,
        predefined_basis=False,
        core_correlation="full",
        description="Coupled Cluster Singles, Doubles, and perturbative Triples (all electrons correlated)",
    ),
    # Composite Thermochemical Methods (Predefined basis sets)
    "G1": ComputationalMethodSpec(
        name="G1",
        category=MethodCategory.COMPOSITE,
        nist_id=19,
        predefined_basis=True,
        core_correlation="full",
        description="Gaussian-1 composite thermochemical procedure",
    ),
    "G2": ComputationalMethodSpec(
        name="G2",
        category=MethodCategory.COMPOSITE,
        nist_id=20,
        predefined_basis=True,
        core_correlation="full",
        description="Gaussian-2 composite thermochemical procedure",
    ),
    "G2MP2": ComputationalMethodSpec(
        name="G2MP2",
        category=MethodCategory.COMPOSITE,
        nist_id=21,
        predefined_basis=True,
        core_correlation="full",
        description="Gaussian-2 (MP2 reduced-order variant) composite procedure",
    ),
    "G3": ComputationalMethodSpec(
        name="G3",
        category=MethodCategory.COMPOSITE,
        nist_id=23,
        predefined_basis=True,
        core_correlation="full",
        description="Gaussian-3 composite thermochemical procedure",
    ),
    "G3B3": ComputationalMethodSpec(
        name="G3B3",
        category=MethodCategory.COMPOSITE,
        nist_id=28,
        predefined_basis=True,
        core_correlation="full",
        description="Gaussian-3 using B3LYP geometries and zero-point energies",
    ),
    "G3MP2": ComputationalMethodSpec(
        name="G3MP2",
        category=MethodCategory.COMPOSITE,
        nist_id=29,
        predefined_basis=True,
        core_correlation="frozen",
        description="Gaussian-3 (MP2 variant) composite procedure",
    ),
    "G3MP2B3": ComputationalMethodSpec(
        name="G3MP2B3",
        category=MethodCategory.COMPOSITE,
        nist_id=30,
        predefined_basis=True,
        core_correlation="frozen",
        description="G3MP2 with B3LYP geometries and frequencies",
    ),
    "G4": ComputationalMethodSpec(
        name="G4",
        category=MethodCategory.COMPOSITE,
        nist_id=53,
        predefined_basis=True,
        core_correlation="full",
        description="Gaussian-4 composite thermochemical theory",
    ),
    "G4MP2": ComputationalMethodSpec(
        name="G4MP2",
        category=MethodCategory.COMPOSITE,
        nist_id=57,
        predefined_basis=True,
        core_correlation="frozen",
        description="Gaussian-4 (MP2 reduced-order variant)",
    ),
    "CBS-Q": ComputationalMethodSpec(
        name="CBS-Q",
        category=MethodCategory.COMPOSITE,
        nist_id=26,
        predefined_basis=True,
        core_correlation="full",
        description="Petersson Complete Basis Set quadratic CI composite method",
    ),
    "CBS-QB3": ComputationalMethodSpec(
        name="CBS-QB3",
        category=MethodCategory.COMPOSITE,
        nist_id=27,
        predefined_basis=True,
        core_correlation="full",
        description="Complete Basis Set CBS-QB3 composite method using B3LYP geometries",
    ),
    "W1U": ComputationalMethodSpec(
        name="W1U",
        category=MethodCategory.COMPOSITE,
        nist_id=70,
        predefined_basis=True,
        core_correlation="full",
        description="Weizmann-1 unrestricted thermochemical protocol",
    ),
    # Semi-Empirical
    "AM1": ComputationalMethodSpec(
        name="AM1",
        category=MethodCategory.SEMI_EMPIRICAL,
        nist_id=17,
        predefined_basis=True,
        core_correlation="n/a",
        description="Austin Model 1 semi-empirical NDDO molecular orbital method",
    ),
    "PM3": ComputationalMethodSpec(
        name="PM3",
        category=MethodCategory.SEMI_EMPIRICAL,
        nist_id=18,
        predefined_basis=True,
        core_correlation="n/a",
        description="Parameterized Model 3 semi-empirical method",
    ),
    "PM6": ComputationalMethodSpec(
        name="PM6",
        category=MethodCategory.SEMI_EMPIRICAL,
        nist_id=50,
        predefined_basis=True,
        core_correlation="n/a",
        description="Parameterized Model 6 semi-empirical method",
    ),
    "MNDOd": ComputationalMethodSpec(
        name="MNDOd",
        category=MethodCategory.SEMI_EMPIRICAL,
        nist_id=40,
        predefined_basis=True,
        core_correlation="n/a",
        description="Modified Neglect of Diatomic Overlap with d-orbitals",
    ),
}

_METHOD_ALIASES: Dict[str, str] = {
    "rhf": "HF",
    "uhf": "HF",
    "hartree-fock": "HF",
    "scf": "HF",
    "svwn": "LSDA",
    "svwn5": "LSDA",
    "lda": "LSDA",
    "pbe": "PBEPBE",
    "pbe0": "PBE1PBE",
    "hse06": "HSEh1PBE",
    "wb97xd": "wB97X-D",
    "ωb97x-d": "wB97X-D",
    "b97-d3": "B97D3",
    "mp2(fc)": "MP2",
    "mp2=fc": "MP2",
    "mp2(full)": "MP2=FULL",
    "mp2=fu": "MP2=FULL",
    "mp3(fc)": "MP3",
    "mp3(full)": "MP3=FULL",
    "mp4(fc)": "MP4",
    "mp4(full)": "MP4=FULL",
    "ccsd(full)": "CCSD=FULL",
    "ccsd(t)(fc)": "CCSD(T)",
    "ccsd(t)(full)": "CCSD(T)=FULL",
    "qcisd(t)(full)": "QCISD(T)=FULL",
    "g2(mp2)": "G2MP2",
    "g3(mp2)": "G3MP2",
}

_LOWER_TO_CANONICAL: Dict[str, str] = {k.lower(): k for k in METHODS_CATALOG}


def resolve_method(name: str) -> ComputationalMethodSpec:
    """Resolve a user-supplied method name or alias to its :class:`ComputationalMethodSpec`."""
    cleaned = name.strip()
    if cleaned in METHODS_CATALOG:
        return METHODS_CATALOG[cleaned]
    lower = cleaned.lower()
    if lower in _LOWER_TO_CANONICAL:
        return METHODS_CATALOG[_LOWER_TO_CANONICAL[lower]]
    if lower in _METHOD_ALIASES:
        return METHODS_CATALOG[_METHOD_ALIASES[lower]]
    raise ValueError(
        f"Unknown computational method '{name}'. Use list_methods() to inspect supported CCCBDB methods."
    )


def get_method_by_nist_id(nist_id: int) -> Optional[ComputationalMethodSpec]:
    """Look up a computational method by its numeric NIST CCCBDB ID."""
    for spec in METHODS_CATALOG.values():
        if spec.nist_id == nist_id:
            return spec
    return None


def list_methods(category: Optional[MethodCategory] = None) -> List[ComputationalMethodSpec]:
    """Return all registered computational methods, optionally filtered by category."""
    specs = list(METHODS_CATALOG.values())
    if category is not None:
        specs = [s for s in specs if s.category == category]
    return specs
