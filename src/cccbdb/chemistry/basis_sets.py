"""Registry of Gaussian basis sets, ECPs, and vibrational scaling factors in NIST CCCBDB.

Contains the exact 45 numeric basis set indices from NIST CCCBDB (`selectbasis1.asp`),
classification by family (Pople, Dunning correlation-consistent, Ahlrichs def2,
Sadlej, Effective Core Potentials), zeta level, polarization/diffuse/core-valence
attributes, and empirical harmonic vibrational scaling factors (`vsfx.asp`).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple


class BasisSetFamily(str, Enum):
    """Basis set family classification in NIST CCCBDB."""

    POPLE = "Pople"
    DUNNING = "Dunning"
    AHLRICHS = "Ahlrichs"
    SADLEJ = "Sadlej"
    ECP = "ECP"
    PREDEFINED = "Predefined"


@dataclass(frozen=True)
class BasisSetSpec:
    """Specification of an atomic orbital basis set in NIST CCCBDB."""

    name: str  # Canonical CCCBDB basis set name (e.g., "6-31G*", "aug-cc-pVTZ")
    nist_id: int  # Exact numeric basis ID from NIST CCCBDB selectbasis1.asp (0 for predefined)
    family: BasisSetFamily
    size_class: str  # "small", "medium", "large", or "predefined"
    zeta_quality: str  # "minimal", "DZ", "TZ", "QZ", "5Z", "predefined"
    has_polarization: bool
    has_diffuse: bool
    has_core_valence: bool = False
    has_tight_d: bool = False
    uses_ecp: bool = False
    description: str = ""


# Exact 45 basis sets from https://cccbdb.nist.gov/selectbasis1.asp + 0 (predefined)
BASIS_SETS_CATALOG: Dict[str, BasisSetSpec] = {
    "Predefined": BasisSetSpec(
        name="Predefined",
        nist_id=0,
        family=BasisSetFamily.PREDEFINED,
        size_class="predefined",
        zeta_quality="predefined",
        has_polarization=True,
        has_diffuse=True,
        description="Predefined basis set intrinsic to composite or semi-empirical method",
    ),
    "STO-3G": BasisSetSpec(
        name="STO-3G",
        nist_id=20,
        family=BasisSetFamily.POPLE,
        size_class="small",
        zeta_quality="minimal",
        has_polarization=False,
        has_diffuse=False,
        description="Minimal Slater-type orbital basis approximated by 3 Gaussian primitives",
    ),
    "3-21G": BasisSetSpec(
        name="3-21G",
        nist_id=2,
        family=BasisSetFamily.POPLE,
        size_class="small",
        zeta_quality="DZ",
        has_polarization=False,
        has_diffuse=False,
        description="Pople split-valence double-zeta basis set",
    ),
    "3-21G*": BasisSetSpec(
        name="3-21G*",
        nist_id=8,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        description="3-21G with d polarization functions on second-row (Na-Ar) heavy atoms",
    ),
    "6-31G": BasisSetSpec(
        name="6-31G",
        nist_id=9,
        family=BasisSetFamily.POPLE,
        size_class="small",
        zeta_quality="DZ",
        has_polarization=False,
        has_diffuse=False,
        description="Pople split-valence double-zeta basis set",
    ),
    "6-31G*": BasisSetSpec(
        name="6-31G*",
        nist_id=1,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        description="Pople 6-31G(d) polarized valence double-zeta basis set (6D on heavy atoms)",
    ),
    "6-31G**": BasisSetSpec(
        name="6-31G**",
        nist_id=3,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        description="Pople 6-31G(d,p) polarized double-zeta basis set (d on heavy, p on H)",
    ),
    "6-31+G**": BasisSetSpec(
        name="6-31+G**",
        nist_id=4,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=True,
        description="Pople 6-31+G(d,p) with diffuse sp functions on heavy atoms",
    ),
    "6-311G*": BasisSetSpec(
        name="6-311G*",
        nist_id=7,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        description="Pople 6-311G(d) polarized triple-split-valence basis set (5D)",
    ),
    "6-311G**": BasisSetSpec(
        name="6-311G**",
        nist_id=21,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        description="Pople 6-311G(d,p) polarized triple-split-valence basis set",
    ),
    "6-31G(2df,p)": BasisSetSpec(
        name="6-31G(2df,p)",
        nist_id=19,
        family=BasisSetFamily.POPLE,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        description="Pople 6-31G with 2d,1f polarization on heavy atoms and 1p on H (used in G3B3/G4)",
    ),
    "6-311+G(3df,2p)": BasisSetSpec(
        name="6-311+G(3df,2p)",
        nist_id=10,
        family=BasisSetFamily.POPLE,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        description="Large Pople triple-zeta basis with diffuse sp and 3d1f/2p polarization",
    ),
    "6-311+G(3df,2pd)": BasisSetSpec(
        name="6-311+G(3df,2pd)",
        nist_id=22,
        family=BasisSetFamily.POPLE,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        description="Large Pople triple-zeta basis with diffuse sp and 3d1f/2p1d polarization",
    ),
    # Dunning correlation-consistent basis sets
    "cc-pVDZ": BasisSetSpec(
        name="cc-pVDZ",
        nist_id=5,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        description="Dunning correlation-consistent polarized valence double-zeta",
    ),
    "cc-pVTZ": BasisSetSpec(
        name="cc-pVTZ",
        nist_id=6,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        description="Dunning correlation-consistent polarized valence triple-zeta",
    ),
    "cc-pVQZ": BasisSetSpec(
        name="cc-pVQZ",
        nist_id=25,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=False,
        description="Dunning correlation-consistent polarized valence quadruple-zeta",
    ),
    "cc-pV5Z": BasisSetSpec(
        name="cc-pV5Z",
        nist_id=39,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="5Z",
        has_polarization=True,
        has_diffuse=False,
        description="Dunning correlation-consistent polarized valence quintuple-zeta",
    ),
    "aug-cc-pVDZ": BasisSetSpec(
        name="aug-cc-pVDZ",
        nist_id=17,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=True,
        description="Augmented Dunning correlation-consistent polarized valence double-zeta",
    ),
    "aug-cc-pVTZ": BasisSetSpec(
        name="aug-cc-pVTZ",
        nist_id=18,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        description="Augmented Dunning correlation-consistent polarized valence triple-zeta",
    ),
    "aug-cc-pVQZ": BasisSetSpec(
        name="aug-cc-pVQZ",
        nist_id=26,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=True,
        description="Augmented Dunning correlation-consistent polarized valence quadruple-zeta",
    ),
    "aug-cc-pV5Z": BasisSetSpec(
        name="aug-cc-pV5Z",
        nist_id=40,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="5Z",
        has_polarization=True,
        has_diffuse=True,
        description="Augmented Dunning correlation-consistent polarized valence quintuple-zeta",
    ),
    "cc-pV(D+d)Z": BasisSetSpec(
        name="cc-pV(D+d)Z",
        nist_id=31,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        has_tight_d=True,
        description="cc-pVDZ with tight d function on second-row elements (Z > 12)",
    ),
    "cc-pV(T+d)Z": BasisSetSpec(
        name="cc-pV(T+d)Z",
        nist_id=27,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        has_tight_d=True,
        description="cc-pVTZ with tight d function on second-row elements",
    ),
    "cc-pV(Q+d)Z": BasisSetSpec(
        name="cc-pV(Q+d)Z",
        nist_id=34,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=False,
        has_tight_d=True,
        description="cc-pVQZ with tight d function on second-row elements",
    ),
    "aug-cc-pV(T+d)Z": BasisSetSpec(
        name="aug-cc-pV(T+d)Z",
        nist_id=28,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        has_tight_d=True,
        description="aug-cc-pVTZ with tight d function on second-row elements",
    ),
    "aug-cc-p(Q+d)Z": BasisSetSpec(
        name="aug-cc-p(Q+d)Z",
        nist_id=45,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=True,
        has_tight_d=True,
        description="aug-cc-pVQZ with tight d function on second-row elements",
    ),
    "cc-pCVDZ": BasisSetSpec(
        name="cc-pCVDZ",
        nist_id=23,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        has_core_valence=True,
        description="Core-valence correlation-consistent polarized double-zeta",
    ),
    "cc-pCVTZ": BasisSetSpec(
        name="cc-pCVTZ",
        nist_id=24,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        has_core_valence=True,
        description="Core-valence correlation-consistent polarized triple-zeta",
    ),
    "cc-pCVQZ": BasisSetSpec(
        name="cc-pCVQZ",
        nist_id=30,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=False,
        has_core_valence=True,
        description="Core-valence correlation-consistent polarized quadruple-zeta",
    ),
    "aug-cc-pCVTZ": BasisSetSpec(
        name="aug-cc-pCVTZ",
        nist_id=33,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        has_core_valence=True,
        description="Augmented core-valence correlation-consistent triple-zeta",
    ),
    "aug-cc-pCVQZ": BasisSetSpec(
        name="aug-cc-pCVQZ",
        nist_id=44,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=True,
        has_core_valence=True,
        description="Augmented core-valence correlation-consistent quadruple-zeta",
    ),
    "aug-cc-pCV5Z": BasisSetSpec(
        name="aug-cc-pCV5Z",
        nist_id=43,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="5Z",
        has_polarization=True,
        has_diffuse=True,
        has_core_valence=True,
        description="Augmented core-valence correlation-consistent quintuple-zeta",
    ),
    "daug-cc-pVDZ": BasisSetSpec(
        name="daug-cc-pVDZ",
        nist_id=36,
        family=BasisSetFamily.DUNNING,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=True,
        description="Doubly-augmented correlation-consistent polarized valence double-zeta",
    ),
    "daug-cc-pVTZ": BasisSetSpec(
        name="daug-cc-pVTZ",
        nist_id=35,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        description="Doubly-augmented correlation-consistent polarized valence triple-zeta",
    ),
    "daug-cc-pVQZ": BasisSetSpec(
        name="daug-cc-pVQZ",
        nist_id=37,
        family=BasisSetFamily.DUNNING,
        size_class="large",
        zeta_quality="QZ",
        has_polarization=True,
        has_diffuse=True,
        description="Doubly-augmented correlation-consistent polarized valence quadruple-zeta",
    ),
    # Ahlrichs & Sadlej
    "TZVP": BasisSetSpec(
        name="TZVP",
        nist_id=32,
        family=BasisSetFamily.AHLRICHS,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        description="Ahlrichs valence triple-zeta polarized basis set",
    ),
    "Def2TZVPP": BasisSetSpec(
        name="Def2TZVPP",
        nist_id=38,
        family=BasisSetFamily.AHLRICHS,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        description="Weigend-Ahlrichs def2-TZVPP triple-zeta doubly-polarized basis set",
    ),
    "Sadlej_pVTZ": BasisSetSpec(
        name="Sadlej_pVTZ",
        nist_id=29,
        family=BasisSetFamily.SADLEJ,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        description="Sadlej electric-property polarized triple-zeta basis set",
    ),
    # Effective Core Potential (ECP) basis sets
    "CEP-31G": BasisSetSpec(
        name="CEP-31G",
        nist_id=11,
        family=BasisSetFamily.ECP,
        size_class="small",
        zeta_quality="DZ",
        has_polarization=False,
        has_diffuse=False,
        uses_ecp=True,
        description="Stevens-Basch-Krauss Compact Effective Potential split-valence basis",
    ),
    "CEP-31G*": BasisSetSpec(
        name="CEP-31G*",
        nist_id=12,
        family=BasisSetFamily.ECP,
        size_class="small",
        zeta_quality="DZ",
        has_polarization=True,
        has_diffuse=False,
        uses_ecp=True,
        description="CEP-31G with d polarization functions",
    ),
    "CEP-121G": BasisSetSpec(
        name="CEP-121G",
        nist_id=13,
        family=BasisSetFamily.ECP,
        size_class="small",
        zeta_quality="TZ",
        has_polarization=False,
        has_diffuse=False,
        uses_ecp=True,
        description="Stevens-Basch-Krauss Compact Effective Potential triple-zeta basis",
    ),
    "CEP-121G*": BasisSetSpec(
        name="CEP-121G*",
        nist_id=14,
        family=BasisSetFamily.ECP,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        uses_ecp=True,
        description="CEP-121G with polarization functions",
    ),
    "LANL2DZ": BasisSetSpec(
        name="LANL2DZ",
        nist_id=15,
        family=BasisSetFamily.ECP,
        size_class="small",
        zeta_quality="DZ",
        has_polarization=False,
        has_diffuse=False,
        uses_ecp=True,
        description="Los Alamos National Laboratory double-zeta ECP basis set",
    ),
    "SDD": BasisSetSpec(
        name="SDD",
        nist_id=16,
        family=BasisSetFamily.ECP,
        size_class="medium",
        zeta_quality="DZ",
        has_polarization=False,
        has_diffuse=False,
        uses_ecp=True,
        description="Stuttgart/Dresden relativistic effective core potential basis set",
    ),
    "cc-pVTZ-PP": BasisSetSpec(
        name="cc-pVTZ-PP",
        nist_id=41,
        family=BasisSetFamily.ECP,
        size_class="medium",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=False,
        uses_ecp=True,
        description="Peterson correlation-consistent triple-zeta basis with relativistic pseudopotential",
    ),
    "aug-cc-pVTZ-PP": BasisSetSpec(
        name="aug-cc-pVTZ-PP",
        nist_id=42,
        family=BasisSetFamily.ECP,
        size_class="large",
        zeta_quality="TZ",
        has_polarization=True,
        has_diffuse=True,
        uses_ecp=True,
        description="Augmented correlation-consistent triple-zeta basis with relativistic pseudopotential",
    ),
}

_BASIS_ALIASES: Dict[str, str] = {
    "6-31g(d)": "6-31G*",
    "6-31g(d,p)": "6-31G**",
    "6-31+g(d,p)": "6-31+G**",
    "6-311g(d)": "6-311G*",
    "6-311g(d,p)": "6-311G**",
    "3-21g(d)": "3-21G*",
    "3-21g(*)": "3-21G*",
    "def2-tzvpp": "Def2TZVPP",
    "def2tzvpp": "Def2TZVPP",
    "sadlej-pvtz": "Sadlej_pVTZ",
    "sadlej pvtz": "Sadlej_pVTZ",
    "aug-cc-pv(q+d)z": "aug-cc-p(Q+d)Z",
    "none": "Predefined",
    "": "Predefined",
}

_LOWER_BASIS_MAP: Dict[str, str] = {k.lower(): k for k in BASIS_SETS_CATALOG}

# Verified NIST CCCBDB Vibrational Frequency Scaling Factors (from vsfx.asp)
VIBRATIONAL_SCALING_FACTORS: Dict[Tuple[str, str], float] = {
    ("HF", "STO-3G"): 0.8165,
    ("HF", "3-21G"): 0.9056,
    ("HF", "6-31G*"): 0.8985,
    ("HF", "6-31G**"): 0.9026,
    ("HF", "cc-pVDZ"): 0.9080,
    ("HF", "cc-pVTZ"): 0.9101,
    ("HF", "cc-pVQZ"): 0.9084,
    ("HF", "aug-cc-pVTZ"): 0.9104,
    ("B3LYP", "STO-3G"): 0.8924,
    ("B3LYP", "3-21G"): 0.9649,
    ("B3LYP", "6-31G*"): 0.9603,
    ("B3LYP", "6-31G**"): 0.9608,
    ("B3LYP", "6-31+G**"): 0.9642,
    ("B3LYP", "6-311G*"): 0.9663,
    ("B3LYP", "6-311G**"): 0.9668,
    ("B3LYP", "6-31G(2df,p)"): 0.9650,
    ("B3LYP", "6-311+G(3df,2p)"): 0.9670,
    ("B3LYP", "TZVP"): 0.9654,
    ("B3LYP", "Def2TZVPP"): 0.9626,
    ("B3LYP", "SDD"): 0.9613,
    ("B3LYP", "cc-pVDZ"): 0.9700,
    ("B3LYP", "cc-pVTZ"): 0.9666,
    ("B3LYP", "cc-pVQZ"): 0.9688,
    ("B3LYP", "aug-cc-pVDZ"): 0.9704,
    ("B3LYP", "aug-cc-pVTZ"): 0.9675,
    ("B3LYP", "aug-cc-pVQZ"): 0.9689,
    ("BLYP", "6-31G*"): 0.9919,
    ("BLYP", "cc-pVTZ"): 0.9970,
    ("BLYP", "aug-cc-pVTZ"): 0.9966,
    ("M06-2X", "6-31G*"): 0.9470,
    ("M06-2X", "cc-pVTZ"): 0.9551,
    ("M06-2X", "aug-cc-pVTZ"): 0.9557,
    ("wB97X-D", "6-31G*"): 0.9485,
    ("wB97X-D", "cc-pVTZ"): 0.9560,
    ("wB97X-D", "aug-cc-pVTZ"): 0.9566,
    ("wB97X-D", "Def2TZVPP"): 0.9552,
    ("PBEPBE", "6-31G*"): 0.9857,
    ("PBEPBE", "cc-pVTZ"): 0.9931,
    ("PBE1PBE", "6-31G*"): 0.9503,
    ("PBE1PBE", "cc-pVTZ"): 0.9610,
    ("MP2", "6-31G*"): 0.9430,
    ("MP2", "6-31G**"): 0.9365,
    ("MP2", "6-311+G(3df,2p)"): 0.9434,
    ("MP2", "cc-pVDZ"): 0.9525,
    ("MP2", "cc-pVTZ"): 0.9495,
    ("MP2", "cc-pVQZ"): 0.9484,
    ("MP2", "aug-cc-pVTZ"): 0.9529,
    ("MP2=FULL", "6-31G*"): 0.9422,
    ("MP2=FULL", "cc-pVTZ"): 0.9494,
    ("MP4", "6-31G"): 0.9438,
    ("MP4", "6-31G*"): 0.9548,
    ("MP4", "cc-pVTZ"): 0.9692,
    ("CCD", "6-311+G(3df,2p)"): 0.9431,
    ("CCSD", "6-31G*"): 0.9443,
    ("CCSD", "cc-pVTZ"): 0.9412,
    ("CCSD", "cc-pVQZ"): 0.9505,
    ("CCSD(T)", "6-31G*"): 0.9621,
    ("CCSD(T)", "cc-pVDZ"): 0.9788,
    ("CCSD(T)", "cc-pVTZ"): 0.9748,
    ("CCSD(T)", "cc-pVQZ"): 0.9697,
    ("CCSD(T)", "aug-cc-pVTZ"): 0.9699,
    ("CCSD(T)=FULL", "cc-pVTZ"): 0.9580,
    ("CCSD(T)=FULL", "aug-cc-pVQZ"): 0.9580,
    ("QCISD", "6-31G*"): 0.9523,
    ("QCISD", "cc-pVTZ"): 0.9555,
    ("QCISD(T)", "6-31G*"): 0.9593,
    ("QCISD(T)", "cc-pVTZ"): 0.9486,
    ("AM1", "Predefined"): 0.9535,
    ("PM3", "Predefined"): 0.9740,
    ("PM6", "Predefined"): 1.0624,
}


def resolve_basis_set(name: Optional[str]) -> BasisSetSpec:
    """Resolve a basis set name or alias (e.g. '6-31G(d)' -> '6-31G*') to its spec."""
    if name is None:
        return BASIS_SETS_CATALOG["Predefined"]
    cleaned = name.strip()
    if cleaned in BASIS_SETS_CATALOG:
        return BASIS_SETS_CATALOG[cleaned]
    lower = cleaned.lower()
    if lower in _LOWER_BASIS_MAP:
        return BASIS_SETS_CATALOG[_LOWER_BASIS_MAP[lower]]
    if lower in _BASIS_ALIASES:
        return BASIS_SETS_CATALOG[_BASIS_ALIASES[lower]]
    raise ValueError(
        f"Unknown basis set '{name}'. Use list_basis_sets() to inspect supported CCCBDB basis sets."
    )


def get_basis_set_by_nist_id(nist_id: int) -> Optional[BasisSetSpec]:
    """Look up a basis set by its numeric NIST CCCBDB ID."""
    for spec in BASIS_SETS_CATALOG.values():
        if spec.nist_id == nist_id:
            return spec
    return None


def list_basis_sets(family: Optional[BasisSetFamily] = None) -> List[BasisSetSpec]:
    """Return all registered basis sets, optionally filtered by family."""
    specs = [s for s in BASIS_SETS_CATALOG.values() if s.nist_id != 0]
    if family is not None:
        specs = [s for s in specs if s.family == family]
    return specs


def get_vibrational_scaling_factor(method: str, basis_set: Optional[str] = None) -> float:
    """Return the NIST CCCBDB harmonic vibrational frequency scaling factor (default 1.0)."""
    from cccbdb.chemistry.methods import resolve_method

    m_spec = resolve_method(method)
    b_spec = resolve_basis_set(basis_set) if not m_spec.predefined_basis else BASIS_SETS_CATALOG["Predefined"]
    return VIBRATIONAL_SCALING_FACTORS.get((m_spec.name, b_spec.name), 1.0)
