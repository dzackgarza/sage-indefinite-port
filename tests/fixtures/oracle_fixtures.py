"""
Typed loader module for oracle corpus and acceptance criteria fixtures.

Provides access to:
- Oracle Manifest (Section 6)
- 8,821 Isotropic Subspace Decision Cases (CI_tests/DATA/IsotropicCases)
- 8,821 Reflective Lorentzian Forms (CI_tests/20_Reflective/ListReflect)
- 145 Lorentzian Equivalence Matrix Pairs with Exact Transporters (CI_tests/13_LorentzianStabEqui)
- 3 Lorentzian Stabilizer Group Cases (CI_tests/13_LorentzianStabEqui)
- 103 Root Systems of Reflective Forms (CI_tests/02A_IntegralAutomorphy)
- Classified Simplices in Dimensions 5, 6, 7 (CI_tests/DATA/ClassificationSimplices)
- Unpolarized Enriques Lattice and Stabilizers (Section 3.2, 3.3)
- All 87 Polarized Enriques Arithmetic Groups (Section 3.4)
- K3 Lattice and Polarized Modular Strata (Sections 3.1, 3.5, 3.6, 3.7)
- Matthew Dawes Independent Tits Buildings (Section 4)
- 40 Lorentzian Perfect Domain Enumeration Forms (Section 5.3)
- CI_tests / 19_IndefiniteComp Metamorphic Suite (Section 5.1)
- 18 Double Coset Finite Group Cases (Section 5.2)
- Indefinite.jl Leaf and Lorentzian Matrix Pairs (Section 5.5)
- Centralizer Involution Rules (Section Phase 8)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, TypedDict

FIXTURES_DIR = Path(__file__).parent


class Enriques87Case(TypedDict):
    case: int
    minimal_degree: int | None
    polarization_orbit_length: int | None
    line_orbits: int
    plane_orbits: int
    flag_orbits: int
    source: dict[str, Any]


class LorentzianPerfectCase(TypedDict):
    id: str
    dimension: int
    gram: list[list[int]]
    has_isotropic: bool
    isotropic_count: int | None
    total_count: int | None


class DoubleCosetCase(TypedDict):
    id: str
    degree: int
    big_group_order: int
    small_group_order: int
    num_vectors: int
    index: int
    group_h: dict[str, Any]
    group_g: dict[str, Any]
    vectors: list[list[int]]
    source: dict[str, Any]


class IndefiniteJlLorPair(TypedDict):
    id: str
    name: str
    source_type: str
    matrix1: list[list[int]]
    matrix2: list[list[int]]
    dimension: int
    signature: list[int]
    is_equivalent: bool
    source: dict[str, Any]


class IndefiniteJlDefiniteLeaf(TypedDict):
    id: str
    name: str
    source_type: str
    dimension: int
    source: dict[str, Any]


class IsotropicDecisionCase(TypedDict):
    id: str
    dimension: int
    gram: list[list[int]]
    has_isotropic: bool


class ReflectiveFormCase(TypedDict):
    id: str
    dimension: int
    gram: list[list[int]]
    num_simple_roots: int


class LorentzianEquivalenceCase(TypedDict):
    id: str
    dimension: int
    mat1: list[list[int]]
    mat2: list[list[int]]
    transporter_witness: list[list[int]] | None


def load_oracle_manifest() -> dict[str, Any]:
    """Load the authoritative oracle manifest."""
    with open(FIXTURES_DIR / "oracle_manifest.json", "r") as f:
        data: dict[str, Any] = json.load(f)
        return data


def load_unpolarized_enriques() -> dict[str, Any]:
    """Load unpolarized Enriques lattice N = U + U(2) + E8(-2) fixture."""
    with open(FIXTURES_DIR / "unpolarized_enriques.json", "r") as f:
        data: dict[str, Any] = json.load(f)
        return data


def load_enriques_87_polarizations() -> list[Enriques87Case]:
    """Load the full 87-row Enriques polarization dataset (Dutour Sikiric-Hulek)."""
    with open(FIXTURES_DIR / "enriques_87_polarizations.json", "r") as f:
        data: list[Enriques87Case] = json.load(f)
        return data


def load_k3_modular_strata() -> dict[str, Any]:
    """Load K3 lattice and polarized modular boundary strata (Jones, Scattone, Attwell-Duval)."""
    with open(FIXTURES_DIR / "k3_modular_strata.json", "r") as f:
        data: dict[str, Any] = json.load(f)
        return data


def load_dawes_buildings() -> list[dict[str, Any]]:
    """Load Matthew Dawes independent Tits-building test cases."""
    with open(FIXTURES_DIR / "dawes_buildings.json", "r") as f:
        data: list[dict[str, Any]] = json.load(f)
        return data


def load_lorentzian_perfect_domains() -> list[LorentzianPerfectCase]:
    """Load 40 Lorentzian perfect-domain forms with isotropic and total counts."""
    with open(FIXTURES_DIR / "lorentzian_perfect_domains.json", "r") as f:
        data: list[LorentzianPerfectCase] = json.load(f)
        return data


def load_ci_indefinite_comp() -> list[dict[str, Any]]:
    """Load upstream CI_tests/19_IndefiniteComp metamorphic forms."""
    with open(FIXTURES_DIR / "ci_indefinite_comp.json", "r") as f:
        data: list[dict[str, Any]] = json.load(f)
        return data


def load_double_coset_cases() -> list[DoubleCosetCase]:
    """Load finite double-coset test cases from CI_tests/DoubleCosets/DBL."""
    with open(FIXTURES_DIR / "double_coset_cases.json", "r") as f:
        data: list[DoubleCosetCase] = json.load(f)
        return data


def load_indefinite_jl_cases() -> dict[str, Any]:
    """Load Indefinite.jl definite leaves and Lorentzian test pairs."""
    with open(FIXTURES_DIR / "indefinite_jl_cases.json", "r") as f:
        data: dict[str, Any] = json.load(f)
        return data


def load_centralizer_involutions() -> dict[str, Any]:
    """Load centralizer and involution structural test cases."""
    with open(FIXTURES_DIR / "centralizer_involutions.json", "r") as f:
        data: dict[str, Any] = json.load(f)
        return data


def load_isotropic_cases() -> list[IsotropicDecisionCase]:
    """Load 8,821 concrete isotropic decision test cases from DATA/IsotropicCases."""
    with open(FIXTURES_DIR / "isotropic_cases_8821.json", "r") as f:
        data: list[IsotropicDecisionCase] = json.load(f)
        return data


def load_reflective_forms() -> list[ReflectiveFormCase]:
    """Load 8,821 concrete reflective forms from 20_Reflective/ListReflect."""
    with open(FIXTURES_DIR / "reflective_forms_8821.json", "r") as f:
        data: list[ReflectiveFormCase] = json.load(f)
        return data


def load_lorentzian_equivalence_cases() -> list[LorentzianEquivalenceCase]:
    """Load 145 Lorentzian equivalence matrix pairs with exact transporter witnesses."""
    with open(FIXTURES_DIR / "lorentzian_equivalence_146.json", "r") as f:
        data: list[LorentzianEquivalenceCase] = json.load(f)
        return data


def load_lorentzian_stabilizer_cases() -> list[dict[str, Any]]:
    """Load Lorentzian stabilizer group cases from TestCasesStab."""
    with open(FIXTURES_DIR / "lorentzian_stabilizers_cases.json", "r") as f:
        data: list[dict[str, Any]] = json.load(f)
        return data


def load_root_systems() -> list[dict[str, Any]]:
    """Load 103 root systems of reflective forms from 02A_IntegralAutomorphy."""
    with open(FIXTURES_DIR / "root_systems_56.json", "r") as f:
        data: list[dict[str, Any]] = json.load(f)
        return data


def load_classification_simplices() -> dict[str, Any]:
    """Load classified simplices in dimensions 5, 6, 7 from DATA/ClassificationSimplices."""
    with open(FIXTURES_DIR / "classification_simplices.json", "r") as f:
        data: dict[str, Any] = json.load(f)
        return data
