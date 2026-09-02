"""Typed loaders for the oracle corpus in this directory.

Provenance and extraction method for every file are recorded in ``DATA_SOURCES.md``; the
acceptance criteria binding fixtures to plan phases are in ``oracle_manifest.yaml``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Required, TypedDict

import yaml

FIXTURES_DIR = Path(__file__).parent

type Vector = list[int]
type Gram = list[list[int]]


class Provenance(TypedDict, total=False):
    kind: Required[str]
    paper: str
    arxiv: str
    table: int
    row: int
    section: str
    sections: list[str]
    theorem: str
    author: str
    code: str
    path: str
    repo: str


class GroupSpec(TypedDict, total=False):
    internal_name: Required[str]
    construction: Required[str]
    source_notation: str


class OrbitCounts(TypedDict, total=False):
    line_orbits: int
    plane_orbits: int
    flag_orbits: int


class Enriques87Case(TypedDict):
    case: int
    minimal_degree: int | None
    polarization_orbit_length: int | None
    line_orbits: int
    plane_orbits: int
    flag_orbits: int
    source: Provenance


class LorentzianPerfectCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    has_isotropic: bool
    isotropic_count: int | None
    total_count: int | None


class PermutationGroupData(TypedDict):
    num_points: int
    generators: list[list[int]]


class DoubleCosetCase(TypedDict):
    id: str
    degree: int
    big_group_order: int
    small_group_order: int
    num_vectors: int
    index: int
    group_h: PermutationGroupData
    group_g: PermutationGroupData
    vectors: list[Vector]
    source: Provenance


class IndefiniteJlDefiniteLeaf(TypedDict):
    id: str
    name: str
    source_type: str
    gram: Gram
    dimension: int
    expected_aut_order: int
    source: Provenance


class IndefiniteJlLorentzianPair(TypedDict):
    id: str
    name: str
    source_type: str
    matrix1: Gram
    matrix2: Gram
    dimension: int
    signature: list[int]
    is_equivalent: bool
    source: Provenance


class IndefiniteJlCases(TypedDict):
    definite_leaves: list[IndefiniteJlDefiniteLeaf]
    lorentzian_pairs: list[IndefiniteJlLorentzianPair]


class IsotropicDecisionCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    has_isotropic: bool


class ReflectiveFormCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    num_simple_roots: int


class LorentzianEquivalenceCase(TypedDict):
    id: str
    dimension: int
    mat1: Gram
    mat2: Gram
    transporter_witness: list[list[int]] | None


class LorentzianStabilizerCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    generators: list[list[list[int]]]
    num_generators: int


class RootSystemCase(TypedDict):
    id: str
    num_roots: int
    dimension: int
    roots: list[Vector]


class ClassifiedSimplices(TypedDict):
    dimension: int
    count: int
    simplices: list[list[Vector]]


class ClassificationSimplices(TypedDict):
    dim5: ClassifiedSimplices
    dim6: ClassifiedSimplices
    dim7: ClassifiedSimplices


class CiIndefiniteCompCase(TypedDict):
    id: str
    components: list[str]
    signature: list[int]
    rank: int
    k_dim: int
    source: Provenance


class DawesBuildingCase(TypedDict):
    id: str
    name: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    gram: Gram
    group: GroupSpec
    source: Provenance


class LatticeSpec(TypedDict):
    name: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    is_even: bool
    gram: Gram


class IsotropicRepresentative(TypedDict):
    id: str
    basis: list[Vector]
    description: str
    discriminant_stabilizer_index: int


class ComponentPreservingGroupData(TypedDict):
    internal_name: str
    construction: str
    source_notation: str
    counts: OrbitCounts
    line_representatives: list[IsotropicRepresentative]
    plane_representatives: list[IsotropicRepresentative]


class OrbitDecompositionNotes(TypedDict):
    lines: str
    planes: str


class StableGroupData(TypedDict):
    internal_name: str
    construction: str
    source_notation: str
    counts: OrbitCounts
    orbit_decompositions: OrbitDecompositionNotes


class UnpolarizedEnriques(TypedDict):
    id: str
    lattice: LatticeSpec
    component_preserving_group: ComponentPreservingGroupData
    stable_component_preserving_group: StableGroupData
    source: Provenance


class PrimitiveVectorTransitivity(TypedDict):
    tested_represented_norms: list[int]
    expected_orbits_per_norm: int


class K3UnimodularLattice(TypedDict):
    id: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    is_even: bool
    is_unimodular: bool
    discriminant_group_order: int
    stable_group_identity: str
    primitive_vector_transitivity: PrimitiveVectorTransitivity
    source: Provenance


class TypeIIComponent(TypedDict):
    index: int
    quotient_root_type: str
    cartan_type: list[str]


class PolarizedK3Stratum(TypedDict):
    id: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    is_even: bool
    group: GroupSpec
    baily_borel_counts: OrbitCounts
    type_ii_boundary_components: list[TypeIIComponent]
    source: Provenance


class AttwellDuvalCase(TypedDict):
    d: int
    omega_d: int
    expected_cusps: int


class AttwellDuvalFamily(TypedDict):
    id: str
    formula: str
    test_cases: list[AttwellDuvalCase]
    source: Provenance


class K3ModularStrata(TypedDict):
    k3_unimodular_lattice: K3UnimodularLattice
    degree_two_polarized_k3: PolarizedK3Stratum
    degree_four_polarized_k3: PolarizedK3Stratum
    attwell_duval_cusp_family: AttwellDuvalFamily


class StructuralRule(TypedDict):
    id: str
    theorem: str
    kind: str


class EnriquesInvolutionCase(TypedDict):
    id: str
    ambient_lattice: str
    fixed_sublattice: str
    anti_fixed_sublattice: str
    gluing_subgroup: str
    source: Provenance


class CentralizerInvolutions(TypedDict):
    orthogonal_direct_sum_rule: StructuralRule
    gluing_overlattice_rule: StructuralRule
    stable_intersection_rule: StructuralRule
    enriques_involution_k3: EnriquesInvolutionCase


class ManifestLattice(TypedDict, total=False):
    construction: Required[str]
    signature: list[int]
    rank: int
    det: int
    is_even: bool
    is_unimodular: bool


class ManifestGroup(TypedDict, total=False):
    construction: Required[str | dict[str, list[dict[str, str]]]]
    internal_name: str
    source_notation: str


class OracleEntry(TypedDict, total=False):
    id: Required[str]
    source: Required[Provenance]
    authoritative_outputs: Required[list[str]]
    lattice: ManifestLattice
    group: ManifestGroup
    group_family: dict[str, str | int]
    domain: dict[str, str]
    expected: dict[str, int | str | list[int] | list[str]]
    fixture_file: str
    formula: str


class OracleManifest(TypedDict):
    version: str
    description: str
    oracles: list[OracleEntry]


def load_oracle_manifest() -> OracleManifest:
    """The acceptance manifest: 13 oracle suites bound to fixtures and plan phases."""
    with open(FIXTURES_DIR / "oracle_manifest.yaml", encoding="utf-8") as f:
        data: OracleManifest = yaml.safe_load(f)
        return data


def load_unpolarized_enriques() -> UnpolarizedEnriques:
    """N = U + U(2) + E8(-2) with its published line and plane orbits and stabilizer indices."""
    with open(FIXTURES_DIR / "unpolarized_enriques.json", encoding="utf-8") as f:
        data: UnpolarizedEnriques = json.load(f)
        return data


def load_enriques_87_polarizations() -> list[Enriques87Case]:
    """The 87 numerical Enriques polarizations (Dutour Sikirić--Hulek, Tables 1 and 2)."""
    with open(FIXTURES_DIR / "enriques_87_polarizations.json", encoding="utf-8") as f:
        data: list[Enriques87Case] = json.load(f)
        return data


def load_k3_modular_strata() -> K3ModularStrata:
    """The K3 lattice and the degree 2 and 4 polarized boundary strata (Scattone, Jones, Attwell-Duval)."""
    with open(FIXTURES_DIR / "k3_modular_strata.json", encoding="utf-8") as f:
        data: K3ModularStrata = json.load(f)
        return data


def load_dawes_buildings() -> list[DawesBuildingCase]:
    """Independent Tits-building computations of Dawes for three lattices."""
    with open(FIXTURES_DIR / "dawes_buildings.json", encoding="utf-8") as f:
        data: list[DawesBuildingCase] = json.load(f)
        return data


def load_lorentzian_perfect_domains() -> list[LorentzianPerfectCase]:
    """40 Lorentzian forms with frozen perfect-domain cell counts (28B_LorentzianPerfStabEqui)."""
    with open(FIXTURES_DIR / "lorentzian_perfect_domains.json", encoding="utf-8") as f:
        data: list[LorentzianPerfectCase] = json.load(f)
        return data


def load_ci_indefinite_comp() -> list[CiIndefiniteCompCase]:
    """The six lattice families of the upstream metamorphic suite (19_IndefiniteComp)."""
    with open(FIXTURES_DIR / "ci_indefinite_comp.json", encoding="utf-8") as f:
        data: list[CiIndefiniteCompCase] = json.load(f)
        return data


def load_double_coset_cases() -> list[DoubleCosetCase]:
    """18 finite double-coset instances (DoubleCosets/DBL)."""
    with open(FIXTURES_DIR / "double_coset_cases.json", encoding="utf-8") as f:
        data: list[DoubleCosetCase] = json.load(f)
        return data


def load_indefinite_jl_cases() -> IndefiniteJlCases:
    """Indefinite.jl definite leaves with automorphism orders and Lorentzian pairs."""
    with open(FIXTURES_DIR / "indefinite_jl_cases.json", encoding="utf-8") as f:
        data: IndefiniteJlCases = json.load(f)
        return data


def load_centralizer_involutions() -> CentralizerInvolutions:
    """Structural identities for involution centralizers and the Enriques involution."""
    with open(FIXTURES_DIR / "centralizer_involutions.json", encoding="utf-8") as f:
        data: CentralizerInvolutions = json.load(f)
        return data


def load_isotropic_cases() -> list[IsotropicDecisionCase]:
    """8,821 Gram matrices with the existence of a nonzero integral isotropic vector (DATA/IsotropicCases)."""
    with open(FIXTURES_DIR / "isotropic_cases_8821.json", encoding="utf-8") as f:
        data: list[IsotropicDecisionCase] = json.load(f)
        return data


def load_reflective_forms() -> list[ReflectiveFormCase]:
    """8,821 reflective Lorentzian forms with their simple-root counts (20_Reflective/ListReflect)."""
    with open(FIXTURES_DIR / "reflective_forms_8821.json", encoding="utf-8") as f:
        data: list[ReflectiveFormCase] = json.load(f)
        return data


def load_lorentzian_equivalence_cases() -> list[LorentzianEquivalenceCase]:
    """145 pairs of equivalent Lorentzian forms with exact transporters (28B_LorentzianPerfStabEqui)."""
    with open(FIXTURES_DIR / "lorentzian_equivalence_146.json", encoding="utf-8") as f:
        data: list[LorentzianEquivalenceCase] = json.load(f)
        return data


def load_lorentzian_stabilizer_cases() -> list[LorentzianStabilizerCase]:
    """Three Lorentzian forms with full orthogonal-group generator sets (28B_LorentzianPerfStabEqui)."""
    with open(FIXTURES_DIR / "lorentzian_stabilizers_cases.json", encoding="utf-8") as f:
        data: list[LorentzianStabilizerCase] = json.load(f)
        return data


def load_root_systems() -> list[RootSystemCase]:
    """103 simple root systems of reflective forms (01_RatIntAutomorphy)."""
    with open(FIXTURES_DIR / "root_systems_56.json", encoding="utf-8") as f:
        data: list[RootSystemCase] = json.load(f)
        return data


def load_classification_simplices() -> ClassificationSimplices:
    """Classified simplices in dimensions 5, 6, and 7 (DATA/ClassificationSimplices)."""
    with open(FIXTURES_DIR / "classification_simplices.json", encoding="utf-8") as f:
        data: ClassificationSimplices = json.load(f)
        return data
