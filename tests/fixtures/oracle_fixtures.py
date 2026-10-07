"""Typed loaders for the oracle corpus in this directory.

Every file is produced by a committed script in ``references/extract`` from sources vendored
under ``references/vendor``; ``oracle_manifest.yaml`` indexes file, script and sources, and
``DATA_SOURCES.md`` describes each dataset and its checks.

Upstream stores each lattice as the Gram matrix of its bilinear form in the standard basis
of ``Z^n``; the preamble's category constructor takes that Gram directly,
``Lattices(ZZ)(record["gram"])``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Required, TypedDict

import yaml

FIXTURES_DIR = Path(__file__).parent

type Vector = list[int]
type Gram = list[list[int]]
type PermutationImages = tuple[int, ...]


class OrbitCounts(TypedDict, total=False):
    line_orbits: int
    plane_orbits: int
    flag_orbits: int


class TableRowSource(TypedDict, total=False):
    kind: Required[str]
    citation: Required[str]
    file: Required[str]
    table: Required[str]
    line: int
    lines: str


class TextSource(TypedDict):
    kind: str
    file: str
    lines: str


class CitedTextSource(TypedDict):
    kind: str
    citation: str
    file: str
    lines: str


class LineSource(TypedDict):
    kind: str
    file: str
    line: int


class Enriques87Case(TypedDict):
    case: int
    S: str
    polarization_orbit_length: int
    group_order: int
    line_orbits: int
    plane_orbits: int
    flag_orbits: int
    minimal_degree: int
    phi_h_min: int
    source: TableRowSource


class PolarizationOrbitCount(TypedDict):
    lattice: str
    two_d: int
    g: int
    primitive_vector_orbits: int
    gamma_conjugacy_classes: int
    primitive_vector_orbits_up_to_2d: int
    gamma_conjugacy_classes_up_to_2d: int
    source: TableRowSource


class IsometryPair(TypedDict, total=False):
    id: Required[str]
    gram1: Required[Gram]
    gram2: Required[Gram]
    signature: Required[list[int]]
    isometric: Required[bool]
    certificate: Required[str]
    witness: list[list[int]]
    witness_convention: str
    source: Required[dict[str, object]]


class AllcockRank3Lattice(TypedDict):
    id: str
    gram: Gram
    elementary_divisors: list[int]
    simple_roots: list[Vector]
    num_simple_roots: int
    weyl_group_id: int
    lattice_id: int
    convention: str
    source: FileSource


class DefiniteOrthogonalGroupOrder(TypedDict):
    id: str
    gram: Gram
    signature: list[int]
    orthogonal_group_order: int
    run_in_oscar_ci: bool
    source: list[LineSource]


class DefiniteIsometryTest(TypedDict):
    id: str
    gram_1: Gram
    gram_2: Gram
    signature_1: list[int]
    signature_2: list[int]
    isometric: bool
    source: list[LineSource]


class GroupIsomorphismType(TypedDict):
    isomorphism_type: str
    order: int


class RootLatticeGroups(TypedDict):
    id: str
    gram: Gram
    signature: list[int]
    stable_orthogonal_group: GroupIsomorphismType
    special_orthogonal_group: GroupIsomorphismType
    special_stable_orthogonal_group: GroupIsomorphismType
    source: list[LineSource]


class HeckeReducedAutomorphismGroupOrder(TypedDict):
    lattice: str
    hecke_reduced_automorphism_group_order: int
    source: LineSource


class OscarLatticeOracles(TypedDict):
    vinberg: list[dict[str, object]]
    discriminant_images: list[dict[str, object]]
    isometry_groups: list[dict[str, object]]
    definite_orthogonal_group_orders: list[DefiniteOrthogonalGroupOrder]
    definite_isometry_tests: list[DefiniteIsometryTest]
    root_lattice_groups: list[RootLatticeGroups]
    hecke_reduced_automorphism_group_orders: list[HeckeReducedAutomorphismGroupOrder]


class DTowerBoundary(TypedDict):
    id: str
    N: int
    lattice: str
    gram: Gram
    signature: list[int]
    group: str
    type_iii_points: int
    type_ii_curves: int
    type_ii_labels: list[str]
    incidences: list[list[str]]
    source: dict[str, TextSource]


class OrbitInvariant(TypedDict):
    invariant: str
    isomorphism_class: str


class AllcockEnriquesPeriodLattice(TypedDict):
    id: str
    lattice: str
    gram: Gram
    signature: list[int]
    group: str
    norm_minus_one_vector_orbits: int
    primitive_isotropic_vector_orbits: list[OrbitInvariant]
    isotropic_plane_orbits: list[OrbitInvariant]
    source: dict[str, TextSource]


class E10VectorOrbit(TypedDict):
    lattice: str
    h_squared: int
    phi: int | None
    stabilizer_image_index_in_O_E10_F2: int
    source: TableRowSource


class E10VectorOrbits(TypedDict):
    orbits: list[E10VectorOrbit]
    source: dict[str, dict[str, object]]


class LabelledRoot(TypedDict):
    label: int
    vector: Vector


class E10OrbitRepresentative(TypedDict):
    degree: int
    gamma_class: int
    g_coefficients: list[int]
    vector: Vector
    source_line: int


class E10FundamentalDomain(TypedDict):
    lattice: str
    gram: Gram
    simple_roots: list[LabelledRoot]
    extreme_rays: list[Vector]
    extreme_ray_gram: Gram
    orthogonal_group_equals_weyl_group: bool
    weyl_vector_norm: int
    orbit_representatives: list[E10OrbitRepresentative]
    source: dict[str, dict[str, object]]


class NamedMatrix(TypedDict):
    name: str
    matrix: list[list[int]]


class MertensPoint(TypedDict):
    point: Vector
    neighbours: int


class MertensGenerators(TypedDict):
    id: str
    gram: Gram
    signature: list[int]
    isometry_convention: str
    omega: str
    omega_generators: list[NamedMatrix]
    mertens_d_perfect_points: list[MertensPoint]
    source: dict[str, object]


class BinaryFormGenerators(TypedDict):
    id: str
    form: list[int]
    gram: Gram
    discriminant: int
    automorphism_group_generators: list[list[list[int]]]
    isometry_convention: str
    source: FileSource


class BinaryFormGeneratorCount(TypedDict):
    """How many generators Hecke's automorphism_group_generators returns; Hecke's output, not a group invariant."""

    id: str
    form: list[int]
    gram: Gram
    discriminant: int
    hecke_generator_count: int
    source: FileSource


class BinaryFormImproper(TypedDict):
    id: str
    form: list[int]
    gram: Gram
    discriminant: int
    has_improper_automorphism: bool
    source: FileSource


class BinaryFormAutomorphisms(TypedDict):
    explicit_generators: list[BinaryFormGenerators]
    generator_counts: list[BinaryFormGeneratorCount]
    improper_automorphisms: list[BinaryFormImproper]


class DTowerTypeCount(TypedDict):
    N: int
    type_ii_components: int
    type_ii_from_genus_of_d: int
    type_ii_from_even_unimodular: int
    type_iii_components: int
    group_note: str
    derivation: str
    source: dict[str, TextSource]


class GenusOfD(TypedDict):
    n: int
    root_sublattices: list[str]


class FHatComponent(TypedDict):
    label: str
    type: str
    dimension_in_F_hat: int
    geometric_meaning: str
    quartic_case: str
    source: TableRowSource


class DTower(TypedDict):
    boundary_pictures: list[DTowerBoundary]
    type_counts: list[DTowerTypeCount]
    genus_of_d: list[GenusOfD]
    f18_hat_type_ii: list[FHatComponent]


class LorentzianPerfectCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    has_isotropic: bool
    isotropic_count: int
    total_count: int
    source: FileSource


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
    group_g: PermutationGroupData
    group_h: PermutationGroupData
    vectors: list[Vector]
    source: FileSource


class Citation(TypedDict, total=False):
    kind: Required[str]
    citation: Required[str]
    chapter: int
    section: str
    statements: list[str]
    quote: str


class CitedGram(TypedDict):
    equation: str
    gram: Gram


class SpinorGenusPair(TypedDict):
    id: str
    genus: str
    determinant: int
    signature: list[int]
    form_a: CitedGram
    form_b: CitedGram
    same_genus: bool
    spinor_genera_in_genus: int
    classes_in_genus: int
    integrally_equivalent: bool
    source: Citation


class ConwaySloaneCases(TypedDict):
    spinor_genus_pair_determinant_minus_128: SpinorGenusPair


class FileSource(TypedDict, total=False):
    kind: Required[str]
    file: Required[str]
    index: int
    line: int
    point_numbering: str
    gram_constructor: str


class IsotropicDecisionCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    has_isotropic: bool
    source: FileSource


class ReflectiveFormCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    num_simple_roots: int
    source: FileSource


class LorentzianEquivalenceCase(TypedDict):
    id: str
    dimension: int
    mat1: Gram
    mat2: Gram
    transporter_witness: list[list[int]]
    witness_convention: str
    source: FileSource


class LorentzianStabilizerCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    generators: list[list[list[int]]]
    num_generators: int
    source: FileSource


class RootSystemCase(TypedDict):
    id: str
    dimension: int
    gram: Gram
    gram_source: FileSource
    num_roots: int
    roots: list[Vector]
    source: FileSource


class ClassifiedSimplices(TypedDict):
    dimension: int
    count: int
    simplices: list[list[Vector]]
    source: FileSource


class ClassificationSimplices(TypedDict):
    dim5: ClassifiedSimplices
    dim6: ClassifiedSimplices
    dim7: ClassifiedSimplices


class CiIndefiniteCompCase(TypedDict):
    id: str
    components: list[str]
    k_dim: int
    gram: Gram
    rank: int
    signature: list[int]
    source: FileSource


class TitsBuildingCounts(TypedDict):
    points: int
    curves: int
    edges: int


class DawesBuildingCase(TypedDict):
    id: str
    lattice: str
    gram: Gram
    signature: list[int]
    group: str
    building: TitsBuildingCounts
    source: dict[str, TextSource]


class FigureSource(TypedDict):
    kind: str
    file: str
    lines: str
    label: str


class DawesIncidenceGraph(TypedDict):
    """Point-curve incidence graph of a building, keyed by DawesBuildingCase.id."""

    building_id: str
    points: list[int]
    curves: list[int]
    point_curve_incidences: list[list[int]]
    source: FigureSource


class DawesIndexChain(TypedDict):
    id: str
    chain: list[str]
    indices: list[int]
    total_index: int
    source: TextSource


class DawesBuildings(TypedDict):
    buildings: list[DawesBuildingCase]
    incidence_graphs: list[DawesIncidenceGraph]
    index_chains: list[DawesIndexChain]


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
    source: CitedTextSource


class StableGroupData(TypedDict):
    internal_name: str
    construction: str
    source_notation: str
    counts: OrbitCounts
    source: CitedTextSource


class UnpolarizedEnriques(TypedDict):
    id: str
    lattice: LatticeSpec
    component_preserving_group: ComponentPreservingGroupData
    stable_component_preserving_group: StableGroupData


class K3UnimodularLattice(TypedDict):
    id: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    is_even: bool
    is_unimodular: bool
    primitive_vector_orbits: str
    tested_represented_norms: list[int]
    source: TextSource


class BailyBorelBoundary(TypedDict, total=False):
    type_iii_points: int
    type_ii_curves: Required[int]
    point_curve_incidences: int


class DegreeTwoSources(TypedDict):
    boundary: TextSource
    root_types: TextSource


class DegreeTwoK3(TypedDict):
    id: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    baily_borel_boundary: BailyBorelBoundary
    type_ii_root_types: list[str]
    source: DegreeTwoSources


class DegreeFourK3(TypedDict):
    id: str
    construction: str
    signature: list[int]
    rank: int
    det: int
    baily_borel_boundary: BailyBorelBoundary
    type_ii_generalised_types: list[str]
    source: TextSource


class K3ModularStrata(TypedDict):
    k3_unimodular_lattice: K3UnimodularLattice
    degree_two_polarized_k3: DegreeTwoK3
    degree_four_polarized_k3: DegreeFourK3


class CentralizerCase(TypedDict, total=False):
    id: Required[str]
    gram: Required[Gram]
    signature: Required[list[int]]
    isometry: Required[list[list[int]]]
    isometry_order: Required[int]
    isometry_convention: Required[str]
    centralizer_image_order: int
    centralizer_image_is_all_of_O_qL: bool
    source: Required[TextSource]


class InvolutionClassCount(TypedDict):
    id: str
    genus_representative_gram: Gram
    signature: list[int]
    characteristic_polynomial: str
    isometry_order: int
    classes_in_genus: int
    local_classes: int
    meaning: str
    source: TextSource


class LatticeClassCount(TypedDict):
    """How many classes of lattices with isometry (L, f) OSCAR's call returns, for an explicit lattice."""

    id: str
    gram: Gram
    signature: list[int]
    isometry_order: int
    oscar_call: str
    class_count: int
    recorded_properties: list[str]
    source: TextSource


class DefiniteClassCount(TypedDict):
    id: str
    gram: Gram
    signature: list[int]
    oscar_call: str
    class_count: int
    recorded_properties: list[str]
    source: TextSource


class HermitianGenusClassCount(TypedDict):
    """As LatticeClassCount, for hermitian-type classes specified by signature pairs and determinant."""

    id: str
    oscar_call: str
    class_count: int
    recorded_properties: list[str]
    source: TextSource


class IsometryCentralizers(TypedDict):
    centralizer_cases: list[CentralizerCase]
    definite_centralizer_cases: list[CentralizerCase]
    involution_classes: list[InvolutionClassCount]
    lattice_class_counts: list[LatticeClassCount]
    hermitian_genus_class_counts: list[HermitianGenusClassCount]
    definite_class_counts: list[DefiniteClassCount]


class OracleEntry(TypedDict, total=False):
    fixture_file: Required[str]
    extraction: Required[str | None]
    runtime: Required[str | None]
    kind: Required[str]
    sources: Required[list[str]]
    vendored: bool
    note: str


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


def load_dawes_buildings() -> DawesBuildings:
    """Dawes's published Tits buildings and index chain (arXiv:2205.10601, 2108.06236)."""
    with open(FIXTURES_DIR / "dawes_buildings.json", encoding="utf-8") as f:
        data: DawesBuildings = json.load(f)
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


def load_conway_sloane_cases() -> ConwaySloaneCases:
    """Conway--Sloane: the same-genus, non-equivalent ternary pair (51a)/(51b) of Chapter 15 §11,
    and the root-lattice automorphism orders of Chapter 4 (A2, A3, D4, E6, E7, E8)."""
    with open(FIXTURES_DIR / "conway_sloane_cases.json", encoding="utf-8") as f:
        data: ConwaySloaneCases = json.load(f)
        return data


def load_isometry_centralizers() -> IsometryCentralizers:
    """Centralizer images and involution-class counts recorded in OSCAR's tests."""
    with open(FIXTURES_DIR / "isometry_centralizers.json", encoding="utf-8") as f:
        data: IsometryCentralizers = json.load(f)
        return data


def load_enriques_polarization_orbits() -> list[PolarizationOrbitCount]:
    """O(U + E8(-1))-orbit counts of primitive vectors of norm 2..72 (Dutour Sikirić--Hulek)."""
    with open(FIXTURES_DIR / "enriques_polarization_orbits.json", encoding="utf-8") as f:
        data: list[PolarizationOrbitCount] = json.load(f)
        return data


def load_indefinite_isometry_pairs() -> list[IsometryPair]:
    """Indefinite pairs with certified isometry verdicts from Indefinite.jl and Hecke tests."""
    with open(FIXTURES_DIR / "indefinite_isometry_pairs.json", encoding="utf-8") as f:
        data: list[IsometryPair] = json.load(f)
        return data


def load_allcock_rank3_reflective() -> list[AllcockRank3Lattice]:
    """Allcock's 8,595 reflective Lorentzian lattices of rank 3 with their simple roots."""
    with open(FIXTURES_DIR / "allcock_rank3_reflective.json", encoding="utf-8") as f:
        data: list[AllcockRank3Lattice] = json.load(f)
        return data


def load_oscar_lattice_oracles() -> OscarLatticeOracles:
    """Vinberg roots, discriminant images and isometry-group orders recorded in OSCAR's tests."""
    with open(FIXTURES_DIR / "oscar_lattice_oracles.json", encoding="utf-8") as f:
        data: OscarLatticeOracles = json.load(f)
        return data


def load_dtower_boundaries() -> DTower:
    """Baily--Borel boundaries of F(N) = Gamma(N) \\ D(Lambda_N), Lambda_N = U^2 + D_{N-2} (Laza--O'Grady)."""
    with open(FIXTURES_DIR / "dtower_boundaries.json", encoding="utf-8") as f:
        data: DTower = json.load(f)
        return data


def load_allcock_i_2_10_orbits() -> AllcockEnriquesPeriodLattice:
    """Orbits of norm -1 vectors, isotropic lines and planes in I_{2,10} (Allcock, Cor. 3-4)."""
    with open(FIXTURES_DIR / "allcock_i_2_10_orbits.json", encoding="utf-8") as f:
        data: AllcockEnriquesPeriodLattice = json.load(f)
        return data


def load_e10_vector_orbits() -> E10VectorOrbits:
    """O(E10)-orbits of primitive vectors of norm 0..10 with stabilizer indices (Brandhorst--Gonzalez-Alonso)."""
    with open(FIXTURES_DIR / "e10_vector_orbits.json", encoding="utf-8") as f:
        data: E10VectorOrbits = json.load(f)
        return data


def load_e10_fundamental_domain() -> E10FundamentalDomain:
    """U + E8(-1): simple roots, chamber rays, and orbit representatives of norm <= 30 (DH section 3)."""
    with open(FIXTURES_DIR / "e10_fundamental_domain.json", encoding="utf-8") as f:
        data: E10FundamentalDomain = json.load(f)
        return data


def load_mertens_generators() -> MertensGenerators:
    """Generators of O^+(L) for Mertens's det -155 example (arXiv:1303.3478)."""
    with open(FIXTURES_DIR / "mertens_generators.json", encoding="utf-8") as f:
        data: MertensGenerators = json.load(f)
        return data


def load_binary_form_automorphisms() -> BinaryFormAutomorphisms:
    """Automorphism-group facts for indefinite binary forms recorded by Hecke's QuadBin tests."""
    with open(FIXTURES_DIR / "binary_form_automorphisms.json", encoding="utf-8") as f:
        data: BinaryFormAutomorphisms = json.load(f)
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
    with open(FIXTURES_DIR / "lorentzian_equivalence_145.json", encoding="utf-8") as f:
        data: list[LorentzianEquivalenceCase] = json.load(f)
        return data


def load_lorentzian_stabilizer_cases() -> list[LorentzianStabilizerCase]:
    """Three Lorentzian forms with full orthogonal-group generator sets (28B_LorentzianPerfStabEqui)."""
    with open(FIXTURES_DIR / "lorentzian_stabilizers_cases.json", encoding="utf-8") as f:
        data: list[LorentzianStabilizerCase] = json.load(f)
        return data


def load_root_systems() -> list[RootSystemCase]:
    """103 simple root systems of reflective forms (01_RatIntAutomorphy)."""
    with open(FIXTURES_DIR / "root_systems_103.json", encoding="utf-8") as f:
        data: list[RootSystemCase] = json.load(f)
        return data


def load_classification_simplices() -> ClassificationSimplices:
    """Classified simplices in dimensions 5, 6, and 7 (DATA/ClassificationSimplices)."""
    with open(FIXTURES_DIR / "classification_simplices.json", encoding="utf-8") as f:
        data: ClassificationSimplices = json.load(f)
        return data
