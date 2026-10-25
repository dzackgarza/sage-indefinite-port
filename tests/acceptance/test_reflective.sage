"""Fundamental domains of reflective Lorentzian lattices (Allcock's edgewalk).

Simple roots are compared up to isometry: two root systems agree when their Gram matrices
agree up to a simultaneous permutation of the roots. Recorded data: polyhedral_common's
8,821 reflective forms (root counts) and 103 explicit root systems, Allcock's 8,595
rank-3 lattices with explicit roots, OSCAR's Vinberg tests, and Dutour Sikirić--Hulek's
ten simple roots of U + E8(-1).
"""

import pytest
from dzack_research.preamble.all import ZZ as PreambleZZ
from dzack_research.preamble.all import HyperbolicLattices

from tests.acceptance.consumer import EDGEWALK, coordinates, lattice, require
from sage_indefinite_port.indefinite.edgewalk import _Edgewalk
from sage_indefinite_port.indefinite.edgewalk_rank2 import _OrientedBinaryCycle
from tests.fixtures.oracle_fixtures import (
    load_allcock_rank3_reflective,
    load_e10_fundamental_domain,
    load_oscar_lattice_oracles,
    load_reflective_forms,
    load_root_systems,
)


def _root_gram(gram, roots):
    form = matrix(ZZ, gram)
    vectors = [vector(ZZ, root) for root in roots]
    return [[int(left * form * right) for right in vectors] for left in vectors]


def _canonical_root_graph(root_gram):
    """The labelled root graph up to simultaneous permutation of the roots."""
    size = len(root_gram)
    graph = Graph([list(range(size)), [(i, j, root_gram[i][j]) for i in range(size) for j in range(i + 1, size) if root_gram[i][j]]], format="vertices_and_edges")
    norms = sorted({root_gram[i][i] for i in range(size)})
    partition = [[i for i in range(size) if root_gram[i][i] == norm] for norm in norms]
    return norms, [len(cell) for cell in partition], graph.canonical_label(partition=partition, edge_labels=True).copy(immutable=True)


def _hyperbolic(gram):
    """The lattice in HyperbolicLattices' signature (n, 1); a (1, n) Gram matrix is negated, which keeps coordinates."""
    negative = sum(1 for value in matrix(QQ, gram).eigenvalues() if value < 0)
    form = matrix(ZZ, gram)
    assert negative == 1 or negative == form.nrows() - 1, f"{gram} is not hyperbolic"
    oriented = form if negative == 1 else -form
    return HyperbolicLattices(PreambleZZ)(lattice([list(row) for row in oriented.rows()]))


def _computed_roots(gram):
    return [coordinates(root) for root in _hyperbolic(gram).allcock_edgewalk().simple_roots()]


def _assert_same_root_system(gram, recorded_roots) -> None:
    computed = _computed_roots(gram)
    assert len(computed) == len(recorded_roots)
    assert _canonical_root_graph(_root_gram(gram, computed)) == _canonical_root_graph(_root_gram(gram, recorded_roots))


def test_u_plus_46_nonreflective_certificate_terminates_after_three_vertices() -> None:
    """Pinned Allcock nonreflectivity branch: U + <46> has three vertex orbits."""
    require(EDGEWALK)

    report = _hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 46]]).allcock_edgewalk()

    assert not report.is_reflective()
    assert len(report.vertices()) == 3
    assert len(report.simple_roots()) == 0


def test_promised_go_left_preserves_left_supplement() -> None:
    """Allcock's rank-two Go-left branch must retain the opposite endpoint."""
    cycle = _OrientedBinaryCycle(lattice([[1, 0], [0, -2]]))

    assert cycle.promised_step(1, (1, 0), (0, 1)) == ((3, 2), (-2, -1))


def test_edgewalk_branch_representative_sample(monkeypatch) -> None:
    """Exercise cusp, finite-corner, next-wall and oriented binary-cycle routes.

    Frozen IDs ``allcock_rk3_4`` and ``allcock_rk3_2`` cover the isotropic
    cusp and rationally anisotropic branches; ``U+<2>`` is the independently
    recorded three-wall example with both finite and ideal corners, and
    ``U+<46>`` covers the nonreflective stop after three vertex orbits. This
    sample does not claim direct coverage of the rank-two public ``binary_walk``
    entry point or higher-rank starting-corner cases; the frozen full suite
    remains unchanged and still runs every recorded case.
    """
    require(EDGEWALK)
    observed = {}

    def instrument(cls, name):
        original = getattr(cls, name)

        def wrapped(self, *args, **kwargs):
            observed[f"{cls.__name__}.{name}"] = True
            return original(self, *args, **kwargs)

        monkeypatch.setattr(cls, name, wrapped)

    for name in (
        "initial_isotropic_vertex", "initial_vertex", "ordinary_vertex_roots",
        "cusp_root", "cusp_roots", "isotropic_next_vertex", "finite_extensions",
        "edge_step", "edge_directions",
    ):
        instrument(_Edgewalk, name)
    for name in ("cycle", "promised_step", "canonical_companion", "shorter_pair", "reduced_start_pair"):
        instrument(_OrientedBinaryCycle, name)

    cases = load_allcock_rank3_reflective()
    # The historical table reuses IDs; the exact Gram disambiguates the two
    # intended branch representatives without depending on list position.
    isotropic_case = next(
        case for case in cases
        if case["id"] == "allcock_rk3_4" and case["gram"] == [[0, -1, 0], [-1, 0, 0], [0, 0, 1]]
    )
    anisotropic_case = next(
        case for case in cases
        if case["id"] == "allcock_rk3_2" and case["gram"] == [[-3, 0, 0], [0, 1, 0], [0, 0, 1]]
    )
    assert _hyperbolic(isotropic_case["gram"]).allcock_edgewalk().is_reflective()
    assert _hyperbolic(anisotropic_case["gram"]).allcock_edgewalk().is_reflective()

    three_wall = _hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 2]]).allcock_edgewalk()
    assert len(three_wall.simple_roots()) == 3
    vertex_squares = tuple(
        int(ray.q())
        for vertex in three_wall.vertices()
        for ray in vertex.primitive_rays()
    )
    assert any(square == 0 for square in vertex_squares)
    assert any(square < 0 for square in vertex_squares)
    nonreflective = _hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 46]]).allcock_edgewalk()
    assert not nonreflective.is_reflective()
    assert len(nonreflective.vertices()) == 3

    required = {
        "_Edgewalk.initial_isotropic_vertex", "_Edgewalk.initial_vertex",
        "_Edgewalk.ordinary_vertex_roots", "_Edgewalk.cusp_root", "_Edgewalk.cusp_roots",
        "_Edgewalk.isotropic_next_vertex", "_Edgewalk.finite_extensions",
        "_Edgewalk.edge_step", "_Edgewalk.edge_directions",
        "_OrientedBinaryCycle.cycle", "_OrientedBinaryCycle.promised_step",
        "_OrientedBinaryCycle.canonical_companion", "_OrientedBinaryCycle.shorter_pair",
        "_OrientedBinaryCycle.reduced_start_pair",
    }
    missing = sorted(name for name in required if name not in observed)
    assert not missing, f"unexercised pinned edgewalk routes: {missing}"


@pytest.mark.parametrize("case", load_reflective_forms(), ids=lambda case: case["id"])
def test_simple_root_counts_match_polyhedral_common(case) -> None:
    require(EDGEWALK)

    report = _hyperbolic(case["gram"]).allcock_edgewalk()

    assert report.is_reflective()
    assert len(report.simple_roots()) == case["num_simple_roots"]


@pytest.mark.parametrize("case", load_root_systems(), ids=lambda case: case["id"])
def test_root_systems_match_polyhedral_common(case) -> None:
    require(EDGEWALK)

    _assert_same_root_system(case["gram"], case["roots"])


@pytest.mark.parametrize("case", load_allcock_rank3_reflective(), ids=lambda case: case["id"])
def test_root_systems_match_allcock(case) -> None:
    require(EDGEWALK)

    _assert_same_root_system(case["gram"], case["simple_roots"])


_VINBERG = {case["id"]: case for case in load_oscar_lattice_oracles()["vinberg"]}


@pytest.mark.parametrize("case", list(_VINBERG.values()), ids=lambda case: case["id"])
def test_simple_root_counts_match_oscar_vinberg(case) -> None:
    require(EDGEWALK)

    assert len(_computed_roots(case["gram"])) == case["num_simple_roots"]


def test_vinberg_coxeter_diagram_has_oscars_edges() -> None:
    require(EDGEWALK)

    case = _VINBERG["vinberg_scharlau_U_A2"]
    roots = _computed_roots(case["gram"])
    root_gram = _root_gram(case["gram"], roots)

    # Oscar.jl/test/NumberTheory/vinberg.jl:43-47 counts intersecting walls,
    # including orthogonal ones, by k^2 / (q_i q_j) <= 1.
    edges = sum(
        1
        for i in range(len(roots))
        for j in range(i + 1, len(roots))
        if root_gram[i][j] ** 2 <= root_gram[i][i] * root_gram[j][j]
    )

    assert edges == case["coxeter_diagram_edges"]


def test_vinberg_roots_match_oscars_chamber_up_to_isometry() -> None:
    require(EDGEWALK)

    case = _VINBERG["vinberg_diag_1_-1_-3"]

    _assert_same_root_system(case["gram"], case["simple_roots_of_test_chamber"])


def test_e10_simple_roots_match_dutour_sikiric_hulek() -> None:
    require(EDGEWALK)

    domain = load_e10_fundamental_domain()

    _assert_same_root_system(domain["gram"], [root["vector"] for root in domain["simple_roots"]])
