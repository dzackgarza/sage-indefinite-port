"""Fundamental domains of reflective Lorentzian lattices (Allcock's edgewalk).

Simple roots are compared up to isometry: two root systems agree when their Gram matrices
agree up to a simultaneous permutation of the roots. Recorded data: polyhedral_common's
8,821 reflective forms (root counts) and 103 explicit root systems, Allcock's 8,595
rank-3 lattices with explicit roots, OSCAR's Vinberg tests, and Dutour Sikirić--Hulek's
ten simple roots of U + E8(-1).
"""

import pytest
from dzack_research.preamble.all import ZZ as PreambleZZ
from dzack_research.preamble.all import HyperbolicLattices, Lattices
from sage.all import QQ, ZZ, Graph, matrix, vector

import sage_indefinite_port.indefinite.edgewalk as edgewalk_module
from sage_indefinite_port.indefinite.edgewalk import _Edgewalk
from sage_indefinite_port.indefinite.edgewalk_rank2 import _OrientedBinaryCycle
from tests.acceptance.consumer import EDGEWALK, coordinates, lattice, require
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
    return norms, [len(cell) for cell in partition], graph.canonical_label(partition=partition, edge_labels=True, algorithm="sage").copy(immutable=True)


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


def test_edgewalk_period_action_matrix_uses_coordinate_label_order() -> None:
    """The induced binary-period action is expressed in coordinate-label order."""
    from dzack_research.preamble.all import Lattices
    from sage.all import QQ as SageQQ
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    from sage_indefinite_port.indefinite.edgewalk import (
        _isometry_row_matrix_in_coordinate_order,
    )

    lattice_ = Lattices(SageZZ)([[2, 1], [1, 2]], module_generators=("right", "left"))
    labels = tuple(lattice_.module_generating_set())
    row_matrix = matrix(SageZZ, [[-1, 1], [0, 1]])
    images = {
        label: lattice_.linear_combination(
            {target_label: int(row_matrix[row, column]) for column, target_label in enumerate(labels) if row_matrix[row, column]}
        )
        for row, label in enumerate(labels)
    }
    isometry = lattice_.O()(images)

    assert _isometry_row_matrix_in_coordinate_order(lattice_, isometry) == row_matrix.change_ring(SageQQ)

    permuted_labels = ("negative", "positive")
    permuted_gram = [[-1, 0], [0, 2]]
    permuted = Lattices(SageZZ)(permuted_gram, module_generators=permuted_labels)
    assert tuple(permuted.module_generating_set()) == permuted_labels
    permuted_row_matrix = matrix(SageZZ, [[1, 0], [0, -1]])
    permuted_images = {
        label: permuted.linear_combination(
            {
                target_label: int(permuted_row_matrix[row, column])
                for column, target_label in enumerate(permuted_labels)
                if permuted_row_matrix[row, column]
            }
        )
        for row, label in enumerate(permuted_labels)
    }
    permuted_isometry = permuted.O()(permuted_images)
    assert _isometry_row_matrix_in_coordinate_order(permuted, permuted_isometry) == permuted_row_matrix.change_ring(
        SageQQ
    )


def test_binary_cycle_coordinate_order_and_period_certificate() -> None:
    """The rank-two period is an isometry in the owned lattice's generator order."""
    require(EDGEWALK)
    for gram, bound, expected_representatives, expected_period in (
        ([[1, 1], [1, -1]], 1, ((1, 2),), ((1, 2), (2, 5))),
        ([[2, 0], [0, -1]], 2, ((1, 1), (3, 4)), ((3, 4), (2, 3))),
    ):
        lattice_ = lattice(gram)
        cycle = _OrientedBinaryCycle(lattice_)
        first, second = tuple(lattice_.module_generator(label) for label in cycle.labels)
        assert cycle.form == (int(first.q()), int(lattice_.b(first, second)), int(second.q()))

        for coordinates_ in ((1, 0), (0, 1), (2, -3), (-5, 2)):
            point = cycle.point(coordinates_)
            assert cycle.coordinates(point) == coordinates_
            assert cycle.pairing(coordinates_, coordinates_) == int(point.q())

        result = cycle.cycle(bound, cycle.point((1, 0)))
        assert result is not None
        period, representatives = result
        assert tuple(cycle.coordinates(point) for point in representatives) == expected_representatives
        assert all(0 < int(point.q()) <= bound and int(point.content()) == 1 for point in representatives)
        assert tuple(cycle.coordinates(period(cycle.point((1, 0))))) == expected_period[0]
        assert tuple(cycle.coordinates(period(cycle.point((0, 1))))) == expected_period[1]
        from sage_indefinite_port.indefinite.edgewalk import _isometry_row_matrix_in_coordinate_order

        assert _isometry_row_matrix_in_coordinate_order(lattice_, period).det() == 1
        basis = tuple(lattice_.module_generators())
        assert all(lattice_.b(period(left), period(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_binary_cycle_period_respects_permuted_coordinate_basis() -> None:
    """The oriented binary cycle returns its period in the declared basis order."""
    from dzack_research.preamble.all import ZZ as PreambleZZ
    from dzack_research.preamble.all import Lattices

    from sage_indefinite_port.indefinite.edgewalk_rank2 import _OrientedBinaryCycle

    require(EDGEWALK)
    labels = ("negative", "positive")
    gram = [[-1, 0], [0, 2]]
    lattice_ = Lattices(PreambleZZ)(gram, module_generators=labels)
    assert tuple(lattice_.module_generating_set()) == labels
    cycle = _OrientedBinaryCycle(lattice_)
    result = cycle.cycle(2, cycle.point((0, 1)))

    default_lattice = Lattices(PreambleZZ)([[2, 0], [0, -1]])
    default_cycle = _OrientedBinaryCycle(default_lattice)
    default_result = default_cycle.cycle(2, default_cycle.point((1, 0)))
    assert result is not None and default_result is not None
    period, representatives = result
    _default_period, default_representatives = default_result

    # The custom basis swaps the two default coordinates. Compare the complete
    # representative sets in the same frame, while allowing the oriented
    # cycle's canonical companion to choose the inverse traversal direction.
    assert {
        (coordinates[1], coordinates[0])
        for coordinates in (cycle.coordinates(point) for point in representatives)
    } == {default_cycle.coordinates(point) for point in default_representatives}
    assert period.determinant() == 1
    assert all(
        lattice_.b(period(left), period(right)) == lattice_.b(left, right)
        for left in tuple(lattice_.module_generators())
        for right in tuple(lattice_.module_generators())
    )
    for coordinates in ((1, 0), (0, 1), (2, -3), (-5, 2)):
        assert cycle.coordinates(cycle.point(coordinates)) == coordinates


def test_edgewalk_branch_representative_sample(monkeypatch) -> None:
    """Exercise cusp, finite-corner, next-wall and oriented binary-cycle routes.

    Frozen IDs ``allcock_rk3_4`` and ``allcock_rk3_2`` cover the isotropic
    cusp and rationally anisotropic branches; ``U+<2>`` is the independently
    recorded three-wall example with both finite and ideal corners, and
    ``U+<46>`` covers the nonreflective stop after three vertex orbits. This
    sample does not claim direct coverage of the rank-two public ``binary_walk``
    entry point; the frozen full suite remains unchanged and still runs every
    recorded case.
    """
    require(EDGEWALK)
    observed = {}

    def instrument(cls, name):
        original = getattr(cls, name)

        def wrapped(self, *args, **kwargs):
            result = original(self, *args, **kwargs)
            observed[f"{cls.__name__}.{name}"] = True
            return result

        monkeypatch.setattr(cls, name, wrapped)

    for name in (
        "initial_isotropic_vertex", "initial_vertex", "ordinary_vertex_roots",
        "cusp_root", "cusp_roots", "extension_roots", "vertex_from_root",
        "isotropic_next_vertex", "finite_extensions",
        "edge_step", "edge_directions", "binary_walk",
    ):
        instrument(_Edgewalk, name)
    for name in ("cycle", "promised_step", "canonical_companion", "shorter_pair", "reduced_start_pair"):
        instrument(_OrientedBinaryCycle, name)
    original_anisotropic_cycle = edgewalk_module.anisotropic_cycle

    def wrapped_anisotropic_cycle(*args, **kwargs):
        result = original_anisotropic_cycle(*args, **kwargs)
        observed["anisotropic_cycle"] = True
        return result

    monkeypatch.setattr(edgewalk_module, "anisotropic_cycle", wrapped_anisotropic_cycle)

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

    # Keep one independent polyhedral_common root-system oracle in this branch
    # sample rather than relying only on the full-file parameterizers. The
    # Julia-backed OSCAR and E10 comparisons are isolated below because the
    # local Sage runtime may lack Julia.
    polyhedral_case = load_root_systems()[0]
    _assert_same_root_system(polyhedral_case["gram"], polyhedral_case["roots"])
    three_wall = _hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 2]]).allcock_edgewalk()
    assert len(three_wall.simple_roots()) == 3
    vertex_squares = tuple(
        int(ray.q())
        for vertex in three_wall.vertices()
        for ray in vertex.primitive_rays()
    )
    assert any(square == 0 for square in vertex_squares)
    assert any(square < 0 for square in vertex_squares)
    # Provider revisions may discover a finite extension before the whole
    # walk reaches its isotropic fallback. Exercise that exact routine on the
    # recorded U+<2> initial edge independently of whether edge_step needs it.
    fallback_walk = _Edgewalk(_hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 2]]))
    fallback_start = fallback_walk.initial_isotropic_vertex() or fallback_walk.initial_vertex()
    assert fallback_start is not None
    fallback_roots, fallback_discarded = fallback_walk.edge_directions(fallback_start)[0]
    fallback_vertex = fallback_walk.isotropic_next_vertex(
        fallback_roots, fallback_discarded, fallback_start.generator
    )
    assert fallback_vertex.generator.q() == 0
    nonreflective = _hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 46]]).allcock_edgewalk()
    assert not nonreflective.is_reflective()
    assert len(nonreflective.vertices()) == 3
    # The compact hexagon has a rootless initial perpendicular and exercises
    # the increasing-prefix starting-corner search rather than the direct
    # full-rank finite-root corner route.
    compact_hexagon = _hyperbolic([[-1, 0, 0], [0, 3, 1], [0, 1, 7]]).allcock_edgewalk()
    assert compact_hexagon.is_reflective()
    assert len(compact_hexagon.simple_roots()) == 6
    labelled_gram = [[7, 0, 1], [0, -1, 0], [1, 0, 3]]
    labelled_labels = ("right", "negative", "left")
    labelled_hexagon_lattice = Lattices(PreambleZZ)(labelled_gram, module_generators=labelled_labels)
    assert tuple(labelled_hexagon_lattice.module_generating_set()) == labelled_labels
    labelled_hexagon = HyperbolicLattices(PreambleZZ)(labelled_hexagon_lattice).allcock_edgewalk()
    assert labelled_hexagon.is_reflective()
    labelled_roots = [
        tuple(int(root.to_vector()(label)) for label in labelled_labels) for root in labelled_hexagon.simple_roots()
    ]
    gram = [[-1, 0, 0], [0, 3, 1], [0, 1, 7]]
    assert _canonical_root_graph(_root_gram(labelled_gram, labelled_roots)) == _canonical_root_graph(
        _root_gram(gram, [coordinates(root) for root in compact_hexagon.simple_roots()])
    )
    binary = _hyperbolic([[1, 0], [0, -2]]).allcock_edgewalk()
    assert binary.is_reflective()
    assert len(binary.simple_roots()) == 2
    binary_isotropic = _hyperbolic([[0, 1], [1, 0]]).allcock_edgewalk()
    assert not binary_isotropic.is_reflective()
    assert len(binary_isotropic.simple_roots()) == 0
    assert len(binary_isotropic.vertices()) == 0

    required = {
        "_Edgewalk.initial_isotropic_vertex", "_Edgewalk.initial_vertex",
        "_Edgewalk.ordinary_vertex_roots", "_Edgewalk.cusp_root", "_Edgewalk.cusp_roots",
        "_Edgewalk.extension_roots", "_Edgewalk.vertex_from_root",
        "_Edgewalk.isotropic_next_vertex", "_Edgewalk.finite_extensions",
        "_Edgewalk.edge_step", "_Edgewalk.edge_directions", "_Edgewalk.binary_walk",
        "anisotropic_cycle",
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
