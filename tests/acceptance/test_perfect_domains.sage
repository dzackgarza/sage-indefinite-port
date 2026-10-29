"""Lorentzian perfect-domain traversal against polyhedral_common's recorded cell counts.

28B_LorentzianPerfStabEqui records, for 40 Lorentzian lattices of rank 3 to 6, the number
of perfect-domain orbits in the "total" and "isotropic" modes. The total mode is reached
through the preamble's reduction complex, whose generators must generate O(L).
"""

from dataclasses import replace

import pytest

from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    ZERO_SQUARES,
    LorentzianPerfectLocalBackend,
    _search_initial_vectors,
    perfect_complex,
    perfect_domain_traversal,
)
from tests.acceptance.consumer import PERFECT_DOMAINS, lattice, require
from tests.fixtures.oracle_fixtures import load_lorentzian_perfect_domains


def test_non_diagonal_positive_shell_uses_the_pinned_cvp_center() -> None:
    """The first shell in a non-diagonal Lorentzian form matches exact enumeration."""
    from sage.all import ZZ as SageZZ
    from sage.all import matrix, vector

    require(PERFECT_DOMAINS)
    lattice_ = lattice([[5, 1, 0], [1, -1, 0], [0, 0, -1]])
    rows = _search_initial_vectors(lattice_, (1, 0, 0), NONNEGATIVE_SQUARES)

    assert set(rows) == {(1, -1, -1), (1, -1, 0), (1, -1, 1)}
    gram = matrix(SageZZ, [[5, 1, 0], [1, -1, 0], [0, 0, -1]])
    for row in rows:
        point = vector(SageZZ, row)
        assert (point * gram * point.column())[0] in {1, 2}
        assert (point * gram * vector(SageZZ, (1, 0, 0)).column())[0] == 4


def test_non_diagonal_positive_shell_respects_explicit_coordinate_labels() -> None:
    """The pinned positive shell is exact when the lattice uses string coordinate labels."""
    from dzack_research.preamble.all import ZZ as PreambleZZ
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ
    from sage.all import matrix, vector

    require(PERFECT_DOMAINS)
    gram_data = [[5, 1, 0], [1, -1, 0], [0, 0, -1]]
    labelled = Lattices(PreambleZZ)(gram_data, module_generators=("time", "space", "auxiliary"))
    rows = _search_initial_vectors(labelled, (1, 0, 0), NONNEGATIVE_SQUARES)

    assert set(rows) == {(1, -1, -1), (1, -1, 0), (1, -1, 1)}
    gram = matrix(SageZZ, gram_data)
    for row in rows:
        point = vector(SageZZ, row)
        assert (point * gram * point.column())[0] in {1, 2}
        assert (point * gram * vector(SageZZ, (1, 0, 0)).column())[0] == 4


def test_non_diagonal_isotropic_shell_uses_the_exact_sphere_center() -> None:
    """The exact-sphere branch finds every isotropic vector in the first shell."""
    from sage.all import QQ, matrix, vector
    from sage.all import ZZ as SageZZ

    from sage_indefinite_port.indefinite.lorentzian_cells import (
        _PositiveVectorEnumerator,
    )

    require(PERFECT_DOMAINS)
    gram = matrix(SageZZ, [[1, 1, 0], [1, 0, 0], [0, 0, -1]])
    lattice_ = lattice([[1, 1, 0], [1, 0, 0], [0, 0, -1]])
    enumerator = _PositiveVectorEnumerator(lattice_, (QQ(1), QQ(0), QQ(0)), QQ(1), ZERO_SQUARES)
    assert enumerator.base_target_engine_coordinates != (QQ(0), QQ(0))
    first = enumerator.first_shell()
    assert first is not None and first[0] == 1
    rows = first[1]

    assert set(rows) == {(0, 1, 0), (1, 0, -1), (1, 0, 1), (2, -1, 0)}
    for row in rows:
        point = vector(SageZZ, row)
        assert (point * gram * point.column())[0] == 0
        assert (point * gram * vector(SageZZ, (1, 0, 0)).column())[0] == 1


def test_positive_shell_scales_rational_direction_pairing_bound() -> None:
    """The first shell is unchanged when a rational direction bound is scaled consistently."""
    from dzack_research.preamble.all import Lattices, ZZ as PreambleZZ
    from sage.all import QQ
    from sage_indefinite_port.indefinite.lorentzian_cells import _PositiveVectorEnumerator

    require(PERFECT_DOMAINS)
    lattice_ = Lattices(PreambleZZ)([[5, 1, 0], [1, -1, 0], [0, 0, -1]])
    integral = _PositiveVectorEnumerator(lattice_, (QQ(1), QQ(0), QQ(0)), QQ(5), NONNEGATIVE_SQUARES)
    rational = _PositiveVectorEnumerator(lattice_, (QQ(1) / 2, QQ(0), QQ(0)), QQ(5) / 2, NONNEGATIVE_SQUARES)

    assert rational.scale == 2
    assert rational.scaled_max == integral.scaled_max == 5
    assert rational.first_shell() == integral.first_shell() == (4, ((1, -1, -1), (1, -1, 0), (1, -1, 1)))


def test_non_diagonal_isotropic_shell_respects_explicit_coordinate_labels() -> None:
    """The exact-sphere CVP branch preserves all isotropic rows in a labelled basis."""
    from dzack_research.preamble.all import ZZ as PreambleZZ
    from dzack_research.preamble.all import Lattices
    from sage.all import QQ, matrix, vector
    from sage.all import ZZ as SageZZ

    from sage_indefinite_port.indefinite.lorentzian_cells import (
        _PositiveVectorEnumerator,
    )

    require(PERFECT_DOMAINS)
    gram_data = [[1, 1, 0], [1, 0, 0], [0, 0, -1]]
    lattice_ = Lattices(PreambleZZ)(gram_data, module_generators=("e", "f", "radical"))
    enumerator = _PositiveVectorEnumerator(lattice_, (QQ(1), QQ(0), QQ(0)), QQ(1), ZERO_SQUARES)
    first = enumerator.first_shell()

    assert first is not None and first[0] == 1
    rows = first[1]
    assert set(rows) == {(0, 1, 0), (1, 0, -1), (1, 0, 1), (2, -1, 0)}
    gram = matrix(SageZZ, gram_data)
    for row in rows:
        point = vector(SageZZ, row)
        assert (point * gram * point.column())[0] == 0
        assert (point * gram * vector(SageZZ, (1, 0, 0)).column())[0] == 1


@pytest.mark.parametrize("case", load_lorentzian_perfect_domains(), ids=lambda case: case["id"])
def test_perfect_domain_orbit_counts_match_polyhedral_common(case) -> None:
    require(PERFECT_DOMAINS)

    traversal = lattice(case["gram"]).lorentzian_reduction_complex()

    assert traversal.is_complete()
    assert len(traversal.cells()) == case["total_count"]
    assert traversal.generates_orthogonal_group()


@pytest.mark.parametrize("case", load_lorentzian_perfect_domains(), ids=lambda case: case["id"])
def test_isotropic_mode_counts_match_polyhedral_common(case) -> None:
    """The "isotropic" mode is a reference-implementation option with no preamble consumer,
    so the case calls the port's entry point directly, on the preamble lattice."""
    require(PERFECT_DOMAINS)

    assert len(perfect_domain_traversal(lattice(case["gram"]), "isotropic").cells()) == case["isotropic_count"]


def test_lazy_neighbor_cones_retain_exact_transformed_rays() -> None:
    """Materializing a deferred side-pairing cone checks its exact ray certificate."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    traversal = lattice(case["gram"]).lorentzian_reduction_complex()
    adjacencies = tuple(traversal.adjacencies())
    assert adjacencies
    for adjacency in adjacencies:
        neighbor = adjacency.neighbor()
        common_face = adjacency.common_face()
        assert neighbor.ambient_lattice() is traversal.lattice()
        assert common_face.ambient_lattice() is traversal.lattice()
        assert int(common_face.dimension()) == int(traversal.lattice().rank()) - 1
        assert common_face.is_face_of(adjacency.source())
        assert common_face.is_face_of(neighbor)
        assert adjacency.target_representative().transport(adjacency.target_to_neighbor()).is_equal_to(neighbor)


def test_traversal_retains_only_adjacency_facet_orbits() -> None:
    """The traversal retains Schreier words only for facets with side pairings."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    complex_ = perfect_complex(lattice(case["gram"]), NONNEGATIVE_SQUARES)
    adjacencies = tuple(complex_.adjacencies())
    assert adjacencies
    for source_id, stored_orbits in complex_._traversal_facet_orbits.items():
        expected_facets = {adjacency.facet for adjacency in adjacencies if id(adjacency.source) == source_id}
        assert set(stored_orbits) == expected_facets


def test_facet_orbit_permutation_view_materializes_exact_isometries() -> None:
    """The public owned-isometry view agrees with the internal Schreier permutations."""
    require(PERFECT_DOMAINS)

    from dzack_research.preamble.all import Lattices, ZZ as PreambleZZ

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    labels = tuple(f"basis_{index}" for index in range(len(case["gram"])))
    labelled_lattice = Lattices(PreambleZZ)(case["gram"], module_generators=labels)
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(labelled_lattice, NONNEGATIVE_SQUARES)
    permutation_orbits = backend.facet_orbit_permutations(cell)
    isometry_orbits = backend.facet_orbits(cell)

    assert tuple(tuple(facet for facet, _ in orbit) for orbit in permutation_orbits) == tuple(
        tuple(facet for facet, _ in orbit) for orbit in isometry_orbits
    )
    for permutation_orbit, isometry_orbit in zip(permutation_orbits, isometry_orbits, strict=True):
        for (facet, permutation), (isometry_facet, isometry) in zip(permutation_orbit, isometry_orbit, strict=True):
            assert facet == isometry_facet
            for position, image_position in enumerate(permutation):
                assert isometry(cell.vector_configuration[position]) == cell.vector_configuration[image_position]


def test_large_cell_cone_uses_exact_normaliz_h_presentation(monkeypatch) -> None:
    """The large-configuration backend preserves the owned exact cone interface."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(lattice(case["gram"]), NONNEGATIVE_SQUARES)
    # Repeating the same rays keeps the cone unchanged while forcing the
    # large-configuration Normaliz path without a costly high-rank fixture.
    large_cell = replace(cell, rows=cell.rows * 4_000)
    assert len(large_cell.rows) >= 10_000
    small_cone = backend.cell_cone(cell)
    # Exercise the direct Normaliz accessors before constructing the owned
    # large cone: these are the hot paths whose purpose is to avoid the PPL
    # reconstruction, not merely to reuse an engine cached by cell_cone().
    with monkeypatch.context() as patch:
        patch.setattr(
            backend,
            "cell_cone",
            lambda _cell: pytest.fail("a direct large-cell Normaliz accessor rebuilt the owned cone"),
        )
        large_covectors = backend.cell_halfspace_covectors(large_cell)
        assert tuple(backend._iter_cell_halfspace_covectors(large_cell)) == large_covectors
        large_primitive_rays = backend.cell_primitive_rays(large_cell)
        large_null_rays = backend.cell_primitive_rays(large_cell, square=0)
    large_cone = backend.cell_cone(large_cell)
    assert backend._normaliz_cell_polyhedron(large_cell)._normaliz_cone is backend._normaliz_cell_cone(large_cell)
    labels = tuple(cell.lattice.module_generating_set())

    def covector_rows(cone):
        return {
            tuple(int(covector.to_vector()(label)) for label in labels)
            for covector in cone.halfspace_covectors()
        }

    def incidence_rows(cone, ray_rows):
        normals = tuple(cone.halfspace_covectors())
        return {
            frozenset(
                position
                for position, ray in enumerate(ray_rows)
                if all(
                    sum(
                        covector.to_vector()(label) * entry
                        for label, entry in zip(labels, ray, strict=True)
                    ) == 0
                    for covector in (normal,)
                )
            )
            for normal in normals
        }

    assert covector_rows(large_cone) == covector_rows(small_cone)
    assert {
        tuple(int(covector.to_vector()(label)) for label in labels)
        for covector in large_covectors
    } == covector_rows(small_cone)
    assert large_cone.dimension() == small_cone.dimension()
    assert set(large_cone.primitive_rays()) == set(small_cone.primitive_rays())
    assert set(large_primitive_rays) == set(small_cone.primitive_rays())
    assert set(large_null_rays) == {ray for ray in small_cone.primitive_rays() if ray.q() == 0}
    assert len(large_cone.faces(int(large_cone.dimension()) - 1)) == len(small_cone.faces(int(small_cone.dimension()) - 1))
    large_incidence = incidence_rows(large_cone, large_cell.rows)
    assert {
        frozenset(position % len(cell.rows) for position in incidence)
        for incidence in large_incidence
    } == incidence_rows(small_cone, cell.rows)
