"""Lorentzian perfect-domain traversal against polyhedral_common's recorded cell counts.

28B_LorentzianPerfStabEqui records, for 40 Lorentzian lattices of rank 3 to 6, the number
of perfect-domain orbits in the "total" and "isotropic" modes. The total mode is reached
through the preamble's reduction complex, whose generators must generate O(L).
"""

import pytest

from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectLocalBackend,
    NONNEGATIVE_SQUARES,
    _dual_lll_enumeration_key,
    _pinned_integer_covector_kernel,
    _pinned_positive_gcd_dot,
    _pinned_sublattice_basis_reduction,
    perfect_complex,
    perfect_domain_traversal,
)
from tests.acceptance.consumer import PERFECT_DOMAINS, lattice, require
from tests.fixtures.oracle_fixtures import load_lorentzian_perfect_domains


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


def test_gap_lll_transformation_uses_row_congruence_orientation() -> None:
    """The GAP provider's transformation acts by P G P^T, not P^T G P."""
    from sage.all import QQ, ZZ, matrix, libgap

    gram = matrix(ZZ, [[100, 31, 17], [31, 11, 6], [17, 6, 4]])
    record = libgap.LLLReducedGramMat(
        [[int(gram[i, j]) for j in range(gram.ncols())] for i in range(gram.nrows())],
        libgap(99) / libgap(100),
    )
    transform = matrix(ZZ, [[ZZ(int(entry)) for entry in row] for row in record["transformation"]])
    reduced = matrix(ZZ, [[ZZ(int(entry)) for entry in row] for row in record["remainder"]])

    assert transform == matrix(ZZ, [[0, -1, 1], [0, -1, 2], [1, -3, 0]])
    assert abs(transform.det()) == 1
    assert transform * gram * transform.transpose() == reduced
    assert transform.transpose() * gram * transform != reduced

    dual_record = libgap.LLLReducedGramMat(
        [[int(gram.adjugate()[i, j]) for j in range(gram.ncols())] for i in range(gram.nrows())],
        libgap(99) / libgap(100),
    )
    dual_transform = matrix(ZZ, [[ZZ(int(entry)) for entry in row] for row in dual_record["transformation"]])
    dual_reduced = matrix(ZZ, [[ZZ(int(entry)) for entry in row] for row in dual_record["remainder"]])
    assert dual_transform == matrix(ZZ, [[1, 0, 0], [3, 1, 1], [6, 2, 1]])
    assert dual_transform * gram.adjugate() * dual_transform.transpose() == dual_reduced
    enumeration_transform = dual_transform.inverse().transpose()
    assert enumeration_transform * gram * enumeration_transform.transpose() == gram.det() * dual_reduced.change_ring(QQ).inverse()
    # Exercise the production ordering key against the actual pinned dual-LLL
    # transform, not only a separately chosen nonsymmetric unimodular matrix.
    candidates = (vector(ZZ, (4, 2, 1)), vector(ZZ, (3, 0, 2)), vector(ZZ, (1, 5, 0)), vector(ZZ, (0, 1, 3)))
    assert tuple(sorted(candidates, key=lambda x: _dual_lll_enumeration_key(tuple(x), dual_transform))) == (
        candidates[3], candidates[2], candidates[1], candidates[0]
    )
    assert tuple(sorted(candidates, key=lambda x: _dual_lll_enumeration_key(tuple(x), dual_transform.transpose()))) != (
        candidates[3], candidates[2], candidates[1], candidates[0]
    )


def test_dual_lll_candidates_sort_in_reduced_column_coordinates() -> None:
    """CVPSolver maps reduced columns y to original columns x=T^{-1}y."""
    from sage.all import ZZ, matrix, vector

    transform = matrix(ZZ, [[1, 2, 0], [0, 1, 1], [0, 0, 1]])
    reduced_columns = (vector(ZZ, (3, 2, 2)), vector(ZZ, (11, 5, 0)), vector(ZZ, (8, 3, 1)))
    # CVPSolver reconstructs x=P^{-1}y. Recover the original columns from
    # deliberately out-of-order recursive coordinates and check the inverse.
    original_columns = tuple(transform.inverse() * reduced.column() for reduced in reduced_columns)
    assert tuple(tuple(column[i, 0] for i in range(3)) for column in original_columns) == (
        (3, 0, 2), (1, 5, 0), (4, 2, 1)
    )
    assert tuple(tuple((transform * column)[i, 0] for i in range(3)) for column in original_columns) == tuple(
        tuple(column) for column in reduced_columns
    )
    expected_order = (original_columns[1], original_columns[2], original_columns[0])
    assert tuple(sorted(original_columns, key=lambda column: _dual_lll_enumeration_key(tuple(column[i, 0] for i in range(3)), transform))) == expected_order
    # The transpose gives a different order on this nonsymmetric transform;
    # this guards against accidentally using Q^{-1}=P^T instead of P.
    assert tuple(sorted(original_columns, key=lambda column: _dual_lll_enumeration_key(tuple(column[i, 0] for i in range(3)), transform.transpose()))) != expected_order


def test_pinned_positive_gcd_dot_matches_exact_coefficients() -> None:
    """The source GCD-dot recurrence fixes the affine closest-vector center."""
    cases = (
        ((6, 10, 15), (1, (-14, 7, 1))),
        ((-12, 18, 0, -30, 42), (6, (1, 1, 0, 0, 0))),
        ((0, 0, 9), (9, (0, 0, 1))),
        ((-5, 0, 0), (5, (-1, 0, 0))),
        ((0, 0, 0), (0, (0, 0, 0))),
        ((0,), (0, (-1,))),
    )
    for values, expected in cases:
        divisor, coefficients = _pinned_positive_gcd_dot(values)
        assert (divisor, coefficients) == expected
        assert sum(a * b for a, b in zip(values, coefficients, strict=True)) == divisor


def test_pinned_sublattice_basis_reduction_preserves_row_lattice() -> None:
    """Alternating direct/dual LLL decreases L1 without changing the lattice."""
    from sage.all import ZZ, matrix

    basis = matrix(ZZ, [[10, 1, 0], [9, 1, 1], [1, 0, 8]])
    reduced = _pinned_sublattice_basis_reduction(basis)
    assert reduced.base_ring() is ZZ
    assert reduced == matrix(ZZ, [[-1, 0, 1], [-1, -1, 0], [3, -3, 3]])
    change = reduced * basis.inverse()
    assert all(entry.denominator() == 1 for entry in change.list())
    assert abs(change.det()) == 1
    assert sum(abs(entry) for entry in reduced.list()) <= sum(abs(entry) for entry in basis.list())
    assert change * basis == reduced


def test_pinned_integer_covector_kernel_matches_exact_contract() -> None:
    """The source-ordered integral kernel is annihilated and saturated."""
    from sage.all import ZZ, matrix, vector

    cases = (
        ((-7,), ()),
        ((7,), ()),
        ((2, 3), ((-3, 2),)),
        ((0, 6, 9), ((1, 0, 0), (0, -3, 2))),
        ((6, 10, 15), ((-5, 3, 0), (-30, 15, 2))),
        ((0, 0, 1, 2), ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, -2, 1))),
        ((7, 0, 0), ((0, 1, 0), (0, 0, 1))),
    )
    for covector, expected_rows in cases:
        kernel = _pinned_integer_covector_kernel(covector)
        assert tuple(tuple(int(entry) for entry in row) for row in kernel.rows()) == expected_rows
        assert kernel.nrows() == len(covector) - 1
        assert kernel * vector(ZZ, covector).column() == 0
        assert all(factor == 1 for factor in kernel.smith_form()[0].diagonal())

    zero_covector = (0, 0, 0)
    assert _pinned_integer_covector_kernel(zero_covector) == matrix.identity(ZZ, 3)
    assert _pinned_integer_covector_kernel((0,)) == matrix.identity(ZZ, 1)


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


def test_facet_orbit_permutation_view_materializes_exact_isometries() -> None:
    """The public owned-isometry view agrees with the internal Schreier permutations."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(lattice(case["gram"]), NONNEGATIVE_SQUARES)
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
