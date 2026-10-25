"""Targeted source-contract regressions for the pinned Lorentzian CVP port.

These probes are deliberately separate from the frozen acceptance corpus.
"""

from sage_indefinite_port.indefinite.lorentzian_cells import (
    _dual_lll_enumeration_key,
    _pinned_integer_covector_kernel,
    _pinned_positive_gcd_dot,
    _pinned_sublattice_basis_reduction,
)

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

    # The kernel supplied by NullspaceIntVect is rectangular in ambient
    # coordinates; exercise the same reducer on a rank-three sublattice of Z^4.
    rectangular = matrix(ZZ, [[10, 1, 0, 2], [9, 1, 1, 0], [1, 0, 8, 1]])
    rectangular_reduced = _pinned_sublattice_basis_reduction(rectangular)
    assert rectangular_reduced == matrix(ZZ, [[-1, 0, 1, -2], [2, 0, 7, 3], [8, 1, 2, -2]])
    rectangular_change = rectangular_reduced * rectangular.transpose() * (rectangular * rectangular.transpose()).inverse()
    assert all(entry.denominator() == 1 for entry in rectangular_change.list())
    assert abs(rectangular_change.det()) == 1
    assert rectangular_change * rectangular == rectangular_reduced


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

