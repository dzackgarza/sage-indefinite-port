"""Targeted source-contract regressions for the pinned Lorentzian CVP port.

These probes are deliberately separate from the frozen acceptance corpus.
"""

from sage_indefinite_port.indefinite.lorentzian_cells import (
    _PositiveVectorEnumerator,
    _dual_lll_enumeration_key,
    _integral_direction_pairing_bound,
    _next_integral_pairing_bound,
    _pinned_integer_covector_kernel,
    _pinned_positive_gcd_dot,
    _pinned_sublattice_basis_reduction,
    _primitive_integral_direction,
    _source_mid_value,
    _source_some_negative_vector,
    _pinned_nearest_integer,
    _source_permutation_sign_matrix,
    _source_elementary_perturbation,
    _source_random_matrix_perturbation,
    _source_graver_k_basis,
    _source_direct_graver_improvement,
    _source_should_reset_positive_perturbation,
    _source_positive_search_best,
    _source_indefinite_short_positive_vector,
    _pinned_nullspace_matrix,
    _pinned_diagonalize_symmetric_matrix,
    _integral_positive_vector_allmeth,
    _source_certify_positive_integral_vector,
    _source_integral_vector_allmeth,
    _source_validate_integral_reduction,
    _source_remove_fraction_vector,
    _source_integral_vector_allmeth_v2,
    _source_reconstruct_integral_vector,
    _pinned_positive_eigen_directions,
    _pinned_positive_diagonal_directions,
    _pinned_approximation_iterator,
    _advance_pinned_approximation_iterator,
    _pinned_independent_column_indices,
    _outside_ray_directions,
    _quadratic_descent_prefix,
    _positive_ldl_factorization,
    _fixed_distance_coordinate_solutions,
    _fixed_distance_dual_lll_solutions,
    _validate_fixed_distance_ldl_factorization,
    _validate_fixed_distance_basis_inverse,
    _upper_bound,
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
    assert enumeration_transform * dual_transform.transpose() == matrix.identity(ZZ, 3)
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
    # The pinned recursion orders y=P*x+s. A target-dependent affine shift s
    # is constant across one result list, so it must not alter the order.
    affine_shift = vector(ZZ, (7, -4, 2))
    shifted_keys = lambda column: tuple(
        (transform * column)[i, 0] + affine_shift[i] for i in (2, 1, 0)
    )
    assert tuple(sorted(original_columns, key=shifted_keys)) == expected_order


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


def test_outside_ray_span_normalizes_each_source_direction_independently() -> None:
    """GetOneOutsideRay receives individually primitive direction rows."""
    from sage.all import QQ, matrix, vector

    gram = matrix(QQ, [[1, 0], [0, 2]])
    affine_nullspace = matrix(QQ, [[0, 2, 0], [0, 0, 3]])
    directions = _outside_ray_directions(affine_nullspace, gram)
    assert directions == matrix(QQ, [[1, 0], [0, 1]])
    assert all(any(entry for entry in row) for row in directions.rows())


def test_initial_cell_nullspace_uses_pinned_gap_row_order() -> None:
    """The GetOnePerfect affine nullspace is the pinned GAP row basis."""
    from sage.all import QQ, matrix

    expanded = matrix(QQ, [[1, 1, 0], [1, 0, 1]])
    assert _pinned_nullspace_matrix(expanded.transpose()) == matrix(QQ, [[-1, 1, 1]])


def test_pinned_column_reduction_matches_greedy_rank_increase() -> None:
    """ColumnReduction scans source columns and retains each rank increase."""
    from random import Random

    from sage.all import QQ, matrix

    random = Random(1618033)
    for rows in range(1, 7):
        for columns in range(1, 7):
            for _ in range(8):
                for denominator in (1, 2, 3):
                    value = matrix(
                        QQ,
                        rows,
                        columns,
                        lambda i, j: QQ(random.randint(-3, 3)) / denominator,
                    )
                    selected = _pinned_independent_column_indices(value)
                    assert len(selected) == value.rank()
                    for index, column in enumerate(selected):
                        assert value.matrix_from_columns(selected[: index + 1]).rank() == index + 1
                        assert all(previous < column for previous in selected[:index])
                    assert selected == tuple(value.echelon_form().pivots())


def test_pinned_positive_diagonal_directions_follow_source_scan_and_l1_norm() -> None:
    """GetPositiveDirections_diag keeps positive rows in diagonal order."""
    from sage.all import QQ, matrix

    transform = matrix(QQ, [[1, 2, 0, 0], [0, 1, 1, 0], [2, 0, -2, 1], [1, 0, 0, 0]])
    diagonal = matrix(QQ, [[2, 0, 0, 0], [0, 0, 0, 0], [0, 0, -3, 0], [0, 0, 0, 5]])
    assert _pinned_positive_diagonal_directions(transform, diagonal) == (
        (QQ(1) / 3, QQ(2) / 3, QQ.zero(), QQ.zero()),
        (QQ(1), QQ.zero(), QQ.zero(), QQ.zero()),
    )
    try:
        _pinned_positive_diagonal_directions(transform, matrix(QQ, [[1, 1], [1, 1]]))
    except ValueError:
        pass
    else:
        raise AssertionError("nondiagonal Gram matrix accepted")


def test_pinned_positive_eigen_directions_preserve_source_row_orientation() -> None:
    """The Eigen consumer reads rows of its eigenvector matrix, not columns."""
    from sage.all import RDF, matrix

    gram = matrix(RDF, [[2.0, 1.0], [1.0, 2.0]])
    data = gram.eigenvectors_right()
    eigenvalues = []
    columns = []
    for eigenvalue, vectors, multiplicity in data:
        assert len(vectors) == multiplicity
        for eigenvector in vectors:
            eigenvalues.append(float(eigenvalue.real()))
            columns.append(tuple(float(entry.real()) for entry in eigenvector))
    expected = tuple(
        tuple(columns[column][row] for column in range(len(columns)))
        for row, eigenvalue in enumerate(eigenvalues)
        if eigenvalue > 0
    )
    directions = _pinned_positive_eigen_directions(gram)
    assert directions == expected
    assert all(len(direction) == gram.nrows() for direction in directions)
    assert any(direction != columns[index] for index, direction in enumerate(directions))


def test_pinned_positive_eigen_directions_rejects_invalid_double_gram() -> None:
    """The optional Eigen path requires finite symmetric input before solving."""
    from math import nan

    from sage.all import RDF, matrix

    for gram in (
        matrix(RDF, [[1.0, 2.0], [0.0, 1.0]]),
        matrix(RDF, [[nan, 0.0], [0.0, 1.0]]),
    ):
        try:
            _pinned_positive_eigen_directions(gram)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid double Gram matrix accepted: {gram}")


def test_pinned_approximation_iterator_rejects_nonfinite_float_directions() -> None:
    """Optional Eigen acceleration cannot poison iterator ordering with NaN/Inf."""
    from math import inf, nan

    for direction in ((nan, 1.0), (inf, 1.0), (1e-320, 1.0)):
        try:
            _pinned_approximation_iterator(direction)
        except ValueError:
            pass
        else:
            raise AssertionError(f"nonfinite ApproxIterator state accepted: {direction}")


def test_pinned_approximation_iterator_advances_tied_steps_in_source_order() -> None:
    """ApproxIterator increments every group tied at its minimum next-time."""
    from sage.all import QQ

    ordered, nexts, approximation = _pinned_approximation_iterator((QQ(1) / 2, QQ(1) / 3, -QQ(1) / 2))
    assert ordered == ((QQ(2), ((0, 1), (2, -1))), (QQ(3), ((1, 1),)))
    assert nexts == [QQ(1), QQ(3) / 2]
    assert _advance_pinned_approximation_iterator(ordered, nexts, approximation) == [1, 0, -1]
    assert nexts == [QQ(3), QQ(3) / 2]
    assert _advance_pinned_approximation_iterator(ordered, nexts, approximation) == [1, 1, -1]
    assert _advance_pinned_approximation_iterator(ordered, nexts, approximation) == [2, 1, -2]

    float_ordered, float_nexts, float_approximation = _pinned_approximation_iterator((1.0, -0.5))
    assert _advance_pinned_approximation_iterator(float_ordered, float_nexts, float_approximation) == [1, 0]
    assert _advance_pinned_approximation_iterator(float_ordered, float_nexts, float_approximation) == [1, -1]

    try:
        _pinned_approximation_iterator((QQ.zero(), QQ.zero()))
    except ValueError:
        pass
    else:
        raise AssertionError("zero direction unexpectedly created an ApproxIterator")


def test_source_positive_candidate_lift_uses_transposed_unimodular_basis() -> None:
    """A V2 candidate is mapped by B^T and keeps its exact quadratic value."""
    from sage.all import ZZ, matrix, vector

    gram = matrix(ZZ, [[2, 1], [1, -1]])
    basis = matrix(ZZ, [[1, 2], [0, 1]])
    reduced = basis * gram * basis.transpose()
    candidate = vector(ZZ, [1, 0])
    expected = basis.transpose() * candidate.column()
    result = _source_reconstruct_integral_vector(gram, basis, reduced, candidate)
    assert result == tuple(expected[i, 0] for i in range(expected.nrows()))
    assert (vector(ZZ, result) * gram * vector(ZZ, result).column())[0] == (candidate * reduced * candidate.column())[0] > 0
    assert _source_reconstruct_integral_vector(gram, basis, reduced, vector(ZZ, [0, 1])) is None

    invalid_reduced = reduced + matrix(ZZ, [[1, 0], [0, 0]])
    try:
        _source_reconstruct_integral_vector(gram, basis, invalid_reduced, candidate)
    except ArithmeticError:
        pass
    else:
        raise AssertionError("a false reduction norm certificate was accepted")


def test_definite_dual_lll_uses_inverse_of_sage_transform() -> None:
    """Dual lattice reduction uses Q=P^{-T} for source row coordinates P=U^T."""
    from sage.all import ZZ, identity_matrix, matrix
    from sage_indefinite_port.indefinite.indefinite_reduction import _definite_lll

    gram = identity_matrix(ZZ, 2)
    sage_transform = matrix(ZZ, [[1, 1], [0, 1]])
    dual, reduced = _definite_lll(gram, dual=True, lll_gram=lambda _form: sage_transform)
    assert dual == sage_transform.inverse()
    assert reduced == dual * gram * dual.transpose()
    direct, direct_reduced = _definite_lll(gram, dual=False, lll_gram=lambda _form: sage_transform)
    assert direct == sage_transform.transpose()
    assert direct_reduced == direct * gram * direct.transpose()


def test_simple_indefinite_reduction_certifies_each_l1_descent(monkeypatch) -> None:
    """The pinned shear search strictly decreases L1 and preserves congruence."""
    from sage.all import ZZ, matrix
    import sage_indefinite_port.indefinite.indefinite_reduction as reduction

    monkeypatch.setattr(reduction, "randrange", lambda _n: 0)
    gram = matrix(ZZ, [[1, 10, 0], [10, 1, 1], [0, 1, -2]])
    basis, reduced = reduction._simple_indefinite_reduction(gram)
    assert abs(basis.det()) == 1
    assert reduced == basis * gram * basis.transpose()
    assert sum(abs(x) for row in reduced.rows() for x in row) <= sum(
        abs(x) for row in gram.rows() for x in row
    )


def test_compute_indefinite_reduction_rejects_false_lll_certificate(monkeypatch) -> None:
    """The outer isotropic/restart loop validates each inner reduction result."""
    from sage.all import ZZ, identity_matrix, matrix
    import sage_indefinite_port.indefinite.indefinite_reduction as reduction

    gram = matrix(ZZ, [[0, 1], [1, 0]])
    identity = identity_matrix(ZZ, 2)
    monkeypatch.setattr(
        reduction,
        "_indefinite_lll",
        lambda current: (identity, current + identity, None),
    )
    try:
        reduction._compute_indefinite_reduction(gram)
    except ArithmeticError as error:
        assert "false Gram certificate" in str(error)
    else:
        raise AssertionError("a false inner LLL result entered the restart loop")


def test_nondegenerate_reduction_rejects_false_operation_certificate(monkeypatch) -> None:
    """Every candidate reduction is checked before its norm can affect the loop."""
    from sage.all import ZZ, identity_matrix, matrix
    import sage_indefinite_port.indefinite.indefinite_reduction as reduction

    gram = matrix(ZZ, [[0, 1], [1, 0]])
    identity = identity_matrix(ZZ, 2)
    monkeypatch.setattr(
        reduction,
        "_compute_indefinite_reduction",
        lambda current: (identity, current + identity, None),
    )
    try:
        reduction._nondegenerate_reduction(gram)
    except ArithmeticError as error:
        assert "false Gram certificate" in str(error)
    else:
        raise AssertionError("a false reduction candidate entered the descent loop")


def test_indefinite_reduction_integrates_source_radical_and_completion() -> None:
    """The reduction orchestration preserves congruence with degenerate forms."""
    from sage.all import ZZ, matrix
    from sage_indefinite_port.indefinite.indefinite_reduction import _indefinite_reduction

    forms = (
        matrix(ZZ, [[1, 0, 1], [0, 1, 1], [1, 1, 2]]),
        matrix(ZZ, [[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
        matrix(ZZ, [[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        matrix(ZZ, 3, 3, 0),
    )
    for gram in forms:
        basis, reduced = _indefinite_reduction(gram)
        assert abs(basis.det()) == 1
        assert reduced == basis * gram * basis.transpose()
        assert reduced.rank() == gram.rank()
        radical_rank = gram.nrows() - gram.rank()
        assert not reduced.matrix_from_rows(range(gram.nrows() - radical_rank, gram.nrows()))


def test_indefinite_gram_schmidt_rejects_invalid_form_or_basis() -> None:
    """The indefinite Gram-Schmidt source operation requires a lattice Gram basis."""
    from sage.all import ZZ, identity_matrix, matrix
    from sage_indefinite_port.indefinite.indefinite_reduction import _indefinite_gram_schmidt

    gram = matrix(ZZ, [[1, 1], [1, -1]])
    assert _indefinite_gram_schmidt(gram, identity_matrix(ZZ, 2))[0] is not None
    invalid = (
        (matrix(ZZ, [[1, 2], [0, 1]]), identity_matrix(ZZ, 2)),
        (gram, matrix(ZZ, [[1, 0], [0, 0]])),
        (gram, matrix(ZZ, 2, 3, 0)),
    )
    for form, basis in invalid:
        try:
            _indefinite_gram_schmidt(form, basis)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid Gram-Schmidt input was accepted")


def test_indefinite_lll_nearest_integer_uses_source_half_tie_rule() -> None:
    """The signed LLL size-reduction operation uses GAP's exact nearest integer."""
    from sage.all import QQ
    from sage_indefinite_port.indefinite.indefinite_reduction import _source_nearest_integer

    cases = {
        QQ(1) / 2: 0,
        -QQ(1) / 2: -1,
        QQ(3) / 2: 1,
        -QQ(3) / 2: -2,
        QQ(7) / 5: 1,
        -QQ(7) / 5: -1,
    }
    assert {value: _source_nearest_integer(value) for value in cases} == cases


def test_source_integer_nullspace_reduces_and_completes_radical() -> None:
    """The production nullspace seam returns the pinned reduced primitive row basis."""
    from sage.all import ZZ, matrix
    from sage_indefinite_port.indefinite.indefinite_reduction import _source_integer_nullspace

    forms = (
        matrix(ZZ, [[1, 0, 1], [0, 1, 1], [1, 1, 2]]),
        matrix(ZZ, [[1, 0, 0], [0, 0, 0], [0, 0, 0]]),
        matrix(ZZ, 3, 3, 0),
    )
    for gram in forms:
        kernel = _source_integer_nullspace(gram)
        assert kernel.nrows() == gram.nrows() - gram.rank()
        assert kernel.ncols() == gram.ncols()
        assert not kernel * gram
        if kernel.nrows():
            assert kernel.rank() == kernel.nrows()


def test_source_subspace_completion_int_returns_unimodular_complements() -> None:
    """The pinned Smith-transform completion extends primitive row lattices."""
    from sage.all import ZZ, identity_matrix, matrix
    from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int

    cases = (
        (matrix(ZZ, [[-1, -1, 1]]), 3),
        (matrix(ZZ, [[1, 0, 1], [0, 1, 1]]), 3),
        (matrix(ZZ, [], ncols=3), 3),
        (matrix(ZZ, [[1, 1, 0, 0], [0, 0, 1, 1]]), 4),
    )
    for kernel, dimension in cases:
        completion = _source_subspace_completion_int(kernel, dimension)
        full_basis = completion.stack(kernel)
        assert full_basis.nrows() == dimension
        assert abs(full_basis.det()) == 1
        if kernel.nrows():
            diagonal, row_transform, column_transform = kernel.smith_form(transformation=True)
            assert row_transform * kernel * column_transform == diagonal
            block = identity_matrix(ZZ, dimension)
            for i in range(kernel.nrows()):
                for j in range(kernel.nrows()):
                    block[i, j] = row_transform[i, j]
            source_basis = block.inverse() * column_transform.inverse()
            assert full_basis == source_basis

    try:
        _source_subspace_completion_int(matrix(ZZ, [[2, 0, 0]]), 3)
    except ValueError:
        pass
    else:
        raise AssertionError("nonprimitive sublattice unexpectedly completed")


def test_source_v2_honors_critical_norm_and_strictness() -> None:
    """V2 returns the first integral approximant satisfying the requested bound."""
    from sage.all import QQ, ZZ, identity_matrix, matrix, vector

    gram = matrix(ZZ, [[1, 0], [0, -1]])
    basis = identity_matrix(ZZ, 2)
    def diagonalize(_gram):
        return identity_matrix(QQ, 2), matrix(QQ, [[1, 0], [0, -1]])

    def no_eigen_directions(_gram):
        return ()

    inclusive = vector(ZZ, _source_integral_vector_allmeth_v2(
        gram, basis, gram, QQ.one(), False,
        diagonalize=diagonalize, eigen_directions=no_eigen_directions,
    ))
    strict = vector(ZZ, _source_integral_vector_allmeth_v2(
        gram, basis, gram, QQ.one(), True,
        diagonalize=diagonalize, eigen_directions=no_eigen_directions,
    ))
    assert inclusive == vector(ZZ, [1, 0])
    assert strict == vector(ZZ, [2, 0])
    assert (inclusive * gram * inclusive.column())[0] == 1
    assert (strict * gram * strict.column())[0] == 4

    def false_diagonalization(_gram):
        return identity_matrix(QQ, 2), matrix(QQ, [[2, 0], [0, -1]])
    try:
        _source_integral_vector_allmeth_v2(
            gram, basis, gram, QQ.one(), False,
            diagonalize=false_diagonalization, eigen_directions=no_eigen_directions,
        )
    except ArithmeticError:
        pass
    else:
        raise AssertionError("V2 accepted an uncertified exact diagonalization")

    negative = matrix(ZZ, [[-2, 0], [0, -1]])
    try:
        _source_integral_vector_allmeth_v2(
            negative, identity_matrix(ZZ, 2), negative,
            diagonalize=lambda _gram: (identity_matrix(QQ, 2), negative.change_ring(QQ)),
            eigen_directions=no_eigen_directions,
        )
    except ValueError as error:
        assert "no positive direction" in str(error)
    else:
        raise AssertionError("V2 entered an unbounded iterator search without a positive direction")


def test_source_integral_vector_allmeth_v2_returns_exact_positive_witness() -> None:
    """Exercise the complete V2 iterator consumer with an identity reduction."""
    from sage.all import QQ, RDF, ZZ, diagonal_matrix, identity_matrix, matrix, vector

    forms = (
        matrix(ZZ, [[1, 0], [0, -1]]),
        matrix(ZZ, [[0, 1], [1, 0]]),
        matrix(ZZ, [[0, 0, 0], [0, 1, 0], [0, 0, -2]]),
        diagonal_matrix(ZZ, [3, -2, 0, 5]),
    )
    for gram in forms:
        basis = identity_matrix(ZZ, gram.nrows())
        candidate = vector(ZZ, _source_integral_vector_allmeth_v2(gram, basis, gram))
        assert (candidate * gram * candidate.column())[0] > 0

    def no_eigen_directions(_gram):
        return ()

    rational_basis = matrix(QQ, [[1, QQ(1) / 2], [0, 1]])
    rational_gram = matrix(QQ, [[1, 0], [0, -1]])
    rational_reduced = rational_basis * rational_gram * rational_basis.transpose()
    try:
        _source_integral_vector_allmeth_v2(
            rational_gram, rational_basis, rational_reduced,
            eigen_directions=no_eigen_directions,
        )
    except ValueError as error:
        assert "integral entries" in str(error)
    else:
        raise AssertionError("a rational nonintegral reduction basis was accepted")

    inexact_gram = rational_gram.change_ring(RDF)
    try:
        _source_integral_vector_allmeth_v2(
            inexact_gram, identity_matrix(ZZ, 2), inexact_gram,
            eigen_directions=no_eigen_directions,
        )
    except ValueError as error:
        assert "exact ring" in str(error)
    else:
        raise AssertionError("an inexact floating-point Gram matrix was accepted")

    inexact_basis = matrix(RDF, [[1.0, 0.0], [0.0, 1.0]])
    try:
        _source_integral_vector_allmeth_v2(
            rational_gram, inexact_basis, rational_gram,
            eigen_directions=no_eigen_directions,
        )
    except ValueError as error:
        assert "exact ring" in str(error)
    else:
        raise AssertionError("an inexact floating-point reduction basis was accepted")

    original = matrix(ZZ, [[2, 1], [1, -1]])
    basis = matrix(ZZ, [[1, 2], [0, 1]])
    reduced = basis * original * basis.transpose()
    candidate = vector(ZZ, _source_integral_vector_allmeth_v2(original, basis, reduced))
    assert (candidate * original * candidate.column())[0] > 0

    negative = diagonal_matrix(ZZ, [-1, -2])
    try:
        _source_integral_vector_allmeth_v2(negative, identity_matrix(ZZ, 2), negative)
    except ValueError as error:
        assert "no positive direction" in str(error)
    else:
        raise AssertionError("V2 iterator accepted a form with no positive direction")


def test_source_remove_fraction_vector_returns_primitive_integral_coordinates() -> None:
    """Clearing rational denominators also removes the integral content."""
    from sage.all import QQ

    assert _source_remove_fraction_vector((QQ(1) / 2, -QQ(1) / 2)) == (1, -1)
    assert _source_remove_fraction_vector((QQ.zero(), QQ.zero())) == (0, 0)


def test_source_integral_reduction_rejects_false_basis_certificates() -> None:
    """The V2 consumer receives only exact unimodular Gram congruences."""
    from sage.all import RDF, ZZ, identity_matrix, matrix

    gram = matrix(ZZ, [[2, 1], [1, -1]])
    identity = identity_matrix(ZZ, 2)
    _source_validate_integral_reduction(gram, identity, gram)
    for bad_gram, bad_basis, bad_reduced in (
        (gram.change_ring(RDF), identity, gram),
        (gram, identity, gram.change_ring(RDF)),
    ):
        try:
            _source_validate_integral_reduction(bad_gram, bad_basis, bad_reduced)
        except ValueError as error:
            assert "exact ring" in str(error)
        else:
            raise AssertionError("an inexact Gram certificate was accepted")
    for basis, reduced in (
        (2 * identity, 4 * gram),
        (identity, gram + identity),
        (matrix(ZZ, [[1, 0], [0, 0]]), gram),
    ):
        try:
            _source_validate_integral_reduction(gram, basis, reduced)
        except (ValueError, ArithmeticError):
            pass
        else:
            raise AssertionError("invalid integral reduction certificate was accepted")


def test_source_positive_integral_vector_contract_certifies_exact_norm() -> None:
    """The public consumer rejects zero, nonpositive, and malformed witnesses."""
    from sage.all import ZZ, matrix

    gram = matrix(ZZ, [[1, 0], [0, -1]])
    assert _source_certify_positive_integral_vector(gram, (ZZ(2), ZZ(1))) == (ZZ(2), ZZ(1))
    for candidate in ((ZZ(0), ZZ(0)), (ZZ(1), ZZ(1)), (ZZ(1),), (0.5, 1)):
        try:
            _source_certify_positive_integral_vector(gram, candidate)
        except (ValueError, ArithmeticError):
            pass
        else:
            raise AssertionError(f"invalid positive-vector witness accepted: {candidate}")


def test_source_integral_vector_allmeth_rejects_invalid_gram_before_reduction() -> None:
    """The public lattice-vector search accepts only nonempty symmetric Gram matrices."""
    from sage.all import QQ, ZZ, matrix

    for gram in (matrix(ZZ, 0, 0), matrix(ZZ, [[1, 2], [0, 1]])):
        try:
            _source_integral_vector_allmeth(
                gram, QQ.zero(), True,
                reduction=lambda _gram: (_ for _ in ()).throw(AssertionError("reduction called")),
                isotropic_reduction=lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("isotropic reduction called")),
            )
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid Gram matrix accepted: {gram}")


def test_source_integral_vector_allmeth_selects_isotropic_shortcut_exactly() -> None:
    """The generic source wrapper chooses the isotropic branch from inequality semantics."""
    from sage.all import QQ, ZZ, identity_matrix, matrix

    gram = matrix(ZZ, [[0, 0], [0, 1]])
    calls = []
    def isotropic_reduction(_gram, *, look_for_isotropic):
        assert look_for_isotropic
        calls.append("isotropic")
        return identity_matrix(ZZ, 2), gram, (QQ(1) / 2, QQ.zero())
    def reduction(_gram):
        calls.append("reduction")
        return identity_matrix(ZZ, 2), gram
    def find_vector(_gram, _basis, _reduced, critical_norm, strict_inequality):
        calls.append(("v2", critical_norm, strict_inequality))
        return (ZZ(0), ZZ(1))

    isotropic = _source_integral_vector_allmeth(
        gram, QQ.zero(), False, reduction=reduction,
        isotropic_reduction=isotropic_reduction, find_vector=find_vector,
    )
    assert isotropic == (ZZ(1), ZZ(0))
    assert calls == ["isotropic"]

    calls.clear()
    positive = _source_integral_vector_allmeth(
        gram, QQ.zero(), True, reduction=reduction,
        isotropic_reduction=isotropic_reduction, find_vector=find_vector,
    )
    assert positive == (ZZ(0), ZZ(1))
    assert calls == ["reduction", ("v2", QQ.zero(), True)]

    calls.clear()
    negative_critical = _source_integral_vector_allmeth(
        gram, -QQ.one(), True, reduction=reduction,
        isotropic_reduction=isotropic_reduction, find_vector=find_vector,
    )
    assert negative_critical == (ZZ(1), ZZ(0))
    assert calls == ["isotropic"]

    calls.clear()
    def no_isotropic(_gram, *, look_for_isotropic):
        assert look_for_isotropic
        calls.append("isotropic-no-witness")
        return identity_matrix(ZZ, 2), gram, None
    fallback = _source_integral_vector_allmeth(
        gram, -QQ.one(), True, reduction=reduction,
        isotropic_reduction=no_isotropic, find_vector=find_vector,
    )
    assert fallback == (ZZ(0), ZZ(1))
    assert calls == ["isotropic-no-witness", ("v2", -QQ.one(), True)]


def test_integral_positive_vector_consumer_uses_pinned_diagonalization() -> None:
    """The real positive-vector consumer returns exact positive-norm witnesses."""
    from sage.all import ZZ, diagonal_matrix, matrix, vector

    forms = (
        matrix(ZZ, [[1, 0], [0, -1]]),
        matrix(ZZ, [[0, 1], [1, 0]]),
        matrix(ZZ, [[0, 0, 0], [0, 1, 0], [0, 0, -2]]),
        diagonal_matrix(ZZ, [3, -2, 0, 5]),
    )
    for gram in forms:
        candidate = vector(ZZ, _integral_positive_vector_allmeth(gram))
        assert candidate in ZZ**gram.nrows()
        assert (candidate * gram * candidate.column())[0] > 0

    for gram in (matrix(ZZ, [[-1, 0], [0, -2]]), matrix(ZZ, 2, 2, 0)):
        try:
            _integral_positive_vector_allmeth(gram)
        except ValueError as error:
            assert "no positive direction" in str(error)
        else:
            raise AssertionError(f"form without positive direction accepted: {gram}")


def test_pinned_symmetric_diagonalization_returns_exact_congruence_witness() -> None:
    """Exercise the production source-ordered diagonalizer on degenerate forms."""
    from sage.all import QQ, diagonal_matrix, matrix

    forms = (
        matrix(QQ, [[1, 2, 3], [2, 4, 6], [3, 6, 9]]),
        matrix(QQ, [[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
        matrix(QQ, [[0, 0, 0], [0, -2, 0], [0, 0, 3]]),
        matrix(QQ, [[0, 1], [1, 0]]),
        diagonal_matrix(QQ, [3, -2, 0, 5]),
        matrix(QQ, 3, 3, 0),
        matrix(QQ, 0, 0),
    )
    expected_transforms = (
        matrix(QQ, [[-2, 1, 0], [-3, 0, 1], [0, 0, 1]]),
        matrix(QQ, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        matrix(QQ, [[1, 1], [-1, 1]]),
        matrix(QQ, [[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
    )
    for gram, expected_transform in zip((forms[0], forms[2], forms[3], forms[4]), expected_transforms, strict=True):
        transform, _diagonal = _pinned_diagonalize_symmetric_matrix(gram)
        assert transform == expected_transform

    for gram in forms:
        transform, diagonal = _pinned_diagonalize_symmetric_matrix(gram)
        assert transform.nrows() == gram.nrows()
        assert transform.ncols() == gram.ncols()
        assert transform.det() != 0
        assert transform * gram * transform.transpose() == diagonal
        assert all(diagonal[i, j] == 0 for i in range(diagonal.nrows()) for j in range(diagonal.ncols()) if i != j)


def test_pinned_diagonalizer_handles_deterministic_rational_symmetric_samples() -> None:
    """Check exact congruence and nonsingularity on degenerate rational samples."""
    from random import Random

    from sage.all import QQ, matrix

    random = Random(314159265)
    for size in range(1, 7):
        for _ in range(10):
            gram = matrix(
                QQ,
                size,
                size,
                lambda i, j: QQ(random.randint(-4, 4)) / random.choice((1, 2, 3)) if i <= j else 0,
            )
            for i in range(size):
                for j in range(i):
                    gram[i, j] = gram[j, i]
            transform, diagonal = _pinned_diagonalize_symmetric_matrix(gram)
            assert transform.det() != 0
            assert transform * gram * transform.transpose() == diagonal
            assert all(diagonal[i, j] == 0 for i in range(size) for j in range(size) if i != j)


def test_pinned_nullspace_preserves_gap_row_basis_order() -> None:
    """The source nullspace rows differ from Sage's canonical kernel basis."""
    from random import Random

    from sage.all import QQ, matrix

    gram = matrix(QQ, [[1, 2, 3], [2, 4, 6], [3, 6, 9]])
    pinned_basis = _pinned_nullspace_matrix(gram)
    assert pinned_basis == matrix(QQ, [[-2, 1, 0], [-3, 0, 1]])
    assert pinned_basis != gram.left_kernel().basis_matrix()

    nondegenerate = matrix(QQ, [[1, 2, 0], [2, 5, 1], [0, 1, 3]])
    previous = matrix(QQ, [[1, 0, 0]])
    orthogonal_constraints = nondegenerate * previous.transpose()
    pinned_complement = _pinned_nullspace_matrix(orthogonal_constraints)
    assert pinned_complement * orthogonal_constraints == 0
    assert pinned_complement != orthogonal_constraints.left_kernel().basis_matrix()
    assert _pinned_nullspace_matrix(matrix(QQ, [[1, 0], [0, 1]])).nrows() == 0
    assert _pinned_nullspace_matrix(matrix(QQ, 0, 3)).dimensions() == (0, 0)
    rectangular_full_row_rank = matrix(QQ, [[1, 0, 0], [0, 1, 0]])
    empty_left_kernel = _pinned_nullspace_matrix(rectangular_full_row_rank)
    assert empty_left_kernel.dimensions() == (0, 2)

    rational_gram = matrix(QQ, [[QQ(1) / 2, 1, QQ(3) / 2], [1, 2, 3]])
    rational_kernel = _pinned_nullspace_matrix(rational_gram)
    assert rational_kernel * rational_gram == 0
    assert rational_kernel.nrows() == rational_gram.nrows() - rational_gram.rank()

    rank_one = matrix(QQ, [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    rank_one_kernel = _pinned_nullspace_matrix(rank_one)
    assert rank_one_kernel * rank_one == 0
    assert rank_one_kernel.nrows() == 2

    random = Random(481516)
    for size in range(1, 6):
        for _ in range(12):
            symmetric = matrix(
                QQ,
                size,
                size,
                lambda i, j: random.randint(-4, 4) if i <= j else 0,
            )
            for i in range(size):
                for j in range(i):
                    symmetric[i, j] = symmetric[j, i]
            kernel = _pinned_nullspace_matrix(symmetric)
            assert kernel * symmetric == 0
            assert kernel.nrows() == size - symmetric.rank()


def test_source_random_perturbations_match_pinned_matrices() -> None:
    """Check both GetRandomMatrixPerturbation branches, including diagonal shear."""
    from sage.all import QQ, matrix

    signed_permutation = _source_permutation_sign_matrix([2, 0, 1], (-1, 1, -1))
    assert signed_permutation == matrix(QQ, [[0, 0, -1], [1, 0, 0], [0, -1, 0]])
    assert abs(signed_permutation.det()) == 1

    for row, column, sign in ((0, 1, 1), (2, 0, -1), (1, 1, 1), (1, 1, -1)):
        elementary = _source_elementary_perturbation(3, row, column, sign)
        expected = matrix(QQ, 3, 3, lambda i, j: int(i == j))
        expected[row, column] = sign
        assert elementary == expected
        assert abs(elementary.det()) == 1

    for invalid in ((3, 0, 1), (0, 3, 1), (0, 0, 0)):
        try:
            _source_elementary_perturbation(3, *invalid)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid perturbation accepted: {invalid}")


def test_source_indefinite_short_positive_vector_integrates_pinned_search_loop() -> None:
    """Run the complete Graver/reframing/termination loop with deterministic seams."""
    from sage.all import QQ, RDF, ZZ, identity_matrix, matrix

    gram = matrix(ZZ, [[10, 0], [0, 10]])
    calls = []
    def find_positive(perturbed):
        calls.append(perturbed)
        assert perturbed == gram
        return (ZZ(1), ZZ(0))
    result = _source_indefinite_short_positive_vector(
        gram, find_positive, lambda dimension: identity_matrix(ZZ, dimension)
    )
    assert result == (ZZ(1), ZZ(0))
    assert len(calls) == 501

    calls.clear()
    try:
        _source_indefinite_short_positive_vector(
            gram, find_positive,
            lambda dimension: matrix(QQ, [[1, QQ(1) / 2], [0, 1]]),
        )
    except ArithmeticError as error:
        assert "integral lattice automorphism" in str(error)
    else:
        raise AssertionError("a rational nonintegral perturbation was accepted")
    assert len(calls) == 1

    calls.clear()
    try:
        _source_indefinite_short_positive_vector(
            gram, find_positive,
            lambda dimension: matrix(RDF, [[1.0, 0.0], [0.0, 1.0]]),
        )
    except ArithmeticError as error:
        assert "exact integral matrix" in str(error)
    else:
        raise AssertionError("an inexact floating-point perturbation was accepted")
    assert len(calls) == 1

    short_gram = matrix(ZZ, [[1, 0], [0, -1]])
    short_calls = []
    def find_short(_perturbed):
        short_calls.append(None)
        return (ZZ(1), ZZ(0))
    short_result = _source_indefinite_short_positive_vector(
        short_gram, find_short, lambda dimension: identity_matrix(ZZ, dimension)
    )
    assert short_result == (ZZ(1), ZZ(0))
    assert len(short_calls) == 1

    # The source always requests degree two, even for a one-dimensional form;
    # BuildSet(1, 2) simply has no support of size two.
    rank_one = matrix(ZZ, [[1]])
    rank_one_result = _source_indefinite_short_positive_vector(
        rank_one,
        lambda _perturbed: (ZZ(1),),
        lambda dimension: identity_matrix(ZZ, dimension),
    )
    assert rank_one_result == (ZZ(1),)

    # Repeated source elementary perturbations grow a shear until the strict
    # 1000x L1 threshold resets the accumulated frame.
    reset_gram = matrix(ZZ, [[10, 0], [0, 10]])
    seen_growth = False
    reset_observed = []
    def find_reset(perturbed):
        nonlocal seen_growth
        if perturbed[0, 0] > 10000:
            seen_growth = True
        elif seen_growth and perturbed[0, 0] == 20:
            reset_observed.append(True)
        return (ZZ(1), ZZ(0))
    def shear(_dimension):
        from sage.all import matrix
        return matrix(ZZ, [[1, 1], [0, 1]])
    reset_result = _source_indefinite_short_positive_vector(reset_gram, find_reset, shear)
    assert reset_result == (ZZ(1), ZZ(0))
    assert reset_observed


def test_source_indefinite_short_positive_vector_certifies_search_witnesses() -> None:
    """The public heuristic rejects malformed witnesses before perturbing again."""
    from sage.all import ZZ, identity_matrix, matrix

    gram = matrix(ZZ, [[1]])
    result = _source_indefinite_short_positive_vector(gram, lambda _gram: (ZZ(1),))
    assert result == (ZZ(1),)
    assert all(type(entry) is type(ZZ(0)) for entry in result)
    for bad_gram, witness in (
        (matrix(ZZ, 0, 0), ()),
        (matrix(ZZ, [[1, 2], [0, 1]]), (ZZ(1), ZZ(1))),
        (gram, (0.5,)),
        (gram, (ZZ(0),)),
        (gram, (ZZ(-1),)),
    ):
        try:
            _source_indefinite_short_positive_vector(bad_gram, lambda _gram, w=witness: w)
        except (ValueError, ArithmeticError):
            pass
        else:
            raise AssertionError(f"invalid positive-search input or witness accepted: {bad_gram}, {witness}")


def test_source_positive_search_best_matches_pinned_termination_state() -> None:
    """Best norm, no-improvement count, and both source stop conditions agree."""
    from sage.all import ZZ, matrix, vector

    gram = matrix(ZZ, [[1, 0], [0, -1]])
    first = vector(ZZ, [4, 1])
    best, norm, failures, stop = _source_positive_search_best(gram, first, None, None, 0)
    assert best == first and norm == 15 and failures == 0 and not stop

    worse = vector(ZZ, [5, 1])
    best, norm, failures, stop = _source_positive_search_best(gram, worse, best, norm, failures)
    assert best == first and norm == 15 and failures == 1 and not stop

    short = vector(ZZ, [2, 1])
    best, norm, failures, stop = _source_positive_search_best(gram, short, best, norm, failures)
    assert best == short and norm == 3 and failures == 0 and stop

    best, norm, failures, stop = _source_positive_search_best(gram, worse, best, norm, 499)
    assert best == short and norm == 3 and failures == 500 and stop

    for state in ((first, None, 0), (None, 15, 0), (first, 14, 0), (first, 15, 500), (first, 15, 501)):
        try:
            _source_positive_search_best(gram, worse, state[0], state[1], state[2])
        except (ValueError, ArithmeticError):
            pass
        else:
            raise AssertionError(f"inconsistent positive-search state accepted: {state}")

    for invalid_gram in (matrix(ZZ, 0, 0), matrix(ZZ, [[1, 2], [0, 1]])):
        try:
            _source_positive_search_best(invalid_gram, worse, None, None, 0)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid positive-search Gram accepted: {invalid_gram}")


def test_source_positive_perturbation_reset_uses_strict_factor_1000() -> None:
    """Match Positivity.h:660, including equality at the threshold."""
    assert not _source_should_reset_positive_perturbation(0, 0)
    assert not _source_should_reset_positive_perturbation(1000, 1)
    assert _source_should_reset_positive_perturbation(1001, 1)
    assert not _source_should_reset_positive_perturbation(10000, 10)
    assert _source_should_reset_positive_perturbation(10001, 10)
    for norms in ((-1, 1), (1, -1)):
        try:
            _source_should_reset_positive_perturbation(*norms)
        except ValueError:
            pass
        else:
            raise AssertionError(f"negative L1 norm accepted: {norms}")


def test_source_direct_graver_improvement_matches_literal_ordered_passes() -> None:
    """The production descent performs repeated full Graver-basis passes."""
    from sage.all import ZZ, matrix, vector

    gram = matrix(ZZ, [[2, 1, 0], [1, -3, 1], [0, 1, 2]])
    directions = tuple(vector(ZZ, row) for row in _source_graver_k_basis(3, 2))
    for initial in (vector(ZZ, [4, 3, 1]), vector(ZZ, [-3, 2, 4]), vector(ZZ, [1, 1, 1])):
        expected = vector(ZZ, initial)
        while True:
            changes = 0
            for direction in directions:
                before = int((expected * gram * expected.column())[0])
                mixed = int((expected * gram * direction.column())[0])
                direction_norm = int((direction * gram * direction.column())[0])
                alpha = 0
                source_norm = before
                while True:
                    candidate = before + 2 * mixed * (alpha + 1) + direction_norm * (alpha + 1) ** 2
                    if candidate >= source_norm or candidate <= 0:
                        break
                    source_norm = candidate
                    alpha += 1
                expected += alpha * direction
                changes += alpha
            if changes == 0:
                break
        actual = _source_direct_graver_improvement(gram, initial, directions)
        assert actual == expected
        assert (actual * gram * actual.column())[0] > 0


def test_source_direct_graver_improvement_rejects_false_descent_prefix(monkeypatch) -> None:
    """The update boundary verifies every claimed GetAlpha step."""
    from sage.all import ZZ, vector
    import sage_indefinite_port.indefinite.lorentzian_cells as cells

    gram = matrix(ZZ, [[1]])
    monkeypatch.setattr(cells, "_quadratic_descent_prefix", lambda *_args: 1)
    try:
        _source_direct_graver_improvement(gram, vector(ZZ, [1]), (vector(ZZ, [1]),))
    except ArithmeticError as error:
        assert "nonpositive or nondecreasing" in str(error)
    else:
        raise AssertionError("a false Graver descent prefix was accepted")


def test_source_graver_k_basis_matches_support_and_sign_enumeration() -> None:
    """GetGraverKbasis enumerates every signed support through its requested degree."""
    from itertools import combinations, product

    for dimension in range(1, 6):
        for degree in range(1, dimension + 3):
            basis = _source_graver_k_basis(dimension, degree)
            expected = []
            for size in range(1, min(degree, dimension) + 1):
                for support in combinations(range(dimension), size):
                    for signs in product((-1, 1), repeat=size):
                        row = [0] * dimension
                        for coordinate, sign in zip(support, signs, strict=True):
                            row[coordinate] = sign
                        expected.append(tuple(row))
            assert basis == tuple(expected)
            assert len(basis) == sum(2**size * __import__('math').comb(dimension, size) for size in range(1, min(degree, dimension) + 1))
            assert len(set(basis)) == len(basis)

    for dimension, degree in ((0, 1), (2, 0)):
        try:
            _source_graver_k_basis(dimension, degree)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid Graver K-basis parameters accepted: {(dimension, degree)}")


def test_source_random_matrix_perturbation_dispatches_both_pinned_branches() -> None:
    """The public random perturbation operation composes the pinned constructors."""
    from sage.all import ZZ, matrix

    choices = iter((0, -1, 1, -1))
    def choose(options):
        value = next(choices)
        assert value in options
        return value
    def reverse(values):
        values.reverse()

    signed = _source_random_matrix_perturbation(3, choose=choose, shuffle_values=reverse)
    assert signed == matrix(ZZ, [[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

    choices = iter((1, -1))
    def choose_elementary(options):
        value = next(choices)
        assert value in options
        return value
    indices = iter((1, 2))
    elementary = _source_random_matrix_perturbation(
        3, choose=choose_elementary, choose_index=lambda _dimension: next(indices)
    )
    assert elementary == matrix(ZZ, [[1, 0, 0], [0, 1, -1], [0, 0, 1]])
    assert abs(signed.det()) == abs(elementary.det()) == 1


def test_source_perturbation_matrices_are_unimodular_in_small_dimensions() -> None:
    """Exhaust signed permutations and all single-entry perturbations through n=4."""
    from itertools import permutations, product

    from sage.all import QQ

    for dimension in range(1, 5):
        for permutation in permutations(range(dimension)):
            for signs in product((-1, 1), repeat=dimension):
                signed = _source_permutation_sign_matrix(list(permutation), signs)
                assert abs(signed.det()) == 1
                assert all(sum(entry != 0 for entry in row) == 1 for row in signed.rows())
                assert all(sum(entry != 0 for entry in signed.column(j)) == 1 for j in range(dimension))
        for row in range(dimension):
            for column in range(dimension):
                for sign in (-1, 1):
                    elementary = _source_elementary_perturbation(dimension, row, column, sign)
                    assert abs(elementary.det()) == 1
                    assert elementary.base_ring() is QQ


def test_pinned_nearest_integer_rounds_half_ties_toward_positive_infinity() -> None:
    """Match Fractions.g:138-148, including negative half-integers."""
    from sage.all import QQ

    cases = (
        (QQ(-5) / 2, -2),
        (QQ(-3) / 2, -1),
        (QQ(-1) / 2, 0),
        (QQ(1) / 2, 1),
        (QQ(3) / 2, 2),
        (QQ(5) / 2, 3),
        (QQ(7) / 3, 2),
        (QQ(-7) / 3, -2),
    )
    for value, expected in cases:
        assert _pinned_nearest_integer(value) == expected


def test_source_negative_vector_respects_radical_and_diagonal_order() -> None:
    """The pinned search selects the first negative-or-zero diagonal row."""
    from sage.all import QQ, matrix, vector

    first_gram = matrix(QQ, [[-1, 0], [0, 1]])
    first = _source_some_negative_vector(first_gram)
    assert first == (1, 0)
    assert (vector(QQ, first) * first_gram * vector(QQ, first).column())[0] < 0

    # The radical contributes a zero candidate first; the returned vector
    # must still have negative norm and remain in the non-radical coordinates.
    radical_gram = matrix(QQ, [[0, 0, 0], [0, -1, 0], [0, 0, 1]])
    result = _source_some_negative_vector(radical_gram)
    assert result == (0, 1, 0)
    result_vector = vector(QQ, result)
    assert (result_vector * radical_gram * result_vector.column())[0] < 0

    hyperbolic_gram = matrix(QQ, [[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    assert _source_some_negative_vector(hyperbolic_gram) == (-1, 1, 0)


def test_source_negative_vector_certifies_congruence_random_indefinite_forms() -> None:
    """Every returned source candidate has negative norm after integral congruence."""
    from random import Random

    from sage.all import QQ, ZZ, diagonal_matrix, identity_matrix, matrix, vector

    random = Random(20261026)
    for size in range(2, 6):
        for _ in range(8):
            change = identity_matrix(ZZ, size)
            for _ in range(3 * size):
                i, j = random.sample(range(size), 2)
                change[i, :] = change[i, :] + random.choice((-2, -1, 1, 2)) * change[j, :]
            diagonal = diagonal_matrix(QQ, [-1] + [1] * (size - 1))
            gram = matrix(QQ, change) * diagonal * matrix(QQ, change).transpose()
            candidate = vector(QQ, _source_some_negative_vector(gram))
            assert (candidate * gram * candidate.column())[0] < 0


def test_source_negative_vector_certifies_rational_congruence_forms() -> None:
    """Exercise exact rational congruences, including a degenerate form."""
    from random import Random

    from sage.all import QQ, diagonal_matrix, identity_matrix, matrix, vector

    random = Random(8675309)
    for size in range(2, 6):
        for radical_rank in (0, 1):
            diagonal = [-1] + [1] * (size - radical_rank - 1) + [0] * radical_rank
            for _ in range(4):
                change = identity_matrix(QQ, size)
                for _ in range(2 * size):
                    i, j = random.sample(range(size), 2)
                    change[i, :] = change[i, :] + QQ(random.choice((-3, -1, 1, 3)))/2 * change[j, :]
                gram = change * diagonal_matrix(QQ, diagonal) * change.transpose()
                candidate = vector(QQ, _source_some_negative_vector(gram))
                assert (candidate * gram * candidate.column())[0] < 0


def test_source_mid_value_matches_pinned_continued_fraction_sequence() -> None:
    """Compare the selected middle value with every pinned CF convergent."""
    from fractions import Fraction

    from sage.all import QQ

    def source_continued_fraction(value: Fraction, level: int) -> Fraction:
        if value.denominator == 1:
            return value
        residue = value.numerator % value.denominator
        integer = value - Fraction(residue, value.denominator)
        if level == 0:
            return integer
        return integer + 1 / source_continued_fraction(
            Fraction(value.denominator, residue), level - 1
        )

    for lower_num in range(-5, 6):
        for lower_den in range(1, 6):
            lower = Fraction(lower_num, lower_den)
            for width_num in range(1, 8):
                for width_den in range(1, 6):
                    upper = lower + Fraction(width_num, width_den)
                    midpoint = (lower + upper) / 2
                    target_low = (2 * lower + upper) / 3
                    target_high = (lower + 2 * upper) / 3
                    expected = None
                    for level in range(32):
                        candidate = source_continued_fraction(midpoint, level)
                        if target_low <= candidate <= target_high:
                            expected = candidate
                            break
                        if candidate == midpoint:
                            break
                    assert expected is not None
                    actual = _source_mid_value(QQ(lower.numerator) / lower.denominator,
                                               QQ(upper.numerator) / upper.denominator)
                    assert actual == QQ(expected.numerator) / expected.denominator, (lower, upper, expected, actual)


def test_quadratic_descent_prefix_matches_pinned_get_alpha() -> None:
    """Compare the accelerated prefix with the literal source loop."""
    for base_norm in range(1, 13):
        for mixed in range(-10, 11):
            for direction_norm in range(-5, 6):
                source_norm = base_norm
                expected = 0
                while True:
                    step = expected + 1
                    candidate_norm = (
                        base_norm + 2 * mixed * step + direction_norm * step * step
                    )
                    if candidate_norm >= source_norm or candidate_norm <= 0:
                        break
                    source_norm = candidate_norm
                    expected += 1
                assert _quadratic_descent_prefix(
                    base_norm, mixed, direction_norm
                ) == expected, (base_norm, mixed, direction_norm, expected)


def test_flip_upper_bound_uses_rational_isotropic_exit_and_rejects_coincident_exit() -> None:
    """Match GetUpperBound's isotropic root and constant-bound collision."""
    from sage.all import QQ, matrix, vector

    gram = matrix(QQ, [[1, 0], [0, -1]])
    base = vector(QQ, [0, 1, 0])
    direction = vector(QQ, [0, 0, 1])
    # q((1,0) + t(0,1)) = 1 - t^2, so the first positive isotropic
    # parameter is exactly one; the doubling search's bound is two.
    assert _upper_bound(gram, base, direction) == QQ.one()

    # The affine constant also vanishes at t=1. The pinned source returns
    # fail in this coincident-bound case rather than accepting that boundary.
    base_with_constant = vector(QQ, [-1, 1, 0])
    direction_with_constant = vector(QQ, [1, 0, 1])
    assert _upper_bound(gram, base_with_constant, direction_with_constant) is None


def test_pinned_sublattice_basis_reduction_preserves_row_lattice() -> None:
    """Alternating direct/dual LLL decreases L1 without changing the lattice."""
    from sage.all import ZZ, matrix

    # Zero- and one-row bases are fixed points of the pinned loop: there is
    # no pair of rows on which either direct or dual LLL can improve the L1 norm.
    assert _pinned_sublattice_basis_reduction(matrix(ZZ, 0, 3)) == matrix(ZZ, 0, 3)
    singleton = matrix(ZZ, [[-7, 2, 0]])
    assert _pinned_sublattice_basis_reduction(singleton) == singleton

    basis = matrix(ZZ, [[10, 1, 0], [9, 1, 1], [1, 0, 8]])
    reduced = _pinned_sublattice_basis_reduction(basis)
    assert reduced.base_ring() is ZZ
    assert reduced == matrix(ZZ, [[-1, 0, 1], [-1, -1, 0], [3, -3, 3]])
    change = reduced * basis.inverse()
    assert all(entry.denominator() == 1 for entry in change.list())
    assert abs(change.det()) == 1
    assert sum(abs(entry) for entry in reduced.list()) <= sum(abs(entry) for entry in basis.list())
    assert change * basis == reduced
    assert _pinned_sublattice_basis_reduction(reduced) == reduced

    # The kernel supplied by NullspaceIntVect is rectangular in ambient
    # coordinates; exercise the same reducer on a rank-three sublattice of Z^4.
    rectangular = matrix(ZZ, [[10, 1, 0, 2], [9, 1, 1, 0], [1, 0, 8, 1]])
    rectangular_reduced = _pinned_sublattice_basis_reduction(rectangular)
    assert rectangular_reduced == matrix(ZZ, [[-1, 0, 1, -2], [2, 0, 7, 3], [8, 1, 2, -2]])
    rectangular_change = rectangular_reduced * rectangular.transpose() * (rectangular * rectangular.transpose()).inverse()
    assert all(entry.denominator() == 1 for entry in rectangular_change.list())
    assert abs(rectangular_change.det()) == 1
    assert rectangular_change * rectangular == rectangular_reduced
    assert _pinned_sublattice_basis_reduction(rectangular_reduced) == rectangular_reduced


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



def test_pinned_sublattice_reduction_on_deterministic_full_rank_samples() -> None:
    """Random-looking small bases preserve their exact row lattices and stop."""
    import random

    from sage.all import ZZ, matrix

    randomizer = random.Random(150)
    for case in range(12):
        rows = 2 + case % 3
        columns = rows + case % 2
        while True:
            basis = matrix(
                ZZ,
                rows,
                columns,
                [randomizer.randint(-8, 8) for _ in range(rows * columns)],
            )
            if basis.rank() == rows:
                break
        reduced = _pinned_sublattice_basis_reduction(basis)
        change = reduced * basis.transpose() * (basis * basis.transpose()).inverse()
        assert all(entry.denominator() == 1 for entry in change.list())
        assert abs(change.det()) == 1
        assert change * basis == reduced
        assert sum(abs(entry) for entry in reduced.list()) <= sum(
            abs(entry) for entry in basis.list()
        )
        assert _pinned_sublattice_basis_reduction(reduced) == reduced


def test_primitive_integral_direction_matches_fraction_multiplier() -> None:
    """The scale is RemoveFractionPlusCoef's TheMult, not its reciprocal."""
    from sage.all import QQ, ZZ, vector

    for coordinates, expected_scale, expected_primitive in (
        ((QQ(2) / 3, QQ(4) / 3), QQ(3) / 2, (1, 2)),
        ((QQ(2), QQ(0)), QQ(1) / 2, (1, 0)),
        ((-QQ(3) / 2, QQ(9) / 4), QQ(4) / 3, (-2, 3)),
        ((QQ(0), QQ(0), QQ(5) / 7), QQ(7) / 5, (0, 0, 1)),
    ):
        scale, primitive = _primitive_integral_direction(coordinates)
        assert scale == expected_scale
        assert primitive == expected_primitive
        assert vector(QQ, coordinates) == vector(ZZ, primitive).change_ring(QQ) / scale
    for invalid in ((), (QQ(0),), (QQ(0), QQ(0))):
        try:
            _primitive_integral_direction(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError(f"zero-dimensional or zero direction accepted: {invalid}")
    try:
        _primitive_integral_direction((QQ(1), QQ(2)), expected_rank=3)
    except ValueError:
        pass
    else:
        raise AssertionError("a direction with the wrong lattice rank was accepted")


def test_ambient_vectors_reconstructs_exact_shell_and_filters_squares() -> None:
    """The kernel-coordinate candidates reconstruct to the prescribed fibre."""
    from types import SimpleNamespace

    from sage.all import ZZ, matrix, vector

    lorentzian = matrix(ZZ, [[1, 0, 0], [0, -1, 0], [0, 0, -1]])
    enumerator = SimpleNamespace(
        partner_row=vector(ZZ, (1, 0, 0)),
        inclusion_rows=matrix(ZZ, [[0, 1, 0], [0, 0, 1]]),
        gram=lorentzian,
        timelike_covector=vector(ZZ, (1, 0, 0)).column(),
        d=ZZ.one(),
        base_bound=ZZ(5),
        scaled_max=ZZ(2),
        allowed_squares={ZZ.zero()},
    )
    close_coordinates = ((1, 0), (0, 1), (0, 0))
    result = _PositiveVectorEnumerator._ambient_vectors(enumerator, 1, close_coordinates, (5, 5, 4))
    assert result == ((1, 1, 0), (1, 0, 1))
    for row in result:
        ambient = vector(ZZ, row)
        assert ambient * lorentzian * ambient.column() == 0
        assert (ambient * enumerator.timelike_covector)[0] == 1

    # At multiplier two, retain both the isotropic shell and an explicitly
    # allowed positive square, while discarding the square-three candidate.
    enumerator.scaled_max = ZZ(2)
    enumerator.allowed_squares = {ZZ.zero(), ZZ(4)}
    result = _PositiveVectorEnumerator._ambient_vectors(
        enumerator, 2, ((2, 0), (0, 2), (0, 0), (1, 0)), (20, 20, 16, 17)
    )
    assert result == ((2, 2, 0), (2, 0, 2), (2, 0, 0))
    assert all(int(vector(ZZ, row) * lorentzian * vector(ZZ, row)) in {0, 4} for row in result)

    try:
        _PositiveVectorEnumerator._ambient_vectors(enumerator, 2, ((0, 0),), ())
    except ValueError as error:
        assert "exactly one squared distance" in str(error)
    else:
        raise AssertionError("a CVP result without its paired squared distance was accepted")

    try:
        _PositiveVectorEnumerator._ambient_vectors(enumerator, 2, ((QQ(1) / 2, 0),), (QQ.zero(),))
    except ArithmeticError as error:
        assert "nonintegral lattice coordinates" in str(error)
    else:
        raise AssertionError("nonintegral CVP lattice coordinates were silently truncated")

    try:
        _PositiveVectorEnumerator._ambient_vectors(enumerator, 2, ((3, 0),), (20,))
    except ArithmeticError as error:
        assert "disagrees with its exact CVP squared distance" in str(error)
    else:
        raise AssertionError("a reconstructed candidate inconsistent with its CVP distance was accepted")


def test_rational_direction_pairing_bound_scales_linearly() -> None:
    """A bilinear pairing bound scales linearly with the direction."""
    from sage.all import QQ

    assert _integral_direction_pairing_bound(QQ(2) / 3, QQ(3) / 2) == QQ(1)
    assert _integral_direction_pairing_bound(QQ(7), QQ(1)) == QQ(7)
    assert _integral_direction_pairing_bound(QQ(5) / 4, QQ(4) / 3) == QQ(5) / 3
    # After lowering, MaxScal belongs to the integral primitive vector and
    # the source advances it by one, independently of the initial scale.
    assert _next_integral_pairing_bound(QQ(5) / 3) == QQ(8) / 3
    assert _next_integral_pairing_bound(QQ(3) / 2) == QQ(5) / 2


def test_initial_vector_search_advances_integral_bound_without_rescanning() -> None:
    """MaxScal advances by one and each multiplier is visited at most once."""
    from sage.all import QQ, ZZ, matrix
    from sage_indefinite_port.indefinite import lorentzian_cells as cells

    visited = []
    bounds = []

    class Enumerator:
        def __init__(self, lattice, direction, max_scal, allowed_squares):
            assert tuple(direction) == (ZZ(1), ZZ(1))
            assert max_scal == 2
            self.scale = QQ(1)
            self.d = ZZ(2)
            self.scaled_max = cells._integral_direction_pairing_bound(max_scal, self.scale)

        def max_multiplier(self):
            bounds.append(self.scaled_max)
            return int(self.scaled_max / self.d)

        def vectors_at_multiplier(self, multiplier):
            visited.append((self.scaled_max, multiplier))
            return ((1, 0),) if multiplier == 2 else ()

    class TimelikeLattice:
        def module_rank(self):
            return 2

        def is_isotropic(self):
            return True

    original_gram = cells._gram_matrix
    original_enumerator = cells._PositiveVectorEnumerator
    try:
        cells._gram_matrix = lambda lattice: matrix(ZZ, [[1, 0], [0, 1]])
        cells._PositiveVectorEnumerator = Enumerator
        result = cells._search_initial_vectors(TimelikeLattice(), (QQ(1) / 2, QQ(1) / 2), {ZZ.zero()})
    finally:
        cells._gram_matrix = original_gram
        cells._PositiveVectorEnumerator = original_enumerator

    assert result == ((1, 0),)
    # The source initializes MaxScal to q of the integral primitive
    # direction: q((1,1))=2, not q((1/2,1/2))=1/2. The direction and bound
    # are lowered together before constructing the enumerator; its scale is
    # then one. The bound advances by one and multiplier one is not rescanned.
    assert visited == [(QQ(2), 1), (QQ(4), 2)]
    assert bounds == [QQ(2), QQ(3), QQ(4)]


def test_isotropic_multiplier_uses_dual_lll_fixed_distance_solver() -> None:
    """The production isotropic multiplier path calls exact sphere enumeration."""
    from sage.all import QQ, ZZ, matrix, vector
    from sage_indefinite_port.indefinite import lorentzian_cells as cells

    class CVPEngine:
        engine_gram_q = matrix(QQ, [[2, 1], [1, 2]])

        def _engine_target_coordinates(self, target):
            return tuple(target)

        def _engine_scalar(self, scalar):
            return QQ(scalar)

        def close_vector_coordinates(self, target, bound):
            raise AssertionError("isotropic mode must not enumerate the full ellipsoid")

    enumerator = object.__new__(cells._PositiveVectorEnumerator)
    enumerator.rational_field = object()
    enumerator.lattice = object()
    enumerator.base_target_engine_coordinates = (QQ(1) / 2, QQ(0))
    enumerator.base_bound = QQ(1) / 2
    enumerator.cvp_engine = CVPEngine()
    enumerator.enumeration_transform = matrix(ZZ, [[1, 0], [-1, 1]])
    enumerator.enumeration_order_transform = matrix(ZZ, [[1, 1], [0, 1]])
    enumerator.partner_row = vector(ZZ, (1, 1, 0))
    enumerator.inclusion_rows = matrix(ZZ, [[0, 1, 0], [0, 0, 1]])
    enumerator.gram = matrix(ZZ, [[2, 0, 0], [0, -2, 0], [0, 0, -2]])
    enumerator.timelike_covector = vector(ZZ, (2, 0, 0)).column()
    enumerator.d = ZZ(2)
    enumerator.scaled_max = ZZ(4)
    enumerator.allowed_squares = cells.ZERO_SQUARES
    calls = []
    original_sphere = cells._fixed_distance_dual_lll_solutions
    original_element = cells._owned_engine_element
    original_rational = cells._owned_rational
    original_ambient = enumerator._ambient_vectors
    try:
        cells._fixed_distance_dual_lll_solutions = lambda *args: (calls.append(args) or (((0, 0), QQ(2)),))
        cells._owned_engine_element = lambda field, value: value
        cells._owned_rational = lambda lattice, value, field: value
        enumerator._ambient_vectors = lambda multiplier, coordinates, squares: (multiplier, coordinates, squares)
        assert enumerator.vectors_at_multiplier(2) == (2, ((0, 0),), (QQ(2),))
    finally:
        cells._fixed_distance_dual_lll_solutions = original_sphere
        cells._owned_engine_element = original_element
        cells._owned_rational = original_rational
        enumerator._ambient_vectors = original_ambient
    assert calls == [(CVPEngine.engine_gram_q, (QQ(1), QQ(0)), QQ(2), enumerator.enumeration_transform, enumerator.enumeration_order_transform)]


def test_fixed_distance_coordinate_solutions_match_exact_brute_force() -> None:
    """The isotropic shell enumerator returns the exact rational-center sphere."""
    from itertools import product
    from sage.all import QQ, RDF, ZZ, matrix

    assert _fixed_distance_coordinate_solutions(matrix(QQ, 0, 0), (), QQ(0)) == (((), QQ(0)),)
    assert _fixed_distance_coordinate_solutions(matrix(QQ, 0, 0), (), QQ(1)) == ()
    # A nonsquare rational squared radius has no rational square root, hence
    # no integer endpoint in rank one.
    assert _fixed_distance_coordinate_solutions(matrix(QQ, [[1]]), (QQ(0),), QQ(2)) == ()

    cases = (
        (matrix(QQ, [[1]]), (QQ(1) / 2,), QQ(1) / 4),
        (matrix(QQ, [[1]]), (QQ(1) / 3,), QQ(4) / 9),
        (matrix(QQ, [[2, 1], [1, 3]]), (QQ(1) / 2, QQ(1) / 3), QQ(1) / 2),
        (matrix(QQ, [[2, 1], [1, 2]]), (QQ(1) / 2, QQ(1) / 2), QQ(1) / 2),
        (matrix(QQ, [[1, 0], [0, 2]]), (QQ(0), QQ(0)), QQ(1)),
    )
    for gram, target, bound in cases:
        actual = _fixed_distance_coordinate_solutions(gram, target, bound)
        actual_coordinates = {coordinates for coordinates, distance in actual if distance == bound}
        assert len(actual_coordinates) == len(actual)
        assert all(distance == bound for _coordinates, distance in actual)
        expected = set()
        for coordinates in product(range(-5, 6), repeat=gram.nrows()):
            delta = tuple(QQ(coordinates[i]) - target[i] for i in range(len(target)))
            square = sum((delta[i] * gram[i, j] * delta[j] for i in range(len(delta)) for j in range(len(delta))), QQ(0))
            if square == bound:
                expected.add(tuple(coordinates))
        assert actual_coordinates == expected

    # The innermost coordinate is solved directly, so this huge radius must
    # still produce only its two exact endpoints rather than scanning 2e10
    # integers in a one-dimensional enclosing interval.
    huge_radius = QQ(10) ** 10
    huge_sphere = _fixed_distance_coordinate_solutions(
        matrix(QQ, [[1]]), (QQ(0),), huge_radius ** 2
    )
    assert {coordinates for coordinates, _distance in huge_sphere} == {
        (ZZ(-huge_radius),), (ZZ(huge_radius),)
    }

    for gram, target, bound in (
        (matrix(QQ, [[1, 2], [0, 1]]), (QQ(0), QQ(0)), QQ(1)),
        (matrix(QQ, [[1, 0], [0, 1]]), (QQ(0),), QQ(1)),
        (matrix(QQ, [[1]]), (0.5,), QQ(1)),
        (matrix(QQ, [[1]]).change_ring(RDF), (RDF(0),), QQ(1)),
        (matrix(QQ, [[1]]), (QQ(0),), 1.0),
        (matrix(QQ, [[1]]), (QQ(0),), QQ(-1)),
    ):
        try:
            _fixed_distance_coordinate_solutions(gram, target, bound)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid fixed-distance sphere input was accepted")


def test_fixed_distance_dual_lll_frame_reconstructs_exact_sphere() -> None:
    """Exact sphere pruning uses the reduced dual-LLL Gram and reconstructs the original frame."""
    from itertools import product
    from sage.all import QQ, RDF, ZZ, diagonal_matrix, identity_matrix, matrix, vector

    gram = diagonal_matrix(QQ, [1, 2, 3])
    # The shear changes this center by an integral vector in the first
    # reduced coordinate, so a mistaken untransformed target still fails.
    target = (QQ(0), QQ(1) / 2, QQ(0))
    bound = QQ(1) / 2
    order = matrix(ZZ, [[1, 2, 0], [0, 1, 1], [0, 0, 1]])
    dual = order.inverse().transpose()
    actual = _fixed_distance_dual_lll_solutions(gram, target, bound, dual, order)
    got = {coordinates for coordinates, distance in actual if distance == bound}
    assert len(got) == len(actual)
    expected = set()
    for coordinates in product(range(-4, 5), repeat=3):
        delta = tuple(QQ(coordinates[i]) - target[i] for i in range(3))
        if sum((delta[i] * gram[i, j] * delta[j] for i in range(3) for j in range(3)), QQ(0)) == bound:
            expected.add(tuple(ZZ(value) for value in coordinates))
    assert got == expected == {
        (ZZ(0), ZZ(0), ZZ(0)),
        (ZZ(0), ZZ(1), ZZ(0)),
    }
    assert tuple(actual) == tuple(sorted(actual, key=lambda item: tuple((order * vector(ZZ, item[0]).column())[i, 0] for i in reversed(range(3)))))

    # Deterministic nontrivial Gram/target/frame combinations check the
    # congruence, target transport and reconstruction together.
    for gram, target, bound, order in (
        (matrix(QQ, [[2, 1], [1, 3]]), (QQ(1) / 2, QQ(1) / 3), QQ(1) / 2, matrix(ZZ, [[1, 2], [0, 1]])),
        (diagonal_matrix(QQ, [1, 2, 3]), (QQ(1) / 2, QQ(0), QQ(0)), QQ(1) / 4, matrix(ZZ, [[1, 0, 0], [1, 1, 0], [0, 0, 1]])),
        (matrix(QQ, [[2, 1], [1, 2]]), (QQ(0), QQ(1) / 2), QQ(1) / 2, matrix(ZZ, [[1, -1], [0, 1]])),
    ):
        dual = order.inverse().transpose()
        actual = _fixed_distance_dual_lll_solutions(gram, target, bound, dual, order)
        got = {coordinates for coordinates, distance in actual if distance == bound}
        expected = set()
        for coordinates in product(range(-6, 7), repeat=gram.nrows()):
            delta = tuple(QQ(coordinates[i]) - target[i] for i in range(len(target)))
            if sum((delta[i] * gram[i, j] * delta[j] for i in range(len(delta)) for j in range(len(delta))), QQ(0)) == bound:
                expected.add(tuple(ZZ(value) for value in coordinates))
        assert got == expected
        assert len(actual) == len(got) and all(distance == bound for _, distance in actual)

    import random

    rng = random.Random(1207)
    for rank in range(1, 4):
        for _ in range(10):
            raw = matrix(ZZ, rank, rank, [rng.randint(-2, 2) for _ in range(rank * rank)])
            gram = (raw * raw.transpose() + identity_matrix(ZZ, rank)).change_ring(QQ)
            target = tuple(QQ(rng.randint(-2, 2)) / rng.randint(1, 3) for _ in range(rank))
            bound = QQ(rng.randint(0, 5)) / rng.randint(1, 3)
            order = identity_matrix(ZZ, rank)
            if rank > 1:
                order[0, 1] = rng.choice((-2, -1, 1, 2))
            if rank > 2:
                order[1, 2] = rng.choice((-1, 1))
            dual = order.inverse().transpose()
            actual = _fixed_distance_dual_lll_solutions(gram, target, bound, dual, order)
            got = {coordinates for coordinates, distance in actual if distance == bound}
            expected = set()
            for coordinates in product(range(-8, 9), repeat=rank):
                delta = tuple(QQ(coordinates[i]) - target[i] for i in range(rank))
                square = sum((delta[i] * gram[i, j] * delta[j] for i in range(rank) for j in range(rank)), QQ(0))
                if square == bound:
                    expected.add(tuple(ZZ(value) for value in coordinates))
            assert expected <= got
            if all(all(abs(value) <= 8 for value in coordinates) for coordinates in got):
                assert got == expected
            assert len(actual) == len(got) and all(distance == bound for _, distance in actual)

    for bad_gram, bad_target, bad_bound, bad_dual, bad_order, exception in (
        (gram, target, bound, identity_matrix(ZZ, 2), order, ValueError),
        (gram, target, bound, identity_matrix(ZZ, 3), order, ArithmeticError),
        (gram, (0.5, QQ(0), QQ(0)), bound, dual, order, ValueError),
        (gram.change_ring(RDF), target, bound, dual, order, ValueError),
        (gram, target, bound, dual.change_ring(RDF), order, ValueError),
        (gram, target, bound, matrix(QQ, [[QQ(1) / 2, 0, 0], [0, 1, 0], [0, 0, 1]]), order, ValueError),
        (gram, target, 0.5, dual, order, ValueError),
        (gram, target, QQ(-1), dual, order, ValueError),
    ):
        try:
            _fixed_distance_dual_lll_solutions(bad_gram, bad_target, bad_bound, bad_dual, bad_order)
        except exception:
            pass
        else:
            raise AssertionError("incompatible dual-LLL transforms were accepted")


def test_fixed_distance_basis_inverse_is_exact_and_two_sided() -> None:
    """A cached inverse is certified before it transports the CVP center."""
    from sage.all import QQ, RDF, identity_matrix, matrix

    basis = matrix(QQ, [[1, 2], [0, 1]])
    inverse = matrix(QQ, [[1, -2], [0, 1]])
    _validate_fixed_distance_basis_inverse(basis, inverse)
    for bad_basis, bad_inverse, exception in (
        (basis, identity_matrix(QQ, 1), ValueError),
        (basis, matrix(QQ, [[1, 0], [0, 1]]), ArithmeticError),
        (basis.change_ring(RDF), inverse, ValueError),
    ):
        try:
            _validate_fixed_distance_basis_inverse(bad_basis, bad_inverse)
        except exception:
            pass
        else:
            raise AssertionError("an invalid fixed-distance basis inverse was accepted")


def test_positive_ldl_factorization_is_exact_and_rejects_nonpositive_pivots() -> None:
    """The fixed-distance sphere uses an exact positive-definite LDL certificate."""
    from sage.all import QQ, matrix, diagonal_matrix

    forms = (
        matrix(QQ, [[1]]),
        matrix(QQ, [[2, 1], [1, 3]]),
        matrix(QQ, [[QQ(5) / 2, QQ(1) / 3], [QQ(1) / 3, 2]]),
        matrix(QQ, [[4, 1, 2], [1, 3, 0], [2, 0, 5]]),
        diagonal_matrix(QQ, [QQ(1) / 2, 3, 7]),
    )
    for gram in forms:
        lower, diagonal = _positive_ldl_factorization(gram)
        lower_matrix = matrix(QQ, lower)
        diagonal_matrix_ = diagonal_matrix(QQ, diagonal)
        assert lower_matrix.is_lower_triangular()
        assert all(pivot > 0 for pivot in diagonal)
        assert lower_matrix * diagonal_matrix_ * lower_matrix.transpose() == gram
        _validate_fixed_distance_ldl_factorization(gram, (lower, diagonal))

    gram = matrix(QQ, [[2, 1], [1, 3]])
    valid_lower, valid_diagonal = _positive_ldl_factorization(gram)
    for bad_factorization in (
        ([[1, 1], [0, 1]], valid_diagonal),
        (valid_lower, [valid_diagonal[0], -valid_diagonal[1]]),
        (valid_lower, [valid_diagonal[0], valid_diagonal[1] + 1]),
        ([[1]], [QQ(1)]),
        ([[1.0, 0.0], [0.0, 1.0]], [1.0, 1.0]),
    ):
        try:
            _validate_fixed_distance_ldl_factorization(gram, bad_factorization)
        except (ValueError, ArithmeticError):
            pass
        else:
            raise AssertionError(f"invalid cached LDL factorization accepted: {bad_factorization}")

    for gram in (
        matrix(QQ, 0, 0),
        matrix(QQ, [[1, 2], [0, 1]]),
    ):
        try:
            _positive_ldl_factorization(gram)
        except ValueError as error:
            assert "nonempty symmetric Gram matrix" in str(error)
        else:
            raise AssertionError(f"invalid Gram matrix accepted: {gram}")

    for gram in (
        matrix(QQ, [[0]]),
        matrix(QQ, [[1, 2], [2, 1]]),
        matrix(QQ, [[1, 0], [0, -1]]),
    ):
        try:
            _positive_ldl_factorization(gram)
        except ArithmeticError as error:
            assert "positive-definite Gram matrix" in str(error)
        else:
            raise AssertionError(f"non-positive-definite Gram matrix accepted: {gram}")


def test_multiplier_scales_cvp_center_and_squared_distance_bound() -> None:
    """The pinned CVP receives m*w and the quadratic radius m^2*base_bound."""
    from sage.all import QQ, ZZ, identity_matrix, vector
    from sage_indefinite_port.indefinite import lorentzian_cells as cells

    class CVPEngine:
        def __init__(self):
            self.call = None

        def close_vector_coordinates(self, target, bound):
            self.call = (target, bound)
            # Deliberately return reverse source order. A nontrivial unimodular
            # transform orders (-2,1) before (-3,0), and each exact squared
            # distance must stay attached to its own coordinate row.
            return (((-3, 0), QQ.zero()), ((-2, 1), QQ(2)))

    enumerator = object.__new__(cells._PositiveVectorEnumerator)
    engine = CVPEngine()
    enumerator.rational_field = object()
    enumerator.lattice = object()
    enumerator.base_target_engine_coordinates = (-QQ.one(), QQ.zero())
    enumerator.base_bound = QQ(5)
    enumerator.cvp_engine = engine
    enumerator.enumeration_order_transform = matrix(ZZ, [[1, 0], [-2, 1]])
    enumerator.partner_row = vector(ZZ, (1, 1, 0))
    enumerator.inclusion_rows = matrix(ZZ, [[0, 1, 0], [0, 0, 1]])
    enumerator.gram = matrix(ZZ, [[5, 0, 0], [0, -1, 0], [0, 0, -1]])
    enumerator.timelike_covector = vector(ZZ, (5, 0, 0)).column()
    enumerator.d = ZZ(5)
    enumerator.scaled_max = ZZ(15)
    enumerator.allowed_squares = {ZZ(43), ZZ(45)}
    original_element = cells._owned_engine_element
    original_rational = cells._owned_rational
    try:
        cells._owned_engine_element = lambda field, value: value
        cells._owned_rational = lambda lattice, value, field: value
        assert enumerator.vectors_at_multiplier(3) == ((3, 1, 1), (3, 0, 0))
    finally:
        cells._owned_engine_element = original_element
        cells._owned_rational = original_rational

    assert engine.call == ((-QQ(3), QQ.zero()), QQ(45))
    engine.close_vector_coordinates = lambda target, bound: (((0, 0), QQ.one()),)
    try:
        enumerator.vectors_at_multiplier(2)
    except ArithmeticError as error:
        assert "disagrees with its exact CVP squared distance" in str(error)
    else:
        raise AssertionError("an inconsistent CVP distance was accepted")
    for invalid_distance in (QQ(-1), QQ(46)):
        engine.close_vector_coordinates = lambda target, bound, distance=invalid_distance: (
            ((0, 0), distance),
        )
        try:
            enumerator.vectors_at_multiplier(2)
        except ArithmeticError as error:
            assert "outside its exact squared-distance bound" in str(error)
        else:
            raise AssertionError(f"an out-of-bound CVP distance was accepted: {invalid_distance}")
