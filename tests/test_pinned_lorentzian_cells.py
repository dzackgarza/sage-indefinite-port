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
    _pinned_nullspace_matrix,
    _pinned_diagonalize_symmetric_matrix,
    _pinned_independent_column_indices,
    _outside_ray_directions,
    _quadratic_descent_prefix,
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


def test_source_random_perturbations_are_the_pinned_unimodular_matrices() -> None:
    """Check both GetRandomMatrixPerturbation branches, including diagonal shear."""
    from sage.all import QQ, matrix

    signed_permutation = _source_permutation_sign_matrix([2, 0, 1], (-1, 1, -1))
    assert signed_permutation == matrix(QQ, [[0, 0, -1], [1, 0, 0], [0, -1, 0]])
    assert abs(signed_permutation.det()) == 1

    for row, column, sign in ((0, 1, 1), (2, 0, -1), (1, 1, 1)):
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
