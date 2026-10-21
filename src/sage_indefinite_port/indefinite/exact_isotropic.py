r"""FindIsotropicExact of polyhedral_common src_isotropy (pinned 1592b246), over Sage number theory.

An isotropic vector of a nondegenerate rational quadratic form, of every rank:

- rank 1 and 2: Isotropic.h:203-249, directly;
- rank 3: the Legendre equation, Legendre_equation.h (reduction to coprime
  square-free coefficients, Legendre's solvability criterion, Lagrange descent);
- rank 4: Simon and Pilnikova's construction, Quaternary.h (local solvability
  from Hilbert symbols, an auxiliary t = p0 q splitting the form into two
  ternary equations);
- rank 5 and more: FifthAndHigherOrder.h (a five-dimensional subspace of
  mixed signature, then the fifth-order construction from one ternary and one
  quaternary equation).

The generic number theory those routines call is Sage's: factorization,
primality, the Hilbert symbol (Quaternary.h ``hilbert_symbol`` and
``Padic_isotropy_ternary`` decide the p-adic isotropy of the ternary form
(a, b, -1), which is (a, b)_p = 1), p-adic square classes
(NumberTheoryPadic.h ``Padic_is_square``), square roots modulo prime powers
and the Chinese remainder theorem (quadratic_residue.h).  Where the source
returns a specific root, the same root is taken: ``find_quadratic_residue``
returns the least root modulo each prime power, combined by CRT.
"""

from math import gcd, lcm

import numpy
from sage.arith.misc import CRT_list, factor, hilbert_symbol, is_prime
from sage.arith.misc import jacobi_symbol as sage_jacobi_symbol
from sage.matrix.constructor import identity_matrix, matrix
from sage.modules.free_module_element import vector
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage_indefinite_port.indefinite.indefinite_reduction import _indefinite_reduction

# ---------------------------------------------------------------------------
# Matrix helpers (SignatureSymmetric.h, MAT_MatrixInt.h)
# ---------------------------------------------------------------------------


def _remove_fraction_vector(values):
    r"""RemoveFractionVector (MAT_MatrixInt.h:360-437): clear denominators, divide by the positive content."""
    rationals = [QQ(value) for value in values]
    common = lcm(*(int(value.denominator()) for value in rationals))
    integers = [int(value * common) for value in rationals]
    content = gcd(*integers)
    if content == 0:
        return integers
    return [entry // content for entry in integers]


def _scaled_inverse(transformation):
    r"""ScaledInverse (MAT_MatrixInt.h:446-456): the inverse with its fractions removed."""
    inverse = transformation.change_ring(QQ).inverse()
    rows = _remove_fraction_vector(inverse.list())
    return matrix(ZZ, inverse.nrows(), inverse.ncols(), rows)


def _find_non_isotrope_vector(form):
    r"""FindNonIsotropeVector (SignatureSymmetric.h:22-42)."""
    size = form.nrows()
    for i in range(size):
        if form[i, i] != 0:
            return vector(QQ, [1 if j == i else 0 for j in range(size)])
    for i in range(size - 1):
        for j in range(i + 1, size):
            if form[i, j] != 0:
                return vector(QQ, [1 if k in (i, j) else 0 for k in range(size)])
    raise ArithmeticError("a nonzero symmetric form has a non-isotropic basis or sum vector")


def _diagonalize_nondegenerate(form):
    r"""DiagonalizeNonDegenerateSymmetricMatrix (SignatureSymmetric.h:60-125): an orthogonal basis by successive perpendiculars."""
    form = form.change_ring(QQ)
    size = form.nrows()
    vectors = []
    for i in range(size):
        if i == 0:
            orthogonal = identity_matrix(QQ, size)
        else:
            orthogonal = (form * matrix(QQ, vectors).transpose()).left_kernel().basis_matrix()
        restricted = orthogonal * form * orthogonal.transpose()
        vectors.append(orthogonal.transpose() * _find_non_isotrope_vector(restricted))
    basis = matrix(QQ, vectors)
    return basis, basis * form * basis.transpose()


def _diagonalize(form):
    r"""DiagonalizeSymmetricMatrix (SignatureSymmetric.h:141-202): the radical first, then the nondegenerate part.

    NullspaceReduction completes a basis of the radical by the unit vectors of
    the non-pivot columns of its echelon basis (ColumnReductionSet).
    """
    form = form.change_ring(QQ)
    size = form.nrows()
    kernel = form.left_kernel().basis_matrix()
    pivots = set(kernel.echelon_form().pivots())
    complement = [vector(QQ, [1 if k == j else 0 for k in range(size)]) for j in range(size) if j not in pivots]
    first = matrix(QQ, [*kernel.rows(), *complement]) if size else identity_matrix(QQ, 0)
    radical_rank = kernel.nrows()
    nondegenerate = (first * form * first.transpose()).matrix_from_rows_and_columns(range(radical_rank, size), range(radical_rank, size))
    inner, _ = _diagonalize_nondegenerate(nondegenerate)
    second = identity_matrix(QQ, size)
    for i in range(size - radical_rank):
        for j in range(size - radical_rank):
            second[radical_rank + i, radical_rank + j] = inner[i, j]
    transform = second * first
    return transform, transform * form * transform.transpose()


def _signature_counts(form):
    r"""DiagonalizeSymmetricMatrixSymbolic (SignatureSymmetric.h:203-240): the numbers of positive and negative squares."""
    positive, negative, _zero = QuadraticForm(QQ, form.change_ring(QQ)).signature_vector()
    return int(positive), int(negative)


def _reduced_diagonal(form):
    r"""get_reduced_diagonal (Legendre_equation.h:390-411): the orthogonal basis and its primitive integral diagonal."""
    transform, reduced = _diagonalize_nondegenerate(form)
    diagonal = _remove_fraction_vector([reduced[i, i] for i in range(reduced.nrows())])
    if any(entry == 0 for entry in diagonal):
        raise ArithmeticError("a nondegenerate form has a nonzero diagonal")
    return transform, diagonal


# ---------------------------------------------------------------------------
# Number theory (quadratic_residue.h, NumberTheoryPadic.h)
# ---------------------------------------------------------------------------


def _factor_map(value):
    return {int(prime): int(exponent) for prime, exponent in factor(abs(int(value)))} if abs(int(value)) > 1 else {}


def _jacobi_filter(value, modulus):
    r"""compute_jacobi_symbol (quadratic_residue.h:256-329): the Jacobi symbol of value over the odd part of modulus is +1.

    The source drops the powers of two of the modulus and evaluates the
    Jacobi symbol over its odd part; a residue modulo that part is a residue
    modulo the whole, so -1 excludes a square root.
    """
    odd = int(modulus)
    while odd % 2 == 0:
        odd //= 2
    if odd == 1 or int(value) % odd == 0:
        return True
    return sage_jacobi_symbol(int(value) % odd, odd) == 1


def _least_square_root(value, modulus):
    r"""find_quadratic_residue_exhaustive (quadratic_residue.h:72-147): the least x >= 0 with x^2 = value mod modulus."""
    roots = Mod(int(value), int(modulus)).sqrt(all=True)
    return min(int(root) for root in roots) if roots else None


def _find_quadratic_residue(value, modulus):
    r"""find_quadratic_residue (quadratic_residue.h:149-209): the least root modulo each prime power, combined by CRT."""
    modulus = abs(int(modulus))
    if modulus == 1:
        return 0
    roots, moduli = [], []
    for prime, exponent in _factor_map(modulus).items():
        power = prime**exponent
        root = _least_square_root(value, power)
        if root is None:
            return None
        roots.append(root)
        moduli.append(power)
    return int(CRT_list(roots, moduli)) % modulus


def _is_quadratic_residue(value, modulus):
    r"""is_quadratic_residue (quadratic_residue.h:366-372)."""
    return Mod(int(value), abs(int(modulus))).is_square()


def _is_quadratic_residue_map(value, factor_map):
    r"""is_quadratic_residue_map (quadratic_residue.h:331-338): a square root modulo every listed prime power."""
    return all(Mod(int(value), prime**exponent).is_square() for prime, exponent in factor_map.items())


def _is_square_free(value):
    return all(exponent < 2 for exponent in _factor_map(value).values())


def _separate_square_factor(value):
    r"""separate_square_factor (Legendre_equation.h:433-453): value = sign * sqr^2 * nosqr with nosqr square-free."""
    square, free = 1, (1 if value > 0 else -1)
    for prime, exponent in _factor_map(value).items():
        free *= prime ** (exponent % 2)
        square *= prime ** (exponent // 2)
    return square, free


# ---------------------------------------------------------------------------
# Ternary forms: the Legendre equation (Legendre_equation.h)
# ---------------------------------------------------------------------------


_SHARED_PRIME_SLOT = {
    (True, False, False): "a",
    (False, True, True): "bc",
    (False, True, False): "b",
    (True, False, True): "ac",
    (False, False, True): "c",
    (True, True, False): "ab",
}
"""reduction_information (:304-332): where a prime of the odd parts goes, by which coefficients it divides."""


def _reduction_information(coefficients):
    r"""reduction_information (Legendre_equation.h:208-388).

    Returns (TransMat, aReduced, factorizations) with B Diag(aReduced) B = u Diag(a):
    even prime powers are moved into the other coordinates, then each prime
    shared by two coefficients is moved to the third.
    """
    a, b, c = (int(value) for value in coefficients)
    maps = [_factor_map(a), _factor_map(b), _factor_map(c)]
    transformation = [1, 1, 1]
    odd_sets = []
    for index, factor_map in enumerate(maps):
        odd = set()
        for prime, exponent in factor_map.items():
            if exponent % 2:
                odd.add(prime)
            if exponent // 2:
                for other in range(3):
                    if other != index:
                        transformation[other] *= prime ** (exponent // 2)
        odd_sets.append(odd)
    a_set, b_set, c_set = odd_sets
    primes = set(maps[0]) | set(maps[1]) | set(maps[2])
    products = {"a": 1, "b": 1, "c": 1, "ab": 1, "ac": 1, "bc": 1}
    product_maps = {key: {} for key in products}
    for prime in sorted(primes):
        membership = (prime in a_set, prime in b_set, prime in c_set)
        key = _SHARED_PRIME_SLOT.get(membership)
        if key is not None:
            products[key] *= prime
            product_maps[key][prime] = product_maps[key].get(prime, 0) + 1
    transformation[0] *= products["bc"]
    transformation[1] *= products["ac"]
    transformation[2] *= products["ab"]

    def sign(value):
        return 1 if value > 0 else -1

    reduced = [sign(a) * products["a"] * products["bc"], sign(b) * products["b"] * products["ac"], sign(c) * products["c"] * products["ab"]]
    factorizations = [{**product_maps["a"], **product_maps["bc"]}, {**product_maps["b"], **product_maps["ac"]}, {**product_maps["c"], **product_maps["ab"]}]
    return matrix(ZZ, 3, 3, [transformation[0], 0, 0, 0, transformation[1], 0, 0, 0, transformation[2]]), reduced, factorizations


def _ternary_is_solvable(reduced, factorizations):
    r"""determine_solvability_dim3 (Legendre_equation.h:97-206): Legendre's criterion for coprime square-free a, b, c."""
    if all(value > 0 for value in reduced) or all(value < 0 for value in reduced):
        return False
    a, b, c = reduced
    counterparts = (-b * c, -a * c, -a * b)
    if not all(_jacobi_filter(counterpart, abs(value)) for counterpart, value in zip(counterparts, reduced, strict=True)):
        return False
    return all(_is_quadratic_residue_map(counterpart, factor_map) for counterpart, factor_map in zip(counterparts, factorizations, strict=True))


def _lagrange_normal(diagonal):
    r"""get_lagrange_normal (Legendre_equation.h:504-571): Z^2 = a X^2 + b Y^2 from the coefficient of least absolute value."""
    least = min(range(3), key=lambda index: (abs(diagonal[index]), index))
    index_z, index_x, index_y = least, (least + 1) % 3, (least + 2) % 3
    value = diagonal[index_z]
    transformation = matrix(ZZ, 3, 3)
    transformation[index_x, 0] = value
    transformation[index_y, 1] = value
    transformation[index_z, 2] = 1
    return transformation, (-diagonal[index_x] * value, -diagonal[index_y] * value)


def _descent_operation(pair):
    r"""descent_operation (Legendre_equation.h:573-629): Lemma 6 of Lagrange's descent, |a| <= |b|."""
    a, b = pair
    root = _find_quadratic_residue(a, b)
    if root is None:
        raise ArithmeticError("Lagrange descent requires a to be a square modulo b")
    if 2 * root >= abs(b):
        root -= abs(b)
    square, free = _separate_square_factor((root * root - a) // b)
    transformation = matrix(ZZ, 3, 3)
    transformation[0, 0] = root
    transformation[0, 2] = 1
    transformation[1, 1] = b * square
    transformation[2, 0] = a
    transformation[2, 2] = root
    return transformation, (a, free)


def _ternary_via_lagrange(diagonal):
    r"""TernaryIsotropicViaLagrange (Legendre_equation.h:666-779)."""
    normal, pair = _lagrange_normal(diagonal)
    steps = []
    while True:
        a, b = pair
        if a == 1 or b == 1:
            # get_trivial_solution (:643-664).
            solution = vector(ZZ, [1, 0, 1]) if a == 1 else vector(ZZ, [0, 1, 1])
            for transformation in reversed(steps):
                solution = _scaled_inverse(transformation) * solution
            return normal * solution
        if abs(a) > abs(b):
            # switch_xy (:631-641).
            swap = matrix(ZZ, 3, 3, [0, 1, 0, 1, 0, 0, 0, 0, 1])
            pair = (b, a)
            steps.append(swap)
        transformation, pair = _descent_operation(pair)
        steps.append(transformation)


def ternary_isotropic_vector_diagonal(coefficients):
    r"""TernaryIsotropicVectorDiagonal (Legendre_equation.h:781-821): an isotropic vector of a1 x^2 + a2 y^2 + a3 z^2, or None."""
    transformation, reduced, factorizations = _reduction_information(coefficients)
    if not _ternary_is_solvable(reduced, factorizations):
        return None
    return transformation * _ternary_via_lagrange(reduced)


def ternary_isotropic_vector(form):
    r"""TernaryIsotropicVector (Legendre_equation.h:823-873)."""
    transform, diagonal = _reduced_diagonal(form)
    solution = ternary_isotropic_vector_diagonal(diagonal)
    return None if solution is None else transform.transpose() * solution.change_ring(QQ)


# ---------------------------------------------------------------------------
# Quaternary forms (Quaternary.h)
# ---------------------------------------------------------------------------


def _local_primes(coefficients):
    r"""get_local_primes (Quaternary.h:127-150): 2 and the primes dividing a coefficient."""
    return sorted({2} | {prime for value in coefficients for prime in _factor_map(value)})


def _quaternary_is_solvable(coefficients, primes):
    r"""determine_solvability_dim4 (Quaternary.h:152-284): Lemma 13 of Simon and Pilnikova at every local prime.

    The form is anisotropic at p exactly when -a1a2 / -a3a4 is a p-adic square
    and (a1, a2)_p = -(-a3, -a4)_p.
    """
    if all(value > 0 for value in coefficients) or all(value < 0 for value in coefficients):
        return False
    a1, a2, a3, a4 = coefficients
    ratio = QQ(-a1 * a2) / QQ(-a3 * a4)
    for prime in primes:
        if ratio.is_padic_square(prime) and hilbert_symbol(a1, a2, prime) == -hilbert_symbol(-a3, -a4, prime):
            return False
    return True


def _residue_classes(prime):
    r"""Padic_get_residue_classes (NumberTheoryPadic.h:386-417): representatives of Q_p^* / (Q_p^*)^2, in source order."""
    if prime == 2:
        return [entry for value in (1, 3, 5, 7) for entry in (value, 2 * value)]
    nonresidue = next(value for value in range(2, prime) if not _is_quadratic_residue(value, prime))
    return [entry for value in (1, nonresidue) for entry in (value, prime * value)]


def _possible_signs(values):
    r"""get_poss_signs (Quaternary.h:352-369): the signs of t for which sum a_i x_i^2 + t z^2 can be indefinite."""
    positive = sum(1 for value in values if value > 0)
    if positive == 0:
        return [1]
    if positive == len(values):
        return [-1]
    return [-1, 1]


def _quaternary_solution(coefficients, primes):
    r"""dim4_compute_solution (Quaternary.h:301-530): t = p0 q splits the form into two solvable ternary equations."""
    a1, a2, a3, a4 = coefficients

    # get_tp_classes (:301-350): t_p with (t, -a b)_p = (a, b)_p for (a1, a2) and (-a3, -a4).
    def matches(t, prime):
        return hilbert_symbol(t, -a1 * a2, prime) == hilbert_symbol(a1, a2, prime) and hilbert_symbol(t, -a3 * a4, prime) == hilbert_symbol(-a3, -a4, prime)

    classes = [next(t for t in _residue_classes(prime) if matches(t, prime)) for prime in primes]
    # get_q_val (:371-402).
    signs = [sign for sign in _possible_signs([-a1, -a2]) if sign in _possible_signs([a3, a4])]
    q = signs[0]
    for prime, class_representative in zip(primes, classes, strict=True):
        q *= prime ** int(ZZ(class_representative).valuation(prime))
    # dim4_pair_legendre_iterate_solution (:404-511).
    symbols12 = [hilbert_symbol(a1, a2, prime) for prime in primes]
    symbols34 = [hilbert_symbol(-a3, -a4, prime) for prime in primes]
    candidate = 2
    while True:
        if is_prime(candidate):
            t = candidate * q
            if all(
                hilbert_symbol(t, -a1 * a2, prime) == symbol12 and hilbert_symbol(t, -a3 * a4, prime) == symbol34
                for prime, symbol12, symbol34 in zip(primes, symbols12, symbols34, strict=True)
            ):
                first = ternary_isotropic_vector_diagonal([a1, a2, -t])
                if first is None:
                    raise ArithmeticError("the first ternary equation of the Simon-Pilnikova splitting must be solvable")
                if first[2] == 0:
                    return vector(ZZ, [first[0], first[1], 0, 0])
                second = ternary_isotropic_vector_diagonal([a3, a4, t])
                if second is None:
                    raise ArithmeticError("the second ternary equation of the Simon-Pilnikova splitting must be solvable")
                return vector(ZZ, [first[0] * second[2], first[1] * second[2], second[0] * first[2], second[1] * first[2]])
        candidate += 1


def quaternary_isotropic_vector_diagonal(coefficients):
    r"""QuaternaryIsotropicVectorDiagonal (Quaternary.h:547-566)."""
    coefficients = [int(value) for value in coefficients]
    primes = _local_primes(coefficients)
    if not _quaternary_is_solvable(coefficients, primes):
        return None
    return _quaternary_solution(coefficients, primes)


def quaternary_isotropic_vector(form):
    r"""QuaternaryIsotropicVector (Quaternary.h:568-606)."""
    transform, diagonal = _reduced_diagonal(form)
    solution = quaternary_isotropic_vector_diagonal(diagonal)
    return None if solution is None else transform.transpose() * solution.change_ring(QQ)


# ---------------------------------------------------------------------------
# Rank five and more (FifthAndHigherOrder.h, Positivity.h)
# ---------------------------------------------------------------------------


def _solve_fifth_equation(coefficients):
    r"""solve_fifth_equation (FifthAndHigherOrder.h:41-125): one ternary and one quaternary equation sharing t."""
    a = [int(value) for value in coefficients]
    signs = [sign for sign in _possible_signs([-a[0], -a[1]]) if sign in _possible_signs([a[2], a[3], a[4]])]
    sign = signs[0]
    t = sign
    while True:
        ternary = ternary_isotropic_vector_diagonal([a[0], a[1], -t])
        if ternary is not None:
            if ternary[2] == 0:
                return vector(ZZ, [ternary[0], ternary[1], 0, 0, 0])
            quaternary = quaternary_isotropic_vector_diagonal([a[2], a[3], a[4], t])
            if quaternary is not None:
                return vector(ZZ, [ternary[0] * quaternary[3], ternary[1] * quaternary[3], quaternary[0] * ternary[2], quaternary[1] * ternary[2], quaternary[2] * ternary[2]])
        t += sign


class _ApproxIterator:
    r"""ApproxIterator (Positivity.h:257-310): the integral approximants round(s V), s increasing."""

    def __init__(self, direction):
        groups = {}
        for index, value in enumerate(direction):
            if value != 0:
                groups.setdefault(1 / abs(value), []).append((index, 1 if value > 0 else -1))
        self.approximation = [0] * len(direction)
        self.steps = list(groups)
        self.next = [step / 2 for step in self.steps]
        self.coefficients = [groups[step] for step in self.steps]

    def increment(self):
        least = min(self.next)
        for position, value in enumerate(self.next):
            if value == least:
                self.next[position] += self.steps[position]
                for index, sign in self.coefficients[position]:
                    self.approximation[index] += sign
        return vector(ZZ, self.approximation)


def _positive_directions_diagonal(form):
    r"""GetPositiveDirections_diag (Positivity.h:138-158): the l1-normalized basis vectors of positive square."""
    transform, reduced = _diagonalize(form)
    return [transform.row(i) / sum(abs(entry) for entry in transform.row(i)) for i in range(form.nrows()) if reduced[i, i] > 0]


def _positive_directions_eigen(form):
    r"""GetPositiveDirections_eigen (Positivity.h:169-189).

    As the source does, row i of Eigen's eigenvector matrix is taken for each
    positive eigenvalue i.
    """
    values, vectors = numpy.linalg.eigh(numpy.array(form, dtype=float))
    return [list(vectors[i, :]) for i in range(form.nrows()) if values[i] > 0]


def _positive_integral_vector(form):
    r"""GetIntegralVector_allmeth with CritNorm 0 and a strict inequality (Positivity.h:399-480).

    A positive CritNorm-free search is not an isotropic one, so the form is
    reduced by IndefiniteReduction and the approximants of the reduced form's
    diagonal and eigen directions are tried in turn (allmeth_V2, :312-397).
    """
    denominator = lcm(*(int(QQ(entry).denominator()) for entry in form.list()))
    integral = (form * denominator).change_ring(ZZ)
    basis, reduced = _indefinite_reduction(integral)
    iterators = [_ApproxIterator(direction) for direction in _positive_directions_diagonal(reduced)]
    iterators += [_ApproxIterator(direction) for direction in _positive_directions_eigen(reduced)]
    while True:
        for iterator in iterators:
            candidate = iterator.increment()
            if candidate.dot_product(reduced * candidate) > 0:
                return basis.transpose() * candidate


def _fifth_basis(form):
    r"""compute_fifth_basis (FifthAndHigherOrder.h:127-218): five vectors spanning a subspace with at most four squares of each sign."""
    form = form.change_ring(QQ)
    size = form.nrows()
    largest = 4
    basis = []

    def acceptable(candidate):
        if len(basis) == 5:
            return False
        if basis and matrix(QQ, [*basis, candidate]).rank() == len(basis):
            return False
        trial = matrix(QQ, [*basis, candidate])
        positive, negative = _signature_counts(trial * form * trial.transpose())
        return positive <= largest and negative <= largest

    for i in range(size):
        candidate = vector(QQ, [1 if k == i else 0 for k in range(size)])
        if acceptable(candidate):
            basis.append(candidate)
    for i in range(size):
        for j in range(i + 1, size):
            candidate = vector(QQ, [1 if k in (i, j) else 0 for k in range(size)])
            if acceptable(candidate):
                basis.append(candidate)

    def sign_vector(sign):
        spanned = matrix(QQ, basis)
        perpendicular = (spanned * form).right_kernel_matrix()
        restricted = sign * perpendicular * form * perpendicular.transpose()
        return perpendicular.transpose() * _positive_integral_vector(restricted)

    while len(basis) < 5:
        spanned = matrix(QQ, basis)
        positive, negative = _signature_counts(spanned * form * spanned.transpose())
        if positive < largest:
            basis.append(sign_vector(1))
        elif negative < largest:
            basis.append(sign_vector(-1))
    return matrix(QQ, basis)


def fifth_and_higher_order_isotropic_vector(form):
    r"""FifthAndHigherOrderIsotropicVector (FifthAndHigherOrder.h:250-315)."""
    form = form.change_ring(QQ)
    basis = identity_matrix(QQ, 5) if form.nrows() == 5 else _fifth_basis(form)
    restricted = basis * form * basis.transpose()
    if restricted.rank() < 5:
        return basis.transpose() * restricted.left_kernel().basis_matrix().row(0)
    # FifthOrderIsotropicVector (:220-248).
    transform, diagonal = _reduced_diagonal(restricted)
    solution = transform.transpose() * _solve_fifth_equation(diagonal).change_ring(QQ)
    return basis.transpose() * solution


# ---------------------------------------------------------------------------
# Dispatcher (Isotropic.h)
# ---------------------------------------------------------------------------


def find_isotropic_exact(form):
    r"""FindIsotropicExact (Isotropic.h:279-306): an isotropic rational vector of a nondegenerate form, or None."""
    form = form.change_ring(QQ)
    size = form.nrows()
    if size == 1:
        # FindIsotropicRankOne (:203-211).
        return vector(QQ, [1]) if form[0, 0] == 0 else None
    if size == 2:
        # FindIsotropicRankTwo (:213-249).
        a, b, c = form[0, 0], form[0, 1], form[1, 1]
        if a == 0:
            return vector(QQ, [1, 0])
        delta = b * b - a * c
        if not delta.is_square():
            return None
        return vector(QQ, [-b + delta.sqrt(), a])
    if size == 3:
        return ternary_isotropic_vector(form)
    if size == 4:
        return quaternary_isotropic_vector(form)
    return fifth_and_higher_order_isotropic_vector(form)


__all__ = ["find_isotropic_exact"]
