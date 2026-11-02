"""Exact certificates for the orbit-representative facet route.

The perfect-cell facet route enumerates facets modulo the cell stabilizer instead of asking an engine for the
complete H-description of the cell cone. It needs three exact primitives, each with a certificate that can be
checked by integer arithmetic alone:

* :func:`rank_lower_certificate` -- certify `rank(rows) >= wanted` by an exact nonzero minor. Modular
  elimination only *finds* the candidate rows (a nonzero reduction never proves a lower bound); the exact
  Bareiss determinant certifies it.
* :func:`facet_normal_from_incidence` -- recover a facet's primitive integer normal from its incidence set
  alone. Cofactor construction over `dimension - 1` independent rows, so the cost does not depend on how many
  rays the facet contains.
* :func:`is_facet_certificate` -- decide "this primitive vector is the facet with this incidence".

Two facts make the tests cheap rather than expensive:

* :meth:`LorentzianPerfectCell.__post_init__` (``lorentzian_cells.py:614-619``) raises unless
  ``rank[1 | A] = n``. So the cell's rows lie in a single affine hyperplane not through the origin, and the
  cone is the cone over an ``(n-1)``-dimensional polytope: **pointed and full-dimensional a priori**, not by
  inference from ``rank A = n`` (which gives full-dimensionality only).
* Consequently the upper bound in both rank tests is free. A vector that is non-negative on every ray and
  vanishes on the incidence set puts that set inside an ``(n-1)``-space, so only ``rank >= n-1`` needs
  certifying; and ``I cap I'`` is annihilated by two independent facet normals, so ``rank <= n-2`` is free and
  adjacency is exactly ``rank(I cap I') = n-2``.

Reference: the upstream route is `polyhedral_common`'s adjacency decomposition
(``src_lorentzian/lorentzian_perfect.h:1433-1472``, ``f_adj`` -> ``DualDescriptionRecordFullDim(EXT_T,
GRPperm, rddo)``). This module supplies the arithmetic those steps consume; it does not walk the facet graph.
"""

from __future__ import annotations

from math import gcd

__all__ = [
    "adjacent_incidences",
    "facet_normal_from_incidence",
    "is_facet_certificate",
    "orient_facet_normal",
    "rank_lower_certificate",
]


def _hadamard_bound(dimension: int, max_entry: int) -> int:
    r"""An entry of a `dimension x dimension` integer minor is at most this in absolute value."""
    return int(max_entry) ** int(dimension) * _factorial_root(int(dimension))


def _factorial_root(dimension: int) -> int:
    """A cheap integer upper bound for `sqrt(dimension!)`, used only to size a prime."""
    value = 1
    for factor in range(2, int(dimension) + 1):
        value *= factor
    root = 1
    while root * root <= value:
        root += 1
    return root


def _prime_above(limit: int) -> int:
    """The least prime strictly above ``limit``."""
    candidate = max(3, int(limit) + 1)
    if candidate % 2 == 0:
        candidate += 1
    while True:
        if all(candidate % divisor for divisor in range(3, int(candidate**0.5) + 1, 2)):
            return candidate
        candidate += 2


def _determinant(matrix_rows) -> int:
    """Exact integer determinant by fraction-free Bareiss elimination."""
    size = len(matrix_rows)
    if size == 0:
        return 1
    working = [[int(entry) for entry in row] for row in matrix_rows]
    if any(len(row) != size for row in working):
        raise ValueError("the determinant requires a square matrix")
    sign, previous = 1, 1
    for column in range(size - 1):
        if working[column][column] == 0:
            swap = next((row for row in range(column + 1, size) if working[row][column]), None)
            if swap is None:
                return 0
            working[column], working[swap] = working[swap], working[column]
            sign = -sign
        pivot = working[column][column]
        for row in range(column + 1, size):
            for position in range(column + 1, size):
                working[row][position] = (
                    working[row][position] * pivot - working[row][column] * working[column][position]
                ) // previous
            working[row][column] = 0
        previous = pivot
    return sign * working[size - 1][size - 1]


def exact_minor_is_nonzero(matrix_rows) -> bool:
    """Whether a square integer minor is exactly nonzero."""
    return _determinant(matrix_rows) != 0


def rank_lower_certificate(rows, wanted: int):
    r"""Certify ``rank(rows) >= wanted``, returning ``(indices, ok)``.

    ``indices`` selects exactly ``wanted`` rows whose ``wanted x wanted`` minor is exactly nonzero, so the
    bound holds unconditionally over ``Q``; ``ok`` is False when the rows cannot certify it. Modular
    elimination over a prime above the Hadamard bound only finds the candidate pivot rows -- a nonzero
    reduction never proves a lower bound -- and the exact Bareiss determinant is what certifies.
    """
    rows = [tuple(int(entry) for entry in row) for row in rows]
    wanted = int(wanted)
    if wanted <= 0:
        return (), True
    if len(rows) < wanted or not rows:
        return (), False
    dimension = len(rows[0])
    max_entry = max((abs(entry) for row in rows for entry in row), default=1)
    prime = _prime_above(_hadamard_bound(max(wanted, 2), max(max_entry, 2)))
    basis: list[tuple[int, list[int]]] = []
    chosen: list[int] = []
    pivot_columns: list[int] = []
    for index, row in enumerate(rows):
        working = [value % prime for value in row]
        for pivot_column, pivot_row in basis:
            factor = working[pivot_column]
            if factor:
                working = [(a - factor * b) % prime for a, b in zip(working, pivot_row, strict=True)]
        pivot_column = next((column for column, value in enumerate(working) if value), None)
        if pivot_column is None:
            continue
        inverse = pow(working[pivot_column], -1, prime)
        basis.append((pivot_column, [(value * inverse) % prime for value in working]))
        chosen.append(index)
        pivot_columns.append(pivot_column)
        if len(basis) == wanted:
            break
    if len(chosen) < wanted:
        return (), False
    square = [[rows[index][column] for column in pivot_columns] for index in chosen]
    if len(square) != wanted or len(pivot_columns) != wanted:
        return (), False
    if not exact_minor_is_nonzero(square):
        return (), False
    return tuple(chosen), True


def independent_rows(rows, wanted: int):
    """A certified independent subfamily of size ``wanted``, or None."""
    indices, ok = rank_lower_certificate(rows, wanted)
    if not ok:
        return None
    return tuple(tuple(int(entry) for entry in rows[index]) for index in indices)


def facet_normal_from_incidence(incidence_rows, dimension: int) -> tuple[int, ...]:
    r"""The primitive integer normal of the hyperplane spanned by ``incidence_rows``.

    Cofactor construction ``y_i = (-1)^i det(M_i)`` over ``dimension - 1`` certified independent rows, so the
    cost is independent of how many rows the incidence has. The result annihilates every given row (checked,
    not assumed) and is primitive with its first nonzero entry positive; orientation is the caller's business.
    """
    dimension = int(dimension)
    basis = independent_rows(tuple(incidence_rows), dimension - 1)
    if basis is None:
        raise ArithmeticError("an incidence set must span a hyperplane of the ray space")
    cofactors = []
    for column in range(dimension):
        minor = [[row[index] for index in range(dimension) if index != column] for row in basis]
        cofactors.append(((-1) ** column) * _determinant(minor))
    divisor = 0
    for value in cofactors:
        divisor = gcd(divisor, value)
    if divisor == 0:
        raise ArithmeticError("the cofactor vector vanished on a spanning incidence set")
    normal = tuple(value // divisor for value in cofactors)
    first = next((value for value in normal if value), 0)
    if first < 0:
        normal = tuple(-value for value in normal)
    for row in incidence_rows:
        if sum(a * b for a, b in zip(normal, row, strict=True)) != 0:
            raise ArithmeticError("the cofactor normal does not annihilate its incidence set")
    return normal


def orient_facet_normal(normal, ray_rows) -> tuple[tuple[int, ...], list[int]]:
    r"""Orient a hyperplane normal so that it is **non-negative on every ray**.

    The cofactor construction fixes only the sign of the first nonzero entry, which is a fact about coordinates
    and not about the cone: the valid orientation is the one with ``<normal, r> >= 0`` for every ray. Exactly one
    of the two signs is non-negative on all rays when the hyperplane is a facet's; both are when the normal
    annihilates every ray, which cannot happen for a full-dimensional cone and is refused here.

    Returns ``(oriented_normal, pairings)`` so the caller can certify without recomputing them.
    """
    normal = tuple(int(entry) for entry in normal)
    values = [sum(a * int(b) for a, b in zip(normal, row, strict=True)) for row in ray_rows]
    if all(value >= 0 for value in values):
        return normal, values
    flipped = tuple(-entry for entry in normal)
    flipped_values = [-value for value in values]
    if all(value >= 0 for value in flipped_values):
        return flipped, flipped_values
    raise ArithmeticError(
        "neither orientation of the normal is non-negative on every ray, so it is not a supporting hyperplane"
    )


def is_facet_certificate(normal, incidence, ray_rows) -> bool:
    r"""Whether ``normal`` is the facet whose incidence is ``incidence``, by exact arithmetic.

    Checks primitivity and nontriviality, that the pairings vanish exactly on ``incidence``, that they are
    strictly positive on every other ray, and that the incidence has rank ``n - 1``. The rank's upper bound is
    free (the normal annihilates the set), so only the certified lower bound enters.
    """
    normal = tuple(int(entry) for entry in normal)
    dimension = len(normal)
    values = [sum(a * int(b) for a, b in zip(normal, row, strict=True)) for row in ray_rows]
    zeros = frozenset(index for index, value in enumerate(values) if value == 0)
    if zeros != frozenset(int(index) for index in incidence):
        return False
    if any(value < 0 for value in values):
        return False
    divisor = 0
    for entry in normal:
        divisor = gcd(divisor, entry)
    if divisor != 1:
        return False
    _, ok = rank_lower_certificate([ray_rows[index] for index in sorted(zeros)], dimension - 1)
    return ok


def adjacent_incidences(first_incidence, second_incidence, ray_rows, dimension: int) -> bool:
    r"""Whether two facets sharing no common interior meet in a ridge: ``rank(I cap I') = n - 2``.

    The upper bound is free: the joint incidence is annihilated by two independent facet normals, so it lies
    in an ``(n-2)``-space. Adjacency therefore needs only the certified lower bound, and a partial rank test
    can exit as soon as ``n - 2`` independent joint rows are found.
    """
    joint = sorted(set(int(index) for index in first_incidence) & set(int(index) for index in second_incidence))
    wanted = int(dimension) - 2
    if wanted <= 0:
        return True
    if len(joint) < wanted:
        return False
    _, ok = rank_lower_certificate([ray_rows[index] for index in joint], wanted)
    return ok
