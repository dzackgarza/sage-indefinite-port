"""Exact primitives for orbit-level facet enumeration of a perfect-cell cone.

The enumeration (the reference's orbit-wise adjacency decomposition) needs three
geometric decisions per facet orbit, and nothing else:

  (P1) is a candidate covector a *facet* normal -- primitive, nonnegative on every
       ray of the cell, tight exactly on a rank-(n-1) incidence set?
  (P2) are two facet normals *adjacent* -- do the two facets share a ridge?
  (P3) what is the exact facet normal belonging to a *given* incidence set?

This module answers (P1)-(P3) in integer arithmetic, with a certificate attached
to every answer.  Two facts make that cheap; both are proved in solution.md:

  * Upper bounds on incidence rank are free.  A facet normal y annihilates its
    incidence set I, so rank(I) <= n-1 with no computation; two distinct normals
    y, y' annihilate I cap I', so rank(I cap I') <= n-2.  Only the *lower* bound
    ever needs a certificate, and a lower bound is certified by exhibiting k
    independent rows together with one nonzero k x k minor.

  * Adjacency (P2) is decided by the reference's joint-rank criterion,
    rank{a_i : <y,a_i> = <y',a_i> = 0} = n-2, and that rank has a *free upper
    bound*: both normals annihilate the joint set and are independent, so the
    rank is at most n-2 always.  Deciding adjacency therefore only needs a
    *lower* bound, i.e. n-2 independent rows out of the joint incidence, which an
    early-exit modular search usually finds in the first few rows.

    A tempting shortcut is *vacuous* and is recorded here so nobody repeats it:
    "some ray is strictly negative on both normals" cannot hold when both normals
    are facet normals, because facet normals are nonnegative on every ray.  The
    condition is necessary for adjacency but is never satisfiable, so it filters
    nothing; the exact joint-rank test (P2) is mandatory, not optional.

Nothing here uses floating point, a tolerance, or a probabilistic rank test: the
one modular reduction is used only to *find* candidate pivots, and every pivot
choice is re-certified by an exact Bareiss minor, while the upper bound comes from
the annihilator, not from a prime.
"""

from __future__ import annotations

from math import gcd, isqrt

import numpy as np

# Largest value a sum of n products may reach and still be exact in int64.
_INT64_SAFE = 2**63 - 1


def hadamard_bound(dimension: int, max_entry: int) -> int:
    """|det| of any `dimension` x `dimension` integer matrix with entries <= max_entry.

    Hadamard: |det| <= prod of row euclidean norms <= (sqrt(d) * max_entry)**d.
    Used to size certificate data and to bound the cofactor vector of (P3).
    """
    return (isqrt(dimension) + 1) ** dimension * max_entry**dimension


def _is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin for n < 3.3e24."""
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_above(limit: int) -> int:
    """The least prime strictly greater than `limit`."""
    candidate = max(3, int(limit) + 1) | 1
    while not _is_prime(candidate):
        candidate += 2
    return candidate


def rank_lower_certificate(rows, wanted: int, *, full_rank: bool = False):
    """Certify rank(rows) >= `wanted`, or report that this row set cannot.

    Returns (indices, rank_mod_p).  `indices` selects `wanted` rows of `rows`
    whose `wanted` x `wanted` minor is *exactly* nonzero (`exact_minor_is_nonzero`),
    so rank >= wanted holds unconditionally over Q.  `rank_mod_p` is the rank of
    the whole row set modulo a prime above the Hadamard bound of the row space
    when `full_rank` is requested (a hint for the caller, never the certificate);
    otherwise it is None and the search stops at `wanted` pivots.

    Cost O(len(rows) * wanted * n) small-integer operations, or O(wanted^2 * n)
    when the early exit succeeds, with no rational arithmetic anywhere.
    """
    rows = [tuple(int(entry) for entry in row) for row in rows]
    if not rows:
        return (), (0 if full_rank else None)
    n = len(rows[0])
    max_entry = max((abs(entry) for row in rows for entry in row), default=1)
    prime = prime_above(hadamard_bound(max(n, 2), max(max_entry, 2)))
    basis: list[tuple[int, list[int]]] = []
    chosen: list[int] = []
    pivot_columns: list[int] = []
    for index, row in enumerate(rows):
        working = [value % prime for value in row]
        for pivot_column, pivot_row in basis:
            factor = working[pivot_column]
            if factor:
                working = [(a - factor * b) % prime for a, b in zip(working, pivot_row, strict=True)]
        try:
            pivot_column = next(column for column, value in enumerate(working) if value)
        except StopIteration:
            continue
        inverse = pow(working[pivot_column], -1, prime)
        basis.append((pivot_column, [(value * inverse) % prime for value in working]))
        chosen.append(index)
        pivot_columns.append(pivot_column)
        if len(basis) == wanted:
            break
    rank_mod_p = _rank_mod_p(rows, prime, n) if full_rank else None
    if len(chosen) < wanted:
        return (), rank_mod_p
    # The chosen rows form an upper-triangular matrix on their own pivot columns
    # modulo the prime, so that square minor cannot vanish mod p; certify it over
    # Z by an exact Bareiss determinant rather than trusting the reduction.
    square = [[rows[index][column] for column in pivot_columns] for index in chosen]
    if not exact_minor_is_nonzero(square):
        # A modular pivot that fails over Z is not a certificate; refuse rather
        # than return an uncertified lower bound.
        return (), rank_mod_p
    return tuple(chosen), rank_mod_p


def _rank_mod_p(rows, prime: int, n: int) -> int:
    basis: list[tuple[int, list[int]]] = []
    for row in rows:
        working = [value % prime for value in row]
        for pivot_column, pivot_row in basis:
            factor = working[pivot_column]
            if factor:
                working = [(a - factor * b) % prime for a, b in zip(working, pivot_row, strict=True)]
        try:
            pivot_column = next(column for column, value in enumerate(working) if value)
        except StopIteration:
            continue
        inverse = pow(working[pivot_column], -1, prime)
        basis.append((pivot_column, [(value * inverse) % prime for value in working]))
    return len(basis)


def exact_minor_is_nonzero(rows) -> bool:
    """Exact decision det(rows) != 0 by fraction-free (Bareiss) elimination.

    `rows` must be a *square* integer matrix: a rank certificate is a minor, and
    a non-square input is a programming error, not a rank statement.
    """
    matrix = [[int(entry) for entry in row] for row in rows]
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("a minor certificate must be square")
    if size == 0:
        return True
    previous = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if matrix[row][column]), None)
        if pivot is None:
            return False
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        pivot_value = matrix[column][column]
        for row in range(column + 1, size):
            for inner in range(column + 1, size):
                matrix[row][inner] = (
                    matrix[row][inner] * pivot_value
                    - matrix[row][column] * matrix[column][inner]
                ) // previous
            matrix[row][column] = 0
        previous = pivot_value
    return matrix[size - 1][size - 1] != 0


def _determinant(matrix) -> int:
    """Exact determinant of a small integer matrix (Bareiss)."""
    values = [[int(entry) for entry in row] for row in matrix]
    size = len(values)
    if size == 0:
        return 1
    sign = 1
    previous = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if values[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            sign = -sign
        pivot_value = values[column][column]
        for row in range(column + 1, size):
            for inner in range(column + 1, size):
                values[row][inner] = (
                    values[row][inner] * pivot_value
                    - values[row][column] * values[column][inner]
                ) // previous
            values[row][column] = 0
        previous = pivot_value
    return sign * values[size - 1][size - 1]


def independent_rows(rows, wanted: int):
    """`wanted` rows independent over Q, or None if there are not that many."""
    indices, _rank_mod_p = rank_lower_certificate(rows, wanted)
    if len(indices) < wanted:
        return None
    return [tuple(int(entry) for entry in rows[index]) for index in indices]


def primitive_normal_from_incidence(incidence_rows, dimension: int) -> tuple[int, ...]:
    """(P3) The primitive integer normal of the hyperplane spanned by `rows`.

    `incidence_rows` has to span a hyperplane of R^dimension (rank dimension-1).
    The result y annihilates every given row, is primitive (entry gcd 1) and is
    sign-normalized by its first nonzero entry.  The facet *orientation* is fixed
    separately by `orient_normal`, as `flip_across` does.

    Cofactor construction: y_i = (-1)^i det(M_i), M_i deleting column i, using
    dimension-1 independent rows only -- O(dimension^4) integer operations,
    independent of |I|.  The annihilation is then certified on *every* given row,
    not just the chosen basis.
    """
    basis = independent_rows(incidence_rows, dimension - 1)
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


def orient_normal(normal, outside_ray) -> tuple[int, ...]:
    """Flip `normal` so that it is strictly positive on `outside_ray`."""
    value = sum(a * b for a, b in zip(normal, outside_ray, strict=True))
    if value == 0:
        raise ArithmeticError("an outside ray must not lie on the facet")
    return tuple(-entry for entry in normal) if value < 0 else tuple(normal)


def pairings(normal, ray_rows) -> np.ndarray:
    """Exact integer pairings of one covector against every ray (int64 path).

    int64 is used only after certifying that no sum of products can overflow:
    the bound is |normal|_1 * max|ray entry| < 2**63.  The check is a proof, not
    a guard against silent wraparound.
    """
    rays = np.asarray(ray_rows, dtype=np.int64)
    vector = np.asarray(normal, dtype=np.int64)
    bound = int(np.abs(vector).sum()) * int(np.abs(rays).max(initial=0))
    if bound > _INT64_SAFE:
        raise ArithmeticError("pairings would exceed int64; use pairings_exact")
    return rays @ vector


def pairings_exact(normal, ray_rows) -> list[int]:
    """The same pairings with unbounded Python integers."""
    return [sum(a * b for a, b in zip(normal, row, strict=True)) for row in ray_rows]


def is_facet_certificate(normal, incidence, ray_rows) -> bool:
    """(P1) Certificate that `normal` is the facet with this incidence set.

    Checks exactly: primitivity and nontriviality; vanishing exactly on
    `incidence`; strict positivity on every other ray; and rank(incidence) = n-1,
    where the upper bound is free (the normal annihilates the incidence set) and
    the lower bound is a certified independent subset with a nonzero minor.
    """
    dimension = len(normal)
    values = pairings_exact(normal, ray_rows)
    zeros = frozenset(index for index, value in enumerate(values) if value == 0)
    if zeros != frozenset(incidence):
        return False
    if any(value < 0 for value in values):
        return False
    divisor = 0
    for entry in normal:
        divisor = gcd(divisor, entry)
    if divisor != 1:
        return False
    indices, _rank_mod_p = rank_lower_certificate(
        [ray_rows[index] for index in sorted(zeros)], dimension - 1
    )
    return len(indices) == dimension - 1


def zero_support(normal, ray_rows) -> int:
    """The incidence set of `normal` as a bitset: bit i iff <normal, a_i> = 0."""
    values = pairings(normal, ray_rows)
    packed = np.packbits(values == 0, bitorder="little")
    return int.from_bytes(packed.tobytes(), "little")


def adjacent_from_zero_sets(first_zeros: int, second_zeros: int, ray_rows, dimension: int):
    """(P2) Exact adjacency from the two incidence bitsets, with an early exit.

    Decides rank(J) = n-2 for J = I cap I'.  The upper bound n-2 is free (both
    normals annihilate J and are independent), so only the lower bound is
    certified: reduce rows of J modulo a prime until n-2 pivots appear, then
    re-certify that minor exactly with Bareiss.

    Returns (adjacent, rows_examined).  For adjacent pairs the loop usually stops
    after n-2 rows, so the cost is a few modular reductions whatever |J| is; for
    non-adjacent pairs it has to consume all of J, which is the honest worst case.
    """
    joint = first_zeros & second_zeros
    rows = []
    positions = []
    while joint and len(rows) < len(ray_rows):
        lowest = joint & -joint
        position = lowest.bit_length() - 1
        joint ^= lowest
        rows.append(tuple(ray_rows[position]))
        positions.append(position)
        if len(rows) >= dimension - 2:
            indices, _rank_mod_p = rank_lower_certificate(rows, dimension - 2)
            if len(indices) == dimension - 2:
                return True, len(rows)
            # keep going: an independent subset may need more of J
    indices, _rank_mod_p = rank_lower_certificate(rows, dimension - 2)
    return len(indices) == dimension - 2, len(rows)


def adjacent_rank_test(normal, other_normal, ray_rows) -> bool:
    """(P2) The reference criterion, straight from the definitions."""
    return adjacent_from_zero_sets(
        zero_support(normal, ray_rows), zero_support(other_normal, ray_rows), ray_rows, len(normal)
    )[0]


def negative_support(normal, ray_rows) -> int:
    """(P2) state of one facet orbit: bit i set iff <normal, a_i> < 0.

    This is the *only* per-facet data the closure test needs for its cheap
    sanity query, because a ray negative on both normals would rule adjacency out
    (Lemma 5 of solution.md: necessary, and vacuous in practice, so it decides
    nothing). (P2) works on the *incident* bitset; this one is the negative side.
    """
    values = pairings(normal, ray_rows)
    packed = np.packbits(values < 0, bitorder="little")
    return int.from_bytes(packed.tobytes(), "little")


def supports_are_disjoint(first: int, second: int) -> bool:
    """(P2) in bitset form: two facet orbits are adjacent iff these are disjoint."""
    return not (first & second)
