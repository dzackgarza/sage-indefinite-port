"""Exact certificates for the orbit-representative facet route, on maintained linear algebra.

The perfect-cell facet route enumerates facets modulo the cell stabilizer rather than asking an engine for the
complete H-description of the cell cone. Each step needs a decision that can be checked exactly, and every
decision here is expressed through the linear algebra the port and Sage already maintain -- ``matrix(...).rank()``,
``.pivots()``, ``.right_kernel()``, ``.det()`` -- not through a private elimination of its own.

Rank certificates.  A nonzero ``k x k`` minor proves ``rank >= k`` over `Q`, and it does so whether the
non-vanishing is established over `Z` or modulo a prime: ``det(M) != 0`` in `Z/pZ` implies ``det(M) != 0`` in
`Z`, hence in `Q`. (I wrote the opposite in an earlier revision; that was wrong. What modular arithmetic cannot
do is certify an *upper* bound -- ``det == 0`` mod `p` says nothing -- which is why the upper bounds below are
argued structurally, not computed.) The certificate produced here is the determinant itself, evaluated by Sage.

Free upper bounds, which is what makes the tests cheap:

* ``LorentzianPerfectCell.__post_init__`` (``lorentzian_cells.py:614-619``) raises unless ``rank[1 | A] = n``,
  which says the rows' *affine* span has dimension ``n-1``. It does **not** check ``rank A = n``, and the two
  are different: with ``rank A = n-1`` the affine hyperplane can pass through the origin (for example the rank-2
  rows ``(1,0), (2,0)`` have ``rank[1|A] = 2`` and ``rank A = 1``), and the cone then lies in a linear hyperplane,
  so it is **not full-dimensional**. Pointedness is a separate matter and does not fail in general -- that same
  example is a ray, which is pointed. The refusal is because the facet theory here needs a cone of dimension `n`.
  So the a-priori statement needs *both* ranks, and :func:`pointedness_functional` verifies both instead of
  assuming the constructor supplies them.
* A vector that is non-negative on every ray and vanishes on an incidence set `I` puts `I` inside an
  `(n-1)`-space, so only ``rank(I) >= n-1`` needs certifying.
* `I cap I'` is annihilated by two independent facet normals, so ``rank(I cap I') <= n-2`` is free and
  adjacency is exactly ``rank(I cap I') = n-2``.

Reference: the upstream route is `polyhedral_common`'s adjacency decomposition
(``src_lorentzian/lorentzian_perfect.h:1433-1472``). This module supplies the arithmetic those steps consume; it
does not walk the facet graph.
"""

from __future__ import annotations

from math import gcd

from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ

__all__ = [
    "adjacent_incidences",
    "facet_normal_from_incidence",
    "incidence_rank_certificate",
    "is_facet_certificate",
    "orient_facet_normal",
    "pointedness_functional",
]

type Normal = tuple[int, ...]


def _integer_rows(rows) -> list[tuple[int, ...]]:
    return [tuple(int(entry) for entry in row) for row in rows]


def incidence_rank_certificate(rows, wanted: int) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    r"""Certify ``rank(rows) >= wanted`` with a nonzero minor, or return None.

    Returns ``(row_indices, column_indices)`` selecting the square submatrix whose determinant Sage evaluates as
    nonzero, so the bound holds over `Q`. The rank itself is Sage's exact ``rank()``; the certificate is the
    determinant, so no private elimination is involved.
    """
    rows = _integer_rows(rows)
    wanted = int(wanted)
    if wanted <= 0:
        return (), ()
    if len(rows) < wanted or not rows:
        return None
    integer_matrix = matrix(ZZ, rows)
    if integer_matrix.rank() < wanted:
        return None
    pivot_columns = tuple(int(column) for column in integer_matrix.pivots())[:wanted]
    if len(pivot_columns) < wanted:
        return None
    restricted = integer_matrix.matrix_from_columns(pivot_columns)
    pivot_rows = tuple(int(row) for row in restricted.transpose().pivots())[:wanted]
    if len(pivot_rows) < wanted:
        return None
    square = integer_matrix.matrix_from_rows_and_columns(pivot_rows, pivot_columns)
    if square.det() == 0:
        return None
    return pivot_rows, pivot_columns


def facet_normal_from_incidence(incidence_rows, dimension: int) -> Normal:
    r"""The primitive integer normal of the hyperplane spanned by ``incidence_rows``.

    Sage's ``right_kernel`` over `Q` supplies the one-dimensional kernel of a rank-``n-1`` row set; the result is
    cleared to integers, made primitive, and then verified to annihilate *every* given row rather than only a
    spanning subset. Which sign is valid for the cone is decided by :func:`orient_facet_normal`, because the sign
    of a kernel basis is a fact about coordinates, not about the cone.
    """
    dimension = int(dimension)
    rows = _integer_rows(incidence_rows)
    if len(rows) < dimension - 1:
        raise ArithmeticError(
            f"an incidence set must have at least {dimension - 1} rows to span a hyperplane, got {len(rows)}"
        )
    kernel = matrix(QQ, rows).right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError(
            f"an incidence set must span a hyperplane, but its kernel has dimension {kernel.nrows()}"
        )
    rational = [QQ(entry) for entry in kernel.row(0)]
    denominator = 1
    for entry in rational:
        denominator = denominator * int(entry.denominator()) // gcd(denominator, int(entry.denominator()))
    integral = [int(entry * denominator) for entry in rational]
    divisor = 0
    for entry in integral:
        divisor = gcd(divisor, entry)
    if divisor == 0:
        raise ArithmeticError("the kernel vector vanished on a spanning incidence set")
    normal = tuple(entry // divisor for entry in integral)
    for row in rows:
        if sum(a * b for a, b in zip(normal, row, strict=True)) != 0:
            raise ArithmeticError("the kernel normal does not annihilate its incidence set")
    return normal


def orient_facet_normal(normal, ray_rows) -> tuple[Normal, list[int]]:
    r"""Orient a hyperplane normal to be **non-negative on every ray**, returning ``(normal, pairings)``.

    The valid orientation is the one with ``<normal, r> >= 0`` for every ray. A normal that annihilates every ray
    is not a supporting hyperplane of a full-dimensional cone: both orientations look non-negative, and the
    earlier version of this function returned that case as if it were a certificate. It raises now, and the test
    suite pins it.
    """
    normal = tuple(int(entry) for entry in normal)
    values = [sum(a * int(b) for a, b in zip(normal, row, strict=True)) for row in ray_rows]
    if all(value == 0 for value in values):
        raise ArithmeticError(
            "the normal annihilates every ray, so neither orientation is a supporting hyperplane of a "
            "full-dimensional cone"
        )
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
    non-negative on every other ray with at least one strict inequality, and that the incidence has rank
    ``n - 1`` (the upper bound is free, so only the certified lower bound enters).
    """
    normal = tuple(int(entry) for entry in normal)
    # Retain the ray sequence once. Pairings below convert every entry to
    # an exact integer, and the incidence rank check converts its own rows;
    # rebuilding every integer row first duplicates the full cone scan.
    ray_rows = tuple(ray_rows)
    values = [sum(a * int(b) for a, b in zip(normal, row, strict=True)) for row in ray_rows]
    return _is_facet_certificate_from_pairings(normal, incidence, ray_rows, values)


def _is_facet_certificate_from_pairings(normal, incidence, ray_rows, values) -> bool:
    r"""Certify using the exact pairings computed for these rows and this normal.

    Both callers obtain the values by evaluating the normal on every ray.
    The walk retains the result of ``orient_facet_normal``; the independent
    public certificate computes it directly. Rank, content and the complete
    zero set are checked here in both cases.
    """
    dimension = len(normal)
    zeros = frozenset(index for index, value in enumerate(values) if value == 0)
    if zeros != frozenset(int(index) for index in incidence):
        return False
    if any(value < 0 for value in values):
        return False
    if all(value == 0 for value in values):
        return False
    divisor = 0
    for entry in normal:
        divisor = gcd(divisor, entry)
    if divisor != 1:
        return False
    return incidence_rank_certificate([ray_rows[index] for index in sorted(zeros)], dimension - 1) is not None


def adjacent_incidences(first_incidence, second_incidence, ray_rows, dimension: int) -> bool:
    r"""Whether two facets meet in a ridge: ``rank(I cap I') = n - 2``, exact.

    The upper bound is free because the joint incidence is annihilated by two independent facet normals; the
    lower bound is Sage's exact rank on the joint rows.
    """
    joint = sorted(
        set(int(index) for index in first_incidence) & set(int(index) for index in second_incidence)
    )
    wanted = int(dimension) - 2
    if wanted <= 0:
        return True
    if len(joint) < wanted:
        return False
    ray_rows = _integer_rows(ray_rows)
    return matrix(ZZ, [ray_rows[index] for index in joint]).rank() >= wanted


def pointedness_functional(cell) -> tuple[Normal, list[int]]:
    r"""The strictly positive functional that the cell's hypotheses supply, computed exactly.

    Two hypotheses are needed and both are verified here, because the cell constructor supplies only the first:

    * ``rank[1 | A] = n`` -- the augmented rows span an `n`-dimensional subspace of `Q^{n+1}`, so their kernel is
      one-dimensional, spanned by ``w = (w_0, w')`` with ``w_0 + w'*a_i = 0`` for every row. This is what
      ``LorentzianPerfectCell.__post_init__`` already enforces.
    * ``rank A = n`` -- **not** enforced by the constructor, and not implied by the first. It is what puts the
      affine hyperplane off the origin: if ``w_0 = 0`` then ``w'`` annihilates ``span A = Q^n``, so ``w' = 0`` and
      ``w = 0``. With ``rank A = n-1`` the hyperplane can pass through the origin and the cone then lies in it, so
      it is **not full-dimensional**. (Pointedness is not what fails there: the ray generated by ``(1,0), (2,0)``
      is pointed. The refusal is that the facet theory needs a cone of dimension `n`.)

    Given both, ``w'*a_i = -w_0`` is the *same nonzero rational value on every ray*, and ``w'`` (sign-corrected,
    then cleared to a primitive integral vector) is strictly positive on every ray.

    Everything is exact rational arithmetic: the kernel is taken over `Q`, ``w_0`` is compared as a rational
    (never truncated to an integer, which would read ``1/2`` as ``0``), and the integer functional comes from the
    port's own ``_primitive_integral_direction`` -- the ``RemoveFractionPlusCoef`` contract, which returns a
    positive multiplier and a primitive integral vector satisfying ``x = v / m`` exactly.

    Returns ``(functional, pairings)``; raises with the failing hypothesis named.
    """
    rows = _integer_rows(cell.rows)
    dimension = len(rows[0])
    row_matrix = matrix(ZZ, rows)
    row_rank = int(row_matrix.rank())
    if row_rank != dimension:
        raise ArithmeticError(
            f"the rows have rank {row_rank} but the ambient dimension is {dimension}: with rank < dimension the "
            "affine hyperplane can pass through the origin, the cone is then not full-dimensional, and the facet "
            "theory here does not apply (the constructor does not check this hypothesis)"
        )
    augmented = matrix(ZZ, [[1, *row] for row in rows])
    augmented_rank = int(augmented.rank())
    if augmented_rank != dimension:
        raise ArithmeticError(
            f"rank[1|A] = {augmented_rank} but the lattice rank is {dimension}, so the rows do not lie on a "
            "single affine hyperplane"
        )
    kernel = matrix(QQ, augmented).right_kernel().basis_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError(
            f"the augmented kernel has dimension {kernel.nrows()}, expected 1 for rank[1|A] = n"
        )
    rational_kernel = [QQ(entry) for entry in kernel.row(0)]
    if rational_kernel[0] == 0:
        raise ArithmeticError(
            "the affine hyperplane passes through the origin (w_0 = 0), so the rows span a linear hyperplane and "
            "the cone is not full-dimensional"
        )
    rational_functional = tuple(rational_kernel[1:])
    if not any(entry != 0 for entry in rational_functional):
        raise ArithmeticError("the affine functional is identically zero, so the rows do not span a hyperplane")

    from sage_indefinite_port.indefinite.lorentzian_cells import _primitive_integral_direction

    multiplier, primitive = _primitive_integral_direction(rational_functional, dimension)
    functional = tuple(int(entry) for entry in primitive)
    pairings = [sum(a * int(b) for a, b in zip(functional, row, strict=True)) for row in rows]
    if any(value == 0 for value in pairings):
        raise ArithmeticError(
            "the affine functional vanishes on a ray, which cannot happen for a cone over a polytope"
        )
    if len(set(pairings)) != 1:
        raise ArithmeticError(
            f"the affine functional is not constant on the rays (values {sorted(set(pairings))[:4]}), so the rows "
            "do not lie in one affine hyperplane"
        )
    if pairings[0] < 0:
        functional = tuple(-entry for entry in functional)
        pairings = [-value for value in pairings]
    if multiplier <= 0:
        raise ArithmeticError("the denominator-clearing multiplier must be positive")
    return functional, pairings
