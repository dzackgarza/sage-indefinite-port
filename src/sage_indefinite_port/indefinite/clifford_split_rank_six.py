r"""The degree-four model of the split rank-six lattice ``3U``.

The Fock representation of ``C(3U)`` on the exterior algebra of ``QQ^3`` sends ``e_i`` to
wedging and ``f_i`` to twice contraction, so ``C_0`` maps isomorphically onto the two
parity blocks ``M_4 x M_4``. For the second model, ``x ^ y = b(x, y) e_1234`` on
``Lambda^2 QQ^4`` is ``3U`` in the basis ``(e12, e34), (e13, -e24), (e14, e23)``, and
``SL_4`` acts by the second compound matrix, with kernel ``{1, -1}`` on rational
points (Knus, Chapter 11; issue #57).

Transfer to ``O(3U)``. Write ``b`` for the bilinear form and ``Q = b / 2``; the spinor
norm ``theta`` sends a reflection ``s_v`` to the class of ``Q(v)``. If ``Lambda^2 g = 1``
then ``g`` fixes every plane, hence every line, so it is a scalar ``c`` with ``c^2 = 1``:
the kernel on ``SL_4(ZZ)`` is ``{1, -1}``. ``SL_4 -> SO(3U)`` is the spin covering, so its
image on rational points is ``ker theta`` on ``SO(V)``. If ``Lambda^2 g`` preserves
``Lambda^2 ZZ_p^4``, the Cartan form ``diag(p^{a_i})`` of ``g`` has every ``a_i + a_j = 0`` for
``i < j``, so every ``a_i`` is ``0``; thus ``g`` is integral at every ``p``. Therefore the image of
``SL_4(ZZ)`` is exactly ``SO(3U) & ker theta``. ``O(3U)`` is even unimodular, so ``theta`` takes
values in ``{1, -1}`` (the Kummer sequence of the smooth group scheme ``SO(Q)`` over ``ZZ``, where
``H^1(ZZ, mu_2) = {1, -1}``), and ``(det, theta)`` has at most four values on ``O(3U)``.
"""

from itertools import combinations

from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage_extensions.all import ZZ as OwnedZZ
from sage_extensions.all import Lattices
from sage_extensions.all import QQ as OwnedQQ

_PAIRS = ((0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2))
_SIGNS = (1, 1, 1, -1, 1, 1)


def fock_generators():
    r"""Return the matrices of ``e_1, f_1, e_2, f_2, e_3, f_3`` on ``Lambda QQ^3``."""
    subsets = [frozenset(c) for r in range(4) for c in combinations(range(3), r)]
    position = {subset: index for index, subset in enumerate(subsets)}
    result = []
    for i in range(3):
        wedge = matrix(QQ, 8, 8)
        contract = matrix(QQ, 8, 8)
        for subset in subsets:
            sign = (-1) ** sum(1 for element in subset if element < i)
            if i in subset:
                contract[position[subset - {i}], position[subset]] = sign
            else:
                wedge[position[subset | {i}], position[subset]] = sign
        result.extend((wedge, 2 * contract))
    return tuple(result)


def even_part_dimension():
    r"""Return the dimension of the span of the even monomials in the Fock generators."""
    gammas = fock_generators()
    monomials = []
    for size in (0, 2, 4, 6):
        for subset in combinations(range(6), size):
            product = matrix.identity(QQ, 8)
            for index in subset:
                product = product * gammas[index]
            monomials.append(product.list())
    return matrix(QQ, monomials).rank()


def second_compound(g):
    r"""Return the action of ``g`` in ``SL_4`` on the basis of ``3U`` in ``Lambda^2``."""
    if g.nrows() != 4 or g.det() != 1:
        raise ValueError("second_compound requires an element of SL_4")
    return matrix(QQ, 6, 6, [
        _SIGNS[row] * _SIGNS[column] * g[list(_PAIRS[row]), list(_PAIRS[column])].det()
        for row in range(6) for column in range(6)
    ])


def three_hyperbolic_planes():
    hyperbolic = Lattices(OwnedZZ)("U")
    return hyperbolic + hyperbolic + hyperbolic


def isometry_from_matrix(action, lattice):
    r"""Return the element of ``O(lattice)`` whose columns are the images of the basis."""
    if action.denominator() != 1:
        raise ValueError(f"{action} does not act integrally on the lattice")
    points = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    images = (sum((int(entry) * point for entry, point in zip(column, points, strict=True)), lattice.zero()) for column in action.columns())
    return lattice.O()(dict(zip(labels, images, strict=True)))


def isometry_of_special_linear(g, lattice):
    r"""Return the element of ``O(lattice)`` given by integral ``g`` in ``SL_4``."""
    return isometry_from_matrix(second_compound(g), lattice)


def elementary_generators():
    r"""Return the elementary matrices ``E_ij(1)``, ``i != j``, which generate ``SL_4(ZZ)``."""
    return tuple(
        matrix.identity(QQ, 4) + matrix(QQ, 4, 4, {(i, j): 1})
        for i in range(4) for j in range(4) if i != j
    )


def spinor_kernel_generators(lattice):
    r"""Return generators of ``SO(3U) & ker theta``, the image of ``SL_4(ZZ)`` modulo ``{1, -1}``."""
    return tuple(isometry_of_special_linear(g, lattice) for g in elementary_generators())


def reflection_matrix(vector_entries, gram):
    r"""Return the matrix of ``x -> x - 2 b(x, v) / b(v, v) v`` in the basis of ``gram``."""
    v = matrix(QQ, len(vector_entries), 1, list(vector_entries))
    return matrix.identity(QQ, gram.nrows()) - 2 * v * v.transpose() * gram / (v.transpose() * gram * v)[0, 0]


def _gram_of(lattice):
    points = tuple(lattice.module_generators())
    return matrix(QQ, [[int(lattice.b(x, y)) for y in points] for x in points])


def hyperbolic_coset_representatives(lattice):
    r"""Return representatives of the nontrivial ``(det, theta)`` classes, keyed by value.

    The first two basis vectors span a hyperbolic plane ``(e, f)``, with ``Q(e + f) = 1`` and
    ``Q(e - f) = -1``. So ``s_{e+f}`` has ``(det, theta) = (-1, 1)``, ``s_{e-f}`` has ``(-1, -1)``, and
    their product, ``-1`` on ``(e, f)`` and ``1`` on the complement, has ``(1, -1)``.
    """
    gram = _gram_of(lattice)
    rank = gram.nrows()
    plus = reflection_matrix([1, 1] + [0] * (rank - 2), gram)
    minus = reflection_matrix([1, -1] + [0] * (rank - 2), gram)
    return {
        (-1, 1): isometry_from_matrix(plus, lattice),
        (-1, -1): isometry_from_matrix(minus, lattice),
        (1, -1): isometry_from_matrix(minus * plus, lattice),
    }


def determinant_and_spinor_values(isometry, lattice):
    r"""Return ``(det, theta)`` of an element of ``O(lattice)``, for ``theta(s_v) = Q(v)``."""
    norm = lattice.spinor_norm(form_multiplier=OwnedQQ(1) / OwnedQQ(2))(isometry).representative()
    if norm.is_square():
        sign = 1
    elif (-norm).is_square():
        sign = -1
    else:
        raise ValueError(f"the spinor norm {norm} of {isometry} is neither 1 nor -1 modulo squares")
    return (int(isometry.determinant()), sign)


def orthogonal_group_generators(lattice):
    r"""Return generators of the full ``O(3U)``: the spinor kernel and its three cosets."""
    return spinor_kernel_generators(lattice) + tuple(hyperbolic_coset_representatives(lattice).values())
