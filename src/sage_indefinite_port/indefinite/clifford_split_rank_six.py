r"""The degree-four model of the split rank-six lattice ``3U``.

The Fock representation of ``C(3U)`` on the exterior algebra of ``QQ^3`` sends ``e_i`` to
wedging and ``f_i`` to twice contraction, so ``C_0`` maps isomorphically onto the two
parity blocks ``M_4 x M_4``. For the second model, ``x ^ y = b(x, y) e_1234`` on
``Lambda^2 QQ^4`` is ``3U`` in the basis ``(e12, e34), (e13, -e24), (e14, e23)``, and
``SL_4`` acts by the second compound matrix, with kernel ``{1, -1}`` on rational
points (Knus, Chapter 11; issue #57).
"""

from itertools import combinations

from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage_extensions.all import ZZ as OwnedZZ
from sage_extensions.all import Lattices

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


def isometry_of_special_linear(g, lattice):
    r"""Return the element of ``O(lattice)`` given by integral ``g`` in ``SL_4``."""
    action = second_compound(g)
    if action.denominator() != 1:
        raise ValueError(f"{g} does not act integrally on Lambda^2")
    points = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    images = (sum((int(entry) * point for entry, point in zip(column, points, strict=True)), lattice.zero()) for column in action.columns())
    return lattice.O()(dict(zip(labels, images, strict=True)))
