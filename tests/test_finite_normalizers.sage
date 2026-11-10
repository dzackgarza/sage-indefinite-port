r"""Normalizers of noncyclic finite subgroups of O(L) for indefinite L (#50).

O(U) = {1, -1, s, -s} with s the swap of e and f.

L = U + U, H = <-1 + 1, 1 + -1>, a Klein four group. An element commuting with H preserves
both eigenspaces, so C(H) = O(U) x O(U) of order 16. An automorphism of H is realized only
if it permutes the eigenspaces by isometries; the swap of the summands does, and an
automorphism sending -1 + 1 to the central -1 kills an eigenspace. So N(H) has order 32.

L = U + U(2): the summands are not isometric, so N(H) = C(H), of order 16.

L = U, H = <s, -1>. Its eigenlattices are Z(e + f) and Z(e - f), of squares 2 and -2, glued
along Z/2. C(H) = O(U), and no automorphism other than the identity is realized.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.finite_normalizers import FiniteSubgroupNormalizer


def _order(lattice, isometries):
    frame = tuple(lattice.module_generators())
    as_matrix = lambda g: matrix(SageQQ, [[int(c) for c in g(x)] for x in frame])
    generators = [as_matrix(g) for g in isometries]
    identity = matrix.identity(SageQQ, len(frame))
    seen = set([tuple(identity.list())])
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            following = current * generator
            if tuple(following.list()) not in seen:
                seen.add(tuple(following.list()))
                frontier.append(following)
    return len(seen)


def _signs(lattice):
    basis = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    half = len(basis) // int(2)
    left = lattice.O()(dict(zip(labels, tuple(-b if i < half else b for i, b in enumerate(basis)), strict=True)))
    right = lattice.O()(dict(zip(labels, tuple(b if i < half else -b for i, b in enumerate(basis)), strict=True)))
    return left, right


def test_klein_four_on_two_hyperbolic_planes():
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    left, right = _signs(lattice)
    normalizer = FiniteSubgroupNormalizer(lattice, (left, right))
    assert _order(lattice, normalizer.centralizer_generators()) == 16
    assert _order(lattice, normalizer.normalizer_generators()) == 32
    swap = normalizer.twisted_transporter((right, left))
    assert swap is not None and swap * left * ~swap == right


def test_klein_four_on_non_isometric_summands():
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)("U").twist(2)
    left, right = _signs(lattice)
    normalizer = FiniteSubgroupNormalizer(lattice, (left, right))
    assert normalizer.twisted_transporter((right, left)) is None
    assert _order(lattice, normalizer.normalizer_generators()) == 16


def test_glued_eigenlattices_of_the_hyperbolic_plane():
    lattice = Lattices(ZZ)("U")
    e, f = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    swap = lattice.O()(dict(zip(labels, (f, e), strict=True)))
    minus = lattice.O()(dict(zip(labels, (-e, -f), strict=True)))
    normalizer = FiniteSubgroupNormalizer(lattice, (swap, minus))
    assert _order(lattice, normalizer.centralizer_generators()) == 4
    assert _order(lattice, normalizer.normalizer_generators()) == 4
    assert normalizer.twisted_transporter((minus * swap, minus)) is None


def test_order_six_rotation_on_an_indefinite_sum_of_hexagonal_planes():
    r"""L = A2 + A2(-1), indefinite of signature (2, 2), F = (rho, -rho) with rho of order 3.

    F has no eigenvalue +-1, so R = L is indefinite. The characteristic polynomials
    x^2 + x + 1 and x^2 - x + 1 differ, so a centralizing isometry preserves both planes;
    in the dihedral group O(A2) of order 12 the centralizer of a rotation of order 3 is the
    six rotations. So C(H) has order 36. Aut(H) = {F -> F, F -> F^-1}, and the inverse is
    realized by a reflection in each plane, so N(H) has order 72.
    """
    lattice = Lattices(ZZ)([[int(2), int(-1), int(0), int(0)], [int(-1), int(2), int(0), int(0)],
                            [int(0), int(0), int(-2), int(1)], [int(0), int(0), int(1), int(-2)]])
    a1, a2, b1, b2 = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    rotation = lattice.O()(dict(zip(labels, (a2, -a1 - a2, -b2, b1 + b2), strict=True)))
    normalizer = FiniteSubgroupNormalizer(lattice, (rotation,))
    assert not normalizer.rest_definite
    assert _order(lattice, normalizer.centralizer_generators()) == 36
    assert _order(lattice, normalizer.normalizer_generators()) == 72
    witness = normalizer.twisted_transporter((~rotation,))
    assert witness is not None and witness * rotation * ~witness == ~rotation
