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
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.finite_normalizers import FiniteSubgroupNormalizer, orthogonal_decomposition_normalizer


def _unreached_transporter(gram, sources, basis, target_basis):
    raise AssertionError("this normalizer does not consult the commutant transporter")


def _from_columns(lattice, columns):
    frame = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    images = {}
    for j, label in enumerate(labels):
        images[label] = sum(
            (lattice.scalar_multiple(int(columns[i][j]), g) for i, g in enumerate(frame)), lattice.zero()
        )
    return lattice.O()(images)


def _reduction_closure(lattice, isometries):
    frame = tuple(lattice.module_generators())
    images = [matrix(GF(int(3)), [[int(c) for c in g(x)] for x in frame]) for g in isometries]
    identity = matrix.identity(GF(int(3)), len(frame))
    identity.set_immutable()
    closure = set([identity])
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for image in images:
            following = current * image
            following.set_immutable()
            if following not in closure:
                closure.add(following)
                frontier.append(following)
    return closure


def test_arithmetic_normalizer_of_the_full_orthogonal_group_of_two_hyperbolic_planes():
    r"""L = U + U, H = <O(L), E> = O(L) with E = E(e1, e2 + f2) unipotent. H contains SO(L) cap Gamma(3) and is
    Zariski dense in O(L) (log E spans an Ad-orbit generating so(4) = sl2 + sl2, whose factors are swapped by the
    summand swap), its congruence kernel is not abelian, and N(H) = O(L). So the returned generators lie in O(L)
    and reduce modulo 3 onto the image of O(L), which contains -1 and the summand swap."""
    from sage_indefinite_port.indefinite.finite_normalizers import normalizer_generators

    lattice = Lattices(ZZ)([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    eichler = _from_columns(lattice, [[1, -1, -1, -1], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]])
    ambient = tuple(lattice.O().group_generators())
    generators = normalizer_generators(lattice, (*ambient, eichler), int(3), orthogonal_decomposition_normalizer, _unreached_transporter)
    target = _reduction_closure(lattice, ambient)
    assert _reduction_closure(lattice, generators) == target
    sigma = _from_columns(lattice, [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
    assert _reduction_closure(lattice, (sigma,)) <= target


def test_orthogonal_pair_owner_gives_the_normalizer_of_sl2_times_sl2_on_two_hyperbolic_planes():
    r"""L = U + U = (M_2(Z), det), SO(L)^0 = SL2 x SL2 / +-1 acting by X -> a X b^{-1}, so g = sl2 + sl2 = so(L_Q),
    one isotypic piece, centralizer c = 0. H = <O(L), E(e1, e2 + f2)> = O(L) contains G(Z) cap Gamma(3). Its
    normalizer is O(L): the owner returns O(L) for N(g), and the returned generators reduce modulo 3 onto the image
    of O(L), which contains the summand swap."""
    from sage_indefinite_port.indefinite.finite_normalizers import arithmetic_normalizer_generators

    lattice = Lattices(ZZ)([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    eichler = _from_columns(lattice, [[1, -1, -1, -1], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]])
    ambient = tuple(lattice.O().group_generators())
    generators = arithmetic_normalizer_generators(lattice, (*ambient, eichler), int(3), orthogonal_decomposition_normalizer)
    target = _reduction_closure(lattice, ambient)
    assert _reduction_closure(lattice, generators) == target
    sigma = _from_columns(lattice, [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
    assert _reduction_closure(lattice, (sigma,)) <= target
