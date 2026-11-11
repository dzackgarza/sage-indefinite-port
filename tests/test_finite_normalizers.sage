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


def _from_columns(lattice, columns):
    frame = tuple(lattice.module_generators())
    labels = tuple(lattice.module_generating_set())
    images = {}
    for j, label in enumerate(labels):
        images[label] = sum(
            (lattice.scalar_multiple(int(columns[i][j]), g) for i, g in enumerate(frame)), lattice.zero()
        )
    return lattice.O()(images)


def test_pell_isometry_normalizer_contains_the_inversion():
    r"""L = <1> + <-2>. f = (3 + 2 sqrt 2) on x + y sqrt 2, columns (3, 2) and (4, 3), has infinite order.
    The reflection y -> -y is conjugation of Z[sqrt 2], so it conjugates f to f^{-1}: N(<f>) = C(f) . 2."""
    from sage_indefinite_port.indefinite.finite_normalizers import has_finite_order, infinite_cyclic_normalizer, normalizer_generators

    lattice = Lattices(ZZ)([[1, 0], [0, -2]])
    f = _from_columns(lattice, [[3, 4], [2, 3]])
    reflection = _from_columns(lattice, [[1, 0], [0, -1]])
    assert reflection * f * ~reflection == ~f
    assert not has_finite_order(matrix(SageZZ, [[3, 4], [2, 3]]))
    centralizer, inversion = infinite_cyclic_normalizer(lattice, f)
    assert inversion is not None and inversion * f * ~inversion == ~f
    assert all(g * f == f * g for g in centralizer)
    assert any(g * f * ~g == ~f for g in normalizer_generators(lattice, (f,), int(3)))


def test_unipotent_eichler_normalizer_contains_the_inversion():
    r"""U + U with basis (e1, f1, e2, f2) and F = E(e1, a), a = e2 + f2. E(e1, a)^{-1} = E(e1, -a), and
    t = 1 + (-1) fixes e1 and sends a to -a, so t F t^{-1} = F^{-1}."""
    from sage_indefinite_port.indefinite.finite_normalizers import has_finite_order, infinite_cyclic_normalizer

    lattice = Lattices(ZZ)([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    eichler = [[1, -1, -1, -1], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]]
    F = _from_columns(lattice, eichler)
    t = _from_columns(lattice, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    assert t * F * ~t == ~F
    assert not has_finite_order(matrix(SageZZ, eichler))
    centralizer, inversion = infinite_cyclic_normalizer(lattice, F)
    assert inversion is not None and inversion * F * ~inversion == ~F
    assert all(g * F == F * g for g in centralizer)


def test_free_abelian_normalizer_of_two_pell_units():
    r"""L = (<1> + <-2>) + (<1> + <-3>), H = <F1, F2> with F1 = (3 + 2 sqrt 2) + 1 and F2 = 1 + (2 + sqrt 3),
    free abelian of rank 2. The joint characters have logarithm vectors (+-log(3 + 2 sqrt 2), 0) and
    (0, +-log(2 + sqrt 3)), of different lengths, so an automorphism F_i -> prod F_j^{Phi_ij} permuting them is
    Phi = diag(+-1, +-1). Each is realized: y -> -y on a summand conjugates its unit to the inverse and fixes the
    other summand. So N(H) / C(H) = (Z/2)^2."""
    from sage_indefinite_port.indefinite.finite_normalizers import normalizer_generators

    lattice = Lattices(ZZ)([[1, 0, 0, 0], [0, -2, 0, 0], [0, 0, 1, 0], [0, 0, 0, -3]])
    F1 = _from_columns(lattice, [[3, 4, 0, 0], [2, 3, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    F2 = _from_columns(lattice, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 2, 3], [0, 0, 1, 2]])
    assert F1 * F2 == F2 * F1
    signs = set()
    for g in normalizer_generators(lattice, (F1, F2), int(3)):
        pair = []
        for F in (F1, F2):
            image = g * F * ~g
            assert image in (F, ~F)
            pair.append(int(1) if image == F else int(-1))
        signs.add(tuple(pair))
    closure = set([(int(1), int(1))])
    frontier = list(closure)
    while frontier:
        current = frontier.pop()
        for s in signs:
            moved = (current[0] * s[0], current[1] * s[1])
            if moved not in closure:
                closure.add(moved)
                frontier.append(moved)
    assert closure == set([(int(1), int(1)), (int(-1), int(1)), (int(1), int(-1)), (int(-1), int(-1))])


def test_infinite_dihedral_normalizer_is_the_orthogonal_group_of_the_pell_plane():
    r"""L = <1> + <-2>, H = <f, r> with f = 3 + 2 sqrt 2 and r the conjugation y -> -y, so r f r = f^{-1} and H is
    infinite dihedral, non-abelian. O(L) = {+-f^k, +-r f^k}, since an isometry of x^2 - 2 y^2 is multiplication by a
    unit of norm one, +-(3 + 2 sqrt 2)^k, possibly composed with conjugation. H has index 2 in O(L), so N(H) = O(L).
    The (1, 1) entry of f^k and of r f^k is positive and that of -f^k and -r f^k is negative, so g lies in H exactly
    when that entry is positive. Hence N(H) normalizes H and contains an element with negative (1, 1) entry."""
    from sage_indefinite_port.indefinite.finite_normalizers import normalizer_generators

    lattice = Lattices(ZZ)([[1, 0], [0, -2]])
    f = _from_columns(lattice, [[3, 4], [2, 3]])
    r = _from_columns(lattice, [[1, 0], [0, -1]])
    assert r * f * ~r == ~f and f * r != r * f
    frame = tuple(lattice.module_generators())
    corner = lambda g: int(list(g(frame[int(0)]))[int(0)])
    generators = normalizer_generators(lattice, (f, r), int(3))
    for g in generators:
        assert g * f * ~g in (f, ~f)
        assert corner(g * r * ~g) > 0
    assert any(corner(g) < 0 for g in generators)


def test_congruence_normalizer_of_the_summand_swap_modulo_three():
    r"""L = U + U with basis (e1, f1, e2, f2), H = <sigma> Gamma_L(3) with sigma the swap of the two summands.
    H is a full preimage of pi(H) = {1, sigma mod 3}, and Gamma_L(3) is normal in O(L), so g normalizes H exactly
    when pi(g) normalizes a group of order 2, i.e. centralizes sigma mod 3. The diagonal (a, a) for a in O(U)
    commutes with sigma, the Eichler cube E(e1, 3(e2 + f2)) lies in Gamma_L(3), and (s, 1) with s the swap e1 <-> f1
    conjugates sigma to the map e1 -> f2, f1 -> e2, which is neither 1 nor sigma modulo 3."""
    from sage_indefinite_port.indefinite.finite_normalizers import congruence_normalizer_generators

    lattice = Lattices(ZZ)([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    sigma = _from_columns(lattice, [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
    diagonal_swap = _from_columns(lattice, [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    diagonal_minus = _from_columns(lattice, [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    eichler = matrix(SageZZ, [[1, -1, -1, -1], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]])
    eichler_cube = _from_columns(lattice, eichler ** int(3))
    left_swap = _from_columns(lattice, [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    frame = tuple(lattice.module_generators())
    reduce = lambda g: matrix(GF(int(3)), [[int(c) for c in g(x)] for x in frame])
    generators = congruence_normalizer_generators(lattice, (sigma,), int(3))
    for g in generators:
        assert reduce(g * sigma * ~g) == reduce(sigma)
    images = [reduce(g) for g in generators]
    identity = matrix.identity(GF(int(3)), int(4))
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
    def contains(g):
        m = reduce(g)
        m.set_immutable()
        return m in closure
    assert all(contains(g) for g in (sigma, diagonal_swap, diagonal_minus, eichler_cube))
    assert not contains(left_swap)


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


def test_torus_of_a_pell_unit_is_the_special_orthogonal_group():
    r"""L = <1> + <-2>, f = 3 + 2 sqrt 2 on columns. Its closure is SO(L)^0, a one-dimensional split torus, whose Lie
    algebra so(L) is spanned by X = [[0, 2], [1, 0]] (X^t B + B X = 0 for B = diag(1, -2)). -1 has finite order and
    trivial torus."""
    from sage_indefinite_port.indefinite.finite_normalizers import _torus_algebra

    algebra = _torus_algebra(matrix(SageZZ, [[3, 4], [2, 3]]))
    assert len(algebra) == int(1)
    X = matrix(SageQQ, [[0, 2], [1, 0]])
    assert matrix(SageQQ, [algebra[int(0)].list(), X.list()]).rank() == int(1)
    assert _torus_algebra(matrix(SageZZ, [[-1, 0], [0, -1]])) == []


def test_torus_of_two_independent_pell_units_is_two_dimensional():
    r"""On (<1> + <-2>) + (<1> + <-3>), F = (3 + 2 sqrt 2) + (2 + sqrt 3). The units 3 + 2 sqrt 2 and 2 + sqrt 3 of
    Q(sqrt 2, sqrt 3) are multiplicatively independent (they generate different quadratic fields), so the eigenvalues
    have only the relations mu * mu' = 1 in each block and the closure of <F> is a two-dimensional torus."""
    from sage_indefinite_port.indefinite.finite_normalizers import _torus_algebra

    F = matrix(SageZZ, [[3, 4, 0, 0], [2, 3, 0, 0], [0, 0, 2, 3], [0, 0, 1, 2]])
    assert len(_torus_algebra(F)) == int(2)


def test_orthogonal_pair_owner_gives_the_normalizer_of_sl2_times_sl2_on_two_hyperbolic_planes():
    r"""L = U + U = (M_2(Z), det), SO(L)^0 = SL2 x SL2 / +-1 acting by X -> a X b^{-1}, so g = sl2 + sl2 = so(L_Q),
    one isotypic piece, centralizer c = 0. H = <O(L), E(e1, e2 + f2)> = O(L) contains G(Z) cap Gamma(3). Its
    normalizer is O(L): the owner returns O(L) for N(g), and the returned generators reduce modulo 3 onto the image
    of O(L), which contains the summand swap."""
    from sage_indefinite_port.indefinite.finite_normalizers import arithmetic_normalizer_generators

    lattice = Lattices(ZZ)([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    eichler = _from_columns(lattice, [[1, -1, -1, -1], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1]])
    ambient = tuple(lattice.O().group_generators())
    generators = arithmetic_normalizer_generators(lattice, (*ambient, eichler), int(3), None)
    target = _reduction_closure(lattice, ambient)
    assert _reduction_closure(lattice, generators) == target
    sigma = _from_columns(lattice, [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
    assert _reduction_closure(lattice, (sigma,)) <= target
