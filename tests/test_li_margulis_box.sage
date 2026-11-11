r"""Box searches, against known orbit counts.

U: the primitive null vectors are +-e, +-f, one O(U)-orbit. U + A1, square -2: by Eichler
two orbits, divisors 1 and 2. A tuple (e, f) in U + U: primitive pairs with Gram U, one orbit.
"""

from sage_extensions.all import *
import pytest
from sage_indefinite_port.indefinite.li_margulis import (
    LiMargulisConstants, embedding_orbit_bound, embedding_orbits_in_box, generator_bound, is_isometric_in_box,
    isometries_in_box, isometry_box_bound, isotropic_pair_grams, isotropic_vector_orbits_in_box, lemma1_constant,
    primitive_tuples_in_box, theorem_one_anisotropic_constant, theorem_one_constant, theorem_two_constant,
)


def test_null_vectors_of_u():
    L = Lattices(ZZ)("U")
    assert len(tuple(primitive_tuples_in_box(L, [[0]], 2))) == int(4)
    assert len(embedding_orbits_in_box(L, [[0]], 2)) == int(1)


def test_two_orbits_of_square_minus_two():
    L = Lattices(ZZ)("U") + Lattices(ZZ)("A1")
    orbits = embedding_orbits_in_box(L, [[-2]], 2)
    assert len(orbits) == int(2)
    assert sorted(int(v[0].divisor()) for v in orbits) == [int(1), int(2)]


def test_rank_two_tuple():
    L = Lattices(ZZ)("U") + Lattices(ZZ)("A1")
    assert len(embedding_orbits_in_box(L, [[0, 1], [1, 0]], 1)) == int(1)


def test_isometry_found_and_absent():
    U = Lattices(ZZ)("U")
    other = Lattices(ZZ)([[0, 1], [1, 0]])
    assert is_isometric_in_box(U, other, 1) is not None
    assert is_isometric_in_box(U, Lattices(ZZ)([[2, 0], [0, -2]]), 3) is None


def test_generator_bound_formula():
    r"""r = 3, B = diag(1, 2, -1): |det|^729 and ||B||^243 with C = 1."""
    gram = [[1, 0, 0], [0, 2, 0], [0, 0, -1]]
    expected = 2 ** int(729) * 2 ** int(243) + 1
    assert int(generator_bound(gram, int(1))) == int(expected)


def test_isometry_bound_exponents():
    r"""n = 3, A = B = diag(1, 2, -1), |det| = 2, norm 2, C = 1.

    (3): exponents 0, (351 - 117 + 156 + 20)/40 = 41/4 and (351 - 117 + 216 - 20)/40 = 43/4,
    with ceilings 11 and 11. (4): (117 + 24 - 26)/30 = 23/6, 1/2, 1, with ceilings 4, 1, 1.
    """
    gram = [[1, 0, 0], [0, 2, 0], [0, 0, -1]]
    assert int(isometry_box_bound(gram, gram, int(1))) == int(2 ** int(22) + 1)
    assert int(isometry_box_bound(gram, gram, int(1), anisotropic=True)) == int(2 ** int(6) + 1)


def test_prescribed_null_gram_on_two_planes():
    r"""Issue 61, item 4. In U + U, A = (0) has one orbit of primitive null vectors in the box."""
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    assert len(embedding_orbits_in_box(L, [[0]], 1)) == int(1)


def test_isotropic_pair_grams_of_u_plus_u_three():
    r"""L = U + U(3), even with exponent 3: d = 1 gives a = 0, and d = 3 gives a = 0, 2, 4."""
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U").twist(3)
    expected = [((int(0), int(d)), (int(d), int(a))) for d, a in ((1, 0), (3, 0), (3, 2), (3, 4))]
    assert list(isotropic_pair_grams(L)) == expected


def _unit_constants(rank):
    from sage.rings.rational_field import QQ as SageQQ

    one = SageQQ(int(1))
    return LiMargulisConstants(int(rank), one, one, one, one, one, one, one)


def test_theorem_one_constant_composition():
    r"""n = 3, every implicit majorant 1. K_24 = ceil(3^{5/2}) = 16 (15.58...).

    e = (39/5)(1 + 3/2) = 39/2, so K_24 enters with exponent n + 2e = 42 and
    C_1 = n^2 (n-1)! 16^42 = 18 * 16^42. Equation (4): 18 * 16^3.
    """
    assert int(lemma1_constant(int(3))) == int(16)
    assert int(lemma1_constant(int(4))) == int(4) ** int(3)
    constants = _unit_constants(3)
    assert theorem_one_constant(constants) == int(18) * int(16) ** int(42)
    assert theorem_one_anisotropic_constant(constants) == int(18) * int(16) ** int(3)


def test_theorem_two_constant_composition():
    r"""n = 3, unit majorants: n^3 (n-1)! (n C_1)^n = 54 * (3 * 18 * 16^42)^3."""
    constants = _unit_constants(3)
    assert theorem_two_constant(constants) == int(54) * (int(54) * int(16) ** int(42)) ** int(3)


def test_embedding_bound_of_a_line():
    r"""B = diag(1, 1, -1), A = (1), unit majorants: ||C|| <= 1 * 1 + 1 = 2, ||B|| = 1.

    The exponent a = 41/4 of ||A|| in (3) rounds to 11: bound C_1 * 2^11 + 1.
    """
    constants = _unit_constants(3)
    expected = int(18) * int(16) ** int(42) * int(2) ** int(11) + int(1)
    assert int(embedding_orbit_bound([[1, 0, 0], [0, 1, 0], [0, 0, -1]], [[1]], constants)) == expected


def test_definite_form_is_rejected():
    with pytest.raises(ValueError):
        embedding_orbit_bound([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1]], _unit_constants(3))


def test_orthogonal_group_of_u():
    r"""O(U) = {+-1, +-swap}: an isometry sends e to a primitive null vector +-e or +-f."""
    assert len(isometries_in_box(Lattices(ZZ)("U"), int(1))) == int(4)


def test_quaternion_unit_of_issue_40():
    r"""Issue 55, D1: g below satisfies g^t B g = B for B = diag(2, 2, -6, -6), entries at most 3.

    The box of size 3 is the smallest containing g; the theorem's box is out of reach.
    """
    from sage.matrix.constructor import matrix as sage_matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    L = Lattices(ZZ)([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    g = sage_matrix(SageZZ, [[int(c) for c in row] for row in [[2, 0, 3, 0], [0, 2, 0, -3], [1, 0, 2, 0], [0, -1, 0, 2]]])
    frame = tuple(L.module_generators())
    columns = set(tuple(int(c) for c in column) for column in g.columns())
    images = [set(tuple(int(c) for c in phi(x)) for x in frame) for phi in isometries_in_box(L, int(3))]
    assert columns in images


def test_isotropic_orbits_of_u_plus_u2():
    r"""Issue 56, E. L = U + U(2), A_L = u(2) with isotropic elements (1,0), (0,1) swapped by O(L).

    L contains U, so by Eichler primitive isotropic orbits are the pairs (div v, v/div v in A_L):
    div 1 one orbit, div 2 one orbit. Representatives e_1 (with f_1) and e_2 (with f_2) lie in
    the box of size 1, which is why that box suffices here.
    """
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U").twist(2)
    orbits = isotropic_vector_orbits_in_box(L, int(1))
    assert sorted(int(v.divisor()) for v in orbits) == [int(1), int(2)]


def test_isotropic_orbits_of_u_plus_u():
    r"""Unimodular L = U + U: by Eichler one orbit of primitive isotropic vectors, represented by e_1."""
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    assert len(isotropic_vector_orbits_in_box(L, int(1))) == int(1)
