r"""Box searches, against known orbit counts.

U: the primitive null vectors are +-e, +-f, one O(U)-orbit. U + A1, square -2: by Eichler
two orbits, divisors 1 and 2. A tuple (e, f) in U + U: primitive pairs with Gram U, one orbit.
"""

from sage_extensions.all import *
from sage_indefinite_port.indefinite.li_margulis import isotropic_pair_grams, generator_bound, isometry_box_bound, embedding_orbits_in_box, is_isometric_in_box, primitive_tuples_in_box


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
