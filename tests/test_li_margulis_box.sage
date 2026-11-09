r"""Box searches, against known orbit counts.

U: the primitive null vectors are +-e, +-f, one O(U)-orbit. U + A1, square -2: by Eichler
two orbits, divisors 1 and 2. A tuple (e, f) in U + U: primitive pairs with Gram U, one orbit.
"""

from sage_extensions.all import *
from sage_indefinite_port.indefinite.li_margulis import embedding_orbits_in_box, is_isometric_in_box, primitive_tuples_in_box


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
