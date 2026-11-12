"""A same-genus non-isometric pair at positive and negative index 2, decided by the port (#68).

L = diag(1, -2, -64, 128) has signature (2, 2). Over Q it is <1, -1> + <-2, 2> (the
squares 64 and 128 = 2 * 64 drop out), two hyperbolic planes, so its rational Witt index
is 2 and the port decides it by the Eichler route, not the low-Witt route. Sage's
Genus(L).representatives() gives two classes in its genus; M below is the second one.
By Eichler-Kneser an indefinite genus of rank at least 3 has one class per spinor genus,
so L and M lie in different spinor genera and are not isometric, while L + U and M + U
are (a unimodular hyperbolic plane makes the 2-adic spinor norm surjective).
"""

from tests.acceptance.consumer import lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion

FIRST = [[1, 0, 0, 0], [0, -2, 0, 0], [0, 0, -64, 0], [0, 0, 0, 128]]
SECOND = [[697, 0, 0, 896], [0, -2, 0, 0], [0, 0, -64, 0], [896, 0, 0, 1152]]


def test_same_genus_pair_at_index_two():
    first, second = lattice(FIRST), lattice(SECOND)
    assert tuple(int(index) for index in first.signature_pair()) == (2, 2)
    assert first.genus() == second.genus()

    assert EichlerRecursion().isometry(first, second) is None
