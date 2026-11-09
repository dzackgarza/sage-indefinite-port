r"""Issue 61, item 4, rank-two sources on U + U: the zero Gram matrix and diag(2, -2) each give one orbit in the box."""

from sage_extensions.all import *
from sage_indefinite_port.indefinite.li_margulis import embedding_orbits_in_box


def test_prescribed_rank_two_gram_fixtures():
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    assert len(embedding_orbits_in_box(L, [[0, 0], [0, 0]], 1)) == int(1)
    assert len(embedding_orbits_in_box(L, [[2, 0], [0, -2]], 1)) == int(1)
