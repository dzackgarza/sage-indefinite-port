"""The Li–Margulis finite box must return actual owned isometries."""

from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.li_margulis import isometries_in_box, is_isometric_in_box


def test_bounded_search_returns_owned_orthogonal_and_isom_arrows():
    r"""A checked Gram identity is not a substitute for the O(L)/Isom parent."""
    source = Lattices(ZZ)("U")
    target = Lattices(ZZ)([[0, 1], [1, 0]])
    found = tuple(isometries_in_box(source, 1))
    assert len(found) == 4
    assert all(g.parent() is source.O() for g in found)
    assert all(g.domain() is source and g.codomain() is source for g in found)
    witness = is_isometric_in_box(source, target, 1)
    assert witness is not None
    assert witness.parent() is source.Isom(target)
    assert witness.domain() is source and witness.codomain() is target
