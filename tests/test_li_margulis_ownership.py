"""The Li–Margulis finite box must return actual owned isometries."""

from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.li_margulis import (
    _vectors_in_box,
    isometries_in_box,
    is_isometric_in_box,
    primitive_tuples_in_box,
)


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


def test_theorem_bounds_allow_unbounded_integer_magnitude_without_materializing_box():
    r"""An arbitrary-precision finite bound must not be converted to a C-sized product pool.

    The Li–Margulis theorem gives numerically finite but enormous upper
    bounds. A finite algorithm is not executable if range/product overflows
    or if it must construct the whole box before emitting its first witness.
    """
    from itertools import islice

    ambient = Lattices(ZZ)("U")
    bound = ZZ(2)**100
    prefix = tuple(islice(_vectors_in_box(ambient, bound), 6))
    assert len(prefix) == 6 and prefix[0] == ambient.zero()
    assert len({tuple(v) for v in prefix}) == 6

    complete_small_box = tuple(_vectors_in_box(ambient, 2))
    assert len(complete_small_box) == 25
    assert {tuple(v) for v in complete_small_box} == {
        (a, b) for a in range(-2, 3) for b in range(-2, 3)
    }

    first_null = next(primitive_tuples_in_box(ambient, ((0,),), bound))
    assert len(first_null) == 1
    assert first_null[0].q() == 0 and first_null[0].content() == 1
    assert len(tuple(primitive_tuples_in_box(ambient, ((0,),), 2))) == 4

    ternary = Lattices(ZZ)([[1, 0, 0], [0, 1, 0], [0, 0, -1]])
    ternary_box = tuple(_vectors_in_box(ternary, 1))
    assert len(ternary_box) == 27
    assert len({tuple(v) for v in ternary_box}) == 27
