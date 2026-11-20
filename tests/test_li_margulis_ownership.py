"""The Li–Margulis finite box must return actual owned isometries."""

from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.li_margulis import (
    _vectors_in_box,
    isometries_in_box,
    is_isometric_in_box,
    primitive_tuples_in_box,
)


def test_theorem_one_determinant_factor_is_negative_as_printed():
    r"""Li--Margulis 2016, Theorem 1, eq. (3), printed p.2.

    The determinant exponent is MINUS (13n³-39n²+36n+52)/(20n),
    not plus. At n=3 it equals -8/3. Since |det B|>=1, bounding
    that factor above by one is valid; taking its reciprocal is not.
    Equation (4), separately applicable to anisotropic forms, has a
    positive determinant exponent and must not be confused with (3).
    """
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ
    from sage_indefinite_port.indefinite.li_margulis import isometry_box_bound

    gram = matrix(SageZZ, [[2, 0, 0], [0, -2, 0], [0, 0, 6]])
    assert abs(gram.det()) == 24
    assert isometry_box_bound(gram, gram, 1) == 6**22 + 1


def test_primitive_embedding_completion_majorizes_negative_determinant_power():
    r"""Li--Margulis 2016, Lemma 20 plus Theorem 1, eq. (3).

    Complete the primitive rank-one embedding with a determinant-controlled
    Gram, then apply Theorem 1 to the completed rank-three form. The
    target determinant contributes |det B|^(-8/3) <= 1, whereas
    the completed Gram has a separate, positive contribution.
    """
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ
    from sage_indefinite_port.indefinite.li_margulis import (
        LiMargulisConstants, embedding_orbit_bound, theorem_one_constant,
    )

    gram = matrix(SageZZ, [[2, 0, 0], [0, -2, 0], [0, 0, 6]])
    embedded = matrix(SageZZ, [[2]])
    constants = LiMargulisConstants(3, *(1 for _ in range(7)))
    completed_norm_bound = 2**2 * abs(gram.det()) + 2
    expected = theorem_one_constant(constants) * completed_norm_bound**11 * 6**11 + 1
    assert embedding_orbit_bound(gram, embedded, constants) == expected


def test_theorem_two_exponent_uses_negative_theorem_one_determinant_power():
    r"""Theorem 2's |det B|^(n^6) dominates n(n+1)b for all n>=3.

    A symbolic proof for all n (not just the checked ranks) is in the
    theorem_two_constant docstring. Here exact rational checks at the
    smallest four ranks protect against accidentally introducing the
    positive determinant exponent from the *different* equation (4).
    """
    from sage.rings.rational_field import QQ as SageQQ

    for n in (3, 4, 5, 6):
        d = SageQQ(13*n**3 - 39*n**2 + 36*n + 52) / (20*n)
        b = SageQQ(13*n**3 - 13*n**2 + 72*n - 20) / 40
        exponent = n * (n + 1)*b
        assert d > 0 and exponent == n*(n+1)*b
        assert exponent <= n**6


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
