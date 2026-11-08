"""Exact witnesses for Dawes 2108.06236, section 3.5.

The discriminant calculation and the two plane-stabilizer transvections
support the orbit proof in references/dawes-L2-boundary.md.
"""

from sage.all import QQ

from tests.acceptance.consumer import basis, lattice


def test_L2_isotropic_discriminant_classes():
    isotropic = [
        (a, b)
        for a in range(2)
        for b in range(6)
        if (-QQ(a * a) / 4 - QQ(b * b) / 12).denominator() == 1
    ]
    assert isotropic == [(0, 0), (1, 3)]


def test_L2_plane_stabilizer_transvections():
    ambient = lattice([
        [0, 1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, -2, 0],
        [0, 0, 0, 0, 0, -6],
    ])
    e1, f1, e2, f2, v, w = basis(ambient)
    null = e2 + e2 + f2 + f2 + v + w
    assert null.q() == 0
    assert null.is_primitive()
    assert int(null.div()) == 2
    assert int(e1.div()) == 1
    first_upper = ambient.eichler_transvection(e1, -f2)
    first_lower = ambient.eichler_transvection(e2, -f1)
    assert first_upper(e1) == e1
    assert first_upper(e2) == e2 + e1
    assert first_lower(e1) == e1 + e2
    assert first_lower(e2) == e2
    upper = ambient.eichler_transvection(e1, -e2)
    lower = ambient.eichler_transvection(null, -f1)
    assert upper(e1) == e1
    assert upper(null) == null + e1 + e1
    assert lower(e1) == e1 + null
    assert lower(null) == null
    for transformation in (first_upper, first_lower, upper, lower):
        for x in basis(ambient):
            assert (~transformation)(transformation(x)) == x
            for y in basis(ambient):
                assert transformation(x).b(transformation(y)) == x.b(y)
