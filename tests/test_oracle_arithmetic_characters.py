"""Arithmetic orthogonal-group oracles, independently checked.

Pinned source: references/vendor/Oscar.jl@d135b70b/test/Groups/spinor_norms.jl,
``@testset \"spinor norms\"``, ``from_sage`` table and the final ``U + <-3>``
example. The Hecke spinor-*genus* tests in
references/vendor/Hecke.jl@e2ab5716/test/QuadForm/Quad/ZGenus.jl,
lines 432--479, concern local genus classes, not an O(L) character;
none of their numerical expectations are transplanted here.
"""

from tests.acceptance.consumer import lattice, isometry_from_rows


def test_oscar_rank_three_discriminant_image_is_proper():
    """OSCAR ``from_sage`` first row: |O(q)|=8, |im O(L)|=4.

    The matrix is the symmetric Gram [[2,-1,0],[-1,-2,0],[0,0,-10]].
    The strict inequality proves the image is a proper subgroup, not merely
    that a collection of listed generators has four elements.
    """
    L = lattice(((2, -1, 0), (-1, -2, 0), (0, 0, -10)))
    assert int(L.discriminant_group().orthogonal_group().cardinality()) == 8
    assert int(L.O().discriminant_image().cardinality()) == 4


def test_hyperbolic_plane_real_spinor_character_has_negative_value():
    r"""Direct Cartan--Dieudonne calculation in U, independent of port output.

    Let e,f have (e,f)=1. The swap s(e)=f,s(f)=e is reflection
    in e-f, of norm -2; the normalized real spinor value -v²/2
    is +1. The negative swap t(e)=-f,t(f)=-e is reflection in
    e+f, of norm +2, so its normalized real spinor value is -1.
    Thus t is NOT in O^+(U), although both have determinant -1;
    in particular the real-spinor kernel is a proper subgroup of O(U).
    """
    L = lattice(((0, 1), (1, 0)))
    positive = isometry_from_rows(L, ((0, 1), (1, 0)))
    negative = isometry_from_rows(L, ((0, -1), (-1, 0)))
    character = L.O_plus().preimage_morphism()
    assert character(positive).is_one()
    assert not character(negative).is_one()
    assert int(positive.determinant()) == int(negative.determinant()) == -1
