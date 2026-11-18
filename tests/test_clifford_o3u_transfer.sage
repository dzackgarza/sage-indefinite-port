r"""Transfer of SL_4(Z) to O(3U), with b the form, Q = b/2, and theta(s_v) = Q(v).

Derived by hand: Lambda^2 g = 1 forces g scalar, c^2 = 1, so the kernel on SL_4(Z) is {1, -1}.
The image is SO(3U) & ker theta (the module docstring proves integrality by Cartan form).
Characters (det, theta): s_{e+f} is (-1, 1), s_{e-f} is (-1, -1), their product is (1, -1).
-id on 3U is (1, -1): theta(-id) = (1)(-1)(1)(-1)(1)(-1) = -1, so -id is outside the image.
Modulo 3, SL_4(F_3) acts with kernel {1, -1}: the image has order |SL_4(3)| / 2 = 6065280.
Since -1 is a nonsquare mod 3, the three cosets stay distinct: the full image is O_6^+(3),
of order 2 * 3^6 (3^3 - 1)(3^2 - 1)(3^4 - 1) = 24261120 = 4 * 6065280.
"""

from sage_extensions.all import *
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.finite_rings.finite_field_constructor import GF as SageGF
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.clifford_split_rank_six import (
    determinant_and_spinor_values,
    hyperbolic_coset_representatives,
    isometry_from_matrix,
    orthogonal_group_generators,
    second_compound,
    spinor_kernel_generators,
    three_hyperbolic_planes,
)


def _order_mod_3(isometries, lattice):
    points = tuple(lattice.module_generators())
    field = SageGF(int(3))
    images = [
        sage_matrix(field, [[int(entry) for entry in f(point).to_vector()] for point in points]).transpose()
        for f in isometries
    ]
    return int(libgap.Group([libgap(m) for m in images]).Size())


def _minus_identity(size):
    return -sage_matrix.identity(SageQQ, int(size))


def test_elementary_images_are_isometries_in_the_spinor_kernel():
    lattice = three_hyperbolic_planes()
    points = tuple(lattice.module_generators())
    for f in spinor_kernel_generators(lattice):
        for x in points:
            for y in points:
                assert lattice.b(f(x), f(y)) == lattice.b(x, y)
        assert determinant_and_spinor_values(f, lattice) == (int(1), int(1))


def test_the_kernel_on_sl4_is_plus_minus_one():
    assert second_compound(_minus_identity(4)) == sage_matrix.identity(SageQQ, int(6))


def test_coset_representatives_have_their_claimed_characters():
    lattice = three_hyperbolic_planes()
    representatives = hyperbolic_coset_representatives(lattice)
    assert sorted(representatives) == [(int(-1), int(-1)), (int(-1), int(1)), (int(1), int(-1))]
    for value, element in representatives.items():
        assert determinant_and_spinor_values(element, lattice) == value


def test_minus_identity_is_determinant_one_with_nontrivial_spinor_norm():
    lattice = three_hyperbolic_planes()
    assert determinant_and_spinor_values(isometry_from_matrix(_minus_identity(6), lattice), lattice) == (int(1), int(-1))


def test_image_modulo_three_has_the_hand_derived_orders():
    lattice = three_hyperbolic_planes()
    assert _order_mod_3(spinor_kernel_generators(lattice), lattice) == int(6065280)
    assert _order_mod_3(orthogonal_group_generators(lattice), lattice) == int(24261120)
