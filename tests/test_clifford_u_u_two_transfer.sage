r"""Transfer of Sp_4(Z) to O(U + U + <2>), with b the form, Q = b/2, and theta(s_v) = Q(v).

Sp_4 acts on the five-dimensional pfaffian space, the complement of Omega = e12 - e34 in 3U.
Derived by hand: an element acting trivially also fixes Omega, so it acts trivially on
Lambda^2 and is scalar with c^2 = 1: the kernel on Sp_4(Z) is {1, -1}. The image is
SO(L) & ker theta (the module docstring proves integrality by Cartan form).
Characters (det, theta): s_{e+f} is (-1, 1), s_{e-f} is (-1, -1), their product is (1, -1).
-id on L is (-1, 1): theta(-id) = (1)(-1)(1)(-1)(1) = 1, and det = -1 in odd rank.
Modulo 3: Sp_4(3) of order 51840 acts with kernel {1, -1}, so the image has order 25920, and the
three cosets stay distinct since -1 is a nonsquare mod 3: the full image is O_5(3), of order 4 * 25920 = 103680.
"""

from sage_extensions.all import *
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.finite_rings.finite_field_constructor import GF as SageGF
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.clifford_rank_five import (
    characters,
    hyperbolic_hyperbolic_two,
    orthogonal_group_generators,
    pfaffian_action,
    spinor_kernel_generators,
    symplectic_generators,
)
from sage_indefinite_port.indefinite.clifford_split_rank_six import hyperbolic_coset_representatives, isometry_from_matrix


def _order_mod_3(isometries, lattice):
    points = tuple(lattice.module_generators())
    field = SageGF(int(3))
    images = [
        sage_matrix(field, [[int(entry) for entry in f(point).to_vector()] for point in points]).transpose()
        for f in isometries
    ]
    return int(libgap.Group([libgap(m) for m in images]).Size())


def test_the_lattice_has_the_gram_of_u_u_two():
    lattice = hyperbolic_hyperbolic_two()
    points = tuple(lattice.module_generators())
    gram = [[int(lattice.b(x, y)) for y in points] for x in points]
    assert gram == [[0, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 2]]


def test_symplectic_images_are_isometries_in_the_spinor_kernel():
    lattice = hyperbolic_hyperbolic_two()
    points = tuple(lattice.module_generators())
    for f in spinor_kernel_generators(lattice):
        for x in points:
            for y in points:
                assert lattice.b(f(x), f(y)) == lattice.b(x, y)
        assert characters(f, lattice) == (int(1), int(1))


def test_the_kernel_on_sp4_is_plus_minus_one():
    minus = -sage_matrix.identity(SageQQ, int(4))
    assert pfaffian_action(minus) == sage_matrix.identity(SageQQ, int(5))
    assert all(pfaffian_action(g) != sage_matrix.identity(SageQQ, int(5)) for g in symplectic_generators())


def test_coset_representatives_have_their_claimed_characters():
    lattice = hyperbolic_hyperbolic_two()
    representatives = hyperbolic_coset_representatives(lattice)
    assert sorted(representatives) == [(int(-1), int(-1)), (int(-1), int(1)), (int(1), int(-1))]
    for value, element in representatives.items():
        assert characters(element, lattice) == value


def test_minus_identity_has_determinant_minus_one_and_trivial_spinor_norm():
    lattice = hyperbolic_hyperbolic_two()
    assert characters(isometry_from_matrix(-sage_matrix.identity(SageQQ, int(5)), lattice), lattice) == (int(-1), int(1))


def test_image_modulo_three_has_the_hand_derived_orders():
    lattice = hyperbolic_hyperbolic_two()
    assert _order_mod_3(spinor_kernel_generators(lattice), lattice) == int(25920)
    assert _order_mod_3(orthogonal_group_generators(lattice), lattice) == int(103680)
