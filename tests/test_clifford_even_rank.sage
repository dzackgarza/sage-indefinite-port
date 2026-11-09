r"""Centres of even Clifford algebras. By hand z^2 = (-1)^(m(2m-1)) prod q(e_i), an invariant modulo squares.

With orthogonal bases e+f, e-f of U, q = 2, -2. So 3U gives (-1)^15 (-4)^3 = 64,
split with components of dimension 32/2 = 16. diag(2,2,-6,-6) gives 144, split with
components of dimension 8/2 = 4. 2U + <1,-2> gives (-1)^15 (16)(-2) = 32, whose
square class is 2, so the discriminant algebra is the field QQ(sqrt 2).
"""

from sage_extensions.all import *
from sage.matrix.constructor import diagonal_matrix, block_diagonal_matrix, matrix as sage_matrix
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.clifford_even_rank import EvenCliffordCentre


def _u():
    return sage_matrix(SageQQ, [[SageQQ(int(0)), SageQQ(int(1))], [SageQQ(int(1)), SageQQ(int(0))]])


def _diagonal(*entries):
    return diagonal_matrix(SageQQ, [SageQQ(int(x)) for x in entries])


def test_three_hyperbolic_planes_split_into_two_degree_four_components():
    centre = EvenCliffordCentre(block_diagonal_matrix([_u(), _u(), _u()], subdivide=False))
    assert (centre.volume_square() / SageQQ(int(64))).is_square()
    assert centre.is_split()
    first, second = centre.idempotents()
    assert first * first == first and second * second == second
    assert first * second == centre.algebra().zero()
    assert first + second == centre.algebra().one()
    assert centre.volume_is_central_in_even_part()
    assert centre.component_dimension(first) == int(16)
    assert centre.component_dimension(second) == int(16)


def test_the_anisotropic_rank_four_form_splits_into_two_quaternion_components():
    centre = EvenCliffordCentre(_diagonal(2, 2, -6, -6))
    assert (centre.volume_square() / SageQQ(int(144))).is_square()
    first, second = centre.idempotents()
    assert centre.component_dimension(first) == int(4)
    assert centre.component_dimension(second) == int(4)


def test_a_nonsplit_rank_six_form_has_a_quadratic_discriminant_field():
    centre = EvenCliffordCentre(block_diagonal_matrix([_u(), _u(), _diagonal(1, -2)], subdivide=False))
    assert (centre.volume_square() / SageQQ(int(32))).is_square()
    assert not centre.is_split()
    assert centre.discriminant_field().discriminant() == int(8)
