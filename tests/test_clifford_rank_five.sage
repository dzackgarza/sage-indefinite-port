r"""Rank-five Clifford models. By hand, s is det(B) modulo squares and lambda = -s.

U + U + <2> has determinant 2, so -lambda is 2 modulo squares and lambda = 1 fails.
<1> + H + H has determinant 1, so -lambda is a square: the required scalar is -1
modulo squares. Negating B multiplies s by (-1)^5, so the two lambdas differ by
-1 modulo squares although the forms are similar.
C_0 has dimension 16. The reversal fixes grades 0 and 4 and negates grade 2, which
gives 1 + 5 = 6 = 4*3/2, the symmetric dimension of a symplectic involution.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.clifford_rank_five import RankFiveCliffordModel


def _gram(rows):
    return sage_matrix(SageQQ, [[SageQQ(int(x)) for x in row] for row in rows])


def _u_u_two():
    return _gram([[0, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 2]])


def _one_h_h():
    return _gram([[1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0]])


def test_the_scaled_pfaffian_recovers_the_original_gram():
    gram = _u_u_two()
    assert RankFiveCliffordModel(gram).recovered_gram() == gram


def test_the_scalar_of_u_u_two_is_not_one_modulo_squares():
    assert (-RankFiveCliffordModel(_u_u_two()).scalar() / SageQQ(int(2))).is_square()
    assert not RankFiveCliffordModel(_u_u_two()).scalar().is_square()


def test_the_split_form_has_scalar_minus_one():
    model = RankFiveCliffordModel(_one_h_h())
    assert (-model.scalar()).is_square()
    assert model.recovered_gram() == _one_h_h()


def test_negated_forms_are_similar_with_scalars_differing_by_minus_one():
    gram = _u_u_two()
    first, second = RankFiveCliffordModel(gram), RankFiveCliffordModel(-gram)
    assert second.recovered_gram() == -gram
    assert (-second.scalar() / first.scalar()).is_square()


def test_the_even_algebra_has_a_symplectic_reversal():
    model = RankFiveCliffordModel(_u_u_two())
    assert model.even_dimension() == int(16)
    assert model.even_symmetric_dimension() == int(6)
