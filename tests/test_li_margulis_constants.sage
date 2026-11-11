r"""Majorants of the Li-Margulis constants at n = 3, recomputed by hand.

n = 3: d = dim SL_3 - dim SO(2,1) = 8 - 3 = 5, binom(3, 1) = 3, (n-1)! = 2.
"""

from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.li_margulis import (
    appendix_inverse_epsilon, corollary1_constant, default_constants, lemma10_radius, lemma11_constant,
    lemma14_constant, lemma15_constant_power, lemma15_inverse_radius, lemma1_constant, lemma20_constant, matrix_coefficient_constant,
    section6_constant, theorem6_constant, zeta_product_majorant,
)

N = int(3)
EPSILON_INVERSE = SageZZ(int(4) * int(27) * int(216)) ** int(9) * int(30)
R3 = int(27) * EPSILON_INVERSE**int(2) // int(2) + int(1)


def test_lemma10_radius():
    assert lemma10_radius(N) == SageQQ(int(1)) / int(16 * 9 * 3)


def test_lemma11_constant():
    assert lemma11_constant(N) == int(8 * 27 * 2)


def test_lemma15_radius_and_zeta():
    assert lemma15_inverse_radius() == int(3)
    assert zeta_product_majorant(N) == int(3)


def test_lemma15_constant_power():
    r"""4 n^2 d^{ceil(d/2)} = 4 * 9 * 5^3."""
    assert lemma15_constant_power(N) == int(4500)


def test_appendix_epsilon():
    r"""(4 * 27 * 216)^9 * (27 + 3)."""
    assert appendix_inverse_epsilon(N) == EPSILON_INVERSE


def test_lemma14_constant():
    r"""floor(3^3 eps^{-2} / 2) + 1; eps^{-1} is even."""
    assert lemma14_constant(N) == R3
    assert int(27) * EPSILON_INVERSE**int(2) % int(2) == int(0)


def test_theorem6_constant():
    assert theorem6_constant(N) == SageQQ(int(81) * R3**int(2)) ** int(5) * int(4500)


def test_corollary1_and_lemma20():
    assert corollary1_constant(int(1)) == int(1)
    assert corollary1_constant(int(2)) == int(2)
    assert corollary1_constant(int(3)) == int(3) * R3**int(2)
    assert lemma20_constant(N) == int(2)


def test_section6_constant():
    r"""K_24 = ceil(sqrt(3^5)) = 16; e = (39/5)(1 + 3/2) = 39/2; exponent 3 + 39 = 42; 13/5 -> 3."""
    assert lemma1_constant(N) == int(16)
    assert section6_constant(N, int(1), int(1)) == int(9 * 2) * SageZZ(int(16)) ** int(42)
    assert section6_constant(N, int(5), int(2)) == int(9 * 2 * 5) * SageZZ(int(16)) ** int(42) * int(8)


def test_default_constants():
    constants = default_constants(N, int(5), int(7), int(11), int(13))
    assert constants.theorem6 == theorem6_constant(N)
    assert constants.lemma18 == section6_constant(N, int(5), theorem6_constant(N))
    assert constants.lemma20 == int(2)
    assert (constants.theorem5, constants.theorem5_compact, constants.corollary2, constants.corollary3) == (
        int(5), int(7), int(11), int(13)
    )


def test_matrix_coefficient_constant():
    r"""sum_{m in Z} 1/(1+m^2) = pi coth pi; partial sum to |m| <= 2 is 1 + 2(1/2 + 1/5) = 12/5,
    and the tail is below 2 * int_2^oo dx/x^2 = 1, so pi coth pi < 17/5; the sharper 3.1534 < 16/5."""
    assert matrix_coefficient_constant() == SageQQ(int(16)) / int(5)
