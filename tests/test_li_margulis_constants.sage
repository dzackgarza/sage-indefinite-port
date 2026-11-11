r"""Majorants of the Li-Margulis constants at n = 3, recomputed by hand.

n = 3: d = dim SL_3 - dim SO(2,1) = 8 - 3 = 5, binom(3, 1) = 3, (n-1)! = 2.
"""

from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.li_margulis import (
    appendix_inverse_epsilon, corollary1_constant, corollary2_constant, corollary3_constant, default_constants,
    haar_ball_inverse_constant, harish_chandra_decay_constant, lemma10_radius, lemma11_constant, lemma14_constant,
    lemma15_constant_power, lemma15_inverse_radius, lemma1_constant, lemma20_constant, lemma24_constant,
    matrix_coefficient_constant, section6_constant, theorem5_compact_constant, theorem5_constant, theorem6_constant,
    zeta_product_majorant,
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


THEOREM5_BASE = int(21) * SageQQ(int(16)) / int(5) * int(72) * SageZZ(int(432)) ** int(3)
P3 = int(27) * SageZZ(int(2)) ** int(24)
R3_SCALED = int(27 * 8 * 9) * SageZZ(int(2)) ** int(22)


def test_matrix_coefficient_constant():
    r"""pi coth pi = pi (1 + 2x/(1-x)), x = e^{-2 pi} < e^{-6} < 1/300 since 2.7^6 > 300 and e > 2.7.
    With pi < 22/7: pi coth pi < (22/7)(1 + 2/299) = 6622/2093 < 16/5."""
    assert SageZZ(int(27)) ** int(6) > int(300) * SageZZ(int(10)) ** int(6)
    assert SageQQ(int(22)) / int(7) * (int(1) + SageQQ(int(2)) / int(299)) == SageQQ(int(6622)) / int(2093)
    assert SageQQ(int(6622)) / int(2093) < SageQQ(int(16)) / int(5)
    assert matrix_coefficient_constant() == SageQQ(int(16)) / int(5)


def test_harish_chandra_decay_constant():
    r"""max (1+t) e^{-t/130} = 130 e^{-129/130} < 130 (10/27)(66/65) = 1320/27, and (1320/27)^25 < 21^32."""
    assert int(130) * SageQQ(int(10)) / int(27) * SageQQ(int(66)) / int(65) == SageQQ(int(1320)) / int(27)
    assert (SageQQ(int(1320)) / int(27)) ** int(25) < SageZZ(int(21)) ** int(32)
    assert harish_chandra_decay_constant() == int(21)


def test_haar_ball_and_theorem5():
    r"""D = dim SO(2,1) = 3: c^{-1} = 2^3 * 3^2 = 72; kappa = 432 = 1/eta_0; ceil(13/5) = 3, ceil(39/5) = 8."""
    assert haar_ball_inverse_constant(N) == int(72)
    assert lemma11_constant(N) == int(432) and int(1) / lemma10_radius(N) == int(432)
    assert theorem5_constant(N) == int(4 * 81) * SageZZ(int(3)) ** int(8) * THEOREM5_BASE ** int(3)
    assert theorem5_compact_constant(N) == int(4 * 81) * THEOREM5_BASE ** int(3)


def test_lemma24_and_corollaries():
    r"""n^n beta'^n (1+beta)^{n-1} d^{n+1} alpha^{3n^2-2n+1}, 3n^2-2n+1 = 22:
    (2,1,1): 27 * 4 * 2^22 = 27 * 2^24; (2,2,1): 27 * 8 * 9 * 2^22; index 5 multiplies by 5^4."""
    assert lemma24_constant(N, int(2), int(1), int(1)) == P3
    assert lemma24_constant(N, int(2), int(2), int(1)) == R3_SCALED
    assert lemma24_constant(N, int(2), int(1), int(5)) == P3 * int(625)
    assert corollary2_constant(N) == int(6) * P3 ** int(3)
    assert corollary3_constant(N) == int(6) * R3_SCALED * P3 ** int(2)


def test_default_constants():
    constants = default_constants(N)
    assert constants.theorem5 == theorem5_constant(N)
    assert constants.theorem5_compact == theorem5_compact_constant(N)
    assert constants.theorem6 == theorem6_constant(N)
    assert constants.lemma18 == section6_constant(N, theorem5_constant(N), theorem6_constant(N))
    assert constants.lemma20 == int(2)
    assert constants.corollary2 == corollary2_constant(N)
    assert constants.corollary3 == corollary3_constant(N)
