"""Independently proved explicit radii for Li--Margulis, Lemmas 9--11.

Primary source: H. Li and G. A. Margulis, *Effective estimates on
integral quadratic forms*, author PDF ``equivalence.pdf``, section 4.2,
printed pp. 11--12 (Lemmas 9, 10, 11), and section 1.7, printed p. 6.
https://li-han.weebly.com/uploads/5/1/3/3/51338781/equivalence.pdf

The paper does NOT specify these numerical radii. The assertions below
use the following independent derivation in the paper's ENTRYWISE maximum
matrix norm |M| = max_ij |M_ij|, with G=SL_3(R).

For |E| < eta = 1/36, the k-th elementary characteristic coefficient
of u=I+E differs from binomial(3,k) by at most

    binomial(3,k) * sum_{j=1}^k binomial(k,j) j! eta**j
    <= binomial(3,k) * k*eta/(1-k*eta) < 1,   1 <= k <= 3.

The inequality follows by bounding binomial(k,j) j! <= k**j.
Consequently, an integral characteristic polynomial of such a u equals
(t-1)^3, so u is unipotent. The compact-orbit argument in Lemma 10
excludes nonidentity unipotents from the arithmetic stabilizer, using
Godement's criterion exactly as in the source.

For |E_i|<r=1/432, the Neumann series in row-sum norm gives
|(I+E_2)^(-1)|_infty <= (1-3r)^(-1), hence

    |(I+E_2)^(-1)(I+E_1)-I|
       < 6*r/(1-3*r) < 12*r = 1/36.

Thus Lemma 9 applies and Lemma 10 holds with r=1/432.

For Lemma 11, put m=|g|, t=m**3. By the Leibniz determinant
formula, 1=|det g|<=3!m**3=6t, so t>=1/6 (NOT m>=1: the
determinant-one matrix (1/3)[[2,-2,1],[2,1,-2],[1,2,2]] has
entrywise maximum m=2/3). Cramer's rule yields
|g^(-1)| <= 2m**2. If gamma in SL_3(Z) is nonidentity, then
1 <= |gamma-I| <= 9|g^(-1)| |g*gamma*g^(-1)-I| |g|
<= 18t |g*gamma*g^(-1)-I|.
For delta=(432t)^(-1) <=1/72, the same Neumann estimate gives
|g_2^(-1)g_1-I| < 6*delta/(1-3*delta) = 6/(432t-3)
< 1/(18t), because the last inequality is equivalent to
324t>3, which follows from t>=1/6. Thus distinct g_1,g_2 in B_G(delta)
cannot have the same image in G/Gamma_N, proving the stated radius.

No claim about the independently unverified constants of Lemmas 14--15,
Theorems 5--6, or section 6 is made here.
"""

from fractions import Fraction
from math import comb

from sage.rings.rational_field import QQ

from sage_indefinite_port.indefinite.li_margulis import (
    lemma10_radius,
    lemma11_constant,
)


def test_lemma10_explicit_rank_three_injectivity_radius():
    """Source Lemmas 9--10; the proof above certifies this choice."""
    n = 3
    radius = Fraction(1, 432)
    eta = 12 * radius
    assert eta == Fraction(1, 36)
    assert all(
        Fraction(comb(n, k) * k) * eta / (1 - k * eta) < 1
        for k in range(1, n + 1)
    )
    assert Fraction(6) * radius / (1 - 3 * radius) < eta
    assert lemma10_radius(n) == QQ(radius.numerator) / radius.denominator


def test_lemma11_explicit_rank_three_uniform_radius():
    """Source Lemmas 9 and 11; Cramer's lower bound proved above."""
    n = 3
    kappa = 432
    # For g in SL_3(R), t=|g|**3 >= 1/6 by det Leibniz, not t>=1.
    # The relative comparison improves as t grows; its strict inequality
    # for all t>=1/6 is equivalent to 324*t > 3.
    minimum_t = Fraction(1, 6)
    delta = 1 / (kappa * minimum_t)
    assert delta <= Fraction(1, 72)
    assert Fraction(6) * delta / (1 - 3 * delta) < Fraction(1, 18) / minimum_t
    assert 324 * minimum_t > 3
    assert lemma11_constant(n) == kappa
