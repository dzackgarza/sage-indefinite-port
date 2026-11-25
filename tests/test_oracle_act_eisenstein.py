"""Exact Eisenstein Heisenberg translations from ACT, §7.8, (7.8.1).

Primary source: Allcock–Carlson–Toledo, *The complex hyperbolic geometry
of the moduli space of cubic surfaces*, J. Algebraic Geom. 11 (2002),
author PDF https://web.ma.utexas.edu/users/allcock/research/cubic.pdf,
printed p. 35, §7.8, equation (7.8.1), inspected in the rendered PDF.

The paper uses column vectors and <v|w> = w* A v (printed p. 35,
following (7.7.2)); its commutator notation is UNUSUAL:
[T,T'] = T' T T^-1 T'^-1. Thus its commutator is the inverse of
conventional T T' T^-1 T'^-1 for these central commutators.

For lambda=e1 and mu=omega*e1, write theta=omega-bar(omega).
For lambda=e1 and mu=omega*e1, the displayed matrix gives the
cross term (bar(mu)*lambda-bar(lambda)*mu)/2 = -theta/2.
Consequently T(lambda,theta/2)T(mu,theta/2)
= T(lambda+mu,theta/2). This is independently derived from the
matrix, not inferred from the port. CAUTION: the paper's printed
symbolic product law writes +Im(<lambda'|lambda>), which, with
its printed <v|w>=w*Av and the usual Im(x)=(x-bar(x))/2,
would give the opposite sign. The matrix itself is unambiguous;
the test therefore takes the displayed matrix as its oracle.
The sign follows directly from the (5,4) entry: its extra term is
-bar(lambda)*mu and its re-expression as z is
(bar(mu)*lambda-bar(lambda)*mu)/2 = -theta/2.
"""

from sage.matrix.constructor import matrix
from sage_indefinite_port.indefinite.act_eisenstein import (
    FIELD, FORM, OMEGA, THETA, is_isometry, translation,
)


def test_act_781_integral_heisenberg_law_and_commutator():
    zero = FIELD.zero()
    one = FIELD.one()
    lam = (one, zero, zero)
    mu = (OMEGA, zero, zero)
    sum_lam = (one + OMEGA, zero, zero)
    t = translation(lam, THETA / 2)
    u = translation(mu, THETA / 2)
    identity = matrix.identity(FIELD, 5)

    # Norm(lambda)=Norm(mu)=Norm(lambda+mu)=1, and
    # theta/2 - 1/2 = omega is an Eisenstein integer.
    assert THETA == OMEGA - OMEGA.conjugate()
    assert (THETA - one) / 2 == OMEGA
    for g in (t, u, t * u, u * t, ~t, ~u):
        assert is_isometry(g)
        assert g.conjugate_transpose() * FORM * g == FORM

    # ACT (7.8.1), using the column-vector order of its displayed matrix.
    assert t * u == translation(sum_lam, THETA / 2)
    assert ~t == translation(tuple(-x for x in lam), -THETA / 2)
    assert ~u == translation(tuple(-x for x in mu), -THETA / 2)
    assert t * ~t == identity == u * ~u

    # Conventional commutator versus ACT's explicitly reversed notation.
    central_minus = translation((zero, zero, zero), -THETA)
    central_plus = translation((zero, zero, zero), THETA)
    assert t * u * ~t * ~u == central_minus
    assert u * t * ~u * ~t == central_plus
    assert central_minus != identity
    assert central_minus * t == t * central_minus
    assert central_minus * u == u * central_minus
    assert is_isometry(central_minus)
    assert is_isometry(central_plus)
