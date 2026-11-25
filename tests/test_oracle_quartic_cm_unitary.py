"""Exact quartic CM unitary group, with an implementation-independent oracle.

Let K=Q(zeta_5), O=Z[zeta_5], and c=zeta_5+zeta_5^(-1).
Since c != 0, preserving c I_2 is equivalent to preserving I_2.
For every embedding sigma:K->C, the column equations of g^*g=I give
sum_j |sigma(g_ji)|^2=1. Thus each conjugate of each entry a has
modulus <=1. If a != 0, its nonzero integral norm satisfies
1 <= |N_K/Q(a)| = product_sigma |sigma(a)| <=1, so all four
conjugates have modulus exactly one. The remaining entry in that column
vanishes. Kronecker's theorem implies a is a root of unity; the full
root-of-unity group of Q(zeta_5) is mu_10 (Washington, Introduction to
Cyclotomic Fields, 2nd ed., Chapter 1, roots of unity in cyclotomic
fields). Orthogonality forces distinct occupied rows. Conversely all
mu_10-monomial matrices are unitary. Therefore U(c I_2,O) is precisely
mu_10^2 semidirect S_2, of cardinality 10^2*2=200.

The real conjugates of c are 2cos(2pi/5)>0 and 2cos(4pi/5)<0;
each contributes four real dimensions of one sign to the trace form.
This proves trace signature (4,4), not definiteness over Q.
"""

from sage.matrix.constructor import matrix
from sage.rings.number_field.number_field import CyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice


def test_quartic_cm_full_public_unitary_group():
    """Public generators must generate the proved full group of order 200."""
    field = CyclotomicField(5)
    zeta = field.gen()
    c = zeta + zeta**(-1)
    hermitian = HermitianLattice(
        tuple(zeta**k for k in range(4)),
        matrix(field, [[c, 0], [0, c]]),
    )
    trace = hermitian.trace_lattice()
    assert tuple(map(int, trace.signature_pair())) == (4, 4)
    generators = tuple(hermitian.unitary_group_generators())
    assert generators
    assert all(hermitian.is_order_linear(g) for g in generators)
    assert all(g.domain() is trace and g.codomain() is trace for g in generators)
    assert trace.O().subgroup(generators).cardinality() == 200
