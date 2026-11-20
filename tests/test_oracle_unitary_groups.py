"""Exact rank-one positive unitary groups over Gaussian/Eisenstein orders.

Independent proof. An O-linear endomorphism of the free rank-one module O
is multiplication by u in O. For h(x,y)=x*conjugate(y), it is unitary iff
u*conjugate(u)=1. For O=Z[i], writing u=a+bi gives a^2+b^2=1,
whose only integer solutions are (+/-1,0),(0,+/-1): four elements.
For O=Z[omega], omega^2+omega+1=0, the equation is
 a^2-ab+b^2=1. Since (2a-b)^2+3b^2=4, |b|<=1; b=0
forces a=+/-1, while b=+/-1 gives a=0 or a=b: six elements.
Thus the full groups have orders 4 and 6, independently of implementation.
"""

import pytest
from sage.matrix.constructor import matrix
from sage.rings.number_field.number_field import CyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice


@pytest.mark.parametrize("conductor,expected", ((4, 4), (3, 6)))
def test_rank_one_positive_full_unitary_group(conductor, expected):
    field = CyclotomicField(conductor)
    zeta = field.gen()
    lattice = HermitianLattice((field.one(), zeta), matrix(field, [[1]]))
    assert lattice.rank() == 1
    assert lattice.gram_matrix() == matrix(field, [[1]])

    # Public full-group generator API, not a reduction modulo a prime.
    generators = tuple(lattice.definite_unitary_group_generators())
    assert generators
    trace = generators[0].domain()
    assert trace.O().subgroup(generators).cardinality() == expected
    assert all(lattice.is_order_linear(g) for g in generators)

    # Independently enumerated norm-one scalar automorphisms in the same
    # original order; every scalar is tested via the public owned isometry.
    units = ((1, -1, zeta, -zeta) if conductor == 4 else
             (1, -1, zeta, -zeta, zeta**2, -zeta**2))
    scalar_maps = tuple(lattice.unit_isometry(field(u)) for u in units)
    assert len(set(scalar_maps)) == expected
    assert trace.O().subgroup(scalar_maps).cardinality() == expected
