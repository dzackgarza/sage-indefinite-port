"""Exact Hecke rank-one integral/even overlattice oracle.

Independent upstream source: Hecke.jl@e2ab5716,
test/QuadForm/Quad/ZLatticeAutIso.jl:684--686, which asserts exactly two
integral overlattices and one even overlattice for the lattice with Gram [4].

Independent proof, including completeness: L=Z e with b(e,e)=4 has
L^#=(1/4)Z e. An intermediate rank-one lattice M with L subset M subset
L^# is Z(e/d), where d divides 4, because its cyclic subgroup M/L is a
subgroup of L^#/L=Z/4. The norm of its generator is 4/d^2, which is
integral if and only if d=1 or d=2. For d=1 the norm is 4 (even), for
d=2 it is 1 (odd), and for d=4 it is 1/4 (nonintegral). Therefore the
complete original source claim follows, without trusting port output.

The discriminant class 2[e/4] has bilinear self-pairing 1 mod Z (zero)
and quadratic value 1 mod 2Z (nonzero). Consequently bilinear isotropy
allows the index-two odd integral overlattice, while quadratic isotropy
does not. No 'even-only' condition may be used to count ALL integral
overlattices.
"""

import pytest

from tests.acceptance.consumer import lattice


def test_gram_four_has_exactly_two_integral_intermediate_lattices_one_even():
    """Hecke ZLatticeAutIso.jl lines 684--686; proof by all divisors of four."""
    original = lattice([[4]])
    discriminant = original.discriminant_group()
    (generator,) = tuple(discriminant.module_generators())

    admitted = {}
    for multiple in (0, 1, 2, 3):
        try:
            inclusion = original.overlattice(multiple * generator)
        except ValueError:
            assert multiple in (1, 3)
        else:
            assert multiple in (0, 2)
            enlarged = inclusion.codomain()
            (basis_vector,) = tuple(enlarged.module_generators())
            admitted[multiple] = (int(basis_vector.q()), bool(enlarged.is_even()))

    assert admitted == {0: (4, True), 2: (1, False)}
    assert len(admitted) == 2
    assert sum(is_even for _, is_even in admitted.values()) == 1


def test_gram_four_order_two_class_is_bilinearly_but_not_quadratically_isotropic():
    """Hecke lines 684--686: the odd index-two extension has integral b, odd q."""
    original = lattice([[4]])
    bilinear = original.discriminant_module()
    quadratic = original.discriminant_group()
    (b_generator,) = tuple(bilinear.module_generators())
    (q_generator,) = tuple(quadratic.module_generators())
    b_class = 2 * b_generator
    q_class = 2 * q_generator

    assert b_class != bilinear.zero()
    assert bilinear.b(b_class, b_class) == bilinear.bilinear_value_module().zero()
    assert quadratic.q(q_class) != quadratic.quadratic_value_module().zero()
    assert quadratic.q(q_class) == quadratic.quadratic_value_module()(1)
