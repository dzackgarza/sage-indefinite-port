r"""Voight's fundamental domain for the maximal order of the quaternion algebra of discriminant 6 over `Q` (#54).

`B = (-1, 3)_Q` ramifies at `3` (`-1` is not a square modulo 3) and at `2` (`(-1, 3)_2 = (-1)^{1 \cdot 1} = -1`), and splits at infinity
(`3 > 0`): it is the indefinite algebra of discriminant `6`. For a maximal order `O`, Shimizu's formula gives
`\mathrm{area}(O^1 \backslash H) = (8\pi / 4\pi^2)\, \zeta(2)\, (2 - 1)(3 - 1) = 2\pi/3`. The numbers of elliptic points of order 2 and 3 are
`\prod_{p \mid 6} (1 - (-4/p)) = 1 \cdot 2 = 2` and `\prod_{p \mid 6} (1 - (-3/p)) = 2 \cdot 1 = 2` (Voight, *Quaternion Algebras*, Sec. 39.4),
so the signature is `(0; 2, 2, 3, 3)`: `2\pi (2 \cdot 0 - 2 + 1/2 + 1/2 + 2/3 + 2/3) = 2\pi/3`."""

from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.fuchsian_domain import FuchsianDomain


def test_fuchsian_domain_of_the_discriminant_six_maximal_order():
    algebra = QuaternionAlgebra(SageQQ, int(-1), int(3))
    domain = FuchsianDomain(algebra, algebra.maximal_order().basis())
    assert domain.is_maximal()
    assert domain.covolume_over_2pi() == SageQQ((int(1), int(3)))
    fundamental = domain.fundamental_domain()
    assert fundamental.area_over_2pi == SageQQ((int(1), int(3)))
    assert fundamental.signature() == (int(0), (int(2), int(2), int(3), int(3)))
    generators = domain.sl1_generators()
    assert all(domain.in_order(g) and g.reduced_norm() == int(1) for g in generators)
    one = algebra(int(1))
    for g, order in fundamental.elliptic:
        power = one
        for exponent in range(int(1), order + int(1)):
            power = power * g
            assert (power == one or power == -one) == (exponent == order)
