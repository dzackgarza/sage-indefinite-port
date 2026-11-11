r"""Principality of right ideals in the definite quaternion algebra of discriminant 23 (#54).

`B = (-1, -23)_Q` ramifies at `23` (`-1` is not a square modulo `23`) and at infinity, and splits at `2` (`-23 \equiv 1 \bmod 8` is a
`2`-adic square). `O = Z\langle 1, i, (1 + j)/2, (i + k)/2 \rangle` is an order: `((1 + j)/2)^2 = -6 + (1 + j)/2`. For
`x = a + b i + c (1 + j)/2 + d (i + k)/2`, `\mathrm{nrd}(x) = (a + c/2)^2 + (b + d/2)^2 + 23 (c^2 + d^2)/4`, so `\mathrm{nrd}(x) \le 2` forces
`c = d = 0`: the units are `\pm 1, \pm i` and the elements of norm `2` are `\pm 1 \pm i = (1 + i) u`. Hence `(1 + i) O` is the only principal
right ideal of norm `2`. `y = (1 + j)/2` has norm `6`, and `(1 + i)^{-1} y = (1 - i + j - k)/4 \notin O` (its `j`-coefficient is `1/4`), so
`I = 2 O + y O`, of reduced norm `2` (it contains `2 O` and an element of norm `6`, and `y \notin 2 O`), is not principal."""

from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.quaternion_voronoi_number_field import QuaternionOrderVoronoi


def test_definite_principal_and_nonprincipal_right_ideals_of_discriminant_23():
    algebra = QuaternionAlgebra(SageQQ, int(-1), int(-23))
    _one, i, j, k = algebra.basis()
    one = algebra(int(1))
    order = [one, i, (one + j) / int(2), (i + k) / int(2)]
    voronoi = QuaternionOrderVoronoi(algebra, order)
    assert voronoi.is_totally_definite()
    principal = [(one + i) * o for o in order]
    x = voronoi.norm_generator(principal, SageZZ.ideal(int(2)))
    assert x is not None and x.reduced_norm() == int(2)
    y = (one + j) / int(2)
    ideal = [int(2) * o for o in order] + [y * o for o in order]
    assert voronoi.norm_generator(ideal, SageZZ.ideal(int(2))) is None
