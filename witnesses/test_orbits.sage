r"""Orbit witnesses whose group computations need long runs.

Gamma(3) on 2U: the orbits of null lines are the 16 isotropic lines of the split form on F_3^4.

L=3U+<-6>, maximal null plane in 3U: minus one on <-6> is in the stabilizer and
is -1 on the discriminant group Z/6, whose orthogonal group is {+-1}. The radical
action uses adjacent transpositions, a shear and a sign, which generate GL_3(Z).
"""

from sage_extensions.all import *


def _hyperbolic(count):
    ambient = Lattices(ZZ)("U")
    for _ in range(count - 1):
        ambient = ambient + Lattices(ZZ)("U")
    return ambient


def test_gamma_three_splits_null_lines_of_two_hyperbolic_planes_into_sixteen_orbits():
    r"""Gamma(3)-orbits are the 16 isotropic lines of the split form on F3^4."""
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic
    group = ambient.O()
    gamma = group.congruence_preimage(3)
    e, _f, _u, _v = tuple(ambient.module_generators())
    arrows = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
    marked = arrows.object(ambient.primitive_isotropic_subobject(e).inclusion())
    orbits = tuple(group.split_orbit(marked, gamma))
    assert len(orbits) == 16
    first, second = orbits[0], orbits[1]
    assert arrows.Core().Mor(first, second).an_element_in(gamma) is None
    assert arrows.Core().Mor(first, first).an_element_in(gamma) is not None
    assert all(ambient.q(arrow.arrow()(next(iter(arrow.arrow().domain().module_generators())))) == 0 for arrow in orbits)


def test_rank_three_discriminant_image():
    ambient = _hyperbolic(3) + Lattices(ZZ)("A1").twist(3)
    plane = ambient.primitive_isotropic_subobject(*tuple(ambient.module_generators())[:6:2])
    assert plane.stabilizer_discriminant_image().cardinality() == cardinal(2)
