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


def test_gamma_two_splits_null_lines_of_two_hyperbolic_planes_into_nine_orbits():
    r"""O(2U) maps onto O+(4, 2) modulo two, and -1 is trivial there, so the Gamma(2)-orbits of primitive
    null lines are the (2 + 1)^2 = 9 singular points of the hyperbolic quadric over F_2."""
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic
    group = ambient.O()
    gamma = group.congruence_preimage(2)
    e, _f, _u, _v = tuple(ambient.module_generators())
    arrows = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
    marked = arrows.object(ambient.primitive_isotropic_subobject(e).inclusion())
    orbits = tuple(group.split_orbit(marked, gamma))
    assert len(orbits) == 9
    assert arrows.Core().Mor(orbits[0], orbits[1]).an_element_in(gamma) is None


def test_gamma_three_transporters_between_null_lines_of_two_hyperbolic_planes():
    r"""Gamma(3) acts trivially on L/3L, so it preserves the line of a null vector modulo 3.

    Ze and Z(e + 3u) have the same reduction, and the Eichler transvection
    x -> x + b(x, e) 3u - b(x, 3u) e lies in Gamma(3) and moves e to e + 3u up to sign.
    Ze and Zf have different reductions, so no element of Gamma(3) moves one to the other.
    """
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic
    group = ambient.O()
    gamma = group.congruence_preimage(3)
    e, f, u, _v = tuple(ambient.module_generators())
    arrows = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
    line = lambda vector: arrows.object(ambient.primitive_isotropic_subobject(vector).inclusion())
    near = arrows.Core().Mor(line(e), line(e + u + u + u)).an_element_in(gamma)
    assert near is not None
    assert arrows.Core().Mor(line(e), line(f)).an_element_in(gamma) is None


def test_rank_three_null_stabilizer_of_twisted_hyperbolic_sum_fixes_one_isotropic_line():
    r"""L = U + U + U(3), I = <e1, e2, e3>. A_L = (Z/3)^2 with the hyperbolic form.

    Its orthogonal group has order four: +-1 fix both isotropic lines, and +-swap exchange
    them. Stab(I) fixes the isotropic line (I_Q meet L^dual)/I, so its image lies in {+-1},
    and -1 on L lies in Stab(I) and acts by -1, which is not the identity on (Z/3)^2.
    """
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic.twist(3)
    basis = tuple(ambient.module_generators())
    labels = tuple(ambient.module_generating_set())
    space = ambient.primitive_isotropic_subobject(*basis[::2])
    minus = ambient.Aut()(dict(zip(labels, tuple(-b for b in basis), strict=True)))
    assert space.stabilizes(minus)
    rho = ambient.discriminant_representation()
    image = space.stabilizer_discriminant_image()
    assert image.cardinality() == cardinal(2)
    assert rho(minus) in image


def test_gamma_three_line_stabilizer_has_schreier_generators_of_index_72():
    r"""P = Stab(Ze) in O(2U) maps onto the stabilizer of the point <e> among the 16 isotropic
    points of F_3^4: an element fixing <e> modulo 3 moves Ze inside its Gamma(3)-orbit, and a
    Gamma(3) element brings it back. The image of O(2U) has order 1152 and is transitive on the
    16 points, so [P : Gamma(3) meet P] = 1152 / 16 = 72.
    """
    from sage_extensions.categories.lattice_centralizers import finite_preimage_generators

    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic
    group = ambient.O()
    gamma = group.congruence_preimage(3)
    e, _f, _u, _v = tuple(ambient.module_generators())
    line = ambient.primitive_isotropic_subobject(e)
    stabilizer = gamma.stabilizer(line, action="setwise")
    generators = tuple(stabilizer.group_generators())
    assert all(g in gamma and line.stabilizes(g) for g in generators)
    ambient_generators = tuple(group.stabilizer(line, action="setwise").group_generators())
    _schreier, index = finite_preimage_generators(ambient_generators, lambda g: g, lambda g: g in gamma)
    assert index == 72
