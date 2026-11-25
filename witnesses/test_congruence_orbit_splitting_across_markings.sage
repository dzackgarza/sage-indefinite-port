r"""Gamma-orbits of every marking of 2U = U + U, split through the finite image modulo three.

The image F of O(2U) in Aut(2U/3) is O^+(4, 3), of order 1152, and O(2U)
is transitive on each marking type below. For Gamma = rho^-1(S) the
Gamma-orbits in an O(2U)-orbit G/P are the double cosets S \ F / rho(P).
The stabilizer images rho(P) are the full stabilizers in F of the reduction,
since each integral parabolic or orthogonal factor maps onto its finite
counterpart:

- a primitive null vector e: O(e^perp/e) and the Eichler transvections give
  O(U/3) (order 4) times F_3^2, order 36, so Gamma(3) has 1152/36 = 32 orbits;
- the line Ze adds the sign: 72 and 16 orbits;
- the plane <e, u>: GL_2(ZZ) and the skew unipotents map onto GL_2(F_3) and
  F_3, order 144, 8 orbits; the flag Ze < <e, u> has the Borel, order 36, 32 orbits;
- the nondegenerate summand U: O(U) x O(U), order 16, 72 orbits.

Gamma(3) acts trivially on 2U/3, so distinct reductions certify distinct orbits.
-id lies in Gamma(2) and not in Gamma(3), so e and -e are equivalent exactly
under the first. For S = rho(P_e), the non-normal stabilizer of the residue of
e, the Gamma-orbits of null vectors are the S-orbits on the 32 isotropic
vectors: +e, -e, the 12 others orthogonal to e, and b(e, w) = 1 or 2 (9 each).

A based pair (e, f) with prescribed Gram U: every embedding U -> 2U has
complement U, so O(2U) is transitive on them, and the pointwise stabilizer is
O(U) on the complement, of order 4 and faithful modulo three: 1152/4 = 288
Gamma(3)-orbits, with distinct residues.

A complete flag of rank three in 3U, under Gamma(2). Reflections in roots
reduce to the transvections of 3U/2, which generate O^+(6, 2), of order 40320.
Its Borel subgroup GL_3-Borel times the alternating radical has order 8 * 8,
and the integral flag stabilizer maps onto it: upper unitriangular GL_3(ZZ)
onto that of GL_3(F_2), skew matrices onto alternating ones. Hence
40320/64 = 630 orbits, with distinct residue flags.

Level and discriminant jointly: on L = U + <6>, A_L = ZZ/6 and O(A_L) = {+-1},
with -id acting by -1. The 3-part ZZ/3 of A_L embeds in L/3L by z |-> 3z, so
an element trivial modulo three acts trivially on A_L: Gamma(3) lies in the
stable subgroup. The joint image F of O(L) in Aut(L/6L), where both conditions
are read, therefore has [F : H] = |im(O(L) -> Aut(L/3L))| for H the image of
Gamma(3) meet the stable subgroup, and not twice that, as a product of the
factor images would give. Its Gamma-orbits are the Gamma(3)-orbits.
"""

from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.finite_rings.finite_field_constructor import GF as SageGF

from sage_extensions.all import *


def _ambient():
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic
    return ambient, ambient.O(), tuple(ambient.module_generators()), tuple(ambient.module_generating_set())


def _residue(labels, point):
    return tuple(int(point.to_vector()(label)) % 3 for label in labels)


def _span(labels, points):
    return sage_matrix(SageGF(3), [_residue(labels, point) for point in points]).row_space()


def _pointed_category(ambient):
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    generator = source.module_generator(next(iter(source.module_generating_set())))
    marking = lambda vector: category.object(source.Mor(ambient)(lambda _label: vector))
    return category, generator, marking


def _arrow_images(marked):
    arrow = marked.arrow()
    return tuple(arrow(point) for point in arrow.domain().module_generators())


def test_based_null_vectors_split_into_thirty_two_gamma_three_orbits():
    ambient, group, (e, _f, _u, _v), labels = _ambient()
    category, generator, marking = _pointed_category(ambient)
    representatives = tuple(group.split_orbit(marking(e), group.congruence_preimage(3)))
    vectors = [category.arrow(representative)(generator) for representative in representatives]
    assert len(representatives) == 32
    assert len({_residue(labels, vector) for vector in vectors}) == 32


def test_a_based_vector_and_its_negative_under_gamma_three_and_gamma_two():
    ambient, group, (e, _f, _u, _v), _labels = _ambient()
    category, _generator, marking = _pointed_category(ambient)
    plus, minus = marking(e), marking(-e)
    isomorphisms = category.Core().Mor(plus, minus)
    assert isomorphisms.an_element_in(group.congruence_preimage(3)) is None
    assert isomorphisms.an_element_in(group.congruence_preimage(2)) is not None


def test_based_vector_stabilizer_in_gamma_three_fixes_the_vector():
    ambient, group, (e, _f, _u, _v), _labels = _ambient()
    gamma = group.congruence_preimage(3)
    generators = tuple(gamma.stabilizer(e).group_generators())
    assert all(g in gamma and g(e) == e for g in generators)


def test_null_vectors_under_the_preimage_of_a_non_normal_stabilizer():
    ambient, group, (e, _f, _u, _v), labels = _ambient()
    representation = group.congruence_representation(3)
    target = representation.codomain().subgroup(tuple(representation(g) for g in group.stabilizer(e).group_generators()))
    assert target.cardinality() == cardinal(36)
    category, generator, marking = _pointed_category(ambient)
    representatives = tuple(group.split_orbit(marking(e), group.congruence_preimage(3, target)))
    vectors = [category.arrow(representative)(generator) for representative in representatives]

    def kind(vector):
        residue = _residue(labels, vector)
        if residue == _residue(labels, e):
            return "e"
        if residue == _residue(labels, -e):
            return "-e"
        return int(ambient.b(e, vector)) % 3

    assert sorted(map(str, map(kind, vectors))) == ["-e", "0", "1", "2", "e"]


def test_null_lines_split_into_sixteen_gamma_three_orbits():
    ambient, group, (e, f, _u, _v), labels = _ambient()
    arrows = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
    gamma = group.congruence_preimage(3)
    marked = arrows.object(ambient.primitive_isotropic_subobject(e).inclusion())
    representatives = tuple(group.split_orbit(marked, gamma))
    assert len(representatives) == 16
    assert len({_span(labels, _arrow_images(r)) for r in representatives}) == 16
    other = arrows.object(ambient.primitive_isotropic_subobject(f).inclusion())
    assert arrows.Core().Mor(marked, other).an_element_in(gamma) is None


def test_null_planes_split_into_eight_gamma_three_orbits_with_stabilizers():
    ambient, group, (e, f, u, v), labels = _ambient()
    arrows = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
    gamma = group.congruence_preimage(3)
    plane = ambient.primitive_isotropic_subobject(e, u)
    marked = arrows.object(plane.inclusion())
    representatives = tuple(group.split_orbit(marked, gamma))
    assert len(representatives) == 8
    assert len({_span(labels, _arrow_images(r)) for r in representatives}) == 8
    opposite = arrows.object(ambient.primitive_isotropic_subobject(f, v).inclusion())
    assert arrows.Core().Mor(marked, opposite).an_element_in(gamma) is None
    generators = tuple(gamma.stabilizer(plane, action="setwise").group_generators())
    assert all(g in gamma and plane.stabilizes(g) for g in generators)


def test_null_flags_split_into_thirty_two_gamma_three_orbits():
    ambient, group, (e, _f, u, _v), _labels = _ambient()
    flags = Cat().Mor(FiniteOrdinalCategory(3), Lattices(ZZ))
    line = ambient.subobject_on([e])
    plane = ambient.subobject_on([e, u])
    line_into_plane = line.Emb(plane)(line.inclusion().factor_through(plane.inclusion()))
    marked = flags((line_into_plane, plane.inclusion()))
    representatives = tuple(group.split_orbit(marked, group.congruence_preimage(3)))
    assert len(representatives) == 32
    assert all(representative in flags for representative in representatives)


def test_hyperbolic_summands_split_into_seventy_two_gamma_three_orbits():
    ambient, group, (e, f, _u, _v), labels = _ambient()
    arrows = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
    marked = arrows.object(ambient.subobject_on([e, f]).inclusion())
    representatives = tuple(group.split_orbit(marked, group.congruence_preimage(3)))
    assert len(representatives) == 72
    assert len({_span(labels, _arrow_images(r)) for r in representatives}) == 72


def _residue_modulo(labels, point, modulus):
    return tuple(int(point.to_vector()(label)) % modulus for label in labels)


def test_based_pairs_of_prescribed_gram_split_into_288_gamma_three_orbits():
    ambient, group, (e, f, _u, _v), labels = _ambient()
    source = Lattices(ZZ)("U")
    category = Lattices(ZZ).CosliceUnder(source)
    names = tuple(source.module_generating_set())
    images = dict(zip(names, (e, f), strict=True))
    marked = category.object(source.Mor(ambient)(lambda label: images[label]))
    representatives = tuple(group.split_orbit(marked, group.congruence_preimage(int(3))))
    assert len(representatives) == int(288)
    pairs = set([
        tuple(_residue(labels, category.arrow(r)(source.module_generator(name))) for name in names)
        for r in representatives
    ])
    assert len(pairs) == int(288)


def test_complete_flags_of_rank_three_split_into_630_gamma_two_orbits():
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic
    group = ambient.O()
    e1, _f1, e2, _f2, e3, _f3 = tuple(ambient.module_generators())
    labels = tuple(ambient.module_generating_set())
    flags = Cat().Mor(FiniteOrdinalCategory(int(4)), Lattices(ZZ))
    line = ambient.subobject_on([e1])
    plane = ambient.subobject_on([e1, e2])
    space = ambient.subobject_on([e1, e2, e3])
    line_into_plane = line.Emb(plane)(line.inclusion().factor_through(plane.inclusion()))
    plane_into_space = plane.Emb(space)(plane.inclusion().factor_through(space.inclusion()))
    marked = flags((line_into_plane, plane_into_space, space.inclusion()))
    representatives = tuple(group.split_orbit(marked, group.congruence_preimage(int(2))))
    assert len(representatives) == int(630)
    assert all(representative in flags for representative in representatives)

    def residue_flag(representative):
        diagram = representative.functor()
        shape = diagram.domain()
        objects = tuple(shape.objects())
        spans = []
        for source in objects[:-int(1)]:
            arrow = diagram(shape.Mor(source, objects[-int(1)]).unique())
            points = [_residue_modulo(labels, arrow(p), int(2)) for p in arrow.domain().module_generators()]
            spans.append(sage_matrix(SageGF(int(2)), points).row_space())
        assert spans[int(0)].dimension() == int(1) and spans[int(1)].dimension() == int(2) and spans[int(2)].dimension() == int(3)
        return tuple(spans)

    assert len(set([residue_flag(r) for r in representatives])) == int(630)


def test_joint_level_and_discriminant_image_on_u_plus_six():
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[int(6)]])
    group = ambient.O()
    e, f, w = tuple(ambient.module_generators())
    labels = tuple(ambient.module_generating_set())
    stable = group.stable_subgroup()
    level_three = group.congruence_preimage(int(3))
    negation = ambient.Aut()(dict(zip(labels, (-e, -f, -w), strict=True)))
    assert negation not in stable and negation not in level_three
    joint = group.congruence_representation(int(6))
    image = joint.codomain()
    section = joint.image_section()
    elements = list(image)
    members = tuple(x for x in elements if section(x) in level_three and section(x) in stable)
    three = list(group.congruence_representation(int(3)).codomain())
    assert len(elements) == len(members) * len(three)
    gamma = group.congruence_preimage(int(6), image.subgroup(members))
    category, _generator, marking = _pointed_category(ambient)
    assert len(tuple(group.split_orbit(marking(e), gamma))) == len(tuple(group.split_orbit(marking(e), level_three)))
    arrows = Cat().Mor(FiniteOrdinalCategory(int(2)), Lattices(ZZ))
    line = arrows.object(ambient.primitive_isotropic_subobject(e).inclusion())
    assert len(tuple(group.split_orbit(line, gamma))) == len(tuple(group.split_orbit(line, level_three)))
