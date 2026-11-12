r"""Null-plane classes and transporters of 2U under a nontrivial congruence subgroup.

Let L = 2U with hyperbolic basis e1, f1, e2, f2, let rho: O(L) -> O(L/2L) be the
reduction and Gamma the preimage of {+-Id}. Since -Id = Id modulo two, Gamma = Gamma(2),
the kernel of rho; it is nontrivial and of finite index. L/2L with q(x) = x^2/2 mod 2 is
the split quadratic space of dimension four over F_2, and |O^+_4(2)| = 2 q^2 (q^2 - 1)^2
= 72, which is the recorded image order of rho.

O(L) is transitive on primitive null planes (hyperbolic completion of a primitive null
basis in the even unimodular L). The stabilizer P_I of I = <e1, e2> is GL(I) extended by
the skew shears Hom(I^*, I)^skew = Z. Its image is GL(2, 2) extended by the F_2 of skew
shears, of order 6 * 2 = 12, so the O(L)-orbit of the reduction of I has 72/12 = 6
points: all 2 (q + 1) = 6 totally singular planes of L/2L. The Gamma-classes are the
double cosets rho(Gamma) \ F / rho(P_I) with rho(Gamma) = {Id}, so there are exactly 6,
and -Id merges none of them. Two primitive null planes are Gamma-equivalent exactly when
their reductions modulo two agree.
"""
from sage.all import GF, matrix

from sage_indefinite_port.indefinite.marked_lattices import ARROWS
from tests.acceptance.consumer import coordinates, lattice


def _two_u():
    return lattice([[int(i // 2 == j // 2 and i != j) for j in range(4)] for i in range(4)])


def _plus_minus_identity_preimage(ambient):
    group = ambient.O()
    representation = group.congruence_representation(2)
    image = representation.codomain()
    labels = tuple(ambient.module_generating_set())
    negation = ambient.Aut()(dict(zip(labels, tuple(-point for point in ambient.module_generators()), strict=True)))
    target = image.subgroup((representation(negation),))
    assert int(image.cardinality()) == 72
    assert int(target.cardinality()) == 1
    return group.congruence_preimage(2, target)


def _plane(ambient, *vectors):
    return ARROWS.object(ambient.primitive_isotropic_subobject(*vectors).inclusion())


def _reduction(marked):
    arrow = marked.arrow()
    rows = [coordinates(arrow(point)) for point in arrow.domain().module_generators()]
    return tuple(tuple(row) for row in matrix(GF(2), rows).echelon_form().rows())


def test_null_planes_split_into_six_classes_under_the_plus_minus_identity_preimage():
    ambient = _two_u()
    e1, _f1, e2, _f2 = tuple(ambient.module_generators())
    gamma = _plus_minus_identity_preimage(ambient)

    representatives = tuple(ambient.O().split_orbit(_plane(ambient, e1, e2), gamma))

    assert len(representatives) == 6
    assert all(representative in ARROWS for representative in representatives)
    assert all(representative.arrow().codomain() is ambient for representative in representatives)
    assert len({_reduction(representative) for representative in representatives}) == 6
    for position, first in enumerate(representatives):
        for second in representatives[position + 1:]:
            assert ARROWS.Core().Mor(first, second).an_element_in(gamma) is None


def test_gamma_transporters_join_exactly_the_planes_with_equal_reductions():
    r"""J = <e1 + 2 f2, e2 - 2 f1> is primitive and null (b(e1 + 2 f2, e2 - 2 f1) = -2 + 2 = 0) and reduces to I = <e1, e2>; K = <e1, f2> does not."""
    ambient = _two_u()
    e1, f1, e2, f2 = tuple(ambient.module_generators())
    gamma = _plus_minus_identity_preimage(ambient)
    source = _plane(ambient, e1, e2)
    congruent_basis = (e1 + ambient.scalar_multiple(2, f2), e2 - ambient.scalar_multiple(2, f1))
    congruent = _plane(ambient, *congruent_basis)
    assert _reduction(source) == _reduction(congruent)

    witness = ARROWS.Core().Mor(source, congruent).an_element_in(gamma)

    assert witness is not None
    apex = witness.forward().right()
    assert apex in gamma
    target = ambient.primitive_isotropic_subobject(*congruent_basis)
    assert all(target.inclusion().is_in_image(apex(point)) for point in (e1, e2))
    basis = (e1, f1, e2, f2)
    assert all(ambient.b(apex(left), apex(right)) == ambient.b(left, right) for left in basis for right in basis)
    assert all(
        (image - entry) % 2 == 0
        for point in basis
        for image, entry in zip(coordinates(apex(point)), coordinates(point), strict=True)
    )

    other = _plane(ambient, e1, f2)
    assert ARROWS.Core().Mor(source, other).is_empty() is False
    assert ARROWS.Core().Mor(source, other).an_element_in(gamma) is None
