"""Finite-index subgroup actions retain the supplied marking category."""
import pytest
from sage_extensions.all import Lattices, ZZ
from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
from tests.acceptance.consumer import lattice


def test_fixed_null_vector_transporter_respects_special_orthogonal_subgroup():
    ambient = lattice([[0, 1], [1, 0]])
    e, f = tuple(ambient.module_generators())
    source = lattice([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    first = category.object(source.hom([e], codomain=ambient))
    second = category.object(source.hom([f], codomain=ambient))
    hom = category.Core().Mor(first, second)
    witness = hom.an_element().forward().right()
    assert witness.domain() is ambient and witness.codomain() is ambient
    assert witness(e) == f
    assert witness.determinant() == -1
    # O(U)={I,-I,swap,-swap}; SO(U)={I,-I} cannot exchange its null lines.
    special = ambient.SO()
    assert hom.an_element_in(special) is None
    fixed = category.Core().Mor(first, first).an_element_in(special)
    assert fixed is not None
    apex = fixed.forward().right()
    assert apex(e) == e and apex.determinant() == 1


def test_composite_level_six_preimage_uses_actual_integral_image():
    """O(U) has four elements; -I is in Gamma(2) but not Gamma(3) or Gamma(6)."""
    ambient = lattice([[0, 1], [1, 0]])
    group = ambient.O()
    minus = ambient.Aut()({label: -ambient.module_generator(label)
                           for label in ambient.module_generating_set()})
    assert minus in group.congruence_preimage(2)
    assert minus not in group.congruence_preimage(3)
    image = group.congruence_representation(6)
    kernel = group.congruence_preimage(6)
    assert int(image.codomain().cardinality()) == 4
    assert minus not in kernel
    identity_image = image(group.identity())
    assert all(image(generator) == identity_image for generator in kernel.group_generators())
    assert kernel.index() == 4


@pytest.mark.parametrize("modulus, expected_order", ((2, 72), (3, 1152), (4, 4608)))
def test_two_hyperbolic_planes_congruence_representation_has_exact_image(modulus, expected_order):
    r"""The split rank-four orthogonal action has the two SL_2 factors
    and their factor-exchange involution. Over Z/m the corresponding
    orders are 2*6^2, 2*24^2, and 2*48^2 for m=2,3,4. Check the
    *actual generated image* instead of assuming independent characters.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    group = ambient.O()
    representation = group.congruence_representation(modulus)
    assert int(representation.codomain().cardinality()) == expected_order
    basis = tuple(ambient.module_generators())
    for generator in group.group_generators():
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for x in basis for y in basis)
        assert representation(generator) in representation.codomain()


@pytest.mark.parametrize("kind, count", (("based", 4), ("moving", 2), ("flag", 2)))
def test_joint_determinant_and_discriminant_conditions_use_the_actual_image(kind, count):
    r"""O(U(3)) consists of +/-I and +/-swap, all distinct modulo three.

    The joint determinant/discriminant image has order four, although
    its two separate images have orders two and four. Its kernel is
    trivial. That kernel has four orbits on the primitive null vectors
    in the orbit of e, and two on the corresponding moving null lines.
    """
    from sage_extensions.all import Cat
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS

    ambient = lattice([[0, 3], [3, 0]])
    e, _f = tuple(ambient.module_generators())
    source = lattice([[0]])
    arrow = source.hom([e], codomain=ambient)
    if kind == "based":
        category = CosliceCategory(Lattices(ZZ), source)
        marked = category(arrow)
    elif kind == "moving":
        category = ARROWS
        marked = category(arrow)
    else:
        category = Cat().Mor(FiniteOrdinalCategory(2), Lattices(ZZ))
        # Fun([1], Lattices) is the arrow category itself.
        marked = category(arrow)
    group = ambient.O()
    stable = ambient.stable_orthogonal_group()
    subgroup = group.intersection(ambient.SO(), stable)
    quotient = subgroup.finite_character_quotient()
    assert int(quotient.image_keys().cardinality()) == 4
    assert int(quotient.subgroup_image_keys().cardinality()) == 1
    assert subgroup.index() == 4
    representatives = tuple(group.split_orbit(marked, subgroup))
    assert len(representatives) == count
    assert all(representative in category for representative in representatives)
    for position, first in enumerate(representatives):
        for second in representatives[position + 1:]:
            assert category.Core().Mor(first, second).an_element_in(subgroup) is None
