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
