"""The determinant-one subgroup preserves the two null rulings of U."""
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
