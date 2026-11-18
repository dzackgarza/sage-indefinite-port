"""Exact finite rational/integral coset data for isotropic markings."""

import pytest
from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isotropic_stabilizer_right_cosets
from tests.acceptance.consumer import lattice


@pytest.mark.parametrize("gram", ([], [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]), ids=("empty_ambient", "two_hyperbolic_planes"))
def test_zero_isotropic_sublattice_has_one_right_coset(gram):
    r"""At I=0, its rational stabilizer is O(L) and the integral
    stabilizer is the same O(L); their right-coset space is a singleton.
    """
    ambient = lattice(gram) if gram else Lattices(ZZ)(0)
    zero = ambient.primitive_isotropic_subobject()
    classes = tuple(isotropic_stabilizer_right_cosets(ARROWS, ARROWS.object(zero.inclusion())))
    assert len(classes) == 1
    (witness,) = classes
    assert witness.domain() is ambient and witness.codomain() is ambient
    assert all(witness(point) == point for point in ambient.module_generators())


def test_glued_line_right_cosets_are_distinct_integral_lattice_images():
    r"""A glued line has a nontrivial finite commensurability orbit.

    The rational line group acts on four integral presentations of the
    invariant perpendicular lattice. Distinct right coset representatives
    differ by a map that does not preserve the selected lattice.
    """
    from itertools import combinations
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[0, 2, 0, 0], [2, 0, 1, 0], [0, 1, 0, 2], [0, 0, 2, 0]])
    point = ambient.module_generators()[0]
    marked = ARROWS.object(ambient.primitive_isotropic_subobject(point).inclusion())
    action = RECURSION.isotropic_stabilizer_action((point,), ((-point,),))
    representatives = tuple(isotropic_stabilizer_right_cosets(ARROWS, marked))
    assert len(representatives) == 4
    inclusion = action.lattice_inclusion()
    assert all(representative.domain() is representative.codomain() for representative in representatives)
    assert any(action._carries_lattice(representative, inclusion, inclusion)
               for representative in representatives)
    assert all(
        not action._carries_lattice((~first) * second, inclusion, inclusion)
        for first, second in combinations(representatives, 2)
    )


def test_split_maximal_null_flag_has_one_right_coset():
    r"""The integral block Borel of a maximal null plane in 2U extends
    integrally to the dual hyperbolic summand; its rational extension
    group therefore coincides with the integral flag stabilizer.
    """
    from sage_extensions.all import Cat
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import _flag_marking_from_basis

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(3), Lattices(ZZ))
    marked = category(_flag_marking_from_basis(ambient, (e, u), (1, 2)))
    cosets = tuple(isotropic_stabilizer_right_cosets(category, marked))
    assert len(cosets) == 1
    assert cosets[0].domain() is cosets[0].codomain()
    assert all(ambient.b(cosets[0](point), cosets[0](other)) == ambient.b(point, other)
               for point in (e, f, u, v) for other in (e, f, u, v))
