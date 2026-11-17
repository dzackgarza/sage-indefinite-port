"""Exact content reduction for nonprimitive pointed lattice isometries (#48)."""
from tests.acceptance.consumer import lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion


def test_imprimitive_isotropic_transport_cancels_content():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    recursion = EichlerRecursion()
    witness = recursion.pointed_isometry(2 * e, 2 * f)
    assert witness is not None
    assert witness(2 * e) == 2 * f
    assert recursion.pointed_isometry(2 * e, 4 * f) is None


def test_imprimitive_nonzero_norm_transport_cancels_content():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    recursion = EichlerRecursion()
    witness = recursion.pointed_isometry(2 * (e + 2 * f), 2 * (2 * e + f))
    assert witness is not None
    assert witness(2 * (e + 2 * f)) == 2 * (2 * e + f)
    assert recursion.pointed_isometry(2 * (e + f), e + 4 * f) is None


def test_imprimitive_isotropic_stabilizer_fixes_original_mark():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    family = EichlerRecursion().pointed_automorphism_generators(3 * e)
    assert all(g(3 * e) == 3 * e and g(e) == e for g in family)
    assert all(plane.b(g(x), g(y)) == plane.b(x, y) for g in family for x in (e, f) for y in (e, f))


def test_imprimitive_nonzero_norm_stabilizer_is_primitive_stabilizer():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    recursion = EichlerRecursion()
    primitive = e + f
    family = tuple(recursion.pointed_automorphism_generators(4 * primitive))
    assert all(g(primitive) == primitive for g in family)
    assert all(g(4 * primitive) == 4 * primitive for g in family)
    assert all(plane.b(g(x), g(y)) == plane.b(x, y) for g in family for x in (e, f) for y in (e, f))


def test_public_imprimitive_pointed_coslice_uses_content_reduction():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    source = Lattices(ZZ)([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    first = category(source.hom([3 * e], codomain=plane))
    second = category(source.hom([3 * f], codomain=plane))
    witness = isometry(category, first, second)
    assert witness is not None and witness(3 * e) == 3 * f
    assert all(g(3 * e) == 3 * e for g in automorphism_group_generators(category, first))
