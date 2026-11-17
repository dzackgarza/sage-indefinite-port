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


def test_fixed_content_orbit_strata_recover_distinct_square_eight_classes():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes, isometry

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[8]])
    category = CosliceCategory(Lattices(ZZ), source)
    primitive = pointed_content_classes(category, plane, 8, 1)
    doubled = pointed_content_classes(category, plane, 8, 2)
    impossible = pointed_content_classes(category, plane, 8, 3)
    assert primitive.cardinality() == cardinal(1)
    assert doubled.cardinality() == cardinal(1)
    assert impossible.cardinality() == cardinal(0)
    assert isometry(category, tuple(primitive)[0], tuple(doubled)[0]) is None
    assert all(int(marking.arrow()(source.module_generators()[0]).content()) == c
               for c, classes in ((1, primitive), (2, doubled)) for marking in classes)


def test_fixed_isotropic_content_strata_do_not_conflate_contents():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    for content in (1, 2, 4):
        classes = pointed_content_classes(category, plane, 0, content)
        assert classes.cardinality() == cardinal(1)
        assert all(int(marking.arrow()(source.module_generators()[0]).content()) == content for marking in classes)


def test_zero_form_pointed_content_strata_use_integral_gl_orbits():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes, isometry

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    source = Lattices(ZZ)([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    for content in (1, 2, 5):
        classes = pointed_content_classes(category, ambient, 0, content)
        assert classes.cardinality() == cardinal(1)
        representative = next(iter(classes))
        assert int(representative.arrow()(source.module_generators()[0]).content()) == content
    first = next(iter(pointed_content_classes(category, ambient, 0, 2)))
    second = next(iter(pointed_content_classes(category, ambient, 0, 5)))
    assert isometry(category, first, second) is None


def test_all_nonzero_square_eight_orbits_have_exact_content_partition():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes, isometry

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[8]])
    category = CosliceCategory(Lattices(ZZ), source)
    classes = pointed_nonzero_norm_classes(category, plane, 8)
    assert classes.cardinality() == cardinal(2)
    first, second = tuple(classes)
    assert {int(marking.arrow()(source.module_generators()[0]).content()) for marking in classes} == {1, 2}
    assert isometry(category, first, second) is None
    assert pointed_nonzero_norm_classes(category, plane, 8).cardinality() == cardinal(2)


def test_unbounded_isotropic_content_cannot_be_materialized_as_finite_orbits():
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    plane = lattice([[0, 1], [1, 0]])
    category = CosliceCategory(Lattices(ZZ), Lattices(ZZ)([[0]]))
    with pytest.raises(ValueError):
        pointed_nonzero_norm_classes(category, plane, 0)


def test_negative_square_content_strata_are_not_lost():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[-8]])
    category = CosliceCategory(Lattices(ZZ), source)
    classes = pointed_nonzero_norm_classes(category, plane, -8)
    assert classes.cardinality() == cardinal(2)
    assert {int(marking.arrow()(source.module_generators()[0]).content()) for marking in classes} == {1, 2}
    assert all(marking.arrow()(source.module_generators()[0]).q() == ZZ(-8) for marking in classes)


def test_lazy_isotropic_orbits_visit_unbounded_content_without_truncation():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import iter_pointed_isotropic_classes

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    iterator = iter_pointed_isotropic_classes(category, plane)
    first = tuple(islice(iterator, 6))
    assert tuple(int(marking.arrow()(source.module_generators()[0]).content()) for marking in first) == (1, 2, 3, 4, 5, 6)
    assert int(next(iterator).arrow()(source.module_generators()[0]).content()) == 7


def test_anisotropic_isotropic_iterator_is_empty():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import iter_pointed_isotropic_classes

    plane = lattice([[1, 0], [0, -2]])
    category = CosliceCategory(Lattices(ZZ), Lattices(ZZ)([[0]]))
    assert tuple(iter_pointed_isotropic_classes(category, plane)) == ()


def test_zero_is_separate_from_all_nonzero_isotropic_orbits():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import iter_pointed_isotropic_classes, isometry

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    zero, first, second = tuple(islice(iter_pointed_isotropic_classes(category, plane, include_zero=True), 3))
    assert zero.arrow()(source.module_generators()[0]) == plane.zero()
    assert first.arrow()(source.module_generators()[0]).content() == ZZ(1)
    assert second.arrow()(source.module_generators()[0]).content() == ZZ(2)
    assert isometry(category, zero, first) is None

    anisotropic = lattice([[1, 0], [0, -2]])
    assert tuple(iter_pointed_isotropic_classes(category, anisotropic, include_zero=True))[0].arrow()(source.module_generators()[0]) == anisotropic.zero()


def test_imprimitive_content_strata_use_the_actual_congruence_subgroup():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes

    plane = lattice([[0, 1], [1, 0]])
    category = CosliceCategory(Lattices(ZZ), Lattices(ZZ)([[0]]))
    gamma = plane.O().congruence_preimage(3)
    classes = pointed_content_classes(category, plane, 0, 2, acting_group=gamma)
    assert classes.cardinality() == cardinal(4)
    assert all(int(marking.arrow()(category.base_object().module_generators()[0]).content()) == 2 for marking in classes)


def test_full_square_eight_strata_split_under_actual_gamma_three():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    plane = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[8]])
    category = CosliceCategory(Lattices(ZZ), source)
    gamma = plane.O().congruence_preimage(3)
    classes = pointed_nonzero_norm_classes(category, plane, 8, acting_group=gamma)
    assert classes.cardinality() == cardinal(6)
    assert {int(marking.arrow()(source.module_generators()[0]).content()) for marking in classes} == {1, 2}


def test_nonisotropic_line_orbits_split_under_gamma_three_without_fixing_sign():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonisotropic_line_classes

    plane = lattice([[0, 1], [1, 0]])
    gamma = plane.O().congruence_preimage(3)
    classes = nonisotropic_line_classes(plane, 8, acting_group=gamma)
    assert classes.cardinality() == cardinal(3)
    assert {int(marked.arrow()(marked.arrow().domain().module_generators()[0]).content()) for marked in classes} == {1, 2}


def test_primitive_isotropic_line_orbits_split_by_actual_so_subgroup():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)([[0]])
    full = tuple(primitive_marking_classes(ARROWS, source, ambient))
    special = tuple(primitive_marking_classes(ARROWS, source, ambient, acting_group=ambient.SO()))
    assert len(full) == 1
    assert len(special) == 2
