"""Saturated rank-one isotropic stabilizers for nonprimitive sublattices."""
from tests.acceptance.consumer import lattice


def test_nonprimitive_line_stabilizer_equals_saturated_line_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    source = Lattices(ZZ)([[0]])
    primitive = ARROWS(source.hom([e], codomain=plane))
    imprimitive = ARROWS(source.hom([3 * e], codomain=plane))
    first = tuple(automorphism_group_generators(ARROWS, primitive))
    second = tuple(automorphism_group_generators(ARROWS, imprimitive))
    assert all(g(e) in (e, -e) for g in second)
    assert all(g(3 * e) in (3 * e, -3 * e) for g in second)
    assert all(plane.b(g(x), g(y)) == plane.b(x, y) for g in second for x in (e, f) for y in (e, f))
    assert {tuple(g(x) for x in (e, f)) for g in first} == {tuple(g(x) for x in (e, f)) for g in second}


def test_nonprimitive_line_in_two_u_has_full_setwise_shear_action():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0]])
    marked = ARROWS(source.hom([2 * e], codomain=ambient))
    family = tuple(automorphism_group_generators(ARROWS, marked))
    assert all(g(e) in (e, -e) for g in family)
    assert all(g(2 * e) in (2 * e, -2 * e) for g in family)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y) for g in family for x in (e, f, u, v) for y in (e, f, u, v))
    assert any(g(u) != u or g(v) != v for g in family)


def test_nonprimitive_moving_line_transporter_cancels_common_content():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    source = Lattices(ZZ)([[0]])
    first = ARROWS(source.hom([3 * e], codomain=plane))
    second = ARROWS(source.hom([3 * f], codomain=plane))
    third = ARROWS(source.hom([6 * f], codomain=plane))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert witness(e) in (f, -f)
    assert isometry(ARROWS, first, third) is None


def test_moving_line_orbit_classes_retain_saturation_index():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_line_classes, ARROWS, isometry

    plane = lattice([[0, 1], [1, 0]])
    first = nonprimitive_isotropic_line_classes(plane, 2)
    second = nonprimitive_isotropic_line_classes(plane, 4)
    assert first.cardinality() == cardinal(1)
    assert second.cardinality() == cardinal(1)
    a, b = tuple(first)[0], tuple(second)[0]
    assert a.arrow().saturation().inclusion().is_primitive()
    assert int(a.arrow()(a.arrow().domain().module_generators()[0]).content()) == 2
    assert int(b.arrow()(b.arrow().domain().module_generators()[0]).content()) == 4
    assert isometry(ARROWS, a, b) is None


def test_nonprimitive_moving_lines_split_under_actual_level_three_kernel():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_line_classes, ARROWS, isometry

    ambient = lattice([[0, 1], [1, 0]])
    congruence = ambient.O().congruence_preimage(3)
    classes = tuple(nonprimitive_isotropic_line_classes(ambient, 2, acting_group=congruence))
    assert len(classes) == 2
    assert all(int(marked.arrow()(next(iter(marked.arrow().domain().module_generators()))).content()) == 2
               for marked in classes)
    assert isometry(ARROWS, classes[0], classes[1]) is not None
    assert len(tuple(ambient.O().split_orbit(classes[0], congruence))) == 2


def test_unrestricted_null_line_stream_dovetails_all_contents_and_subgroup_orbits():
    from itertools import islice
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_line_classes

    ambient = lattice([[0, 1], [1, 0]])
    full = tuple(islice(iter_all_isotropic_line_classes(ambient), 5))
    assert [int(marked.arrow()(next(iter(marked.arrow().domain().module_generators()))).content())
            for marked in full] == [1, 2, 3, 4, 5]

    subgroup = ambient.O().congruence_preimage(3)
    split = tuple(islice(iter_all_isotropic_line_classes(ambient, acting_group=subgroup), 9))
    assert [int(marked.arrow()(next(iter(marked.arrow().domain().module_generators()))).content())
            for marked in split] == [1, 1, 2, 2, 3, 3, 4, 4, 5]
    assert len({(int(marked.arrow()(next(iter(marked.arrow().domain().module_generators()))).content()),
                 tuple(int(v) for v in marked.arrow()(next(iter(marked.arrow().domain().module_generators()))).to_vector()))
                for marked in split}) == 9


def test_unrestricted_null_line_stream_terminates_for_anisotropic_lattice():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_line_classes

    ambient = Lattices(ZZ)([[1, 0], [0, 1]])
    assert tuple(iter_all_isotropic_line_classes(ambient)) == ()


def test_unrestricted_null_line_stream_covers_degenerate_radical_contents():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_line_classes

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    markings = tuple(islice(iter_all_isotropic_line_classes(ambient), 4))
    assert [int(marked.arrow()(next(iter(marked.arrow().domain().module_generators()))).content())
            for marked in markings] == [1, 2, 3, 4]


def test_unrestricted_rank_two_null_sublattices_enumerate_distinct_smith_types():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_sublattice_classes, isometry, ARROWS

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    classes = tuple(islice(iter_all_isotropic_sublattice_classes(ambient, 2), 3))
    assert len(classes) == 3
    indices = []
    for marked in classes:
        arrow = marked.arrow()
        saturation = arrow.saturation()
        embedded = saturation.inclusion().domain().subobject_on(tuple(
            saturation.inclusion().lift(arrow(x)) for x in arrow.domain().module_generators()
        )).inclusion()
        indices.append(int(embedded.cokernel().cardinality()))
    assert indices == [1, 2, 4]
    assert all(isometry(ARROWS, classes[i], classes[j]) is None
               for i in range(3) for j in range(i + 1, 3))
    assert tuple(iter_all_isotropic_sublattice_classes(ambient, 3)) == ()


def test_unrestricted_rank_three_null_sublattices_preserve_smith_index():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_sublattice_classes

    plane = Lattices(ZZ)("U")
    ambient = plane + plane + plane
    classes = tuple(islice(iter_all_isotropic_sublattice_classes(ambient, 3), 2))
    assert len(classes) == 2
    indices = []
    for marked in classes:
        arrow = marked.arrow()
        saturation = arrow.saturation()
        embedded = saturation.inclusion().domain().subobject_on(tuple(
            saturation.inclusion().lift(arrow(x)) for x in arrow.domain().module_generators()
        )).inclusion()
        indices.append(int(embedded.cokernel().cardinality()))
    assert indices == [1, 2]


def test_unrestricted_isotropic_sublattice_iterator_terminates_above_ambient_rank():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_sublattice_classes

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    assert tuple(iter_all_isotropic_sublattice_classes(ambient, 3)) == ()


def test_unrestricted_isotropic_sublattice_iterator_terminates_above_degenerate_witt_bound():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_sublattice_classes

    ambient = Lattices(ZZ)([[0, 0, 0], [0, 1, 0], [0, 0, 1]])
    # The one-dimensional radical is the entire isotropic locus over QQ.
    assert tuple(iter_all_isotropic_sublattice_classes(ambient, 2)) == ()


def test_zero_form_rank_two_isotropic_sublattices_have_complete_smith_strata():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        iter_all_isotropic_sublattice_classes, nonprimitive_isotropic_sublattice_classes,
    )

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    classes = tuple(islice(iter_all_isotropic_sublattice_classes(ambient, 2), 3))
    indices = []
    for marked in classes:
        arrow = marked.arrow()
        saturated = arrow.saturation()
        embedded = saturated.inclusion().domain().subobject_on(tuple(
            saturated.inclusion().lift(arrow(x)) for x in arrow.domain().module_generators()
        )).inclusion()
        indices.append(int(embedded.cokernel().cardinality()))
    assert indices == [1, 2, 4]
    assert len(tuple(nonprimitive_isotropic_sublattice_classes(ambient, 2, (1, 2), acting_group=ambient.O()))) == 1


def test_zero_form_primitive_moving_sublattices_have_one_class_at_each_attainable_rank():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    for rank, expected in ((1, 1), (2, 1), (3, 0)):
        source = Lattices(ZZ)([[0] * rank for _ in range(rank)])
        classes = tuple(primitive_marking_classes(ARROWS, source, ambient))
        assert len(classes) == expected
        assert all(marked.arrow().is_primitive() for marked in classes)


def test_zero_form_index_two_sublattices_split_into_three_level_two_orbits():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    full = tuple(nonprimitive_isotropic_sublattice_classes(ambient, 2, (1, 2)))
    kernel = ambient.O().congruence_preimage(2)
    split = tuple(nonprimitive_isotropic_sublattice_classes(ambient, 2, (1, 2), acting_group=kernel))
    assert len(full) == 1 and len(split) == 3
    # Index-two submodules correspond to the three lines in F_2^2.
    # Gamma(2) fixes every residue class, hence fixes each of the three.
    residue_lines = set()
    for marked in split:
        arrow = marked.arrow()
        residue_lines.add(frozenset(tuple(int(x.to_vector()(label)) % 2
                                          for label in ambient.module_generating_set())
                                    for x in (arrow(point) for point in arrow.domain().module_generators())))
    assert len(residue_lines) == 3


def test_zero_form_moving_finite_index_submodule_transporter_decides_smith_type():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    first = ARROWS.object(source.hom([e, 2 * f], codomain=ambient))
    second = ARROWS.object(source.hom([2 * e, f], codomain=ambient))
    third = ARROWS.object(source.hom([2 * e, 2 * f], codomain=ambient))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    target_subobject = ambient.subobject_on((2 * e, f))
    assert all(target_subobject.inclusion().is_in_image(witness(x)) for x in (e, 2 * f))
    assert isometry(ARROWS, first, third) is None


def test_zero_form_primitive_moving_line_and_plane_stabilizers_are_full_gl_parabolics():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    e, f = tuple(ambient.module_generators())
    line_source = Lattices(ZZ)([[0]])
    plane_source = Lattices(ZZ)([[0, 0], [0, 0]])
    line = ARROWS.object(line_source.hom([e], codomain=ambient))
    plane = ARROWS.object(plane_source.hom([e, f], codomain=ambient))
    line_generators = tuple(automorphism_group_generators(ARROWS, line))
    plane_generators = tuple(automorphism_group_generators(ARROWS, plane))
    assert line_generators and plane_generators
    assert all(g in ambient.O() for g in (*line_generators, *plane_generators))
    assert all(g(e) in (e, -e) for g in line_generators)
    assert any(g(e) == f for g in plane_generators)


def test_zero_form_full_orthogonal_generators_are_owned_isometries():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    e, f = tuple(ambient.module_generators())
    generators = tuple(automorphism_group_generators(Lattices(ZZ), ambient))
    assert generators and all(g in ambient.O() for g in generators)
    assert any(g(e) == f for g in generators)


def test_zero_form_based_primitive_frames_have_one_class_at_each_attainable_rank():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    for rank, expected in ((1, 1), (2, 1), (3, 0)):
        source = Lattices(ZZ)([[0] * rank for _ in range(rank)])
        category = Lattices(ZZ).CosliceUnder(source)
        classes = tuple(primitive_marking_classes(category, source, ambient))
        assert len(classes) == expected
        assert all(marked.arrow().is_primitive() for marked in classes)


def test_zero_form_based_primitive_frame_transporter_preserves_each_marked_generator():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category.object(source.hom([e, f], codomain=ambient))
    second = category.object(source.hom([f, -e], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(e) == f and witness(f) == -e
    assert all((~witness)(witness(x)) == x for x in (e, f))


def test_zero_form_based_primitive_frame_stabilizer_fixes_source_pointwise():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = Lattices(ZZ)([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, f, u = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category.object(source.hom([e, f], codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators
    assert all(g(e) == e and g(f) == f for g in generators)
    assert any(g(u) != u for g in generators)


def test_zero_form_full_based_frame_has_trivial_pointwise_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category.object(source.hom(tuple(ambient.module_generators()), codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators == (ambient.O().identity(),)


def test_zero_form_changed_basis_pointwise_stabilizer_conjugates_correctly():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = Lattices(ZZ)([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, f, u = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category.object(source.hom([e + u, f], codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators
    assert all(g(e + u) == e + u and g(f) == f for g in generators)
    assert any(g(u) != u for g in generators)


def test_unrestricted_rank_two_null_sublattices_split_under_special_orthogonal_group():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_sublattice_classes, ARROWS, isometry

    plane = Lattices(ZZ)("U")
    ambient = plane + plane
    classes = tuple(islice(iter_all_isotropic_sublattice_classes(ambient, 2, acting_group=ambient.SO()), 2))
    assert len(classes) == 2
    assert all(marked.arrow().saturation().inclusion().domain().rank() == 2 for marked in classes)
    assert isometry(ARROWS, classes[0], classes[1]) is not None
    assert ARROWS.Core().Mor(classes[0], classes[1]).an_element_in(ambient.SO()) is None


def test_canonical_moving_null_line_is_sign_invariant_and_retains_content():
    import pytest
    from sage_indefinite_port.indefinite.marked_lattices import canonical_moving_isotropic_line

    ambient = lattice([[0, 1], [1, 0]])
    e, f = tuple(ambient.module_generators())
    outputs = [canonical_moving_isotropic_line(point) for point in (3 * e, -3 * e, 3 * f)]
    canonical = [marked.arrow()(next(iter(marked.arrow().domain().module_generators())))
                 for marked, _witness in outputs]
    assert canonical[0] == canonical[1] == canonical[2]
    assert int(canonical[0].content()) == 3
    for (marked, witness), target in zip(outputs, (3 * e, -3 * e, 3 * f), strict=True):
        generator = marked.arrow()(next(iter(marked.arrow().domain().module_generators())))
        assert witness(generator) in (target, -target)
    with pytest.raises(ValueError, match="nonzero null"):
        canonical_moving_isotropic_line(ambient.zero())


def test_canonical_moving_null_line_on_degenerate_radical_has_exact_witness():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import canonical_moving_isotropic_line

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    e, f = tuple(ambient.module_generators())
    first, first_witness = canonical_moving_isotropic_line(2 * e)
    second, second_witness = canonical_moving_isotropic_line(2 * f)
    x = first.arrow()(next(iter(first.arrow().domain().module_generators())))
    y = second.arrow()(next(iter(second.arrow().domain().module_generators())))
    assert x == y and int(x.content()) == 2
    assert first_witness(x) in (2 * e, -2 * e)
    assert second_witness(y) in (2 * f, -2 * f)


def test_canonical_moving_null_line_under_level_three_subgroup_separates_lines():
    from sage_indefinite_port.indefinite.marked_lattices import canonical_moving_isotropic_line

    ambient = lattice([[0, 1], [1, 0]])
    e, f = tuple(ambient.module_generators())
    subgroup = ambient.O().congruence_preimage(3)
    first, first_witness = canonical_moving_isotropic_line(2 * e, acting_group=subgroup)
    negative, negative_witness = canonical_moving_isotropic_line(-2 * e, acting_group=subgroup)
    second, second_witness = canonical_moving_isotropic_line(2 * f, acting_group=subgroup)
    generator = next(iter(first.arrow().domain().module_generators()))
    assert first.arrow()(generator) == negative.arrow()(generator)
    assert first.arrow()(generator) != second.arrow()(generator)
    for (marked, witness), target in zip(((first, first_witness), (negative, negative_witness), (second, second_witness)),
                                         (2 * e, -2 * e, 2 * f), strict=True):
        assert witness in subgroup
        assert witness(marked.arrow()(generator)) in (target, -target)


def test_moving_nonisotropic_rank_one_markings_use_signed_pointed_transport():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    source = Lattices(ZZ)([[16]])
    first = ARROWS(source.hom([2 * (e + 2 * f)], codomain=plane))
    second = ARROWS(source.hom([-2 * (2 * e + f)], codomain=plane))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert witness(2 * (e + 2 * f)) in (2 * (2 * e + f), -2 * (2 * e + f))
    assert isometry(ARROWS, first, ARROWS(source.hom([e + 8 * f], codomain=plane))) is None


def test_nonisotropic_rank_one_setwise_stabilizer_contains_sign_reversal():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    source = Lattices(ZZ)([[16]])
    marked = ARROWS(source.hom([2 * (e + 2 * f)], codomain=plane))
    family = tuple(automorphism_group_generators(ARROWS, marked))
    assert all(g(2 * (e + 2 * f)) in (2 * (e + 2 * f), -2 * (e + 2 * f)) for g in family)
    assert any(g(e + 2 * f) == -(e + 2 * f) for g in family)
    assert all(plane.b(g(x), g(y)) == plane.b(x, y) for g in family for x in (e, f) for y in (e, f))


def test_nonisotropic_moving_line_classes_retain_all_contents():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonisotropic_line_classes, ARROWS, isometry

    plane = lattice([[0, 1], [1, 0]])
    for square in (8, -8):
        classes = nonisotropic_line_classes(plane, square)
        assert classes.cardinality() == cardinal(2)
        assert {int(marking.arrow()(marking.arrow().domain().module_generators()[0]).content()) for marking in classes} == {1, 2}
        first, second = tuple(classes)
        assert isometry(ARROWS, first, second) is None


def test_rank_two_homothetic_saturation_strata_have_prescribed_index():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import homothetic_isotropic_sublattice_classes

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    classes = homothetic_isotropic_sublattice_classes(ambient, 2, 2)
    assert classes.cardinality() == cardinal(1)
    marking = tuple(classes)[0]
    arrow = marking.arrow()
    saturation = arrow.saturation()
    assert saturation.inclusion().is_primitive()
    lifted = tuple(saturation.inclusion().lift(arrow(point)) for point in arrow.domain().module_generators())
    assert saturation.inclusion().domain().subobject_on(lifted).inclusion().cokernel().cardinality() == cardinal(4)


def test_rank_two_homothetic_stabilizer_is_saturation_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    primitive = ARROWS(source.hom([e, u], codomain=ambient))
    doubled = ARROWS(source.hom([2 * e, 2 * u], codomain=ambient))
    generators = tuple(automorphism_group_generators(ARROWS, doubled))
    assert all(doubled.arrow().is_in_image(g(2 * e)) for g in generators)
    assert all(doubled.arrow().is_in_image(g(2 * u)) for g in generators)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y) for g in generators for x in (e, f, u, v) for y in (e, f, u, v))
    primitive_generators = tuple(automorphism_group_generators(ARROWS, primitive))
    assert {tuple(g(x) for x in (e, f, u, v)) for g in generators} == {tuple(g(x) for x in (e, f, u, v)) for g in primitive_generators}


def test_nonhomothetic_rank_two_is_not_misclassified_as_scalar_saturation():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, _f, u, _v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    marked = ARROWS(source.hom([2 * e, u], codomain=ambient))
    family = tuple(automorphism_group_generators(ARROWS, marked))
    assert all(marked.arrow().is_in_image(g(2 * e)) and marked.arrow().is_in_image(g(u)) for g in family)
    assert all(marked.arrow().is_in_image((~g)(2 * e)) and marked.arrow().is_in_image((~g)(u)) for g in family)


def test_homothetic_rank_two_moving_transporter_uses_saturated_planes():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    first = ARROWS(source.hom([2 * e, 2 * u], codomain=ambient))
    second = ARROWS(source.hom([2 * f, 2 * v], codomain=ambient))
    third = ARROWS(source.hom([3 * f, 3 * v], codomain=ambient))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert second.arrow().saturation().inclusion().is_in_image(witness(e))
    assert second.arrow().saturation().inclusion().is_in_image(witness(u))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f, u, v) for y in (e, f, u, v))
    assert isometry(ARROWS, first, third) is None


def test_rank_two_index_two_submodules_are_quotiented_by_actual_levi_image():
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, _f, u, _v = tuple(ambient.module_generators())
    saturated = ambient.primitive_isotropic_subobject(e, u)
    classes = nonprimitive_sublattice_classes(saturated, (1, 2))
    assert len(classes) == 1
    arrow = classes[0]
    assert arrow.saturation().inclusion().is_primitive()
    assert arrow.codomain() is ambient
    assert arrow.saturation().inclusion().domain().subobject_on(tuple(arrow.saturation().inclusion().lift(arrow(x)) for x in arrow.domain().module_generators())).inclusion().cokernel().cardinality() == 2


def test_index_two_sublattices_respect_proper_integral_levi_image():
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    e, _f, u, _v = tuple(ambient.module_generators())
    saturated = ambient.primitive_isotropic_subobject(e, u)
    group = ambient.O().stabilizer(saturated)
    classes = nonprimitive_sublattice_classes(saturated, (1, 2), stabilizer_group=group)
    assert len(classes) == 2
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_composite_index_four_submodules_keep_smith_type():
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, _f, u, _v = tuple(ambient.module_generators())
    saturated = ambient.primitive_isotropic_subobject(e, u)
    cyclic = nonprimitive_sublattice_classes(saturated, (1, 4))
    homothetic = nonprimitive_sublattice_classes(saturated, (2, 2))
    assert len(cyclic) == len(homothetic) == 1
    assert cyclic[0].saturation().inclusion().is_primitive()
    assert homothetic[0].saturation().inclusion().is_primitive()


def test_nonhomothetic_rank_two_transport_uses_actual_finite_levi_orbit():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    first = ARROWS(source.hom([e, 2 * u], codomain=ambient))
    same = ARROWS(source.hom([e + u, 2 * u], codomain=ambient))
    other = ARROWS(source.hom([2 * e, u], codomain=ambient))
    witness = isometry(ARROWS, first, same)
    assert witness is not None
    assert all(same.arrow().is_in_image(witness(first.arrow()(x))) for x in first.arrow().domain().module_generators())
    assert isometry(ARROWS, first, other) is None


def test_nonhomothetic_rank_two_stabilizer_uses_actual_finite_orbit():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    arrow = source.hom([2 * e, u], codomain=ambient)
    family = tuple(automorphism_group_generators(ARROWS, ARROWS(arrow)))
    assert all(arrow.is_in_image(g(2 * e)) and arrow.is_in_image(g(u)) for g in family)
    assert all(arrow.is_in_image((~g)(2 * e)) and arrow.is_in_image((~g)(u)) for g in family)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y) for g in family for x in (e, f, u, v) for y in (e, f, u, v))


def test_public_orthogonal_group_nonprimitive_sublattice_classes():
    from sage_extensions.all import cardinal

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    e, _f, u, _v = tuple(ambient.module_generators())
    saturated = ambient.primitive_isotropic_subobject(e, u)
    classes = ambient.O().nonprimitive_sublattice_classes(saturated, (1, 2))
    assert classes.cardinality() == cardinal(2)
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_rank_three_index_two_finite_submodule_orbit():
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    e, _f, u, _v, w, _z = tuple(ambient.module_generators())
    saturated = ambient.primitive_isotropic_subobject(e, u, w)
    classes = nonprimitive_sublattice_classes(saturated, (1, 1, 2))
    assert len(classes) == 1
    assert classes[0].codomain() is ambient


def test_homothetic_rank_two_based_transport_and_pointwise_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([2 * e, 2 * u], codomain=ambient))
    second = category(source.hom([2 * f, 2 * v], codomain=ambient))
    third = category(source.hom([3 * f, 3 * v], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(2 * e) == 2 * f and witness(2 * u) == 2 * v
    assert isometry(category, first, third) is None
    family = tuple(automorphism_group_generators(category, first))
    assert all(g(2 * e) == 2 * e and g(2 * u) == 2 * u for g in family)


def test_supplied_trivial_stabilizer_keeps_all_index_two_submodules():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    e, _f, u, _v = tuple(ambient.module_generators())
    saturated = ambient.primitive_isotropic_subobject(e, u)
    trivial = ambient.O().subgroup(())
    classes = nonprimitive_sublattice_classes(saturated, (1, 2), stabilizer_group=trivial)
    assert len(classes) == 3
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_supplied_group_must_preserve_the_saturation():
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    labels = tuple(ambient.module_generating_set())
    exchange = ambient.O()(dict(zip(labels, (f, e, u, v), strict=True)))
    saturated = ambient.primitive_isotropic_subobject(e, u)
    with pytest.raises(ValueError):
        nonprimitive_sublattice_classes(saturated, (1, 2), stabilizer_group=ambient.O().subgroup((exchange,)))


def test_supplied_order_two_stabilizer_has_two_index_two_submodule_orbits():
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    labels = tuple(ambient.module_generating_set())
    exchange = ambient.O()(dict(zip(labels, (u, v, e, f), strict=True)))
    saturated = ambient.primitive_isotropic_subobject(e, u)
    subgroup = ambient.O().subgroup((exchange,))
    classes = nonprimitive_sublattice_classes(saturated, (1, 2), stabilizer_group=subgroup)
    assert len(classes) == 2
    assert all(arrow.codomain() is ambient and arrow.cokernel().is_finite() is False for arrow in classes)


def test_isotropic_lines_of_content_two_split_under_gamma_three():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_line_classes

    plane = lattice([[0, 1], [1, 0]])
    gamma = plane.O().congruence_preimage(3)
    classes = nonprimitive_isotropic_line_classes(plane, 2, acting_group=gamma)
    assert classes.cardinality() == cardinal(2)
    assert all(int(marking.arrow()(marking.arrow().domain().module_generators()[0]).content()) == 2
               for marking in classes)


def test_homothetic_rank_one_sublattices_split_under_gamma_three():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import homothetic_isotropic_sublattice_classes

    plane = lattice([[0, 1], [1, 0]])
    gamma = plane.O().congruence_preimage(3)
    classes = homothetic_isotropic_sublattice_classes(plane, 1, 2, acting_group=gamma)
    assert classes.cardinality() == cardinal(2)
    assert all(int(marking.arrow()(marking.arrow().domain().module_generators()[0]).content()) == 2
               for marking in classes)


def test_public_nonprimitive_rank_two_sublattice_classes_use_actual_levi_image():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    unimodular = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    mixed = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    full_classes = nonprimitive_isotropic_sublattice_classes(unimodular, 2, (1, 2))
    mixed_classes = nonprimitive_isotropic_sublattice_classes(mixed, 2, (1, 2))
    assert full_classes.cardinality() == cardinal(1)
    assert mixed_classes.cardinality() == cardinal(2)
    assert all(marked.arrow().codomain() is mixed for marked in mixed_classes)


def test_public_nonprimitive_sublattice_classes_retain_gamma_three_action():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    plane = lattice([[0, 1], [1, 0]])
    gamma = plane.O().congruence_preimage(3)
    classes = nonprimitive_isotropic_sublattice_classes(plane, 1, (2,), acting_group=gamma)
    assert classes.cardinality() == cardinal(2)
    assert all(int(marked.arrow()(marked.arrow().domain().module_generators()[0]).content()) == 2
               for marked in classes)


def test_nonsplit_nondegenerate_saturation_uses_glue_compatible_source_image():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    saturated = ambient.subobject_on([e + f, u + v])
    classes = nonprimitive_sublattice_classes(saturated, (1, 2))
    assert len(classes) == 2
    assert all(arrow.codomain() is ambient and arrow.saturation().inclusion().is_primitive()
               for arrow in classes)


def test_proper_glue_levi_image_distinguishes_equal_norm_source_summands():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = Lattices(ZZ)([[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 1]])
    e, f, u, v = tuple(ambient.module_generators())
    saturated = ambient.subobject_on([e + f, u])
    classes = nonprimitive_sublattice_classes(saturated, (1, 2))
    # The first norm-three vector has ambient divisibility one, the
    # second divisibility three.  O(S)=O(<3,3>) swaps them, but no
    # ambient isometry does.  Modulo two, the actual source image acts
    # trivially, leaving all three index-two submodules inequivalent.
    assert len(classes) == 3
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_mixed_degenerate_saturation_index_two_submodules_under_actual_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    saturated = ambient.subobject_on([e, z])
    classes = nonprimitive_sublattice_classes(saturated, (1, 2))
    assert len(classes) == 2
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_split_radical_saturation_full_hyperbolic_quotient_and_radical_line():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    saturated = ambient.subobject_on([e, f, z])
    classes = nonprimitive_sublattice_classes(saturated, (1, 1, 2))
    assert len(classes) == 3
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_nonprimitive_mixed_degenerate_sublattice_setwise_schreier_generators():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[8, 0], [0, 0]])
    marked = ARROWS(source.hom([2 * e, z], codomain=ambient))
    generators = tuple(automorphism_group_generators(ARROWS, marked))
    assert all(marked.arrow().is_in_image(g(2 * e)) and marked.arrow().is_in_image(g(z))
               for g in generators)
    assert any(g(e) != e for g in generators)


def test_split_radical_saturation_stabilizer_handles_nonstandard_graph():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    saturated = ambient.subobject_on([e + w, z + w])
    classes = nonprimitive_sublattice_classes(saturated, (1, 2))
    assert len(classes) == 2
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_split_radical_saturation_rank_two_radical_summand_complete_action():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[2, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    e, z, w, t = tuple(ambient.module_generators())
    saturated = ambient.subobject_on([e + t, z + t, w])
    classes = nonprimitive_sublattice_classes(saturated, (1, 1, 2))
    assert len(classes) == 2
    assert all(arrow.codomain() is ambient for arrow in classes)


def test_split_radical_full_saturation_recovers_complete_ambient_group_image():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_sublattice_classes

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    saturated = ambient.subobject_on(tuple(ambient.module_generators()))
    classes = nonprimitive_sublattice_classes(saturated, (1, 1, 2))
    assert len(classes) == 2
