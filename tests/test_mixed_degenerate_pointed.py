"""Exact split-radical mixed pointed transporters and congruence obstructions."""
from tests.acceptance.consumer import lattice


def test_mixed_radical_shear_positive_and_negative_transporters():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0], [0, 0]])
    e, z = tuple(ambient.module_generators())
    for source, target in ((e, e + z), (5 * e + z, 5 * e + 6 * z)):
        witness = RECURSION.pointed_isometry(source, target)
        assert witness is not None and witness(source) == target
        assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
                   for x in (e, z) for y in (e, z))
        assert (~witness)(target) == source
    assert RECURSION.pointed_isometry(5 * e + z, 5 * e + 2 * z) is None


def test_rank_two_radical_uses_actual_gl_two_mod_five_and_shear():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    for source, target in ((e + z, e + w), (5 * e + z, 5 * e + 2 * z)):
        witness = RECURSION.pointed_isometry(source, target)
        assert witness is not None and witness(source) == target
        assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
                   for x in (e, z, w) for y in (e, z, w))
    assert RECURSION.pointed_isometry(5 * e + z, 5 * e + 5 * z) is None


def test_cross_lattice_mixed_radical_transport_uses_ambient_isometry_first():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    source_lattice = lattice([[2, 0], [0, 0]])
    target_lattice = lattice([[0, 0], [0, 2]])
    e, z = tuple(source_lattice.module_generators())
    w, f = tuple(target_lattice.module_generators())
    witness = RECURSION.pointed_isometry(e, f + w)
    assert witness is not None
    assert witness(e) == f + w
    assert (~witness)(f + w) == e
    assert all(target_lattice.b(witness(x), witness(y)) == source_lattice.b(x, y)
               for x in (e, z) for y in (e, z))


def test_public_coslice_mixed_degenerate_marking_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0], [0, 0]])
    e, z = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e], codomain=ambient))
    second = category(source.hom([e + z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(e) == e + z


def test_mixed_degenerate_pointwise_stabilizer_has_exact_radical_sign_group():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0], [0, 0]])
    e, z = tuple(ambient.module_generators())
    family = tuple(RECURSION.pointed_automorphism_generators(e))
    assert all(g(e) == e for g in family)
    assert ambient.O().subgroup(family).cardinality() == cardinal(2)


def test_mixed_rank_two_radical_pointwise_stabilizer_preserves_marking():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    vector = 5 * e + z
    family = tuple(RECURSION.pointed_automorphism_generators(vector))
    assert all(g(vector) == vector for g in family)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in family for x in (e, z, w) for y in (e, z, w))
    assert any(g(w) != w for g in family)


def test_mixed_isotropic_nonradical_vector_uses_shear_not_parabolic_reduction():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    witness = RECURSION.pointed_isometry(e + z, e + 2 * z)
    assert witness is not None and witness(e + z) == e + 2 * z
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, z) for y in (e, f, z))


def test_public_degenerate_nonzero_norm_orbit_classes_split_by_radical_residue():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    ambient = lattice([[2, 0], [0, 0]])
    source = Lattices(ZZ)([[50]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = pointed_nonzero_norm_classes(category, ambient, 50)
    assert classes.cardinality() == cardinal(3)
    assert all(marked.arrow()(source.module_generators()[0]).q() == ZZ(50) for marked in classes)
    assert {int(marked.arrow()(source.module_generators()[0]).content()) for marked in classes} == {1, 5}


def test_degenerate_fixed_content_fibres_separate_primitive_and_imprimitive():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes

    ambient = lattice([[2, 0], [0, 0]])
    source = Lattices(ZZ)([[50]])
    category = Lattices(ZZ).CosliceUnder(source)
    assert pointed_content_classes(category, ambient, 50, 1).cardinality() == cardinal(2)
    assert pointed_content_classes(category, ambient, 50, 5).cardinality() == cardinal(1)
    assert pointed_content_classes(category, ambient, 50, 2).is_empty() is True


def test_degenerate_isotropic_orbit_iterator_is_fair_across_radical_and_mixed_contents():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_pointed_isotropic_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = tuple(islice(iter_pointed_isotropic_classes(category, ambient, include_zero=True), 7))
    vectors = tuple(marked.arrow()(source.module_generators()[0]) for marked in classes)
    assert vectors[0] == ambient.zero()
    assert all(point.q() == ZZ(0) for point in vectors)
    assert any(point in ambient.radical() for point in vectors[1:])
    assert any(point not in ambient.radical() for point in vectors[1:])
    assert {int(point.content()) for point in vectors[1:]} >= {1, 2}


def test_mixed_degenerate_primitive_null_content_has_infinitely_many_orbits():
    from itertools import islice
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        iter_pointed_content_classes, pointed_content_classes, isometry,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    with pytest.raises(ValueError, match="may be infinite"):
        pointed_content_classes(category, ambient, 0, 1)
    classes = tuple(islice(iter_pointed_content_classes(category, ambient, 0, 1), 6))
    assert len(classes) == 6
    vectors = tuple(marked.arrow()(source.module_generators()[0]) for marked in classes)
    assert all(point.q() == 0 and int(point.content()) == 1 for point in vectors)
    assert all(isometry(category, first, second) is None
               for i, first in enumerate(classes) for second in classes[i + 1:])


def test_mixed_degenerate_anisotropic_quotient_has_one_radical_null_orbit_per_content():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import (
        pointed_content_classes, iter_pointed_content_classes,
    )

    ambient = lattice([[2, 0], [0, 0]])
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    for content in (1, 2, 5):
        classes = pointed_content_classes(category, ambient, 0, content)
        assert classes.cardinality() == cardinal(1)
        marked = next(iter(classes))
        point = marked.arrow()(source.module_generators()[0])
        assert point in ambient.radical()
        assert int(point.content()) == content
        assert len(tuple(iter_pointed_content_classes(category, ambient, 0, content))) == 1


def test_mixed_degenerate_indefinite_anisotropic_quotient_has_only_radical_null_vectors():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes

    # x^2 - 3y^2 = 0 has no nonzero rational solution.
    ambient = lattice([[2, 0, 0], [0, -6, 0], [0, 0, 0]])
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = pointed_content_classes(category, ambient, 0, 3)
    assert classes.cardinality() == cardinal(1)
    point = next(iter(classes)).arrow()(source.module_generators()[0])
    assert point in ambient.radical() and int(point.content()) == 3


def test_mixed_degenerate_nonprimitive_fixed_content_stream_is_exhaustive_across_projection_contents():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        iter_pointed_content_classes, isometry,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    representatives = tuple(islice(iter_pointed_content_classes(category, ambient, 0, 2), 5))
    assert len(representatives) == 5
    points = tuple(marked.arrow()(source.module_generators()[0]) for marked in representatives)
    assert all(int(point.content()) == 2 and point.q() == 0 for point in points)
    assert any(point in ambient.radical() for point in points)
    assert any(point not in ambient.radical() for point in points)
    assert all(isometry(category, left, right) is None
               for i, left in enumerate(representatives) for right in representatives[i + 1:])


def test_public_primitive_pointed_classifier_dispatches_degenerate_quotients():
    import pytest
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[2, 0], [0, 0]])
    for square, expected in ((0, 1), (2, 1), (50, 2)):
        source = Lattices(ZZ)([[square]])
        category = Lattices(ZZ).CosliceUnder(source)
        classes = primitive_marking_classes(category, source, ambient)
        assert classes.cardinality() == cardinal(expected)
        assert all(int(marked.arrow()(source.module_generators()[0]).content()) == 1
                   for marked in classes)
    isotropic_quotient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    with pytest.raises(ValueError, match="may be infinite"):
        primitive_marking_classes(category, source, isotropic_quotient)


def test_public_primitive_moving_null_line_classifier_dispatches_degenerate_quotients():
    import pytest
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import (
        ARROWS, primitive_marking_classes,
    )

    source = Lattices(ZZ)([[0]])
    ambient = lattice([[2, 0], [0, 0]])
    classes = primitive_marking_classes(ARROWS, source, ambient)
    assert classes.cardinality() == cardinal(1)
    arrow = next(iter(classes)).arrow()
    point = arrow(next(iter(arrow.domain().module_generators())))
    assert point in ambient.radical() and int(point.content()) == 1
    isotropic_quotient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    with pytest.raises(ValueError, match="may be infinite"):
        primitive_marking_classes(ARROWS, source, isotropic_quotient)


def test_public_primitive_radical_plane_in_anisotropic_degenerate_lattice():
    import pytest
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    classes = primitive_marking_classes(ARROWS, source, ambient)
    assert classes.cardinality() == cardinal(1)
    arrow = next(iter(classes)).arrow()
    assert all(arrow(point) in ambient.radical() for point in arrow.domain().module_generators())
    assert primitive_marking_classes(ARROWS, source, lattice([[2, 0], [0, 0]])).is_empty()
    isotropic_quotient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    with pytest.raises(ValueError, match="may be infinite"):
        primitive_marking_classes(ARROWS, source, isotropic_quotient)


def test_public_primitive_based_radical_plane_in_anisotropic_degenerate_lattice():
    import pytest
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = primitive_marking_classes(category, source, ambient)
    assert classes.cardinality() == cardinal(1)
    arrow = next(iter(classes)).arrow()
    assert all(arrow(point) in ambient.radical() for point in source.module_generators())
    assert primitive_marking_classes(category, source, lattice([[2, 0], [0, 0]])).is_empty()
    isotropic_quotient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    with pytest.raises(ValueError, match="may be infinite"):
        primitive_marking_classes(category, source, isotropic_quotient)


def test_public_primitive_radical_flag_in_anisotropic_degenerate_lattice():
    import pytest
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, primitive_marking_classes

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    plane = Lattices(ZZ)([[0, 0], [0, 0]])
    first = plane.primitive_isotropic_subobject(next(iter(plane.module_generators())))
    classes = primitive_marking_classes(FLAGS, first.inclusion(), ambient)
    assert classes.cardinality() == cardinal(1)
    arrows = tuple(next(iter(classes)).arrows())
    assert tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2)
    assert arrows[-1].codomain() is ambient
    assert primitive_marking_classes(FLAGS, first.inclusion(), lattice([[2, 0], [0, 0]])).is_empty()
    isotropic_quotient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    with pytest.raises(ValueError, match="may be infinite"):
        primitive_marking_classes(FLAGS, first.inclusion(), isotropic_quotient)


def test_primitive_radical_markings_are_covariant_under_reordered_coordinate_labels():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import (
        ARROWS, primitive_marking_classes, iter_pointed_isotropic_classes,
        automorphism_group_generators, FLAGS,
    )

    ambient = Lattices(ZZ)([[0, 0, 0], [0, 2, 0], [0, 0, 0]],
                           module_generators=("radical_first", "positive", "radical_last"))
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    for category in (ARROWS, Lattices(ZZ).CosliceUnder(source)):
        classes = primitive_marking_classes(category, source, ambient)
        assert classes.cardinality() == cardinal(1)
        arrow = next(iter(classes)).arrow()
        assert all(arrow(point) in ambient.radical() for point in arrow.domain().module_generators())
        if category is not ARROWS:
            generators = tuple(automorphism_group_generators(category, next(iter(classes))))
            assert generators
            assert all(all(g(arrow(point)) == arrow(point)
                           for point in arrow.domain().module_generators()) for g in generators)
    pointed_source = Lattices(ZZ)([[0]])
    pointed_category = Lattices(ZZ).CosliceUnder(pointed_source)
    radical_classes = tuple(islice(iter_pointed_isotropic_classes(pointed_category, ambient), 3))
    assert tuple(int(marked.arrow()(pointed_source.module_generators()[0]).content())
                 for marked in radical_classes) == (1, 2, 3)
    assert all(marked.arrow()(pointed_source.module_generators()[0]) in ambient.radical()
               for marked in radical_classes)
    plane = Lattices(ZZ)([[0, 0], [0, 0]])
    line = plane.primitive_isotropic_subobject(next(iter(plane.module_generators())))
    flag = next(iter(primitive_marking_classes(FLAGS, line.inclusion(), ambient)))
    flag_arrows = tuple(flag.arrows())
    first = flag_arrows[-1](flag_arrows[0](next(iter(flag_arrows[0].domain().module_generators()))))
    assert all(g(first) in (first, -first) for g in automorphism_group_generators(FLAGS, flag))


def test_full_radical_plane_moving_stabilizer_is_entire_orthogonal_group():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        ARROWS, primitive_marking_classes, automorphism_group_generators,
    )

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    marked = next(iter(primitive_marking_classes(ARROWS, source, ambient)))
    generators = tuple(automorphism_group_generators(ARROWS, marked))
    assert generators
    radical = ambient.radical()
    assert all(all(radical.inclusion().is_in_image(g(radical.inclusion()(point)))
                   for point in radical.module_generators()) for g in generators)
    assert any(g(tuple(ambient.module_generators())[1]) != tuple(ambient.module_generators())[1]
               for g in generators)


def test_full_radical_plane_based_stabilizer_fixes_each_radical_generator():
    from sage.all import ZZ as SageZZ, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        primitive_marking_classes, automorphism_group_generators,
    )

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = next(iter(primitive_marking_classes(category, source, ambient)))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators
    fixed = tuple(marked.arrow()(point) for point in source.module_generators())
    assert all(all(g(point) == point for point in fixed) for g in generators)
    assert any(g(tuple(ambient.module_generators())[0]) != tuple(ambient.module_generators())[0]
               for g in generators)
    e = tuple(ambient.module_generators())[0]
    labels = tuple(ambient.module_generating_set())
    shear_rows = [tuple(int(g(e).to_vector()(label)) for label in labels[1:])
                  for g in generators]
    assert matrix(SageZZ, shear_rows).rank() == 2
    assert any(int(g(e).to_vector()(labels[0])) == -1 for g in generators)


def test_full_radical_plane_flag_stabilizer_preserves_line_and_plane():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, primitive_marking_classes, automorphism_group_generators,
    )

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    plane = Lattices(ZZ)([[0, 0], [0, 0]])
    line = plane.primitive_isotropic_subobject(next(iter(plane.module_generators())))
    marked = next(iter(primitive_marking_classes(FLAGS, line.inclusion(), ambient)))
    generators = tuple(automorphism_group_generators(FLAGS, marked))
    assert generators
    arrows = tuple(marked.arrows())
    first = arrows[-1](arrows[0](next(iter(arrows[0].domain().module_generators()))))
    assert all(g(first) in (first, -first) for g in generators)
    assert any(g(first) == -first for g in generators)


def test_proper_radical_line_moving_stabilizer_preserves_its_saturation():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        ARROWS, automorphism_group_generators,
    )

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    for radical_vector in (z, z + w, 2 * z + 3 * w):
        line = ambient.primitive_isotropic_subobject(radical_vector)
        marked = ARROWS.object(line.inclusion())
        generators = tuple(automorphism_group_generators(ARROWS, marked))
        assert generators
        assert all(g(radical_vector) in (radical_vector, -radical_vector) for g in generators)
        assert any(g(radical_vector) == -radical_vector for g in generators)


def test_noncoordinate_radical_flag_stabilizer_preserves_primitive_line():
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, automorphism_group_generators,
    )

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    for vector, companion in ((z + w, z), (2 * z + 3 * w, z + w)):
        plane = ambient.primitive_isotropic_subobject(vector, companion)
        inclusion = plane.inclusion()
        line = inclusion.domain().primitive_isotropic_subobject(inclusion.lift(vector))
        marked = FLAGS((line.inclusion(), inclusion))
        generators = tuple(automorphism_group_generators(FLAGS, marked))
        assert generators
        assert all(g(vector) in (vector, -vector) for g in generators)
        assert any(g(vector) == -vector for g in generators)


def test_radical_plane_flags_have_exact_line_transporter():
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    plane = ambient.primitive_isotropic_subobject(z, w)
    inclusion = plane.inclusion()

    def flag(vector):
        line = inclusion.domain().primitive_isotropic_subobject(inclusion.lift(vector))
        return FLAGS((line.inclusion(), inclusion))

    source, target = flag(z), flag(2 * z + 3 * w)
    witness = isometry(FLAGS, source, target)
    assert witness is not None
    assert witness(z) in (2 * z + 3 * w, -2 * z - 3 * w)
    target_ambient = lattice([[0, 0, 0], [0, 2, 0], [0, 0, 0]])
    target_z, _, target_w = tuple(target_ambient.module_generators())
    target_plane = target_ambient.primitive_isotropic_subobject(target_z, target_w)
    target_inclusion = target_plane.inclusion()
    target_line = target_inclusion.domain().primitive_isotropic_subobject(
        target_inclusion.lift(2 * target_z + 3 * target_w)
    )
    cross = FLAGS((target_line.inclusion(), target_inclusion))
    cross_witness = isometry(FLAGS, source, cross)
    assert cross_witness is not None
    assert cross_witness(z) in (2 * target_z + 3 * target_w, -2 * target_z - 3 * target_w)


def test_moving_radical_lines_have_exact_integral_transporters():
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    first = ARROWS.object(ambient.primitive_isotropic_subobject(z).inclusion())
    second = ARROWS.object(ambient.primitive_isotropic_subobject(2 * z + 3 * w).inclusion())
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert witness(z) in (2 * z + 3 * w, -2 * z - 3 * w)


def test_moving_full_radical_planes_have_exact_integral_transporters():
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    first = ARROWS.object(ambient.primitive_isotropic_subobject(z, w).inclusion())
    second = ARROWS.object(ambient.primitive_isotropic_subobject(z + w, z).inclusion())
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert all(witness(point) in ambient.radical() for point in (z, w))
    target_ambient = lattice([[0, 0, 0], [0, 2, 0], [0, 0, 0]])
    target_z, _, target_w = tuple(target_ambient.module_generators())
    cross = ARROWS.object(target_ambient.primitive_isotropic_subobject(target_z, target_w).inclusion())
    cross_witness = isometry(ARROWS, first, cross)
    assert cross_witness is not None
    assert all(cross_witness(point) in target_ambient.radical() for point in (z, w))


def test_proper_radical_plane_stabilizer_in_rank_three_radical():
    from sage_indefinite_port.indefinite.marked_lattices import (
        ARROWS, automorphism_group_generators,
    )

    ambient = lattice([[2, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    z, w, u = tuple(ambient.module_generators())[1:]
    plane = ambient.primitive_isotropic_subobject(z + w, z + u)
    marked = ARROWS.object(plane.inclusion())
    generators = tuple(automorphism_group_generators(ARROWS, marked))
    assert generators
    assert all(all(plane.inclusion().is_in_image(g(point)) for point in (z + w, z + u))
               for g in generators)


def test_proper_radical_based_plane_stabilizer_fixes_frame_pointwise():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = lattice([[2, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    z, w, u = tuple(ambient.module_generators())[1:]
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category.object(source.hom([z + w, z + u], codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators
    assert all(g(z + w) == z + w and g(z + u) == z + u for g in generators)
    assert any(g(u) != u for g in generators)


def test_nonprimitive_based_radical_plane_stabilizer_fixes_saturation_pointwise():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category.object(source.hom([2 * z, 3 * w], codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators
    assert all(g(z) == z and g(w) == w for g in generators)


def test_nonhomothetic_based_radical_frames_have_exact_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    z, w = tuple(ambient.module_generators())[1:]
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category.object(source.hom([2 * z, 3 * w], codomain=ambient))
    second = category.object(source.hom([2 * (z + w), 3 * w], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(2 * z) == 2 * (z + w) and witness(3 * w) == 3 * w

    # Equal Smith invariants do not imply based equivalence: the source
    # frame's column lattice, rather than its unlabelled index, is fixed.
    incompatible = category.object(source.hom([3 * z, 2 * w], codomain=ambient))
    assert isometry(category, first, incompatible) is None
    target_ambient = lattice([[0, 0, 0], [0, 2, 0], [0, 0, 0]])
    target_z, _, target_w = tuple(target_ambient.module_generators())
    cross = category.object(source.hom([2 * target_z, 3 * target_w], codomain=target_ambient))
    cross_witness = isometry(category, first, cross)
    assert cross_witness is not None
    assert cross_witness(2 * z) == 2 * target_z
    assert cross_witness(3 * w) == 3 * target_w


def test_mixed_degenerate_isotropic_rank_above_rational_witt_index_is_empty():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, FLAGS, primitive_marking_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    assert primitive_marking_classes(category, source, ambient).is_empty()
    assert primitive_marking_classes(ARROWS, source, ambient).is_empty()
    first = source.primitive_isotropic_subobject(*tuple(source.module_generators())[:2])
    assert primitive_marking_classes(FLAGS, first.inclusion(), ambient).is_empty()


def test_moving_maximal_null_plane_in_u_plus_radical_is_one_orbit():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    classes = primitive_marking_classes(ARROWS, source, ambient)
    assert classes.cardinality() == cardinal(1)
    arrow = next(iter(classes)).arrow()
    assert arrow.is_primitive()
    assert all(ambient.b(arrow(x), arrow(y)) == 0
               for x in arrow.domain().module_generators()
               for y in arrow.domain().module_generators())
    radical = ambient.radical()
    assert all(arrow.saturation().inclusion().is_in_image(radical.inclusion()(point))
               for point in radical.module_generators())


def test_moving_maximal_null_three_plane_in_u_plus_rank_two_radical():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    source = Lattices(ZZ)([[0] * 3 for _ in range(3)])
    classes = primitive_marking_classes(ARROWS, source, ambient)
    assert classes.cardinality() == cardinal(1)
    arrow = next(iter(classes)).arrow()
    assert arrow.is_primitive()
    radical = ambient.radical()
    assert all(arrow.saturation().inclusion().is_in_image(radical.inclusion()(point))
               for point in radical.module_generators())


def test_moving_maximal_null_planes_in_u_plus_radical_have_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    first = ARROWS.object(ambient.primitive_isotropic_subobject(e, z).inclusion())
    second = ARROWS.object(ambient.primitive_isotropic_subobject(f, z).inclusion())
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert all(second.arrow().is_in_image(witness(first.arrow()(point)))
               for point in first.arrow().domain().module_generators())


def test_moving_maximal_null_plane_in_u_plus_radical_has_exact_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    marked = ARROWS.object(ambient.primitive_isotropic_subobject(e, z).inclusion())
    generators = tuple(automorphism_group_generators(ARROWS, marked))
    assert generators
    failures = tuple((g(e), g(z)) for g in generators
                     if not (marked.arrow().is_in_image(g(e)) and marked.arrow().is_in_image(g(z))))
    assert not failures, failures
    assert any(g(e) == e + z for g in generators)
    assert any(g(f) == f + z for g in generators)
    assert any(g(z) == -z for g in generators)


def test_rectangular_based_radical_frame_transporter_in_rank_three_radical():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    z, w, u = tuple(ambient.module_generators())[1:]
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category.object(source.hom([2 * (z + w), 3 * (z + u)], codomain=ambient))
    second = category.object(source.hom([2 * (z + w + u), 3 * (z + u)], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(2 * (z + w)) == 2 * (z + w + u)
    assert witness(3 * (z + u)) == 3 * (z + u)


def test_mixed_radical_rank_one_stabilizers_respect_residue_and_shear():
    from sage_extensions.all import cardinal
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0], [0, 0]])
    e, z = tuple(ambient.module_generators())
    for point, expected in ((e + z, 2), (5 * e + z, 1), (5 * e, 2)):
        family = tuple(RECURSION.pointed_automorphism_generators(point))
        assert all(g(point) == point for g in family)
        assert ambient.O().subgroup(family).cardinality() == cardinal(expected)


def test_rank_two_radical_norm_fifty_orbits_use_actual_gl_two_mod_five():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes, isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    source = Lattices(ZZ)([[50]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = tuple(pointed_nonzero_norm_classes(category, ambient, 50))
    assert len(classes) == 2
    assert {int(marked.arrow()(source.module_generators()[0]).content()) for marked in classes} == {1, 5}
    assert all(isometry(category, first, second) is None
               for i, first in enumerate(classes) for second in classes[i + 1:])


def test_rank_two_based_mixed_degenerate_gram_embedding_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e, z], codomain=ambient))
    second = category(source.hom([e + w, z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(e) == e + w and witness(z) == z


def test_rank_two_based_mixed_gram_radical_generator_transport_and_obstruction():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 2]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([z, e], codomain=ambient))
    second = category(source.hom([w, e + 2 * z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(z) == w and witness(e) == e + 2 * z
    assert (~witness)(w) == z
    obstructed = category(source.hom([2 * w, e + 2 * z], codomain=ambient))
    assert isometry(category, first, obstructed) is None


def test_rank_two_based_mixed_gram_nonprimitive_quotient_radical_sign_and_shear():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0], [0, 0]])
    e, z = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[8, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([2 * e + z, z], codomain=ambient))
    second = category(source.hom([2 * e + 3 * z, z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(2 * e + z) == 2 * e + 3 * z and witness(z) == z
    impossible = category(source.hom([2 * e + 2 * z, z], codomain=ambient))
    assert isometry(category, first, impossible) is None


def test_rank_two_based_mixed_gram_nonprimitive_radical_negation():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0], [0, 0]])
    e, z = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[8, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([2 * e + z, z], codomain=ambient))
    second = category(source.hom([2 * e + 3 * z, -z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(2 * e + z) == 2 * e + 3 * z and witness(z) == -z


def test_rank_two_based_mixed_gram_rank_two_radical_constrained_congruence():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[50, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([5 * e + w, z], codomain=ambient))
    second = category(source.hom([5 * e + w + 5 * z, z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(5 * e + w) == 5 * e + w + 5 * z and witness(z) == z
    # Fixing z forces the radical matrix first row to be (1,0).
    # The w-coordinate of the first marking is therefore a unit modulo 5.
    third = category(source.hom([5 * e + 2 * w, z], codomain=ambient))
    assert isometry(category, first, third) is None


def test_rank_two_radical_stabilizer_transporter_lifts_nontrivial_mod_five_word():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[50, 0], [0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([5 * e + w, z], codomain=ambient))
    second = category(source.hom([5 * e + w + z, z], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(5 * e + w) == 5 * e + w + z and witness(z) == z
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, z, w) for y in (e, z, w))


def test_full_quotient_based_mixed_rank_three_gram_marking_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e, z, w], codomain=ambient))
    second = category(source.hom([e + z + w, z + w, w], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(e) == e + z + w and witness(z) == z + w and witness(w) == w
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, z, w) for y in (e, z, w))
    nonprimitive = category(source.hom([e, 2 * z, w], codomain=ambient))
    assert isometry(category, first, nonprimitive) is None


def test_full_quotient_based_mixed_rank_four_gram_with_two_radical_generators():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, f, z, w, t = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e, f, z, w], codomain=ambient))
    second = category(source.hom([e + z, f + w, w, z + t], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert all(witness(first.arrow()(x)) == second.arrow()(x) for x in source.module_generators())
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())


def test_full_quotient_based_mixed_gram_across_distinct_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    source_lattice = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    target_lattice = lattice([[0, 0, 0], [0, 2, 0], [0, 0, 0]])
    e, z, w = tuple(source_lattice.module_generators())
    u, f, v = tuple(target_lattice.module_generators())
    source = Lattices(ZZ)([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e, z, w], codomain=source_lattice))
    second = category(source.hom([f + u, u, v], codomain=target_lattice))
    witness = isometry(category, first, second)
    assert witness is not None
    assert all(witness(first.arrow()(x)) == second.arrow()(x) for x in source.module_generators())
    assert all(target_lattice.b(witness(x), witness(y)) == source_lattice.b(x, y)
               for x in source_lattice.module_generators() for y in source_lattice.module_generators())


def test_full_quotient_based_mixed_gram_nonprimitive_radical_frame_obstruction():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e, 2 * z, w], codomain=ambient))
    second = category(source.hom([e + z, 2 * z + 2 * w, w], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert all(witness(first.arrow()(x)) == second.arrow()(x) for x in source.module_generators())
    impossible = category(source.hom([e, 2 * z + w, w], codomain=ambient))
    assert isometry(category, first, impossible) is None


def test_full_rank_based_mixed_gram_nonprimitive_quotient_exact_shear_integrality():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    e, z, w = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[8, 0, 0], [0, 0, 0], [0, 0, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([2 * e + z, 2 * z, w], codomain=ambient))
    second = category(source.hom([2 * e + 3 * z, 2 * z, w], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert all(witness(first.arrow()(x)) == second.arrow()(x) for x in source.module_generators())
    impossible = category(source.hom([2 * e + 2 * z, 2 * z, w], codomain=ambient))
    assert isometry(category, first, impossible) is None


def test_public_degenerate_primitive_radical_line_orbit_stream():
    from sage_indefinite_port.indefinite.marked_lattices import iter_primitive_isotropic_line_classes
    from sage_extensions.all import ZZ

    ambient = lattice([[2, 0], [0, 0]])
    representative = next(iter_primitive_isotropic_line_classes(ambient))
    inclusion = representative.arrow()
    vector = inclusion(next(iter(inclusion.domain().module_generators())))
    assert vector.q() == ZZ(0) and vector.content() == ZZ(1)
    assert vector in ambient.radical()


def test_public_degenerate_rank_two_radical_line_orbit_stream():
    from sage_indefinite_port.indefinite.marked_lattices import iter_primitive_isotropic_line_classes
    from sage_extensions.all import ZZ

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    representative = next(iter_primitive_isotropic_line_classes(ambient))
    inclusion = representative.arrow()
    vector = inclusion(next(iter(inclusion.domain().module_generators())))
    assert vector.q() == ZZ(0) and vector.content() == ZZ(1)
    assert vector in ambient.radical()


def test_public_pure_radical_primitive_null_line_stream():
    from sage_indefinite_port.indefinite.marked_lattices import iter_primitive_isotropic_line_classes
    from sage_extensions.all import Lattices, ZZ

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    representative = next(iter_primitive_isotropic_line_classes(ambient))
    arrow = representative.arrow()
    vector = arrow(next(iter(arrow.domain().module_generators())))
    assert vector.q() == ZZ(0) and vector.content() == ZZ(1)
