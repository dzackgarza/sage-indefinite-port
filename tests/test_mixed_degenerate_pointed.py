"""Exact split-radical mixed pointed transporters and congruence obstructions."""
import pytest
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


def test_rank_zero_pointed_null_stream_has_only_optional_zero_orbit():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_pointed_isotropic_classes

    ambient = Lattices(ZZ)(0)
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    assert tuple(iter_pointed_isotropic_classes(category, ambient)) == ()
    marked = tuple(iter_pointed_isotropic_classes(category, ambient, include_zero=True))
    assert len(marked) == 1
    assert marked[0].arrow()(source.module_generators()[0]) == ambient.zero()


def test_rank_zero_fixed_pointed_content_has_no_nonzero_vectors():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import pointed_content_classes, iter_pointed_content_classes

    ambient = Lattices(ZZ)(0)
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    assert tuple(pointed_content_classes(category, ambient, 0, 1)) == ()
    assert tuple(iter_pointed_content_classes(category, ambient, 0, 2)) == ()


def test_rank_zero_nonzero_norm_pointed_classes_are_empty():
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    ambient = Lattices(ZZ)(0)
    source = Lattices(ZZ)([[2]])
    category = Lattices(ZZ).CosliceUnder(source)
    assert tuple(pointed_nonzero_norm_classes(category, ambient, 2)) == ()
    with pytest.raises(ValueError, match="source form"):
        pointed_nonzero_norm_classes(category, ambient, 0)


def test_degenerate_nonzero_norm_orbit_classifier_rejects_incompatible_source():
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    ambient = Lattices(ZZ)([[2, 0], [0, 0]])
    wrong_source = Lattices(ZZ)([[8]])
    category = Lattices(ZZ).CosliceUnder(wrong_source)
    with pytest.raises(ValueError, match="source form"):
        pointed_nonzero_norm_classes(category, ambient, 50)


def test_infinite_fixed_content_iterator_rejects_incompatible_pointed_source():
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_pointed_content_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    wrong_source = Lattices(ZZ)([[2]])
    wrong_category = Lattices(ZZ).CosliceUnder(wrong_source)
    with pytest.raises(ValueError, match="source form"):
        next(iter_pointed_content_classes(wrong_category, ambient, 0, 1))


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


def test_odd_rank_radical_flag_so_orbits_equal_full_orthogonal_orbits():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, primitive_marking_classes

    ambient = Lattices(ZZ)([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    plane = Lattices(ZZ)([[0, 0], [0, 0]])
    line = plane.primitive_isotropic_subobject(next(iter(plane.module_generators())))
    full = tuple(primitive_marking_classes(FLAGS, line.inclusion(), ambient))
    special = tuple(primitive_marking_classes(FLAGS, line.inclusion(), ambient, acting_group=ambient.SO()))
    assert len(full) == len(special) == 1


def test_mixed_degenerate_maximal_moving_null_plane_is_finite_despite_infinite_lines():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    classes = tuple(primitive_marking_classes(ARROWS, source, ambient))
    assert len(classes) == 1
    plane = classes[0].arrow()
    assert plane.is_primitive()
    assert plane.domain().rank() == 2
    assert all(ambient.b(plane(x), plane(y)) == 0
               for x in plane.domain().module_generators() for y in plane.domain().module_generators())
    radical = ambient.radical().inclusion()
    assert all(plane.is_in_image(radical(point)) for point in radical.domain().module_generators())


def test_mixed_degenerate_full_radical_maximal_isotropic_stabilizer_generators():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, automorphism_group_generators

    for radical_rank in (1, 2):
        ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * radical_rank for _ in range(radical_rank)])
        basis = tuple(ambient.module_generators())
        e, f = basis[:2]
        radical = basis[2:]
        plane = ambient.primitive_isotropic_subobject(e, *radical)
        generators = INDUCTION.plane_automorphism_generators(plane)
        assert generators
        assert all(plane.inclusion().is_in_image(generator(point))
                   for generator in generators for point in (e, *radical))
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for generator in generators for x in basis for y in basis)
        marked = ARROWS.object(plane.inclusion())
        public = tuple(automorphism_group_generators(ARROWS, marked))
        assert public
        assert all(plane.inclusion().is_in_image(generator(point))
                   for generator in public for point in (e, *radical))


def test_mixed_degenerate_full_radical_plane_transporter_lifts_quotient_isometry():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    for k in (1, 2):
        ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * k for _ in range(k)])
        e, f, *radical = tuple(ambient.module_generators())
        first = ambient.primitive_isotropic_subobject(e, *radical)
        second = ambient.primitive_isotropic_subobject(f, *radical)
        witness = INDUCTION.plane_isometry(first, second)
        assert witness is not None
        assert all(second.inclusion().is_in_image(witness(point)) for point in (e, *radical))
        assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
                   for x in (e, f, *radical) for y in (e, f, *radical))


def test_mixed_degenerate_full_radical_plane_transporter_between_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    source = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    target = Lattices(ZZ)([[0, -1, 0], [-1, 0, 0], [0, 0, 0]])
    e, f, z = tuple(source.module_generators())
    a, b, r = tuple(target.module_generators())
    first = source.primitive_isotropic_subobject(e, z)
    second = target.primitive_isotropic_subobject(b, r)
    witness = INDUCTION.plane_isometry(first, second)
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert all(second.inclusion().is_in_image(witness(point)) for point in (e, z))
    assert all(target.b(witness(x), witness(y)) == source.b(x, y)
               for x in (e, f, z) for y in (e, f, z))


def test_mixed_degenerate_plane_transporter_rejects_distinct_radical_intersection_ranks():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    radical_plane = ambient.primitive_isotropic_subobject(z, w)
    mixed_plane = ambient.primitive_isotropic_subobject(e, z)
    assert INDUCTION.plane_isometry(radical_plane, mixed_plane) is None
    assert isometry(ARROWS, ARROWS.object(radical_plane.inclusion()),
                    ARROWS.object(mixed_plane.inclusion())) is None


def test_moving_primitive_radical_planes_transport_by_unimodular_radical_action():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    first = ambient.primitive_isotropic_subobject(z, w)
    second = ambient.primitive_isotropic_subobject(w, t)
    witness = INDUCTION.plane_isometry(first, second)
    assert witness is not None
    assert all(second.inclusion().is_in_image(witness(point)) for point in (z, w))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, z, w, t) for y in (e, f, z, w, t))
    public = isometry(ARROWS, ARROWS.object(first.inclusion()), ARROWS.object(second.inclusion()))
    assert public is not None
    assert all(second.inclusion().is_in_image(public(point)) for point in (z, w))


def test_mixed_degenerate_primitive_projection_partial_radical_plane_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = ambient.primitive_isotropic_subobject(e + z, w)
    second = ambient.primitive_isotropic_subobject(f + w, z)
    witness = INDUCTION.plane_isometry(first, second)
    assert witness is not None
    assert all(second.inclusion().is_in_image(witness(point)) for point in (e + z, w))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, z, w) for y in (e, f, z, w))
    public = isometry(ARROWS, ARROWS.object(first.inclusion()), ARROWS.object(second.inclusion()))
    assert public is not None
    assert all(second.inclusion().is_in_image(public(point)) for point in (e + z, w))


def test_mixed_degenerate_primitive_projection_partial_radical_plane_stabilizer():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e + z, w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (e + z, w))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, z, w) for y in (e, f, z, w))
    labels = tuple(ambient.module_generating_set())[-2:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w)
    ]) for generator in generators]
    assert MatrixGroup(actions).order() == 2
    public = tuple(automorphism_group_generators(ARROWS, ARROWS.object(plane.inclusion())))
    assert public and all(plane.inclusion().is_in_image(generator(point))
                          for generator in public for point in (e + z, w))


def test_mixed_degenerate_plane_transporter_rejects_distinct_projection_smith_types():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = ambient.primitive_isotropic_subobject(e, w)
    second = ambient.primitive_isotropic_subobject(2 * e + z, w)
    assert INDUCTION.plane_isometry(first, second) is None
    assert isometry(ARROWS, ARROWS.object(first.inclusion()),
                    ARROWS.object(second.inclusion())) is None


def test_mixed_degenerate_rank_three_primitive_projection_and_proper_radical_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    first = ambient.primitive_isotropic_subobject(e + z, u, w)
    second = ambient.primitive_isotropic_subobject(f + w, v, z)
    witness = INDUCTION.plane_isometry(first, second)
    assert witness is not None
    assert all(second.inclusion().is_in_image(witness(point)) for point in (e + z, u, w))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, u, v, z, w) for y in (e, f, u, v, z, w))


def test_mixed_degenerate_rank_three_primitive_projection_partial_radical_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e + z, u, w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (e + z, u, w))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z, w)
               for y in (e, f, u, v, z, w))


def test_mixed_degenerate_primitive_projection_graph_plane_stabilizer_without_radical():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e + z, u + w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (e + z, u + w))
    labels = tuple(ambient.module_generating_set())[-2:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w)
    ]) for generator in generators]
    assert MatrixGroup(actions).order() == 6


def test_mixed_degenerate_primitive_projection_transporter_across_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    source = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    target = Lattices(ZZ)([
        [0, -1, 0, 0], [-1, 0, 0, 0],
        [0, 0, 0, 0], [0, 0, 0, 0],
    ])
    e, f, z, w = tuple(source.module_generators())
    a, b, r, s = tuple(target.module_generators())
    first = source.primitive_isotropic_subobject(e + z, w)
    second = target.primitive_isotropic_subobject(b + r, s)
    witness = INDUCTION.plane_isometry(first, second)
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert all(second.inclusion().is_in_image(witness(point)) for point in (e + z, w))
    assert all(target.b(witness(x), witness(y)) == source.b(x, y)
               for x in (e, f, z, w) for y in (e, f, z, w))


def test_primitive_radical_plane_stabilizer_is_full_radical_parabolic():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, w))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, z, w, t) for y in (e, f, z, w, t))
    public = tuple(automorphism_group_generators(ARROWS, ARROWS.object(plane.inclusion())))
    assert public
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in public for point in (z, w))
    labels = tuple(ambient.module_generating_set())[-3:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w, t)
    ]) for generator in generators]
    assert MatrixGroup(actions).order() == 24


def test_zero_form_rank_three_primitive_plane_transporter_uses_unimodular_completion():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, g = tuple(ambient.module_generators())
    first = ambient.primitive_isotropic_subobject(e, f)
    second = ambient.primitive_isotropic_subobject(f, g)
    witness = INDUCTION.plane_isometry(first, second)
    assert witness is not None
    assert all(second.inclusion().is_in_image(witness(point)) for point in (e, f))
    generators = INDUCTION.plane_automorphism_generators(first)
    assert generators
    assert all(first.inclusion().is_in_image(generator(point))
               for generator in generators for point in (e, f))


def test_full_radical_plane_self_transporter_uses_full_rank_completion():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    radical = ambient.primitive_isotropic_subobject(z, w)
    witness = INDUCTION.plane_isometry(radical, radical)
    assert witness is not None
    assert all(radical.inclusion().is_in_image(witness(point)) for point in (z, w))


def test_isotropic_plane_transporter_rejects_nondegenerate_to_degenerate_ambient():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    nondegenerate = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    degenerate = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v = tuple(nondegenerate.module_generators())
    a, b, z, w = tuple(degenerate.module_generators())
    first = nondegenerate.primitive_isotropic_subobject(e, u)
    second = degenerate.primitive_isotropic_subobject(a, z)
    assert INDUCTION.plane_isometry(first, second) is None
    assert INDUCTION.plane_isometry(second, first) is None


def test_public_moving_full_radical_plane_transporter_between_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    source = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    target = Lattices(ZZ)([[0, -1, 0], [-1, 0, 0], [0, 0, 0]])
    e, f, z = tuple(source.module_generators())
    a, b, r = tuple(target.module_generators())
    first = ARROWS.object(source.primitive_isotropic_subobject(e, z).inclusion())
    second = ARROWS.object(target.primitive_isotropic_subobject(b, r).inclusion())
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert all(second.arrow().is_in_image(witness(point)) for point in (e, z))


def test_full_radical_maximal_isotropic_stabilizer_has_full_levi_and_shear_images():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e, z, w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    labels = tuple(ambient.module_generating_set())[-2:]
    radical_actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w)
    ]) for generator in generators]
    assert MatrixGroup(radical_actions).order() == 6
    assert any(generator(e) == e + z for generator in generators)
    assert any(generator(f) == f + w for generator in generators)
    assert any(generator(e) == -e for generator in generators)


def test_full_radical_rank_three_isotropic_stabilizer_over_two_hyperbolic_planes():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, u, v, z = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e, u, z)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (e, u, z))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z) for y in (e, f, u, v, z))
    assert any(generator(e) == e + z for generator in generators)
    assert any(generator(u) == u + z for generator in generators)


def test_full_radical_nonmaximal_isotropic_plane_stabilizer_over_three_u():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic + Lattices(ZZ)([[0]])
    e, f, u, v, a, b, z = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e, u, z)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (e, u, z))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, a, b, z)
               for y in (e, f, u, v, a, b, z))


def test_mixed_degenerate_radical_to_isotropic_plane_flag_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, INDUCTION, automorphism_group_generators, _flag_marking_from_basis,
    )

    for k in (1, 2):
        ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * k for _ in range(k)])
        e, f, *radical = tuple(ambient.module_generators())
        basis = (*radical, e)
        flag = FLAGS(_flag_marking_from_basis(ambient, basis, (k, k + 1)))
        generators = tuple(automorphism_group_generators(FLAGS, flag))
        expected = INDUCTION.plane_automorphism_generators(ambient.primitive_isotropic_subobject(*basis))
        assert generators == expected
        assert all(ambient.primitive_isotropic_subobject(*basis).inclusion().is_in_image(generator(point))
                   for generator in generators for point in basis)


def test_proper_radical_intersection_first_flag_stabilizer_is_top_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, INDUCTION, automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    for top_vector in (e, e + w, e + 3 * w):
        flag = FLAGS(_flag_marking_from_basis(ambient, (z, top_vector), (1, 2)))
        generators = tuple(automorphism_group_generators(FLAGS, flag))
        expected = INDUCTION.plane_automorphism_generators(
            ambient.primitive_isotropic_subobject(z, top_vector)
        )
        assert generators and generators == expected
        assert all(ambient.primitive_isotropic_subobject(z).inclusion().is_in_image(generator(z))
                   for generator in generators)


def test_proper_radical_intersection_first_flag_transporter_uses_top_subspace():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = FLAGS(_flag_marking_from_basis(ambient, (z, e + w), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(ambient, (w, f + z), (1, 2)))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    assert ambient.primitive_isotropic_subobject(w).inclusion().is_in_image(witness(z))
    top = ambient.primitive_isotropic_subobject(w, f + z)
    assert all(top.inclusion().is_in_image(witness(point)) for point in (z, e + w))


def test_proper_radical_intersection_flag_transporter_rejects_projection_smith_content():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, INDUCTION, isometry, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first_plane = ambient.primitive_isotropic_subobject(z, e + w)
    second_plane = ambient.primitive_isotropic_subobject(z, 2 * e + w)
    assert INDUCTION.plane_isometry(first_plane, second_plane) is None
    first = FLAGS(_flag_marking_from_basis(ambient, (z, e + w), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(ambient, (z, 2 * e + w), (1, 2)))
    assert isometry(FLAGS, first, second) is None


def test_nonprimitive_projection_proper_radical_flag_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = FLAGS(_flag_marking_from_basis(ambient, (z, 2 * e + w), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(ambient, (w, 2 * f + z), (1, 2)))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    assert ambient.primitive_isotropic_subobject(w).inclusion().is_in_image(witness(z))
    target = ambient.primitive_isotropic_subobject(w, 2 * f + z)
    assert all(target.inclusion().is_in_image(witness(point)) for point in (z, 2 * e + w))


def test_nonprimitive_projection_radical_rank_three_quotient_line_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    for c in (2, 3, 5):
        source = ambient.primitive_isotropic_subobject(z, c * e + w)
        target = ambient.primitive_isotropic_subobject(t, c * f + w + z)
        witness = INDUCTION.plane_isometry(source, target)
        assert witness is not None
        assert all(target.inclusion().is_in_image(witness(point))
                   for point in (z, c * e + w))
        public = isometry(ARROWS, ARROWS.object(source.inclusion()),
                          ARROWS.object(target.inclusion()))
        assert public is not None


def test_nonprimitive_projection_rank_two_radical_intersection_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    source = ambient.primitive_isotropic_subobject(z, w, 3 * e + t)
    target = ambient.primitive_isotropic_subobject(w, t, 3 * f + z)
    witness = INDUCTION.plane_isometry(source, target)
    assert witness is not None
    assert all(target.inclusion().is_in_image(witness(point))
               for point in (z, w, 3 * e + t))


def test_nonprimitive_projection_rank_two_radical_intersection_flag_public_api():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, automorphism_group_generators, isometry, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    first = FLAGS(_flag_marking_from_basis(ambient, (z, w, 3 * e + t), (2, 3)))
    second = FLAGS(_flag_marking_from_basis(ambient, (w, t, 3 * f + z), (2, 3)))
    generators = tuple(automorphism_group_generators(FLAGS, first))
    assert generators
    radical = ambient.primitive_isotropic_subobject(z, w)
    top = ambient.primitive_isotropic_subobject(z, w, 3 * e + t)
    assert all(radical.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, w))
    assert all(top.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, w, 3 * e + t))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    target = ambient.primitive_isotropic_subobject(w, t, 3 * f + z)
    assert all(target.inclusion().is_in_image(witness(point))
               for point in (z, w, 3 * e + t))


def test_nonprimitive_projection_proper_radical_flag_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, automorphism_group_generators, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    for c in (2, 3, 5):
        flag = FLAGS(_flag_marking_from_basis(ambient, (z, c * e + w), (1, 2)))
        generators = tuple(automorphism_group_generators(FLAGS, flag))
        assert generators
        first = ambient.primitive_isotropic_subobject(z)
        top = ambient.primitive_isotropic_subobject(z, c * e + w)
        assert all(first.inclusion().is_in_image(generator(z)) for generator in generators)
        assert all(top.inclusion().is_in_image(generator(point))
                   for generator in generators for point in (z, c * e + w))
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for generator in generators for x in (e, f, z, w) for y in (e, f, z, w))
        assert any(generator(z) == -z for generator in generators)


def test_nonprimitive_projection_parabolic_contains_full_integral_shear_kernel():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, 5 * e + w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert any(generator(e) == e + z for generator in generators)
    assert any(generator(f) == f + z for generator in generators)
    assert any(generator(f) == f + w for generator in generators)
    assert any(generator(w) == w + z for generator in generators)
    assert any(generator(z) == -z for generator in generators)


def test_nonprimitive_projection_rank_three_radical_saturation_stabilizer():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    for basis in ((z, 3 * e + w), (z, w, 3 * e + t)):
        plane = ambient.primitive_isotropic_subobject(*basis)
        generators = INDUCTION.plane_automorphism_generators(plane)
        assert generators
        assert all(plane.inclusion().is_in_image(generator(point))
                   for generator in generators for point in basis)
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for generator in generators for x in (e, f, z, w, t)
                   for y in (e, f, z, w, t))
        if len(basis) == 3:
            labels = tuple(ambient.module_generating_set())[-3:-1]
            actions = [matrix(GF(2), [
                [int(generator(point).to_vector()(label)) for label in labels]
                for point in (z, w)
            ]) for generator in generators]
            assert MatrixGroup(actions).order() == 6


def test_nonprimitive_projection_proper_radical_stabilizer_over_two_u():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, 3 * e + w)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, 3 * e + w))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z, w)
               for y in (e, f, u, v, z, w))


def test_nonprimitive_rank_two_projection_mixed_isotropic_subspace_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, 2 * e + w, u)
    generators = INDUCTION.plane_automorphism_generators(plane)
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, 2 * e + w, u))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z, w)
               for y in (e, f, u, v, z, w))
    # The returned group is not required to contain a particular shear as
    # one of its generators: Schreier reduction changes generating families.
    assert any(generator(z) != z for generator in generators)


def test_nonprimitive_rank_two_projection_graph_without_radical_intersection():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    graph = ambient.primitive_isotropic_subobject(2 * e + z, u + w)
    assert graph.inclusion().is_primitive()
    generators = INDUCTION.plane_automorphism_generators(graph)
    assert generators
    assert all(graph.inclusion().is_in_image(generator(point))
               for generator in generators for point in (2 * e + z, u + w))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z, w)
               for y in (e, f, u, v, z, w))


def test_nonprimitive_rank_two_projection_radical_first_flag_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        INDUCTION, FLAGS, automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, 2 * e + w, u)
    generators = INDUCTION.flag_automorphism_generators((z, plane))
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, 2 * e + w, u))
    assert all(generator(z) in (z, -z) for generator in generators)
    flag = FLAGS(_flag_marking_from_basis(ambient, (z, 2 * e + w, u), (1, 3)))
    public = tuple(automorphism_group_generators(FLAGS, flag))
    assert public
    assert all(generator(z) in (z, -z) for generator in public)
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in public for point in (z, 2 * e + w, u))


def test_nonprimitive_rank_two_projection_rank_two_radical_first_flag():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    flag = FLAGS(_flag_marking_from_basis(ambient, (z, w, 2 * e + t, u), (2, 4)))
    generators = tuple(automorphism_group_generators(FLAGS, flag))
    assert generators
    radical = ambient.primitive_isotropic_subobject(z, w)
    top = ambient.primitive_isotropic_subobject(z, w, 2 * e + t, u)
    assert all(radical.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, w))
    assert all(top.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, w, 2 * e + t, u))


def test_nonprimitive_projection_nested_radical_flag_public_stabilizer():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    flag = category(_flag_marking_from_basis(ambient, (z, w, 2 * e + t, u), (1, 2, 4)))
    generators = tuple(automorphism_group_generators(category, flag))
    assert generators
    for rank in (1, 2, 4):
        points = (z, w, 2 * e + t, u)[:rank]
        term = ambient.primitive_isotropic_subobject(*points)
        assert all(term.inclusion().is_in_image(generator(point))
                   for generator in generators for point in points)
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in ambient.module_generators()
               for y in ambient.module_generators())
    # Independently, the image on R/2R is the full upper parabolic
    # preserving the complete radical flag <z> < <z,w> < R.
    # Its three independent off-diagonal entries give order 2^3.
    labels = tuple(ambient.module_generating_set())
    radical_actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(labels[j])) for j in (4, 5, 6)]
        for point in (z, w, t)
    ]) for generator in generators]
    assert MatrixGroup(radical_actions).order() == 8


def test_nonprimitive_projection_nested_radical_flag_public_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    source = category(_flag_marking_from_basis(ambient, (z, w, 2 * e + t, u), (1, 2, 4)))
    target = category(_flag_marking_from_basis(ambient, (w, z, 2 * e + t, u), (1, 2, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank, source_points, target_points in (
        (1, (z,), (w,)),
        (2, (z, w), (w, z)),
        (4, (z, w, 2 * e + t, u), (w, z, 2 * e + t, u)),
    ):
        term = ambient.primitive_isotropic_subobject(*target_points)
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_points)
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())
    # A non-permutation change of the first radical line must also lift.
    shear_target = category(_flag_marking_from_basis(
        ambient, (z + w, w, 2 * e + t, u), (1, 2, 4)))
    shear_witness = isometry(category, source, shear_target)
    assert shear_witness is not None
    assert ambient.primitive_isotropic_subobject(z + w).inclusion().is_in_image(shear_witness(z))
    assert ambient.primitive_isotropic_subobject(z, w).inclusion().is_in_image(shear_witness(w))
    assert all(ambient.b(shear_witness(x), shear_witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())
    # A change in the top adapted basis, without changing its integral
    # sublattice or any intermediate term, must not introduce a false
    # obstruction to transport.
    top_basis_target = category(_flag_marking_from_basis(
        ambient, (w, z, 2 * e + t + z, u + w), (1, 2, 4)))
    top_basis_witness = isometry(category, source, top_basis_target)
    assert top_basis_witness is not None
    for rank in (1, 2, 4):
        target_term = ambient.primitive_isotropic_subobject(
            *(w, z, 2 * e + t + z, u + w)[:rank])
        assert all(target_term.inclusion().is_in_image(top_basis_witness(point))
                   for point in (z, w, 2 * e + t, u)[:rank])


def test_projected_finite_index_isotropic_transporter_and_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        _finite_index_sublattice_transporter,
        _finite_index_sublattice_stabilizer_generators,
    )

    complement = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(complement.module_generators())
    source = complement.subobject_on((2 * e, u)).inclusion()
    target = complement.subobject_on((2 * e, u + e)).inclusion()
    witness = _finite_index_sublattice_transporter(source, target)
    assert witness is not None
    assert all(target.is_in_image(witness(source(point)))
               for point in source.domain().module_generators())
    generators = tuple(_finite_index_sublattice_stabilizer_generators(target))
    assert generators
    assert all(target.is_in_image(generator(target(point)))
               for generator in generators for point in target.domain().module_generators())


def test_nonprimitive_projection_nested_radical_flag_moving_top_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    source_basis = (z, w, 2 * e + t, u)
    target_basis = (z, w, 2 * e + t, u + e)
    source = category(_flag_marking_from_basis(ambient, source_basis, (1, 2, 4)))
    target = category(_flag_marking_from_basis(ambient, target_basis, (1, 2, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank in (1, 2, 4):
        term = ambient.primitive_isotropic_subobject(*target_basis[:rank])
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_basis[:rank])


def test_nonprimitive_projection_moving_top_and_radical_flag_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    source_basis = (z, w, 2 * e + t, u)
    target_basis = (w, z, 2 * e + t, u + e)
    source = category(_flag_marking_from_basis(ambient, source_basis, (1, 2, 4)))
    target = category(_flag_marking_from_basis(ambient, target_basis, (1, 2, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank in (1, 2, 4):
        term = ambient.primitive_isotropic_subobject(*target_basis[:rank])
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_basis[:rank])
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())


    # The characteristic radical intersection itself may move inside R;
    # it is not an invariant submodule, only its rank is invariant.
    target_basis = (z, t, 2 * e + w, u + e)
    target = category(_flag_marking_from_basis(ambient, target_basis, (1, 2, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank in (1, 2, 4):
        term = ambient.primitive_isotropic_subobject(*target_basis[:rank])
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_basis[:rank])
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())
    source_basis = (z, w, 4 * e + t, u)
    target_basis = (z, w, 4 * e + t, u + e)
    source = category(_flag_marking_from_basis(ambient, source_basis, (1, 2, 4)))
    target = category(_flag_marking_from_basis(ambient, target_basis, (1, 2, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank in (1, 2, 4):
        term = ambient.primitive_isotropic_subobject(*target_basis[:rank])
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_basis[:rank])
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())
    # The same graph-shear construction must work for a non-binary
    # projection index; a single arbitrary projected transporter need not
    # satisfy the graph congruence, so its finite stabilizer coset is used.
    source_basis = (z, w, 3 * e + t, u)
    target_basis = (z, w, 3 * e + t, u + e)
    source = category(_flag_marking_from_basis(ambient, source_basis, (1, 2, 4)))
    target = category(_flag_marking_from_basis(ambient, target_basis, (1, 2, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank in (1, 2, 4):
        term = ambient.primitive_isotropic_subobject(*target_basis[:rank])
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_basis[:rank])


def test_rank_three_radical_intersection_complete_flag_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(5), Lattices(ZZ))
    source_basis = (z, w, t, 2 * e + u)
    target_basis = (w, z, t, 2 * e + u)
    source = category(_flag_marking_from_basis(ambient, source_basis, (1, 2, 3, 4)))
    target = category(_flag_marking_from_basis(ambient, target_basis, (1, 2, 3, 4)))
    witness = isometry(category, source, target)
    assert witness is not None
    for rank in (1, 2, 3, 4):
        term = ambient.primitive_isotropic_subobject(*target_basis[:rank])
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_basis[:rank])
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in ambient.module_generators() for y in ambient.module_generators())


def test_mixed_three_stage_flag_action_data_uses_sage_integer_matrices():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, 2 * e + w, u)
    domain = plane.inclusion().domain()
    first = domain.primitive_isotropic_subobject(plane.inclusion().lift(z))
    second = domain.primitive_isotropic_subobject(
        plane.inclusion().lift(z), plane.inclusion().lift(2 * e + w))
    basis, images = INDUCTION.flag_action_data(plane, (first, second))
    assert len(basis) == 3
    assert images
    assert all(len(image) == 3 for image in images)


def test_mixed_flag_transporter_rejects_intermediate_projection_smith_difference():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    first = category(_flag_marking_from_basis(ambient, (z, e + w, u), (1, 2, 3)))
    second = category(_flag_marking_from_basis(ambient, (z, 2 * e + w, u), (1, 2, 3)))
    assert isometry(category, first, second) is None


def test_induction_legacy_flag_transporter_with_proper_radical_line():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = (z, ambient.primitive_isotropic_subobject(z, e + w))
    second = (w, ambient.primitive_isotropic_subobject(w, f + z))
    witness = INDUCTION.flag_isometry(first, second)
    assert witness is not None
    assert witness(z) in (w, -w)
    assert all(second[1].inclusion().is_in_image(witness(point)) for point in (z, e + w))


def test_induction_legacy_flag_transporter_with_nonprimitive_projection():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = (z, ambient.primitive_isotropic_subobject(z, 5 * e + w))
    second = (w, ambient.primitive_isotropic_subobject(w, 5 * f + z))
    witness = INDUCTION.flag_isometry(first, second)
    assert witness is not None
    assert witness(z) in (w, -w)
    assert all(second[1].inclusion().is_in_image(witness(point))
               for point in (z, 5 * e + w))


def test_induction_legacy_flag_stabilizer_with_proper_radical_line():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, e + w)
    generators = INDUCTION.flag_automorphism_generators((z, plane))
    assert generators
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, e + w))
    assert all(ambient.primitive_isotropic_subobject(z).inclusion().is_in_image(generator(z))
               for generator in generators)


def test_proper_radical_intersection_three_stage_flag_stabilizer():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    flag = category(_flag_marking_from_basis(ambient, (z, e + w, u), (1, 2, 3)))
    generators = tuple(automorphism_group_generators(category, flag))
    assert generators
    for points in ((z,), (z, e + w), (z, e + w, u)):
        term = ambient.primitive_isotropic_subobject(*points)
        assert all(term.inclusion().is_in_image(generator(point))
                   for generator in generators for point in points)
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z, w)
               for y in (e, f, u, v, z, w))
    labels = tuple(ambient.module_generating_set())
    quotient_actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(labels[j])) for j in (0, 2)]
        for point in (e + w, u)
    ]) for generator in generators]
    assert MatrixGroup(quotient_actions).order() == 2
    radical_actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(labels[j])) for j in (4, 5)]
        for point in (z, w)
    ]) for generator in generators]
    assert MatrixGroup(radical_actions).order() == 2


def test_proper_radical_intersection_three_stage_flag_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    first = category(_flag_marking_from_basis(ambient, (z, e + w, u), (1, 2, 3)))
    second = category(_flag_marking_from_basis(ambient, (w, v + z, f), (1, 2, 3)))
    witness = isometry(category, first, second)
    assert witness is not None
    for left, right in (
        ((z,), (w,)), ((z, e + w), (w, v + z)),
        ((z, e + w, u), (w, v + z, f)),
    ):
        target = ambient.primitive_isotropic_subobject(*right)
        assert all(target.inclusion().is_in_image(witness(point)) for point in left)


def test_proper_radical_intersection_four_stage_flag_stabilizer_and_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        automorphism_group_generators, isometry, _flag_marking_from_basis,
    )

    ambient = (Lattices(ZZ)("U") + Lattices(ZZ)("U")
               + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]]))
    e, f, u, v, a, b, z, w = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(5), Lattices(ZZ))
    first = category(_flag_marking_from_basis(ambient, (z, e + w, u, a), (1, 2, 3, 4)))
    second = category(_flag_marking_from_basis(ambient, (w, a + z, e, u), (1, 2, 3, 4)))
    generators = tuple(automorphism_group_generators(category, first))
    assert generators
    for rank in (1, 2, 3, 4):
        term = ambient.primitive_isotropic_subobject(*(z, e + w, u, a)[:rank])
        assert all(term.inclusion().is_in_image(generator(point))
                   for generator in generators for point in (z, e + w, u, a)[:rank])
    witness = isometry(category, first, second)
    assert witness is not None
    for rank in (1, 2, 3, 4):
        target = ambient.primitive_isotropic_subobject(*(w, a + z, e, u)[:rank])
        assert all(target.inclusion().is_in_image(witness(point))
                   for point in (z, e + w, u, a)[:rank])


def test_rank_two_proper_radical_intersection_partial_flag_parabolic_and_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        automorphism_group_generators, isometry, _flag_marking_from_basis,
    )

    ambient = (Lattices(ZZ)("U") + Lattices(ZZ)("U")
               + Lattices(ZZ)([[0] * 3 for _ in range(3)]))
    e, f, u, v, z, w, t = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    first = category(_flag_marking_from_basis(ambient, (z, w, e + t, u), (2, 3, 4)))
    second = category(_flag_marking_from_basis(ambient, (w, t, u + z, e), (2, 3, 4)))
    generators = tuple(automorphism_group_generators(category, first))
    assert generators
    for rank in (2, 3, 4):
        term = ambient.primitive_isotropic_subobject(*(z, w, e + t, u)[:rank])
        assert all(term.inclusion().is_in_image(generator(point))
                   for generator in generators for point in (z, w, e + t, u)[:rank])
    witness = isometry(category, first, second)
    assert witness is not None
    for rank in (2, 3, 4):
        term = ambient.primitive_isotropic_subobject(*(w, t, u + z, e)[:rank])
        assert all(term.inclusion().is_in_image(witness(point))
                   for point in (z, w, e + t, u)[:rank])


def test_provider_proper_radical_intersection_flag_stabilizer_matches_marked_api():
    from sage_extensions.all import Cat, IsotropicFlag, FlagType, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        INDUCTION, automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, u, v, z, w = tuple(ambient.module_generators())
    basis = (z, e + w, u)
    provider_flag = IsotropicFlag(ambient, basis, flag_type=FlagType((1, 2, 3)))
    provider = tuple(INDUCTION.flag_automorphism_generators(provider_flag))
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    marked = category(_flag_marking_from_basis(ambient, basis, (1, 2, 3)))
    public = tuple(automorphism_group_generators(category, marked))
    assert provider and provider == public


def test_provider_proper_radical_four_stage_flag_stabilizer_matches_marked_api():
    from sage_extensions.all import Cat, IsotropicFlag, FlagType, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        INDUCTION, automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = (Lattices(ZZ)("U") + Lattices(ZZ)("U")
               + Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]]))
    e, f, u, v, a, b, z, w = tuple(ambient.module_generators())
    basis = (z, e + w, u, a)
    provider_flag = IsotropicFlag(ambient, basis, flag_type=FlagType((1, 2, 3, 4)))
    provider = tuple(INDUCTION.flag_automorphism_generators(provider_flag))
    category = Cat().Mor(FiniteOrdinalCategory(5), Lattices(ZZ))
    marked = category(_flag_marking_from_basis(ambient, basis, (1, 2, 3, 4)))
    public = tuple(automorphism_group_generators(category, marked))
    assert provider and provider == public


def test_mixed_degenerate_proper_radical_partial_flag_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    first = FLAGS(_flag_marking_from_basis(ambient, (z, w), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(ambient, (t, w), (1, 2)))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    assert ambient.primitive_isotropic_subobject(t).inclusion().is_in_image(witness(z))
    target_plane = ambient.primitive_isotropic_subobject(t, w)
    assert all(target_plane.inclusion().is_in_image(witness(point)) for point in (z, w))


def test_mixed_degenerate_flag_transporter_rejects_radical_intersection_profile():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    first = FLAGS(_flag_marking_from_basis(ambient, (z, e), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(ambient, (e, z), (1, 2)))
    assert isometry(FLAGS, first, second) is None


def test_mixed_degenerate_radical_flag_transporter_between_integral_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    source = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    target = Lattices(ZZ)([
        [0, -1, 0, 0, 0], [-1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0],
    ])
    e, f, z, w, t = tuple(source.module_generators())
    a, b, r, s, q = tuple(target.module_generators())
    first = FLAGS(_flag_marking_from_basis(source, (z, w), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(target, (q, s), (1, 2)))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert target.primitive_isotropic_subobject(q).inclusion().is_in_image(witness(z))
    plane = target.primitive_isotropic_subobject(q, s)
    assert all(plane.inclusion().is_in_image(witness(point)) for point in (z, w))


def test_mixed_degenerate_proper_radical_three_stage_flag_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 4 for _ in range(4)])
    e, f, z, w, t, s = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    first = category(_flag_marking_from_basis(ambient, (z, w, t), (1, 2, 3)))
    second = category(_flag_marking_from_basis(ambient, (s, t, w), (1, 2, 3)))
    witness = isometry(category, first, second)
    assert witness is not None
    for source_points, target_points in (
        ((z,), (s,)), ((z, w), (s, t)), ((z, w, t), (s, t, w))
    ):
        term = ambient.primitive_isotropic_subobject(*target_points)
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_points)


def test_provider_pure_radical_three_stage_flag_stabilizer_is_full_parabolic():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import IsotropicFlag, FlagType, Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 4 for _ in range(4)])
    e, f, z, w, t, s = tuple(ambient.module_generators())
    flag = IsotropicFlag(ambient, (z, w, t), flag_type=FlagType((1, 2, 3)))
    generators = INDUCTION.flag_automorphism_generators(flag)
    assert generators
    for points in ((z,), (z, w), (z, w, t)):
        term = ambient.primitive_isotropic_subobject(*points)
        assert all(term.inclusion().is_in_image(generator(point))
                   for generator in generators for point in points)
    labels = tuple(ambient.module_generating_set())[-4:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w, t, s)
    ]) for generator in generators]
    assert MatrixGroup(actions).order() == 64


def test_legacy_pure_radical_line_plane_flag_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, w)
    generators = INDUCTION.flag_automorphism_generators((z, plane))
    assert generators
    assert all(generator(z) in (z, -z) for generator in generators)
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, w))


def test_pure_radical_flag_pointwise_stabilizer_fixes_top_frame_exactly():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import IsotropicFlag, FlagType, Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    flag = IsotropicFlag(ambient, (z, w), flag_type=FlagType((1, 2)))
    generators = INDUCTION.flag_pointwise_automorphism_generators(flag)
    assert generators
    assert all(generator(z) == z and generator(w) == w for generator in generators)
    labels = tuple(ambient.module_generating_set())[-3:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w, t)
    ]) for generator in generators]
    assert MatrixGroup(actions).order() == 4


def test_full_radical_flag_pointwise_stabilizer_has_unrestricted_shears():
    from sage_extensions.all import IsotropicFlag, FlagType, Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    flag = IsotropicFlag(ambient, (z, w), flag_type=FlagType((1, 2)))
    generators = INDUCTION.flag_pointwise_automorphism_generators(flag)
    assert generators
    assert all(generator(z) == z and generator(w) == w for generator in generators)
    assert any(generator(e) == e + z for generator in generators)
    assert any(generator(f) == f + w for generator in generators)


def test_full_radical_quotient_plane_pointwise_flag_stabilizer_fixes_top_frame():
    from sage_extensions.all import IsotropicFlag, FlagType, Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, u, v, z = tuple(ambient.module_generators())
    flag = IsotropicFlag(ambient, (z, e, u), flag_type=FlagType((2, 3)))
    generators = INDUCTION.flag_pointwise_automorphism_generators(flag)
    assert generators
    assert all(generator(z) == z and generator(e) == e and generator(u) == u
               for generator in generators)
    assert any(generator(f) == f + z for generator in generators)
    assert any(generator(v) == v + z for generator in generators)
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z)
               for y in (e, f, u, v, z))


def test_full_radical_rank_two_pointwise_flag_stabilizer_over_u():
    from sage_extensions.all import IsotropicFlag, FlagType, Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    flag = IsotropicFlag(ambient, (z, e), flag_type=FlagType((1, 2)))
    generators = INDUCTION.flag_pointwise_automorphism_generators(flag)
    assert generators
    assert all(generator(e) == e and generator(z) == z for generator in generators)
    assert any(generator(f) == f + z for generator in generators)


def test_mixed_degenerate_full_radical_partial_quotient_flag_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, automorphism_group_generators, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, u, v, z = tuple(ambient.module_generators())
    flag = FLAGS(_flag_marking_from_basis(ambient, (z, e, u), (2, 3)))
    generators = tuple(automorphism_group_generators(FLAGS, flag))
    assert generators
    first = ambient.primitive_isotropic_subobject(z, e)
    top = ambient.primitive_isotropic_subobject(z, e, u)
    assert all(first.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, e))
    assert all(top.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, e, u))
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, u, v, z) for y in (e, f, u, v, z))


def test_provider_isotropic_flag_full_radical_stabilizer_uses_quotient_parabolic():
    from sage_extensions.all import IsotropicFlag, FlagType, Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, u, v, z = tuple(ambient.module_generators())
    flag = IsotropicFlag(ambient, (z, e, u), flag_type=FlagType((2, 3)))
    generators = INDUCTION.flag_automorphism_generators(flag)
    assert generators
    first = ambient.primitive_isotropic_subobject(z, e)
    top = ambient.primitive_isotropic_subobject(z, e, u)
    assert all(first.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, e))
    assert all(top.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, e, u))


def test_legacy_mixed_radical_first_flag_stabilizer_uses_top_parabolic():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(z, e)
    generators = INDUCTION.flag_automorphism_generators((z, plane))
    assert generators
    assert all(generator(z) in (z, -z) for generator in generators)
    assert all(plane.inclusion().is_in_image(generator(point))
               for generator in generators for point in (z, e))


def test_legacy_mixed_radical_first_flag_transporter_across_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    source = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    target = Lattices(ZZ)([[0, -1, 0], [-1, 0, 0], [0, 0, 0]])
    e, f, z = tuple(source.module_generators())
    a, b, r = tuple(target.module_generators())
    first = source.primitive_isotropic_subobject(z, e)
    second = target.primitive_isotropic_subobject(r, b)
    witness = INDUCTION.flag_isometry((z, first), (r, second))
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert target.primitive_isotropic_subobject(r).inclusion().is_in_image(witness(z))
    assert all(second.inclusion().is_in_image(witness(point)) for point in (z, e))


def test_legacy_mixed_flag_transporter_rejects_radical_to_nonradical_line():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e, z)
    assert INDUCTION.flag_isometry((z, plane), (e, plane)) is None


def test_mixed_degenerate_full_radical_partial_quotient_flag_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        FLAGS, isometry, _flag_marking_from_basis,
    )

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, u, v, z = tuple(ambient.module_generators())
    first = FLAGS(_flag_marking_from_basis(ambient, (z, e, u), (2, 3)))
    second = FLAGS(_flag_marking_from_basis(ambient, (z, u, e), (2, 3)))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    first_target = ambient.primitive_isotropic_subobject(z, u)
    top_target = ambient.primitive_isotropic_subobject(z, u, e)
    assert all(first_target.inclusion().is_in_image(witness(point)) for point in (z, e))
    assert all(top_target.inclusion().is_in_image(witness(point)) for point in (z, e, u))


def test_full_radical_quotient_flag_transporter_rejects_distinct_line_divisibilities():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    ambient = Lattices(ZZ)([
        [0, 2, 0, 0, 0], [2, 0, 0, 0, 0],
        [0, 0, 0, 3, 0], [0, 0, 3, 0, 0],
        [0, 0, 0, 0, 0],
    ])
    e, f, u, v, z = tuple(ambient.module_generators())
    assert e.divisor() != u.divisor()
    first = FLAGS(_flag_marking_from_basis(ambient, (z, e, u), (2, 3)))
    second = FLAGS(_flag_marking_from_basis(ambient, (z, u, e), (2, 3)))
    assert isometry(FLAGS, first, second) is None


def test_mixed_degenerate_radical_first_flag_transporter_between_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, isometry, _flag_marking_from_basis

    source = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    target = Lattices(ZZ)([[0, -1, 0], [-1, 0, 0], [0, 0, 0]])
    e, f, z = tuple(source.module_generators())
    a, b, r = tuple(target.module_generators())
    first = FLAGS(_flag_marking_from_basis(source, (z, e), (1, 2)))
    second = FLAGS(_flag_marking_from_basis(target, (r, b), (1, 2)))
    witness = isometry(FLAGS, first, second)
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert target.primitive_isotropic_subobject(r).inclusion().is_in_image(witness(z))
    plane = target.primitive_isotropic_subobject(r, b)
    assert all(plane.inclusion().is_in_image(witness(point)) for point in (z, e))


def test_mixed_degenerate_full_radical_three_stage_flag_stabilizer():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        automorphism_group_generators, _flag_marking_from_basis,
    )

    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic + Lattices(ZZ)([[0]])
    e, f, u, v, a, b, z = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    flag = category(_flag_marking_from_basis(ambient, (z, e, u, a), (2, 3, 4)))
    generators = tuple(automorphism_group_generators(category, flag))
    assert generators
    for points in ((z, e), (z, e, u), (z, e, u, a)):
        term = ambient.primitive_isotropic_subobject(*points)
        assert all(term.inclusion().is_in_image(generator(point))
                   for generator in generators for point in points)


def test_mixed_degenerate_full_radical_three_stage_flag_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry, _flag_marking_from_basis

    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic + Lattices(ZZ)([[0]])
    e, f, u, v, a, b, z = tuple(ambient.module_generators())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    first = category(_flag_marking_from_basis(ambient, (z, e, u, a), (2, 3, 4)))
    second = category(_flag_marking_from_basis(ambient, (z, a, e, u), (2, 3, 4)))
    witness = isometry(category, first, second)
    assert witness is not None
    for source_points, target_points in (
        ((z, e), (z, a)), ((z, e, u), (z, a, e)), ((z, e, u, a), (z, a, e, u))
    ):
        term = ambient.primitive_isotropic_subobject(*target_points)
        assert all(term.inclusion().is_in_image(witness(point)) for point in source_points)


def test_mixed_degenerate_maximal_moving_rank_three_contains_full_radical():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    source = Lattices(ZZ)([[0] * 3 for _ in range(3)])
    classes = tuple(primitive_marking_classes(ARROWS, source, ambient))
    assert len(classes) == 1
    embedding = classes[0].arrow()
    assert embedding.is_primitive()
    assert int(embedding.domain().rank()) == 3
    radical = ambient.radical().inclusion()
    assert all(embedding.is_in_image(radical(point)) for point in radical.domain().module_generators())


def test_mixed_degenerate_maximal_plane_index_two_uses_integral_radical_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    classes = tuple(nonprimitive_isotropic_sublattice_classes(ambient, 2, (1, 2)))
    # The saturation stabilizer preserves the radical line. Its mod-2
    # image has two orbits on the three index-two sublattices of Z^2.
    assert len(classes) == 2
    assert all(marked.arrow().codomain() is ambient for marked in classes)
    assert all(marked.arrow().saturation().inclusion().domain().rank() == 2 for marked in classes)


def test_mixed_degenerate_maximal_plane_index_two_preserves_requested_so_action():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    full = tuple(nonprimitive_isotropic_sublattice_classes(ambient, 2, (1, 2)))
    special = tuple(nonprimitive_isotropic_sublattice_classes(
        ambient, 2, (1, 2), acting_group=ambient.SO()
    ))
    assert len(full) == 2
    assert len(special) == len(full)
    assert all(marked.arrow().codomain() is ambient for marked in special)


def test_odd_rank_moving_maximal_null_plane_has_identical_o_and_so_orbits():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    full = tuple(primitive_marking_classes(ARROWS, source, ambient))
    special = tuple(primitive_marking_classes(ARROWS, source, ambient, acting_group=ambient.SO()))
    assert len(full) == len(special) == 1


def test_even_rank_mixed_degenerate_maximal_plane_smith_two_so_orbits():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    full = tuple(nonprimitive_isotropic_sublattice_classes(ambient, 3, (1, 1, 2)))
    special = tuple(nonprimitive_isotropic_sublattice_classes(
        ambient, 3, (1, 1, 2), acting_group=ambient.SO()
    ))
    assert len(full) == len(special) == 2


def test_odd_rank_moving_null_line_stream_has_same_so_and_o_orbits():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_isotropic_line_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    full = tuple(islice(iter_all_isotropic_line_classes(ambient), 3))
    special = tuple(islice(iter_all_isotropic_line_classes(ambient, acting_group=ambient.SO()), 3))
    assert len(full) == len(special) == 3
    assert all(left.arrow().codomain() is ambient and right.arrow().codomain() is ambient
               for left, right in zip(full, special, strict=True))
    assert all(left.arrow()(next(iter(left.arrow().domain().module_generators()))) ==
               right.arrow()(next(iter(right.arrow().domain().module_generators())))
               for left, right in zip(full, special, strict=True))


def test_odd_rank_homothetic_maximal_null_plane_so_orbits_equal_full_orbits():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import homothetic_isotropic_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    full = tuple(homothetic_isotropic_sublattice_classes(ambient, 2, 3))
    special = tuple(homothetic_isotropic_sublattice_classes(ambient, 2, 3, acting_group=ambient.SO()))
    assert len(full) == len(special) == 1


def test_odd_rank_moving_nonisotropic_lines_do_not_split_under_so():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonisotropic_line_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    full = tuple(nonisotropic_line_classes(ambient, 2))
    special = tuple(nonisotropic_line_classes(ambient, 2, acting_group=ambient.SO()))
    assert len(full) == len(special) == 1


def test_mixed_degenerate_imprimitive_null_line_stabilizer_equals_primitive():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e = next(iter(ambient.module_generators()))
    primitive = RECURSION.line_stabilizer_generators(e)
    imprimitive = RECURSION.line_stabilizer_generators(3 * e)
    assert primitive == imprimitive
    assert all(generator(3 * e) in (3 * e, -3 * e) for generator in imprimitive)


def test_mixed_degenerate_primitive_projection_line_stabilizer_conjugates_radical_shear():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    for vector in (e, e + z, e + 3 * z):
        generators = RECURSION.line_stabilizer_generators(vector)
        assert generators
        assert all(generator(vector) in (vector, -vector) for generator in generators)
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for generator in generators for x in (e, f, z) for y in (e, f, z))


def test_mixed_degenerate_rank_two_radical_line_stabilizer_preserves_full_form():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    vector = e + 2 * z - 3 * w
    generators = RECURSION.line_stabilizer_generators(vector)
    assert generators
    assert all(generator(vector) in (vector, -vector) for generator in generators)
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in generators for x in (e, f, z, w) for y in (e, f, z, w))
    labels = tuple(ambient.module_generating_set())[-2:]
    radical_images = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w)
    ]) for generator in generators]
    assert MatrixGroup(radical_images).order() == 6


def test_mixed_degenerate_pointwise_null_stabilizer_is_sign_kernel_of_line_stabilizer():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    vector = e + 2 * z - 3 * w
    line = RECURSION.line_stabilizer_generators(vector)
    fixed = RECURSION.pointed_automorphism_generators(vector)
    assert len(line) == len(fixed)
    assert all(generator(vector) == vector for generator in fixed)
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in fixed for x in (e, f, z, w) for y in (e, f, z, w))
    labels = tuple(ambient.module_generating_set())[-2:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w)
    ]) for generator in fixed]
    assert MatrixGroup(actions).order() == 6


def test_mixed_degenerate_scaled_hyperbolic_quotient_line_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)([[0, 2, 0], [2, 0, 0], [0, 0, 0]])
    e, f, z = tuple(ambient.module_generators())
    vector = e + z
    line = RECURSION.line_stabilizer_generators(vector)
    fixed = RECURSION.pointed_automorphism_generators(vector)
    assert line and fixed
    assert all(generator(vector) in (vector, -vector) for generator in line)
    assert all(generator(vector) == vector for generator in fixed)


def test_mixed_degenerate_nonprimitive_projection_primitive_null_line_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    for c in (2, 3, 4):
        vector = c * e + z
        line = RECURSION.line_stabilizer_generators(vector)
        fixed = RECURSION.pointed_automorphism_generators(vector)
        assert line and fixed
        assert all(generator(vector) in (vector, -vector) for generator in line)
        assert all(generator(vector) == vector for generator in fixed)
        assert any(generator(z) == -z for generator in fixed) == (c == 2)
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for generator in line for x in (e, f, z) for y in (e, f, z))


def test_mixed_degenerate_rank_two_radical_nonprimitive_projection_exact_schreier():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    for c, vector, expected_order in (
        (2, 2 * e + z, 2), (3, 3 * e + z + w, 6), (5, 5 * e + z + w, 10)
    ):
        line = RECURSION.line_stabilizer_generators(vector)
        fixed = RECURSION.pointed_automorphism_generators(vector)
        assert line and fixed
        assert all(generator(vector) in (vector, -vector) for generator in line)
        assert all(generator(vector) == vector for generator in fixed)
        assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
                   for generator in line for x in (e, f, z, w) for y in (e, f, z, w))
        labels = tuple(ambient.module_generating_set())[-2:]
        actions = [matrix(GF(c), [
            [int(generator(point).to_vector()(label)) for label in labels]
            for point in (z, w)
        ]) for generator in fixed]
        assert MatrixGroup(actions).order() == expected_order


def test_mixed_degenerate_nonprimitive_projection_line_transporters_respect_residue_orbits():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    for c in (2, 3, 5):
        source = c * e + z
        target = c * e + z + w
        witness = INDUCTION.line_isometry(source, target)
        assert witness is not None
        assert witness(source) in (target, -target)
        assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
                   for x in (e, f, z, w) for y in (e, f, z, w))


def test_mixed_degenerate_rank_three_radical_congruence_stabilizer_image_is_exact():
    from sage.all import GF, MatrixGroup, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, z, w, t = tuple(ambient.module_generators())
    vector = 2 * e + z + w
    fixed = RECURSION.pointed_automorphism_generators(vector)
    assert fixed and all(generator(vector) == vector for generator in fixed)
    labels = tuple(ambient.module_generating_set())[-3:]
    actions = [matrix(GF(2), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w, t)
    ]) for generator in fixed]
    assert MatrixGroup(actions).order() == 24


def test_mixed_degenerate_composite_modulus_radical_stabilizer_image():
    from sage.all import MatrixGroup, Zmod, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    vector = 4 * e + z + w
    fixed = RECURSION.pointed_automorphism_generators(vector)
    assert fixed and all(generator(vector) == vector for generator in fixed)
    labels = tuple(ambient.module_generating_set())[-2:]
    actions = [matrix(Zmod(4), [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in (z, w)
    ]) for generator in fixed]
    assert MatrixGroup(actions).order() == 8


def test_mixed_degenerate_line_stabilizer_on_non_split_integral_presentation():
    from sage.all import ZZ as SageZZ, matrix
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    standard = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    gram = matrix(SageZZ, [[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    change = matrix(SageZZ, [[1, 1, 0], [0, 1, 1], [1, 0, 0]])
    assert abs(change.det()) == 1
    twisted = Lattices(ZZ)([[int(entry) for entry in row]
                            for row in (change * gram * change.transpose()).rows()])
    labels = tuple(twisted.module_generating_set())
    coordinates = matrix(SageZZ, [[2, 0, 1]]) * change.inverse()
    vector = twisted.linear_combination({
        label: twisted.base_ring()(int(coordinates[0, j]))
        for j, label in enumerate(labels) if coordinates[0, j]
    })
    generators = RECURSION.line_stabilizer_generators(vector)
    assert generators
    assert all(generator(vector) in (vector, -vector) for generator in generators)
    assert all(twisted.b(generator(x), generator(y)) == twisted.b(x, y)
               for generator in generators for x in twisted.module_generators()
               for y in twisted.module_generators())


def test_mixed_degenerate_two_hyperbolic_planes_with_radical_line_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + Lattices(ZZ)([[0]])
    e, f, u, v, z = tuple(ambient.module_generators())
    vector = e + u + z
    line = RECURSION.line_stabilizer_generators(vector)
    fixed = RECURSION.pointed_automorphism_generators(vector)
    assert line and fixed
    assert all(generator(vector) in (vector, -vector) for generator in line)
    assert all(generator(vector) == vector for generator in fixed)
    assert all(ambient.b(generator(x), generator(y)) == ambient.b(x, y)
               for generator in line for x in (e, f, u, v, z) for y in (e, f, u, v, z))


def test_public_induction_mixed_degenerate_null_line_stabilizer_uses_exact_split():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import INDUCTION, RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    vector = e + z
    expected = RECURSION.line_stabilizer_generators(vector)
    assert INDUCTION.line_automorphism_generators(ambient, vector) == expected
    assert INDUCTION.line_automorphism_generators(ambient, 3 * vector) == expected


def test_public_moving_arrow_mixed_degenerate_null_line_stabilizer():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    line = ambient.primitive_isotropic_subobject(e + z)
    marked = ARROWS.object(line.inclusion())
    generators = tuple(automorphism_group_generators(ARROWS, marked))
    assert generators
    assert all(generator(e + z) in (e + z, -e - z) for generator in generators)


def test_public_based_mixed_degenerate_null_vector_stabilizer_is_pointwise():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, z = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category(source.hom([e + z], codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert generators
    assert all(generator(e + z) == e + z for generator in generators)


def test_public_mixed_degenerate_projection_content_three_line_and_based_stabilizers():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    vector = 3 * e + z + w
    moving = ARROWS.object(ambient.primitive_isotropic_subobject(vector).inclusion())
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    based = category(source.hom([vector], codomain=ambient))
    moving_generators = tuple(automorphism_group_generators(ARROWS, moving))
    based_generators = tuple(automorphism_group_generators(category, based))
    assert moving_generators and based_generators
    assert all(generator(vector) in (vector, -vector) for generator in moving_generators)
    assert all(generator(vector) == vector for generator in based_generators)


def test_mixed_degenerate_radical_null_line_stabilizer_adds_exact_sign():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION, INDUCTION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    fixed = RECURSION.pointed_automorphism_generators(z)
    line = RECURSION.line_stabilizer_generators(z)
    assert len(line) == len(fixed) + 1
    assert all(generator(z) == z for generator in fixed)
    assert line[-1](z) == -z
    assert all(generator(z) in (z, -z) for generator in line)
    assert RECURSION.line_stabilizer_generators(3 * z) == line
    assert INDUCTION.line_automorphism_generators(ambient, z) == line


def test_zero_form_rank_two_moving_line_stabilizer_preserves_primitive_line():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)([[0, 0], [0, 0]])
    e, f = tuple(ambient.module_generators())
    generators = RECURSION.line_stabilizer_generators(e)
    assert generators
    assert all(generator(e) in (e, -e) for generator in generators)
    assert all(generator.domain() is ambient and generator.codomain() is ambient
               for generator in generators)
    assert all((~generator)(generator(point)) == point
               for generator in generators for point in (e, f))


def test_zero_form_rank_one_line_stabilizer_is_exact_sign_group():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)([[0]])
    e = next(iter(ambient.module_generators()))
    generators = RECURSION.line_stabilizer_generators(e)
    assert generators
    assert any(generator(e) == -e for generator in generators)
    assert all(generator(e) in (e, -e) for generator in generators)


def test_public_moving_radical_line_stabilizer_contains_sign_reversal():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, z, w = tuple(ambient.module_generators())
    marked = ARROWS.object(ambient.primitive_isotropic_subobject(z).inclusion())
    generators = tuple(automorphism_group_generators(ARROWS, marked))
    assert generators
    assert all(generator(z) in (z, -z) for generator in generators)
    assert any(generator(z) == -z for generator in generators)


@pytest.mark.parametrize("smith", ((1, 1, 2), (1, 2, 2)))
def test_mixed_degenerate_maximal_rank_three_index_two_respects_radical_plane(smith):
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_isotropic_sublattice_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    classes = tuple(nonprimitive_isotropic_sublattice_classes(ambient, 3, smith))
    # The stabilizer of the full radical plane in the maximal null rank-three
    # space has two mod-two orbits on lines and on hyperplanes: those
    # distinguished by their incidence with the radical plane.
    assert len(classes) == 2
    assert all(marked.arrow().codomain() is ambient for marked in classes)


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
