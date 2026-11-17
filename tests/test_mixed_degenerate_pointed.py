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
