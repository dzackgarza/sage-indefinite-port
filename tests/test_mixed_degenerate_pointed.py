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
