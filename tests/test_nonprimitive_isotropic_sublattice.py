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
    import pytest
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    e, _f, u, _v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    marked = ARROWS(source.hom([2 * e, u], codomain=ambient))
    with pytest.raises(ValueError):
        tuple(automorphism_group_generators(ARROWS, marked))
