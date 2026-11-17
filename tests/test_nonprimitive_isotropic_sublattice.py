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
