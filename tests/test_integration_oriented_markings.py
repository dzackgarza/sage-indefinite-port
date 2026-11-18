"""Regression for the integral oriented-marking route retained across integration."""
from sage_extensions.all import Lattices, ZZ
from sage_extensions.categories.isotropic_orbits import OrientedIsotropicFlag
from sage_extensions.categories.lattices import OrientedIsotropicFlags
from sage_indefinite_port.indefinite.marked_lattices import (
    automorphism_group_generators,
    isometry,
    primitive_marking_classes,
)


def test_two_u_oriented_plane_transporter_and_full_integral_stabilizer():
    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U")
    e, f, u, v = tuple(ambient.module_generators())
    category = OrientedIsotropicFlags(lattices, 1)
    positive = category(OrientedIsotropicFlag(ambient, (e, u)))
    negative = category(OrientedIsotropicFlag(ambient, (u, e)))
    plane = ambient.primitive_isotropic_subobject(e, u)
    inclusion = plane.inclusion()
    labels = tuple(inclusion.domain().module_generating_set())

    def determinant(g):
        left, right = (inclusion.lift(g(point)).to_vector() for point in (e, u))
        return left(labels[0]) * right(labels[1]) - left(labels[1]) * right(labels[0])

    witness = isometry(category, positive, negative)
    assert witness is not None
    assert determinant(witness) == -1
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, u, v) for y in (e, f, u, v))
    family = tuple(automorphism_group_generators(category, positive))
    assert family
    assert all(determinant(generator) == 1 for generator in family)
    assert len(primitive_marking_classes(category, lattices([[0, 0], [0, 0]]), ambient)) == 1
