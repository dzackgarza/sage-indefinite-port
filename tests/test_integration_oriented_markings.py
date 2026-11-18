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


def test_zero_form_full_flag_shape_from_distinct_ambient_parent():
    """The terminal zero-form ambient is not an additional isotropic term."""
    from sage_extensions.all import Cat
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import (
        _flag_arrows,
        _flag_marking_from_basis,
    )

    lattices = Lattices(ZZ)
    source_ambient = lattices([[0] * 3 for _ in range(3)])
    target = lattices([[0] * 3 for _ in range(3)], module_generators=("x", "y", "z"))
    assert source_ambient is not target
    category = Cat().Mor(FiniteOrdinalCategory(3), lattices)
    source_basis = tuple(source_ambient.module_generators())
    shape = category(_flag_marking_from_basis(source_ambient, source_basis[:2], (1, 2)))
    classes = tuple(primitive_marking_classes(category, shape, target))
    assert len(classes) == 1
    arrows = _flag_arrows(classes[0])
    assert tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2)
    assert all(arrow.is_primitive() for arrow in arrows)
