"""Fixed isotropic bases split according to the actual integral Levi image."""
from sage.all import ZZ as SageZZ, matrix
from sage_indefinite_port.indefinite.marked_lattices import INDUCTION
from tests.acceptance.consumer import lattice


def test_bases_of_a_maximal_plane_in_u_plus_two_u_have_three_classes():
    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    space = ambient.primitive_isotropic_subobject(basis[0], basis[2])
    inclusion = space.inclusion()
    labels = tuple(inclusion.domain().module_generating_set())
    representatives = tuple(INDUCTION.basis_marking_classes(ambient, (basis[0], basis[2])))
    assert len(representatives) == 3
    divisibility_lines = set()
    for images in representatives:
        coordinates = [inclusion.lift(point).to_vector() for point in images]
        assert abs(matrix(SageZZ, [[int(row(label)) for label in labels] for row in coordinates]).det()) == 1
        assert all(ambient.b(left, right) == 0 for left in images for right in images)
        # The divisor-two vectors form one line in I/2I. GL2(Z) has
        # three such lines, and its stabilizer has index three.
        even = tuple((a, b) for a, b in ((1, 0), (0, 1), (1, 1)) if int((a * images[0] + b * images[1]).divisor()) % 2 == 0)
        assert len(even) == 1
        divisibility_lines.add(even[0])
    assert len(divisibility_lines) == 3


def test_public_fixed_maximal_plane_classes_in_two_u():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    source = lattice([[0, 0], [0, 0]])
    category = CosliceCategory(Lattices(ZZ), source)
    classes = tuple(primitive_marking_classes(category, source, ambient))
    assert len(classes) == 1
    arrow = classes[0].arrow()
    assert arrow.domain() is source and arrow.codomain() is ambient
    images = tuple(arrow(point) for point in source.module_generators())
    assert all(ambient.b(left, right) == 0 for left in images for right in images)
    assert arrow.is_primitive()
