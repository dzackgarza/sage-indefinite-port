"""Fixed isotropic bases split according to the actual integral Levi image."""
from sage.all import ZZ as SageZZ, matrix
from sage_indefinite_port.indefinite.marked_lattices import INDUCTION
from tests.acceptance.consumer import lattice
import pytest


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


def test_public_moving_maximal_plane_class_in_two_u():
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    source = lattice([[0, 0], [0, 0]])
    classes = tuple(primitive_marking_classes(ARROWS, source, ambient))
    assert len(classes) == 1
    arrow = classes[0].arrow()
    assert int(arrow.domain().rank()) == 2 and arrow.codomain() is ambient
    images = tuple(arrow(point) for point in arrow.domain().module_generators())
    assert all(ambient.b(left, right) == 0 for left in images for right in images)
    assert arrow.is_primitive()


@pytest.mark.parametrize("exchange", (False, True), ids=("integral-levi", "divisibility-obstruction"))
def test_fixed_plane_transporters_retain_the_basis_divisibilities(exchange):
    r"""In U+2U, e and u have divisibilities one and two.

    Swapping them cannot extend to an integral ambient isometry. The
    change (e,u) -> (e+u,u) does extend: fix f and send v to v-2f.
    Both changes preserve the same unbased primitive null plane.
    """
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS

    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    source = lattice([[0, 0], [0, 0]])
    category = CosliceCategory(Lattices(ZZ), source)
    first_arrow = source.hom([e, u], codomain=ambient)
    target_images = (u, e) if exchange else (e + u, u)
    second_arrow = source.hom(list(target_images), codomain=ambient)
    first, second = category(first_arrow), category(second_arrow)
    homset = category.Core().Mor(first, second)
    assert homset.is_empty() is exchange
    if not exchange:
        witness = homset.an_element()
        apex = witness.forward().right()
        inverse = (~witness).forward().right()
        assert apex(e) == e + u and apex(u) == u
        basis = (e, f, u, v)
        assert all(inverse(apex(point)) == point for point in basis)
        assert all(ambient.b(apex(x), apex(y)) == ambient.b(x, y) for x in basis for y in basis)
    moving_homset = ARROWS.Core().Mor(ARROWS(first_arrow), ARROWS(second_arrow))
    assert moving_homset.is_empty() is False
    moving = moving_homset.an_element()
    forward, backward = moving.forward(), (~moving).forward()
    for point in source.module_generators():
        assert second_arrow(forward.left()(point)) == forward.right()(first_arrow(point))
        assert backward.left()(forward.left()(point)) == point
    assert all(backward.right()(forward.right()(point)) == point for point in (e, f, u, v))
