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


def test_rank_three_fixed_basis_transport_preserves_each_marked_generator():
    """The based rank-three marking requires the identity action on its source."""
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    source = Lattices(ZZ)([[0] * 3 for _ in range(3)])
    e, f, u, v, w, z = tuple(ambient.module_generators())
    category = CosliceCategory(Lattices(ZZ), source)
    first = category(source.hom([e, u, w], codomain=ambient))
    second = category(source.hom([e + u, u, w], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert (witness(e), witness(u), witness(w)) == (e + u, u, w)
    assert all((~witness)(witness(point)) == point for point in (e, f, u, v, w, z))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f, u, v, w, z) for y in (e, f, u, v, w, z))


def test_rank_three_fixed_basis_rejects_divisibility_exchange():
    """A based lift cannot exchange primitive null vectors of divisibilities one and two."""
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0],
                       [0, 0, 0, 2, 0, 0], [0, 0, 2, 0, 0, 0],
                       [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0]])
    e, _f, u, _v, w, _z = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0] * 3 for _ in range(3)])
    category = CosliceCategory(Lattices(ZZ), source)
    first = category(source.hom([e, u, w], codomain=ambient))
    second = category(source.hom([u, e, w], codomain=ambient))
    assert int(e.divisor()) == 1 and int(u.divisor()) == 2
    assert isometry(category, first, second) is None


def test_rank_three_fixed_basis_stabilizer_fixes_each_source_vector():
    """The full based stabilizer fixes each of three isotropic generators pointwise."""
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    e, f, u, v, w, z = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0] * 3 for _ in range(3)])
    category = CosliceCategory(Lattices(ZZ), source)
    marked = category(source.hom([e, u, w], codomain=ambient))
    family = tuple(automorphism_group_generators(category, marked))
    assert all((g(e), g(u), g(w)) == (e, u, w) for g in family)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y) for g in family for x in (e, f, u, v, w, z) for y in (e, f, u, v, w, z))


def test_rank_four_fixed_isotropic_basis_has_integral_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(8)] for i in range(8)])
    basis = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0] * 4 for _ in range(4)])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom(list(basis[::2]), codomain=ambient))
    second = category(source.hom(list(basis[1::2]), codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert all(witness(basis[2 * i]) == basis[2 * i + 1] for i in range(4))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in basis for y in basis)


def test_rank_four_pointwise_stabilizer_generates_all_six_skew_shears():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(8)] for i in range(8)])
    basis = tuple(ambient.module_generators())
    null = basis[::2]
    dual = basis[1::2]
    source = Lattices(ZZ)([[0] * 4 for _ in range(4)])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category(source.hom(list(null), codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    labels = tuple(ambient.module_generating_set())
    parameter_space = ambient.base_ring().free_module(6)
    parameter_basis = tuple(parameter_space.module_generators())
    parameters = []
    for g in generators:
        assert all(g(point) == point for point in null)
        entries = []
        for i in range(4):
            assert g(dual[i]) - dual[i] == sum((ambient.scalar_multiple(g(dual[i]).to_vector()(labels[2 * j]), null[j]) for j in range(4)), ambient.zero())
            for j in range(i + 1, 4):
                a = g(dual[i]).to_vector()(labels[2 * j])
                assert a == -g(dual[j]).to_vector()(labels[2 * i])
                entries.append(a)
        parameters.append(sum((parameter_space.scalar_multiple(a, point) for a, point in zip(entries, parameter_basis, strict=True)), parameter_space.zero()))
    assert parameter_space.subobject_on(parameters).inclusion().cokernel().cardinality() == cardinal(1)
