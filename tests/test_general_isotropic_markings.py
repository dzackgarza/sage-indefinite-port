"""Maximal isotropic sublattices of 3U admit ambient transporters."""
import pytest

from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry
from tests.acceptance.consumer import lattice


def test_extension_shortening_preserves_the_integral_flag_orbit():
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(4)] for i in range(4)])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    e, f, x, y = (ambient.module_generator(label) for label in labels)
    integral_shear = ambient.Isom(ambient)(dict(zip(labels, (e, f + x, x, y - e), strict=True)))
    shear = integral_shear.base_change(ambient.base_ring().fraction_field_map())
    identity = rational.Isom(rational).identity()
    integral_basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    radical, extension = integral_basis[0], integral_basis[3]
    for transporter, point in ((shear * shear * shear, extension), (identity, extension + radical + radical + radical)):
        shortened = IsotropicInduction.shortened_extension(ambient, transporter, point, lambda: (shear,))
        assert shortened == extension
        assert ambient.b(shortened, radical) == 0
        assert shortened.q() == ambient.base_ring().zero()


@pytest.mark.parametrize("dimensions", ((3,), (1, 3), (1, 2, 3)), ids=("plane", "partial", "complete"))
def test_general_flag_invariants_are_preserved_by_hyperbolic_pair_exchange(dimensions):
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    e, f, x, y, u, v = tuple(ambient.module_generators())
    induction = IsotropicInduction(EichlerRecursion())
    source = induction.flag_invariant(ambient, (e + x, x, u), dimensions)
    target = induction.flag_invariant(ambient, (f + y, y, v), dimensions)
    assert source == target
    assert source.dimensions == dimensions
    assert source.ambient.rank == 6
    assert source.quotient.rank == 0
    assert source.quotient.discriminant == 1
    assert source.denominator_primes == ()


def test_glued_line_invariant_retains_the_rational_denominator_prime():
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[0, 2, 0, 0], [2, 0, 1, 0], [0, 1, 0, 2], [0, 0, 2, 0]])
    line = next(iter(ambient.module_generators()))
    invariant = IsotropicInduction(EichlerRecursion()).flag_invariant(ambient, (line,), (1,))
    assert invariant.dimensions == (1,)
    assert invariant.ambient.rank == 4
    assert invariant.quotient.rank == 2
    assert invariant.quotient.discriminant == 4
    assert invariant.denominator_primes == (2,)


def test_public_pointwise_plane_stabilizer_in_two_u_is_the_integral_skew_shear_group():
    r"""Fixing e,u in 2U leaves f -> f+t*u, v -> v-t*e, with t in ZZ.

    The pairing equations force this form. A unit parameter ideal proves
    the returned family generates every integral value of t.
    """
    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    e, f, u, v = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e, u)
    family = tuple(ambient.O().pointwise_stabilizer(plane.inclusion()).group_generators())
    labels = tuple(ambient.module_generating_set())
    parameters = []
    for generator in family:
        parameter = generator(f).to_vector()(labels[2])
        assert generator(e) == e and generator(u) == u
        assert generator(f) == f + parameter * u
        assert generator(v) == v - parameter * e
        parameters.append(parameter)
    assert ambient.base_ring().ideal(parameters) == ambient.base_ring().ideal(1)


def test_public_pointwise_maximal_stabilizer_in_three_u_has_all_three_skew_parameters():
    r"""Fixing e,u,w in 3U leaves the full three-parameter integral skew group."""
    from sage_extensions.all import cardinal

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    e, f, u, v, w, z = tuple(ambient.module_generators())
    subspace = ambient.primitive_isotropic_subobject(e, u, w)
    family = tuple(ambient.O().pointwise_stabilizer(subspace.inclusion()).group_generators())
    labels = tuple(ambient.module_generating_set())
    parameter_space = ambient.base_ring().free_module(3)
    x, y, s = tuple(parameter_space.module_generators())
    parameters = []
    for generator in family:
        first = generator(f).to_vector()
        second = generator(v).to_vector()
        a, b, c = first(labels[2]), first(labels[4]), second(labels[4])
        assert generator(e) == e and generator(u) == u and generator(w) == w
        assert generator(f) == f + a * u + b * w
        assert generator(v) == v - a * e + c * w
        assert generator(z) == z - b * e - c * u
        parameters.append(a * x + b * y + c * s)
    assert parameter_space.subobject_on(parameters).inclusion().cokernel().cardinality() == cardinal(1)


def test_isotropic_orbit_bases_retain_the_original_lattice_presentation():
    r"""Every primitive maximal isotropic plane in 2U has an integral hyperbolic completion.

    Thus O(2U) has one such orbit. The first block here is U written in
    the basis (e+f,f), so reduction must transport its result back.
    """
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[2, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    bases = IsotropicInduction(EichlerRecursion()).isotropic_flag_bases(ambient, (2,))
    assert len(bases) == 1
    basis = bases[0]
    assert len(basis) == 2
    assert all(point.parent() is ambient for point in basis)
    assert all(ambient.b(left, right) == 0 for left in basis for right in basis)
    subobject = ambient.subobject_on(basis)
    assert subobject.inclusion().domain().rank() == 2
    assert subobject.inclusion().is_primitive()


def test_three_u_maximal_isotropic_sublattices_have_a_transporter():
    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    ambient = lattice(gram)
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    source = ambient.primitive_isotropic_subobject(*basis[::2])
    target = ambient.primitive_isotropic_subobject(*basis[1::2])
    witness = isometry(ARROWS, ARROWS.object(source.inclusion()), ARROWS.object(target.inclusion()))
    assert witness is not None
    for point in basis[::2]:
        assert target.inclusion().is_in_image(witness(point))
    for point in basis[1::2]:
        assert source.inclusion().is_in_image((~witness)(point))
    for left in basis:
        for right in basis:
            assert ambient.b(witness(left), witness(right)) == ambient.b(left, right)


def test_three_u_partial_flags_have_a_transporter():
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS

    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    ambient = lattice(gram)
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    source = ambient.primitive_isotropic_subobject(*basis[::2])
    target = ambient.primitive_isotropic_subobject(*basis[1::2])
    source_line = basis[0] + basis[2]
    target_line = basis[1] + basis[3]
    source_term = source.inclusion().domain().primitive_isotropic_subobject(source.inclusion().lift(source_line))
    target_term = target.inclusion().domain().primitive_isotropic_subobject(target.inclusion().lift(target_line))
    source_marking = FLAGS((source_term.inclusion(), source.inclusion()))
    target_marking = FLAGS((target_term.inclusion(), target.inclusion()))
    witness = isometry(FLAGS, source_marking, target_marking)
    assert witness is not None
    target_ambient_line = ambient.primitive_isotropic_subobject(target_line)
    assert target_ambient_line.inclusion().is_in_image(witness(source_line))
    for point in basis[::2]:
        assert target.inclusion().is_in_image(witness(point))
    for left in basis:
        for right in basis:
            assert ambient.b(witness(left), witness(right)) == ambient.b(left, right)


def test_three_u_complete_flags_have_a_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory

    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    ambient = lattice(gram)
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))

    def marking(vectors):
        space = ambient.primitive_isotropic_subobject(*vectors)
        domain = space.inclusion().domain()
        plane = domain.primitive_isotropic_subobject(*(space.inclusion().lift(point) for point in vectors[:2]))
        line = plane.inclusion().domain().primitive_isotropic_subobject(plane.inclusion().lift(space.inclusion().lift(vectors[0] + vectors[1])))
        return category((line.inclusion(), plane.inclusion(), space.inclusion()))

    witness = isometry(category, marking(basis[::2]), marking(basis[1::2]))
    assert witness is not None
    for size in (1, 2, 3):
        source_vectors = (basis[0] + basis[2],) if size == 1 else basis[::2][:size]
        target_vectors = (basis[1] + basis[3],) if size == 1 else basis[1::2][:size]
        target_term = ambient.primitive_isotropic_subobject(*target_vectors)
        assert all(target_term.inclusion().is_in_image(witness(point)) for point in source_vectors)
    assert all(ambient.b(witness(left), witness(right)) == ambient.b(left, right) for left in basis for right in basis)


def test_complete_isotropic_flags_in_three_u_split_into_two_special_orthogonal_orbits():
    r"""A maximal isotropic stabilizer has determinant det(A) det(A^-T)=1.

    Thus its O(3U)-orbit has two SO(3U)-orbits. Exchanging a hyperbolic
    pair has determinant -1, so both determinant classes occur in O(3U).
    """
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import _flag_marking_from_basis

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    basis = tuple(ambient.module_generators())[::2]
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    marked = category(_flag_marking_from_basis(ambient, basis, (1, 2, 3)))

    representatives = tuple(ambient.O().split_orbit(marked, ambient.SO()))

    assert len(representatives) == 2
    assert all(representative in category for representative in representatives)
    assert all(
        tuple(int(arrow.domain().rank()) for arrow in representative.arrows()) == (1, 2, 3)
        for representative in representatives
    )


def test_isotropic_rank_above_witt_index_has_no_markings():
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    shape = lattice([[0 for _ in range(4)] for _ in range(4)])
    assert tuple(primitive_marking_classes(ARROWS, shape, ambient)) == ()


def test_three_u_general_induction_classifies_maximal_isotropic_spaces():
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    shape = lattice([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    markings = tuple(primitive_marking_classes(ARROWS, shape, ambient))
    representatives = tuple(
        tuple(marking.arrow()(marking.arrow().domain().module_generator(label)) for label in marking.arrow().domain().module_generating_set())
        for marking in markings
    )
    assert len(representatives) == 1
    for basis in representatives:
        assert len(basis) == 3
        assert all(ambient.b(left, right) == 0 for left in basis for right in basis)
        assert ambient.primitive_isotropic_subobject(*basis).inclusion().domain().rank() == 3


def test_complete_flag_marking_retains_each_primitive_term():
    from sage_indefinite_port.indefinite.marked_lattices import _flag_marking_from_basis

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())[::2]
    arrows = _flag_marking_from_basis(ambient, basis, (1, 2, 3))
    assert tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2, 3)
    assert arrows[-1].codomain() is ambient
    assert all(arrow.is_primitive() for arrow in arrows)
    assert all(left.codomain() is right.domain() for left, right in zip(arrows, arrows[1:]))


def test_partial_flag_above_witt_index_has_no_markings():
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, primitive_marking_classes

    ambient = lattice([[0, 1], [1, 0]])
    plane = lattice([[0, 0], [0, 0]])
    line = plane.primitive_isotropic_subobject(next(iter(plane.module_generators())))
    assert tuple(primitive_marking_classes(FLAGS, line.inclusion(), ambient)) == ()


def test_three_u_general_induction_classifies_complete_and_partial_flags():
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    induction = IsotropicInduction(EichlerRecursion())
    for dimensions in ((1, 3), (1, 2, 3)):
        representatives = induction.isotropic_flag_bases(ambient, dimensions)
        assert len(representatives) == 1
        for basis in representatives:
            assert len(basis) == 3
            assert all(ambient.b(left, right) == 0 for left in basis for right in basis)
            for rank in dimensions:
                assert ambient.primitive_isotropic_subobject(*basis[:rank]).inclusion().domain().rank() == rank


def test_public_partial_two_three_flag_classes_in_three_u():
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, primitive_marking_classes

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    shape = lattice([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    shape_basis = tuple(shape.module_generator(label) for label in shape.module_generating_set())
    first_term = shape.primitive_isotropic_subobject(*shape_basis[:2])
    classes = tuple(primitive_marking_classes(FLAGS, first_term.inclusion(), ambient))
    assert len(classes) == 1
    for flag in classes:
        arrows = tuple(flag.arrows())
        assert tuple(int(arrow.domain().rank()) for arrow in arrows) == (2, 3)
        assert arrows[-1].codomain() is ambient
        assert arrows[0].codomain() is arrows[1].domain()
        assert all(arrow.is_primitive() for arrow in arrows)


def test_public_complete_flag_classes_in_three_u():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    shape = lattice([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    shape_basis = tuple(shape.module_generator(label) for label in shape.module_generating_set())
    plane = shape.primitive_isotropic_subobject(*shape_basis[:2])
    line = plane.inclusion().domain().primitive_isotropic_subobject(plane.inclusion().lift(shape_basis[0] + shape_basis[1]))
    shape_category = Cat().Mor(FiniteOrdinalCategory(3), Lattices(ZZ))
    source = shape_category((line.inclusion(), plane.inclusion()))
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))
    classes = tuple(primitive_marking_classes(category, source, ambient))
    assert len(classes) == 1
    for flag in classes:
        arrows = tuple(flag.arrows())
        assert tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2, 3)
        assert arrows[-1].codomain() is ambient
        assert all(arrow.is_primitive() for arrow in arrows)
        assert all(left.codomain() is right.domain() for left, right in zip(arrows, arrows[1:]))


def test_core_partial_flags_have_a_natural_isomorphism_across_presentations():
    from itertools import pairwise
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS

    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    source = lattice(gram)
    # Replace e_0 by e_0 + f_0 in the target presentation of 3U.
    target_gram = [row[:] for row in gram]
    target_gram[0][0] = 2
    target = lattice(target_gram)
    e = tuple(source.module_generators())
    f = tuple(target.module_generators())

    def marking(ambient, vectors):
        space = ambient.primitive_isotropic_subobject(*vectors)
        line = space.inclusion().domain().primitive_isotropic_subobject(
            space.inclusion().lift(vectors[0])
        )
        return FLAGS((line.inclusion(), space.inclusion()))

    first = marking(source, (e[0] + e[2], e[2], e[4]))
    second = marking(target, (f[0] - f[1] + f[2], f[2], f[4]))
    homset = FLAGS.Core().Mor(first, second)
    assert homset.is_empty() is False
    witness = homset.an_element()
    assert witness.domain().functor() is first.functor()
    assert witness.codomain().functor() is second.functor()
    forward = witness.forward()
    inverse = (~witness).forward()
    shape = FLAGS.domain_category()
    objects = tuple(shape.objects())
    apex = forward.component(objects[-1])
    assert apex.domain() is source and apex.codomain() is target
    assert all(target.b(apex(x), apex(y)) == source.b(x, y) for x in e for y in e)
    for obj in objects:
        assert all(
            inverse.component(obj)(forward.component(obj)(point)) == point
            for point in first.functor()(obj).module_generators()
        )
    for left, right in pairwise(objects):
        arrow = shape.Mor(left, right).unique()
        assert all(
            second.functor()(arrow)(forward.component(left)(point))
            == forward.component(right)(first.functor()(arrow)(point))
            for point in first.functor()(left).module_generators()
        )


def test_zero_isotropic_marking_is_the_unique_empty_embedding():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[0, 1], [1, 0]])
    source = Lattices(ZZ)(0)
    fixed = CosliceCategory(Lattices(ZZ), source)
    for category in (ARROWS, fixed):
        classes = tuple(primitive_marking_classes(category, source, ambient))
        assert len(classes) == 1
        arrow = classes[0].arrow()
        assert int(arrow.domain().rank()) == 0
        assert arrow.codomain() is ambient
        witness = isometry(category, classes[0], classes[0])
        assert witness is not None
        basis = tuple(ambient.module_generators())
        assert all((~witness)(witness(point)) == point for point in basis)
        assert all(ambient.b(witness(left), witness(right)) == ambient.b(left, right) for left in basis for right in basis)
