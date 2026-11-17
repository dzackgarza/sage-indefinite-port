"""The radical action on an unflagged rank-three isotropic space is GL(3,Z).

Its reduction modulo two is GL(3,F_2), of order 168.
"""
from sage.all import GF, MatrixGroup, matrix

from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
from tests.acceptance.consumer import lattice


def test_rank_three_radical_action_has_full_mod_two_image():
    radical = lattice([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    subobject = radical.primitive_isotropic_subobject(*tuple(radical.module_generators()))
    images = IsotropicInduction.sublattice_action_images(subobject)
    labels = tuple(radical.module_generating_set())
    generators = [matrix(GF(2), [[int(point.to_vector()(label)) for label in labels] for point in family]) for family in images]
    assert MatrixGroup(generators).order() == 168


def test_rank_three_partial_flag_action_is_a_proper_subgroup():
    radical = lattice([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    subobject = radical.primitive_isotropic_subobject(*tuple(radical.module_generators()))
    inclusion = subobject.inclusion()
    domain = inclusion.domain()
    basis = tuple(domain.module_generator(label) for label in domain.module_generating_set())
    line = domain.primitive_isotropic_subobject(basis[0] + basis[1])
    adapted, images = IsotropicInduction.flag_action_data(subobject, (line,))
    assert line.inclusion().is_in_image(inclusion.lift(adapted[0]))
    labels = tuple(radical.module_generating_set())
    basis_matrix = matrix(GF(2), [[int(point.to_vector()(label)) for label in labels] for point in adapted])
    generators = [~basis_matrix * matrix(GF(2), [[int(point.to_vector()(label)) for label in labels] for point in family]) for family in images]
    assert MatrixGroup(generators).order() == 24


def test_rank_two_induction_extension_stabilizer_has_expected_mod_two_image():
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    vectors = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    space = ambient.primitive_isotropic_subobject(vectors[0], vectors[2])
    induction = IsotropicInduction(EichlerRecursion())
    reduction = space.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    isotropic = next(iter(quotient.module_generators()))
    basis, images = induction.extension_stabilizer_images(reduction, induction.plane_vectors(space), induction.sublattice_action_images(space), isotropic)
    labels = tuple(ambient.module_generating_set())
    basis_matrix = matrix(GF(2), [[int(point.to_vector()(label)) for label in labels] for point in basis])
    generators = []
    for family in images:
        rows = matrix(GF(2), [[int(point.to_vector()(label)) for label in labels] for point in family])
        generators.append(basis_matrix.transpose().solve_right(rows.transpose()).transpose())
    assert MatrixGroup(generators).order() == 96


def test_extension_stabilizer_does_not_negate_a_pointwise_marked_plane():
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    vectors = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(vectors[0], vectors[2])
    induction = IsotropicInduction(EichlerRecursion())
    reduction = plane.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    radical = induction.plane_vectors(plane)
    basis, families = induction.extension_stabilizer_images(
        reduction, radical, (), next(iter(quotient.module_generators()))
    )
    assert families
    assert all(images[:len(radical)] == radical for images in families)
    assert any(images != basis for images in families)
    extension = ambient.primitive_isotropic_subobject(*radical, basis[len(radical)])
    for images in families:
        assert all(ambient.b(left, right) == ambient.b(image_left, image_right)
                   for left, image_left in zip(basis, images, strict=True)
                   for right, image_right in zip(basis, images, strict=True))
        # The chosen quotient null line is the first quotient-frame vector.
        # The entire primitive extension must be preserved setwise.
        assert all(extension.inclusion().is_in_image(point)
                   for point in images[:len(radical) + 1])


def test_extension_stabilizer_preserves_a_proper_flag_on_the_radical():
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(8)] for i in range(8)])
    vectors = tuple(ambient.module_generators())
    space = ambient.primitive_isotropic_subobject(vectors[0], vectors[2], vectors[4])
    inclusion = space.inclusion()
    domain = inclusion.domain()
    first = domain.primitive_isotropic_subobject(next(iter(domain.module_generators())))
    induction = IsotropicInduction(EichlerRecursion())
    radical, action = induction.flag_action_data(space, (first,))
    reduction = space.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    basis, families = induction.extension_stabilizer_images(
        reduction, radical, action, next(iter(quotient.module_generators()))
    )
    line = ambient.primitive_isotropic_subobject(radical[0])
    extension = ambient.primitive_isotropic_subobject(*radical, basis[len(radical)])
    assert families
    assert all(line.inclusion().is_in_image(images[0]) for images in families)
    assert all(extension.inclusion().is_in_image(point)
               for images in families for point in images[:len(radical) + 1])


def test_pointwise_extension_stabilizer_has_quotient_sign_without_radical_sign():
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    e, _, f, _, _, _ = tuple(ambient.module_generators())
    plane = ambient.primitive_isotropic_subobject(e, f)
    induction = IsotropicInduction(EichlerRecursion())
    reduction = plane.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    radical = induction.plane_vectors(plane)
    basis, families = induction.extension_stabilizer_images(
        reduction, radical, (), next(iter(quotient.module_generators()))
    )
    sign = families[0]
    assert sign[:len(radical)] == radical
    assert sign[len(radical):] == tuple(-point for point in basis[len(radical):])
    assert all(images[:len(radical)] == radical for images in families)
