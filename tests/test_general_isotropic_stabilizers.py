"""Maximal isotropic stabilizers in 3U have their split orthogonal parabolic images.

Over F3 the stabilizer of a maximal isotropic 3-space is
Alt^2(F3^3) semidirect GL3(F3). Its order is 3^3 times
(3^3-1)(3^3-3)(3^3-3^2). Integral GL3 generators map onto GL3(F3)
because the only nonzero determinants in F3 are +1 and -1.
"""
from sage.all import GF, MatrixGroup, matrix

from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators
from tests.acceptance.consumer import lattice


def test_maximal_isotropic_stabilizer_in_three_u_has_full_parabolic_image():
    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    labels = tuple(ambient.module_generating_set())
    basis = tuple(ambient.module_generator(label) for label in labels)
    space = ambient.primitive_isotropic_subobject(*basis[::2])
    generators = tuple(automorphism_group_generators(ARROWS, ARROWS.object(space.inclusion())))
    for generator in generators:
        assert all(space.inclusion().is_in_image(generator(point)) for point in basis[::2])
        assert all(ambient.b(generator(left), generator(right)) == ambient.b(left, right) for left in basis for right in basis)
    images = [matrix(GF(3), [[int(generator(point).to_vector()(label)) for label in labels] for point in basis]) for generator in generators]
    assert MatrixGroup(images).order() == 3**3 * (27 - 1) * (27 - 3) * (27 - 9)


def test_partial_and_complete_flag_stabilizers_have_their_parabolic_images():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import _flag_marking_from_basis

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    labels = tuple(ambient.module_generating_set())
    basis = tuple(ambient.module_generator(label) for label in labels)
    adapted = (basis[0] + basis[2], basis[2], basis[4])
    for dimensions, levi_order in (((1, 3), 2 * 48 * 3**2), ((2, 3), 48 * 2 * 3**2), ((1, 2, 3), 2**3 * 3**3)):
        category = Cat().Mor(FiniteOrdinalCategory(len(dimensions) + 1), Lattices(ZZ))
        marking = category(_flag_marking_from_basis(ambient, adapted, dimensions))
        generators = tuple(automorphism_group_generators(category, marking))
        terms = tuple(ambient.primitive_isotropic_subobject(*adapted[:rank]) for rank in dimensions)
        for generator in generators:
            for rank, term in zip(dimensions, terms, strict=True):
                assert all(term.inclusion().is_in_image(generator(point)) for point in adapted[:rank])
            assert all(ambient.b(generator(left), generator(right)) == ambient.b(left, right) for left in basis for right in basis)
        images = [matrix(GF(3), [[int(generator(point).to_vector()(label)) for label in labels] for point in basis]) for generator in generators]
        assert MatrixGroup(images).order() == 3**3 * levi_order
