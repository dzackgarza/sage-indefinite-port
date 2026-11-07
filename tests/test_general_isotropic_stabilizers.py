"""Isotropic flag stabilizers have their split orthogonal parabolic images.

Over F3 the stabilizer of a maximal isotropic 3-space is
Alt^2(F3^3) semidirect GL3(F3). Its order is 3^3 times
(3^3-1)(3^3-3)(3^3-3^2). Integral GL3 generators map onto GL3(F3)
because the only nonzero determinants in F3 are +1 and -1.
In 4U the flag (1,3,4) has GL4 blocks (1,2,1), with five
unipotent block entries and a six-dimensional orthogonal radical.
"""
import pytest
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


@pytest.mark.parametrize(
    "hyperbolic_rank, flags",
    (
        (3, (((1, 3), 2 * 48 * 3**2), ((2, 3), 48 * 2 * 3**2), ((1, 2, 3), 2**3 * 3**3))),
        (4, (((1, 3, 4), 2 * 48 * 2 * 3**5),)),
    ),
    ids=("three_u", "four_u"),
)
def test_partial_and_complete_flag_stabilizers_have_their_parabolic_images(hyperbolic_rank, flags):
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
    from sage_indefinite_port.indefinite.marked_lattices import _flag_marking_from_basis

    size = 2 * hyperbolic_rank
    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(size)] for i in range(size)])
    labels = tuple(ambient.module_generating_set())
    basis = tuple(ambient.module_generator(label) for label in labels)
    adapted = (basis[0] + basis[2], *basis[2::2])
    for dimensions, levi_order in flags:
        category = Cat().Mor(FiniteOrdinalCategory(len(dimensions) + 1), Lattices(ZZ))
        marking = category(_flag_marking_from_basis(ambient, adapted, dimensions))
        generators = tuple(automorphism_group_generators(category, marking))
        terms = tuple(ambient.primitive_isotropic_subobject(*adapted[:rank]) for rank in dimensions)
        for generator in generators:
            for rank, term in zip(dimensions, terms, strict=True):
                assert all(term.inclusion().is_in_image(generator(point)) for point in adapted[:rank])
            assert all(ambient.b(generator(left), generator(right)) == ambient.b(left, right) for left in basis for right in basis)
        images = [matrix(GF(3), [[int(generator(point).to_vector()(label)) for label in labels] for point in basis]) for generator in generators]
        # The maximal orthogonal radical is Alt^2(F3^r). The Levi factor
        # is the indicated flag parabolic in GL_r(F3); its block-upper
        # entries contribute one F3 factor for each pair of distinct blocks.
        assert MatrixGroup(images).order() == 3**(hyperbolic_rank * (hyperbolic_rank - 1) // 2) * levi_order


def test_fixed_rank_two_embedding_has_pointwise_parabolic_image():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(4)] for i in range(4)])
    labels = tuple(ambient.module_generating_set())
    basis = tuple(ambient.module_generator(label) for label in labels)
    source = lattice([[0, 0], [0, 0]])
    category = CosliceCategory(Lattices(ZZ), source)
    marked = category.object(source.hom([basis[0], basis[2]], codomain=ambient))
    generators = tuple(automorphism_group_generators(category, marked))
    assert all(generator(point) == point for generator in generators for point in basis[::2])
    images = [matrix(GF(3), [[int(generator(point).to_vector()(label)) for label in labels] for point in basis]) for generator in generators]
    # The pointwise maximal parabolic in O(2U) is Alt^2(Z^2).
    assert MatrixGroup(images).order() == 3


def test_fixed_embedding_transporter_respects_every_source_generator():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(4)] for i in range(4)])
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    source = lattice([[0, 0], [0, 0]])
    category = CosliceCategory(Lattices(ZZ), source)
    first = category.object(source.hom([basis[0], basis[2]], codomain=ambient))
    second = category.object(source.hom([basis[0], basis[2] + basis[0]], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(basis[0]) == basis[0]
    assert witness(basis[2]) == basis[2] + basis[0]
    assert (~witness)(basis[2] + basis[0]) == basis[2]
    assert all(ambient.b(witness(left), witness(right)) == ambient.b(left, right) for left in basis for right in basis)
