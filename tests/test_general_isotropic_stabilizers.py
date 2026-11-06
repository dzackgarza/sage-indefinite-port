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
