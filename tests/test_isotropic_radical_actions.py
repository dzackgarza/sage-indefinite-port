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
