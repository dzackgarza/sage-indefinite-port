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
