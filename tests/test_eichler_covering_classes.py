"""Eichler covering classes under the full signed coordinate action."""

from itertools import combinations

from sage_extensions.categories.eichler_criterion import TwoUEichlerModel

from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from tests.acceptance.consumer import lattice


def test_covering_classes_of_two_u_plus_four_negative_roots():
    """Weight and norm determine the discriminant classes modulo coordinate permutations."""
    complement = lattice([[-2 * int(i == j) for j in range(4)] for i in range(4)])
    model = TwoUEichlerModel(complement)
    ambient = model.lattice()
    ring = ambient.base_ring()
    e, f, _u, _v = model.hyperbolic_basis()
    inclusion = model.complement_inclusion()
    roots = tuple(complement.module_generators())
    for square, weight in ((0, 4), (2, 3)):
        representatives = tuple(OrbitCoverModel.split_covering_vectors(model, ring(square)))
        assert len(representatives) == 2
        assert all(int(point.q()) == square and point.is_primitive() for point in representatives)
        assert sorted(int(point.div()) for point in representatives) == [1, 2]
        expected = {
            ambient.divided_discriminant_class(
                2 * e + 2 * f + inclusion(sum(selected, complement.zero()))
            )
            for selected in combinations(roots, weight)
        }
        actual = {
            ambient.divided_discriminant_class(point)
            for point in representatives if int(point.div()) == 2
        }
        assert len(actual) == 1 and actual <= expected
        assert all(int(point.additive_order()) == 2 for point in expected)
    assert not tuple(OrbitCoverModel.split_covering_vectors(model, ring(1)))
