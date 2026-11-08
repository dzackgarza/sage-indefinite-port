"""The pointwise maximal isotropic stabilizer in 4U is the full skew group.

Write the hyperbolic basis as e_i,f_i with b(e_i,f_j)=delta_ij.
An isometry fixing every e_i must send f_i to f_i+sum_j B_ij e_j.
The pairings of its f-images give B_ij+B_ji=0. Conversely every integral
skew B defines such an isometry, and composition adds B. Thus the returned
parameter vectors must generate the entire integral module of skew matrices.
"""

from itertools import combinations

from sage_extensions.all import Lattices, ZZ, cardinal


def test_pointwise_maximal_stabilizer_in_four_u_has_the_full_skew_parameter_module():
    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic + hyperbolic
    basis = tuple(ambient.module_generators())
    null, dual = basis[::2], basis[1::2]
    subspace = ambient.primitive_isotropic_subobject(*null)
    family = tuple(ambient.O().pointwise_stabilizer(subspace.inclusion()).group_generators())
    labels = tuple(ambient.module_generating_set())
    pairs = tuple(combinations(range(len(null)), 2))
    parameters = ZZ.free_module(len(pairs))
    parameter_labels = tuple(parameters.module_generating_set())
    images = []
    for generator in family:
        assert all(generator(point) == point for point in null)
        coefficients = tuple(
            tuple(generator(point).to_vector()(labels[2 * j]) for j in range(len(null)))
            for point in dual
        )
        for i, point in enumerate(dual):
            assert coefficients[i][i] == ZZ.zero()
            assert generator(point) == point + sum(
                (coefficients[i][j] * null[j] for j in range(len(null))), ambient.zero()
            )
            assert (~generator)(point) == point - sum(
                (coefficients[i][j] * null[j] for j in range(len(null))), ambient.zero()
            )
            for j in range(len(null)):
                assert coefficients[i][j] == -coefficients[j][i]
        images.append(parameters.linear_combination({
            label: coefficients[i][j]
            for label, (i, j) in zip(parameter_labels, pairs, strict=True)
        }))
    assert parameters.subobject_on(images).inclusion().cokernel().cardinality() == cardinal(1)
