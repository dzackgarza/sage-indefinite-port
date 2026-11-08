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


def test_special_orthogonal_group_has_two_moving_maximal_isotropic_orbits_in_four_u():
    r"""The maximal parabolic in O(4U) has determinant one.

    Unimodularity supplies integral duals to a primitive null basis;
    evenness permits null corrections, giving a hyperbolic completion.
    Its diagonal blocks are A and A^(-T), and its unipotent radical has
    determinant one. Swapping one hyperbolic pair has determinant minus
    one, so the single O(4U) orbit splits into two SO(4U) orbits.
    """
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS

    hyperbolic = Lattices(ZZ)("U")
    ambient = hyperbolic + hyperbolic + hyperbolic + hyperbolic
    basis = tuple(ambient.module_generators())
    null = ambient.primitive_isotropic_subobject(*basis[::2])
    marking = ARROWS(null.inclusion())
    group, special = ambient.O(), ambient.SO()
    representatives = tuple(group.split_orbit(marking, special))
    assert len(representatives) == 2
    for representative in representatives:
        arrow = representative.arrow()
        assert arrow.codomain() is ambient
        assert arrow.is_primitive()
        assert arrow.domain().rank() == 4
        null_basis = tuple(arrow(point) for point in arrow.domain().module_generators())
        assert all(ambient.b(x, y) == ZZ.zero() for x in null_basis for y in null_basis)
        self_witness = ARROWS.Core().Mor(representative, representative).an_element_in(special)
        assert self_witness is not None
        apex = self_witness.forward().right()
        assert apex in special
        assert all(arrow.is_in_image(apex(point)) for point in null_basis)
    homset = ARROWS.Core().Mor(*representatives)
    assert homset.an_element_in(special) is None
    witness = homset.an_element().forward().right()
    assert witness.determinant() == -ZZ.one()
    first, second = (representative.arrow() for representative in representatives)
    assert all(second.is_in_image(witness(first(point))) for point in first.domain().module_generators())
    assert all(first.is_in_image((~witness)(second(point))) for point in second.domain().module_generators())
    assert all((~witness)(witness(point)) == point for point in basis)
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in basis for y in basis)
