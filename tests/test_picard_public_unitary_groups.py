"""Full integral unitary groups as commutants of CM-order actions on trace lattices."""

import pytest

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.rings.integer_ring import ZZ as SageZZ

from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite import euclidean_picard, gaussian_picard


@pytest.mark.parametrize("d", (1, 2, 7, 11))
def test_picard_full_integral_unitary_generators_over_imaginary_quadratic_orders(d):
    r"""Full U(2,1;O_{-d}) through a skew-adjoint order generator, not a finite scalar group.

    In the three source-supported non-Gaussian Euclidean imaginary quadratic
    orders the only norm-one scalar units are ±1. Hence centralizing a scalar
    isometry cannot recover O-linear unitary maps: -1 is central in all
    O(Tr h). A primitive field element instead acts as a skew-adjoint
    *endomorphism*, whose rational commutant is the full K-linear algebra.
    The method works for the Gaussian Picard source as well.
    """
    hermitian = gaussian_picard.lattice() if d == 1 else euclidean_picard.lattice(d)
    trace = hermitian.trace_lattice()
    field = hermitian.field()
    action = hermitian.order_action(field.gen())
    frame = tuple(trace.module_generators())
    labels = tuple(trace.module_generating_set())
    gram = matrix(SageZZ, [[int(trace.b(x, y)) for y in frame] for x in frame])
    columns = matrix(SageZZ, [
        [int(action(v).to_vector()(label)) for v in frame] for label in labels
    ])
    assert columns.transpose() * gram == -gram * columns
    assert IsometryCommutant(gram, (columns,)).is_compact() is False

    generators = tuple(hermitian.unitary_group_generators())
    assert generators
    assert all(g.domain() is trace and g.codomain() is trace for g in generators)
    assert all(hermitian.is_order_linear(g) for g in generators)
    matrices = [matrix(SageZZ, [
        [int(g(v).to_vector()(label)) for v in frame] for label in labels
    ]) for g in generators]
    one = identity_matrix(SageZZ, int(trace.rank()))
    # A nontrivial unipotent is an explicit infinite-order integral unitary.
    assert any(M != one and (M - one)**3 == 0 for M in matrices)
