from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.indefinite.edgewalk_rank2 import (
    anisotropic_cycle,
    canonical_companion,
    first_next_vector,
    fixed_norm_vectors_isotropic,
    isotropic_factorization,
    oriented_complement,
    oriented_determinant,
    primitive_isotropic_vectors,
    quadratic_eval,
    scalar_eval,
)


def test_hyperbolic_plane_factorization_and_fixed_norm_enumeration():
    G = matrix(ZZ, [[0, 1], [1, 0]])
    F = isotropic_factorization(G)
    x = vector(ZZ, [3, -2])
    assert F is not None and quadratic_eval(G, x) == -12
    assert scalar_eval(G, x, vector(ZZ, [1, 4])) == 10
    assert (F[0, 0] * x[0] + F[0, 1] * x[1]) * (F[1, 0] * x[0] + F[1, 1] * x[1]) == -12
    assert {tuple(v) for v in primitive_isotropic_vectors(G)} == {(-1, 0), (0, 1)}
    assert {tuple(v) for v in fixed_norm_vectors_isotropic(G, 12)} == {
        (-6, -1),
        (-3, -2),
        (-2, -3),
        (-1, -6),
        (1, 6),
        (2, 3),
        (3, 2),
        (6, 1),
    }


def test_canonical_companion_preserves_orientation_and_bound():
    G = matrix(ZZ, [[2, 0], [0, -3]])
    r = vector(ZZ, [1, 0])
    l = canonical_companion(G, 2, r, oriented_complement(r))
    assert oriented_determinant(r, l) == 1 and quadratic_eval(G, l) <= 2


def test_anisotropic_cycle_closes_by_integral_isometry():
    G = matrix(ZZ, [[2, 0], [0, -3]])
    result = anisotropic_cycle(G, 2)
    assert result is not None
    automorphism, cycle = result
    assert cycle and automorphism.det() in (-1, 1)
    assert automorphism * G * automorphism.transpose() == G
    assert all(quadratic_eval(G, v) > 0 for v in cycle)


def test_first_next_vector_obeys_source_orientation_conditions():
    G = matrix(ZZ, [[0, 1], [1, 0]])
    r = vector(ZZ, [1, 1])
    nxt = first_next_vector(G, r, 12)
    assert nxt is not None and quadratic_eval(G, nxt) == 12
    assert scalar_eval(G, r, nxt) > 0 and oriented_determinant(r, nxt) > 0
