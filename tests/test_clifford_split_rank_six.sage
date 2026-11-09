r"""The split rank-six model. C_0 has dimension 32 = dim(M_4 x M_4).

By hand the transvection e2 -> e1 + e2 sends e23 to e13 + e23 and e24 to e14 + e24, and
fixes e12, e13, e14, e34. So in the basis (e12,e34,e13,-e24,e14,e23) it fixes
all vectors except f2 = -e24, sent to -e14 + f2, and f3 = e23, sent to e13 + f3.
-1 acts trivially on Lambda^2, and the transvection does not.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.clifford_split_rank_six import (
    even_part_dimension, fock_generators, isometry_of_special_linear, second_compound, three_hyperbolic_planes,
)


def _integers(rows):
    return sage_matrix(SageQQ, [[SageQQ(int(x)) for x in row] for row in rows])


def _transvection():
    return _integers([[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def test_the_fock_generators_satisfy_the_clifford_relations_of_3u():
    gammas = fock_generators()
    pairing = [[1 if {a, b} in ({0, 1}, {2, 3}, {4, 5}) else 0 for b in range(6)] for a in range(6)]
    identity = gammas[0].parent().identity_matrix()
    for a in range(int(6)):
        for b in range(int(6)):
            assert gammas[a] * gammas[b] + gammas[b] * gammas[a] == int(2) * int(pairing[a][b]) * identity


def test_the_even_part_is_all_of_two_matrix_blocks():
    assert even_part_dimension() == int(32)


def test_the_transvection_has_the_hand_action():
    action = second_compound(_transvection())
    expected = _integers([
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, -1, 1, 0],
        [0, 0, 0, 0, 0, 1],
    ])
    assert action == expected


def test_the_scalar_kernel_contains_minus_one_and_not_the_transvection():
    minus = _integers([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
    identity = minus.parent().identity_matrix()[:int(4), :int(4)]
    assert second_compound(minus) == second_compound(identity)
    assert second_compound(_transvection()) != second_compound(identity)


def test_the_compound_is_multiplicative_and_preserves_3u_integrally():
    lattice = three_hyperbolic_planes()
    g = _transvection()
    h = _integers([[1, 0, 0, 0], [2, 1, 0, 0], [0, 0, 1, 3], [0, 0, 0, 1]])
    assert second_compound(g * h) == second_compound(g) * second_compound(h)
    forward = isometry_of_special_linear(g * h, lattice)
    points = tuple(lattice.module_generators())
    for x in points:
        for y in points:
            assert lattice.b(forward(x), forward(y)) == lattice.b(x, y)
