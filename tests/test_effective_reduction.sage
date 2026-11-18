r"""Proved transporter decisions.

U + U is unimodular and swapping the two hyperbolic planes carries e1 to e3, an isometry with entries
at most one. In U + U(2) the null vector e1 has divisor 1 and e3 has divisor 2, so no isometry joins them.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage_indefinite_port.indefinite.effective_reduction import (
    centralizer_transporter_in_box, orthogonal_transporter, separated_by_invariant,
)

U_U = [[int(0), int(1), int(0), int(0)], [int(1), int(0), int(0), int(0)],
       [int(0), int(0), int(0), int(1)], [int(0), int(0), int(1), int(0)]]
U_U2 = [[int(0), int(1), int(0), int(0)], [int(1), int(0), int(0), int(0)],
        [int(0), int(0), int(0), int(2)], [int(0), int(0), int(2), int(0)]]


def test_null_vectors_of_u_plus_u_are_joined():
    x = sage_matrix(SageZZ, 4, 1, [int(1), int(0), int(0), int(0)])
    y = sage_matrix(SageZZ, 4, 1, [int(0), int(0), int(1), int(0)])
    assert not separated_by_invariant(U_U, x, y, [])
    g = centralizer_transporter_in_box(U_U, [], x, y, int(1))
    gram = sage_matrix(SageZZ, U_U)
    assert g is not None and g * x == y and g.transpose() * gram * g == gram


def test_different_divisors_have_no_transporter():
    x = sage_matrix(SageZZ, 4, 1, [int(1), int(0), int(0), int(0)])
    y = sage_matrix(SageZZ, 4, 1, [int(0), int(0), int(1), int(0)])
    assert separated_by_invariant(U_U2, x, y, [])
    assert orthogonal_transporter(U_U2, x, y) is None
