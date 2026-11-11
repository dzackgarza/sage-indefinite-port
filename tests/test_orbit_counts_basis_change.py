"""Orbit counts of norm-0 and norm-2 vectors, planes and flags do not depend on the basis (#18).

A unimodular change of basis P gives the Gram matrix P G P^T of the same lattice, so every
O(L)-orbit count must agree. The expected counts are derived by hand:

2U + A2 (A2 negative definite): D(L) = Z/3 with q = -2/3 on both nonzero classes, and
L contains 2U, so by the Eichler criterion the O(L)-orbit of a primitive v is fixed by
v^2 and v/div(v) in D(L). For v^2 = 0 only div 1 occurs (D has no nonzero isotropic class):
one orbit. For v^2 = 2, div 3 would need q(v/3) = 2/9, which is not a value of q, and
div 2 is impossible since D has no 2-torsion: one orbit. Dawes arXiv:2205.10601 (fixture
dawes_2U_A2_stable) gives one plane and one flag orbit under the stable O^+(L) <= O(L),
so O(L) has one of each.

2U + <-2> + <-6>: lines Ze1 and Zl, l = 2e2 + 2f2 + v + w, have divisors 1 and 2; planes
E1 = <e1, e2> and E2 = <e1, l> have E^perp/E of determinant 12 and 3; flags e1 < E1,
e1 < E2 and l < E2 differ in these invariants. These lower bounds meet the O^+(L)
counts (2, 2, 3) of Dawes (fixture dawes_2U_minus2_minus6), which bound the O(L) counts
from above. Primitive isotropic vectors: v and -v lie in one O(L)-orbit, so 2.
"""

from sage.all import ZZ, matrix

from tests.acceptance.consumer import flag_classes, isotropic_classes, lattice, vector_classes

TWO_U_A2 = [
    [0, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, -2, 1],
    [0, 0, 0, 0, 1, -2],
]
TWO_U_MINUS2_MINUS6 = [
    [0, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, -2, 0],
    [0, 0, 0, 0, 0, -6],
]


def _changed(gram):
    """P G P^T for P the upper unitriangular all-ones matrix, which has determinant 1."""
    rank = len(gram)
    change = matrix(ZZ, [[1 if j >= i else 0 for j in range(rank)] for i in range(rank)])
    assert change.determinant() == 1
    return [[int(entry) for entry in row] for row in (change * matrix(ZZ, gram) * change.transpose()).rows()]


def _isotropic_counts(gram):
    lattice_ = lattice(gram)
    return (
        len(vector_classes(lattice_, 0)),
        len(isotropic_classes(lattice_, 1)),
        len(isotropic_classes(lattice_, 2)),
        len(flag_classes(lattice_)),
    )


def test_2U_A2_counts_are_basis_independent():
    for gram in (TWO_U_A2, _changed(TWO_U_A2)):
        assert _isotropic_counts(gram) == (1, 1, 1, 1)
        assert len(vector_classes(lattice(gram), 2)) == 1


def test_2U_minus2_minus6_counts_are_basis_independent():
    for gram in (TWO_U_MINUS2_MINUS6, _changed(TWO_U_MINUS2_MINUS6)):
        assert _isotropic_counts(gram) == (2, 2, 2, 3)
