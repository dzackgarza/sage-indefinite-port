"""Orbit counts of norm-0 and norm-2 vectors, null planes and flags do not depend on the basis (#18).

A unimodular change of basis P gives the Gram matrix P G P^T of the same lattice, so every
O(L)-orbit count must agree. The expected counts for L = 2U are derived by hand:

L is even unimodular, so D(L) = 0 and L contains 2U. By the Eichler criterion the
O(L)-orbit of a primitive v is fixed by v^2 and v/div(v) in D(L) = 0, hence by v^2 alone:
one orbit of primitive vectors of norm 0 and one of norm 2 (e1 + f1 has norm 2). A
primitive null plane E has a null complement E' with E + E' = L (hyperbolic completion in
a unimodular lattice), so any two primitive null planes are carried to each other by an
isometry; the stabilizer of E induces all of GL(E), which is transitive on primitive
vectors of E, so there is one complete-flag orbit. Hence one null-line, one null-plane and
one complete-flag orbit.
"""

from sage.all import ZZ, matrix

from tests.acceptance.consumer import flag_classes, isotropic_classes, lattice, vector_classes

TWO_U = [
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0],
]
CHANGE = matrix(ZZ, [[1, 2, 0, 1], [0, 1, 1, 0], [0, 0, 1, 3], [0, 0, 0, 1]])


def _changed(gram):
    assert CHANGE.determinant() == 1
    return [[int(entry) for entry in row] for row in (CHANGE * matrix(ZZ, gram) * CHANGE.transpose()).rows()]


def _counts(gram):
    lattice_ = lattice(gram)
    return (
        len(vector_classes(lattice_, 0)),
        len(vector_classes(lattice_, 2)),
        len(isotropic_classes(lattice_, 1)),
        len(isotropic_classes(lattice_, 2)),
        len(flag_classes(lattice_)),
    )


def test_2U_counts_are_basis_independent():
    assert _counts(TWO_U) == _counts(_changed(TWO_U)) == (1, 1, 1, 1, 1)
