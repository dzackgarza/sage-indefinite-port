"""Spacelike vector classes at positive index 1 (#18, marked_lattices.py spacelike branch).

L = U + <-2> has signature (1, 2), so its vectors of square -2 are spacelike. Identify L
with the integral symmetric matrices M = [[a, c], [c, b]] under q(M) = 2 det M, with
e, f, w the matrices of (a, b, c) = (1, 0, 0), (0, 1, 0), (0, 0, 1): then e^2 = f^2 = 0,
e.f = 1, w^2 = -2. GL_2(Z) acts by M -> g M g^T and preserves q and L. A vector of square
-2 is a symmetric M with det M = -1, a unimodular indefinite binary form. Such forms are
odd (equivalent to x^2 - y^2) or even (equivalent to 2xy) [Serre, A Course in
Arithmetic, Ch. V, section 2.2], so GL_2(Z) has two orbits on them. These are not merged by O(L):
their representatives e - f and w have divisors 1 and 2, an O(L)-invariant. Every vector of
square -2 is primitive. Hence exactly two classes, with divisors {1, 2}.
"""

from tests.acceptance.consumer import marked_vector, lattice, vector_classes

U_PLUS_MINUS_TWO = [[0, 1, 0], [1, 0, 0], [0, 0, -2]]


def test_spacelike_square_minus_two_has_two_classes_with_divisors_one_and_two():
    lattice_ = lattice(U_PLUS_MINUS_TWO)
    assert tuple(int(index) for index in lattice_.signature_pair()) == (1, 2)

    classes = vector_classes(lattice_, -2)

    vectors = [marked_vector(marked) for marked in classes]
    assert all(int(vector.q()) == -2 for vector in vectors)
    assert sorted(int(vector.div()) for vector in vectors) == [1, 2]
