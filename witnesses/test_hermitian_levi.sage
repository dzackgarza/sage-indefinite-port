r"""Rank-two Levi images whose negative lift decisions need long computations.

Split: Z[i], pi = 1 + i, L = Pi O^5 with Pi = diag(1, pi, 1, 1, 1/conj(pi)). The pair
(A, U) lifts exactly when A preserves the line of e_1 in F_2^2 = O^2 / pi O^2, so the
index is |P^1(F_2)| = 3.

Glued: Z[omega], theta = omega - conj(omega), h(w, w) = 9, and the glue vector
(e_1 + w) / theta. A must preserve the line of e_1 in F_3^2, of index |P^1(F_3)| = 4,
and act on it by the unit u modulo theta, of index 2. The index is 8.
"""

from sage.all import *
from sage_extensions.categories.hermitian_lattices import HermitianLattice, HermitianWittFrame


def test_split_gaussian_rank_two_levi_index():
    E = CyclotomicField(4)
    p = 1 + E.gen()
    J = matrix(E, 5, 5, [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0])
    Pi = diagonal_matrix(E, [1, p, 1, 1, 1 / p.conjugate()])
    G = Pi.conjugate_transpose() * J * Pi
    frame = HermitianWittFrame(HermitianLattice((E.one(), E.gen()), G.transpose()), ~Pi, matrix(E, 1, 1, [1]), 2)
    index, generators = frame.levi_image_schreier()
    assert index == 3
    assert all(frame.levi_lift(a, u) is not None for a, u in generators)


def test_glued_eisenstein_rank_two_levi_index():
    E = CyclotomicField(3)
    w = E.gen()
    th = w - w.conjugate()
    H = E(9)
    J = matrix(E, 5, 5, [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, H, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0])
    Pi = matrix(E, 5, 5, [1, 0, 1 / th, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 / th, 0, 0, 0, 0, 0, th, 0, 0, 0, 0, 0, 1])
    G = Pi.conjugate_transpose() * J * Pi
    frame = HermitianWittFrame(HermitianLattice((E.one(), w), G.transpose()), ~Pi, matrix(E, 1, 1, [H]), 2)
    index, generators = frame.levi_image_schreier()
    assert index == 8
    assert all(frame.levi_lift(a, u) is not None for a, u in generators)
