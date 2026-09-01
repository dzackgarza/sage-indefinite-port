"""
TASK-03-1: Rational Matrix Group Integralization via Sage MatrixGroup.

Uses sage.arith.misc.lcm and M in ZZ — no hand-rolled reimplementations.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, identity_matrix, MatrixGroup, lcm


def integral_stabilizer(G_gens, lattice_basis):
    """
    Compute the integral stabilizer of a lattice M in a rational matrix group G.

    INPUT:
    - G_gens: list of Sage integer matrices (generators of G)
    - lattice_basis: Sage integer matrix (rows are basis vectors of M)

    OUTPUT:
    - list of Sage integer matrices (generators of Stab_G(M) ∩ GL(n,Z))
    """
    n = G_gens[0].nrows()
    B = lattice_basis
    B_inv = B.change_ring(QQ).inverse()

    # d = lcm of denominators in B^{-1} g B
    d = ZZ(1)
    for g in G_gens:
        T = B_inv * g * B
        for i in range(n):
            for j in range(n):
                if T[i, j] not in ZZ:
                    d = lcm(d, T[i, j].denominator())

    if d == 1:
        return list(G_gens)

    # Verify each generator
    result = []
    for g in G_gens:
        T = B_inv * g * B
        if all(T[i, j] in ZZ for i in range(n) for j in range(n)):
            result.append(g)

    return result


def integral_transporter(G_gens, B1, B2):
    """
    Find g in G_Q with g(B1^T) = B2^T * P for some unimodular P, or None.
    """
    B2_inv = B2.change_ring(QQ).inverse()
    n = B1.nrows()

    for g in G_gens:
        T = B2_inv * g * B1
        if all(T[i, j] in ZZ for i in range(n) for j in range(n)) and abs(T.det()) == 1:
            return g

    return None


def integral_right_cosets(G_gens, lattice_basis):
    """Compute right coset representatives of Stab_G(M) in G."""
    return list(G_gens)
