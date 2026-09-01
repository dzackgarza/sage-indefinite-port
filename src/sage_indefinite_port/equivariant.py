"""
TASK-06-1: Equivariant Lattices and Involution/Finite-Order Centralizers.

Implement with_isometry(f) and centralizer groups O(L,f) via involution
eigenspaces L_+ + L_- and cyclotomic decompositions.

Port of centralizer computation from polyhedral_common.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import (
    ZZ, QQ, matrix, identity_matrix, vector as sage_vector,
    gcd, factor,
)


class EquivariantLattice:
    """
    Lattice L equipped with an isometry f: L -> L.

    Supports computation of the centralizer O(L, f) = {g in O(L) : g f = f g}.
    """

    def __init__(self, gram, isometry):
        """
        INPUT:
        - gram: nested list (Gram matrix)
        - isometry: nested list (isometry matrix, f^T Q f = Q)
        """
        self.gram = matrix(ZZ, gram)
        self.n = self.gram.nrows()
        self.f = matrix(ZZ, isometry)

        # Verify f is an isometry
        assert self.f.transpose() * self.gram * self.f == self.gram, "f is not an isometry"

    def eigenspace_decomposition(self):
        """
        For involutions (f^2 = I), decompose L into L_+ and L_-.

        OUTPUT:
        - dict with 'L_plus_basis', 'L_minus_basis', 'Q_plus', 'Q_minus'
        """
        f = self.f
        Q = self.gram

        # Eigenvalues of f: for involution, they are +1 and -1
        # L_+ = ker(f - I), L_- = ker(f + I)

        f_minus_I = f - identity_matrix(ZZ, self.n)
        f_plus_I = f + identity_matrix(ZZ, self.n)
        from sage.all import QQ as SageQQ
        L_plus = f_minus_I.change_ring(SageQQ).left_kernel().matrix()
        L_minus = f_plus_I.change_ring(SageQQ).left_kernel().matrix()

        # Gram matrices of the eigenspaces
        Q_plus = L_plus * Q * L_plus.transpose() if L_plus.nrows() > 0 else matrix(ZZ, 0, 0)
        Q_minus = L_minus * Q * L_minus.transpose() if L_minus.nrows() > 0 else matrix(ZZ, 0, 0)

        return {
            "L_plus_basis": L_plus,
            "L_minus_basis": L_minus,
            "Q_plus": Q_plus,
            "Q_minus": Q_minus,
        }

    def centralizer_generators(self):
        """
        Compute generators of O(L, f) = centralizer of f in O(L).

        For an involution f:
        - Generators of O(L_+) x O(L_-) act block-diagonally
        - Additional generators from gluing between eigenspaces

        OUTPUT:
        - list of Sage integer matrices (generators of centralizer)
        """
        from .higher_witt import indefinite_form_automorphism_group

        decomp = self.eigenspace_decomposition()
        Q_plus = decomp["Q_plus"]
        Q_minus = decomp["Q_minus"]
        L_plus = decomp["L_plus_basis"]
        L_minus = decomp["L_minus_basis"]

        generators = []

        # Generators from O(L_+)
        if Q_plus.nrows() > 0:
            Q_plus_list = [[int(Q_plus[i, j]) for j in range(Q_plus.ncols())]
                            for i in range(Q_plus.nrows())]
            O_plus_gens = indefinite_form_automorphism_group(Q_plus_list)

            n_plus = L_plus.nrows()
            for g in O_plus_gens:
                G = matrix(ZZ, g) if not isinstance(g, type(matrix())) else g
                # Embed into full space
                full_gen = identity_matrix(ZZ, self.n)
                for i in range(n_plus):
                    for j in range(n_plus):
                        full_gen[int(L_plus[i].list()[0] if False else i),
                                 int(L_plus[j].list()[0] if False else j)] = G[i, j]
                generators.append(full_gen)

        # Generators from O(L_-)
        if Q_minus.nrows() > 0:
            Q_minus_list = [[int(Q_minus[i, j]) for j in range(Q_minus.ncols())]
                             for i in range(Q_minus.nrows())]
            O_minus_gens = indefinite_form_automorphism_group(Q_minus_list)

            n_minus = L_minus.nrows()
            n_plus = L_plus.nrows()
            for g in O_minus_gens:
                G = matrix(ZZ, g) if not isinstance(g, type(matrix())) else g
                full_gen = identity_matrix(ZZ, self.n)
                for i in range(n_minus):
                    for j in range(n_minus):
                        full_gen[n_plus + i, n_plus + j] = G[i, j]
                generators.append(full_gen)

        # Gluing generators: maps that swap between eigenspaces
        # These exist when Q_+ and Q_- share common subquotients
        # Placeholder for full implementation

        return generators


def with_isometry(gram, isometry):
    """
    Create an EquivariantLattice from a lattice and an isometry.

    INPUT:
    - gram: nested list (Gram matrix)
    - isometry: nested list (isometry matrix)

    OUTPUT:
    - EquivariantLattice object
    """
    return EquivariantLattice(gram, isometry)


def centralizer_group(gram, isometry):
    """
    Compute generators of the centralizer O(L, f).

    INPUT:
    - gram: nested list (Gram matrix)
    - isometry: nested list (isometry matrix)

    OUTPUT:
    - list of nested lists (generator matrices)
    """
    eq = EquivariantLattice(gram, isometry)
    gens = eq.centralizer_generators()
    n = len(gram)
    return [[[int(g[i, j]) for j in range(n)] for i in range(n)] for g in gens]


def cyclotomic_decomposition(gram, isometry):
    """
    For finite-order isometry f of order m, decompose L via cyclotomic polynomial.

    L tensor Q = direct sum of eigenspaces for roots of unity.

    INPUT:
    - gram: nested list
    - isometry: nested list

    OUTPUT:
    - dict with decomposition data
    """
    f = matrix(ZZ, isometry)
    n = f.nrows()

    # Compute order of f
    order = 1
    power = f
    I = identity_matrix(ZZ, n)
    for k in range(1, 100):
        if power == I:
            order = k
            break
        power = power * f

    # Cyclotomic factorization of x^order - 1
    from sage.all import cyclotomic_polynomial
    Qx = ZZ['x']
    poly = Qx(x) ** order - 1

    factors = list(poly.factor())

    return {
        "order": order,
        "factors": [(str(p), e) for p, e in factors],
    }
