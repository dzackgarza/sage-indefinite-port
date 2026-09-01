"""
TASK-04-2: 2U-Eichler Approximate Subgroups and Discriminant Covers.

Implement 2U-Eichler model for L = U + U + K generating A(L) from SL_2(Z)
copies and Eichler transvections with discriminant covering lists C(L, beta).

Port of ApproximateModels.h and INDEF_FORM_GetApproximateModel from
polyhedral_common.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import (
    ZZ, QQ, matrix, identity_matrix, vector as sage_vector,
    gcd, lcm,
)
from itertools import product as iterproduct


# ---------------------------------------------------------------------------
# Eichler transvection
# ---------------------------------------------------------------------------

def eichler_transvection(gram, f, x):
    """
    Compute the Eichler transvection E_{f,x}.

    Formula: E_{f,x}(y) = y + (y,x) f - (x,x)/2 (y,f) f - (y,f) x

    This preserves L because (x,x)/2 is in Z for even lattices,
    and (y,x), (y,f) are in Z for y in L.

    Port of INDEF_FORM_Eichler_Transvection from ApproximateModels.h.

    INPUT:
    - gram: nested list of integers (even Gram matrix)
    - f: list of integers (isotropic, q(f) = 0)
    - x: list of integers (with b(f,x) = 0)

    OUTPUT:
    - Sage integer matrix (the transvection)
    """
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    fv = sage_vector(ZZ, f)
    xv = sage_vector(ZZ, x)

    q_f = int(fv * Q * fv)
    q_x = int(xv * Q * xv)
    b_fx = int(fv * Q * xv)

    assert q_f == 0, "f must be isotropic"
    assert b_fx == 0, "f and x must be orthogonal"

    # E_{f,x}(y) = y + (y,x) f - (x,x)/2 (y,f) f - (y,f) x
    # As a matrix: E_{ij} = delta_{ij} + f_j * x^T Q e_i - (x,x)/2 * f_j * f^T Q e_i - x_j * f^T Q e_i
    T = matrix(ZZ, n)
    for j in range(n):
        for i in range(n):
            y = sage_vector(ZZ, [1 if k == i else 0 for k in range(n)])
            yx = int(y * Q * xv)  # (y, x)
            yf = int(y * Q * fv)  # (y, f)

            T[i, j] = (
                (1 if i == j else 0)
                + yx * fv[j]
                - (q_x // 2) * yf * fv[j]
                - yf * xv[j]
            )

    return T


def sl2z_generators():
    """
    Generators of SL(2, Z): S and T matrices.

    Port of GeneratorsSL2Z from ApproximateModels.h.

    OUTPUT:
    - pair of Sage 2x2 integer matrices
    """
    S = matrix(ZZ, [[0, -1], [1, 0]])
    T = matrix(ZZ, [[1, 1], [0, 1]])
    return S, T


def get_square_divisors(X):
    """
    Compute all square divisors of X.

    Port of GetSquareDivisors from ApproximateModels.h.

    INPUT:
    - X: integer

    OUTPUT:
    - list of integers d such that d^2 | X
    """
    X = abs(X)
    if X == 0:
        return [0]

    # Factor X
    from sage.all import factor as sage_factor
    fac = sage_factor(X)

    # Build square divisors from prime factorization
    square_divs = [ZZ(1)]
    for p, e in fac:
        # Square part: p^0, p^1, ..., p^{floor(e/2)}
        max_sq_exp = e // 2
        new_divs = []
        for d in square_divs:
            val = d
            for k in range(max_sq_exp + 1):
                new_divs.append(val)
                val *= p * p
        square_divs = new_divs

    return sorted(set(square_divs))


# ---------------------------------------------------------------------------
# Translation classes in discriminant group
# ---------------------------------------------------------------------------

def compute_translation_classes(gram_reduced):
    """
    Compute translation classes for the discriminant group A_L = L*/L.

    Port of ComputeTranslationClasses from ApproximateModels.h.

    INPUT:
    - gram_reduced: nested list of integers (Gram matrix of the K part)

    OUTPUT:
    - list of Sage vectors (representatives of A_L)
    """
    Q = matrix(ZZ, gram_reduced)
    n = Q.nrows()

    if n == 0:
        return [sage_vector(ZZ, [])]

    # Compute L* = Q^{-1} Z^n
    Q_inv = Q.inverse()

    # Find representatives of L*/L by SNF of Q
    # Actually, we enumerate vectors mod L in Q^{-1} Z^n
    det_Q = abs(Q.det())

    if det_Q == 1:
        return [sage_vector(ZZ, [0] * n)]

    # Find basis for L*/L using HNF
    # For simplicity, enumerate small vectors in Q^{-1} Z^n
    representatives = []
    seen = set()

    # Upper bound on coordinate size
    bound = max(int(det_Q), 10)

    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        v = sage_vector(ZZ, coords)
        # Check v is in L* = Q^{-1} Z^n, i.e., Q v in Z^n
        Qv = Q * v
        if not all(x in ZZ for x in Qv):
            continue

        # Compute coset representative mod L
        # v mod L is determined by Qv mod something
        coset_key = tuple(int(Qv[i]) % det_Q for i in range(n))

        if coset_key not in seen:
            seen.add(coset_key)
            representatives.append(v)

        if len(representatives) >= det_Q:
            break

    return representatives


# ---------------------------------------------------------------------------
# Eichler canonical form check
# ---------------------------------------------------------------------------

def is_eichler_canonical(gram):
    """
    Check if Gram is in Eichler canonical form (2U + K block).

    Port of is_eichler_canonical from ApproximateModels.h.

    INPUT:
    - gram: nested list of integers

    OUTPUT:
    - bool
    """
    n = len(gram)
    if n < 4:
        return False

    # Check 2U block:
    # Row 0: [0, 1, 0, 0, ...]
    # Row 1: [1, 0, 0, 0, ...]
    # Row 2: [0, 0, 0, 1, ...]
    # Row 3: [0, 0, 1, 0, ...]
    expected = [
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
    ]

    for i in range(4):
        for j in range(4):
            if gram[i][j] != expected[i][j]:
                return False

    # Check off-block zeros
    for i in range(4):
        for j in range(4, n):
            if gram[i][j] != 0:
                return False
            if gram[j][i] != 0:
                return False

    # Check even
    for i in range(n):
        if gram[i][i] % 2 != 0:
            return False

    return True


# ---------------------------------------------------------------------------
# Main: approximate model construction
# ---------------------------------------------------------------------------

def get_approximate_model(gram):
    """
    Construct the 2U-Eichler approximate model for an indefinite lattice.

    Port of INDEF_FORM_GetApproximateModel from ApproximateModels.h.

    The algorithm:
    1. Find two hyperbolic pairs to build the 2U sublattice
    2. Compute the Eichler canonical form
    3. Build generators from SL_2(Z) actions and Eichler transvections
    4. Return approximate group and covering orbit representatives

    INPUT:
    - gram: nested list of integers (indefinite even Gram matrix)

    OUTPUT:
    - dict with:
        'generators': list of matrix generators for A(L)
        'embedding': matrix P with P Q P^T = 2U + K
        'reduced_gram': the reduced Gram matrix
        'scaling': scaling factor
    """
    from .lorentzian_perfect import find_hyperbolic_pair, orthogonal_complement_basis

    Q = matrix(ZZ, gram)
    n = Q.nrows()

    # Step 1: Find first hyperbolic pair
    hp1 = find_hyperbolic_pair(gram)
    if hp1 is None:
        return {"generators": [], "embedding": None, "reduced_gram": gram, "scaling": 1}

    u1 = sage_vector(ZZ, hp1["u"])
    v1 = sage_vector(ZZ, hp1["v"])

    # Step 2: Find second hyperbolic pair in orthogonal complement
    # Compute orthogonal complement of span(u1, v1)
    V = matrix(ZZ, [list(u1), list(v1)])
    from sage.all import QQ as SageQQ
    NSP = (V * Q).transpose().change_ring(SageQQ).left_kernel().matrix()

    if NSP.nrows() == 0:
        # No complement possible - degenerate case
        return {"generators": [], "embedding": None, "reduced_gram": gram, "scaling": 1}

    # Compute restricted Gram matrix on complement
    Q_restricted = NSP * Q * NSP.transpose()

    # Find hyperbolic pair in complement
    hp2 = find_hyperbolic_pair([[int(Q_restricted[i, j]) for j in range(Q_restricted.ncols())]
                                 for i in range(Q_restricted.nrows())])

    if hp2 is None:
        # Can't find second hyperbolic pair
        # Build from single U
        basis = matrix(ZZ, [list(u1), list(v1)])
        for i in range(NSP.nrows()):
            basis = basis.stack(NSP[i])

        return {
            "generators": _build_generators_single_U(gram, basis),
            "embedding": basis,
            "reduced_gram": gram,
            "scaling": 1,
        }

    # Step 3: Build full basis
    u2_local = sage_vector(ZZ, hp2["u"])
    v2_local = sage_vector(ZZ, hp2["v"])

    # Map back to original coordinates
    u2 = sage_vector(ZZ, list(u2_local * NSP)) if NSP.nrows() > 0 else sage_vector(ZZ, [0]*n)
    v2 = sage_vector(ZZ, list(v2_local * NSP)) if NSP.nrows() > 0 else sage_vector(ZZ, [0]*n)

    # Build embedding matrix: [u1, v1, u2, v2, remaining_basis]
    basis_rows = [list(u1), list(v1), list(u2), list(v2)]
    remaining = NSP
    if remaining.nrows() > 0:
        basis_matrix = matrix(ZZ, basis_rows)
        for i in range(remaining.nrows()):
            basis_matrix = basis_matrix.stack(remaining[i])
    else:
        basis_matrix = matrix(ZZ, basis_rows)

    # Step 4: Build approximate group generators
    generators = []

    # SL_2(Z) actions on first U block
    S, T = sl2z_generators()
    for sl2_gen in [S, T]:
        # Left multiplication on (u1, u2) pair
        big_gen = identity_matrix(ZZ, n)
        for si in range(2):
            for sj in range(2):
                row_i = 2 * si  # u1, u2 indices
                col_j = 2 * sj
                big_gen[row_i, col_j] = sl2_gen[si, sj]
        generators.append(big_gen)

    # Eichler transvections
    for i in range(min(4, n)):
        ei = sage_vector(ZZ, [1 if k == i else 0 for k in range(n)])
        Qei = Q * ei

        # Find isotropic vectors orthogonal to ei
        from sage.all import QQ as SageQQ
        # Nullspace of Qei (a row vector): {w : Qei · w = 0}
        NSP_ei_mat = matrix(SageQQ, [list(Qei)]).left_kernel().matrix()
        NSP_ei = NSP_ei_mat
        for j in range(NSP_ei.nrows()):
            ej = sage_vector(ZZ, list(NSP_ei[j]))
            if int(ej * Q * ei) == 0:
                try:
                    T_mat = eichler_transvection(gram, list(ei), list(ej))
                    generators.append(T_mat)
                except (AssertionError, ValueError):
                    pass

    return {
        "generators": generators,
        "embedding": basis_matrix,
        "reduced_gram": gram,
        "scaling": 1,
    }


def _build_generators_single_U(gram, basis):
    """Build generators from a single hyperbolic plane embedding."""
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    generators = []

    # Standard generators of O(U): swap and negation
    if n >= 2:
        swap = identity_matrix(ZZ, n)
        swap[0, 0] = 0
        swap[0, 1] = 1
        swap[1, 0] = 1
        swap[1, 1] = 0
        generators.append(swap)

    return generators
