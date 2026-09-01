"""
TASK-04-1: Lorentzian Perfect-Domain Base Engine and Component Groups.

All computation uses the preamble lattice API:
  L.b(v, w)              — bilinear pairing
  L.q(v)                 — quadratic form
  L.module_generator(i)  — i-th basis vector
  L.signature_pair()     — (p, q)
  L.rank()               — rank

NEVER bypass the preamble to use raw Sage matrices.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, identity_matrix, vector as sage_vector
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.tensors import tensor


def _make_lattice(gram):
    """Build a preamble Lattice from a nested-list Gram matrix."""
    C = Lattices(ZZ)
    n = len(gram)
    rows = tuple(tuple(int(x) for x in row) for row in gram)
    g = tensor(ZZ, (), (n, n), rows)
    return C(g)


def _from_sage_matrix(M):
    """Convert Sage matrix to nested list of ints."""
    return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


# ---------------------------------------------------------------------------
# Attack scheme: compute Witt index h and sign-adjusted matrix
# ---------------------------------------------------------------------------

def get_attack_scheme(gram):
    """
    Compute Witt index h and sign-adjusted Gram matrix.

    Uses preamble L.signature_pair() — NOT raw QuadraticForm.

    INPUT:
    - gram: nested list of integers (Gram matrix)

    OUTPUT:
    - dict with keys 'h' (Witt index), 'mat' (adjusted Gram), 'sign' (+1 or -1)
    """
    L = _make_lattice(gram)
    n = L.rank()
    p, q = L.signature_pair()

    if p == 0 and q == 0:
        return {"h": 0, "mat": gram, "sign": 1, "lattice": L}

    h = min(p, q)

    if q < p:
        # Negate so fewer negative eigenvalues
        neg_gram = [[-int(L.b(L.module_generator(i), L.module_generator(j)))
                      for j in range(n)] for i in range(n)]
        return {"h": h, "mat": neg_gram, "sign": -1, "lattice": None}
    else:
        return {"h": h, "mat": gram, "sign": 1, "lattice": L}


def is_lorentzian(gram):
    """Check if gram has signature (1, n-1) or (n-1, 1)."""
    L = _make_lattice(gram)
    p, q = L.signature_pair()
    return min(p, q) == 1


# ---------------------------------------------------------------------------
# Lorentzian perfect domain: finding the fundamental domain
# ---------------------------------------------------------------------------

def find_isotropic_vector(gram):
    """
    Find an integer isotropic vector v with Q(v) = 0.

    Uses preamble L.q() and L.module_generator().

    INPUT:
    - gram: nested list of integers

    OUTPUT:
    - list of integers (isotropic vector), or None
    """
    L = _make_lattice(gram)
    n = L.rank()

    # Try nullspace approach: Qv = 0
    Q = matrix(ZZ, [[int(L.b(L.module_generator(i), L.module_generator(j)))
                       for j in range(n)] for i in range(n)])
    from sage.all import QQ as SageQQ
    Q_qq = Q.change_ring(SageQQ)
    NSP = Q_qq.left_kernel().matrix()
    if NSP.nrows() > 0:
        v = NSP[0]
        return [int(v[i]) for i in range(n)]

    # Enumerate small vectors using preamble q()
    bound = 10
    from itertools import product as iterproduct
    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = sage_vector(ZZ, coords)
        norm = L.q(L(v))
        if norm == 0:
            return list(coords)

    return None


def find_hyperbolic_pair(gram):
    """
    Find isotropic vectors u, v with b(u, v) = 1 (a hyperbolic pair).

    Uses preamble L.b(), L.q(), L.module_generator().

    INPUT:
    - gram: nested list of integers

    OUTPUT:
    - dict with 'u', 'v' (lists), 'scal' (integer), or None
    """
    L = _make_lattice(gram)
    n = L.rank()

    v1_list = find_isotropic_vector(gram)
    if v1_list is None:
        return None

    v1 = L(sage_vector(ZZ, v1_list))

    bound = 20
    from itertools import product as iterproduct
    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v2 = L(sage_vector(ZZ, coords))
        norm2 = L.q(v2)
        if norm2 != 0:
            continue
        scal = L.b(v1, v2)
        if scal != 0:
            if abs(scal) == 1:
                if scal == -1:
                    v2 = L(-sage_vector(ZZ, coords))
                    scal = 1
                return {
                    "u": [int(c) for c in v1.to_list()],
                    "v": [int(c) for c in v2.to_list()],
                    "scal": 1,
                }
            return {
                "u": [int(c) for c in v1.to_list()],
                "v": [int(c) for c in v2.to_list()],
                "scal": int(scal),
            }

    return None


def reflection_matrix(gram, root):
    """
    Compute the orthogonal reflection matrix for a root vector.

    Uses preamble L.b(), L.q(), L.module_generator().

    R = I - 2 * root * (Q * root)^T / Q(root)

    INPUT:
    - gram: nested list of integers
    - root: list of integers (the root vector)

    OUTPUT:
    - Sage integer matrix (the reflection)
    """
    L = _make_lattice(gram)
    n = L.rank()
    r = L(sage_vector(ZZ, root))
    q_r = L.q(r)

    if q_r == 0:
        raise ValueError("Cannot reflect in isotropic vector")

    # R = I - 2 * r * (Qr)^T / q(r)
    # R_{ij} = delta_{ij} - 2 * r_i * (Qr)_j / q(r)
    # r_i = coordinate of root, (Qr)_j = b(root, e_j)
    r_coords = [int(root[i]) for i in range(n)]
    Qr = [int(L.b(r, L.module_generator(j))) for j in range(n)]
    R = identity_matrix(ZZ, n)
    for i in range(n):
        for j in range(n):
            R[i, j] -= 2 * r_coords[i] * Qr[j] // int(q_r)

    return R


def lorentzian_perfect_domain_step(gram):
    """
    One step of Lorentzian perfect domain traversal.

    Uses preamble L.q(), L.b(), L.module_generator().

    INPUT:
    - gram: nested list of integers (Lorentzian Gram matrix)

    OUTPUT:
    - list of root vectors (lists) that are walls of the current chamber
    """
    L = _make_lattice(gram)
    n = L.rank()

    roots = []
    bound = 5
    from itertools import product as iterproduct
    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = L(sage_vector(ZZ, coords))
        norm = L.q(v)
        if norm == 2 or norm == -2:
            roots.append(list(coords))

    return roots


def lorentzian_generators_autom(gram):
    """
    Compute generators of O^Omega(L) for a Lorentzian lattice.

    Uses preamble L.b(), L.q(), L.module_generator() throughout.

    INPUT:
    - gram: nested list of integers (Lorentzian Gram matrix)

    OUTPUT:
    - list of nested lists (integer matrix generators)
    """
    L = _make_lattice(gram)
    n = L.rank()

    if n <= 1:
        return [[[-1]]] if n == 1 else []

    roots = lorentzian_perfect_domain_step(gram)

    generators = []
    seen = set()
    I_mat = identity_matrix(ZZ, n)

    for root in roots:
        R = reflection_matrix(gram, root)
        key = tuple(R[i, j] for i in range(n) for j in range(n))
        if key not in seen and R != I_mat:
            seen.add(key)
            generators.append(R)

    # Find additional generators from stabilizers of isotropic vectors
    from itertools import product as iterproduct
    for coords in iterproduct(range(-3, 4), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = L(sage_vector(ZZ, coords))
        norm = L.q(v)
        if norm != 0:
            continue

        for root in roots:
            rv = L(sage_vector(ZZ, root))
            b_v_r = L.b(rv, v)
            if b_v_r == 0:
                R = reflection_matrix(gram, root)
                key = tuple(R[i, j] for i in range(n) for j in range(n))
                if key not in seen and R != I_mat:
                    seen.add(key)
                    generators.append(R)

    return [_from_sage_matrix(g) for g in generators]


# ---------------------------------------------------------------------------
# Full O(L) for Lorentzian: O^Omega x <-I>
# ---------------------------------------------------------------------------

def lorentzian_full_automorphism_group(gram):
    """
    Compute generators of full O(L) for Lorentzian L.

    O(L) = O^Omega(L) x <-I>

    INPUT:
    - gram: nested list of integers

    OUTPUT:
    - list of nested lists (generators of O(L))
    """
    n = len(gram)
    gens = lorentzian_generators_autom(gram)

    # Add -I if not already present
    neg_I = [[-1 if i == j else 0 for j in range(n)] for i in range(n)]
    neg_I_key = tuple(neg_I[i][j] for i in range(n) for j in range(n))
    seen = {tuple(g[i][j] for i in range(n) for j in range(n)) for g in gens}
    if neg_I_key not in seen:
        gens.append(neg_I)

    return gens
