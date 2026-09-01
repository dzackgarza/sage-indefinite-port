"""
TASK-04-1: Lorentzian Perfect-Domain Base Engine and Component Groups.

All lattice operations use the preamble API.
Isometries are preamble morphisms (L.Hom(L) elements).
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, identity_matrix, vector as sage_vector
from . import make_lattice, matrix_to_morphism, unique_morphisms


def get_attack_scheme(gram):
    """
    Compute Witt index h via preamble L.signature_pair().

    OUTPUT: dict with 'h', 'mat', 'sign', 'lattice'
    """
    L = make_lattice(gram)
    n = L.rank()
    p, q = L.signature_pair()

    if p == 0 and q == 0:
        return {"h": 0, "mat": gram, "sign": 1, "lattice": L}

    h = min(p, q)
    if q < p:
        neg_gram = [[-int(L.b(L.module_generator(i), L.module_generator(j)))
                      for j in range(n)] for i in range(n)]
        return {"h": h, "mat": neg_gram, "sign": -1, "lattice": None}
    return {"h": h, "mat": gram, "sign": 1, "lattice": L}


def is_lorentzian(gram):
    L = make_lattice(gram)
    p, q = L.signature_pair()
    return min(p, q) == 1


def find_isotropic_vector(gram):
    """Find isotropic v with L.q(v) = 0."""
    L = make_lattice(gram)
    n = L.rank()
    bound = 10
    from itertools import product as iterproduct
    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = L(sage_vector(ZZ, coords))
        if L.q(v) == 0:
            return list(coords)
    return None


def find_hyperbolic_pair(gram):
    """Find isotropic u, v with L.b(u, v) = 1."""
    L = make_lattice(gram)
    n = L.rank()
    v1_list = find_isotropic_vector(gram)
    if v1_list is None:
        return None
    v1 = L(sage_vector(ZZ, v1_list))

    from itertools import product as iterproduct
    for coords in iterproduct(range(-20, 21), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v2 = L(sage_vector(ZZ, coords))
        if L.q(v2) != 0:
            continue
        scal = L.b(v1, v2)
        if scal == 1:
            return {"u": [int(c) for c in v1.to_list()],
                    "v": [int(c) for c in v2.to_list()], "scal": 1}
        if scal == -1:
            v2 = L(-sage_vector(ZZ, coords))
            return {"u": [int(c) for c in v1.to_list()],
                    "v": [int(c) for c in v2.to_list()], "scal": 1}
    return None


def reflection_morphism(gram, root):
    """
    Build an isometry morphism: reflection in root using preamble b(), q().

    Returns a preamble morphism element.
    """
    L = make_lattice(gram)
    n = L.rank()
    r = L(sage_vector(ZZ, root))
    q_r = L.q(r)
    if q_r == 0:
        return None

    r_coords = [int(root[i]) for i in range(n)]
    Qr = [int(L.b(r, L.module_generator(j))) for j in range(n)]

    # Build reflection matrix, then wrap as morphism
    R = identity_matrix(ZZ, n)
    for i in range(n):
        for j in range(n):
            R[i, j] -= 2 * r_coords[i] * Qr[j] // int(q_r)

    return matrix_to_morphism(L, R)


def lorentzian_generators_autom(gram):
    """
    Compute isometry morphisms generating O^Omega(L) for Lorentzian L.

    Returns list of preamble morphism elements.
    """
    L = make_lattice(gram)
    n = L.rank()
    if n <= 1:
        return [matrix_to_morphism(L, matrix(ZZ, [[-1]]))] if n == 1 else []

    # Find roots (norm ±2 vectors)
    roots = []
    from itertools import product as iterproduct
    for coords in iterproduct(range(-5, 6), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = L(sage_vector(ZZ, coords))
        if L.q(v) in (2, -2):
            roots.append(list(coords))

    # Build reflection morphisms
    mors = []
    for root in roots:
        mor = reflection_morphism(gram, root)
        if mor is not None and verify_isometry(L, mor):
            mors.append(mor)

    # Find additional generators from isotropic vector stabilizers
    for coords in iterproduct(range(-3, 4), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = L(sage_vector(ZZ, coords))
        if L.q(v) != 0:
            continue
        for root in roots:
            rv = L(sage_vector(ZZ, root))
            if L.b(rv, v) == 0:
                mor = reflection_morphism(gram, root)
                if mor is not None and verify_isometry(L, mor):
                    mors.append(mor)

    return unique_morphisms(mors, L)


def verify_isometry(L, mor):
    """Check morphism preserves bilinear form."""
    n = L.rank()
    for i in range(n):
        ei = L.module_generator(i)
        for j in range(n):
            ej = L.module_generator(j)
            if L.b(mor(ei), mor(ej)) != L.b(ei, ej):
                return False
    return True


def lorentzian_full_automorphism_group(gram):
    """Full O(L) = O^Omega(L) x <-I>."""
    L = make_lattice(gram)
    n = L.rank()
    mors = lorentzian_generators_autom(gram)

    # Add -I if not present
    neg_I = matrix_to_morphism(L, matrix(ZZ, [[-1 if i == j else 0 for j in range(n)] for i in range(n)]))
    id_mor = L.identity_morphism()
    seen = set()
    for m in mors:
        M = m.matrix()
        seen.add(tuple(M[i, j] for i in range(n) for j in range(n)))
    neg_key = tuple(neg_I.matrix()[i, j] for i in range(n) for j in range(n))
    if neg_key not in seen:
        mors.append(neg_I)

    return mors
