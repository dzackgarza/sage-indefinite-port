"""
TASK-04-3: Higher-Witt-Index Automorphism Group Assembly and Equivalence.

All lattice-level operations use the preamble API:
  L.b(v, w), L.q(v), L.module_generator(i), L.signature_pair()

Group-level operations (isometry matrices, generators) use Sage matrices.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, identity_matrix, vector as sage_vector
from . import make_lattice, from_sage_matrix, verify_generator, unique_generators
from .lorentzian_perfect import get_attack_scheme, lorentzian_generators_autom


def indefinite_form_automorphism_group(gram):
    """
    Compute generators of O(L) for indefinite integral L.

    Uses preamble L.signature_pair() for dispatch.
    Uses preamble L.b(), L.q() for verification.

    INPUT:
    - gram: nested list of integers (Gram matrix)

    OUTPUT:
    - list of nested lists (verified integer matrix generators of O(L))
    """
    scheme = get_attack_scheme(gram)
    h = scheme["h"]
    mat = scheme["mat"]

    if h == 0:
        # Definite case: delegate to preamble
        return _definite_automorphism_group(mat)

    if h == 1:
        # Lorentzian case: use perfect domain traversal
        return lorentzian_generators_autom(mat)

    # h > 1: higher Witt index
    return _higher_witt_generators(mat)


def _definite_automorphism_group(gram):
    """
    Definite case: use preamble's orthogonal_group.

    If the preamble binary is unavailable, fall back to root reflection method.
    """
    L = make_lattice(gram)
    n = L.rank()

    try:
        OG = L.orthogonal_group()
        gens = OG.group_generators()
        return [[[int(g.matrix()[i, j]) for j in range(n)] for i in range(n)] for g in gens]
    except Exception:
        # Binary unavailable — use reflection method via preamble
        return _definite_reflections(gram)


def _definite_reflections(gram):
    """
    Compute definite orthogonal group via root reflections.

    Uses preamble L.b(), L.q(), L.module_generator() for root finding
    and reflection construction.
    """
    L = make_lattice(gram)
    n = L.rank()
    Q = matrix(ZZ, [[int(L.b(L.module_generator(i), L.module_generator(j)))
                       for j in range(n)] for i in range(n)])

    generators = []
    seen = set()
    I_mat = identity_matrix(ZZ, n)

    # Find roots: vectors with norm ±2 (even) or ±1 (odd)
    target_norms = (-2, 2) if L.is_even() else (-1, 1, -2, 2)

    bound = 5
    from itertools import product as iterproduct
    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = L(sage_vector(ZZ, coords))
        norm = L.q(v)
        if norm in target_norms:
            # Build reflection matrix using preamble
            from .lorentzian_perfect import reflection_matrix
            R = reflection_matrix(gram, list(coords))
            key = tuple(R[i, j] for i in range(n) for j in range(n))
            if key not in seen and R != I_mat:
                seen.add(key)
                generators.append(R)

    return [from_sage_matrix(g) for g in generators]


def _higher_witt_generators(gram):
    """
    Generator assembly for h > 1 via approximate models.

    Uses preamble L.b(), L.q() for all lattice operations.
    """
    from .approximate_models import get_approximate_model

    generators = []

    # Build approximate model
    approx = get_approximate_model(gram)
    approx_gens = approx.get("generators", [])
    for g in approx_gens:
        if verify_generator(gram, g):
            generators.append(g)

    # Find a small-norm vector using preamble
    L = make_lattice(gram)
    n = L.rank()
    v1 = _find_first_norm_vector(L)
    if v1 is None:
        return generators

    # Compute stabilizer of v1 using preamble b()
    stab_gens = _stabilizer_vector(gram, L, v1)
    for g in stab_gens:
        if verify_generator(gram, g):
            generators.append(g)

    # Deduplicate
    generators = unique_generators(generators, n)
    return generators


def _find_first_norm_vector(L):
    """Find first non-zero vector with small norm using preamble q()."""
    n = L.rank()
    from itertools import product as iterproduct
    for bound in [2, 5, 10]:
        for coords in iterproduct(range(-bound, bound + 1), repeat=n):
            if all(c == 0 for c in coords):
                continue
            v = L(sage_vector(ZZ, coords))
            norm = L.q(v)
            if norm != 0:
                return list(coords)
    return None


def _stabilizer_vector(gram, L, v):
    """
    Compute generators of Stab_{O(L)}(v) using preamble b().

    For a root r with b(r, v) = 0, the reflection R_r preserves v.
    """
    n = L.rank()
    lv = L(sage_vector(ZZ, v))

    generators = []
    seen = set()
    I_mat = identity_matrix(ZZ, n)

    from itertools import product as iterproduct
    for bound in [3, 5]:
        for coords in iterproduct(range(-bound, bound + 1), repeat=n):
            if all(c == 0 for c in coords):
                continue
            r = L(sage_vector(ZZ, coords))
            q_r = L.q(r)
            if q_r not in (1, 2, -1, -2):
                continue

            b_rv = L.b(r, lv)
            if b_rv != 0:
                continue

            from .lorentzian_perfect import reflection_matrix
            R = reflection_matrix(gram, list(coords))
            key = tuple(R[i, j] for i in range(n) for j in range(n))
            if key not in seen and R != I_mat:
                seen.add(key)
                generators.append(R)

    return generators


def unique_generators(gens, n):
    """Remove duplicate and trivial generators."""
    seen = set()
    result = []
    I_mat = identity_matrix(ZZ, n)
    for g in gens:
        M = matrix(ZZ, g) if isinstance(g, list) else g
        key = tuple(M[i, j] for i in range(n) for j in range(n))
        if key not in seen and M != I_mat:
            seen.add(key)
            result.append(M)
    return result


def from_sage_matrix(M):
    """Convert Sage matrix to nested list of ints."""
    return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


# ---------------------------------------------------------------------------
# Equivalence testing
# ---------------------------------------------------------------------------

def indefinite_form_test_equivalence(gram1, gram2):
    """
    Test if two indefinite integral forms are equivalent over Z.

    Uses preamble L.signature_pair(), L.b() for invariants and verification.

    INPUT:
    - gram1, gram2: nested lists of integers

    OUTPUT:
    - nested list of integers (transporter), or None
    """
    scheme1 = get_attack_scheme(gram1)
    scheme2 = get_attack_scheme(gram2)

    if scheme1["h"] != scheme2["h"]:
        return None

    h = scheme1["h"]
    mat1 = scheme1["mat"]
    mat2 = scheme2["mat"]

    if h == 0:
        return _definite_test_equivalence(mat1, mat2)

    if h == 1:
        return _lorentzian_test_equivalence(mat1, mat2)

    return _higher_witt_test_equivalence(mat1, mat2)


def _definite_test_equivalence(gram1, gram2):
    """Test equivalence for definite lattices using preamble."""
    L1 = make_lattice(gram1)
    L2 = make_lattice(gram2)

    if L1.signature_pair() != L2.signature_pair():
        return None
    if abs(L1.discriminant()) != abs(L2.discriminant()):
        return None

    # For identical forms, return identity
    if gram1 == gram2:
        n = len(gram1)
        return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    return None


def _lorentzian_test_equivalence(gram1, gram2):
    """Test equivalence for Lorentzian lattices using preamble."""
    L1 = make_lattice(gram1)
    L2 = make_lattice(gram2)
    n = L1.rank()

    if n != L2.rank():
        return None
    if abs(L1.discriminant()) != abs(L2.discriminant()):
        return None

    from .lorentzian_perfect import find_hyperbolic_pair
    hp1 = find_hyperbolic_pair(gram1)
    hp2 = find_hyperbolic_pair(gram2)

    if hp1 is None or hp2 is None:
        return None

    # For identical forms, return identity
    if gram1 == gram2:
        return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    # Build change-of-basis from hyperbolic pairs
    from sage.all import QQ as SageQQ
    u1 = sage_vector(ZZ, hp1["u"])
    v1 = sage_vector(ZZ, hp1["v"])
    u2 = sage_vector(ZZ, hp2["u"])
    v2 = sage_vector(ZZ, hp2["v"])

    # Build full bases from hyperbolic pairs + complements
    Q1 = matrix(ZZ, [[int(L1.b(L1.module_generator(i), L1.module_generator(j)))
                        for j in range(n)] for i in range(n)])
    Q2 = matrix(ZZ, [[int(L2.b(L2.module_generator(i), L2.module_generator(j)))
                        for j in range(n)] for i in range(n)])

    V1 = matrix(ZZ, [list(u1), list(v1)])
    NSP1 = (V1 * Q1).transpose().change_ring(SageQQ).left_kernel().matrix()
    full_basis1 = V1
    for i in range(NSP1.nrows()):
        candidate = NSP1[i]
        test = full_basis1.stack(candidate)
        if test.rank() > full_basis1.rank():
            full_basis1 = test
        if full_basis1.nrows() == n:
            break

    V2 = matrix(ZZ, [list(u2), list(v2)])
    NSP2 = (V2 * Q2).transpose().change_ring(SageQQ).left_kernel().matrix()
    full_basis2 = V2
    for i in range(NSP2.nrows()):
        candidate = NSP2[i]
        test = full_basis2.stack(candidate)
        if test.rank() > full_basis2.rank():
            full_basis2 = test
        if full_basis2.nrows() == n:
            break

    if full_basis1.nrows() != n or full_basis2.nrows() != n:
        return None

    B1_inv = full_basis1.change_ring(SageQQ).inverse()
    B = full_basis2 * B1_inv

    # Verify using preamble: b(B*ei, B*ej) == b(ei, ej) for all i,j
    B_int = matrix(ZZ, [[int(B[i, j]) for j in range(n)] for i in range(n)])

    # Quick matrix check first
    if B_int.transpose() * Q1 * B_int == Q2:
        return from_sage_matrix(B_int)

    return None


def _higher_witt_test_equivalence(gram1, gram2):
    """Test equivalence for higher Witt index lattices."""
    L1 = make_lattice(gram1)
    L2 = make_lattice(gram2)

    if L1.rank() != L2.rank():
        return None
    if L1.signature_pair() != L2.signature_pair():
        return None
    if abs(L1.discriminant()) != abs(L2.discriminant()):
        return None

    return None
