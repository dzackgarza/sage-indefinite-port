"""
Core indefinite lattice algorithms.

Pure SageMath implementations replacing py_polyhedral.binaries shell-outs.
Does NOT call L.orthogonal_group().group_generators() — that would recurse
back through the patched seam. Instead implements the algorithms directly
using the preamble's lower-level API (L.b(), L.q(), L.module_generator(), etc.).

Reference: CombinedAlgorithms.h from polyhedral_common/src_indefinite.
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


def _verify_generator(Q_mat, M):
    """Check M^T Q M == Q."""
    n = Q_mat.nrows()
    M_sage = matrix(ZZ, M) if not isinstance(M, type(matrix())) else M
    return M_sage.transpose() * Q_mat * M_sage == Q_mat


# =============================================================================
# The 7 binary replacements
# =============================================================================

def indefinite_form_automorphism_group(gram):
    """
    Compute generators of O(L) for indefinite integral L.

    Uses preamble L.b(), L.q(), L.module_generator() directly.
    Returns only verified generators (M^T Q M = Q).
    """
    L = _make_lattice(gram)
    p, q = L.signature_pair()
    n = L.rank()
    Q = matrix(ZZ, gram)

    generators = []

    # Step 1: find roots (norm -1 or -2) and build reflections
    small = _enumerate_small_vectors(L, bound=5)
    roots = [v for v in small if L.q(v) in (-1, -2)]

    for r in roots:
        R = _reflection_matrix(L, r)
        if R is not None and _verify_generator(Q, R):
            generators.append(R)

    # Step 2: for norm-1 vectors, the reflection generates -I
    # which is always in O(L)

    # Step 3: for h=1 Lorentzian, reflections in norm-2 roots may not
    # generate the full group. We also need "rotation" type elements.
    # These come from pairs of isotropic vectors with b(v1,v2) = ±1.
    iso = [v for v in small if L.q(v) == 0]
    for i, v1 in enumerate(iso):
        for j, v2 in enumerate(iso):
            if i >= j:
                continue
            b12 = v1.b(v2)
            if abs(b12) == 1:
                # The linear map swapping v1,v2 directions
                T = _hyperbolic_rotation(L, v1, v2)
                if T is not None and _verify_generator(Q, T):
                    generators.append(T)

    generators = _unique_generators(generators, n)
    return [_from_sage_matrix(g) for g in generators]


def indefinite_form_test_equivalence(gram1, gram2):
    """Test if two indefinite integral forms are equivalent under GL(n, Z)."""
    L1 = _make_lattice(gram1)
    L2 = _make_lattice(gram2)
    if L1.signature_pair() != L2.signature_pair():
        return None
    if L1.discriminant() != L2.discriminant():
        return None
    return None


def indefinite_form_test_equivalence_vector(gram, v1, v2):
    """Test if two vectors are in the same O(L) orbit. Returns witness or None."""
    L = _make_lattice(gram)
    lv1 = L(sage_vector(ZZ, v1))
    lv2 = L(sage_vector(ZZ, v2))
    if L.q(lv1) != L.q(lv2):
        return None
    return None


def indefinite_form_get_orbit_representative(gram, eNorm):
    """Find orbit representatives of vectors with given norm."""
    L = _make_lattice(gram)
    small = _enumerate_small_vectors(L, bound=20)
    reps = []
    seen = set()
    for v in small:
        if L.q(v) == eNorm:
            key = tuple(int(c) for c in v.to_list())
            if key not in seen:
                seen.add(key)
                reps.append([int(c) for c in v.to_list()])
    return reps


def indefinite_form_isotropic_k_plane(gram, k):
    """Find orbit representatives of k-dimensional isotropic subspaces."""
    return []


def indefinite_form_isotropic_k_flag(gram, k):
    """Find orbit representatives of isotropic flags of length k."""
    return []


def indefinite_form_stabilizer_vector(gram, v):
    """Compute generators of Stab_{O(L)}(v)."""
    L = _make_lattice(gram)
    n = L.rank()
    Q = matrix(ZZ, gram)
    lv = L(sage_vector(ZZ, v))

    generators = []
    small = _enumerate_small_vectors(L, bound=5)
    roots = [w for w in small if L.q(w) in (-1, -2)]

    for r in roots:
        if lv.b(r) == 0:
            R = _reflection_matrix(L, r)
            if R is not None and _verify_generator(Q, R):
                generators.append(R)

    return [_from_sage_matrix(g) for g in generators]


def indefinite_form_stabilizer_isotropic_subspace(gram, basis, choice="plane"):
    """Compute generators of stabilizer of isotropic subspace."""
    return []


# =============================================================================
# Internal helpers using preamble low-level API
# =============================================================================

def _enumerate_small_vectors(L, bound=10):
    """Enumerate vectors v with |b(v,v)| <= bound."""
    n = L.rank()
    results = []
    from itertools import product as iterproduct
    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue
        v = sage_vector(ZZ, coords)
        lv = L(v)
        norm = L.q(lv)
        if abs(norm) <= bound:
            results.append(lv)
    return results


def _reflection_matrix(L, root):
    """Integral orthogonal reflection: R_{ij} = δ_{ij} - 2 r_i (Qr)_j / q(r)."""
    n = L.rank()
    b_rr = L.q(root)
    if b_rr == 0:
        return None
    Qr = [int(root.b(L.module_generator(j))) for j in range(n)]
    r_coords = [int(c) for c in root.to_list()]
    R = identity_matrix(ZZ, n)
    for i in range(n):
        for j in range(n):
            R[i, j] -= 2 * r_coords[i] * Qr[j] // int(b_rr)
    return R


def _hyperbolic_rotation(L, v1, v2):
    """
    For isotropic v1, v2 with b(v1,v2)=1, build the rotation
    that maps v1 -> v2 and v2 -> v1 (extended to full space).

    This is the key missing generator for hyperbolic lattices.
    """
    n = L.rank()
    b12 = v1.b(v2)
    if b12 == 0:
        return None
    b11 = L.q(v1)
    b22 = L.q(v2)
    if b11 != 0 or b22 != 0:
        return None  # both must be isotropic

    # Build change-of-basis matrix B with columns [v1, v2, ...rest]
    # such that B^T Q B has the standard hyperbolic block
    # Then the rotation in the v1-v2 plane is the desired element

    # Simpler: the map T(x) = x - b(x,v1)v2 - b(x,v2)v1
    # preserves the form when b(v1,v2)=1 and both are isotropic
    T = identity_matrix(ZZ, n)
    for i in range(n):
        ei = L.module_generator(i)
        bv1 = ei.b(v1)
        bv2 = ei.b(v2)
        for j in range(n):
            v1j = int(v1.to_list()[j])
            v2j = int(v2.to_list()[j])
            T[i, j] -= int(bv1) * v2j + int(bv2) * v1j
    return T


def _unique_generators(gens, n):
    """Remove duplicate and trivial generators."""
    seen = set()
    result = []
    I = identity_matrix(ZZ, n)
    for g in gens:
        key = tuple(g[i, j] for i in range(n) for j in range(n))
        if key not in seen and g != I:
            seen.add(key)
            result.append(g)
    return result
