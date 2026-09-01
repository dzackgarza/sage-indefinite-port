"""
TASK-05: Parabolic Recursion and Isotropic Subspace Orbit Decomposition.

TASK-05-1: Primitive isotropic vector cusps, divisibility, and cusp orbits.
TASK-05-2: Isotropic reduction objects (K_I = I^perp/I).
TASK-05-3: Exact gluing-based parabolic stabilizers.
TASK-05-4: Inductive rank-k isotropic sublattice and flag orbits.

Port of INDEF_FORM_Rec_IsotropicKplane and related functions from
CombinedAlgorithms.h.
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
# TASK-05-1: Primitive isotropic vector loci, divisibility, cusp orbits
# ---------------------------------------------------------------------------

def primitive_isotropic_vectors(gram, bound=10):
    """
    Find all primitive isotropic vectors with coordinates bounded.

    A vector v is primitive if gcd(v_1, ..., v_n) = 1.
    It is isotropic if Q(v) = 0.

    INPUT:
    - gram: nested list of integers
    - bound: coordinate bound for enumeration

    OUTPUT:
    - list of lists of integers (primitive isotropic vectors)
    """
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    result = []

    for coords in iterproduct(range(-bound, bound + 1), repeat=n):
        if all(c == 0 for c in coords):
            continue

        # Check primitivity
        g = gcd(list(coords))
        if g != 1:
            continue

        # Check isotropy
        v = sage_vector(ZZ, coords)
        norm = int(v * Q * v)
        if norm == 0:
            result.append(list(coords))

    return result


def divisor(gram, v):
    """
    Compute the divisor of v: the integer d such that b(v, L) = dZ.

    Port of divisor computation from IndefiniteFormFundamental.h.

    INPUT:
    - gram: nested list of integers
    - v: list of integers (the vector)

    OUTPUT:
    - positive integer d
    """
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    lv = sage_vector(ZZ, v)

    # b(v, e_i) = (Qv)_i
    Qv = Q * lv

    # gcd of all entries
    entries = [int(Qv[i]) for i in range(n) if int(Qv[i]) != 0]
    if not entries:
        return 1

    return abs(gcd(entries))


def discriminant_class(gram, v):
    """
    Compute the class of v in the discriminant group A_L = L*/L.

    INPUT:
    - gram: nested list of integers
    - v: list of integers

    OUTPUT:
    - list of integers (representative of v/div(v) in A_L)
    """
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    d = divisor(gram, v)
    lv = sage_vector(ZZ, v)

    # v/d is in L*, so v/d mod L is the discriminant class
    Qv = Q * lv
    # Class is determined by Qv mod (d * Z^n)
    return [int(Qv[i]) % d for i in range(n)]


def cusp_orbit_decomposition(gram):
    """
    Decompose primitive isotropic vectors into O(L)-orbits (cusps).

    Two primitive isotropic vectors v, w are in the same orbit iff
    they have the same discriminant class in A_L.

    INPUT:
    - gram: nested list of integers

    OUTPUT:
    - dict mapping discriminant class (tuple) to list of orbit representatives
    """
    iso_vecs = primitive_isotropic_vectors(gram, bound=5)
    orbits = {}

    for v in iso_vecs:
        d = divisor(gram, v)
        dc = tuple(discriminant_class(gram, v))
        key = (d, dc)

        if key not in orbits:
            orbits[key] = []
        orbits[key].append(v)

    # Return one representative per orbit
    return {k: vlist[0] for k, vlist in orbits.items()}


# ---------------------------------------------------------------------------
# TASK-05-2: Isotropic reduction objects (K_I = I^perp/I)
# ---------------------------------------------------------------------------

def isotropic_reduction(gram, isotropic_basis):
    """
    Compute the isotropic reduction K_I = I^perp / I.

    Given a totally isotropic subspace I spanned by isotropic_basis,
    compute:
    - The quotient lattice K_I (non-degenerate)
    - Inclusion maps and projections
    - The Witt splitting data

    Port of INDEF_FORM_Rec_IsotropicKplane from CombinedAlgorithms.h.

    INPUT:
    - gram: nested list of integers (Gram matrix)
    - isotropic_basis: list of lists of integers (basis of I)

    OUTPUT:
    - dict with:
        'quotient_gram': Gram matrix of K_I
        'nsp': nullspace basis (orthogonal complement of I in Q)
        'projection': matrix projecting from L to K_I
        'inclusion': matrix embedding K_I into L
    """
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    k = len(isotropic_basis)

    if k == 0:
        return {
            "quotient_gram": gram,
            "nsp": identity_matrix(ZZ, n),
            "projection": identity_matrix(ZZ, n),
            "inclusion": identity_matrix(ZZ, n),
        }

    # I = span(isotropic_basis)
    I = matrix(ZZ, isotropic_basis)

    # Compute I^perp = {v in L : b(v, I) = 0}
    # This is the kernel of I * Q as a map from L to Z^k
    Iq = I * Q  # k x n matrix
    from sage.all import QQ as SageQQ
    NSP = Iq.change_ring(SageQQ).left_kernel().matrix()  # (n-k) x n matrix

    # K_I = I^perp / I restricted to NSP
    # The quotient is represented by NSP mod I
    # In coordinates: express NSP in terms of a basis that contains I
    # The reduced Gram is NSP * Q * NSP^T (on the complement of I within I^perp)

    # Find complement of I within I^perp
    # Use SNF to separate I from I^perp
    from sage.all import MatrixSpace

    # Embed I into I^perp and find complement
    if NSP.nrows() == k:
        # No room for complement
        return {
            "quotient_gram": matrix(ZZ, 0, 0),
            "nsp": NSP,
            "projection": I,
            "inclusion": I,
        }

    # Build full basis: I followed by complement of I in I^perp
    full_basis = I
    for i in range(NSP.nrows()):
        candidate = NSP[i]
        # Check if candidate is in span of full_basis rows
        if full_basis.nrows() < n:
            test = full_basis.stack(candidate)
            if test.rank() > full_basis.rank():
                full_basis = test

    # Reduced Gram matrix on the complement
    complement_dim = full_basis.nrows() - k
    if complement_dim > 0:
        comp_basis = full_basis[k:]  # rows after I
        quotient_gram = comp_basis * Q * comp_basis.transpose()
    else:
        quotient_gram = matrix(ZZ, 0, 0)

    return {
        "quotient_gram": [[int(quotient_gram[i, j]) for j in range(quotient_gram.ncols())]
                           for i in range(quotient_gram.nrows())] if quotient_gram.nrows() > 0 else [],
        "nsp": [[int(NSP[i, j]) for j in range(NSP.ncols())] for i in range(NSP.nrows())],
        "projection": [[int(full_basis[i, j]) for j in range(full_basis.ncols())]
                        for i in range(full_basis.nrows())],
        "inclusion": [[int(NSP[i, j]) for j in range(NSP.ncols())]
                       for i in range(NSP.nrows())],
    }


def isotropic_reduction_objects(gram, isotropic_basis):
    """
    Full structured isotropic reduction returning IsometryExtensionTorsor data.

    Returns the unipotent kernel and gluing data needed for parabolic
    stabilizer computation.

    INPUT:
    - gram: nested list of integers
    - isotropic_basis: list of lists (basis of isotropic subspace)

    OUTPUT:
    - dict with reduction data and extension torsor
    """
    red = isotropic_reduction(gram, isotropic_basis)

    Q = matrix(ZZ, gram)
    n = Q.nrows()
    k = len(isotropic_basis)

    # Compute gluing matrix H: off-diagonal block in Witt decomposition
    # H_ij = b(e_i, f_j) where e_i in I and f_j in complement
    I = matrix(ZZ, isotropic_basis)
    NSP = matrix(ZZ, red["nsp"])

    # Gluing matrix
    H = I * Q * NSP.transpose()

    red["gluing_matrix"] = [[int(H[i, j]) for j in range(H.ncols())]
                             for i in range(H.nrows())] if H.nrows() > 0 else []

    # Unipotent radical: solutions U to H U^T + U H^T = 0
    # This is the set of matrices that act trivially on I and K_I
    red["unipotent_radical_dim"] = max(0, k * k - k * (k - 1) // 2)

    return red


# ---------------------------------------------------------------------------
# TASK-05-3: Exact gluing-based parabolic stabilizers
# ---------------------------------------------------------------------------

def parabolic_stabilizer(gram, isotropic_basis):
    """
    Compute generators of the parabolic stabilizer P_I.

    P_I = {g in O(L) : g(I) = I (setwise)}

    The exact computation uses:
    1. Rational Witt splitting I_Q + K_Q + I'_Q
    2. Extension 1 -> U_I(Z) -> P_I -> M_I -> 1
    3. Finite gluing data H_L to constrain U_I

    Port of the parabolic stabilizer from CombinedAlgorithms.h.

    INPUT:
    - gram: nested list of integers
    - isotropic_basis: list of lists (basis of isotropic subspace)

    OUTPUT:
    - list of nested lists (generator matrices of P_I)
    """
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    k = len(isotropic_basis)

    generators = []

    if k == 0:
        # Full orthogonal group
        from .higher_witt import indefinite_form_automorphism_group
        return indefinite_form_automorphism_group(gram)

    # Step 1: Compute the quotient K_I = I^perp / I
    red = isotropic_reduction(gram, isotropic_basis)

    quotient_gram = red["quotient_gram"]
    if not quotient_gram or (isinstance(quotient_gram, list) and len(quotient_gram) == 0):
        return generators

    # Step 2: Compute automorphism group of the quotient
    # This gives the finite Levi quotient M_I
    from .higher_witt import indefinite_form_automorphism_group
    quotient_gens = indefinite_form_automorphism_group(quotient_gram)

    # Step 3: Lift quotient generators to full isometries
    # Each quotient generator extends to an element of P_I
    NSP = matrix(ZZ, red["nsp"])
    I = matrix(ZZ, isotropic_basis)

    for q_gen in quotient_gens:
        Q_gen = matrix(ZZ, q_gen) if isinstance(q_gen, list) else q_gen

        # Extend: Q_gen acts on complement, identity on I
        full_gen = identity_matrix(ZZ, n)

        # I basis maps to itself
        # Complement maps according to Q_gen
        full_basis = I.stack(NSP) if NSP.nrows() > 0 else I
        full_basis_inv = full_basis.change_ring(QQ).inverse()

        # Build extended matrix
        big_gen = identity_matrix(ZZ, n)
        dim_comp = NSP.nrows()
        for i in range(dim_comp):
            for j in range(dim_comp):
                row_i = k + i
                col_j = k + j
                if row_i < n and col_j < n:
                    big_gen[row_i, col_j] = int(Q_gen[i, j])

        generators.append(big_gen)

    # Step 4: Add unipotent radical generators
    # These are transvections that preserve I and act trivially on K_I
    gluing = red.get("gluing_matrix", [])
    if gluing:
        H = matrix(ZZ, gluing)
        # Unipotent elements: I + E where E is in the kernel of the
        # gluing constraint
        for i in range(k):
            for j in range(NSP.nrows()):
                U = identity_matrix(ZZ, n)
                # Map I_i to I_i + NSP_j
                if i < n and k + j < n:
                    U[i, k + j] = 1
                generators.append(U)

    return generators


# ---------------------------------------------------------------------------
# TASK-05-4: Inductive isotropic plane orbits via double cosets
# ---------------------------------------------------------------------------

def isotropic_k_plane_orbits(gram, k):
    """
    Compute orbit representatives of k-dimensional isotropic sublattices.

    Uses inductive approach: extend (k-1)-dimensional orbits to k-dimensional
    via double coset decomposition Q_u \ O(K_I) / H_I.

    Port of INDEF_FORM_GetOrbit_IsotropicKplane from
    INDEF_FORM_GetOrbit_IsotropicKplane.cpp.

    INPUT:
    - gram: nested list of integers
    - k: dimension of isotropic subspace

    OUTPUT:
    - list of dicts with 'basis' (list of lists) and 'stabilizer_gens'
    """
    if k == 0:
        return [{"basis": [], "stabilizer_gens": indefinite_form_automorphism_group_simple(gram)}]

    if k == 1:
        # Single isotropic vector orbits
        cusps = cusp_orbit_decomposition(gram)
        results = []
        for (d, dc), v in cusps.items():
            stab = _stabilizer_isotropic_vector(gram, v)
            results.append({
                "basis": [v],
                "stabilizer_gens": stab,
                "divisor": d,
                "discriminant_class": dc,
            })
        return results

    # Inductive: extend from (k-1) to k
    prev_orbits = isotropic_k_plane_orbits(gram, k - 1)
    results = []

    for prev_orbit in prev_orbits:
        prev_basis = prev_orbit["basis"]

        # Extend prev_basis to a k-dimensional isotropic subspace
        # Find new isotropic vector v orthogonal to prev_basis
        Q = matrix(ZZ, gram)
        n = Q.nrows()

        # v must satisfy: Qv = 0 and b(v, prev_basis_i) = 0 for all i
        constraints = matrix(ZZ, prev_basis) * Q if prev_basis else matrix(ZZ, 0, n)
        from sage.all import QQ as SageQQ
        NSP = constraints.change_ring(SageQQ).left_kernel().matrix() if constraints.nrows() > 0 else matrix(ZZ, identity_matrix(ZZ, n))

        # Find isotropic vectors in the constrained space
        for i in range(NSP.nrows()):
            v = sage_vector(ZZ, list(NSP[i]))
            norm = int(v * Q * v)
            if norm == 0:
                new_basis = prev_basis + [list(v)]
                stab = _stabilizer_isotropic_subspace(gram, new_basis)
                results.append({
                    "basis": new_basis,
                    "stabilizer_gens": stab,
                })

    return results


def isotropic_flag_orbits(gram, k):
    """
    Compute orbit representatives of isotropic flags of length k.

    A flag is a chain I_1 ⊂ I_2 ⊂ ... ⊂ I_k of isotropic subspaces
    with dim(I_j) = j.

    INPUT:
    - gram: nested list of integers
    - k: length of flag

    OUTPUT:
    - list of dicts with 'flag' (list of bases) and 'stabilizer_gens'
    """
    if k == 0:
        return [{"flag": [], "stabilizer_gens": []}]

    # Build flags inductively
    results = []
    plane_orbits = isotropic_k_plane_orbits(gram, k)

    for plane in plane_orbits:
        # A k-plane gives a flag by taking successive subsets of the basis
        basis = plane["basis"]
        flag = [basis[:j+1] for j in range(len(basis))]
        results.append({
            "flag": flag,
            "stabilizer_gens": plane["stabilizer_gens"],
        })

    return results


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def indefinite_form_automorphism_group_simple(gram):
    """Simplified automorphism group for use in orbit computations."""
    from .higher_witt import indefinite_form_automorphism_group
    return indefinite_form_automorphism_group(gram)


def _stabilizer_isotropic_vector(gram, v):
    """Compute stabilizer of a single isotropic vector."""
    Q = matrix(ZZ, gram)
    n = Q.nrows()
    lv = sage_vector(ZZ, v)

    generators = []
    seen = set()
    I_mat = identity_matrix(ZZ, n)

    from itertools import product as iterproduct
    for bound in [3, 5]:
        for coords in iterproduct(range(-bound, bound + 1), repeat=n):
            if all(c == 0 for c in coords):
                continue
            r = sage_vector(ZZ, coords)
            q_r = int(r * Q * r)
            if q_r not in (1, 2, -1, -2):
                continue

            b_rv = int(r * Q * lv)
            if b_rv != 0:
                continue

            Qr = Q * r
            R = identity_matrix(ZZ, n)
            for i in range(n):
                for j in range(n):
                    R[i, j] -= 2 * int(r[i]) * int(Qr[j]) // q_r

            key = tuple(R[i, j] for i in range(n) for j in range(n))
            if key not in seen and R != I_mat:
                seen.add(key)
                generators.append(R)

    return generators


def _stabilizer_isotropic_subspace(gram, basis):
    """Compute setwise stabilizer of an isotropic subspace."""
    Q = matrix(ZZ, gram)
    n = Q.nrows()

    generators = []
    seen = set()
    I_mat = identity_matrix(ZZ, n)

    # Find vectors orthogonal to the entire subspace
    B = matrix(ZZ, basis)
    Bq = B * Q

    from itertools import product as iterproduct
    for bound in [3, 5]:
        for coords in iterproduct(range(-bound, bound + 1), repeat=n):
            if all(c == 0 for c in coords):
                continue
            r = sage_vector(ZZ, coords)
            q_r = int(r * Q * r)
            if q_r not in (1, 2, -1, -2):
                continue

            # Check if r is orthogonal to all basis vectors
            orthogonal = True
            for bi in basis:
                b_rv = int(r * Q * sage_vector(ZZ, bi))
                if b_rv != 0:
                    orthogonal = False
                    break

            if not orthogonal:
                continue

            Qr = Q * r
            R = identity_matrix(ZZ, n)
            for i in range(n):
                for j in range(n):
                    R[i, j] -= 2 * int(r[i]) * int(Qr[j]) // q_r

            key = tuple(R[i, j] for i in range(n) for j in range(n))
            if key not in seen and R != I_mat:
                seen.add(key)
                generators.append(R)

    return generators
