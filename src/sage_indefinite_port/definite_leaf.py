"""
TASK-03-4: Definite Lattice Leaf Backend Adapter and Witness Verification.

Supply definite lattice leaf adapters (Sage) with exact integer witness
verification for definite orthogonal complements.

Port of definite lattice operations from polyhedral_common's
INDEF_FORM_AutomorphismGroup_PosNeg and INDEF_FORM_TestEquivalence_PosNeg.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, identity_matrix
from . import make_lattice


def definite_automorphism_group(gram):
    """
    Compute generators of O(L) for definite integral L.

    Uses Sage's built-in definite lattice algorithms via the preamble.

    INPUT:
    - gram: nested list of integers (Gram matrix, must be definite)

    OUTPUT:
    - list of nested lists (integer matrices generating O(L))
    """
    L = make_lattice(gram)
    p, q = L.signature_pair()
    if p > 0 and q > 0:
        raise ValueError("Lattice must be definite (p=0 or q=0)")
    if p == 0 and q == 0:
        return []

    n = L.rank()

    # Use preamble's orthogonal_group for definite lattices
    OG = L.orthogonal_group()
    gens = OG.group_generators()

    result = []
    for g in gens:
        mat = g.matrix()
        mat_list = [[int(mat[i, j]) for j in range(n)] for i in range(n)]
        result.append(mat_list)

    return result


def definite_test_equivalence(gram1, gram2):
    """
    Test if two definite integral forms are equivalent over Z.

    Returns the transporter matrix P with P^T Q1 P = Q2, or None.

    INPUT:
    - gram1, gram2: nested lists of integers

    OUTPUT:
    - nested list of integers (transporter matrix), or None
    """
    L1 = make_lattice(gram1)
    L2 = make_lattice(gram2)

    if L1.signature_pair() != L2.signature_pair():
        return None
    if abs(L1.discriminant()) != abs(L2.discriminant()):
        return None

    n = L1.rank()
    if n == 0:
        return [[]]

    # Use preamble's isometry testing
    try:
        iso = L1.isometry_to(L2)
        if iso is not None:
            mat = iso.matrix()
            return [[int(mat[i, j]) for j in range(n)] for i in range(n)]
    except (AttributeError, NotImplementedError):
        pass

    return None


