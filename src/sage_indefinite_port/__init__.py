"""
sage_indefinite_port: Pure SageMath indefinite lattice algorithms.

All lattice-level operations go through the preamble API:
  L.b(), L.q(), L.module_generator(), L.signature_pair(), L.discriminant()
"""

__version__ = "0.1.0"

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, matrix, identity_matrix, vector as sage_vector
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.tensors import tensor


def make_lattice(gram):
    """Build a preamble Lattice from a nested-list Gram matrix."""
    n = len(gram)
    rows = tuple(tuple(int(x) for x in row) for row in gram)
    g = tensor(ZZ, (), (n, n), rows)
    return Lattices(ZZ)(g)


def from_sage_matrix(M):
    """Convert Sage matrix to nested list of ints."""
    return [[int(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


def verify_generator(gram, M):
    """Check M^T Q M == Q using preamble b()."""
    L = make_lattice(gram)
    n = L.rank()
    if isinstance(M, list):
        M = matrix(ZZ, M)
    for i in range(n):
        ei = L.module_generator(i)
        for j in range(n):
            ej = L.module_generator(j)
            Mi = L(sage_vector(ZZ, [int(M[i, k]) for k in range(n)]))
            Mj = L(sage_vector(ZZ, [int(M[j, k]) for k in range(n)]))
            if L.b(Mi, Mj) != L.b(ei, ej):
                return False
    return True


def unique_generators(gens, n):
    """Remove duplicate and trivial (identity) generators."""
    seen = set()
    result = []
    I = identity_matrix(ZZ, n)
    for g in gens:
        M = matrix(ZZ, g) if isinstance(g, list) else g
        key = tuple(M[i, j] for i in range(n) for j in range(n))
        if key not in seen and M != I:
            seen.add(key)
            result.append(M)
    return result


# Public API
from .higher_witt import indefinite_form_automorphism_group, indefinite_form_test_equivalence
from .lorentzian_perfect import get_attack_scheme, is_lorentzian, lorentzian_generators_autom
from .pre_sieves import LatticeInvariant, VectorInvariant, lattice_pre_sieve, vector_pre_sieve
from .equivariant import EquivariantLattice, with_isometry, centralizer_group
from .semantic_wrapper import indefinite_lattice, orthogonal_generators, test_isometry, isotropic_orbits
