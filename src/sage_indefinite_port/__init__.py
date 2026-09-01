"""
sage_indefinite_port: Pure SageMath indefinite lattice algorithms.

All lattice-level operations go through the preamble API.
Isometries are returned as preamble morphisms, not raw matrices.
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


def matrix_to_morphism(L, M):
    """Wrap a Sage integer matrix M as a preamble isometry morphism L -> L.

    Returns a morphism element with .matrix(), .domain(), .codomain().
    """
    if isinstance(M, list):
        M = matrix(ZZ, M)
    n = L.rank()
    H = L.Hom(L)
    return H(M)


def morphism_to_matrix(mor):
    """Extract the Sage integer matrix from a preamble morphism."""
    return mor.matrix()


def verify_isometry(L, mor):
    """Check that a morphism is an isometry: b(mor(ei), mor(ej)) == b(ei, ej)."""
    n = L.rank()
    for i in range(n):
        ei = L.module_generator(i)
        for j in range(n):
            ej = L.module_generator(j)
            if L.b(mor(ei), mor(ej)) != L.b(ei, ej):
                return False
    return True


def verify_matrix_isometry(gram, M):
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


def unique_morphisms(mors, L):
    """Remove duplicate and identity morphisms."""
    seen = set()
    result = []
    n = L.rank()
    id_mor = L.identity_morphism()
    for m in mors:
        M = m.matrix()
        key = tuple(M[i, j] for i in range(n) for j in range(n))
        if key not in seen and m != id_mor:
            seen.add(key)
            result.append(m)
    return result


def matrix_list_to_morphism_list(L, mat_list):
    """Convert a list of integer matrices to preamble morphisms."""
    return [matrix_to_morphism(L, M) for M in mat_list]


# Public API
from .higher_witt import indefinite_form_automorphism_group, indefinite_form_test_equivalence
from .lorentzian_perfect import get_attack_scheme, is_lorentzian, lorentzian_generators_autom
from .pre_sieves import LatticeInvariant, VectorInvariant, lattice_pre_sieve, vector_pre_sieve
from .equivariant import EquivariantLattice, with_isometry, centralizer_group
from .semantic_wrapper import indefinite_lattice, orthogonal_generators, test_isometry, isotropic_orbits
