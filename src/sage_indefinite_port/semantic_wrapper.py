"""
TASK-06-2: Coordinate-Free Sage Category Integration.

Package all primitives into coordinate-free Sage category interfaces
(Lattices, OrthogonalGroups, Subobjects, Loci) with zero coordinate
or temporary file leaks.

Designed for future upstream contribution to dzack_research.preamble.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.tensors import tensor
from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Public API: coordinate-free wrappers around low-level algorithms
# ---------------------------------------------------------------------------

@dataclass
class IndefiniteLattice:
    """
    Coordinate-free wrapper for an indefinite integral lattice.

    Wraps a Gram matrix with coordinate-free semantic methods.
    """
    _gram: list  # Gram matrix as nested list

    @classmethod
    def from_gram(cls, gram):
        """Create from a Gram matrix (nested list of integers)."""
        return cls(_gram=[list(row) for row in gram])

    @classmethod
    def from_preamble_lattice(cls, L):
        """Create from a preamble Lattice object (without editing preamble)."""
        n = L.rank()
        gram = [[int(L.b(L.module_generator(i), L.module_generator(j)))
                  for j in range(n)] for i in range(n)]
        return cls(_gram=gram)

    @property
    def rank(self):
        return len(self._gram)

    @property
    def gram_matrix(self):
        return self._gram

    def signature(self):
        """Return (p, q) signature via preamble L.signature_pair()."""
        L = self._preamble_lattice()
        return L.signature_pair()

    def _preamble_lattice(self):
        """Build preamble lattice from stored gram."""
        C = Lattices(ZZ)
        n = self.rank
        rows = tuple(tuple(int(x) for x in row) for row in self._gram)
        g = tensor(ZZ, (), (n, n), rows)
        return C(g)

    def determinant(self):
        """Return |det(Q)|."""
        from sage.all import matrix as sage_matrix
        Q = sage_matrix(ZZ, self._gram)
        return abs(Q.det())

    def is_even(self):
        """Check if all diagonal entries are even."""
        return all(self._gram[i][i] % 2 == 0 for i in range(self.rank))

    def orthogonal_group(self):
        """Compute generators of O(L)."""
        from .higher_witt import indefinite_form_automorphism_group
        gens = indefinite_form_automorphism_group(self._gram)
        return OrthogonalGroup(self, gens)

    def isometry_to(self, other):
        """Test if self is isometric to other. Returns transporter or None."""
        from .higher_witt import indefinite_form_test_equivalence
        result = indefinite_form_test_equivalence(self._gram, other._gram)
        return result

    def invariants(self):
        """Compute fast invariants for pre-sieving."""
        from .pre_sieves import LatticeInvariant
        from dzack_research.preamble.categories.lattices import Lattices
        from dzack_research.preamble.tensors import tensor

        C = Lattices(ZZ)
        n = self.rank
        rows = tuple(tuple(int(x) for x in row) for row in self._gram)
        g = tensor(ZZ, (), (n, n), rows)
        L = C(g)
        return LatticeInvariant.from_lattice(L)


@dataclass
class OrthogonalGroup:
    """
    Wraps generators of an orthogonal group with semantic operations.
    """
    lattice: IndefiniteLattice
    _generators: list  # list of nested lists (matrix generators)

    @property
    def generators(self):
        return self._generators

    def order(self):
        """Attempt to compute group order (finite groups only)."""
        from sage.all import MatrixGroup, ZZ as SageZZ
        gens = [matrix(SageZZ, g) for g in self._generators]
        try:
            G = MatrixGroup(gens)
            return G.order()
        except Exception:
            return None

    def centralizer(self, isometry):
        """Compute centralizer of an isometry in this group."""
        from .equivariant import centralizer_group
        return centralizer_group(self.lattice._gram, isometry)

    def orbit(self, v):
        """Compute orbit of vector v under this group."""
        # Placeholder for orbit enumeration
        return [v]


@dataclass
class IsotropicSubspace:
    """
    Wraps an isotropic subspace with reduction data.
    """
    lattice: IndefiniteLattice
    basis: list  # list of vectors (nested lists)

    @property
    def dimension(self):
        return len(self.basis)

    def quotient_lattice(self):
        """Compute K_I = I^perp / I."""
        from .isotropic_orbits import isotropic_reduction
        return isotropic_reduction(self.lattice._gram, self.basis)

    def stabilizer(self):
        """Compute parabolic stabilizer P_I."""
        from .isotropic_orbits import parabolic_stabilizer
        gens = parabolic_stabilizer(self.lattice._gram, self.basis)
        return OrthogonalGroup(self.lattice, gens)


# ---------------------------------------------------------------------------
# Category-style methods (aligned with preamble interfaces)
# ---------------------------------------------------------------------------

def indefinite_lattice(gram):
    """
    Create an IndefiniteLattice with full computational backend.

    This is the main entry point, analogous to `Lattices(ZZ)(Q)` in preamble
    but for indefinite lattices.
    """
    return IndefiniteLattice.from_gram(gram)


def orthogonal_generators(gram):
    """
    Compute generators of O(L) for indefinite L.

    Main API entry point replacing py_polyhedral.binaries shell-outs.
    """
    L = IndefiniteLattice.from_gram(gram)
    OG = L.orthogonal_group()
    return OG.generators


def test_isometry(gram1, gram2):
    """
    Test if two indefinite integral forms are equivalent.

    Returns transporter matrix or None.
    """
    L1 = IndefiniteLattice.from_gram(gram1)
    L2 = IndefiniteLattice.from_gram(gram2)
    return L1.isometry_to(L2)


def isotropic_orbits(gram, k):
    """
    Compute orbit representatives of k-dimensional isotropic sublattices.
    """
    from .isotropic_orbits import isotropic_k_plane_orbits
    return isotropic_k_plane_orbits(gram, k)
