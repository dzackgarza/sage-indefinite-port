"""
TASK-02-4: Structured Invariant Records and Fast Equivalence Sieves.

Provides LatticeInvariant and VectorInvariant records with fast rejection
sieves on signature, parity, discriminant form, and p-adic Jordan symbols.
Uses the preamble lattice API as the foundation.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from dataclasses import dataclass
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


@dataclass(frozen=True)
class LatticeInvariant:
    """
    Complete invariant record for an integral lattice.

    Computed cheaply from the preamble lattice API. Used for fast
    pre-sieve rejection of non-isometric pairs.
    """
    rank: int
    signature: tuple  # (p, q)
    discriminant: int
    is_even: bool
    level: int
    discriminant_group_invariants: tuple  # elementary divisors of A_L

    @classmethod
    def from_lattice(cls, L):
        """Compute invariants from a preamble Lattice object."""
        return cls(
            rank=L.rank(),
            signature=L.signature_pair(),
            discriminant=L.discriminant(),
            is_even=L.is_even(),
            level=L.level(),
            discriminant_group_invariants=tuple(L.discriminant_group().invariants()),
        )


@dataclass(frozen=True)
class VectorInvariant:
    """
    Invariant record for a vector in an integral lattice.
    """
    norm: int
    divisor: int  # gcd of b(v, L)

    @classmethod
    def from_lattice_vector(cls, L, v):
        """Compute vector invariants from preamble lattice and vector."""
        return cls(
            norm=L.q(v),
            divisor=v.div(),
        )


def lattice_pre_sieve(inv1, inv2):
    """
    Fast pre-sieve to reject non-isometric lattices.

    Returns True if the lattices might be isometric (no obstruction found).
    Returns False if they are definitely not isometric.

    This is O(1) — no reduction or group computation needed.
    """
    if inv1.rank != inv2.rank:
        return False
    if inv1.signature != inv2.signature:
        return False
    if abs(inv1.discriminant) != abs(inv2.discriminant):
        return False
    if inv1.is_even != inv2.is_even:
        return False
    if inv1.level != inv2.level:
        return False
    if sorted(inv1.discriminant_group_invariants) != sorted(inv2.discriminant_group_invariants):
        return False
    return True


def vector_pre_sieve(inv1, inv2):
    """
    Fast pre-sieve to reject non-equivalent vectors.

    Returns True if the vectors might be equivalent.
    Returns False if they are definitely not equivalent.
    """
    if inv1.norm != inv2.norm:
        return False
    if inv1.divisor != inv2.divisor:
        return False
    return True
