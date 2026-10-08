"""Isometry prefilters and the sign normalization of the Eichler recursion.

Replaces the ``size_t`` hash of ``INDEF_FORM_Invariant`` (polyhedral_common
``src_indefinite/IndefiniteFormFundamental.h``) by structured records of preamble invariants.
A prefilter may reject or bucket candidates; it never decides equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from dzack_research.preamble.categories.lattices import Lattices


@dataclass(frozen=True)
class LatticePrefilter:
    """Rank, signature, parity and discriminant: equal for isometric lattices."""

    rank: int
    signature: tuple[int, int]
    even: bool
    discriminant: int

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> LatticePrefilter:
        positive, negative = lattice.signature_pair()
        return cls(int(lattice.rank()), (int(positive), int(negative)), bool(lattice.is_even()), int(lattice.discriminant()))


@dataclass(frozen=True)
class AttackProfile:
    """L or L(-1), whichever has positive index min(n_+, n_-)."""

    lattice: Lattices.ParentMethods
    signed_view: Lattices.ParentMethods
    sign: int
    positive_index: int
    negative_index: int

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> AttackProfile:
        positive, negative = (int(index) for index in lattice.signature_pair())
        if positive <= negative:
            return cls(lattice, lattice, 1, positive, negative)
        return cls(lattice, lattice.twist(-lattice.base_ring().one()), -1, negative, positive)


@dataclass(frozen=True)
class VectorPrefilter:
    r"""O(L)-invariants of a vector v: q(v), content, div(v), and the order and norm of its class v/div(v) in A_L.

    For a primitive isotropic v it also holds the prefilter of v^perp/v.
    """

    norm: int
    content: int
    divisibility: int
    class_order: int
    class_norm: Fraction
    isotropic_reduction: LatticePrefilter | None

    @classmethod
    def from_vector(cls, vector: Lattices.ElementMethods) -> VectorPrefilter:
        norm = int(vector.q())
        content = int(vector.content())
        divisibility = int(vector.div())
        discriminant_class = vector.divided_discriminant_class()
        reduction = LatticePrefilter.from_lattice(vector.isotropic_reduction()) if norm == 0 and content == 1 else None
        return cls(norm, content, divisibility, int(discriminant_class.order()), Fraction(norm, divisibility * divisibility), reduction)
