"""Structured isometry prefilters.

Replaces the ``size_t`` hash of ``INDEF_FORM_Invariant`` (polyhedral_common
``src_indefinite/IndefiniteFormFundamental.h``) by a structured record. A prefilter may
reject or bucket candidates; it never decides equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Literal

from dzack_research.preamble.categories.lattices import FiniteRankLattices, Lattices
from dzack_research.preamble.categories.rings.ring_foundation import _engine_element
from sage.rings.integer import Integer

Parity = Literal["even", "odd"]


@dataclass(frozen=True)
class LatticePrefilter:
    """Isometry invariants that are cheap to compute and to compare."""

    rank: Integer
    signature: tuple[Integer, Integer, Integer]
    parity: Parity
    discriminant: Integer


def lattice_prefilter(lattice: FiniteRankLattices.ParentMethods) -> LatticePrefilter:
    """Return the structured form of the pinned ``INDEF_FORM_Invariant`` key.

    The reference implementation hashes rank, real signature, parity, and the
    determinant of the nondegenerate quotient.  This prefilter keeps those
    invariants as a record instead of hashing them; it deliberately does not
    construct the discriminant group or correlation cokernel.
    """
    base_ring = lattice.base_ring()
    rank = Integer(lattice.rank())
    positive, negative = (Integer(index) for index in lattice.signature_pair())
    return LatticePrefilter(
        rank=rank,
        signature=(positive, negative, rank - positive - negative),
        parity="even" if lattice.is_even() else "odd",
        discriminant=Integer(_engine_element(base_ring, lattice.discriminant())),
    )


@dataclass(frozen=True)
class AttackProfile:
    lattice: Lattices.ParentMethods
    signed_view: Lattices.ParentMethods
    sign: int
    positive_index: int
    negative_index: int

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> AttackProfile:
        positive, negative = lattice.signature_pair()
        positive_index = int(positive)
        negative_index = int(negative)
        if positive_index <= negative_index:
            return cls(lattice, lattice, 1, positive_index, negative_index)
        signed = lattice.twist(-lattice.base_ring().one())
        return cls(lattice, signed, -1, negative_index, positive_index)


def vector_content(vector: Lattices.ElementMethods) -> int:
    coordinates = vector.to_vector()
    content = 0
    for label in vector.parent().module_generating_set():
        content = gcd(content, abs(int(coordinates(label))))
    return content


@dataclass(frozen=True)
class VectorPrefilter:
    norm: int
    content: int
    divisor: int
    discriminant_class_orbit_key: str | None
    orthogonal_reduction_prefilter: LatticePrefilter | None

    @classmethod
    def from_vector(
        cls,
        vector: Lattices.ElementMethods,
        *,
        include_orthogonal_reduction: bool = True,
    ) -> VectorPrefilter:
        norm = int(vector.q())
        content = vector_content(vector)
        divisor = int(vector.div())
        discriminant_key = None
        match divisor:
            case 0:
                pass
            case _:
                positive_divisor = abs(divisor)
                order = positive_divisor // gcd(positive_divisor, content)
                discriminant_key = (
                    "trivial"
                    if order == 1
                    else repr(
                        (
                            order,
                            Fraction(norm, positive_divisor * positive_divisor),
                        )
                    )
                )
        reduction_prefilter = None
        if include_orthogonal_reduction and vector.is_isotropic() and vector.is_primitive():
            reduction_prefilter = lattice_prefilter(vector.isotropic_reduction())
        return cls(
            norm,
            content,
            divisor,
            discriminant_key,
            reduction_prefilter,
        )
