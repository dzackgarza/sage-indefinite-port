"""Structured isometry prefilters.

Replaces the ``size_t`` hash of ``INDEF_FORM_Invariant`` (polyhedral_common
``src_indefinite/IndefiniteFormFundamental.h``) by a structured record. A prefilter may
reject or bucket candidates; it never decides equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass
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
    discriminant_elementary_divisors: tuple[Integer, ...]


def lattice_prefilter(lattice: FiniteRankLattices.ParentMethods) -> LatticePrefilter:
    """Rank, real signature ``(n_+, n_-, n_0)``, parity, signed discriminant, and the
    invariants of the cokernel of the correlation ``L -> Hom(L, Z)``, the finite abelian
    group ``L^vee / L``."""
    base_ring = lattice.base_ring()
    rank = Integer(lattice.rank())
    positive, negative = (Integer(index) for index in lattice.signature_pair())
    divisors = tuple(Integer(_engine_element(base_ring, d)) for d in lattice.correlation_morphism().cokernel().invariant_factors())
    return LatticePrefilter(
        rank=rank,
        signature=(positive, negative, rank - positive - negative),
        parity="even" if lattice.is_even() else "odd",
        discriminant=Integer(_engine_element(base_ring, lattice.discriminant())),
        discriminant_elementary_divisors=divisors,
    )


@dataclass(frozen=True)
class AttackProfile:
    lattice: Lattices.ParentMethods
    signed_view: Lattices.ParentMethods
    sign: int
    positive_index: int
    negative_index: int

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> "AttackProfile":
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
    def from_vector(cls, vector: Lattices.ElementMethods) -> "VectorPrefilter":
        lattice = vector.parent()
        discriminant_key = None
        if vector.is_primitive() and vector.div() != lattice.base_ring().zero():
            discriminant_key = repr(vector.divided_discriminant_class())
        reduction_prefilter = None
        if vector.is_isotropic() and vector.is_primitive():
            reduction_prefilter = lattice_prefilter(vector.isotropic_reduction())
        return cls(int(vector.q()), vector_content(vector), int(vector.div()), discriminant_key, reduction_prefilter)
