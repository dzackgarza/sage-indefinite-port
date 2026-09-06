"""Structured isometry prefilters.

Replaces the ``size_t`` hash of ``INDEF_FORM_Invariant`` (polyhedral_common
``src_indefinite/IndefiniteFormFundamental.h``) by a structured record. A prefilter may
reject or bucket candidates; it never decides equivalence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from dzack_research.preamble.categories.lattices import FiniteRankLattices
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
    rank = Integer(lattice.rank())
    positive, negative = (Integer(index) for index in lattice.signature_pair())
    divisors = tuple(Integer(d) for d in lattice.correlation_morphism().cokernel().invariant_factors())
    return LatticePrefilter(
        rank=rank,
        signature=(positive, negative, rank - positive - negative),
        parity="even" if lattice.is_even() else "odd",
        discriminant=Integer(lattice.discriminant()),
        discriminant_elementary_divisors=divisors,
    )
