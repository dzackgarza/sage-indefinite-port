"""Signed attack profile and structured prefilters.

Ports ``AttackScheme`` and ``INDEF_FORM_GetAttackScheme`` from polyhedral_common
``src_indefinite/IndefiniteFormFundamental.h`` as a signed view of the lattice, and replaces
the ``size_t`` hashes of ``INDEF_FORM_Invariant`` and ``INDEF_FORM_InvariantVector`` by
structured records. A prefilter may reject or bucket candidates; it never decides
equivalence.

The upstream ``divisor``/``index`` scalars of a vector are named here by their meaning:
``content(v) = max{d : v in dL}`` and ``div(v) = gcd b(v, L)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from dzack_research.preamble.categories._lattice import Lattice
from sage.arith.misc import gcd
from sage.matrix.constructor import matrix
from sage.rings.integer import Integer

type LatticeElement = Lattice.Element

Parity = Literal["even", "odd"]


class ZeroVectorError(ValueError):
    """The zero vector has no content."""


@dataclass(frozen=True)
class AttackProfile:
    """The lattice with the sign-normalized view whose positive index is minimal.

    ``signed_view`` is ``lattice`` when ``n_+ <= n_-`` and ``lattice.twist(-1)`` otherwise;
    ``positive_index`` is ``min(n_+, n_-)``, the bound on the Witt index that drives dispatch.
    """

    lattice: Lattice
    signed_view: Lattice
    sign: Integer
    positive_index: Integer
    negative_index: Integer


def attack_profile(lattice: Lattice) -> AttackProfile:
    """Normalize the sign of the form so that the positive index is ``min(n_+, n_-)``."""
    positive, negative = (Integer(index) for index in lattice.signature_pair())
    if positive <= negative:
        return AttackProfile(lattice, lattice, Integer(1), positive, negative)
    return AttackProfile(lattice, lattice.twist(-1), Integer(-1), negative, positive)


@dataclass(frozen=True)
class LatticePrefilter:
    """Isometry invariants that are cheap to compute and to compare."""

    rank: Integer
    signature: tuple[Integer, Integer, Integer]
    parity: Parity
    discriminant: Integer
    discriminant_elementary_divisors: tuple[Integer, ...]


def lattice_prefilter(lattice: Lattice) -> LatticePrefilter:
    """Rank, real signature ``(n_+, n_-, n_0)``, parity, signed discriminant, and the
    elementary divisors of the Gram matrix other than 1 (the invariants of ``L^vee/L``)."""
    rank = Integer(lattice.rank())
    positive, negative = (Integer(index) for index in lattice.signature_pair())
    gram = matrix(lattice.gram_tensor().components())
    divisors = tuple(d for d in gram.elementary_divisors() if d != 1)
    return LatticePrefilter(
        rank=rank,
        signature=(positive, negative, rank - positive - negative),
        parity="even" if lattice.is_even() else "odd",
        discriminant=Integer(lattice.discriminant()),
        discriminant_elementary_divisors=divisors,
    )


def content(vector: LatticeElement) -> Integer:
    """The largest ``d`` with ``v in dL``: the gcd of the coordinates of ``v``."""
    value = Integer(gcd(list(vector.to_tuple())))
    if value == 0:
        raise ZeroVectorError("the zero vector has no content")
    return value


@dataclass(frozen=True)
class VectorPrefilter:
    """Orbit invariants of a nonzero lattice vector."""

    norm: Integer
    content: Integer
    divisor: Integer


def vector_prefilter(lattice: Lattice, vector: LatticeElement) -> VectorPrefilter:
    """Norm ``q(v) = b(v, v)``, content, and divisor ``gcd b(v, L)`` of ``v``."""
    return VectorPrefilter(
        norm=Integer(lattice.q(vector)),
        content=content(vector),
        divisor=Integer(vector.div()),
    )
