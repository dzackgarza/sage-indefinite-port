"""Eichler-cover primitives delegated to the preamble lattice owner."""

from __future__ import annotations

from dataclasses import dataclass

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices


class InfiniteLocusError(ValueError):
    """Raised when a requested orbit cover is not finite."""


@dataclass(frozen=True)
class OrbitCover:
    """A finite covering family, not necessarily a full-orbit decomposition."""

    representatives: tuple[Lattices.ElementMethods, ...]

    def __iter__(self):
        return iter(self.representatives)

    def __len__(self) -> int:
        return len(self.representatives)


@dataclass(frozen=True)
class EichlerOrbitCover:
    """Immutable covering model for the preamble's represented ``2U + K`` lattice."""

    model: object

    def lattice(self):
        return self.model.lattice()

    def subgroup(self):
        return self.model.approximate_generating_family()

    def covering_representatives(self, norm, *, primitive: bool) -> OrbitCover:
        ring = self.lattice().base_ring()
        owned_norm = norm if getattr(norm, "parent", lambda: None)() is ring else ring(int(norm))
        if primitive:
            family = self.model.covering_vector_representatives(owned_norm)
            return OrbitCover(tuple(family[label] for label in family.index_set()))
        if int(owned_norm) == 0:
            raise InfiniteLocusError("nonprimitive isotropic vectors form an infinite locus")
        lattice = self.lattice()
        representatives = []
        for divisor in square_divisors(owned_norm):
            primitive_norm = ring(int(owned_norm) // (divisor * divisor))
            family = self.model.covering_vector_representatives(primitive_norm)
            scalar = lattice.base_ring()(divisor)
            representatives.extend(lattice.scalar_multiple(scalar, family[label]) for label in family.index_set())
        return OrbitCover(tuple(representatives))

    def one_representative(self, norm, *, primitive: bool):
        cover = self.covering_representatives(norm, primitive=primitive)
        if not cover.representatives:
            raise ValueError(f"no covering representative exists for norm {norm}")
        return cover.representatives[0]


def eichler_transvection(
    isotropic: Lattices.ElementMethods,
    orthogonal: Lattices.ElementMethods,
) -> LatticeIsometryMethods:
    r"""Return the exact Eichler transvection ``E_(f,x)`` in the ambient lattice."""
    lattice = isotropic.parent()
    if orthogonal.parent() is not lattice:
        raise ValueError("an Eichler transvection needs two vectors in one lattice")
    if not lattice.is_even():
        raise ValueError("an Eichler transvection in this port requires an even lattice")
    if not isotropic.is_isotropic():
        raise ValueError("the first Eichler-transvection vector must be isotropic")
    if lattice.b(isotropic, orthogonal) != lattice.base_ring().zero():
        raise ValueError("the second Eichler-transvection vector must lie in f^perp")
    return lattice.eichler_transvection(isotropic, orthogonal)


def square_divisors(integer) -> tuple[int, ...]:
    r"""Return positive ``c`` such that ``c^2`` divides the nonzero integer ``integer``."""
    value = abs(int(integer))
    if value == 0:
        raise ValueError("zero has infinitely many square divisors")
    return tuple(divisor for divisor in range(1, int(value**0.5) + 1) if value % (divisor * divisor) == 0)


__all__ = [
    "EichlerOrbitCover",
    "InfiniteLocusError",
    "OrbitCover",
    "eichler_transvection",
    "square_divisors",
]
