"""Eichler-cover primitives delegated to the preamble lattice owner."""

from __future__ import annotations

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices


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
    return tuple(
        divisor
        for divisor in range(1, int(value**0.5) + 1)
        if value % (divisor * divisor) == 0
    )


__all__ = ["eichler_transvection", "square_divisors"]
