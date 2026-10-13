"""Reduced lattice presentations retaining exact isometries."""

from __future__ import annotations

from dataclasses import dataclass

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices

from sage_indefinite_port.backends.canonization import presentation_bucket_key


@dataclass(frozen=True)
class ReducedLatticePresentation:
    """A reduced presentation together with its exact source isometry."""

    source: Lattices.ParentMethods
    reduced: Lattices.ParentMethods
    isometry: LatticeIsometryMethods

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> ReducedLatticePresentation:
        reduction = lattice.lll_reduction()
        return cls(reduction.isometry.domain(), reduction.isometry.codomain(), reduction.isometry)


__all__ = ["ReducedLatticePresentation", "presentation_bucket_key"]
