"""Definite-lattice leaves delegated to the research preamble."""

from __future__ import annotations

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeIsometryMethods,
    LatticeIsometryMor,
)


def definite_orthogonal_group(
    lattice: Lattices.ParentMethods,
) -> LatticeIsometryMor:
    """Return the exact preamble orthogonal group of a definite lattice."""
    if not lattice.is_definite():
        raise ValueError(f"{lattice} is not definite")
    return lattice.O()


def definite_isometry(
    source: Lattices.ParentMethods,
    target: Lattices.ParentMethods,
) -> LatticeIsometryMethods:
    """Return and verify one exact preamble isometry between definite lattices."""
    if not source.is_definite() or not target.is_definite():
        raise ValueError("the definite leaf requires both lattices to be definite")
    witness = source.Isom(target).an_element()
    source_generators = tuple(source.module_generators())
    if any(source.b(left, right) != target.b(witness(left), witness(right)) for left in source_generators for right in source_generators):
        raise ArithmeticError("the preamble definite-isometry backend returned a map that does not preserve the form")
    return witness
