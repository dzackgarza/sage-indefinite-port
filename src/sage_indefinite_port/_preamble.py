"""Boundary between raw Gram data and the research preamble's lattice category.

``lattice_from_gram`` is the single entry point where integer Gram rows (fixtures,
published tables) become preamble ``Lattice`` objects. Every other module consumes
preamble objects and their types directly.
"""

from __future__ import annotations

from collections.abc import Sequence

from dzack_research.preamble.categories._lattice import Lattice
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.tensors import tensor
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


def lattice_from_gram(rows: Sequence[Sequence[int | Integer]]) -> Lattice:
    """The preamble lattice on ``Z^n`` whose Gram matrix has the given integer rows."""
    size = len(rows)
    entries = tuple(tuple(Integer(entry) for entry in row) for row in rows)
    return Lattices(ZZ)(tensor(ZZ, (), (size, size), entries))
