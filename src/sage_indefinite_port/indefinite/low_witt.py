"""Integral isometry witnesses at rational Witt index below two (#40).

The source's approximate model needs two rational hyperbolic planes, so it does
not cover these forms. Eichler's theorem and the spinor-genus decision settle
existence. Once existence is certified, a fair enumeration of GL_r(Z) by growing
coordinate boxes terminates at a witness: the search has no bound to exhaust.
"""

from __future__ import annotations

from collections.abc import Iterator
from itertools import count, product

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage_extensions.categories.lattices import Lattices
from sage_extensions.tensors.tensor import _engine_component_matrix


def _gram(lattice: Lattices.ParentMethods):
    return _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)


def _search_box(source_gram, target_gram, radius: int):
    r"""Yield the rows of each ``T`` in the box with ``T G_target T^t = G_source`` and ``det T = +-1``."""
    rank = source_gram.nrows()
    target = [[int(entry) for entry in row] for row in target_gram.rows()]

    def pairing(left, right):
        return sum(left[i] * target[i][j] * right[j] for i in range(rank) for j in range(rank))

    by_norm: dict[int, list[tuple[int, ...]]] = {}
    for candidate in product(range(-radius, radius + 1), repeat=rank):
        by_norm.setdefault(pairing(candidate, candidate), []).append(candidate)

    def extend(rows):
        depth = len(rows)
        if depth == rank:
            if abs(matrix(SageZZ, rows).det()) == 1:
                yield rows
            return
        for candidate in by_norm.get(int(source_gram[depth, depth]), ()):
            if all(pairing(rows[i], candidate) == int(source_gram[i, depth]) for i in range(depth)):
                yield from extend([*rows, candidate])

    yield from extend([])


def _isometry_from_rows(source, target, rows):
    labels = tuple(target.module_generating_set())
    images = dict(zip(source.module_generating_set(), (
        target.linear_combination({
            label: target.base_ring()(entry) for label, entry in zip(labels, row, strict=True) if entry
        })
        for row in rows
    ), strict=True))
    return source.Isom(target)(images)


def fair_orthogonal_group_elements(lattice: Lattices.ParentMethods) -> Iterator[LatticeIsometryMethods]:
    r"""Enumerate ``O(L)`` without repetition by growing boxes: every element appears at a finite step."""
    gram = _gram(lattice)
    seen: set[tuple[tuple[int, ...], ...]] = set()
    for radius in count(1):
        for rows in _search_box(gram, gram, radius):
            key = tuple(tuple(row) for row in rows)
            if key not in seen:
                seen.add(key)
                yield _isometry_from_rows(lattice, lattice, rows)


def low_witt_isometry_witness(
    source: Lattices.ParentMethods, target: Lattices.ParentMethods
) -> LatticeIsometryMethods:
    r"""An isometry ``source -> target``. The caller has certified that one exists."""
    source_gram, target_gram = _gram(source), _gram(target)
    for radius in count(1):
        rows = next(_search_box(source_gram, target_gram, radius), None)
        if rows is not None:
            return _isometry_from_rows(source, target, rows)
    raise AssertionError("count() does not terminate")


__all__ = ["fair_orthogonal_group_elements", "low_witt_isometry_witness"]
