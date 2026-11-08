"""Duplicate-free orbit iterators over an enumerable set with decidable equivalence (#62)."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from itertools import count, product
from typing import TypeVar

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage_extensions.categories.lattices import Lattices

Point = TypeVar("Point")


def canonical_orbit_representatives(candidates: Iterable[Point], equivalent: Callable[[Point, Point], bool]) -> Iterator[Point]:
    r"""Yield the first member of each orbit, in the order of ``candidates``.

    #62 part B. A candidate is emitted iff it is equivalent to no earlier
    emitted one, so the output is duplicate-free and meets every orbit. The
    end of the output of an infinite walk carries no completion certificate.
    """
    emitted: list[Point] = []
    for candidate in candidates:
        if not any(equivalent(earlier, candidate) for earlier in emitted):
            emitted.append(candidate)
            yield candidate


def _hermite_matrices(rank: int, ambient: int, height: int) -> Iterator[list[list[int]]]:
    r"""Row Hermite normal forms of ``rank`` by ``ambient`` matrices of entries of maximum absolute value ``height``."""
    span = range(-height, height + 1)
    for entries in product(span, repeat=rank * ambient):
        rows = [list(entries[i * ambient : (i + 1) * ambient]) for i in range(rank)]
        candidate = matrix(SageZZ, rows)
        if candidate.rank() == rank and candidate.hermite_form() == candidate:
            yield rows


def primitive_isotropic_hermite_subobjects(lattice: Lattices.ParentMethods, rank: int, height_bound: int | None = None) -> Iterator:
    r"""Primitive isotropic rank-``rank`` subobjects by Hermite normal form, by height, without repetition."""
    labels = tuple(lattice.module_generating_set())
    seen: set[tuple[tuple[int, ...], ...]] = set()
    for height in count(1) if height_bound is None else range(1, height_bound + 1):
        for rows in _hermite_matrices(rank, len(labels), height):
            key = tuple(tuple(row) for row in rows)
            if key in seen:
                continue
            seen.add(key)
            if any(d != 1 for d in matrix(SageZZ, rows).elementary_divisors()[:rank]):
                continue
            vectors = tuple(
                lattice.linear_combination({label: lattice.base_ring()(entry) for label, entry in zip(labels, row, strict=True) if entry})
                for row in rows
            )
            if all(lattice.b(left, right) == 0 for left in vectors for right in vectors):
                yield lattice.primitive_isotropic_subobject(*vectors)


__all__ = ["canonical_orbit_representatives", "primitive_isotropic_hermite_subobjects"]
