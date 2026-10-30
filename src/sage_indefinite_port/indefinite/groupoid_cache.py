"""Isometry-groupoid memoization from IndefiniteForms.g:1173-1299.

The pinned source retains exact equivalence witnesses, proven non-equivalences
between their connected components, and orthogonal groups transported along
isometries. Equal Gram keys give a coordinate-identification isometry;
different keys require witnesses, never a canonical-form assumption.
"""

from __future__ import annotations

from collections import deque

from dzack_research.preamble.all import ZZ, Lattices
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ

type GramKey = tuple[tuple[int, ...], ...]
type Lattice = Lattices.ParentMethods
type Isometry = LatticeIsometryMethods


def _gram_in_label_order(lattice: Lattice):
    """Gram matrix of the declared coordinate generators, in their given order."""
    labels = tuple(lattice.module_generating_set())
    basis = tuple(lattice.module_generator(label) for label in labels)
    if lattice.base_ring() == ZZ:
        entries = [[int(lattice.b(left, right)) for right in basis] for left in basis]
        return matrix(SageZZ, entries)
    entries = [
        [
            SageQQ(int(value.numerator())) / int(value.denominator())
            for value in (lattice.b(left, right) for right in basis)
        ]
        for left in basis
    ]
    return matrix(SageQQ, entries)


def coordinate_identification(source: Lattice, target: Lattice) -> Isometry:
    """The explicit basis-to-basis isometry of identical Gram presentations."""
    if source is target:
        return source.O().identity()
    if source.base_ring() != target.base_ring():
        raise ValueError("coordinate identification requires a common base ring")
    source_labels = tuple(source.module_generating_set())
    target_labels = tuple(target.module_generating_set())
    if len(source_labels) != len(target_labels):
        raise ValueError("coordinate identification requires equal coordinate ranks")
    source_gram = _gram_in_label_order(source)
    target_gram = _gram_in_label_order(target)
    if source_gram != target_gram:
        raise ValueError("coordinate identification requires identical Gram matrices in declared label order")
    images = []
    for generator in source.module_generators():
        coordinates = generator.to_vector()
        images.append(
            target.linear_combination(
                {target_label: coordinates(source_label) for source_label, target_label in zip(source_labels, target_labels, strict=True)}
            )
        )
    return source.Isom(target)(tuple(images))


class IsometryGroupoidCache:
    """Source-backed equivalence graph and exact Aut transport."""

    def __init__(self) -> None:
        self._representatives: dict[GramKey, Lattice] = {}
        self._adjacency: dict[GramKey, list[tuple[GramKey, Isometry]]] = {}
        self._nonisometric: set[frozenset[GramKey]] = set()
        self._groups: dict[GramKey, tuple[Isometry, ...]] = {}

    def _key(self, lattice: Lattice) -> GramKey:
        if lattice.base_ring() != ZZ:
            raise TypeError("the isometry groupoid cache requires integral ZZ-lattices")
        gram = _gram_in_label_order(lattice)
        key = tuple(tuple(int(gram[row, column]) for column in range(gram.ncols())) for row in range(gram.nrows()))
        self._representatives.setdefault(key, lattice)
        return key

    def _component(self, start: GramKey) -> set[GramKey]:
        reached = {start}
        pending = [start]
        while pending:
            for neighbour, _witness in self._adjacency.get(pending.pop(), ()):
                if neighbour not in reached:
                    reached.add(neighbour)
                    pending.append(neighbour)
        return reached

    def _path(self, start: GramKey, target: GramKey) -> Isometry | None:
        pending = deque([(start, self._representatives[start].O().identity())])
        reached = {start}
        while pending:
            current, witness = pending.popleft()
            if current == target:
                return witness
            for neighbour, edge in self._adjacency.get(current, ()):
                if neighbour not in reached:
                    reached.add(neighbour)
                    pending.append((neighbour, edge * witness))
        return None

    def cached_isometry(self, source: Lattice, target: Lattice) -> tuple[bool, Isometry | None]:
        """Return ``(known, witness)``; a known negative has witness ``None``.

        As in ``INDEF_FORM_TestEquivalence_Memoize`` (source lines 1239-1273),
        positive edges compose along paths and negative edges separate entire
        connected components, not merely their originally queried endpoints.
        """
        start, end = self._key(source), self._key(target)
        path = self._path(start, end)
        if path is not None:
            return (
                True,
                coordinate_identification(self._representatives[end], target) * path * coordinate_identification(source, self._representatives[start]),
            )
        left, right = self._component(start), self._component(end)
        if any((negative & left) and (negative & right) for negative in self._nonisometric):
            return True, None
        return False, None

    def record_isometry(self, source: Lattice, target: Lattice, witness: Isometry | None) -> None:
        """Record an exact kernel verdict, retaining its isometry certificate.

        Equal Gram matrices need no extra edge: the coordinate identification
        is already a certificate. Non-equivalence is retained only after an
        exact full isometry decision, never after an invariant prefilter.
        """
        start, end = self._key(source), self._key(target)
        if witness is None:
            if self._path(start, end) is not None:
                raise ArithmeticError("an established isometry contradicts non-equivalence")
            self._nonisometric.add(frozenset((start, end)))
            return
        if witness.domain() is not source or witness.codomain() is not target:
            raise ValueError("the stored isometry must have the supplied endpoints")
        if start == end:
            return
        if any((negative & self._component(start)) and (negative & self._component(end)) for negative in self._nonisometric):
            raise ArithmeticError("an established non-isometry contradicts the supplied witness")
        forward = coordinate_identification(target, self._representatives[end]) * witness * coordinate_identification(self._representatives[start], source)
        self._adjacency.setdefault(start, []).append((end, forward))
        self._adjacency.setdefault(end, []).append((start, forward.inverse()))

    def cached_group(self, lattice: Lattice) -> tuple[Isometry, ...] | None:
        """Conjugate a computed complete O(L) generating family along a path.

        This is ``INDEF_FORM_AutomorphismGroup_Memoize`` (lines 1275-1310):
        the generators are transported by a proved isometry, with no attempt
        to infer isometry from agreement of discriminant or genus.
        """
        destination = self._key(lattice)
        for source, generators in self._groups.items():
            path = self._path(source, destination)
            if path is None:
                continue
            transport = coordinate_identification(self._representatives[destination], lattice) * path
            if self._representatives[source] is lattice:
                return generators
            inverse = transport.inverse()
            return tuple(transport * generator * inverse for generator in generators)
        return None

    def record_group(self, lattice: Lattice, generators: tuple[Isometry, ...]) -> None:
        """Store an accepted full automorphism generating family at its Gram key."""
        key = self._key(lattice)
        if any(generator.domain() is not lattice or generator.codomain() is not lattice for generator in generators):
            raise ValueError("a cached orthogonal generator must be an endomorphism of the given lattice")
        representative = self._representatives[key]
        if representative is lattice:
            self._groups[key] = generators
            return
        transport = coordinate_identification(lattice, representative)
        inverse = transport.inverse()
        self._groups[key] = tuple(transport * generator * inverse for generator in generators)
