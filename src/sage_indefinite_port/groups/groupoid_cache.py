"""Exact isometry-groupoid cache for lattice presentations."""

from __future__ import annotations

from dataclasses import dataclass, field

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices

from sage_indefinite_port.backends.canonization import presentation_bucket_key


type BucketKey = tuple[int, tuple[int, int], bool, int]
type LatticePair = tuple[Lattices.ParentMethods, Lattices.ParentMethods]


@dataclass
class IsometryGroupoidCache:
    """Remember exact arrows and automorphism generators between lattice objects."""

    _isometries: dict[LatticePair, LatticeIsometryMethods] = field(default_factory=dict)
    _nonisometric: set[LatticePair] = field(default_factory=set)
    _orthogonal_groups: dict[Lattices.ParentMethods, tuple[LatticeIsometryMethods, ...]] = field(default_factory=dict)
    _buckets: dict[BucketKey, set[Lattices.ParentMethods]] = field(default_factory=dict)

    def _remember_bucket(self, lattice: Lattices.ParentMethods) -> None:
        self._buckets.setdefault(presentation_bucket_key(lattice), set()).add(lattice)

    def remember_isometry(self, isometry: LatticeIsometryMethods) -> None:
        source = isometry.domain()
        target = isometry.codomain()
        self._remember_bucket(source)
        self._remember_bucket(target)
        self._isometries[(source, target)] = isometry
        self._isometries[(target, source)] = ~isometry
        self._nonisometric.discard((source, target))
        self._nonisometric.discard((target, source))

    def remember_nonisometric(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> None:
        if presentation_bucket_key(source) != presentation_bucket_key(target):
            return
        self._remember_bucket(source)
        self._remember_bucket(target)
        self._nonisometric.add((source, target))
        self._nonisometric.add((target, source))
        self._isometries.pop((source, target), None)
        self._isometries.pop((target, source), None)

    def lookup_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods):
        if presentation_bucket_key(source) != presentation_bucket_key(target):
            return None
        if (source, target) in self._nonisometric:
            return None
        return self._isometries.get((source, target))

    def remember_orthogonal_group(
        self,
        lattice: Lattices.ParentMethods,
        generators: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        automorphisms = lattice.Aut()
        verified = tuple(automorphisms(generator) for generator in generators)
        self._remember_bucket(lattice)
        self._orthogonal_groups[lattice] = verified

    def lookup_orthogonal_group(self, lattice: Lattices.ParentMethods):
        direct = self._orthogonal_groups.get(lattice)
        if direct is not None:
            return direct
        key = presentation_bucket_key(lattice)
        for cached in self._buckets.get(key, ()):
            generators = self._orthogonal_groups.get(cached)
            if generators is None:
                continue
            isometry = self.lookup_isometry(lattice, cached)
            if isometry is None:
                continue
            inverse = ~isometry
            automorphisms = lattice.Aut()
            transported = tuple(automorphisms(inverse * generator * isometry) for generator in generators)
            self._orthogonal_groups[lattice] = transported
            return transported
        return None


__all__ = ["IsometryGroupoidCache"]
