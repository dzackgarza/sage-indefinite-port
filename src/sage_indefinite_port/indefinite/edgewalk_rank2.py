"""Oriented binary cycles for Allcock's edgewalk (arXiv:2110.03784, section 8).

The lattice supplies its form and semantic isometries.  Canonical supplements,
the ordered search and the period stopping rule belong to the edgewalk.
"""

from __future__ import annotations

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ

type _Vector = Lattices.ElementMethods
type _Pair = tuple[_Vector, _Vector]


class _OrientedBinaryCycle:
    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        if int(lattice.rank()) != 2 or tuple(int(index) for index in lattice.signature_pair()) != (1, 1) or lattice.is_isotropic():
            raise ValueError("an oriented binary cycle needs an anisotropic lattice of signature (1,1)")
        self.lattice = lattice

    def coordinates(self, point: _Vector) -> tuple[int, ...]:
        coordinates = point.to_vector()
        return tuple(int(coordinates(label)) for label in self.lattice.module_generating_set())

    def canonical_companion(self, bound: int, right: _Vector, left: _Vector) -> _Vector:
        """The last supplement before the positive branch of the bound hyperbola."""
        rr = int(right.q())
        rl = int(self.lattice.b(right, left))
        ll = int(left.q())
        radicand = rl * rl - rr * (ll - bound)
        root = SageZZ(radicand).isqrt()
        multiple = (-rl + int(root)) // rr
        return left + self.lattice.scalar_multiple(self.lattice.base_ring()(multiple), right)

    def promised_step(self, bound: int, right: _Vector, left: _Vector) -> _Pair:
        """Find the next primitive positive vector of square at most bound.

        The oriented unimodular pair brackets the required vector. Each step
        decreases the sum of its nonnegative coordinates in that pair, as in
        Allcock's Promised algorithm. No iteration cap decides nonexistence.
        """
        while True:
            middle = right + left
            if int(middle.q()) <= bound or self.lattice.b(middle, right) < 0:
                left = middle
                continue
            if left.q() >= 0 and self.lattice.b(right, left) > 0:
                return left, -right
            right = middle

    def characteristic(self, right: _Vector, left: _Vector) -> tuple[int, int, int]:
        return int(right.q()), int(self.lattice.b(right, left)), int(left.q())

    def shorter_pair(self, right: _Vector, left: _Vector) -> _Pair | None:
        """A strictly shorter positive vector, or a full period proving absence."""
        initial = self.characteristic(right, left)
        while True:
            right, left = self.promised_step(int(right.q()), right, left)
            left = self.canonical_companion(int(right.q()), right, left)
            if int(right.q()) < initial[0]:
                return right, left
            if self.characteristic(right, left) == initial:
                return None

    def reduced_start_pair(self, bound: int, right: _Vector, left: _Vector) -> _Pair | None:
        """Allcock's NotPromised, beginning with a canonical q(right)-supplement."""
        if int(right.q()) <= bound:
            right, left = self.promised_step(bound, right, left)
            return right, self.canonical_companion(bound, right, left)
        while True:
            pair = self.shorter_pair(right, left)
            if pair is None:
                return None
            right, left = pair
            if int(right.q()) <= bound:
                return right, self.canonical_companion(bound, right, left)

    def cycle(self, bound: int, start: _Vector) -> tuple[LatticeIsometryMethods, tuple[_Vector, ...]] | None:
        if start.parent() is not self.lattice or int(start.content()) != 1 or start.q() <= 0:
            raise ValueError("the cycle start must be a primitive positive vector of its lattice")
        if bound <= 0:
            return None
        coordinates = self.coordinates(start)
        _gcd, first, second = SageZZ(coordinates[0]).xgcd(SageZZ(coordinates[1]))
        labels = tuple(self.lattice.module_generating_set())
        ring = self.lattice.base_ring()
        companion = self.lattice.linear_combination({labels[0]: ring(-int(second)), labels[1]: ring(int(first))})
        # Shorter requires the q(start)-supplement, even when bound < q(start).
        companion = self.canonical_companion(int(start.q()), start, companion)
        pair = self.reduced_start_pair(bound, start, companion)
        if pair is None:
            return None
        right, left = pair
        initial = self.characteristic(right, left)
        frame = matrix(SageQQ, [self.coordinates(point) for point in pair])
        representatives: list[_Vector] = []
        while True:
            representatives.append(right)
            right, left = self.promised_step(bound, right, left)
            left = self.canonical_companion(bound, right, left)
            if self.characteristic(right, left) == initial:
                final = matrix(SageQQ, [self.coordinates(point) for point in (right, left)])
                period = (frame.inverse() * final).change_ring(SageZZ)
                images = tuple(self.lattice.linear_combination({label: ring(int(entry)) for label, entry in zip(labels, row, strict=True)}) for row in period.rows())
                return self.lattice.Isom(self.lattice)(images), tuple(representatives)


def anisotropic_cycle(lattice: Lattices.ParentMethods, bound: Rational, start: _Vector) -> tuple[LatticeIsometryMethods, tuple[_Vector, ...]] | None:
    """The positive oriented period and all primitive representatives after start."""
    return _OrientedBinaryCycle(lattice).cycle(int(SageQQ(bound).floor()), start)
