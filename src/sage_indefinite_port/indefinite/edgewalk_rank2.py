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
type _Point = tuple[int, int]
type _Pair = tuple[_Point, _Point]


class _OrientedBinaryCycle:
    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        if int(lattice.rank()) != 2 or tuple(int(index) for index in lattice.signature_pair()) != (1, 1) or lattice.is_isotropic():
            raise ValueError("an oriented binary cycle needs an anisotropic lattice of signature (1,1)")
        self.lattice = lattice
        first, second = tuple(lattice.module_generators())
        self.form = int(first.q()), int(lattice.b(first, second)), int(second.q())

    def coordinates(self, point: _Vector) -> tuple[int, ...]:
        coordinates = point.to_vector()
        return tuple(int(coordinates(label)) for label in self.lattice.module_generating_set())

    def point(self, coordinates: _Point) -> _Vector:
        ring = self.lattice.base_ring()
        return self.lattice.linear_combination({label: ring(entry) for label, entry in zip(self.lattice.module_generating_set(), coordinates, strict=True)})

    def pairing(self, left: _Point, right: _Point) -> int:
        a, b, c = self.form
        return a * left[0] * right[0] + b * (left[0] * right[1] + left[1] * right[0]) + c * left[1] * right[1]

    def canonical_companion(self, bound: int, right: _Point, left: _Point) -> _Point:
        """The last supplement before the positive branch of the bound hyperbola."""
        rr = self.pairing(right, right)
        rl = self.pairing(right, left)
        ll = self.pairing(left, left)
        radicand = rl * rl - rr * (ll - bound)
        root = SageZZ(radicand).isqrt()
        multiple = (-rl + int(root)) // rr
        return left[0] + multiple * right[0], left[1] + multiple * right[1]

    def promised_step(self, bound: int, right: _Point, left: _Point) -> _Pair:
        """Find the next primitive positive vector of square at most bound.

        The oriented unimodular pair brackets the required vector. Each step
        decreases the sum of its nonnegative coordinates in that pair, as in
        Allcock's Promised algorithm. No iteration cap decides nonexistence.
        """
        while True:
            middle = right[0] + left[0], right[1] + left[1]
            # Allcock, Algorithm Promised, step [Go right?]: the test is
            # m·l < 0, not m·r < 0. The negative pairing certifies that
            # the left sector contains no eligible next vector.
            if self.pairing(middle, middle) <= bound or self.pairing(middle, left) < 0:
                left = middle
                continue
            if self.pairing(left, left) >= 0 and self.pairing(right, left) > 0:
                return left, (-right[0], -right[1])
            # Allcock, proof following Algorithm Promised and the
            # performance recurrence: the left sector has basis (m,l).
            # The printed (m,r) in step [Go left] contradicts both.
            right = middle

    def characteristic(self, right: _Point, left: _Point) -> tuple[int, int, int]:
        return self.pairing(right, right), self.pairing(right, left), self.pairing(left, left)

    def shorter_pair(self, right: _Point, left: _Point) -> _Pair | None:
        """A strictly shorter positive vector, or a full period proving absence."""
        initial = self.characteristic(right, left)
        while True:
            right, left = self.promised_step(self.pairing(right, right), right, left)
            left = self.canonical_companion(self.pairing(right, right), right, left)
            if self.pairing(right, right) < initial[0]:
                return right, left
            if self.characteristic(right, left) == initial:
                return None

    def reduced_start_pair(self, bound: int, right: _Point, left: _Point) -> _Pair | None:
        """Allcock's NotPromised, beginning with a canonical q(right)-supplement."""
        if self.pairing(right, right) <= bound:
            right, left = self.promised_step(bound, right, left)
            return right, self.canonical_companion(bound, right, left)
        while True:
            pair = self.shorter_pair(right, left)
            if pair is None:
                return None
            right, left = pair
            if self.pairing(right, right) <= bound:
                return right, self.canonical_companion(bound, right, left)

    def cycle(self, bound: int, start: _Vector) -> tuple[LatticeIsometryMethods, tuple[_Vector, ...]] | None:
        if start.parent() is not self.lattice or int(start.content()) != 1 or start.q() <= 0:
            raise ValueError("the cycle start must be a primitive positive vector of its lattice")
        if bound <= 0:
            return None
        coordinates = self.coordinates(start)
        initial_point = coordinates[0], coordinates[1]
        _gcd, first, second = SageZZ(coordinates[0]).xgcd(SageZZ(coordinates[1]))
        companion = -int(second), int(first)
        # Shorter requires the q(start)-supplement, even when bound < q(start).
        companion = self.canonical_companion(int(start.q()), initial_point, companion)
        pair = self.reduced_start_pair(bound, initial_point, companion)
        if pair is None:
            return None
        right, left = pair
        initial = self.characteristic(right, left)
        frame = matrix(SageQQ, pair)
        representatives: list[_Point] = []
        while True:
            representatives.append(right)
            right, left = self.promised_step(bound, right, left)
            left = self.canonical_companion(bound, right, left)
            if self.characteristic(right, left) == initial:
                final = matrix(SageQQ, (right, left))
                period = (frame.inverse() * final).change_ring(SageZZ)
                images = tuple(self.point((int(row[0]), int(row[1]))) for row in period.rows())
                return self.lattice.Isom(self.lattice)(images), tuple(self.point(point) for point in representatives)


def anisotropic_cycle(lattice: Lattices.ParentMethods, bound: Rational, start: _Vector) -> tuple[LatticeIsometryMethods, tuple[_Vector, ...]] | None:
    """The positive oriented period and all primitive representatives after start."""
    return _OrientedBinaryCycle(lattice).cycle(int(SageQQ(bound).floor()), start)
