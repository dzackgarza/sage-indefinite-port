"""Oriented binary cycles for Allcock's edgewalk (arXiv:2110.03784, section 8).

The lattice supplies its form and semantic isometries.  Canonical supplements,
the ordered search and the period stopping rule belong to the edgewalk.
"""

from __future__ import annotations

from sage_extensions.all import Lattices
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ

from .edgewalk_scalars import _rational

type _Vector = Lattices.ElementMethods
type _Point = tuple[int, int]
type _Pair = tuple[_Point, _Point]


class _OrientedBinaryCycle:
    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        if int(lattice.rank()) != 2 or tuple(int(index) for index in lattice.signature_pair()) != (1, 1) or lattice.is_isotropic():
            raise ValueError("an oriented binary cycle needs an anisotropic lattice of signature (1,1)")
        self.lattice = lattice
        # The form coefficients and every coordinate conversion must use the
        # same ordered basis. `module_generators()` returns elements, whereas
        # the coordinate routines below are indexed by `module_generating_set()`;
        # derive the basis elements from those labels rather than assuming the
        # two independently exposed iteration orders coincide.
        self.labels = tuple(lattice.module_generating_set())
        if len(self.labels) != 2:
            raise ArithmeticError("a rank-two binary cycle needs exactly two basis labels")
        first, second = tuple(lattice.module_generator(label) for label in self.labels)
        self.form = int(first.q()), int(lattice.b(first, second)), int(second.q())

    def coordinates(self, point: _Vector) -> tuple[int, ...]:
        coordinates = point.to_vector()
        values = tuple(_rational(coordinates(label)) for label in self.labels)
        if any(value.denominator() != 1 for value in values):
            raise ArithmeticError("an integral binary lattice vector has nonintegral coordinates")
        return tuple(int(value) for value in values)

    def point(self, coordinates: _Point) -> _Vector:
        ring = self.lattice.base_ring()
        return self.lattice.linear_combination({label: ring(int(entry)) for label, entry in zip(self.labels, coordinates, strict=True)})

    def pairing(self, left: _Point, right: _Point) -> int:
        a, b, c = self.form
        return a * left[0] * right[0] + b * (left[0] * right[1] + left[1] * right[0]) + c * left[1] * right[1]

    def canonical_companion(self, bound: int, right: _Point, left: _Point) -> _Point:
        """The last supplement before the positive branch of the bound hyperbola."""
        rr = self.pairing(right, right)
        rl = self.pairing(right, left)
        ll = self.pairing(left, left)
        radicand = rl * rl - rr * (ll - bound)
        if rr <= 0 or radicand < 0:
            raise ArithmeticError("canonical companion requires a positive right norm and real hyperbola intersection")
        root = SageZZ(int(radicand)).isqrt()
        multiple = (-rl + int(root)) // rr
        return left[0] + multiple * right[0], left[1] + multiple * right[1]

    def promised_step(self, bound: int, right: _Point, left: _Point) -> _Pair:
        """Find the next primitive positive vector of square at most bound.

        The oriented unimodular pair brackets the required vector. Each step
        decreases the sum of its nonnegative coordinates in that pair, as in
        Allcock's Promised algorithm. No iteration cap decides nonexistence.
        """
        if self.pairing(right, right) <= 0 or right[0] * left[1] - right[1] * left[0] != 1:
            raise ValueError("Promised requires a positive right vector and oriented unimodular pair")
        # An M-supplement lies outside the *interior* of the future
        # hyperbola's convex hull.  A vector of square > M in the opposite
        # time component is still a valid supplement; only a future vector
        # of square > M lies strictly inside that hull (Allcock, section 8).
        if self.pairing(left, left) > bound and self.pairing(right, left) > 0:
            raise ValueError("Promised requires an M-supplement outside the future hyperbola's convex hull")
        while True:
            middle = right[0] + left[0], right[1] + left[1]
            # The pinned C++ tests the pairing with the *right* endpoint
            # (two_dim_lorentzian.h:124, scal_mr). Allcock's printed step
            # [Go right?] (main.tex:1501) tests m.l instead, and that variant
            # returns a vector above the prescribed bound on q = -3x^2-6xy-y^2,
            # r0 = (-1,1), M = 2: it emits (-7,4) with q = 5. The C++ reading
            # is the one the paper's own proof uses, so follow the source.
            if self.pairing(middle, middle) <= bound or self.pairing(middle, right) < 0:
                left = middle
                continue
            if self.pairing(left, left) >= 0 and self.pairing(right, left) > 0:
                if self.pairing(left, left) > bound:
                    raise ArithmeticError("Promised returned a vector above its norm bound")
                return left, (-right[0], -right[1])
            # Allcock, proof following Algorithm Promised and the
            # performance recurrence: the left sector has basis (m,l).
            # The printed (m,r) in step [Go left] (main.tex:1508) contradicts
            # both the proof and the pinned C++ (two_dim_lorentzian.h:134,
            # Promised(G, M, m, l)); with it the returned period has
            # determinant -1 on the specimen above.
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

    def start_point(self, start: _Vector) -> _Point:
        if start.parent() is not self.lattice or int(start.content()) != 1 or start.q() <= 0:
            raise ValueError("the cycle start must be a primitive positive vector of its lattice")
        coordinates = self.coordinates(start)
        return coordinates[0], coordinates[1]

    def cycle(self, bound: int, start: _Vector) -> tuple[LatticeIsometryMethods, tuple[_Vector, ...]] | None:
        initial_point = self.start_point(start)
        if bound <= 0:
            return None
        _gcd, first, second = SageZZ(initial_point[0]).xgcd(SageZZ(initial_point[1]))
        companion = -int(second), int(first)
        # Shorter requires the q(start)-supplement, even when bound < q(start).
        companion = self.canonical_companion(self.pairing(initial_point, initial_point), initial_point, companion)
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
                rational_period = frame.inverse() * final
                if any(entry.denominator() != 1 for entry in rational_period.list()):
                    raise ArithmeticError("the binary period does not preserve the integral lattice")
                period = rational_period.change_ring(SageZZ)
                gram = matrix(SageZZ, ((self.form[0], self.form[1]), (self.form[1], self.form[2])))
                if period * gram * period.transpose() != gram or period.det() != 1:
                    raise ArithmeticError("the binary period is not an orientation-preserving isometry")
                return self.result(period, tuple(representatives))

    def result(self, period: Matrix_integer_dense, representatives: tuple[_Point, ...]) -> tuple[LatticeIsometryMethods, tuple[_Vector, ...]]:
        """The period as an isometry of the lattice and the representatives as its vectors."""
        # The period rows are the images of the two coordinate basis
        # vectors in `self.labels` order. The Isom constructor, however,
        # consumes images in `module_generators()` order. Translate each
        # actual generator through its coordinates rather than assuming
        # these two iteration orders coincide.
        images = []
        for generator in self.lattice.module_generators():
            coordinates = self.coordinates(generator)
            image_coordinates = tuple(sum(coordinates[row] * int(period[row, column]) for row in range(2)) for column in range(2))
            images.append(self.point(image_coordinates))
        return self.lattice.Isom(self.lattice)(tuple(images)), tuple(self.point(point) for point in representatives)


class _FormBinaryCycle(_OrientedBinaryCycle):
    """The oriented cycle of an edge plane read on its Gram coefficients (edgewalk.h:772, Anisotropic(Gwork, ...)).

    The walk consumes only coordinates: the period as its integral row
    matrix and the representatives as coordinate pairs, the data the
    pinned source returns. The lattice-backed parent class wraps the same
    result into the owned isometry and vectors.
    """

    def __init__(self, form: tuple[int, int, int]) -> None:
        a, b, c = form
        if a * c - b * b >= 0 or SageZZ(b * b - a * c).is_square():
            raise ValueError("an oriented binary cycle needs an anisotropic form of signature (1,1)")
        self.form = form

    def start_point(self, start: _Point) -> _Point:
        if SageZZ(start[0]).gcd(SageZZ(start[1])) != 1 or self.pairing(start, start) <= 0:
            raise ValueError("the cycle start must be a primitive positive vector of its lattice")
        return start

    def result(self, period: Matrix_integer_dense, representatives: tuple[_Point, ...]) -> tuple[Matrix_integer_dense, tuple[_Point, ...]]:
        return period, representatives


def anisotropic_cycle(form: tuple[int, int, int], bound: Rational, start: _Point) -> tuple[Matrix_integer_dense, tuple[_Point, ...]] | None:
    """The positive oriented period (integral row matrix) and all primitive representatives after start."""
    return _FormBinaryCycle(form).cycle(int(_rational(bound).floor()), start)
