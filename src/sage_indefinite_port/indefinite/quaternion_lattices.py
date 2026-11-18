r"""Norm lattices of quaternion orders and the isometries given by norm-one units.

For integers ``a, b`` the order ``ZZ<1, i, j, k>`` of ``(a, b)_QQ`` with the form
``q = s Nrd`` has the diagonal Gram matrix ``s diag(1, -a, -b, a b)``. Left and
right multiplication by a unit of reduced norm one preserve ``Nrd`` and the order,
so they are isometries of this lattice (issues #40 and #57).
"""

from itertools import product

from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage_extensions.all import ZZ as OwnedZZ
from sage_extensions.all import Lattices


class QuaternionNormLattice(SageObject):
    r"""The lattice ``(ZZ<1, i, j, k>, s Nrd)`` of the quaternion algebra ``(a, b)``."""

    def __init__(self, a: int, b: int, scale: int) -> None:
        self._algebra = QuaternionAlgebra(QQ, a, b)
        self._basis = (self._algebra.one(), *self._algebra.gens())
        self._gram = [[int(scale * (x * y.conjugate()).reduced_trace() / 2) for y in self._basis] for x in self._basis]
        self._lattice = Lattices(OwnedZZ)(self._gram)

    def _repr_(self) -> str:
        return f"Norm lattice of the order ZZ<1, i, j, k> of {self._algebra}"

    def algebra(self):
        return self._algebra

    def lattice(self):
        return self._lattice

    def element(self, *coordinates):
        r"""Return the quaternion with the given integer coordinates on ``1, i, j, k``."""
        return sum(int(c) * b for c, b in zip(coordinates, self._basis, strict=True))

    def _multiplication_matrix(self, unit, side: str):
        if unit.reduced_norm() != 1:
            raise ValueError(f"{unit} does not have reduced norm one")
        products = {"left": lambda x: unit * x, "right": lambda x: x * unit}[side]
        columns = matrix(QQ, [list(products(x)) for x in self._basis]).transpose()
        if columns.denominator() != 1:
            raise ValueError(f"{unit} does not lie in the order ZZ<1, i, j, k>")
        return columns

    def left_multiplication_matrix(self, unit):
        return self._multiplication_matrix(unit, "left")

    def right_multiplication_matrix(self, unit):
        return self._multiplication_matrix(unit, "right")

    def isometry(self, unit, side: str):
        r"""Return left or right multiplication by ``unit`` as an element of ``O(L)``."""
        columns = self._multiplication_matrix(unit, side)
        generators = tuple(self._lattice.module_generators())
        labels = tuple(self._lattice.module_generating_set())
        images = (sum((int(entry) * point for entry, point in zip(column, generators, strict=True)), self._lattice.zero()) for column in columns.columns())
        return self._lattice.O()(dict(zip(labels, images, strict=True)))

    def _as_isometry(self, columns):
        if columns.denominator() != 1:
            raise ValueError("the map does not preserve the order ZZ<1, i, j, k>")
        generators = tuple(self._lattice.module_generators())
        labels = tuple(self._lattice.module_generating_set())
        images = (sum((int(entry) * point for entry, point in zip(column, generators, strict=True)), self._lattice.zero()) for column in columns.columns())
        return self._lattice.O()(dict(zip(labels, images, strict=True)))

    def pair_matrix(self, u, v):
        r"""Return the matrix of ``x -> u x conj(v)`` for ``nrd(u) nrd(v) = 1``."""
        if u.reduced_norm() * v.reduced_norm() != 1:
            raise ValueError("the pair does not scale the norm by one")
        return matrix(QQ, [list(u * x * v.conjugate()) for x in self._basis]).transpose()

    def pair_isometry(self, u, v):
        return self._as_isometry(self.pair_matrix(u, v))

    def conjugation_matrix(self):
        return matrix(QQ, [list(x.conjugate()) for x in self._basis]).transpose()

    def conjugation_isometry(self):
        r"""Return quaternion conjugation, which has determinant ``-1``."""
        return self._as_isometry(self.conjugation_matrix())

    def pair_spinor_class(self, u):
        r"""Return the squarefree part of ``nrd(u)``, the spinor class of the pair action of ``u``."""
        return u.reduced_norm().squarefree_part()

    def _coordinates_are_integral(self, element):
        return all(c.denominator() == 1 for c in element)

    def normalizes_order(self, element):
        r"""Return whether conjugation by ``element`` preserves ``ZZ<1, i, j, k>``."""
        inverse = ~element
        return all(self._coordinates_are_integral(element * x * inverse) for x in self._basis)

    def _same_normalizer_class(self, first, second):
        ratio = first * ~second
        content = gcd([c.numerator() for c in ratio]) / lcm([c.denominator() for c in ratio])
        primitive = ratio / content
        return self._coordinates_are_integral(primitive) and abs(primitive.reduced_norm()) == 1

    def normalizer_classes_in_box(self, radius: int):
        r"""Return distinct classes of ``N(Lambda) / QQ^x Lambda^x`` with representatives in a coordinate box."""
        classes = []
        for point in product(range(-radius, radius + 1), repeat=4):
            element = self.element(*point)
            if element == 0 or element.reduced_norm() == 0 or not self.normalizes_order(element):
                continue
            if not any(self._same_normalizer_class(element, known) for known in classes):
                classes.append(element)
        return tuple(classes)

    def conjugation_by_matrix(self, element):
        r"""Return the matrix of ``x -> element x element^{-1}``."""
        inverse = ~element
        return matrix(QQ, [list(element * x * inverse) for x in self._basis]).transpose()

    def conjugation_by_isometry(self, element):
        r"""Return the isometry ``x -> element x element^{-1}`` for ``element`` normalizing the order."""
        return self._as_isometry(self.conjugation_by_matrix(element))

