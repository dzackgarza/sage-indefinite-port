r"""Norm lattices of quaternion orders and the isometries given by norm-one units.

For integers ``a, b`` the order ``ZZ<1, i, j, k>`` of ``(a, b)_QQ`` with the form
``q = s Nrd`` has the diagonal Gram matrix ``s diag(1, -a, -b, a b)``. Left and
right multiplication by a unit of reduced norm one preserve ``Nrd`` and the order,
so they are isometries of this lattice (issues #40 and #57).
"""

from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
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
