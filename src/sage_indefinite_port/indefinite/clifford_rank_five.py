r"""The even Clifford model of a nondegenerate rank-five rational quadratic space.

The quadratic map is ``q(v) = v^t B v``, with Clifford relation ``v^2 = q(v)``.
The volume element ``z`` of an orthogonal basis is central with ``z^2 = s`` for
``s`` the product of the ``q(e_i)``. With ``n_+(x) = -x^2`` on the five-dimensional
space ``z^{-1} V`` of even elements, ``(z^{-1} v)^2 = q(v) / s`` gives
``q(v) = -s n_+(z^{-1} v)``: this is Knus, Quadratic Forms, Clifford Algebras and
Spinors, equation (12.1), with ``lambda = -s`` (issue #57).
"""

from sage.algebras.clifford_algebra import CliffordAlgebra
from sage.matrix.constructor import matrix
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject


class RankFiveCliffordModel(SageObject):
    r"""The algebra ``C_0(V, q)`` with the pfaffian space of Knus (12.1)."""

    def __init__(self, gram) -> None:
        self._gram = matrix(QQ, gram)
        if self._gram.nrows() != 5 or not self._gram.is_symmetric() or self._gram.det() == 0:
            raise ValueError("the Gram matrix must be symmetric, nondegenerate and of rank five")
        form = QuadraticForm(QQ, 2 * self._gram)
        self._algebra = CliffordAlgebra(form)
        _, transform = form.rational_diagonal_form(return_matrix=True)
        generators = self._algebra.gens()
        orthogonal = [sum(entry * g for entry, g in zip(column, generators, strict=True)) for column in transform.columns()]
        self._volume = orthogonal[0] * orthogonal[1] * orthogonal[2] * orthogonal[3] * orthogonal[4]
        self._s = self._scalar(self._volume * self._volume)
        self._volume_inverse = self._volume / self._s

    def _repr_(self) -> str:
        return f"Even Clifford model of the rank-five form with Gram matrix {self._gram.list()}"

    def _scalar(self, element):
        coefficients = element.monomial_coefficients()
        value = next(iter(coefficients.values()), QQ.zero())
        if element != value * self._algebra.one():
            raise ValueError(f"{element} is not a scalar")
        return value

    def algebra(self):
        return self._algebra

    def volume_square(self):
        r"""Return ``s = z^2``."""
        return self._s

    def scalar(self):
        r"""Return ``lambda = -s``."""
        return -self._s

    def pfaffian_basis(self):
        r"""Return the images ``z^{-1} e_i`` of the original basis."""
        return tuple(self._volume_inverse * g for g in self._algebra.gens())

    def pfaffian_gram(self):
        r"""Return the Gram matrix of the polar form of ``n_+`` on :meth:`pfaffian_basis`."""
        basis = self.pfaffian_basis()
        return matrix(QQ, [[-self._scalar(x * y + y * x) / 2 for y in basis] for x in basis])

    def recovered_gram(self):
        r"""Return ``lambda N``, which equals the original Gram matrix."""
        return self.scalar() * self.pfaffian_gram()

    def even_dimension(self):
        return sum(1 for key in self._algebra.basis().keys() if len(key) % 2 == 0)

    def even_symmetric_dimension(self):
        r"""Return the dimension of the fixed space of the reversal on ``C_0``."""
        keys = tuple(key for key in self._algebra.basis().keys() if len(key) % 2 == 0)
        positions = {key: index for index, key in enumerate(keys)}
        difference = matrix(QQ, len(keys), len(keys))
        for column, key in enumerate(keys):
            element = self._algebra.basis()[key]
            for image_key, value in (element.transpose() - element).monomial_coefficients().items():
                difference[positions[image_key], column] = value
        return len(keys) - difference.rank()
