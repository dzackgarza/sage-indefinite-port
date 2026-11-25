r"""The even Clifford model of a nondegenerate rank-five rational quadratic space.

The quadratic map is ``q(v) = v^t B v``, with Clifford relation ``v^2 = q(v)``.
The volume element ``z`` of an orthogonal basis is central with ``z^2 = s`` for
``s`` the product of the ``q(e_i)``. With ``n_+(x) = -x^2`` on the five-dimensional
space ``z^{-1} V`` of even elements, ``(z^{-1} v)^2 = q(v) / s`` gives
``q(v) = -s n_+(z^{-1} v)``: this is Knus, Quadratic Forms, Clifford Algebras and
Spinors, equation (12.1), with ``lambda = -s`` (issue #57).

Integral transfer. In ``Lambda^2 ZZ^4 = 3U`` take ``Omega = e12 - e34``, of norm ``-2``; its
orthogonal complement has the basis ``(e13, -e24), (e14, e23), e12 + e34`` and is ``U + U + <2>``,
with ``Q = b / 2`` the pfaffian. ``Sp_4(ZZ)`` is the stabilizer of ``Omega`` in ``SL_4(ZZ)``. It acts on
``L = U + U + <2>``; the kernel is ``{1, -1}``, since a kernel element also fixes ``Omega`` and so
acts trivially on ``Lambda^2``. ``Sp_4 -> SO(L)`` is the spin covering. The Cartan form
``diag(p^a, p^b, p^-a, p^-b)`` acts on ``L tensor ZZ_p`` with weights ``0, +-a +-b``, so integrality gives
``a = b = 0`` and the image of ``Sp_4(ZZ)`` is exactly ``SO(L) & ker theta``. With ``Q = b / 2``,
``theta(-1) = 1 * -1 * 1 * -1 * 1 = 1`` (the orthogonal basis ``e_i +- f_i, t``), so ``-1`` has
``(det, theta) = (-1, 1)``. ``L`` is regular over ``ZZ`` for ``Q``, and ``theta`` takes values in ``{1, -1}``.
"""

from sage.algebras.clifford_algebra import CliffordAlgebra
from sage.matrix.constructor import matrix
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage_extensions.all import ZZ as OwnedZZ
from sage_extensions.all import Lattices

from sage_indefinite_port.indefinite.clifford_split_rank_six import (
    determinant_and_spinor_values,
    hyperbolic_coset_representatives,
    isometry_from_matrix,
    second_compound,
)

_SYMPLECTIC_BIVECTOR = matrix(QQ, [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
_SYMPLECTIC_FORM = _SYMPLECTIC_BIVECTOR.inverse()


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


def hyperbolic_hyperbolic_two():
    r"""Return ``U + U + <2>`` with its basis in the order of :func:`pfaffian_action`."""
    hyperbolic = Lattices(OwnedZZ)("U")
    return hyperbolic + hyperbolic + Lattices(OwnedZZ)("A1").twist(OwnedZZ(-1))


def pfaffian_action(g):
    r"""Return the action of ``g`` in ``Sp_4`` on the pfaffian basis ``(e13, -e24, e14, e23, e12 + e34)``."""
    action = second_compound(g)
    if action.column(0) - action.column(1) != matrix(QQ, 6, 1, [1, -1, 0, 0, 0, 0]).column(0):
        raise ValueError(f"{g} does not fix the symplectic bivector e12 - e34")
    columns = [action.column(2), action.column(3), action.column(4), action.column(5), action.column(0) + action.column(1)]
    return matrix(QQ, [[column[2], column[3], column[4], column[5], column[0]] for column in columns]).transpose()


def symplectic_generators():
    r"""Return the eight root elements of ``Sp_4(ZZ)`` for the bivector ``e12 - e34``.

    ``x -> x + omega(x, u) u`` is a long root element, and ``x -> x + omega(x, u) v + omega(x, v) u``
    for ``u, v`` in distinct hyperbolic pairs is a short one. These generate ``Sp_4(ZZ)`` (Hua--Reiner 1949).
    """
    basis = matrix.identity(QQ, 4).columns()

    def element(u, v):
        def image(x):
            shift = (x * _SYMPLECTIC_FORM * u) * v
            if u != v:
                shift += (x * _SYMPLECTIC_FORM * v) * u
            return x + shift

        return matrix(QQ, [list(image(x)) for x in basis]).transpose()

    long_roots = tuple(element(u, u) for u in basis)
    short_roots = tuple(element(u, v) for u in basis[:2] for v in basis[2:])
    return long_roots + short_roots


def spinor_kernel_generators(lattice):
    r"""Return generators of ``SO(L) & ker theta``, the image of ``Sp_4(ZZ)`` modulo ``{1, -1}``."""
    return tuple(isometry_from_matrix(pfaffian_action(g), lattice) for g in symplectic_generators())


def orthogonal_group_generators(lattice):
    r"""Return generators of the full ``O(L)`` for ``L = U + U + <2>``."""
    return spinor_kernel_generators(lattice) + tuple(hyperbolic_coset_representatives(lattice).values())


def characters(isometry, lattice):
    r"""Return ``(det, theta)``; see :func:`determinant_and_spinor_values`."""
    return determinant_and_spinor_values(isometry, lattice)
