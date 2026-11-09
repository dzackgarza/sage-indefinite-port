r"""The centre of the even Clifford algebra of a nondegenerate even-rank space.

For rank ``2m`` the volume element ``z`` of an orthogonal basis commutes with
``C_0`` and ``Z(C_0) = QQ[z]``. When ``z^2 = d`` is a square the centre is
``QQ x QQ`` and ``C_0`` splits along the idempotents ``(1 +- z / sqrt(d)) / 2``;
otherwise it is the quadratic field ``QQ(sqrt d)`` (Knus, Chapters 9 and 11; issue #57).
"""

from sage.algebras.clifford_algebra import CliffordAlgebra
from sage.matrix.constructor import matrix
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject


class EvenCliffordCentre(SageObject):
    r"""``Z(C_0(V, q))`` for ``q(v) = v^t B v`` with ``B`` of even rank."""

    def __init__(self, gram) -> None:
        self._gram = matrix(QQ, gram)
        rank = self._gram.nrows()
        if rank % 2 or not self._gram.is_symmetric() or self._gram.det() == 0:
            raise ValueError("the Gram matrix must be symmetric, nondegenerate and of even rank")
        form = QuadraticForm(QQ, 2 * self._gram)
        self._algebra = CliffordAlgebra(form)
        _, transform = form.rational_diagonal_form(return_matrix=True)
        generators = self._algebra.gens()
        self._volume = self._algebra.one()
        for column in transform.columns():
            self._volume = self._volume * sum(entry * g for entry, g in zip(column, generators, strict=True))
        square = self._volume * self._volume
        self._d = next(iter(square.monomial_coefficients().values()))
        if square != self._d * self._algebra.one():
            raise ValueError("the volume element does not square to a scalar")

    def _repr_(self) -> str:
        return f"Centre of the even Clifford algebra with Gram matrix {self._gram.list()}"

    def algebra(self):
        return self._algebra

    def volume(self):
        return self._volume

    def volume_square(self):
        return self._d

    def is_split(self):
        r"""Return whether the discriminant algebra is ``QQ x QQ``."""
        return self._d.is_square()

    def discriminant_field(self):
        r"""Return ``QQ(sqrt d)`` when the discriminant algebra is a field."""
        if self.is_split():
            raise ValueError("the discriminant algebra is split, so it is not a field")
        return QQ["x"]("x^2 - " + str(self._d)).splitting_field("a")

    def idempotents(self):
        r"""Return the two central idempotents of the split case."""
        if not self.is_split():
            raise ValueError("the discriminant algebra is a field, so C_0 has no central idempotents")
        root = self._d.sqrt()
        one = self._algebra.one()
        return ((one + self._volume / root) / 2, (one - self._volume / root) / 2)

    def _even_keys(self):
        return tuple(key for key in self._algebra.basis().keys() if len(key) % 2 == 0)

    def component_dimension(self, idempotent):
        r"""Return the ``QQ``-dimension of the two-sided component ``e C_0``."""
        keys = self._even_keys()
        positions = {key: index for index, key in enumerate(keys)}
        action = matrix(QQ, len(keys), len(keys))
        for column, key in enumerate(keys):
            for image_key, value in (idempotent * self._algebra.basis()[key]).monomial_coefficients().items():
                action[positions[image_key], column] = value
        return action.rank()

    def volume_is_central_in_even_part(self):
        r"""Return whether ``z`` commutes with every even basis element."""
        return all(self._volume * self._algebra.basis()[key] == self._algebra.basis()[key] * self._volume for key in self._even_keys())
