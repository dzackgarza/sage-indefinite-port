"""Normed Dynkin extensions of a root configuration, for the Allcock edgewalk (ComputePossibleExtensions)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ

type CoxeterLabel = Integer | type(Infinity)

_FOUR_COS_SQUARED: dict[CoxeterLabel, Rational] = {Integer(2): QQ(0), Integer(3): QQ(1), Integer(4): QQ(2), Integer(6): QQ(3), Infinity: QQ(4)}
_FINITE_LABELS: tuple[CoxeterLabel, ...] = (Integer(2), Integer(3), Integer(4), Integer(6))


@dataclass(frozen=True)
class NormedDynkinExtension:
    r"""A root a of norm ``norm`` with prescribed Coxeter labels against the roots r_i.

    ``u_component`` is the projection of a onto the span of the r_i, in coordinates of L tensor Q,
    and ``residual_norm`` is q(a) - q(u_component), the norm left for the component orthogonal to them.
    """

    u_component: tuple[Rational, ...]
    residual_norm: Rational
    norm: Rational


def _norm_ratio_allows(label: CoxeterLabel, old_norm: Rational, new_norm: Rational) -> bool:
    if label == 3:
        return new_norm == old_norm
    if label == 4:
        return 2 * new_norm == old_norm or new_norm == 2 * old_norm
    if label == 6:
        return 3 * new_norm == old_norm or new_norm == 3 * old_norm
    return True


def _extensions(
    lattice: Lattices.ParentMethods, roots: Sequence[Lattices.ElementMethods], norms: Sequence[int], labels: tuple[CoxeterLabel, ...]
) -> tuple[NormedDynkinExtension, ...]:
    if not roots:
        return tuple(NormedDynkinExtension((QQ.zero(),) * int(lattice.rank()), QQ(norm), QQ(norm)) for norm in norms)
    coordinates = matrix(QQ, [[QQ(int(root.to_vector()(label))) for label in lattice.module_generating_set()] for root in roots])
    gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(QQ)
    pairings = coordinates * gram * coordinates.transpose()
    root_norms = tuple(pairings[i, i] for i in range(len(roots)))
    if not pairings.is_positive_definite():
        raise ValueError(f"the roots {roots} span a form that is not positive definite: they are not the walls through a point of hyperbolic space")
    inverse = pairings.inverse()
    if any(pairings[i, j] > 0 for i in range(len(roots)) for j in range(i)):
        raise ValueError("Dynkin extensions require simple roots with nonpositive mutual pairings")
    out: list[NormedDynkinExtension] = []

    def extend(target: tuple[Rational, ...], projected_square: Rational, norm: Rational, choices: Sequence[tuple[Rational, ...]]) -> None:
        position = len(target)
        if position == len(roots):
            weights = inverse * vector(QQ, target)
            out.append(NormedDynkinExtension(tuple(weights * coordinates), norm - projected_square, norm))
            return
        for pairing in choices[position]:
            square = projected_square + inverse[position, position] * pairing * pairing
            square += 2 * pairing * sum(inverse[position, j] * target[j] for j in range(position))
            # A positive-definite simple-root Gram matrix is a Stieltjes
            # matrix: its inverse is entrywise nonnegative. All chosen
            # pairings are nonpositive, so later terms cannot lower this
            # partial squared projection.
            if square <= norm:
                extend((*target, pairing), square, norm, choices)

    for new_norm in (QQ(norm) for norm in norms):
        choices = []
        for old_norm in root_norms:
            allowed = []
            for label in labels:
                square = _FOUR_COS_SQUARED[label] * old_norm * new_norm / 4
                if _norm_ratio_allows(label, old_norm, new_norm) and square.is_square():
                    allowed.append(-square.sqrt())
            choices.append(tuple(allowed))

        extend((), QQ.zero(), new_norm, choices)
    return tuple(out)


def finite_type_extensions(lattice: Lattices.ParentMethods, roots: Sequence[Lattices.ElementMethods], norms: Sequence[int]) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root at angles pi/2, pi/3, pi/4, pi/6 with a nonzero orthogonal component."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, _FINITE_LABELS) if extension.residual_norm != 0)


def parabolic_extensions(lattice: Lattices.ParentMethods, roots: Sequence[Lattices.ElementMethods], norms: Sequence[int]) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root lying in the span of the roots (residual norm 0), at any crystallographic angle including parallel walls."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, (*_FINITE_LABELS, Infinity)) if extension.residual_norm == 0)
