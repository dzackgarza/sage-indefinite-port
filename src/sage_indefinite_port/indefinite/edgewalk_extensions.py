"""Normed Dynkin extensions of a root configuration, for the Allcock edgewalk (ComputePossibleExtensions)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from dzack_research.preamble.all import CoxeterDiagrams, Lattices
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ

type CoxeterLabel = Integer | type(Infinity)

_FOUR_COS_SQUARED: dict[CoxeterLabel, Rational] = {
    Integer(2): QQ(0),
    Integer(3): QQ(1),
    Integer(4): QQ(2),
    Integer(6): QQ(3),
    Infinity: QQ(4),
}
_FINITE_LABELS: tuple[CoxeterLabel, ...] = (
    Integer(2),
    Integer(3),
    Integer(4),
    Integer(6),
)


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
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
    labels: tuple[CoxeterLabel, ...],
) -> tuple[NormedDynkinExtension, ...]:
    if not roots:
        return tuple(NormedDynkinExtension((QQ.zero(),) * int(lattice.rank()), QQ(norm), QQ(norm)) for norm in norms)
    coordinates = matrix(
        QQ,
        [[QQ(int(root.to_vector()(label))) for label in lattice.module_generating_set()] for root in roots],
    )
    gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(QQ)
    pairings = coordinates * gram * coordinates.transpose()
    root_norms = tuple(pairings[i, i] for i in range(len(roots)))
    if not pairings.is_positive_definite():
        raise ValueError(f"the roots {roots} span a form that is not positive definite: they are not the walls through a point of hyperbolic space")
    inverse = pairings.inverse()
    if any(pairings[i, j] > 0 for i in range(len(roots)) for j in range(i)):
        raise ValueError("Dynkin extensions require simple roots with nonpositive mutual pairings")
    # ComputePossibleExtensions (coxeter_dynkin.h:2305-2336) iterates diagram
    # extensions before prescribed norms. The source enumerates diagrams
    # by Dynkin type (FindDiagramExtensions_Efficient:715-1671), rather than
    # lexicographic label order; the finite diagram check is indispensable.
    out: list[tuple[tuple[CoxeterLabel, ...], int, NormedDynkinExtension]] = []
    # The preamble owns the Coxeter diagram and its exact Schlaefli
    # classification. The pinned CheckDiagram recognizes the same spherical
    # and Euclidean types before any root-norm test.
    coxeter_rows = [[Integer(1) if i == j else Integer(2) for j in range(len(roots) + 1)] for i in range(len(roots) + 1)]
    for i in range(len(roots)):
        for j in range(i):
            old_pair = pairings[i, j]
            squares = 4 * old_pair * old_pair / (root_norms[i] * root_norms[j])
            label = next((edge for edge, value in _FOUR_COS_SQUARED.items() if squares == value), None)
            if label is None:
                raise ValueError("the supplied roots have noncrystallographic Coxeter pairings")
            coxeter_rows[i][j] = coxeter_rows[j][i] = Integer(-1) if label == Infinity else label

    def diagram_allowed(diagram: tuple[CoxeterLabel, ...]) -> bool:
        for i, label in enumerate(diagram):
            coxeter_rows[i][-1] = coxeter_rows[-1][i] = Integer(-1) if label == Infinity else label
        extended = CoxeterDiagrams()(coxeter_rows)
        return extended.is_elliptic() if len(labels) == len(_FINITE_LABELS) else (extended.is_elliptic() or extended.is_parabolic())

    def extend(
        target: tuple[Rational, ...],
        projected_square: Rational,
        norm: Rational,
        choices: Sequence[tuple[tuple[CoxeterLabel, Rational], ...]],
        diagram: tuple[CoxeterLabel, ...],
        norm_position: int,
    ) -> None:
        position = len(target)
        if position == len(roots):
            if not diagram_allowed(diagram):
                return
            weights = inverse * vector(QQ, target)
            out.append(
                (
                    diagram,
                    norm_position,
                    NormedDynkinExtension(tuple(weights * coordinates), norm - projected_square, norm),
                )
            )
            return
        for label, pairing in choices[position]:
            square = projected_square + inverse[position, position] * pairing * pairing
            square += 2 * pairing * sum(inverse[position, j] * target[j] for j in range(position))
            # A positive-definite simple-root Gram matrix is a Stieltjes
            # matrix: its inverse is entrywise nonnegative. All chosen
            # pairings are nonpositive, so later terms cannot lower this
            # partial squared projection.
            if square <= norm:
                extend(
                    (*target, pairing),
                    square,
                    norm,
                    choices,
                    (*diagram, label),
                    norm_position,
                )

    for norm_position, new_norm in enumerate(QQ(norm) for norm in norms):
        choices = []
        for old_norm in root_norms:
            allowed = []
            for label in labels:
                square = _FOUR_COS_SQUARED[label] * old_norm * new_norm / 4
                if _norm_ratio_allows(label, old_norm, new_norm) and square.is_square():
                    allowed.append((label, -square.sqrt()))
            choices.append(tuple(allowed))

        extend((), QQ.zero(), new_norm, choices, (), norm_position)
    out.sort(key=lambda item: (item[0], item[1]))
    return tuple(item[2] for item in out)


def finite_type_extensions(
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root at angles pi/2, pi/3, pi/4, pi/6 with a nonzero orthogonal component."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, _FINITE_LABELS) if extension.residual_norm != 0)


def parabolic_extensions(
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root lying in the span of the roots (residual norm 0), at any crystallographic angle including parallel walls."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, (*_FINITE_LABELS, Infinity)) if extension.residual_norm == 0)
