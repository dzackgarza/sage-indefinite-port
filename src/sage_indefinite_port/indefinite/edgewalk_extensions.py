"""Normed Dynkin extensions of a root configuration, for the Allcock edgewalk (ComputePossibleExtensions)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations

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


def _diagram_extension_labels(coxeter_rows: list[list[Integer]], spherical: bool) -> tuple[tuple[CoxeterLabel, ...], ...]:
    r"""Admissible diagram extensions (coxeter_dynkin.h:715–1671).

    The reference constructs at most four noncommuting bonds using finite
    and affine Dynkin types. The Coxeter category certifies each resulting
    diagram before the root-norm test.
    """
    n = len(coxeter_rows)
    base: tuple[CoxeterLabel, ...] = (Integer(2),) * n
    candidates: list[tuple[CoxeterLabel, ...]] = [base]
    isolated = tuple(i for i in range(n) if all(coxeter_rows[i][j] == 2 for j in range(n) if i != j))

    def add(entries: tuple[tuple[int, CoxeterLabel], ...]) -> None:
        labels = list(base)
        for index, value in entries:
            labels[index] = value
        candidates.append(tuple(labels))

    # The pinned FindDiagramExtensions_Efficient emits the disconnected A1,
    # then the single-edge Dynkin enlargements, and only subsequently its
    # disconnected A1/G2/parallel cases. Keep diagram-first output ordering.
    for label in (Integer(3), Integer(4), Integer(6)):
        for i in range(n):
            add(((i, label),))
    if not spherical:
        for i in isolated:
            add(((i, Infinity),))

    # Source: finite A/B/D/E/F and affine A/B/C/D/E/F/G extensions.
    for i, j in combinations(range(n), 2):
        add(((i, Integer(3)), (j, Integer(3))))
        add(((i, Integer(3)), (j, Integer(4))))
        add(((i, Integer(4)), (j, Integer(3))))
        if not spherical and i in isolated and j in isolated:
            add(((i, Integer(4)), (j, Integer(4))))
            add(((i, Integer(3)), (j, Integer(6))))
            add(((i, Integer(6)), (j, Integer(3))))
    # Source: D4 and the higher D/E triples, plus affine triples.
    for subset in combinations(range(n), 3):
        add(tuple((i, Integer(3)) for i in subset))
        if not spherical and all(i in isolated for i in subset):
            for singled in subset:
                add(tuple((i, Integer(4) if i == singled else Integer(3)) for i in subset))
    if not spherical:
        for subset in combinations(isolated, 4):
            add(tuple((i, Integer(3)) for i in subset))

    rows = [list(row) + [Integer(2)] for row in coxeter_rows]
    rows.append([Integer(2)] * n + [Integer(1)])
    accepted: list[tuple[CoxeterLabel, ...]] = []
    seen: set[tuple[CoxeterLabel, ...]] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        for i, label in enumerate(candidate):
            rows[i][-1] = rows[-1][i] = Integer(-1) if label == Infinity else Integer(label)
        diagram = CoxeterDiagrams()(rows)
        if diagram.is_elliptic() or (not spherical and diagram.is_parabolic()):
            accepted.append(candidate)
    return tuple(accepted)


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
    # extensions before prescribed norms. FindDiagramExtensions_Efficient
    # (715-1671) lists them in Dynkin construction order.
    out: list[tuple[tuple[CoxeterLabel, ...], int, NormedDynkinExtension]] = []
    coxeter_rows = [[Integer(1) if i == j else Integer(2) for j in range(len(roots) + 1)] for i in range(len(roots) + 1)]
    for i in range(len(roots)):
        for j in range(i):
            old_pair = pairings[i, j]
            squares = 4 * old_pair * old_pair / (root_norms[i] * root_norms[j])
            label = next((edge for edge, value in _FOUR_COS_SQUARED.items() if squares == value), None)
            if label is None:
                raise ValueError("the supplied roots have noncrystallographic Coxeter pairings")
            coxeter_rows[i][j] = coxeter_rows[j][i] = Integer(-1) if label == Infinity else label

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

    for diagram in _diagram_extension_labels(
        [row[:-1] for row in coxeter_rows[:-1]],
        spherical=len(labels) == len(_FINITE_LABELS),
    ):
        for norm_position, new_norm in enumerate(QQ(norm) for norm in norms):
            choices: list[tuple[tuple[CoxeterLabel, Rational], ...]] = []
            for old_norm, label in zip(root_norms, diagram, strict=True):
                square = _FOUR_COS_SQUARED[label] * old_norm * new_norm / 4
                allowed = ((label, -square.sqrt()),) if _norm_ratio_allows(label, old_norm, new_norm) and square.is_square() else ()
                choices.append(allowed)
            if all(choices):
                extend((), QQ.zero(), new_norm, choices, (), norm_position)
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
