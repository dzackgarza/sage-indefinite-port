"""Normed Dynkin extensions of a root configuration, for the Allcock edgewalk (ComputePossibleExtensions)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product

from dzack_research.preamble.all import Lattices
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
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
    labels: tuple[CoxeterLabel, ...],
) -> tuple[NormedDynkinExtension, ...]:
    root_norms = tuple(QQ(int(root.q())) for root in roots)
    pairings = matrix(QQ, [[QQ(int(lattice.b(left, right))) for right in roots] for left in roots])
    if not pairings.is_positive_definite():
        raise ValueError(f"the roots {roots} span a form that is not positive definite: they are not the walls through a point of hyperbolic space")
    inverse = pairings.inverse()
    coordinates = matrix(QQ, [[QQ(int(entry)) for entry in root.to_tuple()] for root in roots])
    out: list[NormedDynkinExtension] = []
    for entries in product(labels, repeat=len(roots)):
        for new_norm in (QQ(norm) for norm in norms):
            if not all(_norm_ratio_allows(entries[i], root_norms[i], new_norm) for i in range(len(roots))):
                continue
            squares = tuple(_FOUR_COS_SQUARED[entries[i]] * root_norms[i] * new_norm / 4 for i in range(len(roots)))
            if not all(square.is_square() for square in squares):
                continue
            target = vector(QQ, [-square.sqrt() for square in squares])
            weights = inverse * target
            residual = new_norm - target.dot_product(weights)
            if residual < 0:
                continue
            out.append(NormedDynkinExtension(tuple(weights * coordinates), residual, new_norm))
    return tuple(out)


def finite_type_extensions(lattice: Lattices.ParentMethods, roots: Sequence[Lattices.ElementMethods], norms: Sequence[int]) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root at angles pi/2, pi/3, pi/4, pi/6 with a nonzero orthogonal component."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, _FINITE_LABELS) if extension.residual_norm != 0)


def parabolic_extensions(lattice: Lattices.ParentMethods, roots: Sequence[Lattices.ElementMethods], norms: Sequence[int]) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root lying in the span of the roots (residual norm 0), at any crystallographic angle including parallel walls."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, (*_FINITE_LABELS, Infinity)) if extension.residual_norm == 0)
