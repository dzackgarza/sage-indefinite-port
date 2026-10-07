"""Exact normed-Dynkin extension candidates for the Allcock edgewalk."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product

from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ

type _RationalVector = FreeModuleElement[Rational]
type _RootRow = Sequence[int | Integer | Rational]

COXETER_INFINITY = QQ(-1)


@dataclass(frozen=True)
class NormedDynkinExtension:
    u_component: tuple[Rational, ...]
    residual_norm: Rational
    norm: Rational
    coxeter_entries: tuple[Rational, ...]


def _as_qq_matrix(gram: Matrix_integer_dense | Matrix_rational_dense) -> Matrix_rational_dense:
    result = matrix(QQ, gram)
    if not result.is_square() or result != result.transpose():
        raise ValueError("the Gram matrix must be square and symmetric")
    return result


def _as_qq_roots(roots: Sequence[_RootRow], dimension: int) -> tuple[_RationalVector, ...]:
    converted = tuple(vector(QQ, root) for root in roots)
    if any(len(root) != dimension for root in converted):
        raise ValueError("every root must have the Gram matrix dimension")
    return converted


def compute_coxeter_matrix(
    gram: Matrix_integer_dense | Matrix_rational_dense,
    roots: Sequence[_RootRow],
) -> tuple[Matrix_rational_dense, Matrix_rational_dense]:
    """Port ComputeCoxeterMatrix from pinned coxeter_dynkin.h."""
    gram_matrix = _as_qq_matrix(gram)
    root_vectors = _as_qq_roots(roots, gram_matrix.nrows())
    scalar = matrix(
        QQ,
        len(root_vectors),
        len(root_vectors),
        lambda i, j: (root_vectors[i] * gram_matrix * root_vectors[j].column())[0],
    )
    coxeter = matrix(QQ, scalar.nrows(), scalar.ncols())
    quotient_to_label = {
        QQ(1) / 4: QQ(3),
        QQ(1) / 2: QQ(4),
        QQ(3) / 4: QQ(6),
        QQ(1): COXETER_INFINITY,
    }
    for i in range(scalar.nrows()):
        if scalar[i, i] <= 0:
            raise ValueError("edgewalk roots must have positive norm")
        for j in range(scalar.ncols()):
            if i == j:
                coxeter[i, j] = scalar[i, j]
            elif scalar[i, j] == 0:
                coxeter[i, j] = 2
            else:
                quotient = scalar[i, j] ** 2 / (scalar[i, i] * scalar[j, j])
                if quotient not in quotient_to_label:
                    raise ValueError(f"root pair ({i}, {j}) has unsupported Coxeter quotient {quotient}")
                coxeter[i, j] = quotient_to_label[quotient]
    return coxeter, scalar


def _exact_negative_pairing(label: Rational, old_norm: Rational, new_norm: Rational) -> Rational | None:
    cs = {
        QQ(2): QQ(0),
        QQ(3): QQ(1) / 4,
        QQ(4): QQ(1) / 2,
        QQ(6): QQ(3) / 4,
        COXETER_INFINITY: QQ(1),
    }[label]
    square = cs * old_norm * new_norm
    return -square.sqrt() if square.is_square() else None


def _norm_ratio_allows(label: Rational, old_norm: Rational, new_norm: Rational) -> bool:
    if label == 3:
        return new_norm == old_norm
    if label == 4:
        return 2 * new_norm == old_norm or new_norm == 2 * old_norm
    if label == 6:
        return 3 * new_norm == old_norm or new_norm == 3 * old_norm
    return True


def compute_possible_extensions(
    gram: Matrix_integer_dense | Matrix_rational_dense,
    roots: Sequence[_RootRow],
    norms: Sequence[int | Integer | Rational],
    *,
    only_spherical: bool = False,
) -> tuple[NormedDynkinExtension, ...]:
    """Port ComputePossibleExtensions using the equivalent exact Gram test."""
    gram_matrix = _as_qq_matrix(gram)
    root_vectors = _as_qq_roots(roots, gram_matrix.nrows())
    coxeter, scalar = compute_coxeter_matrix(gram_matrix, roots)
    if scalar.nrows() and not scalar.is_positive_definite():
        raise ValueError("the existing root Gram matrix must be positive definite")
    inv = scalar.inverse() if scalar.nrows() else matrix(QQ, 0, 0)
    candidate_norms = tuple(QQ(x) for x in norms)
    if any(x <= 0 for x in candidate_norms):
        raise ValueError("candidate root norms must be positive")
    labels = (QQ(2), QQ(3), QQ(4), QQ(6)) + (() if only_spherical else (COXETER_INFINITY,))
    out: list[NormedDynkinExtension] = []
    seen: set[tuple[tuple[Rational, ...], Rational]] = set()
    for entries in product(labels, repeat=len(root_vectors)):
        for new_norm in candidate_norms:
            if any(not _norm_ratio_allows(entries[i], coxeter[i, i], new_norm) for i in range(len(entries))):
                continue
            exact = tuple(_exact_negative_pairing(entries[i], coxeter[i, i], new_norm) for i in range(len(entries)))
            pairings = tuple(x for x in exact if x is not None)
            if len(pairings) != len(exact):
                continue
            pv = vector(QQ, pairings)
            weights = inv * pv if len(entries) else vector(QQ, [])
            residual = new_norm - pv.dot_product(weights)
            if residual < 0 or (only_spherical and residual == 0):
                continue
            key = (entries, new_norm)
            if key in seen:
                continue
            seen.add(key)
            component = vector(QQ, gram_matrix.nrows())
            for i, w in enumerate(weights):
                component += w * root_vectors[i]
            out.append(NormedDynkinExtension(tuple(component), residual, new_norm, entries))
    return tuple(out)
