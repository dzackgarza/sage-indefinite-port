r"""Explicit bounds that turn bounded transporter searches into decisions (#53).

Two tools are proved here, for an integral nondegenerate Gram matrix `B` on `L = \mathbb Z^n`, isometries
`F_1, \dots, F_s` and the centralizer `C = C_{O(L)}(F)`.

Image invariant. For a marking `X : \mathbb Z^k \to L` put `M_X = [X^t B W]_{W \in \{1, F_1, \dots, F_s\}}`,
a `k(s+1) \times n` matrix, and let `I(X) \subseteq \mathbb Z^{k(s+1)}` be its column lattice. For `g \in C`
and `Y = g X`, `g^t B = B g^{-1}` and `g W = W g` give `Y^t B W = X^t B W g^{-1}`, so `M_Y = M_X g^{-1}` and
`I(Y) = I(X)`, since `g^{-1} \in GL_n(\mathbb Z)`. Unequal lattices prove that no transporter exists. For
`s = 0` and a vector `X = v` this is `\mathrm{div}(v) \mathbb Z`.

Transporter box for `C = O(L)`. Let `X` be primitive with `A = X^t B X` nonsingular, `K_X = X^\perp` and
`K_Y = Y^\perp` with Gram matrices `G_X, G_Y` in integral bases, indefinite of rank `r \ge 3`. Every
`g \in O(L)` with `g X = Y` is `g = (Y | K_Y) \,\mathrm{diag}(1, h)\, (X | K_X)^{-1}` for an isometry
`h : G_X \to G_Y`, and the admissible `h` form one coset `h_1 S` of the stabilizer `S \le O(G_X)` of the
glue point (the projections of the standard basis to `K_X^\vee / K_X`). `S` contains the kernel of
`O(G_X) \to \mathrm{Aut}(K_X^\vee / K_X)`, so `[O(G_X) : S] \le |\mathrm{Aut}(D)| \le |D|^d` for
`D = K_X^\vee / K_X` with `d` invariant factors (an automorphism is fixed by the images of `d` generators).
Li-Margulis, *Effective estimates on integral quadratic forms*, GAFA 26 (2016), Theorem 1, eq. (3), gives
`h_0 : G_X \to G_Y` with entries below `N_1` (:func:`li_margulis.isometry_box_bound`), and Theorem 2 gives
generators of `O(G_X)` with entries below `N_2` (:func:`li_margulis.generator_bound`); the constants are
:func:`li_margulis.default_constants` of rank `r`. A spanning tree of the Schreier graph of `O(G_X)` on
`O(G_X)/S` has depth below the index, so `h_0^{-1} h_1 S = w S` for a word `w` of length `\ell < |D|^d`
(Holt-Eick-O'Brien, *Handbook of Computational Group Theory*, Section 4.1). With
`\|PQ\|_\infty \le r \|P\|_\infty \|Q\|_\infty`, `h = h_0 w` has entries at most `r^\ell N_1 N_2^\ell`, and
`g` has entries at most `n^2 \|(Y|K_Y)\|_\infty \max(1, \|h\|_\infty) \|(X|K_X)^{-1}\|_\infty`.
"""

from itertools import product

from sage.matrix.constructor import block_matrix, matrix
from sage.matrix.special import identity_matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.quadratic_forms.quadratic_form import QuadraticForm

from sage_indefinite_port.indefinite.li_margulis import (
    default_constants,
    generator_bound,
    isometry_box_bound,
    theorem_one_constant,
    theorem_two_constant,
)


def _entry_norm(m):
    return max((abs(c) for c in m.list()), default=0)


def image_invariant(gram, marking, actions):
    r"""Return the Hermite normal form of the column lattice `I(X)` of `M_X`, a `C(F)`-orbit invariant."""
    gram, marking = matrix(ZZ, gram), matrix(ZZ, marking)
    words = [identity_matrix(ZZ, gram.nrows()), *(matrix(ZZ, F) for F in actions)]
    stacked = block_matrix(ZZ, [[marking.transpose() * gram * W] for W in words], subdivide=False)
    rows = stacked.transpose().echelon_form()
    return rows.matrix_from_rows(range(rows.rank()))


def separated_by_invariant(gram, marking, target, actions):
    r"""Return ``True`` when `I(X) \ne I(Y)`, which proves that no `g \in C(F)` has `g X = Y`."""
    return image_invariant(gram, marking, actions) != image_invariant(gram, target, actions)


def centralizer_transporter_in_box(gram, actions, marking, target, bound):
    r"""Return `g \in C(F)` with `g X = Y` and entries at most ``bound``, or ``None`` when no such `g` has them.

    The columns `g e_j` are enumerated among the vectors of the box with `c^t B c = B_{jj}`, with the cross
    conditions `c_i^t B c_j = B_{ij}` checked as each column is placed; a full match is an isometry of the
    nondegenerate `B`, hence in `GL_n(\mathbb Z)`.
    """
    gram, marking, target = matrix(ZZ, gram), matrix(ZZ, marking), matrix(ZZ, target)
    actions = [matrix(ZZ, F) for F in actions]
    n = gram.nrows()
    assert gram.det() != 0, "the form must be nondegenerate"
    bound = int(bound)
    span = range(-bound, bound + 1)
    box = [matrix(ZZ, n, 1, list(c)) for c in product(span, repeat=n)]
    norms = {}
    for c in box:
        norms.setdefault(int((c.transpose() * gram * c)[0, 0]), []).append(c)
    candidates = [norms.get(int(gram[j, j]), []) for j in range(n)]

    def extend(columns):
        j = len(columns)
        if j == n:
            g = block_matrix(ZZ, [columns], subdivide=False)
            if g * marking == target and all(g * F == F * g for F in actions):
                return g
            return None
        for c in candidates[j]:
            if all((columns[i].transpose() * gram * c)[0, 0] == gram[i, j] for i in range(j)):
                found = extend([*columns, c])
                if found is not None:
                    return found
        return None

    return extend([])


def _complement(gram, marking):
    basis = (marking.transpose() * gram).right_kernel_matrix().transpose()
    return basis, basis.transpose() * gram * basis


def orthogonal_transporter_bound(gram, marking, target):
    r"""Return an integer `N` such that, if some `g \in O(L)` has `g X = Y`, one has entries at most `N`.

    The hypotheses of the module docstring are asserted: `X` and `Y` primitive with the same nonsingular
    Gram matrix and an indefinite orthogonal complement of rank at least three.
    """
    gram, marking, target = matrix(ZZ, gram), matrix(ZZ, marking), matrix(ZZ, target)
    n = gram.nrows()
    source_gram = marking.transpose() * gram * marking
    assert source_gram.det() != 0, "the source Gram matrix must be nonsingular"
    assert target.transpose() * gram * target == source_gram, "the target has another Gram matrix"
    for m in (marking, target):
        assert all(e == 1 for e in m.elementary_divisors()), f"{m} is not primitive"
    basis_x, gram_x = _complement(gram, marking)
    basis_y, gram_y = _complement(gram, target)
    r = gram_x.nrows()
    assert r >= 3, "Li-Margulis Theorems 1 and 2 need a complement of rank at least three"
    positive, negative, _zero = QuadraticForm(QQ, 2 * gram_x.change_ring(QQ)).signature_vector()
    assert positive and negative, "Li-Margulis section 6 needs an indefinite complement"
    constants = default_constants(r)
    n1 = isometry_box_bound(gram_x, gram_y, theorem_one_constant(constants))
    n2 = generator_bound(gram_x, theorem_two_constant(constants))
    divisors = [e for e in gram_x.elementary_divisors() if abs(e) > 1]
    order = 1
    for e in divisors:
        order *= abs(e)
    length = order ** len(divisors) - 1
    h = r**length * n1 * n2**length
    frame_x = block_matrix(QQ, [[marking, basis_x]], subdivide=False)
    frame_y = block_matrix(QQ, [[target, basis_y]], subdivide=False)
    return int((n**2 * _entry_norm(frame_y) * max(1, h) * _entry_norm(frame_x.inverse())).ceil())


def orthogonal_transporter(gram, marking, target):
    r"""Return `g \in O(L)` with `g X = Y`, or ``None``, which is then proved.

    ``None`` comes from :func:`separated_by_invariant` or from the exhausted box of
    :func:`orthogonal_transporter_bound`; boxes `1, 2, 4, \dots` are searched up to that bound.
    """
    gram, marking, target = matrix(ZZ, gram), matrix(ZZ, marking), matrix(ZZ, target)
    if target.transpose() * gram * target != marking.transpose() * gram * marking:
        return None
    if separated_by_invariant(gram, marking, target, []):
        return None
    limit = orthogonal_transporter_bound(gram, marking, target)
    size = 1
    while True:
        found = centralizer_transporter_in_box(gram, [], marking, target, min(size, limit))
        if found is not None or size >= limit:
            return found
        size *= 2
