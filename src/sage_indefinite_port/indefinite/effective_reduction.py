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
    if bound < 0:
        raise ValueError("the integral coordinate bound must be nonnegative")
    if marking.rank() != target.rank():
        return None
    if marking.transpose().elementary_divisors() != target.transpose().elementary_divisors():
        # A unimodular ambient transformation preserves the Smith
        # invariants of the marked column module, including its index in
        # the primitive saturation. Equal Gram matrices alone do not
        # detect this obstruction for isotropic markings.
        return None
    relations = marking.change_ring(QQ).right_kernel_matrix()
    if target.change_ring(QQ) * relations.transpose() != 0:
        # An invertible integral transporter preserves every rational
        # linear relation among the marked columns. Reject impossible
        # dependent markings before enumerating the coefficient box.
        return None
    if marking.nrows() == n and marking.ncols() == n and marking.det() != 0:
        # A full-rank marking determines the transporter uniquely over Q:
        # g X = Y implies g = Y X^-1. This applies equally to a
        # nonprimitive full-rank marking; its index is irrelevant.
        candidate = target.change_ring(QQ) * marking.change_ring(QQ).inverse()
        if (candidate.denominator() == 1
                and _entry_norm(candidate) <= bound):
            g = candidate.change_ring(ZZ)
            if (abs(g.det()) == 1 and g.transpose() * gram * g == gram
                    and all(g * F == F * g for F in actions)):
                return g
        return None
    if n == 2 and marking.rank():
        # A nonzero binary marking has at most two integral transporters
        # when it contains a nonisotropic vector; for an isotropic vector
        # the split binary orthogonal group has at most one transporter.
        # Check this complete candidate family against the coordinate
        # bound and the prescribed centralizer, avoiding an enormous box.
        for j in range(marking.ncols()):
            v, w = marking.matrix_from_columns([j]), target.matrix_from_columns([j])
            if (v.transpose()*gram*v)[0, 0] == 0:
                continue
            bv, bw = gram*v, gram*w
            source_frame = matrix(ZZ, [[v[0, 0], bv[1, 0]],
                                       [v[1, 0], -bv[0, 0]]])
            for sign in (1, -1):
                target_frame = matrix(ZZ, [[w[0, 0], sign*bw[1, 0]],
                                           [w[1, 0], -sign*bw[0, 0]]])
                candidate = orthogonal_transporter(gram, source_frame, target_frame)
                if (candidate is not None and _entry_norm(candidate) <= bound
                        and candidate*marking == target
                        and all(candidate*F == F*candidate for F in actions)):
                    return candidate
            return None
        candidate = orthogonal_transporter(gram, marking, target)
        if (candidate is not None and _entry_norm(candidate) <= bound
                and all(candidate*F == F*candidate for F in actions)):
            return candidate
        return None

    def candidates(column):
        # The effective Li--Margulis bound can exceed machine-sized range
        # lengths. Enumerate each finite coordinate box with an unbounded
        # integer odometer, without materializing its exponentially many
        # vectors or a dictionary of norm shells.
        coordinates = [-bound] * n
        while True:
            c = matrix(ZZ, n, 1, coordinates)
            if (c.transpose() * gram * c)[0, 0] == gram[column, column]:
                yield c
            for position in range(n - 1, -1, -1):
                if coordinates[position] < bound:
                    coordinates[position] += 1
                    coordinates[position + 1:] = [-bound] * (n - position - 1)
                    break
            else:
                return

    def extend(columns):
        j = len(columns)
        if j == n:
            g = block_matrix(ZZ, [columns], subdivide=False)
            if g * marking == target and all(g * F == F * g for F in actions):
                return g
            return None
        for c in candidates(j):
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
        # For an n-by-k column marking with k<n, Sage includes n-k zero
        # elementary divisors. The k nonzero Smith invariants are those
        # of its transpose; zeros are not failures of primitivity.
        assert all(e == 1 for e in m.transpose().elementary_divisors()), f"{m} is not primitive"
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
    if (gram.nrows() != gram.ncols() or marking.nrows() != gram.nrows()
            or target.dimensions() != marking.dimensions()):
        raise ValueError("marked transporter requires matching ambient and marking dimensions")
    if not marking.ncols():
        return identity_matrix(ZZ, gram.nrows())
    if target.transpose() * gram * target != marking.transpose() * gram * marking:
        return None
    if gram.nrows() <= 1:
        # The zero-dimensional group is trivial; for a nonzero rank-one
        # form O([d],Z)={+1,-1}. Enumerate these exact integral isometries
        # instead of applying Li--Margulis, whose hypotheses require n>=3.
        identity = identity_matrix(ZZ, gram.nrows())
        for candidate in (identity, -identity):
            if candidate * marking == target:
                return candidate
        return None
    if marking.rank() != target.rank():
        # An invertible isometry preserves the rank of every marked frame,
        # even when its restricted Gram form is degenerate.
        return None
    if gram.rank() == 0:
        # For the zero quadratic form, O(0_n,Z)=GL_n(Z). Left integral
        # row operations preserve the row Hermite normal form of the
        # marked matrix, and equality of those forms is sufficient:
        # U_X X=H=U_Y Y implies (U_Y^-1 U_X) X=Y. Both U's are
        # unimodular, giving a complete decision and integral witness
        # for arbitrary rank, nonprimitive and dependent markings.
        source_hnf, source_change = marking.hermite_form(transformation=True)
        target_hnf, target_change = target.hermite_form(transformation=True)
        if source_hnf != target_hnf:
            return None
        return target_change.inverse()*source_change
    if marking.transpose().elementary_divisors() != target.transpose().elementary_divisors():
        # Smith invariants of the marking are unchanged by left
        # multiplication with GL_n(Z). In particular a primitive null
        # vector cannot be transported to a nonprimitive null vector.
        return None
    if marking.rank() < marking.ncols():
        # Dependent markings have a unique rational relation matrix on
        # their independent columns. Every linear isometry preserves those
        # relations, so discard dependent columns only after verifying the
        # target satisfies exactly the same relations.
        pivots = tuple(marking.pivots())
        if not pivots:
            return identity_matrix(ZZ, gram.nrows()) if not target else None
        independent = marking.matrix_from_columns(pivots)
        moved = target.matrix_from_columns(pivots)
        relations = independent.change_ring(QQ).solve_right(marking.change_ring(QQ))
        if moved.change_ring(QQ) * relations != target.change_ring(QQ):
            return None
        return orthogonal_transporter(gram, independent, moved)
    if marking.rank() == marking.ncols() and marking.ncols() < gram.nrows() and any(
        abs(d) != 1 for d in marking.transpose().elementary_divisors()
    ):
        # An isometry maps the primitive saturation of a nonprimitive
        # marking to the saturation of its image. Extend the prescribed
        # map rationally to that saturation, then demand integrality and
        # primitivity before applying the primitive marked decision.
        saturated = marking.transpose().saturation().transpose()
        coefficients = marking.change_ring(QQ).solve_right(saturated.change_ring(QQ))
        moved = target.change_ring(QQ) * coefficients
        if moved.denominator() != 1:
            return None
        moved = moved.change_ring(ZZ)
        if any(abs(d) != 1 for d in moved.transpose().elementary_divisors()):
            return None
        return orthogonal_transporter(gram, saturated, moved)
    if marking.nrows() == gram.nrows() and marking.ncols() == gram.nrows() and marking.det() != 0:
        # No Li--Margulis complement is needed for a full-rank marking:
        # the rational extension is unique and its integrality is exact.
        candidate = target.change_ring(QQ) * marking.change_ring(QQ).inverse()
        if candidate.denominator() != 1:
            return None
        g = candidate.change_ring(ZZ)
        return g if abs(g.det()) == 1 and g.transpose() * gram * g == gram else None
    if gram.is_positive_definite() or (-gram).is_positive_definite():
        # For a definite integral form the entire orthogonal group is
        # finite. Enumerate its exact owned integral automorphisms;
        # neither an indefinite complement nor the Li--Margulis bound
        # applies to this case. This covers every marking rank and
        # nonprimitive content without imposing primitive hypotheses.
        from sage_extensions.all import Lattices, ZZ as OwnedZZ

        lattice = Lattices(OwnedZZ)([[int(c) for c in row] for row in gram.rows()])
        frame = tuple(lattice.module_generators())
        for isometry in lattice.O():
            candidate = matrix(ZZ, [[int(c) for c in isometry(v)] for v in frame]).transpose()
            if candidate*marking == target:
                return candidate
        return None
    if (gram.det() != 0 and marking.ncols() == gram.nrows()-1
            and (marking.transpose()*gram*marking).det() != 0):
        # A nondegenerate codimension-one marking has a rational
        # orthogonal complement of dimension one. An extension is
        # determined by the marked images and the sign of the image of
        # a normal vector. If u,v are nonzero integral normals, that
        # image is +/-sqrt(q(u)/q(v))*v. Test both unique rational
        # candidates for integral unimodularity and the full Gram
        # identity. This is complete in every signature and Witt index.
        normal_x = (marking.transpose()*gram).right_kernel_matrix().transpose()
        normal_y = (target.transpose()*gram).right_kernel_matrix().transpose()
        norm_x = (normal_x.transpose()*gram*normal_x)[0, 0]
        norm_y = (normal_y.transpose()*gram*normal_y)[0, 0]
        ratio = QQ(norm_x)/norm_y
        if ratio < 0 or not ratio.is_square():
            return None
        scale = ratio.sqrt()
        source_frame = marking.augment(normal_x).change_ring(QQ)
        for sign in (1, -1):
            target_frame = target.augment(normal_y).change_ring(QQ)
            target_frame.set_column(target_frame.ncols()-1, sign*scale*normal_y.column(0))
            candidate = target_frame*source_frame.inverse()
            if candidate.denominator() == 1:
                candidate = candidate.change_ring(ZZ)
                if abs(candidate.det()) == 1 and candidate.transpose()*gram*candidate == gram:
                    return candidate
        return None
    if (gram.det() != 0 and 0 < marking.ncols() <= gram.nrows()-2
            and (marking.transpose()*gram*marking).det() != 0):
        # A nondegenerate marking splits the rational space into its
        # marked span and an integral orthogonal complement of rank >=2.
        # If that complement is definite, its entire integral isometry
        # torsor is finite. Enumerate the exact complement isometries
        # and check the induced rational ambient maps for integrality;
        # the latter is precisely the finite-index gluing obstruction.
        basis_x, gram_x = _complement(gram, marking)
        basis_y, gram_y = _complement(gram, target)
        if gram_x.is_positive_definite() or (-gram_x).is_positive_definite():
            from sage_extensions.all import Lattices, ZZ as OwnedZZ

            category = Lattices(OwnedZZ)
            Kx = category([[int(c) for c in row] for row in gram_x.rows()])
            Ky = category([[int(c) for c in row] for row in gram_y.rows()])
            torsor = Kx.Isom(Ky)
            if torsor.is_empty():
                return None
            source_frame = marking.augment(basis_x).change_ring(QQ)
            for h in torsor:
                H = matrix(ZZ, [[int(c) for c in h(v)] for v in Kx.module_generators()]).transpose()
                target_frame = target.augment(basis_y*H).change_ring(QQ)
                candidate = target_frame*source_frame.inverse()
                if candidate.denominator() == 1:
                    candidate = candidate.change_ring(ZZ)
                    if abs(candidate.det()) == 1 and candidate.transpose()*gram*candidate == gram:
                        return candidate
            return None
    if gram.nrows() == 2 and gram.rank() == 1:
        # A primitive radical vector r extends to a unimodular basis
        # (u,r), in which B=diag(a,0), a != 0. Its entire integral
        # isometry group consists of [[eps,0],[b,delta]], where
        # eps,delta in {+1,-1}, b in ZZ. Solve the marked equations
        # exactly for b; no bounded search or nondegenerate hypothesis.
        from sage.arith.misc import xgcd

        radical = gram.right_kernel_matrix().row(0)
        d, s, t = xgcd(ZZ(radical[0]), ZZ(radical[1]))
        radical = radical / d
        frame = matrix(ZZ, [[t, radical[0]], [-s, radical[1]]])
        assert frame.det() == 1
        inverse = frame.inverse()
        x, y = inverse*marking, inverse*target
        for epsilon in (1, -1):
            if any(y[0, j] != epsilon*x[0, j] for j in range(x.ncols())):
                continue
            for delta in (1, -1):
                first = next((j for j in range(x.ncols()) if x[0, j]), None)
                if first is None:
                    b = ZZ(0)
                else:
                    quotient = QQ(y[1, first]-delta*x[1, first]) / x[0, first]
                    if quotient.denominator() != 1:
                        continue
                    b = ZZ(quotient)
                candidate = frame*matrix(ZZ, [[epsilon, 0], [b, delta]])*inverse
                if candidate*marking == target:
                    return candidate
        return None
    if gram.nrows() == 2 and gram.det() < 0 and ZZ(gram[0, 1]**2 - gram[0, 0]*gram[1, 1]).is_square():
        # Every rationally split indefinite binary form has precisely two
        # rational isotropic lines. A unimodular integral isometry permutes
        # these lines and sends their primitive integral generators to
        # primitive generators with signs. Enumerating the eight signed
        # permutations and checking integrality gives the *entire* O(B),
        # even when the two rays do not form an integral basis.
        from math import gcd
        a, b, c = (ZZ(gram[0, 0]), ZZ(gram[0, 1]), ZZ(gram[1, 1]))
        root = ZZ(b*b-a*c).sqrt()
        rays = ((-b+root, a), (-b-root, a)) if a else ((1, 0), (-c, 2*b))
        primitive = tuple((ZZ(x)//g, ZZ(y)//g)
                          for x, y in rays for g in (gcd(int(x), int(y)),))
        frame = matrix(ZZ, primitive).transpose().change_ring(QQ)
        if frame.det() == 0:
            raise ArithmeticError("the two split isotropic rays must be independent")
        for permutation in ((0, 1), (1, 0)):
            for signs in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                images = matrix(ZZ, [[signs[i]*primitive[permutation[i]][j]
                                      for j in range(2)] for i in range(2)]).transpose()
                candidate = images.change_ring(QQ)*frame.inverse()
                if candidate.denominator() != 1:
                    continue
                candidate = candidate.change_ring(ZZ)
                if (abs(candidate.det()) == 1 and candidate.transpose()*gram*candidate == gram
                        and candidate*marking == target):
                    return candidate
        return None
    if gram.nrows() == 2 and gram.det() != 0:
        # For a nonisotropic v, its perpendicular line is generated over Q
        # by t=((Bv)_2,-(Bv)_1). Its square is det(B)*q(v). Therefore
        # exactly two rational isometries can send v to w, obtained by
        # sending t to +/-t_w. Testing integrality and the full marking
        # decides the Pell-type binary transporter without an infinite
        # search or the rank >= 3 Li--Margulis theorem.
        for j in range(marking.ncols()):
            v = marking.matrix_from_columns([j])
            w = target.matrix_from_columns([j])
            if (v.transpose()*gram*v)[0, 0] == 0:
                continue
            bv, bw = gram*v, gram*w
            frame = matrix(ZZ, [[v[0, 0], bv[1, 0]], [v[1, 0], -bv[0, 0]]]).change_ring(QQ)
            for sign in (1, -1):
                image_frame = matrix(ZZ, [[w[0, 0], sign*bw[1, 0]],
                                          [w[1, 0], -sign*bw[0, 0]]]).change_ring(QQ)
                candidate = image_frame*frame.inverse()
                if candidate.denominator() == 1:
                    candidate = candidate.change_ring(ZZ)
                    if (abs(candidate.det()) == 1 and candidate.transpose()*gram*candidate == gram
                            and candidate*marking == target):
                        return candidate
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
