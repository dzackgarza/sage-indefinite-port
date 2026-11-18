r"""Stabilizers of markings in a degenerate ambient lattice (#41).

Split `L = R + C` with `R` the radical. An isometry is `(D, B, A)`, sending `r + c` to
`D r + B c + A c`, with `A` in `O(C)`, `D` in `GL(R)` and `B` in `Hom(C, R)`. A sublattice `S`
meeting `R` trivially is the graph of `\phi` over its projection `\bar S`. Then `g` preserves
`S` exactly when `A` preserves `\bar S` and `B = \phi A - D \phi` on `\bar S`. For
primitive `\bar S`, `B` extends to `C`, and any two extensions differ by `Hom(C / \bar S, R)`.
"""

from copy import copy

from sage.matrix.constructor import block_matrix, matrix
from sage.matrix.special import identity_matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage_extensions.categories.isotropic_parabolics import general_linear_images
from sage_extensions.categories.lattices import Lattices


def _full_projection_marking_generators(lattice, vectors, pointwise):
    r"""Exact stabilizer of a marking complementary to the rational radical.

    Let S be the integral span of the marked vectors and R=rad(L). Here
    S_Q + R_Q = L_Q and S_Q intersect R_Q = 0. The stabilizer of S in O(L)
    restricts rationally to O(S) x GL(R), with no shear kernel: S and R
    together span L_Q. Conversely each pair is an isometry of the rational
    form, and it preserves the integral lattice precisely when it belongs
    to the integral-structure stabilizer of L inside this rational group.

    The finite-index integral frame M=S+R is stable under the entire group,
    and L lies in (1/d)M for d=exponent(L/M). Thus the existing bounded
    integral-structure action gives complete Schreier generators of the
    integral intersection. No primitiveness of the projection S->L/R is
    assumed. Pointwise marking fixes S and removes its O(S) factor.
    """
    from sage_extensions.all import Modules
    from sage_extensions.categories.modules.pure.modules import ModulesOverIntegralDomains
    from sage_extensions.categories.rational_integral_stabilizers import IntegralStructureAction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _integral_generator_set

    sub = lattice.subobject_on(tuple(vectors))
    sub_inc = sub.inclusion()
    radical = lattice.radical()
    rad_inc = radical.inclusion()
    s = sub_inc.domain()
    r = rad_inc.domain()
    s_basis = tuple(sub_inc(point) for point in s.module_generators())
    r_basis = tuple(rad_inc(point) for point in r.module_generators())
    n = int(lattice.rank())
    if len(s_basis) + len(r_basis) != n:
        raise ValueError("the marked subspace must complement the rational radical")
    labels = tuple(lattice.module_generating_set())
    frame = matrix(QQ, [
        [int(point.to_vector()(label)) for label in labels]
        for point in (*s_basis, *r_basis)
    ])
    if frame.det() == 0:
        raise ValueError("the marked subspace must intersect the radical trivially")
    s_labels = tuple(s.module_generating_set())
    r_labels = tuple(r.module_generating_set())
    actions = []
    if not pointwise:
        for generator in EichlerRecursion().orthogonal_group_generators(s):
            block = identity_matrix(ZZ, n)
            for i, point in enumerate(s.module_generators()):
                coordinates = generator(point).to_vector()
                for j, label in enumerate(s_labels):
                    block[i, j] = int(coordinates(label))
            actions.append(block)
    for generator in r.general_linear_group().group_generators():
        block = identity_matrix(ZZ, n)
        for i, point in enumerate(r.module_generators()):
            coordinates = generator(point).to_vector()
            for j, label in enumerate(r_labels):
                block[len(s_basis) + i, len(s_basis) + j] = int(coordinates(label))
        actions.append(block)

    rational = lattice.base_change(lattice.base_ring().fraction_field_map())
    rational_labels = tuple(rational.module_generating_set())
    extensions = []
    for action in actions:
        rows = frame.inverse() * action * frame
        extensions.append(rational.Isom(rational)({
            label: rational.linear_combination({
                target: rational.base_ring()(rows[i, j])
                for j, target in enumerate(rational_labels) if rows[i, j]
            }) for i, label in enumerate(rational_labels)
        }))
    group = rational.O().subgroup(tuple(extensions) or (rational.O().identity(),))
    inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
    split = lattice.subobject_on((*s_basis, *r_basis))
    divisor = split.inclusion().cokernel().exponent()
    containing = Modules(lattice.base_ring()).Mono(split, inclusion.codomain())(
        lambda label: inclusion.codomain().wrap(
            rational.scalar_multiple(
                rational.base_ring().one() / rational.base_ring()(divisor),
                inclusion(split.inclusion()(split.module_generator(label))).underlying_element()
            )
        )
    )
    integral = IntegralStructureAction(group, inclusion).with_stable_containing_bound(containing)
    generators = _integral_generator_set(lattice, integral)
    basis = tuple(lattice.module_generator(label) for label in labels)
    return tuple(matrix(ZZ, [[int(generator(point).to_vector()(label)) for point in basis]
                             for label in labels]) for generator in generators)


def degenerate_marking_generators(lattice, vectors, pointwise):
    r"""Return generators of the stabilizer in `O(L)` of the span of ``vectors``, as matrices on the frame.

    With `S_0 = S \cap R`, the radical action is block triangular on `S_0` and a complement
    `R_1`. A pointwise stabilizer fixes `S_0` pointwise. A setwise one acts on it by `GL`, and
    `B` may also add `S_0` on `\bar S`.
    """
    frame = tuple(lattice.module_generators())
    size = len(frame)
    gram = matrix(ZZ, [[int(lattice.b(x, y)) for y in frame] for x in frame])
    base = gram.right_kernel().basis_matrix()
    rank_radical = base.nrows()
    assert 0 < rank_radical < size, f"{lattice} is not properly degenerate"
    _diagonal, _left, right = base.smith_form()
    complement = right.inverse()[rank_radical:]
    span = matrix(ZZ, [[int(c) for c in v] for v in vectors])
    trial = span.change_ring(QQ) * base.stack(complement).change_ring(QQ).inverse()
    kernel = trial[:, rank_radical:].left_kernel().basis_matrix()
    scale = kernel.denominator() if kernel.nrows() else 1
    inner = matrix(ZZ, scale * kernel).saturation() if kernel.nrows() else matrix(ZZ, 0, span.nrows())
    in_radical = matrix(ZZ, (inner.change_ring(QQ) * trial)[:, :rank_radical]) if inner.nrows() else matrix(ZZ, 0, rank_radical)
    if in_radical.nrows():
        in_radical = in_radical.saturation()
    k0 = in_radical.nrows()
    if k0:
        _d, _l, extra = in_radical.smith_form()
        radical_rows = in_radical.stack(extra.inverse()[k0:])
    else:
        radical_rows = identity_matrix(ZZ, rank_radical)
    radical = radical_rows * base
    change = radical.stack(complement)
    assert abs(change.det()) == 1
    inverse = change.change_ring(QQ).inverse()
    rows = span.change_ring(QQ) * inverse
    pivots = rows[:, rank_radical:].pivot_rows()
    lifts = rows[list(pivots), :]
    radical_part, projected = lifts[:, :rank_radical], lifts[:, rank_radical:]
    rho = projected.nrows()
    quotient_rank = size - rank_radical
    if (k0 == 0 and rho == quotient_rank and projected.denominator() == 1
            and any(d != 1 for d in matrix(ZZ, projected).elementary_divisors())):
        return _full_projection_marking_generators(lattice, vectors, pointwise)
    assert projected.denominator() == 1 and all(d == 1 for d in matrix(ZZ, projected).elementary_divisors()), (
        f"the projection of {vectors} to the nondegenerate quotient is not primitive"
    )
    quotient_gram = complement * gram * complement.transpose()
    quotient = Lattices(lattice.base_ring())([[int(c) for c in row] for row in quotient_gram.rows()])
    quotient_frame = tuple(quotient.module_generators())

    def in_quotient(row):
        return sum(
            (quotient.scalar_multiple(int(c), g) for c, g in zip(row, quotient_frame, strict=True)), quotient.zero()
        )

    group = quotient.O()
    if rho:
        embedding = quotient.subobject_on(tuple(in_quotient(r) for r in projected.rows())).inclusion()
        stabilizer = group.pointwise_stabilizer(embedding) if pointwise else group.setwise_stabilizer(embedding)
    else:
        stabilizer = group
    if rho:
        _d, _l, completion_right = matrix(ZZ, projected).smith_form()
        basis_rows = matrix(ZZ, projected).stack(completion_right.inverse()[rho:])
    else:
        basis_rows = identity_matrix(ZZ, quotient_rank)
    assert abs(basis_rows.det()) == 1
    basis_inverse = basis_rows.change_ring(QQ).inverse()
    phi = radical_part
    identity_radical = identity_matrix(ZZ, rank_radical)
    identity_quotient = identity_matrix(ZZ, quotient_rank)

    def assemble(action, unipotent, shear_on_basis):
        shear = basis_inverse * matrix(QQ, quotient_rank, rank_radical, shear_on_basis)
        assert shear.denominator() == 1
        block = block_matrix(
            ZZ, [[unipotent, matrix(ZZ, shear).transpose()], [matrix(ZZ, quotient_rank, rank_radical), action]], subdivide=False
        )
        columns = change.transpose()
        return matrix(ZZ, columns * block.change_ring(QQ) * columns.inverse())

    def padded(values):
        return [list(row) for row in values.rows()] + [[0] * rank_radical for _ in range(quotient_rank - rho)]

    results = []
    for generator in stabilizer.group_generators():
        action = matrix(ZZ, [[int(c) for c in generator(v)] for v in quotient_frame]).transpose()
        shear_values = phi * 0
        if rho:
            moved = matrix(QQ, projected) * action.transpose()
            shear_values = matrix(QQ, projected).solve_left(moved) * phi - phi
        results.append(assemble(action, identity_radical, padded(shear_values) if rho else [[0] * rank_radical] * quotient_rank))

    def radical_generator(unipotent):
        shear_values = phi - phi * unipotent.transpose()
        results.append(assemble(identity_quotient, unipotent, padded(shear_values) if rho else [[0] * rank_radical] * quotient_rank))

    def diagonal_block(images, start, length):
        block = identity_matrix(ZZ, rank_radical)
        for i, row in enumerate(images):
            for j, value in enumerate(row):
                block[start + j, start + i] = value
        return block

    tail = rank_radical - k0
    if tail:
        for images in general_linear_images(tuple(identity_matrix(ZZ, tail).rows())):
            radical_generator(diagonal_block(images, k0, tail))
    if k0 and not pointwise:
        for images in general_linear_images(tuple(identity_matrix(ZZ, k0).rows())):
            radical_generator(diagonal_block(images, 0, k0))
    for i in range(k0):
        for j in range(k0, rank_radical):
            shear = copy(identity_radical)
            shear[i, j] = 1
            radical_generator(shear)
    for row in range(rho, quotient_rank):
        for column in range(rank_radical):
            values = [[0] * rank_radical for _ in range(quotient_rank)]
            values[row][column] = 1
            results.append(assemble(identity_quotient, identity_radical, values))
    if not pointwise:
        for row in range(rho):
            for column in range(k0):
                values = [[0] * rank_radical for _ in range(quotient_rank)]
                values[row][column] = 1
                results.append(assemble(identity_quotient, identity_radical, values))
    return tuple(results)
