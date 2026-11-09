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
    assert projected.denominator() == 1 and all(d == 1 for d in matrix(ZZ, projected).elementary_divisors()), (
        f"the projection of {vectors} to the nondegenerate quotient is not primitive"
    )
    quotient_rank = size - rank_radical
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
