r"""Stabilizers of markings in a degenerate ambient lattice (#41).

Split `L = R + C` with `R` the radical. An isometry is `(D, B, A)`, sending `r + c` to
`D r + B c + A c`, with `A` in `O(C)`, `D` in `GL(R)` and `B` in `Hom(C, R)`. A sublattice `S`
meeting `R` trivially is the graph of `\phi` over its projection `\bar S`. Then `g` preserves
`S` exactly when `A` preserves `\bar S` and `B = \phi A - D \phi` on `\bar S`. For
primitive `\bar S`, `B` extends to `C`, and any two extensions differ by `Hom(C / \bar S, R)`.
"""

from sage.matrix.constructor import block_matrix, matrix
from sage.matrix.special import identity_matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage_extensions.categories.isotropic_parabolics import general_linear_images
from sage_extensions.categories.lattices import Lattices


def degenerate_marking_generators(lattice, vectors, pointwise):
    r"""Return generators of the stabilizer in `O(L)` of the span of ``vectors``, as matrices on the frame."""
    frame = tuple(lattice.module_generators())
    size = len(frame)
    gram = matrix(ZZ, [[int(lattice.b(x, y)) for y in frame] for x in frame])
    radical = gram.right_kernel().basis_matrix()
    rank_radical = radical.nrows()
    assert 0 < rank_radical < size, f"{lattice} is not properly degenerate"
    _diagonal, _left, right = radical.smith_form()
    complement = right.inverse()[rank_radical:]
    change = radical.stack(complement)
    assert abs(change.det()) == 1
    inverse = change.change_ring(QQ).inverse()
    rows = matrix(QQ, [[int(c) for c in v] for v in vectors]) * inverse
    radical_part, projected = rows[:, :rank_radical], rows[:, rank_radical:]
    assert projected.rank() == len(vectors), f"{vectors} meets the radical of {lattice}"
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

    embedding = quotient.subobject_on(tuple(in_quotient(r) for r in projected.rows())).inclusion()
    group = quotient.O()
    stabilizer = group.pointwise_stabilizer(embedding) if pointwise else group.setwise_stabilizer(embedding)
    _d, _l, completion_right = matrix(ZZ, projected).smith_form()
    basis_rows = matrix(ZZ, projected).stack(completion_right.inverse()[projected.nrows() :])
    assert abs(basis_rows.det()) == 1
    basis_inverse = basis_rows.change_ring(QQ).inverse()
    phi = radical_part
    identity_radical = identity_matrix(ZZ, rank_radical)

    def assemble(action, unipotent, shear_on_basis):
        shear = basis_inverse * matrix(QQ, shear_on_basis)
        assert shear.denominator() == 1
        block = block_matrix(
            ZZ, [[unipotent, matrix(ZZ, shear).transpose()], [matrix(ZZ, quotient_rank, rank_radical), action]], subdivide=False
        )
        columns = change.transpose()
        return matrix(ZZ, columns * block.change_ring(QQ) * columns.inverse())

    results = []
    zero_shear = lambda: [[0] * rank_radical for _ in range(quotient_rank)]
    for generator in stabilizer.group_generators():
        action = matrix(ZZ, [[int(c) for c in generator(v)] for v in quotient_frame]).transpose()
        moved = matrix(QQ, projected) * action.transpose()
        coordinates = matrix(QQ, projected).solve_left(moved)
        shear_values = coordinates * phi - phi
        rows_on_basis = [list(row) for row in shear_values.rows()]
        rows_on_basis += [[0] * rank_radical for _ in range(quotient_rank - len(vectors))]
        results.append(assemble(action, identity_radical, rows_on_basis))
    unit_rows = tuple(identity_radical.rows())
    for images in general_linear_images(unit_rows):
        unipotent = matrix(ZZ, [list(row) for row in images]).transpose()
        shear_values = phi - phi * unipotent.transpose()
        rows_on_basis = [list(row) for row in shear_values.rows()]
        rows_on_basis += [[0] * rank_radical for _ in range(quotient_rank - len(vectors))]
        results.append(assemble(identity_matrix(ZZ, quotient_rank), unipotent, rows_on_basis))
    for row in range(len(vectors), quotient_rank):
        for column in range(rank_radical):
            rows_on_basis = zero_shear()
            rows_on_basis[row][column] = 1
            results.append(assemble(identity_matrix(ZZ, quotient_rank), identity_radical, rows_on_basis))
    return tuple(results)
