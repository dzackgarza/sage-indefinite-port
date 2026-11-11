r"""GL-type summand transporters of rank two over the non-maximal order `Z[\sqrt 5]` (#54).

`\varepsilon = 2 + \sqrt 5` is a unit of `Z[\sqrt 5]` (`(2 + \sqrt 5)(\sqrt 5 - 2) = 1`) with `Z[\varepsilon] = Z[\sqrt 5]`. On `Z[\sqrt 5]` (basis
`1, \sqrt 5`) it acts by `A = [[2, 5], [1, 2]]`, on `Z[\varphi]` (basis `1, \varphi`, `\sqrt 5 = 2\varphi - 1`) by `A_2 = [[1, 2], [2, 3]]`; both
have determinant `-1`. Let `X = Z[\sqrt 5]^2` (`\varepsilon` by `A \oplus A`) and `X_2 = Z[\sqrt 5] \oplus Z[\varphi]` (`A \oplus A_2`) in
`L = Z^8`, `B_L = [[0, I_4], [I_4, 0]]`, `F = B \oplus B^{-t}`. Both have multiplier ring `Z[\sqrt 5]` and
`O_K X = O_K X_2 = Z[\varphi]^2`, but `[O_K X : X] = 4` while `[O_K X_2 : X_2] = [Z[\varphi] : Z[\sqrt 5]] = 2`; an isomorphism
`X \to X_2` is `K`-linear, so it carries `O_K X` onto `O_K X_2` and preserves this index: no transporter exists. A conjugate by `P \oplus P^{-t}`, `P \in GL_4(Z)`, has one."""

from sage.matrix.constructor import block_diagonal_matrix, block_matrix, matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.integer_ring import ZZ as SageZZ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import equivariant_transporter


def _isometry(B):
    return block_diagonal_matrix(B, B.inverse().transpose()).change_ring(SageZZ)


def test_gl_type_transporter_of_rank_two_over_a_nonmaximal_order():
    zero, one = sage_matrix(SageZZ, int(4), int(4)), sage_identity(SageZZ, int(4))
    gram = block_matrix(SageZZ, [[zero, one], [one, zero]], subdivide=False)
    A = sage_matrix(SageZZ, [[int(2), int(5)], [int(1), int(2)]])
    A2 = sage_matrix(SageZZ, [[int(1), int(2)], [int(2), int(3)]])
    F = _isometry(block_diagonal_matrix(A, A))
    F2 = _isometry(block_diagonal_matrix(A, A2))
    assert F.transpose() * gram * F == gram and F2.transpose() * gram * F2 == gram
    commutant = IsometryCommutant(gram, [F])
    assert equivariant_transporter(commutant, gram, [F2]) is None
    P = sage_matrix(SageZZ, [[int(1), int(1), int(0), int(2)], [int(0), int(1), int(1), int(0)], [int(0), int(0), int(1), int(1)], [int(0), int(0), int(0), int(1)]])
    t = _isometry(P)
    moved = t * F * t.inverse()
    found = equivariant_transporter(commutant, gram, [moved])
    assert found is not None
    assert found.transpose() * gram * found == gram and found * F == moved * found
