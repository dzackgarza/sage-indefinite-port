r"""GL-type summand transporters with `X \cong Z[\varphi]^2`, of rank two over `Q(\sqrt 5)` (#54).

`A = [[0, 1], [1, 1]]` is multiplication by `\varphi` on `Z[\varphi] = Z 1 + Z \varphi` (`\varphi \cdot 1 = \varphi`, `\varphi^2 = 1 + \varphi`), of
determinant `-1`. On `X = Z^4` the isometry acts by `B = A \oplus A`, so `X \cong Z[\varphi]^2` and `Q[F]|_X = Q(\sqrt 5)` has dimension
`2 < 4`: the rank-two case of the Morita reduction. `L = Z^8` with `B_L = [[0, I_4], [I_4, 0]]` and `F = B \oplus B^{-t}`; the eigenvalue
polynomials `x^2 - x - 1` and `x^2 + x - 1` differ, so the commutant `M_2(K) \times M_2(K)` is not commutative and its centre is not a field.
For `P \in GL_4(Z)`, `t = P \oplus P^{-t}` is an isometry carrying `F` to `t F t^{-1}`, so a transporter exists; the target summand
`P X` is again a `Z[\varphi]`-module of rank two, with Steinitz class trivial as `Z[\varphi]` has class number one."""

from sage.matrix.constructor import block_diagonal_matrix, block_matrix, matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.integer_ring import ZZ as SageZZ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import equivariant_transporter


def test_gl_type_transporter_of_rank_two_over_the_golden_order():
    zero, one = sage_matrix(SageZZ, int(4), int(4)), sage_identity(SageZZ, int(4))
    gram = block_matrix(SageZZ, [[zero, one], [one, zero]], subdivide=False)
    A = sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(1)]])
    B = block_diagonal_matrix(A, A)
    F = block_diagonal_matrix(B, B.inverse().transpose()).change_ring(SageZZ)
    assert F.transpose() * gram * F == gram
    P = sage_matrix(SageZZ, [[int(1), int(1), int(0), int(2)], [int(0), int(1), int(1), int(0)], [int(0), int(0), int(1), int(1)], [int(0), int(0), int(0), int(1)]])
    t = block_diagonal_matrix(P, P.inverse().transpose()).change_ring(SageZZ)
    moved = t * F * t.inverse()
    commutant = IsometryCommutant(gram, [F])
    found = equivariant_transporter(commutant, gram, [moved])
    assert found is not None
    assert found.transpose() * gram * found == gram and found * F == moved * found
