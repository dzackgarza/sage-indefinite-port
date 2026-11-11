r"""GL-type summand transporters with a non-maximal multiplier ring `Z[\sqrt 5] \subset Z[\varphi]` (#54).

`\varepsilon = 2 + \sqrt 5` has `\varepsilon (\sqrt 5 - 2) = 1`, minimal polynomial `x^2 - 4x - 1`, and `Z[\varepsilon] = Z[\sqrt 5]`. On `X = Z^2`
let `\varepsilon` act by `A = [[2, 5], [1, 2]]`: from `e_1 = 1`, `A e_1 = (2, 1)` gives `e_2 = \varepsilon - 2 = \sqrt 5`, so `X \cong Z[\sqrt 5]`. By
`A_2 = [[1, 2], [2, 3]]`: `A_2 e_1 = (1, 2)` gives `e_2 = (\varepsilon - 1)/2 = \varphi`, so `X_2 \cong Z[\varphi]`. Both determinants are `-1`.
`L = Z^6` with `B = [[0, I], [I, 0]] \oplus I_2` and `F = A \oplus A^{-t} \oplus 1`; the eigenvalue polynomials `x^2 - 4x - 1`,
`x^2 + 4x - 1` and `x - 1` differ, so the commutant has three central idempotents, the first two swapped by the adjoint, and it is not
commutative (it contains `M_2(Q)` on the last summand). An integral transporter to `F_2 = A_2 \oplus A_2^{-t} \oplus 1` would restrict to a
`Z[\varepsilon]`-isomorphism `Z[\sqrt 5] \to Z[\varphi]`, impossible since their multiplier rings `Z[\sqrt 5]` and `Z[\varphi]` differ. The
isometry `t = P \oplus P^{-t} \oplus 1`, `P = [[1, 1], [0, 1]]`, carries `F` to `t F t^{-1}`, so a transporter to it exists."""

from sage.matrix.constructor import block_diagonal_matrix, matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.integer_ring import ZZ as SageZZ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import equivariant_transporter


def _isometry(A):
    return block_diagonal_matrix(A, A.inverse().transpose(), sage_identity(SageZZ, int(2))).change_ring(SageZZ)


def test_gl_type_transporter_over_a_nonmaximal_order():
    hyperbolic = sage_matrix(SageZZ, [[int(0), int(0), int(1), int(0)], [int(0), int(0), int(0), int(1)], [int(1), int(0), int(0), int(0)], [int(0), int(1), int(0), int(0)]])
    gram = block_diagonal_matrix(hyperbolic, sage_identity(SageZZ, int(2)))
    A = sage_matrix(SageZZ, [[int(2), int(5)], [int(1), int(2)]])
    A2 = sage_matrix(SageZZ, [[int(1), int(2)], [int(2), int(3)]])
    F, F2 = _isometry(A), _isometry(A2)
    assert F.transpose() * gram * F == gram and F2.transpose() * gram * F2 == gram
    commutant = IsometryCommutant(gram, [F])
    assert equivariant_transporter(commutant, gram, [F2]) is None
    P = sage_matrix(SageZZ, [[int(1), int(1)], [int(0), int(1)]])
    t = _isometry(P)
    moved = t * F * t.inverse()
    found = equivariant_transporter(commutant, gram, [moved])
    assert found is not None
    assert found.transpose() * gram * found == gram and found * F == moved * found
