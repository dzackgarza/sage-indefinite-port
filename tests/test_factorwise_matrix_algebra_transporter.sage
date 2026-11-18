r"""GL-type summand transporters when `Q[F]|_X = M_2(Q)` is not commutative (#54).

`S = [[0, 1], [1, 0]]` and `T = [[1, 2], [0, 1]]` generate `M_2(Q)` (`T - 1 = 2 E_{12}`, `S E_{12} S = E_{21}`), and the ring
`\Lambda = Z[S, T] \subsetneq M_2(Z)`. On `X = Z^2 \oplus Z^2` they act by `S \oplus S`, `T \oplus T`. The sublattice
`L' = \{(x, y) : x \equiv y \bmod 2\}` of `Z^2` is `\Lambda`-stable; in its basis `(1, 1), (1, -1)` they act by `S' = \mathrm{diag}(1, -1)` and
`T' = [[2, -1], [1, 0]]`. Let `X_2 = Z^2 \oplus L'` (`S \oplus S'`, `T \oplus T'`), `L = Z^8`, `B_L = [[0, I_4], [I_4, 0]]`, `F = E \oplus E^{-t}`.
A `\Lambda`-isomorphism `X \to X_2` is `M_2(Q)`-linear, so it carries `M_2(Z) X = X` onto `M_2(Z) X_2 = Z^2 \oplus Z^2 \supsetneq X_2`
(`E_{11} (1, 1) = (1, 0)`), but `[M_2(Z) X : X] = 1 \ne 2 = [M_2(Z) X_2 : X_2]`: no transporter exists. A conjugate by `P \oplus P^{-t}`,
`P \in GL_4(Z)`, has one."""

from sage.matrix.constructor import block_diagonal_matrix, block_matrix, matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.integer_ring import ZZ as SageZZ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import equivariant_transporter


def _isometry(B):
    return block_diagonal_matrix(B, B.inverse().transpose()).change_ring(SageZZ)


def test_gl_type_transporter_over_a_matrix_algebra():
    zero, one = sage_matrix(SageZZ, int(4), int(4)), sage_identity(SageZZ, int(4))
    gram = block_matrix(SageZZ, [[zero, one], [one, zero]], subdivide=False)
    S = sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(0)]])
    T = sage_matrix(SageZZ, [[int(1), int(2)], [int(0), int(1)]])
    S2 = sage_matrix(SageZZ, [[int(1), int(0)], [int(0), int(-1)]])
    T2 = sage_matrix(SageZZ, [[int(2), int(-1)], [int(1), int(0)]])
    basis = sage_matrix(SageZZ, [[int(1), int(1)], [int(1), int(-1)]])
    assert S * basis == basis * S2 and T * basis == basis * T2
    F = [_isometry(block_diagonal_matrix(S, S)), _isometry(block_diagonal_matrix(T, T))]
    F2 = [_isometry(block_diagonal_matrix(S, S2)), _isometry(block_diagonal_matrix(T, T2))]
    assert all(G.transpose() * gram * G == gram for G in F + F2)
    commutant = IsometryCommutant(gram, F)
    assert equivariant_transporter(commutant, gram, F2) is None
    P = sage_matrix(SageZZ, [[int(1), int(1), int(0), int(2)], [int(0), int(1), int(1), int(0)], [int(0), int(0), int(1), int(1)], [int(0), int(0), int(0), int(1)]])
    t = _isometry(P)
    moved = [t * G * t.inverse() for G in F]
    found = equivariant_transporter(commutant, gram, moved)
    assert found is not None
    assert found.transpose() * gram * found == gram and all(found * G == H * found for G, H in zip(F, moved))
