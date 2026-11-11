r"""Noncompact commutant centralizers and transporters (#54).

Unipotent commutant. L = U + <2> with basis e, f, g, e.f = 1, g.g = 2. The Eichler transvection
E = E_{e,g} fixes e, sends g to g - 2e and f to f + g - e; it is one unipotent Jordan block, so
A = Q[E - 1] and A / rad = Q, and C(Q) = {+-E^s}. The reflection s in g gives s E s = E_{e,-g} = E^{-1}.
The rational isometry D: e -> 2e, f -> f/2, g -> g gives D E D^{-1} = E_{2e,g} = E^2, so every rational
t with t E = E^2 t is +-D E^s, and D E^s f = f/2 + s g - 2 s^2 e is never integral: E and E^2 are not
conjugate in O(L).

Symplectic quaternion commutant. D = (-1, 3)_Q is a division algebra (ramified at 2 and 3), L = Z<1,i,j,k>
with b(x, y) = trd(x conj(y)), and F_1, F_2 the left multiplications by i and 2 + j, of norm one. They
generate L_D, so A is the right multiplications R_D, with adjoint R_v -> R_{conj(v)}. Hence
C(Z) = {R_v : v in Z<1,i,j,k>, n(v) = 1}. Conjugation by q = 1 + i maps i, j, k to i, k, -j, so it is an
integral isometry carrying (L_i, L_{2+j}) to (L_i, L_{2+k}).

Orthogonal quaternion commutant. L = Z<1,i,j,k>^2 with b(x, y) = trd((conj(x_1) y_2 - conj(x_2) y_1) j). The
matrices [[1,1],[0,1]], [[1,0],[1,1]], i 1 and (2 + j) 1 acting on the left preserve b and generate
M_2(D), so A = R_D with adjoint R_v -> R_{j conj(v) j^{-1}}. A unitary R_v with n(v) = 1 commutes with j,
so v = +-(2 + j)^n, the norm-one units of Z[sqrt 3]. One with n(v) = -1 anticommutes with j, v = x i + y k
with x^2 - 3 y^2 = -1, which has no solution modulo 3. So C(Z) = <-1> x <R_{2+j}>.

A transporter to the Gram matrix 3B has determinant squared 1/3^r, so there is none.

Levi commutant. On L = M + M, M = U + <2>, F = E + E has commutant M_2(Q[E - 1]), with radical
M_2((E - 1)) and Levi quotient M_2(Q). With s the reflection in g of the second summand,
t = 1 + s carries F to E + E^{-1}, so the Levi transporter exists.

Gaussian hyperbolic plane. L = Z[i]^2 with b = Re h, h(x, y) = x_1 conj(y_2) + x_2 conj(y_1), and F
the multiplication by i. The commutant is M_2(Q(i)), centre Q(i), with the second-kind adjoint of h,
and C(Z) = U(1,1; Z[i]). Its scalar units u I, u in mu_4, are the generators D_u of the model
(D_u = diag(u, conj(u)^{-1}) = u I), and they are F, -1, -F in every basis; T is unipotent.

Real quadratic symplectic. L = O^2 (x) O^2, O = Z[sqrt 2], with b = Tr(w (x) w), w the standard
alternating form, and F_i = 1 (x) g for g in E_12(1), E_21(1), diag(u, u^{-1}), u = 1 + sqrt 2. These
generate 1 (x) M_2(K), so A = M_2(K) (x) 1 with the symplectic adjoint of w and C(Z) = SL_2(O) (x) 1.
The model generator diag(u, u^{-1}) has eigenvalues u, u^{-1} and their conjugates, each twice on Q^8,
so characteristic polynomial (x^2 - 2x - 1)^2 (x^2 + 2x - 1)^2 for every choice of fundamental unit.

Gaussian hyperbolic transporters. On the basis (1,0), (i,0), (0,1), (0,i) of Z[i]^2, complex conjugation c is
diag(1, -1, 1, -1); it preserves b = Re h, since h(cx, cy) = conj h(x, y), and c F c^{-1} = -F. A Z[i]-linear
P in GL_2(Z[i]) commutes with F, so P is a transporter from (B, F) to (P^{-t} B P^{-1}, F). The Gram 3B has
determinant 81 != 1, so no transporter to it exists.

Anisotropic Gaussian plane. L = Z[i]^2 with h = diag(1, -3), b = Re h = diag(1, 1, -3, -3), F the
multiplication by i. h is anisotropic: |x|^2 = 3 |y|^2 has no nonzero solution, since 3 is not a sum of
two rational squares. So C(Z) = U(h) over Z[i], from the Bestvina-Savin ocean. It contains [[2, 3], [1, 2]]
(columns of h-values 1, -3 and orthogonal), of infinite order, so the generators include an element of
infinite order; a finite-order element of GL_4(Z) has order at most 12. The vectors v = (1 + i, 1), w = (3, 1 - i)
have h(v) = -1, h(w) = 3, h(v, w) = 0 and det [v w] = -1, so [v w] is an isometry of (Z[i]^2, h) onto
(Z[i]^2, -h): a transporter to (Re(-h), F) = (diag(-1, -1, 3, 3), F) exists. As above c carries F to -F.

Real quadratic orthogonal plane. L = O^2 (x) O^2, O = Z[sqrt 2], with b = Tr(w (x) w), w = diag(1, -1), and
F_i = 1 (x) g for g = R = [[r, 1], [1, r]], r = sqrt 2, and D = diag(1, -1); both preserve w. R + R^{-1} = 2r, so
the F_i generate 1 (x) M_2(K) and A = M_2(K) (x) 1 with the orthogonal adjoint of w: C(Z) = O(O^2, w) (x) 1.
D R D = R^{-1} != R, so C(Z) is not abelian, and R has eigenvalues r +- 1, of infinite order.

Real quadratic orthogonal transporters. In the plane above, A = M_2(K) (x) 1 and Lambda = O^2 (x) O^2, so a transporter
between the forms Tr(w (x) w) and Tr(w' (x) w) with the same F_i is P (x) 1 with P in GL_2(O) and P^t w' P = w. The shear
P = [[1, r], [0, 1]] gives one by construction. For w' = diag(u, -u^{-1}), u = 1 + r, of the same determinant, w would represent
u: x^2 - y^2 = u forces x - y = e, x + y = u / e with e = +-u^k, and 2x = e + u/e is congruent to 1 + u = 2 + r modulo 2 O
(u^2 = 3 + 2r is 1 modulo 2), which is not 0. So no transporter exists.

Nonmaximal Levi order. In M + M with coordinates (e_1, f_1, g_1, e_2, f_2, g_2), the lattice L' spanned by e_1, e_2, g_1, g_2,
f_1, 2 f_2 is stable under F = E + E (E e = e, E g = g - 2e, E f = f + g - e). Its order reduces modulo the radical to the
matrices in M_2(Z) that preserve Z + 2Z, the Eichler order of level 2, not M_2(Z). A unimodular t_0 = 1 + E_{16} is a
transporter from (B', F) to (t_0^{-t} B' t_0^{-1}, t_0 F t_0^{-1}) by construction; the Gram 3 B' has another determinant.

Product Levi quotient. On M + <1> with basis e, f, g, h, F = E + 1 has A / rad = Q x Q (one Jordan block of length three on
M and a fixed h). t_0 = 1 + E_{14} is a transporter from (B, F) to (t_0^{-t} B t_0^{-1}, t_0 F t_0^{-1}) by construction;
3B has another determinant, so no transporter to it exists.

Anisotropic Eisenstein plane. L = O^2, O = Z[w], w^2 + w + 1 = 0, with h = diag(1, -2) and b = Tr h; on the basis
(1,0), (w,0), (0,1), (0,w), b = A_2 + A_2(-2), and F, the multiplication by w, is [[0, -1], [1, -1]] on each summand.
h is anisotropic: 2 is inert in Q(w), so x conj(x) = 2 y conj(y) has no nonzero solution (compare 2-adic valuations).
The commutant is M_2(Q(w)) with the adjoint of h, and C(Z) = U(O^2, h), from the quaternion order of SU(h). It contains
X = [[a, 2 conj(c)], [c, conj(a)]] for a = 1 - w, c = 1, since N(a) - 2 N(c) = 3 - 2 = 1; its trace 3 exceeds 2, so X has
infinite order. Coordinatewise conjugation c(x) = conj(x) preserves b and carries F to F^2. The vectors v = (1, 1),
u = (2, 1) have h(v) = -1, h(u) = 2, h(u, v) = 2 - 2 = 0 and det [v u] = -1, so [v u] is an isometry of (O^2, h) onto
(O^2, -h): a transporter to (-b, F) exists. The Gram 3b has another determinant, so no transporter to it exists.

Hermitian quaternion plane. L = Delta^2, Delta = Z<1,i,j,k> in D = (-1, 3)_Q, b(x, y) = trd(conj(y_1) x_2 + conj(y_2) x_1), and F_1, F_2
the right multiplications by i and 2 + j (reduced norm one). They generate R_D, so A = M_2(D) acting on the left, with the adjoint
of the Hermitian H = [[0, 1], [1, 0]], and C(Z) = U(Delta^2, H). It contains n(i) = [[1, i], [0, 1]] (i is pure), unipotent and not 1,
and n(i), W = [[0, 1], [1, 0]] do not commute. The left matrix [[1, 1], [0, 1]] is in GL_2(Delta) and commutes with the F_i, so it is a
transporter from (B, F) to its pushed-forward Gram; 3B has another determinant.

Golden orthogonal rank four. L = O^4 (x) O^2, O = Z[phi], phi^2 = phi + 1, b = Tr(w' (x) w), w' = H + H, w = diag(1, -1), F_i = 1 (x) R,
1 (x) D with R = [[sqrt 5, -2], [-2, sqrt 5]], D = diag(1, -1). Then A = M_4(K) (x) 1 with the orthogonal adjoint of w', C(Z) = O(O^4, H + H)
(x) 1. 2 and 3 are inert in K (residue fields F_4, F_9), so the Kneser rank conditions are void, and H + H is even unimodular. C(Z)
contains the Eichler transvections t(e, a), unipotent and not 1, and t(e, e_1), t(f, e_1) do not commute. P (x) 1 for the shear
P = 1 + phi E_13 in GL_4(O) is a transporter from (B, F) to its pushed-forward Gram; 3B has another determinant.

Split orthogonal M_2(Q). L = U, F = -1: A = M_2(Q), O = M_2(Z), and the adjoint is the adjoint of the form U on the simple
module Z^2, so C(Z) = O(U) = {+-1, +-swap}, a group of order four. On L = <1> + <1> + <-1> with F = -1, C(Z) = O(I_{2,1}),
which contains the reflections s_v, s_{e_1} in v = (1, 1, 1) and e_1 (both of norm 1). s_{e_1} s_v = [[1, 2, -2], [-2, -1, 2],
[-2, -2, 3]] lies in SO(2, 1) and has trace 3, so all its eigenvalues are 1; it is not 1, so it is a nontrivial unipotent
element, of infinite order, and C(Z) is infinite.

Split symplectic M_2(Q). L = Z^2 (x) Z^2 with B = J (x) J, J = [[0, 1], [-1, 0]], symmetric of determinant one, and F_i = 1 (x) g for
g = [[1, 1], [0, 1]], [[1, 0], [1, 1]], which preserve J and generate 1 (x) M_2(Q). So A = M_2(Q) (x) 1 with the adjoint of J,
symplectic, and C(Z) = Sp_2(Z) (x) 1 = SL_2(Z) (x) 1, which contains the unipotent transvections.

Two Pell units. L = (<1> + <-2>) + (<1> + <-3>), F_1 = [[3, 4], [2, 3]] + 1 and F_2 = 1 + [[2, 3], [1, 2]] (rows), the
multiplications by 3 + 2 sqrt 2 and 2 + sqrt 3. The F_i have distinct eigenvalues on the two summands, so
A = Q(sqrt 2) x Q(sqrt 3), each factor fixed by the adjoint, which is the Galois conjugation. y -> -y on the first summand
conjugates F_1 to F_1^{-1} and fixes F_2, so a transporter to (B, F_1^{-1}, F_2) exists; with y -> -y on both summands,
to (B, F_1^{-1}, F_2^{-1}). A transporter to B' = (<1> + <-2>) + (<-1> + <3>) with the same F_i commutes with them, so it
preserves both summands and maps <1> + <-3> isometrically onto <-1> + <3>: x^2 - 3 y^2 = -1 has no solution modulo 3,
so there is none, although the determinants agree. 3B has another determinant.

Swapped factors. L = U + U with basis (e_1, e_2, f_1, f_2), B = [[0, 1], [1, 0]] in 2 x 2 blocks, and F = M + M^{-t}, M = [[1, 1],
[1, 0]], which preserves B. M has eigenvalues phi, -1/phi and M^{-t} = M^{-1} has 1/phi, -phi, so A = Q[M] x Q[M^{-1}] =
Q(sqrt 5) x Q(sqrt 5), and the adjoint swaps the factors. The block swap s preserves B and carries F to s F s = M^{-1} + M,
so a transporter to (B, s F s) exists.

Vector orbits. L = U, F = -1, so C(Z) = O(U) = {+-1, +-swap}. The orbit of v = (1, 1) is {+-v}; its stabilizer
is {1, swap}, and (1, -1) has another norm, so it is not in the orbit.

Split orthogonal transporters. L = <1> + <-1>, F = -1, so A = M_2(Q) and every isometry commutes with F. For g = [[2, 1], [1, 1]]
(det 1) the matrix g is a transporter from (B, F) to (g^{-t} B g^{-1}, F). The hyperbolic plane U has the same determinant -1, but
it is even and B represents 1, so no transporter to (U, F) exists.

Split symplectic transporters. On L = Z^2 (x) Z^2 with B = J (x) J and F_i = 1 (x) g, the matrix t = 1 (x) Q, Q = diag(1, -1), has
Q J Q = -J, so it carries (B, F_i) to (-B, 1 (x) Q g Q).

Definite Hermitian quaternion plane. As the Hermitian quaternion plane, for D = (-1, -1)_Q (definite, ramified at 2) and the Lipschitz
order Delta = Z<1,i,j,k>; the right multiplications by i and j have reduced norm one and generate R_D. Delta^x = {+-1, +-i, +-j, +-k}
is finite, and R^0 = Z i + Z j + Z k with nrd(x i + y j + z k) = x^2 + y^2 + z^2 has covering radius sqrt(3)/2 < 1 for the
sup-nearest translate, so the height reduction applies. n(i) is unipotent and not 1, and does not commute with W. The left matrix
[[1, 1], [0, 1]] commutes with the F_i, so it is a transporter from (B, F) to its pushed-forward Gram; 3B has another determinant.

Ramified dyadic Eichler generators. K = Q(sqrt 2), where 2 = (sqrt 2)^2 ramifies; O_K = Z[sqrt 2] is principal and K totally real.
Lambda = H + H + H + <2 sqrt 2> is even, and rank_P(Lambda) = 6 at P = (sqrt 2), residue field F_2; 3 is inert. The vector
v = sqrt 2 b of the last basis vector b has h(v, v) = 4 sqrt 2 and h(v, Lambda) = 4 O_K, so 2 h(v, Lambda) / h(v, v) = sqrt 2 O_K is
integral and sigma_v is an integral reflection with c = -2 h(v, v) = -8 sqrt 2 of odd P-valuation 7. So the spinor-norm bound
at P is realized and generators are returned; one of them is a reflection whose norm has odd P-valuation.
"""

from sage_extensions.all import *
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import block_matrix, matrix as sage_matrix
from sage.modules.free_module_element import vector as sage_vector
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import VectorOrbit, centralizer_generators, equivariant_transporter

ALGEBRA = QuaternionAlgebra(SageQQ, int(-1), int(3))
ONE, I, J, K = ALGEBRA.basis()
BASIS = (ONE, I, J, K)


def _eichler():
    gram = sage_matrix(SageZZ, [[int(0), int(1), int(0)], [int(1), int(0), int(0)], [int(0), int(0), int(2)]])
    E = sage_matrix(SageZZ, [[int(1), int(-1), int(-2)], [int(0), int(1), int(0)], [int(0), int(1), int(1)]])
    return gram, E


def _is_transporter(t, gram, target_gram, sources, targets):
    return (
        t.det() in (int(1), int(-1))
        and t.transpose() * target_gram * t == gram
        and all(t * F == T * t for F, T in zip(sources, targets))
    )


def test_unipotent_transporter_to_the_inverse_exists():
    gram, E = _eichler()
    commutant = IsometryCommutant(gram, [E])
    target = E.inverse().change_ring(SageZZ)
    t = equivariant_transporter(commutant, gram, [target])
    assert t is not None and _is_transporter(t, gram, gram, [E], [target])


def test_unipotent_transporter_to_the_square_is_absent():
    gram, E = _eichler()
    commutant = IsometryCommutant(gram, [E])
    assert equivariant_transporter(commutant, gram, [E * E]) is None


def _coordinates(q):
    return [SageQQ(c) for c in list(q)]


def _left(u):
    return sage_matrix(SageZZ, [_coordinates(u * b) for b in BASIS]).transpose()


def _norm_lattice():
    gram = sage_matrix(SageZZ, [[int((x * y.conjugate()).reduced_trace()) for y in BASIS] for x in BASIS])
    return gram, [_left(I), _left(int(2) + J)]


def _right_multiplier(g):
    r"""Return ``v`` when ``g`` is the right multiplication by ``v`` on Z<1,i,j,k>."""
    v = sum(c * b for c, b in zip(g.column(int(0)), BASIS))
    assert g == sage_matrix(SageZZ, [_coordinates(b * v) for b in BASIS]).transpose()
    return v


def test_symplectic_quaternion_centralizer_is_the_norm_one_units():
    gram, sources = _norm_lattice()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    values = [_right_multiplier(g) for g in generators]
    assert all(v.reduced_norm() == int(1) and all(c in SageZZ for c in v) for v in values)
    assert any(abs(v.reduced_trace()) > int(2) for v in values)


def test_symplectic_quaternion_transporters():
    gram, sources = _norm_lattice()
    commutant = IsometryCommutant(gram, sources)
    q = int(1) + I
    conjugated = [_left(q * I * ~q), _left(q * (int(2) + J) * ~q)]
    assert conjugated == [_left(I), _left(int(2) + K)]
    t = equivariant_transporter(commutant, gram, conjugated)
    assert t is not None and _is_transporter(t, gram, gram, sources, conjugated)
    P = _left(int(2) + J)
    moved = [P * F * P.inverse() for F in sources]
    t = equivariant_transporter(commutant, gram, moved)
    assert t is not None and _is_transporter(t, gram, gram, sources, moved)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


PAIRS = tuple((b, int(0) * ONE) for b in BASIS) + tuple((int(0) * ONE, b) for b in BASIS)


def _pair_form(x, y):
    return ((x[0].conjugate() * y[1] - x[1].conjugate() * y[0]) * J).reduced_trace()


def _pair_coordinates(x):
    return _coordinates(x[0]) + _coordinates(x[1])


def _act(g, x):
    return (g[0][0] * x[0] + g[0][1] * x[1], g[1][0] * x[0] + g[1][1] * x[1])


def _pair_matrix(g):
    return sage_matrix(SageZZ, [_pair_coordinates(_act(g, x)) for x in PAIRS]).transpose()


def _orthogonal_lattice():
    gram = sage_matrix(SageZZ, [[int(_pair_form(x, y)) for y in PAIRS] for x in PAIRS])
    zero, one, u = int(0) * ONE, ONE, int(2) + J
    actions = [
        ((one, one), (zero, one)),
        ((one, zero), (one, one)),
        ((I, zero), (zero, I)),
        ((u, zero), (zero, u)),
    ]
    return gram, [_pair_matrix(g) for g in actions], _pair_matrix(((u, zero), (zero, u)))


def test_orthogonal_quaternion_centralizer_is_the_real_torus():
    gram, sources, _P = _orthogonal_lattice()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    unit = int(2) + J
    rows = [[int(2), int(0)]]
    for g in generators:
        v = sum(c * b for c, b in zip(g.column(int(0))[: int(4)], BASIS))
        assert list(g.column(int(0))[int(4) :]) == [int(0)] * int(4)
        assert g == sage_matrix(SageZZ, [_pair_coordinates((x[0] * v, x[1] * v)) for x in PAIRS]).transpose()
        exponent = next(
            (n, s) for n in range(int(-20), int(21)) for s in (int(1), int(-1)) if v == s * unit**n
        )
        rows.append([int(0) if exponent[1] == int(1) else int(1), exponent[0]])
    assert sage_matrix(SageZZ, rows).echelon_form().submatrix(int(0), int(0), int(2), int(2)).det() == int(1)


def test_orthogonal_quaternion_transporters():
    gram, sources, P = _orthogonal_lattice()
    commutant = IsometryCommutant(gram, sources)
    moved = [P * F * P.inverse() for F in sources]
    t = equivariant_transporter(commutant, gram, moved)
    assert t is not None and _is_transporter(t, gram, gram, sources, moved)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


def _gaussian_plane():
    entries = [int(c) for c in (0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0)]
    gram = sage_matrix(SageZZ, int(4), int(4), entries)
    rotation = [int(c) for c in (0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0)]
    return gram, sage_matrix(SageZZ, int(4), int(4), rotation)


def test_gaussian_hyperbolic_centralizer_contains_the_scalar_units():
    gram, F = _gaussian_plane()
    generators = centralizer_generators(IsometryCommutant(gram, [F]))
    one = sage_identity(SageZZ, int(4))
    assert all(g * F == F * g and g.transpose() * gram * g == gram for g in generators)
    assert F in generators and -one in generators and -F in generators
    assert any(g != one and (g - one) ** int(2) == int(0) for g in generators)


def _symplectic_tensor_lattice():
    field = QuadraticField(int(2), "s")
    root = field.gen()
    omega = ((int(0), int(1)), (int(-1), int(0)))

    def pairing(x, y):
        return sum(
            x[int(2) * a + b] * y[int(2) * c + d] * omega[a][c] * omega[b][d]
            for a in range(int(2)) for b in range(int(2)) for c in range(int(2)) for d in range(int(2))
        )

    frame = []
    for i in range(int(4)):
        for scalar in (field.one(), root):
            v = [field.zero()] * int(4)
            v[i] = scalar
            frame.append(v)
    gram = sage_matrix(SageZZ, [[SageZZ(pairing(x, y).trace()) for y in frame] for x in frame])

    def on_lattice(block):
        images = []
        for v in frame:
            image = block * sage_matrix(field, int(4), int(1), v)
            images.append([SageQQ(c) for i in range(int(4)) for c in image[i, int(0)].list()])
        return sage_matrix(SageZZ, images).transpose()

    one = sage_identity(field, int(2))
    unit = int(1) + root
    moves = [
        sage_matrix(field, int(2), int(2), [int(1), int(1), int(0), int(1)]),
        sage_matrix(field, int(2), int(2), [int(1), int(0), int(1), int(1)]),
        sage_matrix(field, int(2), int(2), [unit, int(0), int(0), ~unit]),
    ]
    return gram, [on_lattice(one.tensor_product(g)) for g in moves]


def test_real_quadratic_symplectic_centralizer():
    gram, sources = _symplectic_tensor_lattice()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    one = sage_identity(SageZZ, int(8))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert -one in generators
    x = PolynomialRing(SageQQ, "x").gen()
    expected = (x**int(2) - int(2) * x - int(1)) ** int(2) * (x**int(2) + int(2) * x - int(1)) ** int(2)
    assert any(g.charpoly() == expected for g in generators)


GAUSSIAN = QuadraticField(int(-1), "i")


def _gaussian_linear(entries):
    r"""Return the integral matrix of the Z[i]-linear map ``entries`` on the basis (1,0), (i,0), (0,1), (0,i)."""
    g = sage_matrix(GAUSSIAN, int(2), int(2), entries)
    i = GAUSSIAN.gen()
    frame = [(int(1), int(0)), (i, int(0)), (int(0), int(1)), (int(0), i)]
    columns = []
    for x in frame:
        image = g * sage_matrix(GAUSSIAN, int(2), int(1), list(x))
        columns.append([SageZZ(c) for k in range(int(2)) for c in image[k, int(0)].list()])
    return sage_matrix(SageZZ, columns).transpose()


CONJUGATION = sage_matrix(SageZZ, int(4), int(4), [int(c) for c in (1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1)])


def test_gaussian_hyperbolic_transporters():
    gram, F = _gaussian_plane()
    commutant = IsometryCommutant(gram, [F])
    assert CONJUGATION * F * CONJUGATION.inverse() == -F
    t = equivariant_transporter(commutant, gram, [-F])
    assert t is not None and _is_transporter(t, gram, gram, [F], [-F])
    i = GAUSSIAN.gen()
    P = _gaussian_linear([int(1), int(1) + i, int(0), int(1)])
    assert P * F == F * P
    moved = (P.inverse().transpose() * gram * P.inverse()).change_ring(SageZZ)
    t = equivariant_transporter(commutant, moved, [F])
    assert t is not None and _is_transporter(t, gram, moved, [F], [F])
    assert equivariant_transporter(commutant, int(3) * gram, [F]) is None


def _ocean_plane():
    gram = sage_matrix(SageZZ, int(4), int(4), int(0))
    for k, value in enumerate((1, 1, -3, -3)):
        gram[k, k] = int(value)
    _hyperbolic, F = _gaussian_plane()
    return gram, F


def test_anisotropic_gaussian_centralizer_is_infinite():
    gram, F = _ocean_plane()
    generators = centralizer_generators(IsometryCommutant(gram, [F]))
    assert all(g * F == F * g and g.transpose() * gram * g == gram for g in generators)
    # An integral 4 x 4 matrix of finite order has roots of unity as eigenvalues, so |trace| <= 4.
    assert any(abs((g * h).trace()) > int(4) for g in generators for h in generators)


def test_anisotropic_gaussian_transporters():
    gram, F = _ocean_plane()
    commutant = IsometryCommutant(gram, [F])
    t = equivariant_transporter(commutant, gram, [-F])
    assert t is not None and _is_transporter(t, gram, gram, [F], [-F])
    negative = -gram
    i = GAUSSIAN.gen()
    witness = _gaussian_linear([int(1) + i, int(3), int(1), int(1) - i])
    assert witness.det() in (int(1), int(-1)) and witness.transpose() * negative * witness == gram
    t = equivariant_transporter(commutant, negative, [F])
    assert t is not None and _is_transporter(t, gram, negative, [F], [F])


def _orthogonal_tensor_lattice(first):
    r"""Return the Gram of `Tr(w' (x) w)` on `O^2 (x) O^2`, `w' = diag(` ``first(root)`` `)`, the sources `1 (x) R`, `1 (x) D`, and the map `P (x) Q -> GL_8(Z)`."""
    field = QuadraticField(int(2), "s")
    root = field.gen()
    signs = (int(1), int(-1))
    weights = first(root)

    def pairing(x, y):
        return sum(x[int(2) * a + b] * y[int(2) * a + b] * weights[a] * signs[b] for a in range(int(2)) for b in range(int(2)))

    frame = []
    for i in range(int(4)):
        for scalar in (field.one(), root):
            v = [field.zero()] * int(4)
            v[i] = scalar
            frame.append(v)
    gram = sage_matrix(SageZZ, [[SageZZ(pairing(x, y).trace()) for y in frame] for x in frame])

    def on_lattice(block):
        images = []
        for v in frame:
            image = block * sage_matrix(field, int(4), int(1), v)
            images.append([SageQQ(c) for i in range(int(4)) for c in image[i, int(0)].list()])
        return sage_matrix(SageZZ, images).transpose()

    one = sage_identity(field, int(2))
    moves = [
        sage_matrix(field, int(2), int(2), [root, int(1), int(1), root]),
        sage_matrix(field, int(2), int(2), [int(1), int(0), int(0), int(-1)]),
    ]
    return gram, [on_lattice(one.tensor_product(g)) for g in moves], on_lattice


def _hyperbolic_weights(root):
    return (int(1), int(-1))


def test_real_quadratic_orthogonal_centralizer():
    gram, sources, _on_lattice = _orthogonal_tensor_lattice(_hyperbolic_weights)
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert any(g * h != h * g for g in generators for h in generators)
    assert any(not all(f.is_cyclotomic() for f, _m in g.charpoly().factor()) for g in generators)


def test_vector_orbit_in_the_hyperbolic_plane():
    gram = sage_matrix(SageZZ, int(2), int(2), [int(0), int(1), int(1), int(0)])
    F = -sage_identity(SageZZ, int(2))
    orbit = VectorOrbit(gram, [F], [int(1), int(1)])
    swap = sage_matrix(SageZZ, int(2), int(2), [int(0), int(1), int(1), int(0)])
    stabilizer = orbit.stabilizer_generators()
    assert all(g * orbit.vector == orbit.vector and g.transpose() * gram * g == gram for g in stabilizer)
    assert swap in stabilizer
    t = orbit.transporter([int(-1), int(-1)])
    assert t is not None and t.transpose() * gram * t == gram and list(t * orbit.vector) == [int(-1), int(-1)]
    assert orbit.transporter([int(1), int(-1)]) is None


def test_real_quadratic_orthogonal_transporters():
    gram, sources, on_lattice = _orthogonal_tensor_lattice(_hyperbolic_weights)
    commutant = IsometryCommutant(gram, sources)
    field = QuadraticField(int(2), "s")
    root = field.gen()
    one = sage_identity(field, int(2))
    shear = sage_matrix(field, int(2), int(2), [int(1), root, int(0), int(1)])
    witness = on_lattice(shear.tensor_product(one))
    moved = witness.inverse().transpose() * gram * witness.inverse()
    assert witness.det() in (int(1), int(-1)) and moved.denominator() == int(1)
    t = equivariant_transporter(commutant, sage_matrix(SageZZ, moved), sources)
    assert t is not None and _is_transporter(t, gram, moved, sources, sources)
    other, other_sources, _ = _orthogonal_tensor_lattice(lambda r: (int(1) + r, int(1) - r))
    assert other_sources == sources and abs(other.det()) == abs(gram.det())
    assert equivariant_transporter(commutant, other, sources) is None


def test_product_levi_quotient_transporter():
    gram = sage_matrix(SageZZ, int(4), int(4), [int(c) for c in (0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1)])
    F = sage_matrix(SageZZ, int(4), int(4), [int(c) for c in (1, -1, -2, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1)])
    commutant = IsometryCommutant(gram, [F])
    witness = sage_identity(SageZZ, int(4)).__copy__()
    witness[int(0), int(3)] = int(1)
    target = witness * F * witness.inverse()
    target_gram = witness.inverse().transpose() * gram * witness.inverse()
    t = equivariant_transporter(commutant, target_gram, [target])
    assert t is not None and _is_transporter(t, gram, target_gram, [F], [target])
    assert equivariant_transporter(commutant, int(3) * gram, [F]) is None


EISENSTEIN = QuadraticField(int(-3), "s")
OMEGA = (EISENSTEIN.gen() - int(1)) / int(2)


def _eisenstein_linear(entries):
    r"""Return the integral matrix of the O-linear map ``entries`` on the basis (1,0), (w,0), (0,1), (0,w)."""
    g = sage_matrix(EISENSTEIN, int(2), int(2), entries)
    frame = [(int(1), int(0)), (OMEGA, int(0)), (int(0), int(1)), (int(0), OMEGA)]
    columns = []
    for x in frame:
        image = g * sage_matrix(EISENSTEIN, int(2), int(1), list(x))
        for k in range(int(2)):
            entry = image[k, int(0)]
            second = (entry - entry.conjugate()) / (OMEGA - OMEGA.conjugate())
            columns.append(SageZZ(entry - second * OMEGA))
            columns.append(SageZZ(second))
    return sage_matrix(SageZZ, int(4), int(4), columns).transpose()


def _eisenstein_plane():
    entries = (2, -1, 0, 0, -1, 2, 0, 0, 0, 0, -4, 2, 0, 0, 2, -4)
    gram = sage_matrix(SageZZ, int(4), int(4), [int(c) for c in entries])
    F = sage_matrix(SageZZ, int(4), int(4), [int(c) for c in (0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0, -1, 0, 0, 1, -1)])
    assert F == _eisenstein_linear([OMEGA, int(0), int(0), OMEGA])
    return gram, F


def test_anisotropic_eisenstein_centralizer_is_infinite():
    gram, F = _eisenstein_plane()
    a = int(1) - OMEGA
    X = _eisenstein_linear([a, int(2), int(1), a.conjugate()])
    assert X * F == F * X and X.transpose() * gram * X == gram and X.trace() == int(6)
    generators = centralizer_generators(IsometryCommutant(gram, [F]))
    assert all(g * F == F * g and g.transpose() * gram * g == gram for g in generators)
    assert any(abs((g * h).trace()) > int(4) for g in generators for h in generators)


def test_anisotropic_eisenstein_transporters():
    gram, F = _eisenstein_plane()
    commutant = IsometryCommutant(gram, [F])
    conjugation = sage_matrix(SageZZ, int(4), int(4), [int(c) for c in (1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, -1)])
    assert conjugation * F * conjugation.inverse() == F * F and conjugation.transpose() * gram * conjugation == gram
    t = equivariant_transporter(commutant, gram, [F * F])
    assert t is not None and _is_transporter(t, gram, gram, [F], [F * F])
    negative = -gram
    witness = _eisenstein_linear([int(1), int(2), int(1), int(1)])
    assert witness.det() in (int(1), int(-1)) and witness.transpose() * negative * witness == gram
    t = equivariant_transporter(commutant, negative, [F])
    assert t is not None and _is_transporter(t, gram, negative, [F], [F])
    assert equivariant_transporter(commutant, int(3) * gram, [F]) is None


def _hermitian_quaternion_plane():
    r"""Return `\Delta^2`, `\Delta = Z<1, i, j, k>` in `(-1, 3)_Q`, with `b(x, y) = trd(conj(y_1) x_2 + conj(y_2) x_1)` and the right
    multiplications by `i` and `2 + j` (reduced norm one, so they preserve `b`; they generate `R_D`, whose commutant is `M_2(D)`)."""

    def form(x, y):
        return (y[int(0)].conjugate() * x[int(1)] + y[int(1)].conjugate() * x[int(0)]).reduced_trace()

    gram = sage_matrix(SageZZ, [[int(form(x, y)) for y in PAIRS] for x in PAIRS])

    def right(d):
        return sage_matrix(SageZZ, [_pair_coordinates((x[int(0)] * d, x[int(1)] * d)) for x in PAIRS]).transpose()

    return gram, [right(I), right(int(2) + J)]


def test_hermitian_quaternion_centralizer_is_infinite():
    gram, sources = _hermitian_quaternion_plane()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    one = sage_identity(SageZZ, int(8))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert any(g != one and (g - one) ** int(8) == int(0) for g in generators)
    assert any(g * h != h * g for g in generators for h in generators)


def _golden_tensor_lattice():
    r"""Return `O^4 \otimes O^2`, `O = Z[\phi]`, `\phi^2 = \phi + 1`, with `b = Tr(w' \otimes w)`, `w' = H + H`, `w = diag(1, -1)`, the
    sources `1 \otimes R`, `1 \otimes D` for `R = [[\sqrt 5, -2], [-2, \sqrt 5]]` (`5 - 4 = 1`, so `R^t w R = w`), `D = diag(1, -1)`, and the map
    `P \otimes Q \mapsto GL_{16}(Z)`. `R + R^{-1} = 2 \sqrt 5`, and `R - \sqrt 5` is `-2` times the swap, so the sources generate `1 \otimes M_2(K)`."""
    field = QuadraticField(int(int(5)), "s")
    root = field.gen()
    phi = (int(int(1)) + root) / int(int(2))
    outer = sage_matrix(field, int(int(4)), int(int(4)))
    for a in (int(int(0)), int(int(2))):
        outer[a, a + int(int(1))] = outer[a + int(int(1)), a] = int(int(1))
    inner = sage_matrix(field, int(int(2)), int(int(2)), [int(int(1)), int(int(0)), int(int(0)), int(-int(1))])
    weights = outer.tensor_product(inner)
    frame = []
    for i in range(int(int(8))):
        for scalar in (field.one(), phi):
            v = [field.zero()] * int(int(8))
            v[i] = scalar
            frame.append(sage_matrix(field, int(int(8)), int(int(1)), v))
    gram = sage_matrix(SageZZ, [[SageZZ((x.transpose() * weights * y)[int(int(0)), int(int(0))].trace()) for y in frame] for x in frame])
    coordinates = sage_matrix(SageQQ, [[int(int(1)), int(int(0))], [SageQQ(int(1)) / int(int(2)), SageQQ(int(1)) / int(int(2))]])

    def on_lattice(block):
        images = []
        for v in frame:
            image = block * v
            images.append([c for i in range(int(int(8))) for c in coordinates.solve_left(sage_vector(SageQQ, image[i, int(int(0))].list()))])
        return sage_matrix(SageZZ, images).transpose()

    one = sage_identity(field, int(int(4)))
    moves = [sage_matrix(field, int(int(2)), int(int(2)), [root, int(-int(2)), int(-int(2)), root]), inner]
    return gram, [on_lattice(one.tensor_product(g)) for g in moves], on_lattice, field, phi


def test_golden_orthogonal_centralizer_rank_four():
    gram, sources, _on_lattice, _field, _phi = _golden_tensor_lattice()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    one = sage_identity(SageZZ, int(16))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert any(g != one and (g - one) ** int(16) == int(0) for g in generators)
    assert any(g * h != h * g for g in generators for h in generators)


def _closure(generators):
    one = generators[int(0)] ** int(0)
    elements, frontier = [one], [one]
    while frontier:
        current = frontier.pop()
        for g in generators:
            image = current * g
            if image not in elements:
                elements.append(image)
                frontier.append(image)
    return elements


def test_split_orthogonal_centralizer_of_minus_one_on_U():
    gram = sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(0)]])
    generators = centralizer_generators(IsometryCommutant(gram, [-sage_identity(SageZZ, int(2))]))
    one = sage_identity(SageZZ, int(2))
    swap = sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(0)]])
    group = _closure(list(generators))
    assert len(group) == int(4)
    assert all(g in (one, -one, swap, -swap) for g in group)


def test_split_orthogonal_centralizer_rank_three_is_infinite():
    from sage.groups.matrix_gps.finitely_generated import MatrixGroup

    gram = sage_matrix(SageZZ, [[int(1), int(0), int(0)], [int(0), int(1), int(0)], [int(0), int(0), int(-1)]])
    generators = centralizer_generators(IsometryCommutant(gram, [-sage_identity(SageZZ, int(3))]))
    assert all(g.transpose() * gram * g == gram for g in generators)
    assert not MatrixGroup(list(generators)).is_finite()


def test_split_symplectic_centralizer_is_SL2():
    J = sage_matrix(SageZZ, [[int(0), int(1)], [int(-1), int(0)]])
    one = sage_identity(SageZZ, int(2))
    gram = J.tensor_product(J)
    sources = [one.tensor_product(sage_matrix(SageZZ, rows)) for rows in ([[int(1), int(1)], [int(0), int(1)]], [[int(1), int(0)], [int(1), int(1)]])]
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    blocks = []
    for g in generators:
        P = g.matrix_from_rows_and_columns([int(0), int(2)], [int(0), int(2)])
        assert g == P.tensor_product(one) and P.det() == int(1)
        blocks.append(P)
    assert any(P != one and (P - one) ** int(2) == int(0) for P in blocks)


def _two_pell_lattice():
    gram = sage_matrix(SageZZ, [[int(1), int(0), int(0), int(0)], [int(0), int(-2), int(0), int(0)], [int(0), int(0), int(1), int(0)], [int(0), int(0), int(0), int(-3)]])
    F1 = sage_matrix(SageZZ, [[int(3), int(4), int(0), int(0)], [int(2), int(3), int(0), int(0)], [int(0), int(0), int(1), int(0)], [int(0), int(0), int(0), int(1)]])
    F2 = sage_matrix(SageZZ, [[int(1), int(0), int(0), int(0)], [int(0), int(1), int(0), int(0)], [int(0), int(0), int(2), int(3)], [int(0), int(0), int(1), int(2)]])
    return gram, [F1, F2]


def test_two_pell_units_transporters():
    gram, sources = _two_pell_lattice()
    F1, F2 = sources
    commutant = IsometryCommutant(gram, sources)
    for targets in ([F1.inverse(), F2], [F1.inverse(), F2.inverse()]):
        targets = [T.change_ring(SageZZ) for T in targets]
        t = equivariant_transporter(commutant, gram, targets)
        assert t is not None and _is_transporter(t, gram, gram, sources, targets)
    other = sage_matrix(SageZZ, [[int(1), int(0), int(0), int(0)], [int(0), int(-2), int(0), int(0)], [int(0), int(0), int(-1), int(0)], [int(0), int(0), int(0), int(3)]])
    assert all(F.transpose() * other * F == other for F in sources) and other.det() == gram.det()
    assert equivariant_transporter(commutant, other, sources) is None
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


def test_swapped_factor_transporter():
    zero, one = sage_matrix(SageZZ, int(2), int(2)), sage_identity(SageZZ, int(2))
    gram = block_matrix(SageZZ, [[zero, one], [one, zero]], subdivide=False)
    M = sage_matrix(SageZZ, [[int(1), int(1)], [int(1), int(0)]])
    F = block_matrix(SageZZ, [[M, zero], [zero, M.inverse().transpose().change_ring(SageZZ)]], subdivide=False)
    swap = block_matrix(SageZZ, [[zero, one], [one, zero]], subdivide=False)
    assert F.transpose() * gram * F == gram
    commutant = IsometryCommutant(gram, [F])
    assert commutant.dimension() == int(4) and all(x * y == y * x for x in commutant.order_basis for y in commutant.order_basis)
    target = swap * F * swap
    t = equivariant_transporter(commutant, gram, [target])
    assert t is not None and _is_transporter(t, gram, gram, [F], [target])
    assert equivariant_transporter(commutant, int(3) * gram, [F]) is None


def test_split_orthogonal_transporters():
    gram = sage_matrix(SageZZ, [[int(1), int(0)], [int(0), int(-1)]])
    minus = -sage_identity(SageZZ, int(2))
    commutant = IsometryCommutant(gram, [minus])
    g = sage_matrix(SageZZ, [[int(2), int(1)], [int(1), int(1)]])
    target = (g.inverse().transpose() * gram * g.inverse()).change_ring(SageZZ)
    t = equivariant_transporter(commutant, target, [minus])
    assert t is not None and _is_transporter(t, gram, target, [minus], [minus])
    hyperbolic = sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(0)]])
    assert equivariant_transporter(commutant, hyperbolic, [minus]) is None


def test_split_symplectic_transporters():
    J = sage_matrix(SageZZ, [[int(0), int(1)], [int(-1), int(0)]])
    one = sage_identity(SageZZ, int(2))
    Q = sage_matrix(SageZZ, [[int(1), int(0)], [int(0), int(-1)]])
    gram = J.tensor_product(J)
    moves = [sage_matrix(SageZZ, rows) for rows in ([[int(1), int(1)], [int(0), int(1)]], [[int(1), int(0)], [int(1), int(1)]])]
    sources = [one.tensor_product(m) for m in moves]
    targets = [one.tensor_product(Q * m * Q) for m in moves]
    commutant = IsometryCommutant(gram, sources)
    t = equivariant_transporter(commutant, -gram, targets)
    assert t is not None and _is_transporter(t, gram, -gram, sources, targets)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


DEFINITE = QuaternionAlgebra(SageQQ, int(-1), int(-1))
DEFINITE_BASIS = DEFINITE.basis()
DEFINITE_PAIRS = tuple((b, int(0) * DEFINITE_BASIS[int(0)]) for b in DEFINITE_BASIS) + tuple(
    (int(0) * DEFINITE_BASIS[int(0)], b) for b in DEFINITE_BASIS
)


def _definite_coordinates(x):
    return [SageQQ(c) for c in list(x[int(0)])] + [SageQQ(c) for c in list(x[int(1)])]


def _definite_hermitian_plane():
    def form(x, y):
        return (y[int(0)].conjugate() * x[int(1)] + y[int(1)].conjugate() * x[int(0)]).reduced_trace()

    gram = sage_matrix(SageZZ, [[int(form(x, y)) for y in DEFINITE_PAIRS] for x in DEFINITE_PAIRS])

    def right(d):
        return sage_matrix(SageZZ, [_definite_coordinates((x[int(0)] * d, x[int(1)] * d)) for x in DEFINITE_PAIRS]).transpose()

    _one, i, j, _k = DEFINITE_BASIS
    shear = sage_matrix(SageZZ, [_definite_coordinates((x[int(0)] + x[int(1)], x[int(1)])) for x in DEFINITE_PAIRS]).transpose()
    return gram, [right(i), right(j)], shear


def test_definite_hermitian_quaternion_centralizer_is_infinite():
    gram, sources, _shear = _definite_hermitian_plane()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    one = sage_identity(SageZZ, int(8))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert any(g != one and (g - one) ** int(8) == int(0) for g in generators)
    assert any(g * h != h * g for g in generators for h in generators)


def test_ramified_dyadic_eichler_generators():
    from sage_indefinite_port.indefinite import eichler_orthogonal

    field = QuadraticField(int(2), "r")
    root = field.gen()
    gram = sage_matrix(field, int(7), int(7))
    for a in (int(0), int(2), int(4)):
        gram[a, a + int(1)] = gram[a + int(1), a] = int(1)
    gram[int(6), int(6)] = int(2) * root
    generators = eichler_orthogonal.generators(field, gram)
    one = sage_identity(field, int(7))
    prime = field.ideal(root)
    assert all(g.transpose() * gram * g == gram and all(c.is_integral() for c in g.list()) for g in generators)
    odd = []
    for g in generators:
        if (g - one).rank() == int(1) and g.det() == int(-1):
            v = (g - one).column_space().basis()[int(0)]
            norm = (sage_matrix(field, [v]) * gram * sage_matrix(field, [v]).transpose())[int(0), int(0)]
            odd.append(norm.valuation(prime) % int(2) == int(1))
    assert any(odd)


QUADS = tuple((x, (int(0) * ONE, int(0) * ONE)) for x in PAIRS) + tuple(((int(0) * ONE, int(0) * ONE), x) for x in PAIRS)


def _skew(x, y):
    return y[int(0)].conjugate() * x[int(1)] - y[int(1)].conjugate() * x[int(0)]


def _quad_coordinates(X):
    return _pair_coordinates(X[int(0)]) + _pair_coordinates(X[int(1)])


def _quad_matrix(move):
    return sage_matrix(SageZZ, [_quad_coordinates(move(X)) for X in QUADS]).transpose()


def _skew_hermitian_quaternion_lattice():
    r"""Return `\Delta^2 \otimes Z^2` with `b = \mathrm{trd} \circ (h \otimes \omega)`, `h(x, y) = \bar y_1 x_2 - \bar y_2 x_1` skew-Hermitian,
    `\omega` the standard alternating form, and sources generating `D^{op} \otimes M_2(Q)`: the right multiplications by `i`, `2 + j` and the
    two shears of `SL_2(Z)` on the `Z^2` factor."""

    def form(X, Y):
        return (_skew(X[int(0)], Y[int(1)]) - _skew(X[int(1)], Y[int(0)])).reduced_trace()

    gram = sage_matrix(SageZZ, [[int(form(X, Y)) for Y in QUADS] for X in QUADS])

    def right(d):
        return _quad_matrix(lambda X: tuple((x[int(0)] * d, x[int(1)] * d) for x in X))

    def add(x, y):
        return (x[int(0)] + y[int(0)], x[int(1)] + y[int(1)])

    sources = [
        right(I),
        right(int(2) + J),
        _quad_matrix(lambda X: (add(X[int(0)], X[int(1)]), X[int(1)])),
        _quad_matrix(lambda X: (X[int(0)], add(X[int(1)], X[int(0)]))),
    ]
    return gram, sources


def _diagonal_left(g):
    return _quad_matrix(lambda X: (_act(g, X[int(0)]), _act(g, X[int(1)])))


def test_skew_hermitian_quaternion_centralizer_is_infinite():
    gram, sources = _skew_hermitian_quaternion_lattice()
    assert gram == gram.transpose() and gram.det() != int(0)
    assert all(F.transpose() * gram * F == gram for F in sources)
    shear = _diagonal_left(((ONE, ONE), (int(0) * ONE, ONE)))
    assert shear.transpose() * gram * shear == gram and all(shear * F == F * shear for F in sources)
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    one = sage_identity(SageZZ, int(16))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert any(g != one and (g - one) ** int(16) == int(0) for g in generators)
    assert any(g * h != h * g for g in generators for h in generators)


def test_diagonal_null_vector_beyond_ternary_subforms():
    r"""`<1, 1, 1, 1, -7>` over `Q(sqrt 2)` is isotropic (`2^2 + 1 + 1 + 1 = 7`); `<1, 1, 1, -7>` over `Q` is not (`7` is not a sum of three
    rational squares), so a null vector needs all five coordinates or the quaternary reduction."""
    from sage_indefinite_port.indefinite.commutant_groups import _diagonal_null

    field = QuadraticField(int(2), "r")
    values = [int(1), int(1), int(1), int(1), int(-7)]
    null = _diagonal_null(field, values)
    assert null is not None and any(z != int(0) for z in null)
    assert sum(a * z**int(2) for a, z in zip(values, null)) == int(0)


def test_dyadic_spinor_norm_groups_of_binary_two_adic_lattices():
    r"""Over `K = Q(sqrt 17)`, where `2` splits and `K_P = Q_2`. `<1, 3>`: unit vectors of norms `1, 3, 7 = 2^2 + 3, 13 = 1 + 3 \cdot 2^2` give
    integral symmetries, so products of two of them give every unit class; Beli (Def. 3 (III), `R = 0`, `d(-3) = 2 = 2e`) gives
    `O^x F^2`, so `2 \notin \theta`. `<1, 2>`: `R = 1` gives `N(-2)`, the norms of `Q_2(sqrt -2)`: `2 = 0^2 + 2 \cdot 1^2` and `3 = 1 + 2` are
    norms, while `(-2, -1)_2 = (-2, 5)_2 = -1`."""
    from sage_indefinite_port.indefinite.eichler_orthogonal import LocalSquareClasses, spinor_norm_group

    field = QuadraticField(int(17), "t")
    prime = field.ideal(int(2)).factor()[int(0)][int(0)]
    classes = LocalSquareClasses(field, prime)
    odd = spinor_norm_group(field, sage_matrix(field, [[int(1), int(0)], [int(0), int(3)]]), prime, classes)
    assert all(classes(x) in odd for x in (int(-1), int(3), int(5)))
    assert classes(int(2)) not in odd
    even = spinor_norm_group(field, sage_matrix(field, [[int(1), int(0)], [int(0), int(2)]]), prime, classes)
    assert classes(int(2)) in even and classes(int(3)) in even
    assert classes(int(-1)) not in even and classes(int(5)) not in even


def test_dyadic_spinor_norm_group_after_maximal_norm_splitting():
    r"""`<1> + A(0, 0)` over `Q(sqrt 17)_P = Q_2`: the computed blocks `<1>`, `A(0, 0)` share the scale `O` with norms `O`, `2O`, so they are
    not a maximal norm splitting and Beli's exchange (Lemma 1.4.5) is needed. The lattice is `<1, 1, -1>`; the symmetries of `e_1` (`Q = 1`)
    and of `e + f` (`Q = 2`, `2 B(e + f, L) / 2 = Z_2`) are integral, so `2 \\in \\theta(O^+)`, and Lemma 1.7.2(ii) (`ord nH = 1`, `ord n<1> = 0`)
    gives `\\theta(O^+) = Q_2^x`."""
    from sage_indefinite_port.indefinite.eichler_orthogonal import LocalSquareClasses, spinor_norm_group

    field = QuadraticField(int(17), "t")
    prime = field.ideal(int(2)).factor()[int(0)][int(0)]
    classes = LocalSquareClasses(field, prime)
    gram = sage_matrix(field, [[int(1), int(0), int(0)], [int(0), int(0), int(1)], [int(0), int(1), int(0)]])
    theta = spinor_norm_group(field, gram, prime, classes)
    assert classes(int(2)) in theta
    assert theta == classes.space


def _lipschitz_over_sqrt2():
    r"""Return `(-1, -1)` over `Z = Q(sqrt 2)` and the order `Z[sqrt 2]<1, i, j, k>`. The algebra is totally definite (both real places
    give the Hamilton quaternions), so the canonical involution is positive."""
    from sage_indefinite_port.indefinite.quaternion_voronoi_number_field import QuaternionOrderVoronoi

    field = QuadraticField(int(2), "w")
    algebra = QuaternionAlgebra(field, int(-1), int(-1))
    basis = [field.gen() ** p * g for g in algebra.basis() for p in (int(0), int(1))]
    return algebra, QuaternionOrderVoronoi(algebra, basis)


def test_voronoi_perfect_forms_over_real_quadratic_field():
    algebra, voronoi = _lipschitz_over_sqrt2()
    forms = voronoi.perfect_forms()
    assert len(forms) >= int(1)
    for s in forms:
        assert voronoi.involution(s) == s
        assert voronoi.gram(s).is_positive_definite()
        assert voronoi.is_perfect(s)


def test_voronoi_sl1_of_lipschitz_order_over_sqrt2():
    r"""`SL_1` of `Z[sqrt 2]<1, i, j, k>` in the totally definite `(-1, -1)` over `Q(sqrt 2)`: `a^2 + b^2 + c^2 + d^2 = 1` at both real places
    forces every coordinate `p + q sqrt 2` to satisfy `|p| <= 1` and `|q| sqrt 2 <= 1`, so `q = 0` and `SL_1 = Q_8 = \\{\\pm 1, \\pm i, \\pm j, \\pm k\\}`."""
    algebra, voronoi = _lipschitz_over_sqrt2()
    _one, i, j, k = algebra.basis()
    units = voronoi.unit_generators()
    assert all(voronoi.is_unit(u) for u in units)
    generators = voronoi.sl1_generators()
    assert all(voronoi.in_order(u) and u.reduced_norm() == int(1) for u in generators)
    one = algebra(int(1))
    elements, frontier = [one], [one]
    while frontier:
        current = frontier.pop()
        for g in generators:
            image = current * g
            if image not in elements:
                elements.append(image)
                frontier.append(image)
    assert len(elements) == int(8)
    assert i in elements and j in elements and k in elements


def test_anisotropic_skew_hermitian_quaternion_centralizer():
    r"""`\Delta^2 \otimes Z^2` with `b = \mathrm{trd} \circ (h \otimes \omega)`, `h(x, y) = \bar y_1 i x_1 + \bar y_2 (3i + j) x_2` over `(-1, 3)_Q`.
    `nrd(i) = 1`, `nrd(3i + j) = 9 - 3 = 6`, and `6` is not a rational square, so `h` has no null vector (a null vector forces
    `nrd(u)^2 nrd(i) = nrd(3i + j)`). The sources are those of the hyperbolic fixture, so `A = M_2(D)` with the adjoint of `h`. Every generator
    is an integral isometry commuting with the sources."""
    alpha, beta = I, int(3) * I + J

    def form(X, Y):
        def h(x, y):
            return y[int(0)].conjugate() * alpha * x[int(0)] + y[int(1)].conjugate() * beta * x[int(1)]

        return (h(X[int(0)], Y[int(1)]) - h(X[int(1)], Y[int(0)])).reduced_trace()

    gram = sage_matrix(SageZZ, [[int(form(X, Y)) for Y in QUADS] for X in QUADS])
    assert gram == gram.transpose() and gram.det() != int(0)
    _skew_gram, sources = _skew_hermitian_quaternion_lattice()
    assert all(F.transpose() * gram * F == gram for F in sources)
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)


def test_anisotropic_skew_hermitian_quaternion_transporters():
    r"""In the anisotropic plane `h = diag(i, 3i + j)`, the diagonal left action `P` of `[[1, i], [0, 1]] \in GL_2(\Delta)` commutes with the
    sources and is unimodular, so it is a transporter from `(B, F)` to `(P^{-t} B P^{-1}, F)`; `3B` has another determinant."""
    alpha, beta = I, int(3) * I + J

    def form(X, Y):
        def h(x, y):
            return y[int(0)].conjugate() * alpha * x[int(0)] + y[int(1)].conjugate() * beta * x[int(1)]

        return (h(X[int(0)], Y[int(1)]) - h(X[int(1)], Y[int(0)])).reduced_trace()

    gram = sage_matrix(SageZZ, [[int(form(X, Y)) for Y in QUADS] for X in QUADS])
    _skew_gram, sources = _skew_hermitian_quaternion_lattice()
    commutant = IsometryCommutant(gram, sources)
    P = _diagonal_left(((ONE, I), (int(0) * ONE, ONE)))
    assert P.det() in (int(1), int(-1)) and all(P * F == F * P for F in sources)
    moved = (P.inverse().transpose() * gram * P.inverse()).change_ring(SageZZ)
    t = equivariant_transporter(commutant, moved, sources)
    assert t is not None and _is_transporter(t, gram, moved, sources, sources)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


def test_maximal_orders_containing_the_lipschitz_order_over_sqrt2():
    r"""`(-1, -1)` over `Q(sqrt 2)` ramifies exactly at the two real places (even count, Hilbert reciprocity with no finite place), so a
    maximal order has `|\det Tr(o_i o_j)| = d_Z^4 = 8^4`. `Z[sqrt 2]<1, i, j, k>` has a larger discriminant (it misses `(1 + i)/sqrt 2`, which
    has reduced trace `sqrt 2` and reduced norm `1`), so it is not maximal, and every order returned contains it, is closed and integral, and
    has the maximal discriminant."""
    from sage_indefinite_port.indefinite.quaternion_orders_number_field import OrderLattices, maximal_orders_containing

    algebra, voronoi = _lipschitz_over_sqrt2()
    lattices = OrderLattices(algebra)
    order = list(voronoi.order)
    expected = int(8) ** int(4)
    assert lattices.discriminant(order) > expected
    root = algebra.base_ring().gen()
    hurwitz_like = (int(1) + algebra.basis()[int(1)]) / root
    assert lattices.integral(hurwitz_like)
    maximal = maximal_orders_containing(algebra, order, expected)
    assert len(maximal) >= int(1)
    for basis in maximal:
        assert lattices.discriminant(basis) == expected
        assert lattices.closure(basis) is not None and lattices.key(lattices.closure(basis)) == lattices.key(basis)
        assert lattices.key([*basis, *order]) == lattices.key(basis)


def test_witt_one_parabolic_generators_over_sqrt2():
    r"""`U + <2>` over `Q(sqrt 2)`. `<2>` is anisotropic and totally positive, and `O(<2>) = \\{\\pm 1\\}` (`u^2 = 1`). Every generator is an integral
    isometry. `v = e - f + m` is primitive isotropic (`2 \\cdot 1 \\cdot (-1) + 2 = 0`), and the Euclidean reduction returns an integral
    isometry `g` with `g v` a unit multiple of `e`. The transvection `E(e, m)` is unipotent and not `1`, so the group is infinite."""
    from sage_indefinite_port.indefinite.witt_one_number_field import WittOneLattice, eichler

    field = QuadraticField(int(2), "r")
    gram = sage_matrix(field, [[int(0), int(1), int(0)], [int(1), int(0), int(0)], [int(0), int(0), int(2)]])
    model = WittOneLattice(field, gram)
    generators = model.generators()
    assert all(g.transpose() * gram * g == gram and all(c.is_integral() for c in g.list()) for g in generators)
    one = sage_identity(field, int(3))
    transvection = eichler(gram, sage_vector(field, [int(1), int(0), int(0)]), sage_vector(field, [int(0), int(0), int(1)]))
    assert transvection != one and (transvection - one) ** int(3) == int(0)
    v = sage_vector(field, [int(1), int(-1), int(1)])
    g = model.reduce_isotropic(v)
    image = g * v
    assert g.transpose() * gram * g == gram and all(c.is_integral() for c in g.list())
    assert image[int(1)] == int(0) and image[int(2)] == int(0) and image[int(0)].is_unit()


def test_witt_one_parabolic_generators_mixed_sign_rank_one():
    r"""`U + <2 sqrt 5>` over `Q(sqrt 5)`: `2 sqrt 5` is positive at one real place and negative at the other, `M` is anisotropic (rank one)
    and even. `|N(2 sqrt 5 z^2)| = 20 N(z)^2 < 4` needs `|N(z)| < 1/sqrt 5`, which the Euclidean minimum `1/4` of `Q(sqrt 5)` allows,
    so the certificate holds. `v = sqrt 5 e - f + m` is isotropic (`2 sqrt 5 \\cdot (-1) + 2 sqrt 5 = 0`), primitive, and it reduces to a unit multiple of `e`."""
    from sage_indefinite_port.indefinite.witt_one_number_field import WittOneLattice

    field = QuadraticField(int(5), "s")
    root = field.gen()
    gram = sage_matrix(field, [[int(0), int(1), int(0)], [int(1), int(0), int(0)], [int(0), int(0), int(2) * root]])
    model = WittOneLattice(field, gram)
    generators = model.generators()
    assert all(g.transpose() * gram * g == gram and all(c.is_integral() for c in g.list()) for g in generators)
    v = sage_vector(field, [root, int(-1), int(1)])
    assert (sage_matrix([v]) * gram * sage_matrix([v]).transpose())[int(0), int(0)] == int(0)
    g = model.reduce_isotropic(v)
    image = g * v
    assert image[int(1)] == int(0) and image[int(2)] == int(0) and image[int(0)].is_unit()


def test_binary_orthogonal_group_of_one_minus_three_over_sqrt5():
    r"""`M = <1, -3>` over `Q(sqrt 5)`: `3` is not a square in `Q(sqrt 5)`, so `M` is anisotropic, of mixed sign at no place but indefinite at
    both. `SO(M)` is the integral norm-one torus of `L = Q(sqrt 5, sqrt 3)`, totally real of degree `4`, whose relative norm-one units have rank
    `3 - 1 = 2`; the only roots of unity in the totally real `L` are `\\pm 1`, so a proper generator other than `\\pm 1` has infinite order. The reflection in `m_1` (`h(m_1, m_1) = 1`) is integral and improper, so the improper coset is found."""
    from sage_indefinite_port.indefinite.binary_orthogonal_number_field import BinaryOrthogonalGroup

    field = QuadraticField(int(5), "s")
    gram = sage_matrix(field, [[int(1), int(0)], [int(0), int(-3)]])
    group = BinaryOrthogonalGroup(field, gram)
    generators = group.generators()
    assert all(g.transpose() * gram * g == gram and all(c.is_integral() for c in g.list()) for g in generators)
    assert any(g.det() == int(-1) for g in generators)
    one = sage_identity(field, int(2))
    assert any(g.det() == int(1) and g != one and g != -one for g in generators)


def _hyperbolic_plus_pell():
    r"""`L = U + <1> + <-2>` and `F = 1_U + P`, `P = [[3, 4], [2, 3]]` the Pell unit `3 + 2 sqrt 2` on `<1, -2>` (`P^t diag(1, -2) P =
    diag(1, -2)`). The commutant is `M_2(Q) \\times Q(sqrt 2)`, semisimple with a centre `Q \\times Q(sqrt 2)` of dimension 3."""
    gram = block_matrix(SageZZ, [[sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(0)]]), int(0)], [int(0), sage_matrix(SageZZ, [[int(1), int(0)], [int(0), int(-2)]])]], subdivide=False)
    F = block_matrix(SageZZ, [[sage_identity(SageZZ, int(2)), int(0)], [int(0), sage_matrix(SageZZ, [[int(3), int(4)], [int(2), int(3)]])]], subdivide=False)
    return gram, F


def test_factorwise_centralizer_and_transporter():
    gram, F = _hyperbolic_plus_pell()
    assert F.transpose() * gram * F == gram
    commutant = IsometryCommutant(gram, [F])
    generators = centralizer_generators(commutant)
    assert all(g.transpose() * gram * g == gram and g * F == F * g for g in generators)
    one = sage_identity(SageZZ, int(4))
    assert any(g.submatrix(int(2), int(2), int(2), int(2)) not in (sage_identity(SageZZ, int(2)), -sage_identity(SageZZ, int(2))) for g in generators)
    P = block_matrix(SageZZ, [[sage_matrix(SageZZ, [[int(1), int(1)], [int(0), int(1)]]), int(0)], [int(0), sage_identity(SageZZ, int(2))]], subdivide=False)
    assert P * F == F * P
    moved = (P.inverse().transpose() * gram * P.inverse()).change_ring(SageZZ)
    assert moved != gram
    t = equivariant_transporter(commutant, moved, [F])
    assert t is not None and _is_transporter(t, gram, moved, [F], [F])
    assert equivariant_transporter(commutant, int(3) * gram, [F]) is None


def test_factorwise_centralizer_with_a_gl_type_pair():
    r"""`L = (X + Y) + U` with `X = <e_1, e_2>`, `Y = <f_1, f_2>` paired by `B(e_i, f_j) = \\delta_{ij}`, and `F = A + A^{-t} + 1_U` with
    `A = [[1, 1], [1, 0]]` (characteristic polynomial `x^2 - x - 1`, not self-reciprocal, so `X` and `Y` are non-isomorphic simple `Q[F]`-modules).
    The commutant is `Q(sqrt 5) \\times Q(sqrt 5) \\times M_2(Q)`; the adjoint swaps the two field idempotents, a GL-type pair. `F` itself
    centralizes `F` and acts on `X` by `A`, which is not `\\pm 1`."""
    I2 = sage_identity(SageZZ, int(2))
    Z2 = sage_matrix(SageZZ, int(2), int(2))
    U = sage_matrix(SageZZ, [[int(0), int(1)], [int(1), int(0)]])
    gram = block_matrix(SageZZ, [[Z2, I2, Z2], [I2, Z2, Z2], [Z2, Z2, U]], subdivide=False)
    A = sage_matrix(SageZZ, [[int(1), int(1)], [int(1), int(0)]])
    F = block_matrix(SageZZ, [[A, Z2, Z2], [Z2, A.inverse().transpose().change_ring(SageZZ), Z2], [Z2, Z2, I2]], subdivide=False)
    assert F.transpose() * gram * F == gram
    generators = centralizer_generators(IsometryCommutant(gram, [F]))
    assert all(g.transpose() * gram * g == gram and g * F == F * g for g in generators)
    assert any(g.submatrix(int(0), int(0), int(2), int(2)) not in (I2, -I2) for g in generators)
