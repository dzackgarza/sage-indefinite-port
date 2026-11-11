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


def test_hermitian_quaternion_transporters():
    gram, sources = _hermitian_quaternion_plane()
    commutant = IsometryCommutant(gram, sources)
    zero = int(0) * ONE
    witness = _pair_matrix(((ONE, ONE), (zero, ONE)))
    moved = witness.inverse().transpose() * gram * witness.inverse()
    assert witness.det() in (int(1), int(-1)) and moved.denominator() == int(1)
    assert all(witness * F == F * witness for F in sources)
    t = equivariant_transporter(commutant, sage_matrix(SageZZ, moved), sources)
    assert t is not None and _is_transporter(t, gram, moved, sources, sources)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


def _golden_tensor_lattice():
    r"""Return `O^4 \otimes O^2`, `O = Z[\phi]`, `\phi^2 = \phi + 1`, with `b = Tr(w' \otimes w)`, `w' = H + H`, `w = diag(1, -1)`, the
    sources `1 \otimes R`, `1 \otimes D` for `R = [[\sqrt 5, -2], [-2, \sqrt 5]]` (`5 - 4 = 1`, so `R^t w R = w`), `D = diag(1, -1)`, and the map
    `P \otimes Q \mapsto GL_{16}(Z)`. `R + R^{-1} = 2 \sqrt 5`, and `R - \sqrt 5` is `-2` times the swap, so the sources generate `1 \otimes M_2(K)`."""
    field = QuadraticField(int(5), "s")
    root = field.gen()
    phi = (int(1) + root) / int(2)
    outer = sage_matrix(field, int(4), int(4))
    for a in (int(0), int(2)):
        outer[a, a + int(1)] = outer[a + int(1), a] = int(1)
    inner = sage_matrix(field, int(2), int(2), [int(1), int(0), int(0), int(-1)])
    weights = outer.tensor_product(inner)
    frame = []
    for i in range(int(8)):
        for scalar in (field.one(), phi):
            v = [field.zero()] * int(8)
            v[i] = scalar
            frame.append(sage_matrix(field, int(8), int(1), v))
    gram = sage_matrix(SageZZ, [[SageZZ((x.transpose() * weights * y)[int(0), int(0)].trace()) for y in frame] for x in frame])
    coordinates = sage_matrix(SageQQ, [[int(1), int(0)], [SageQQ(1) / int(2), SageQQ(1) / int(2)]])

    def on_lattice(block):
        images = []
        for v in frame:
            image = block * v
            images.append([c for i in range(int(8)) for c in coordinates.solve_left(sage_vector(SageQQ, image[i, int(0)].list()))])
        return sage_matrix(SageZZ, images).transpose()

    one = sage_identity(field, int(4))
    moves = [sage_matrix(field, int(2), int(2), [root, int(-2), int(-2), root]), inner]
    return gram, [on_lattice(one.tensor_product(g)) for g in moves], on_lattice, field, phi


def test_golden_orthogonal_centralizer_rank_four():
    gram, sources, _on_lattice, _field, _phi = _golden_tensor_lattice()
    generators = centralizer_generators(IsometryCommutant(gram, sources))
    one = sage_identity(SageZZ, int(16))
    assert all(g.transpose() * gram * g == gram and all(g * F == F * g for F in sources) for g in generators)
    assert any(g != one and (g - one) ** int(16) == int(0) for g in generators)
    assert any(g * h != h * g for g in generators for h in generators)


def test_golden_orthogonal_transporters_rank_four():
    gram, sources, on_lattice, field, phi = _golden_tensor_lattice()
    commutant = IsometryCommutant(gram, sources)
    shear = sage_identity(field, int(4))
    shear[int(0), int(2)] = phi
    witness = on_lattice(shear.tensor_product(sage_identity(field, int(2))))
    moved = witness.inverse().transpose() * gram * witness.inverse()
    assert witness.det() in (int(1), int(-1)) and moved.denominator() == int(1)
    t = equivariant_transporter(commutant, sage_matrix(SageZZ, moved), sources)
    assert t is not None and _is_transporter(t, gram, moved, sources, sources)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


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
