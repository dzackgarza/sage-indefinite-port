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
"""

from sage_extensions.all import *
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import centralizer_generators, equivariant_transporter

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
