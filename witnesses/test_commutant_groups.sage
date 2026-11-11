r"""Commutant-group witnesses that need long runs (#54)."""

from sage_extensions.all import *
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import centralizer_generators, equivariant_transporter

def _eichler():
    gram = sage_matrix(SageZZ, [[int(0), int(1), int(0)], [int(1), int(0), int(0)], [int(0), int(0), int(2)]])
    E = sage_matrix(SageZZ, [[int(1), int(-1), int(-2)], [int(0), int(1), int(0)], [int(0), int(1), int(1)]])
    return gram, E
from sage.matrix.constructor import block_matrix, matrix as sage_matrix
from sage.modules.free_module_element import vector as sage_vector
from sage_indefinite_port.indefinite.commutant_groups import VectorOrbit, centralizer_generators, equivariant_transporter

PAIRS = tuple((b, int(0) * ONE) for b in BASIS) + tuple((int(0) * ONE, b) for b in BASIS)


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


def _pair_matrix(g):
    return sage_matrix(SageZZ, [_pair_coordinates(_act(g, x)) for x in PAIRS]).transpose()


def _act(g, x):
    return (g[0][0] * x[0] + g[0][1] * x[1], g[1][0] * x[0] + g[1][1] * x[1])


def _pair_coordinates(x):
    return _coordinates(x[0]) + _coordinates(x[1])


def _coordinates(q):
    return [SageQQ(c) for c in list(q)]


def _hermitian_quaternion_plane():
    r"""Return `\Delta^2`, `\Delta = Z<1, i, j, k>` in `(-1, 3)_Q`, with `b(x, y) = trd(conj(y_1) x_2 + conj(y_2) x_1)` and the right
    multiplications by `i` and `2 + j` (reduced norm one, so they preserve `b`; they generate `R_D`, whose commutant is `M_2(D)`)."""

    def form(x, y):
        return (y[int(0)].conjugate() * x[int(1)] + y[int(1)].conjugate() * x[int(0)]).reduced_trace()

    gram = sage_matrix(SageZZ, [[int(form(x, y)) for y in PAIRS] for x in PAIRS])

    def right(d):
        return sage_matrix(SageZZ, [_pair_coordinates((x[int(0)] * d, x[int(1)] * d)) for x in PAIRS]).transpose()

    return gram, [right(I), right(int(2) + J)]


def _is_transporter(t, gram, target_gram, sources, targets):
    return (
        t.det() in (int(1), int(-1))
        and t.transpose() * target_gram * t == gram
        and all(t * F == T * t for F, T in zip(sources, targets))
    )


def test_levi_transporter_to_the_inverse_on_one_summand():
    gram, E = _eichler()
    big = gram.block_sum(gram)
    F = E.block_sum(E)
    commutant = IsometryCommutant(big, [F])
    target = E.block_sum(E.inverse().change_ring(SageZZ))
    t = equivariant_transporter(commutant, big, [target])
    assert t is not None and _is_transporter(t, big, big, [F], [target])
    assert equivariant_transporter(commutant, int(3) * big, [F]) is None


def test_nonmaximal_levi_transporter():
    gram, E = _eichler()
    double = sage_matrix(SageZZ, int(6), int(6), int(0))
    pair = sage_matrix(SageZZ, int(6), int(6), int(0))
    for offset in (int(0), int(3)):
        double.set_block(offset, offset, gram)
        pair.set_block(offset, offset, E)
    basis = sage_matrix(SageZZ, int(6), int(6), int(0))
    for column, (row, scale) in enumerate([(int(0), int(1)), (int(3), int(1)), (int(2), int(1)), (int(5), int(1)), (int(1), int(1)), (int(4), int(2))]):
        basis[row, column] = scale
    F = basis.inverse() * pair * basis
    lattice_gram = basis.transpose() * double * basis
    assert F.denominator() == int(1) and lattice_gram.denominator() == int(1)
    F, lattice_gram = sage_matrix(SageZZ, F), sage_matrix(SageZZ, lattice_gram)
    commutant = IsometryCommutant(lattice_gram, [F])
    witness = sage_identity(SageZZ, int(6)).__copy__()
    witness[int(0), int(5)] = int(1)
    target = witness * F * witness.inverse()
    target_gram = witness.inverse().transpose() * lattice_gram * witness.inverse()
    t = equivariant_transporter(commutant, target_gram, [target])
    assert t is not None and _is_transporter(t, lattice_gram, target_gram, [F], [target])
    assert equivariant_transporter(commutant, int(3) * target_gram, [target]) is None


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


def test_definite_hermitian_quaternion_transporters():
    gram, sources, shear = _definite_hermitian_plane()
    assert all(shear * F == F * shear for F in sources) and shear.det() in (int(1), int(-1))
    commutant = IsometryCommutant(gram, sources)
    target = (shear.inverse().transpose() * gram * shear.inverse()).change_ring(SageZZ)
    t = equivariant_transporter(commutant, target, sources)
    assert t is not None and _is_transporter(t, gram, target, sources, sources)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None


BASIS = (ONE, I, J, K)


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


def test_skew_hermitian_quaternion_transporters():
    gram, sources = _skew_hermitian_quaternion_lattice()
    commutant = IsometryCommutant(gram, sources)
    P = _diagonal_left(((ONE, I), (int(0) * ONE, ONE)))
    assert P.det() in (int(1), int(-1)) and all(P * F == F * P for F in sources)
    moved = (P.inverse().transpose() * gram * P.inverse()).change_ring(SageZZ)
    assert moved != gram
    t = equivariant_transporter(commutant, moved, sources)
    assert t is not None and _is_transporter(t, gram, moved, sources, sources)
    assert equivariant_transporter(commutant, int(3) * gram, sources) is None
