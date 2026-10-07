"""Exact rank-two Lorentzian arithmetic for the Allcock edgewalk.

Translation of two_dim_lorentzian.h from polyhedral_common at
a55fcb7b71af48c88d7abbeab7889e9347916e43.
"""

from collections.abc import Sequence
from math import isqrt

from sage.arith.misc import divisors, gcd, xgcd
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ

type _Gram = Matrix_integer_dense | Matrix_rational_dense
type _Scalar = int | Integer | Rational
type _IntegralVector = FreeModuleElement[Integer]
type _IntegralRow = Sequence[int | Integer] | _IntegralVector
type _RationalRow = Sequence[_Scalar] | FreeModuleElement[Rational] | _IntegralVector
type _Pair = tuple[_IntegralVector, _IntegralVector]


def _gram(G: _Gram) -> Matrix_rational_dense:
    rational = matrix(QQ, G)
    if rational.dimensions() != (2, 2) or not rational.is_symmetric() or rational.det() >= 0:
        raise ValueError("G must be a symmetric Lorentzian 2 by 2 matrix")
    return rational


def _zvector(v: _IntegralRow) -> _IntegralVector:
    v = vector(ZZ, v)
    if len(v) != 2:
        raise ValueError("rank-two vectors must have length 2")
    return v


def quadratic_eval(G: _Gram, v: _IntegralRow) -> Rational:
    G, v = _gram(G), _zvector(v)
    return G[0, 0] * v[0] ** 2 + 2 * G[0, 1] * v[0] * v[1] + G[1, 1] * v[1] ** 2


def scalar_eval(G: _Gram, v: _IntegralRow, w: _IntegralRow) -> Rational:
    G, v, w = _gram(G), _zvector(v), _zvector(w)
    return G[0, 0] * v[0] * w[0] + G[0, 1] * (v[0] * w[1] + v[1] * w[0]) + G[1, 1] * v[1] * w[1]


def oriented_determinant(r: _IntegralRow, l: _IntegralRow) -> Integer:
    r, l = _zvector(r), _zvector(l)
    return r[0] * l[1] - r[1] * l[0]


def _primitive_integral(v: _RationalRow) -> _IntegralVector:
    v = vector(QQ, v)
    den = ZZ.one()
    for x in v:
        den = den.lcm(x.denominator())
    z = vector(ZZ, [den * x for x in v])
    c = gcd(z)
    if c == 0:
        raise ValueError("zero vector has no primitive reduction")
    return vector(ZZ, [x // c for x in z])


def isotropic_factorization(G: _Gram) -> Matrix_rational_dense | None:
    """Rows F_i satisfy Q(x,y)=(F_0.(x,y))(F_1.(x,y))."""
    G = _gram(G)
    a, b, c = G[0, 0], G[0, 1], G[1, 1]
    delta = b * b - a * c
    if not delta.is_square():
        return None
    root = delta.sqrt()
    if a != 0:
        x1, x2 = (-b + root) / a, (-b - root) / a
        return matrix(QQ, [[a, -a * x1], [1, -x2]])
    if c != 0:
        y1, y2 = (-b + root) / c, (-b - root) / c
        return matrix(QQ, [[-c * y1, c], [-y2, 1]])
    return matrix(QQ, [[2 * b, 0], [0, 1]])


def primitive_isotropic_vectors(G: _Gram) -> tuple[_IntegralVector, ...]:
    F = isotropic_factorization(G)
    if F is None:
        raise ValueError("the form is anisotropic over QQ")
    return tuple(_primitive_integral([-row[1], row[0]]) for row in F.rows())


def fixed_norm_vectors_isotropic(G: _Gram, norm: _Scalar) -> tuple[_IntegralVector, ...]:
    F = isotropic_factorization(G)
    if F is None:
        raise ValueError("the form is anisotropic over QQ")
    norm = QQ(norm)
    if norm == 0:
        u, v = primitive_isotropic_vectors(G)
        return (u, -u, v, -v)
    rows: list[_IntegralVector] = []
    mult: list[Rational] = []
    for row in F.rows():
        p = _primitive_integral(row)
        rows.append(p)
        mult.append(next(row[i] / p[i] for i in range(2) if p[i]))
    scaled = norm / (mult[0] * mult[1])
    if scaled.denominator() != 1:
        return ()
    integral_scaled = ZZ(scaled)
    Ainv: Matrix_rational_dense = matrix(QQ, rows).inverse()
    out: set[tuple[Integer, ...]] = set()
    for d in divisors(abs(integral_scaled)):
        for first in (ZZ(d), -ZZ(d)):
            sol = Ainv * vector(QQ, [first, integral_scaled / first])
            if all(x.denominator() == 1 for x in sol):
                out.add(tuple(ZZ(x) for x in sol))
    return tuple(vector(ZZ, x) for x in sorted(out))


def _floor_sqrt(value: _Scalar) -> Integer:
    value = QQ(value)
    if value < 0:
        raise ValueError("negative radicand")
    x = ZZ(isqrt(ZZ(value.numerator()) // ZZ(value.denominator())))
    while QQ((x + 1) ** 2) <= value:
        x += 1
    return x


def canonical_companion(G: _Gram, bound: _Scalar, r: _IntegralRow, l: _IntegralRow) -> _IntegralVector:
    rr, rl, ll = quadratic_eval(G, r), scalar_eval(G, r, l), quadratic_eval(G, l)
    if rr <= 0:
        raise ValueError("r must have positive norm")
    root = _floor_sqrt(rl * rl - rr * (ll - QQ(bound)))
    k = ZZ(((-rl + root) / rr).floor())
    return _zvector(l) + k * _zvector(r)


def promised_step(G: _Gram, bound: _Scalar, r: _IntegralRow, l: _IntegralRow) -> _Pair:
    G, M, r, l = _gram(G), QQ(bound), _zvector(r), _zvector(l)
    if quadratic_eval(G, r) <= 0 or oriented_determinant(r, l) != 1:
        raise ValueError("Promised requires Q(r)>0 and det(r,l)=1")
    if quadratic_eval(G, l) > M:
        raise ValueError("Promised requires Q(l)<=M")
    while True:
        m = r + l
        if quadratic_eval(G, m) <= M or scalar_eval(G, m, r) < 0:
            l = m
            continue
        if quadratic_eval(G, l) >= 0 and scalar_eval(G, r, l) > 0:
            return l, -r
        r = m


def shorter_pair(G: _Gram, r: _IntegralRow, l: _IntegralRow) -> _Pair | None:
    G, r, l = _gram(G), _zvector(r), _zvector(l)
    M = quadratic_eval(G, r)
    char0 = (quadratic_eval(G, r), scalar_eval(G, r, l), quadratic_eval(G, l))
    while True:
        current = quadratic_eval(G, r)
        r, l = promised_step(G, current, r, l)
        current = quadratic_eval(G, r)
        l = canonical_companion(G, current, r, l)
        if current < M:
            return r, l
        if (quadratic_eval(G, r), scalar_eval(G, r, l), quadratic_eval(G, l)) == char0:
            return None


def reduced_start_pair(G: _Gram, bound: _Scalar, r: _IntegralRow, l: _IntegralRow) -> _Pair | None:
    G, M, r, l = _gram(G), QQ(bound), _zvector(r), _zvector(l)
    if quadratic_eval(G, r) <= M:
        r, l = promised_step(G, M, r, l)
        return r, canonical_companion(G, M, r, l)
    while True:
        pair = shorter_pair(G, r, l)
        if pair is None:
            return None
        r, l = pair
        if quadratic_eval(G, r) <= M:
            return r, canonical_companion(G, M, r, l)


def oriented_complement(r: _IntegralRow) -> _IntegralVector:
    r = _zvector(r)
    g, s, t = xgcd(r[0], r[1])
    if g == -1:
        g, s, t = -g, -s, -t
    if g != 1:
        raise ValueError("r must be primitive")
    return vector(ZZ, [-t, s])


def anisotropic_cycle(
    G: _Gram,
    bound: _Scalar,
    r: _IntegralRow,
) -> tuple[Matrix_integer_dense, tuple[_IntegralVector, ...]] | None:
    G = _gram(G)
    if isotropic_factorization(G) is not None:
        raise ValueError("requires a QQ-anisotropic form")
    r0 = _zvector(r)
    l0 = canonical_companion(G, bound, r0, oriented_complement(r0))
    start = reduced_start_pair(G, bound, r0, l0)
    if start is None:
        return None
    r1, l1 = start
    char0 = (quadratic_eval(G, r1), scalar_eval(G, r1, l1), quadratic_eval(G, l1))
    cycle: list[_IntegralVector] = []
    r, l = r1, l1
    while True:
        cycle.append(r)
        r, l = promised_step(G, bound, r, l)
        l = canonical_companion(G, bound, r, l)
        if (quadratic_eval(G, r), scalar_eval(G, r, l), quadratic_eval(G, l)) == char0:
            transform = matrix(QQ, [r1, l1]).inverse() * matrix(QQ, [r, l])
            if any(x.denominator() != 1 for x in transform.list()):
                raise ArithmeticError("nonintegral cycle")
            return matrix(ZZ, transform), tuple(cycle)
