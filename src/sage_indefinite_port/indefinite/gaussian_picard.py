r"""The Gaussian Picard group ``SU(2,1; Z[i])`` (Falbel-Francsics-Lax-Parker, arXiv:0911.1104).

Theorem 1 of the source: four matrices generate the group. Its proof reduces
``|g_13|^2`` by at least half with ``J N``, and writes the cusp stabilizer
through ``N_{(0,i)}``, ``N_{(1+i,1)}`` and ``M_i`` (Proposition 3.2).
Elements are matrices over ``Q(i)`` acting on the form ``C = -E_13 - E_31 + E_22``.
"""

from sage.matrix.constructor import matrix
from sage.rings.number_field.number_field import CyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice

FIELD = CyclotomicField(4)
I = FIELD.gen()
ONE = FIELD.one()
ZERO = FIELD.zero()


def _matrix(rows):
    return matrix(FIELD, 3, 3, [FIELD(entry) for row in rows for entry in row])


FORM = _matrix(((ZERO, ZERO, -ONE), (ZERO, ONE, ZERO), (-ONE, ZERO, ZERO)))
GENERATORS = {
    "N0": _matrix(((ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (I, ZERO, ONE))),
    "N1": _matrix(((ONE, ZERO, ZERO), (ONE + I, ONE, ZERO), (ONE, ONE - I, ONE))),
    "M": _matrix(((I, ZERO, ZERO), (ZERO, -ONE, ZERO), (ZERO, ZERO, I))),
    "J": _matrix(((ZERO, ZERO, -ONE), (ZERO, -ONE, ZERO), (-ONE, ZERO, ZERO))),
}
SCALAR = I * matrix.identity(FIELD, 3)
CUSP_STABILIZER_GENERATORS = tuple(GENERATORS[name] for name in ("N0", "N1", "M"))


def _parts(value):
    real, imaginary = value.vector()
    return real, imaginary


def is_element(g) -> bool:
    r"""Return whether ``g`` lies in ``SU(2,1; Z[i])``."""
    integral = all(part.denominator() == 1 for entry in g.list() for part in _parts(entry))
    return integral and g.det() == ONE and g.conjugate_transpose() * FORM * g == FORM


def evaluate(word):
    r"""Return the product of a word of ``(generator name, exponent)`` pairs."""
    result = matrix.identity(FIELD, 3)
    for name, exponent in word:
        result = result * GENERATORS[name] ** int(exponent)
    return result


def _translation(gamma, height):
    return _matrix(((ONE, ZERO, ZERO), (gamma, ONE, ZERO), (height, gamma.conjugate(), ONE)))


def _stabilizer_word(p):
    r"""Proposition 3.2: ``P = N_{(0,i)}^{r-2kl} N_1^k M N_1^l M^{-1} M_beta``."""
    beta = p[0, 0]
    exponent = next(e for e in range(4) if I**e == beta)
    translation = p * ~(GENERATORS["M"] ** exponent)
    gamma = translation[1, 0]
    m, n = _parts(gamma)
    r = _parts(translation[2, 0])[1]
    k, l = (m + n) / 2, (n - m) / 2
    return [("N0", r - 2 * k * l), ("N1", k), ("M", 1), ("N1", l), ("M", -1), ("M", exponent)]


def _reduction_step(g):
    r"""Return ``(gamma, height)`` for equations (e60)-(e65) of the source."""
    ratio = g[1, 2] / g[0, 2]
    a, b = _parts(-ratio)
    x, y = (a + b) / 2, (b - a) / 2
    m, n = x.round(), y.round()
    gamma = FIELD(m - n) + FIELD(m + n) * I
    shift = _parts(gamma.conjugate() * ratio)[1] + _parts(g[2, 2] / g[0, 2])[1]
    k = (-shift).round()
    return gamma, FIELD(m * m + n * n) + FIELD(k) * I


def word(g):
    r"""Return a checked word in the four generators representing ``g``."""
    if not is_element(g):
        raise ValueError(f"{g} is not an element of SU(2,1; Z[i])")
    prefix = []
    current = g
    while current[0, 2] != ZERO:
        gamma, height = _reduction_step(current)
        translation = _translation(gamma, height)
        following = GENERATORS["J"] * translation * current
        if 2 * following[0, 2].norm() > current[0, 2].norm():
            raise ArithmeticError("the source reduction did not halve |g_13|^2")
        prefix += [*_stabilizer_word(~translation), ("J", 1)]
        current = following
    result = prefix + _stabilizer_word(current)
    if evaluate(result) != g:
        raise ArithmeticError("the word does not evaluate to the given element")
    return tuple(result)


def lattice() -> HermitianLattice:
    r"""Return ``Z[i]^3`` with the source form."""
    return HermitianLattice((ONE, I), FORM.transpose())


def trace_isometry(g):
    r"""Return ``g`` as an isometry of the rank-six trace lattice."""
    hermitian = lattice()
    trace = hermitian.trace_lattice()
    positions = dict(zip(hermitian._labels(), trace.module_generating_set(), strict=True))
    basis = hermitian.order_basis()

    def image(name):
        i, a = next(key for key, value in positions.items() if value == name)
        result = {}
        for row in range(3):
            coordinates = hermitian._coordinates_in_order(g[row, i] * basis[a])
            for b, value in enumerate(coordinates):
                if value:
                    result[positions[(row, b)]] = trace.base_ring()(value)
        return trace.linear_combination(result)

    return trace.Aut()(image)
