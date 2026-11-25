r"""Generators and words for `U(H \oplus I_m)` over `E = Z[\omega]`, `m \le 2` (#51).

Coordinates are `(x; \mu_1, \dots, \mu_m; y)` with `h(v, w) = w^* J v` and
`J = [[0, 0, 1], [0, I_m, 0], [1, 0, 0]]`, so `e = (1; 0; 0)` and `f = (0; 0; 1)` are null with
`h(f, e) = 1`, and `h(v, e) = y`. For `\lambda \in E^m` and `c \in E` with `c + \bar c = -|\lambda|^2`,
`N(\lambda, c)` fixes `e`, sends `b_i \mapsto b_i - \bar\lambda_i e` and `f \mapsto c e + \lambda + f`. Then
`N v = (x - \lambda^* \mu + c y; \mu + y \lambda; y)`. The imaginary elements of `E` are `\theta Z` with
`\theta = \sqrt{-3}`, and `\theta / 2 - 1/2 = \omega`, so `c` exists for every `\lambda`.

Reduction of a primitive null `v` with `y \ne 0`. Choose `\lambda` with `|\mu_i / y + \lambda_i|^2 \le 1/3`,
the covering radius of `E`. Then `\mathrm{Re}(x'/y) = -|\mu'|^2 / (2 |y|^2) \in [-m/6, 0]`, and a choice of
`c` in its coset of `\theta Z` gives `|\mathrm{Im}(x'/y)| \le \sqrt 3 / 2`. So
`|x'/y|^2 \le m^2/36 + 3/4 < 1` for `m \le 2`, and the swap `W` of `e` and `f` lowers `|y|`.
At `y = 0`, `h(\mu, \mu) = 0` forces `\mu = 0`, so `v` is a unit multiple of `e`. An element
fixing `e` is `N(\lambda, c)` times an element fixing `e` and `f`, which is a unitary map of `I_m`:
a monomial matrix in the units. So the translations, `W`, the units on `e` and the
monomial group generate `U(H \oplus I_m)`.
"""

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.number_field.number_field import CyclotomicField

FIELD = CyclotomicField(3)
OMEGA = FIELD.gen()
THETA = OMEGA - OMEGA.conjugate()
ONE, ZERO = FIELD.one(), FIELD.zero()
UNIT = -OMEGA


def form(m: int):
    size = m + 2
    J = matrix(FIELD, size, size)
    J[0, size - 1] = J[size - 1, 0] = ONE
    for i in range(1, m + 1):
        J[i, i] = ONE
    return J


def _check_rank(m: int) -> None:
    if m not in (0, 1, 2):
        raise ValueError(f"the height reduction decreases |y| only for m <= 2, not m = {m}")


def is_element(m: int, g) -> bool:
    integral = all(c.denominator() == 1 for entry in g.list() for c in entry.vector())
    return integral and g.conjugate_transpose() * form(m) * g == form(m)


def translation(m: int, lam, c):
    lam = vector(FIELD, lam)
    if c + c.conjugate() != -lam.dot_product(lam.conjugate()):
        raise ValueError(f"c = {c} does not satisfy c + conj(c) = -|lambda|^2")
    size = m + 2
    g = matrix.identity(FIELD, size)
    for i in range(m):
        g[0, i + 1] = -lam[i].conjugate()
        g[i + 1, size - 1] = lam[i]
    g[0, size - 1] = c
    return g


def _base_c(lam):
    square = lam.dot_product(lam.conjugate())
    return -square / 2 + (THETA / 2 if int(square.vector()[0]) % 2 else ZERO)


def swap(m: int):
    size = m + 2
    g = matrix.identity(FIELD, size)
    g[0, 0] = g[size - 1, size - 1] = ZERO
    g[0, size - 1] = g[size - 1, 0] = ONE
    return g


def unit_on_e(m: int, u):
    size = m + 2
    g = matrix.identity(FIELD, size)
    g[0, 0] = u
    g[size - 1, size - 1] = ~u.conjugate()
    return g


def monomial(m: int, permutation, units):
    size = m + 2
    g = matrix.identity(FIELD, size)
    for i in range(m):
        g[i + 1, i + 1] = ZERO
    for i in range(m):
        g[permutation[i] + 1, i + 1] = units[i]
    return g


def generators(m: int) -> dict:
    r"""Return the named generators of `U(H \oplus I_m)`."""
    _check_rank(m)
    table = {"W": swap(m), "D": unit_on_e(m, UNIT), "Z": translation(m, [ZERO] * m, THETA)}
    for i in range(m):
        for name, scalar in (("", ONE), ("w", OMEGA)):
            lam = vector(FIELD, [scalar if j == i else ZERO for j in range(m)])
            table[f"T{name}{i}"] = translation(m, lam, _base_c(lam))
    if m:
        table["S"] = monomial(m, tuple(range(m)), (UNIT, *([ONE] * (m - 1))))
    if m == 2:
        table["P"] = monomial(m, (1, 0), (ONE, ONE))
    return table


def _nearest(value):
    a, b = value.vector()
    candidates = [FIELD(a.round() + i) + FIELD(b.round() + j) * OMEGA for i in (-1, 0, 1) for j in (-1, 0, 1)]
    return min(candidates, key=lambda t: (value - t).norm())


def _translation_word(m: int, lam, c):
    r"""Return a word for `N(\lambda, c)` in the basis translations and `Z`."""
    table = generators(m)
    word = []
    built = matrix.identity(FIELD, m + 2)
    for i in range(m):
        a, b = lam[i].vector()
        for name, power in ((f"T{i}", int(a)), (f"Tw{i}", int(b))):
            if power:
                word.append((name, power))
                built = built * table[name] ** power
    remainder = (~built * translation(m, lam, c))[0, m + 1]
    multiple = remainder / THETA
    if not multiple.is_rational() or multiple.denominator() != 1:
        raise ArithmeticError(f"the central part {remainder} is not in theta Z")
    if multiple:
        word.append(("Z", int(multiple)))
    return word


def evaluate(m: int, word):
    table = generators(m)
    result = matrix.identity(FIELD, m + 2)
    for name, power in word:
        result = result * table[name] ** int(power)
    return result


def _inverse(word):
    return [(name, -power) for name, power in reversed(word)]


def _unit_power(u):
    for power in range(6):
        if UNIT**power == u:
            return power
    raise ArithmeticError(f"{u} is not a unit of Z[omega]")


def _monomial_word(m: int, g):
    r"""Return a word for a monomial `g` fixing `e` and `f`: `P^a S^{p_1} (P S^{p_2} P)`."""
    word = []
    current = g
    if m == 2 and current[1, 1] == ZERO:
        word.append(("P", 1))
        current = swap_block(m) * current
    word.append(("S", _unit_power(current[1, 1])))
    if m == 2:
        word.extend([("P", 1), ("S", _unit_power(current[2, 2])), ("P", 1)])
    return word


def swap_block(m: int):
    return monomial(m, (1, 0), (ONE, ONE))


def word(m: int, g):
    r"""Return a checked word in :func:`generators` equal to ``g``."""
    _check_rank(m)
    if not is_element(m, g):
        raise ValueError(f"{g} is not in U(H + I_{m}) over Z[omega]")
    table = generators(m)
    size = m + 2
    reducer = []
    current = g
    while True:
        v = current.column(0)
        y = v[size - 1]
        if y == ZERO:
            break
        mu = vector(FIELD, [v[i + 1] for i in range(m)])
        lam = vector(FIELD, [-_nearest(mu[i] / y) for i in range(m)])
        base = _base_c(lam)
        step = translation(m, lam, base)
        moved = step * v
        q = moved[0] / y
        shift = (-(q - q.conjugate()) / (2 * THETA)).vector()[0].round()
        c = base + THETA * FIELD(shift)
        step = translation(m, lam, c)
        moved = step * v
        if moved[0].norm() >= y.norm():
            raise ArithmeticError("the translation did not make |x| < |y|")
        prefix = [("W", 1), *_translation_word(m, lam, c)]
        reducer = prefix + reducer
        current = swap(m) * step * current
    unit = current[0, 0]
    for power in range(6):
        if UNIT**power * unit == ONE:
            break
    else:
        raise ArithmeticError(f"{unit} is not a unit")
    reducer = [("D", power), *reducer]
    current = table["D"] ** power * current
    f = current.column(size - 1)
    lam = vector(FIELD, [f[i + 1] for i in range(m)])
    parabolic = _translation_word(m, lam, f[0])
    levi = ~translation(m, lam, f[0]) * current
    result = _inverse(reducer) + parabolic + (_monomial_word(m, levi) if m else [])
    if evaluate(m, result) != g:
        raise ArithmeticError("the word does not evaluate to g")
    return result


def lattice(m: int):
    r"""Return `E^{m+2}` with the form `J`, whose trace lattice has rank `2m + 4`."""
    from sage_extensions.categories.hermitian_lattices import HermitianLattice

    return HermitianLattice((ONE, OMEGA), form(m).transpose())


def trace_generators(m: int) -> tuple:
    r"""Return :func:`generators` as integral isometries of the trace lattice commuting with `\omega`."""
    hermitian = lattice(m)
    return tuple(hermitian.isometry_from_matrix(g) for g in generators(m).values())
