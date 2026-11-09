r"""Generators of the Euclidean Picard groups ``PU(2,1; O_d)``, ``d = 2, 7, 11`` (Zhao, arXiv:1006.3331).

Theorems 5.2-5.4 of the source give five matrices generating the projective
group. The kernel of ``U(2,1; O_d) -> PU(2,1; O_d)`` is the scalar units
``+-1``, so the five lifts and ``-I`` generate ``U(2,1; O_d)``. The form is
``z_1 conj(w_3) + z_2 conj(w_2) + z_3 conj(w_1)``.
"""

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.rational_field import QQ
from sage_extensions.categories.hermitian_lattices import HermitianLattice


def field(d: int):
    r"""Return ``Q(sqrt(-d))`` for a source value of ``d``."""
    if d not in (2, 7, 11):
        raise ValueError(f"the source covers d = 2, 7, 11, not {d}")
    return QuadraticField(-d, "a")


def form(d: int):
    r"""Return the source Hermitian Gram matrix."""
    E = field(d)
    return matrix(E, 3, 3, [0, 0, 1, 0, 1, 0, 1, 0, 0])


def order_basis(d: int) -> tuple:
    r"""Return a ``ZZ``-basis of the maximal order."""
    E = field(d)
    return (E.one(), E.gen() if d == 2 else (1 + E.gen()) / 2)


def generators(d: int) -> dict:
    r"""Return the five source matrices of Theorems 5.2-5.4."""
    E = field(d)
    a = E.gen()
    omega = (1 + a) / 2
    bar = omega.conjugate()

    def entries(*rows):
        return matrix(E, 3, 3, [entry for row in rows for entry in row])

    I0 = entries((0, 0, 1), (0, -1, 0), (1, 0, 0))
    R1 = entries((1, 0, 0), (0, -1, 0), (0, 0, 1))
    T = entries((1, 0, a), (0, 1, 0), (0, 0, 1))
    match d:
        case 2:
            R2 = entries((1, 2, -2), (0, -1, 2), (0, 0, 1))
            R3 = entries((1, -a, -1), (0, -1, a), (0, 0, 1))
        case 7:
            R2 = entries((1, 1, -bar), (0, -1, 1), (0, 0, 1))
            R3 = entries((1, bar, -1), (0, -1, omega), (0, 0, 1))
        case 11:
            R2 = entries((1, 1, -bar), (0, -1, 1), (0, 0, 1))
            R3 = entries((1, bar, -1 - bar), (0, -1, omega), (0, 0, 1))
    return {"I0": I0, "R1": R1, "R2": R2, "R3": R3, "T": T}


def is_element(d: int, g) -> bool:
    r"""Return whether ``g`` is an integral unitary matrix of the source form."""
    order = field(d).ring_of_integers()
    return all(entry in order for entry in g.list()) and g.conjugate_transpose() * form(d) * g == form(d)


def full_unitary_generators(d: int) -> tuple:
    r"""Return the five lifts with the scalar kernel generator ``-I``."""
    E = field(d)
    return (*generators(d).values(), -matrix.identity(E, 3))


def evaluate(d: int, word) -> "matrix":
    r"""Return the product of a word of ``(name, exponent)`` pairs; ``S`` is ``-I``."""
    E = field(d)
    table = dict(generators(d), S=-matrix.identity(E, 3))
    result = matrix.identity(E, 3)
    for name, exponent in word:
        result = result * table[name] ** int(exponent)
    return result


def _shears(d: int):
    r"""Return ``R2 R1`` and ``R3 R1``, the unipotent cusp elements spanning the ``z`` lattice."""
    g = generators(d)
    return g["R2"] * g["R1"], g["R3"] * g["R1"]


def _shear_coordinates(d: int, value):
    r"""Return rational ``(m, n)`` with ``value = m z_A + n z_B``."""
    A, B = _shears(d)
    basis = matrix(QQ, [list(A[1, 2].vector()), list(B[1, 2].vector())])
    return tuple(vector(QQ, list(value.vector())) * ~basis)


def _cusp_word(d: int, element):
    r"""Return a word for an upper triangular element of ``U(2,1; O_d)``."""
    E = field(d)
    sign = element[0, 0]
    if sign not in (E.one(), -E.one()) or element[2, 2] != sign:
        raise ArithmeticError(f"{element} is not in the source cusp stabilizer")
    scalar = [] if sign == E.one() else [("S", 1)]
    current = sign * element
    rotation = []
    if current[1, 1] == -E.one():
        current = current * generators(d)["R1"]
        rotation = [("R1", 1)]
    m, n = _shear_coordinates(d, current[1, 2])
    if m.denominator() != 1 or n.denominator() != 1:
        raise ArithmeticError(f"{current} has z outside the lattice of the source generators")
    A, B = _shears(d)
    central = ~(A ** int(m) * B ** int(n)) * current
    k = central[0, 2] / E.gen()
    if k not in QQ or QQ(k).denominator() != 1:
        raise ArithmeticError(f"{central} is not a power of the vertical translation")
    return [*scalar, *_shear_power("R2", int(m)), *_shear_power("R3", int(n)), ("T", int(k)), *rotation]


def _shear_power(name: str, exponent: int):
    r"""Return the word of ``(name R1)^exponent``."""
    if exponent >= 0:
        return [(name, 1), ("R1", 1)] * exponent
    return [("R1", -1), (name, -1)] * (-exponent)


def _reduction(d: int, g):
    r"""Return ``(M, k)`` with ``|((T^k M g)_{11} / g_{31})| < 1`` and minimal among the window."""
    E = field(d)
    A, B = _shears(d)
    ratio = g[1, 0] / g[2, 0]
    m0, n0 = _shear_coordinates(d, -ratio)
    best = None
    for dm in range(-2, 3):
        for dn in range(-2, 3):
            M = A ** (m0.round() + dm) * B ** (n0.round() + dn)
            value = (M * g)[0, 0] / g[2, 0]
            coefficient = value.vector()[1] / E.gen().vector()[1]
            for k in (-coefficient.floor() - 1, -coefficient.floor(), -coefficient.floor() + 1):
                norm = (value + int(k) * E.gen()).norm()
                if best is None or norm < best[0]:
                    best = (norm, (m0.round() + dm, n0.round() + dn), int(k))
    return best


def word(d: int, g):
    r"""Return a checked word for ``g`` in the five source matrices and ``-I``.

    Each step applies ``I0 T^k N`` for a cusp translation ``N``. The new
    ``|g_31|^2`` is ``|(T^k N g)_{11}|^2``, a smaller positive integer, so the
    reduction stops at an upper triangular element.
    """
    if not is_element(d, g):
        raise ValueError(f"{g} is not an element of U(2,1; O_{d})")
    table = generators(d)
    A, B = _shears(d)
    prefix = []
    current = g
    while current[2, 0] != 0:
        norm, (m, n), k = _reduction(d, current)
        if norm >= 1:
            raise ArithmeticError(f"the source covering failed to reduce {current}")
        M = A ** int(m) * B ** int(n)
        following = table["I0"] * table["T"] ** k * M * current
        # g = M^-1 T^-k I0^-1 following, and I0 is an involution.
        prefix += [*_shear_power("R3", -int(n)), *_shear_power("R2", -int(m)), ("T", -k), ("I0", 1)]
        current = following
    result = prefix + _cusp_word(d, current)
    if evaluate(d, result) != g:
        raise ArithmeticError("the word does not evaluate to the given element")
    return tuple(result)


def lattice(d: int) -> HermitianLattice:
    r"""Return ``O_d^3`` with the source form."""
    return HermitianLattice(order_basis(d), form(d).transpose())
