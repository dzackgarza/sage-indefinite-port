r"""The Eisenstein lattice ``I_{4,1}`` of Allcock-Carlson-Toledo (math/0007048, section 7).

Column vectors carry ``<v|w> = w^* A v`` with ``A`` of (7.7.1). ``R_1..R_7`` are
the ``(-omega)``-reflections in the roots of (7.7.2), ``B`` is the biflection
of (7.15), ``T_{lambda,z}`` are the translations of (7.8). Theorem 7.21:
``Gamma = <R_i> = Aut^+`` acts transitively on primitive null vectors.
"""

from itertools import product

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.number_field.number_field import CyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice

FIELD = CyclotomicField(3)
OMEGA = FIELD.gen()
BAR = OMEGA.conjugate()
THETA = OMEGA - BAR
ONE, ZERO = FIELD.one(), FIELD.zero()
RHO = vector(FIELD, [ZERO, ZERO, ZERO, ZERO, ONE])

_A = matrix(FIELD, 5, 5)
for _i in range(3):
    _A[_i, _i] = ONE
_A[3, 4] = _A[4, 3] = ONE
FORM = _A


def _root(*coordinates):
    return vector(FIELD, [FIELD(c) for c in coordinates])


ROOTS = {
    1: _root(1, 0, 0, 0, 0),
    2: _root(1, 0, 0, 0, 1),
    3: _root(0, 0, 0, 1, -OMEGA),
    4: _root(0, 1, 0, 0, 1),
    5: _root(0, 1, 0, 0, 0),
    6: _root(0, 0, 1, 0, 1),
    7: _root(0, 0, 1, 0, 0),
}


def h(v, w):
    r"""Return ``<v|w> = w^* A v``."""
    return w.conjugate() * FORM * v


def reflection(root, xi):
    r"""Return ``x -> x - (1 - xi) <x|r> r / <r|r>`` as a matrix."""
    norm = h(root, root)
    images = [
        e - (ONE - xi) * h(e, root) / norm * root
        for e in matrix.identity(FIELD, 5).columns()
    ]
    return matrix(FIELD, images).transpose()


MINUS_OMEGA = -OMEGA
R = {index: reflection(root, MINUS_OMEGA) for index, root in ROOTS.items()}
BIFLECTION = reflection(_root(1, -1, 0, 0, 0), -ONE)
AUT_PLUS_GENERATORS = tuple(R[index] for index in (1, 2, 3, 4, 5, 6, 7))
FULL_GENERATORS = (*AUT_PLUS_GENERATORS, BIFLECTION)


def translation(lam, z):
    r"""Return ``T_{lambda, z}`` of (7.8), with ``lambda`` in ``K^3`` and ``z`` imaginary."""
    lam = vector(FIELD, lam)
    star = lam.conjugate()
    square = lam.dot_product(star)
    result = matrix.identity(FIELD, 5)
    for i in range(3):
        result[i, 3] = lam[i]
        result[4, i] = -star[i]
    result[4, 3] = z - square / 2
    return result


def is_isometry(g) -> bool:
    r"""Return whether ``g`` is unitary with entries in ``Z[omega]``."""
    integral = all(c.denominator() == 1 for entry in g.list() for c in entry.vector())
    return integral and g.conjugate_transpose() * FORM * g == FORM


def is_primitive(v) -> bool:
    r"""Return whether the coordinates generate the unit ideal of ``Z[omega]``."""
    return FIELD.ideal(list(v)).norm() == 1


def _nearest(value):
    r"""Return an element of ``Z[omega]`` within the covering radius of ``value``."""
    a, b = value.vector()
    candidates = [
        FIELD(a.round() + i) + FIELD(b.round() + j) * OMEGA
        for i, j in product((-1, 0, 1), repeat=2)
    ]
    return min(candidates, key=lambda c: (value - c).norm())


def _admissible_translation(lam):
    square = lam.dot_product(lam.conjugate())
    return translation(lam, ZERO if int(square.vector()[0]) % 2 == 0 else THETA / 2)


def _height_root(v):
    r"""Return ``r = (0,0,0;1,n theta - omega)`` that lowers the height or is orthogonal to ``v``."""
    height = v[3]
    best = None
    ratio = v[4] / height
    centre = (ONE + BAR * ratio - BAR * BAR) / (BAR * THETA)
    start = centre.vector()[0].round()
    for n in range(int(start) - 3, int(start) + 4):
        root = _root(0, 0, 0, 1, FIELD(n) * THETA - OMEGA)
        image = reflection(root, MINUS_OMEGA) * v
        if best is None or image[3].norm() < best[0]:
            best = (image[3].norm(), root)
        if h(v, root) == ZERO:
            return root, True
    return best[1], False


def reduce_null_vector(v):
    r"""Return ``(g, unit)`` with ``g`` in ``Gamma`` and ``g v = unit rho`` (7.14, 7.22)."""
    v = vector(FIELD, v)
    if h(v, v) != ZERO or not is_primitive(v):
        raise ValueError(f"{v} is not a primitive null vector")
    g = matrix.identity(FIELD, 5)
    while v[3] != ZERO:
        height = v[3]
        lam = vector(FIELD, [-_nearest(v[i] / height) for i in range(3)])
        step = _admissible_translation(lam)
        v = step * v
        g = step * g
        root, orthogonal = _height_root(v)
        if orthogonal:
            transport = _move_root_to_r7(root)
            v = transport * v
            g = transport * g
            if v[2] != ZERO:
                raise ArithmeticError("the transported vector is not orthogonal to r_7")
            root, orthogonal = _height_root(v)
            if orthogonal:
                raise ArithmeticError("a vector of I_{3,1} is orthogonal to a height-one root")
        step = reflection(root, MINUS_OMEGA)
        next_vector = step * v
        if next_vector[3].norm() >= height.norm():
            raise ArithmeticError("the source reduction did not lower the height")
        v = next_vector
        g = step * g
    unit = v[4]
    if not unit.is_unit() or any(v[i] != ZERO for i in range(4)):
        raise ArithmeticError("the reduction did not end at a multiple of rho")
    return g, unit


def _move_root_to_r7(root):
    r"""Return an element of ``Gamma`` carrying the height-one ``root`` to a multiple of ``r_7``.

    A translation carries it to ``r_3``; braiding ``R_3, R_6`` and ``R_6, R_7`` then moves the root.
    """
    n = (root[4] + OMEGA) / THETA
    shift = translation(_root(0, 0, 0), -FIELD(n) * THETA)
    moved = shift * root
    if moved != ROOTS[3] and moved / moved[3] != ROOTS[3]:
        raise ArithmeticError("the translation did not reach r_3")
    return (R[6] * R[7]) * (R[3] * R[6]) * shift


def lattice() -> HermitianLattice:
    r"""Return ``Z[omega]^5`` with the source form; its trace lattice has rank 10."""
    return HermitianLattice((ONE, OMEGA), FORM.transpose())
