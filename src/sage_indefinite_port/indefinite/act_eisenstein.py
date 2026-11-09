r"""The Eisenstein lattice ``I_{4,1}`` of Allcock-Carlson-Toledo (math/0007048, section 7).

Column vectors carry ``<v|w> = w^* A v`` with ``A`` of (7.7.1). ``R_1..R_7`` are
the ``(-omega)``-reflections in the roots of (7.7.2), ``B`` is the biflection
of (7.15), ``T_{lambda,z}`` are the translations of (7.8). Theorem 7.21:
``Gamma = <R_i> = Aut^+`` acts transitively on primitive null vectors.
"""

from itertools import count, product

from sage.arith.functions import lcm

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
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
    r"""Return ``g`` in ``Gamma`` with ``g v = rho`` (7.14, 7.16, 7.22)."""
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
    if any(v[i] != ZERO for i in range(4)) or v[4].norm() != ONE:
        raise ArithmeticError("the reduction did not end at a unit multiple of rho")
    step = _unit_transport(v[4])
    g = step * g
    if step * v != RHO:
        raise ArithmeticError("the unit transport did not reach rho")
    return g


def null_transporter(v, w):
    r"""Return ``g`` in ``Gamma`` with ``g v = w`` for primitive null vectors (Theorem 7.21)."""
    g = ~reduce_null_vector(w) * reduce_null_vector(v)
    if g * vector(FIELD, v) != vector(FIELD, w) or not is_isometry(g):
        raise ArithmeticError("the transporter does not carry the null vector")
    return g


def _null_vector_orthogonal_to(root):
    r"""Return a primitive null vector orthogonal to ``root``, by increasing coefficient height.

    ``root^perp`` is a Hermitian space of signature ``(3, 1)`` over ``Q(omega)``, so it is isotropic.
    """
    equation = matrix(FIELD, [(root.conjugate() * FORM).list()])
    basis = equation.right_kernel().basis()
    for bound in count(1):
        scalars = [FIELD(a) + FIELD(b) * OMEGA for a in range(-bound, bound + 1) for b in range(-bound, bound + 1)]
        for coefficients in product(scalars, repeat=len(basis)):
            if not any(coefficients):
                continue
            candidate = vector(FIELD, list(sum(c * b for c, b in zip(coefficients, basis))))
            if h(candidate, candidate) != ZERO:
                continue
            denominator = lcm([c.denominator() for entry in candidate for c in entry.vector()])
            integral = denominator * candidate
            ideal = FIELD.ideal(list(integral))
            primitive = integral / ideal.gens_reduced()[0]
            if is_primitive(primitive):
                return primitive
    raise AssertionError("an isotropic Hermitian space has a null vector")


def reduce_short_root(root):
    r"""Return ``g`` in ``Gamma`` with ``g root = r_7``, for a short root (Theorem 7.21)."""
    root = vector(FIELD, root)
    if h(root, root) != ONE:
        raise ValueError(f"{root} is not a short root")
    g = reduce_null_vector(_null_vector_orthogonal_to(root))
    root = g * root
    if root[3] != ZERO:
        raise ArithmeticError("a root orthogonal to rho has nonzero height")
    position = next(i for i in range(3) if root[i] != ZERO)
    unit = root[position]
    if unit.norm() != ONE or any(root[i] != ZERO for i in range(3) if i != position):
        raise ArithmeticError("a height-zero short root is not a unit coordinate vector")
    lam = vector(FIELD, [ZERO, ZERO, ZERO])
    lam[position] = (root[4] / unit).conjugate()
    step = _admissible_translation(lam)
    root = step * root
    g = step * g
    if root[3] != ZERO or root[4] != ZERO:
        raise ArithmeticError("the translation did not clear the root")
    target = ROOTS[7]
    letters = [*AUT_PLUS_GENERATORS, *(~r for r in AUT_PLUS_GENERATORS)]
    frontier = [(root, matrix.identity(FIELD, 5))]
    seen = {str(root)}
    while frontier:
        current, word = frontier.pop(0)
        if current == target:
            return word * g
        for letter in letters:
            image = letter * current
            if str(image) not in seen and all(c.norm() <= 1 for c in image):
                seen.add(str(image))
                frontier.append((image, letter * word))
    raise ArithmeticError(f"no bounded word carries {root} to r_7")


def root_transporter(r, s):
    r"""Return ``g`` in ``Gamma`` with ``g r = s`` for short roots (Theorem 7.21)."""
    g = ~reduce_short_root(s) * reduce_short_root(r)
    if g * vector(FIELD, r) != vector(FIELD, s) or not is_isometry(g):
        raise ArithmeticError("the transporter does not carry the short root")
    return g


def _unit_transport(unit):
    r"""Return an element of ``Gamma`` carrying ``unit rho`` to ``rho`` (7.15, 7.16).

    ``F = R_3 T_{0,0,0;-theta}`` has ``F^2 = omega`` on ``I_{1,1}``, and ``B``
    times the central involution of ``I_{1,1}`` carries ``rho`` to ``-rho``.
    """
    F = R[3] * translation(_root(0, 0, 0), -THETA)
    involution = matrix.identity(FIELD, 5)
    involution[3, 3] = involution[4, 4] = -ONE
    sign = BIFLECTION * involution
    for negate in (False, True):
        for power in range(3):
            candidate = (F * F) ** power * (sign if negate else matrix.identity(FIELD, 5))
            if candidate * (unit * RHO) == RHO:
                return candidate
    raise ArithmeticError(f"{unit} is not a unit of Z[omega]")


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


def _mod_theta(entries):
    return entries.apply_map(lambda c: GF(3)(int(c.vector()[0]) + int(c.vector()[1])), GF(3))


def spinor_character(g) -> int:
    r"""Return the spinor norm ``nu(g)`` of (2.12) on ``Lambda / theta Lambda``, as ``+1`` or ``-1``.

    Wall's form on ``im(1 - g)`` is ``chi((1-g) x, u) = b(x, u)``. A reflection in
    ``v`` has ``chi(v, v) = q(v) / 2``, so ``nu = det(chi) 2^k`` for ``k = rk(1 - g)``.
    """
    if not is_isometry(g):
        raise ValueError(f"{g} is not an isometry of the Eisenstein lattice")
    reduced = _mod_theta(g)
    form = _mod_theta(FORM)
    difference = matrix.identity(GF(3), 5) - reduced
    pivots = difference.pivots()
    images = difference.matrix_from_columns(pivots)
    sources = matrix.identity(GF(3), 5).matrix_from_columns(pivots)
    chi = images.transpose() * form * sources
    value = chi.det() * GF(3)(2) ** len(pivots) if pivots else GF(3)(1)
    return 1 if value == 1 else -1


def is_aut_plus(g) -> bool:
    r"""Return whether ``g`` lies in the index-two subgroup ``Aut^+`` of (2.12)."""
    return spinor_character(g) == 1


def is_gamma_theta(g) -> bool:
    r"""Return whether ``g`` lies in ``Aut^+`` and acts trivially on ``Lambda / theta Lambda`` (7.24)."""
    return is_aut_plus(g) and _mod_theta(g).is_one()


def triflection(root):
    r"""Return the ``omega``-reflection in a short root, the generators of ``Gamma_theta`` in Theorem 7.25."""
    if h(root, root) != ONE:
        raise ValueError(f"{root} is not a short root")
    return reflection(root, OMEGA)


def lattice() -> HermitianLattice:
    r"""Return ``Z[omega]^5`` with the source form; its trace lattice has rank 10."""
    return HermitianLattice((ONE, OMEGA), FORM.transpose())


def trace_generators(generators=FULL_GENERATORS) -> tuple:
    r"""Return the generators as integral isometries of the rank-ten trace lattice."""
    hermitian = lattice()
    return tuple(hermitian.isometry_from_matrix(g) for g in generators)
