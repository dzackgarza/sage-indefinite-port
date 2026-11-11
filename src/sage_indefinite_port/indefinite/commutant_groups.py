r"""Centralizers and equivariant transporters for noncompact isometry commutants (#54).

Let `A` be the rational commutant of isometries `F_i` of a nonsingular integral Gram matrix
`B`, `O = A \cap M_r(Z)` and `a^\dagger = B^{-1} a^t B`. The centralizer `C(Z)` is the set of
units `u \in O` with `u^\dagger u = 1`. A transporter to `(B', F'_i)` is `t \in GL_r(Z)` with
`t F_i = F'_i t` and `t^t B' t = B`; every rational one is `t_0 a` for one invertible rational
intertwiner `t_0` and `a \in A` with `a^\dagger h a = 1`, `h = B^{-1} t_0^t B' t_0 \in A`.

The provider :class:`IsometryCommutant` treats compact, split, commutative and radical
centralizers. This module adds the routes that need the port's arithmetic groups:

- `A` an indefinite division quaternion algebra over `Q`. An involution of the first kind
  is the canonical one (symplectic) or `x \mapsto c^{-1} \bar x c` for a pure `c` (orthogonal)
  (Knus, Merkurjev, Rost, Tignol, *The Book of Involutions*, Prop. 2.21). In the symplectic
  case `C(Z) = O^1`, generated from the Voronoi unit generators of `O`. In the orthogonal case
  `C(Z)` is the norm-one torus of `Q(c)` in `O` and at most one coset of anticommuting `w`.
- `A / rad = Q`: the transporter coset `t_r (\pm N(Q))` is searched for an integral point.

Other noncompact commutants raise :class:`NotImplementedError` naming the missing group.
"""

from itertools import product

from sage.arith.functions import lcm
from sage.arith.misc import binomial, gcd
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.libs.pari import pari
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.quadratic_forms.binary_qf import BinaryQF
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage_extensions.categories.lattice_centralizers import finite_preimage_generators
from sage_indefinite_port.indefinite.quaternion_voronoi import QuaternionVoronoi


def intertwiner_basis(sources, targets):
    r"""Return a `Q`-basis of `\{t : t F_i = T_i t\}`."""
    r = sources[0].nrows()
    rows = []
    for F, T in zip(sources, targets, strict=True):
        for i in range(r):
            for j in range(r):
                row = [0] * r**2
                for k in range(r):
                    row[i * r + k] += F[k, j]
                    row[k * r + j] -= T[i, k]
                rows.append(row)
    return [matrix(QQ, r, r, list(v)) for v in matrix(QQ, rows).right_kernel().basis()]


def rational_intertwiner(sources, targets):
    r"""Return an invertible rational intertwiner, or ``None`` when there is none.

    `\det(\sum x_m t_m)` is a polynomial of degree `r` in the coordinates. If it is not zero,
    it is not zero at a point of `\{0, \ldots, r\}^n` (Schwartz, J. ACM 27 (1980), Lemma 1).
    """
    basis = intertwiner_basis(sources, targets)
    if not basis:
        return None
    r = sources[0].nrows()
    for values in product(range(r + 1), repeat=len(basis)):
        candidate = sum(c * m for c, m in zip(values, basis))
        if candidate.det() != 0:
            return candidate
    return None


def solution_lattice(commutant, start):
    r"""Return a `Z`-basis, as elements of `A`, of `\{a \in A : t_0 a \in M_r(Z)\}`."""
    basis = [m.change_ring(QQ) for m in commutant.order_basis]
    rows = matrix(QQ, [(start * b).list() for b in basis])
    scale = rows.denominator()
    lattice = matrix(ZZ, scale * rows).saturation()
    return [
        sum(c * b for c, b in zip(rows.solve_left(vector(QQ, row)), basis))
        for row in lattice.rows()
    ]


def _verified(commutant, candidate, target_gram, targets):
    if candidate.denominator() != 1 or candidate.det() not in (1, -1):
        return None
    assert candidate.transpose() * target_gram * candidate == commutant.gram
    assert all(candidate * F == T * candidate for F, T in zip(commutant.isometries, targets, strict=True))
    return matrix(ZZ, candidate)


def _rational_square_root(value):
    value = QQ(value)
    return value.sqrt() if value > 0 and value.is_square() else None


def _integral_kernel(elements, equation):
    r"""Return the `Z`-combinations of ``elements`` (matrices) killed by the linear ``equation``."""
    rows = matrix(QQ, [equation(x).list() for x in elements])
    scale = rows.denominator()
    return [sum(c * x for c, x in zip(v, elements)) for v in matrix(ZZ, scale * rows).left_kernel().basis()]


def _binary_representation(norm, first, second, value):
    r"""Return an integral `(x, y)` with `n(x w_1 + y w_2) =` ``value``, or ``None``.

    The binary form is scaled to a primitive integral form and decided by PARI's ``qfbsolve``.
    """
    a, c = norm(first), norm(second)
    b = norm(first + second) - a - c
    scale = lcm([QQ(a).denominator(), QQ(b).denominator(), QQ(c).denominator(), QQ(value).denominator()])
    coefficients = [ZZ(scale * a), ZZ(scale * b), ZZ(scale * c)]
    content = gcd(coefficients)
    target = scale * value
    if target % content != 0:
        return None
    solution = BinaryQF([x // content for x in coefficients]).solve_integer(ZZ(target / content))
    return None if solution is None else (ZZ(solution[0]), ZZ(solution[1]))


class QuaternionCommutant:
    r"""An explicit isomorphism `A \to (a, b)_Q` for a division quaternion commutant `A`.

    `I` is the pure part `x - \mathrm{Tr}(x)/r` of an order basis element, so `I^2 = a` is
    rational, and `J = I y - y I` anticommutes with `I`. When `a, b > 0`, `(IJ, J)` replaces
    `(I, J)` since `(IJ)^2 = -ab`; then `a < 0 < b` whenever `A_R \cong M_2(R)`.
    """

    def __init__(self, commutant):
        self.commutant = commutant
        r = commutant.rank
        one = identity_matrix(QQ, r)
        basis = [m.change_ring(QQ) for m in commutant.order_basis]
        assert len(basis) == 4 and not commutant.radical_basis()
        pure = lambda x: x - x.trace() / r * one
        first = next(pure(m) for m in basis if not pure(m).is_zero())
        first = first * first.denominator()
        second = next(c for c in (first * pure(m) - pure(m) * first for m in basis) if not c.is_zero())
        second = second * second.denominator()
        a, b = (first * first)[0, 0], (second * second)[0, 0]
        if a * b == 0:
            # A pure element of square zero: `A \cong M_2(Q)`.
            self.algebra = None
            return
        assert first * first == a * one and second * second == b * one and first * second == -second * first
        if a > 0 and b > 0:
            first, a = first * second, -a * b
        elif a > 0 > b:
            first, second, a, b = second, first, b, a
        self.a, self.b = ZZ(a), ZZ(b)
        self.algebra = QuaternionAlgebra(QQ, self.a, self.b)
        self.frame = (one, first, second, first * second)
        self._frame_rows = matrix(QQ, [m.list() for m in self.frame])
        self.order_rows = matrix(QQ, [list(self.to_model(m)) for m in basis])
        table = matrix(QQ, commutant.involution_table())
        self._symmetric_dimension = (table - identity_matrix(QQ, 4)).left_kernel().dimension()
        self._skew = (table + identity_matrix(QQ, 4)).left_kernel().basis()

    def is_division(self):
        return self.algebra is not None and self.algebra.is_division_algebra()

    def is_symplectic(self):
        r"""Decide whether the adjoint is the canonical involution: its symmetric space is `Q`."""
        assert self._symmetric_dimension in (1, 3)
        return self._symmetric_dimension == 1

    def to_model(self, element):
        coordinates = self._frame_rows.solve_left(vector(QQ, matrix(QQ, element).list()))
        return self.algebra(list(coordinates))

    def from_model(self, quaternion):
        return sum(QQ(c) * m for c, m in zip(list(quaternion), self.frame))

    def norm(self, element):
        return self.to_model(element).reduced_norm()

    def voronoi(self):
        if not hasattr(self, "_voronoi"):
            assert self.a < 0 < self.b, f"the commutant ({self.a}, {self.b})_Q is definite"
            self._voronoi = QuaternionVoronoi(int(self.a), int(self.b), self.order_rows)
        return self._voronoi

    def _pure_skew(self):
        r"""Return the integral generator `c` of the skew elements, with `x^\dagger = c^{-1} \bar x c`."""
        assert len(self._skew) == 1
        c = sum(QQ(x) * m.change_ring(QQ) for x, m in zip(self._skew[0], self.commutant.order_basis))
        return c * c.denominator()

    def centralizer_generators(self):
        assert self.is_division()
        if self.is_symplectic():
            norm_one, central = self.voronoi().norm_one_generators()
            result = [matrix(ZZ, self.from_model(u)) for u in (*norm_one, *central)]
            assert all(self.commutant.is_centralizer_element(g) for g in result)
            return tuple(result)
        return self._torus_generators()

    def _torus_generators(self):
        r"""Return generators of `C(Z)` for an orthogonal adjoint `x^\dagger = c^{-1} \bar x c`.

        `u^\dagger u = 1` gives `n(u) = \pm 1`. If `n(u) = 1` then `u` commutes with `c`, so `u` lies in
        `Q(c)^1`, and an integral one lies in `O_{Q(c)}^1`. If `n(u) = -1` then `u` anticommutes with `c`,
        and conversely such a `w` satisfies `w^\dagger w = w^2 = -n(w)`. So `C(Z)` is generated by the
        `O`-points of `\langle -1, \varepsilon \rangle` (Schreier) and one `w \in O` with `wc = -cw`,
        `n(w) = -1`, decided as a value of the binary norm form.
        """
        c = self._pure_skew()
        d = (c * c)[0, 0]
        assert d > 0 and not QQ(d).is_square(), f"the torus Q(sqrt({d})) is compact"
        polynomial = PolynomialRing(QQ, "x")([-d, 0, 1])
        field = NumberField(polynomial, "s")
        one = identity_matrix(QQ, self.commutant.rank)
        image = lambda u: sum(QQ(x) * c**k for k, x in enumerate(u.list()))
        fundamental = field(field.unit_group().fundamental_units()[0])
        norm_one = fundamental if fundamental.norm() == 1 else fundamental**2
        moves = [-one, image(norm_one)]
        order = [m.change_ring(QQ) for m in self.commutant.order_basis]
        anticommuting = _integral_kernel(order, lambda x: x * c + c * x)
        assert len(anticommuting) == 2
        solution = _binary_representation(self.norm, *anticommuting, -1)
        if solution is not None:
            moves.append(solution[0] * anticommuting[0] + solution[1] * anticommuting[1])
        result, _index = finite_preimage_generators(moves, lambda x: x, self.commutant.in_order)
        result = [matrix(ZZ, g) for g in result]
        assert all(self.commutant.is_centralizer_element(g) for g in result)
        return tuple(result)

    def equivariant_transporter(self, target_gram, targets):
        r"""Return an integral transporter to ``(target_gram, targets)``, or ``None`` as a proof of absence."""
        commutant = self.commutant
        start = rational_intertwiner(commutant.isometries, targets)
        if start is None:
            return None
        square = commutant.gram.change_ring(QQ).inverse() * start.transpose() * target_gram * start
        assert commutant.contains(square)
        lattice = solution_lattice(commutant, start)
        if self.is_symplectic():
            return self._symplectic_transporter(start, square, lattice, target_gram, targets)
        return self._orthogonal_transporter(start, square, lattice, target_gram, targets)

    def _symplectic_transporter(self, start, square, lattice, target_gram, targets):
        r"""Decide `\bar a h a = 1` on the lattice `M` of `a` with `t_0 a` integral.

        `h` is self-adjoint, hence a rational `\eta`, and the equation is `n(a) = \nu = 1/\eta`. `M` is a
        right `O`-module, so `\bar M` is a left `O`-module. By :meth:`QuaternionVoronoi.cell_trace_bound`,
        every `b \in \bar M` with `n(b) = \nu` has a unit translate `ub` with `\mathrm{Tr}(\beta\beta^t)` at most
        `|\nu|` times the bound, `\beta` its real image, and `n(ub) = \pm\nu`. A unit of norm `-1`, if any, is
        among the unit generators. So the finite enumeration decides existence.
        """
        eta = square[0, 0]
        assert square == eta * identity_matrix(QQ, self.commutant.rank)
        value = 1 / eta
        voronoi = self.voronoi()
        conjugates = [self.to_model(a).conjugate() for a in lattice]
        images = [voronoi.real_matrix(q) for q in conjugates]
        gram = matrix(RDF, [[float((x * y.transpose()).trace()) for y in images] for x in images])
        bound = float(abs(value) * voronoi.cell_trace_bound()) * (1 + 1e-9)
        found = pari(gram).qfminim(bound, None, 2)[2]
        columns = int(pari(found).matsize()[1])
        reversing = next((u for u in voronoi.unit_generators() if u.reduced_norm() == -1), None)
        for column in range(columns):
            coordinates = [int(round(float(found[row, column]))) for row in range(len(conjugates))]
            element = sum(c * q for c, q in zip(coordinates, conjugates))
            norm = element.reduced_norm()
            if norm == -value and reversing is not None:
                element, norm = reversing * element, value
            if norm != value:
                continue
            return _verified(self.commutant, start * self.from_model(element.conjugate()), target_gram, targets)
        return None

    def _orthogonal_transporter(self, start, square, lattice, target_gram, targets):
        r"""Decide `c^{-1} \bar a c h a = 1`, that is `\bar a g a = c` for `g = c h`.

        Norms give `n(a)^2 = n(c)/n(g) = \sigma^2`. With `n(a) = \epsilon\sigma`, `\bar a = \epsilon\sigma a^{-1}`, so
        the equation is the linear condition `g a = a\, \epsilon c / \sigma`, and conversely. Its solutions in `M`
        form a rank-two lattice on which `n` is a multiple of the norm form of `Q(c)`; the value
        `\epsilon\sigma` is decided by :func:`_binary_representation`.
        """
        c = self._pure_skew()
        g = c * square
        sigma = _rational_square_root(self.norm(c) / self.norm(g))
        if sigma is None:
            return None
        for sign in (1, -1):
            right = sign * c / sigma
            solutions = _integral_kernel(lattice, lambda a: g * a - a * right)
            if not solutions:
                continue
            assert len(solutions) == 2
            found = _binary_representation(self.norm, *solutions, sign * sigma)
            if found is None:
                continue
            element = found[0] * solutions[0] + found[1] * solutions[1]
            return _verified(self.commutant, start * element, target_gram, targets)
        return None


def unipotent_transporter(commutant, target_gram, targets):
    r"""Return an integral transporter when `A / rad A = Q`, or ``None`` as a proof of absence.

    Here `C(Q) = \pm N(Q)` with `N = \exp(n)`. Write `h = \lambda (1 + x)`, `x \in rad A`; `x` is
    self-adjoint, so `a = \mu (1 + x)^{-1/2}` with `\mu^2 \lambda = 1` solves `a^\dagger h a = 1`, and no
    solution exists unless `1/\lambda` is a rational square. The integral transporters are the
    integral points of `t_0 a N(Q)`, decided by the provider's Mal'cev-coordinate search
    :meth:`IsometryCommutant._integral_point`, which requires only that the coset be `g N(Q)`.
    """
    start = rational_intertwiner(commutant.isometries, targets)
    if start is None:
        return None
    r = commutant.rank
    one = identity_matrix(QQ, r)
    square = commutant.gram.change_ring(QQ).inverse() * start.transpose() * target_gram * start
    assert commutant.contains(square)
    scalar = square.trace() / r
    nilpotent = square / scalar - one
    nilpotency = commutant._radical_nilpotency()[0]
    assert (nilpotent ** (nilpotency + 1)).is_zero() and commutant.adjoint(nilpotent) == nilpotent
    root = _rational_square_root(1 / scalar)
    if root is None:
        return None
    correction = sum(binomial(QQ((-1, 2)), k) * nilpotent**k for k in range(nilpotency + 1))
    rational = start * root * correction
    assert rational.transpose() * target_gram * rational == commutant.gram
    found = commutant._integral_point(rational, 1)
    return None if found is None else _verified(commutant, matrix(QQ, found), target_gram, targets)


def _missing_group(commutant):
    r"""Return the description of the arithmetic group of a commutant with no implemented route."""
    center = commutant._center_basis()
    if commutant.radical_basis():
        quotient = commutant.dimension() - len(commutant.radical_basis())
        return (
            f"the integral transporter of the commutant of {commutant.isometries} needs reduction theory "
            f"for the unitary group of its Levi quotient of dimension {quotient} over Q (#54)"
        )
    kind = "first" if all(commutant.adjoint(z) == z for z in center) else "second"
    if len(center) == commutant.dimension():
        return (
            f"the integral transporter of the commutative commutant of dimension {commutant.dimension()} of "
            f"{commutant.isometries} needs the norm-one unit orbit decision of its etale algebra (#54)"
        )
    if len(center) == 1:
        return (
            f"the unitary group of the central simple commutant of dimension {commutant.dimension()} over Q of "
            f"{commutant.isometries}, with an involution of the first kind, has no implemented generation or "
            f"transporter route (#54)"
        )
    generic = next(
        z for values in product(range(-3, 4), repeat=len(center))
        for z in (sum(c * m for c, m in zip(values, center)),)
        if z.minimal_polynomial().degree() == len(center)
    )
    polynomial = generic.minimal_polynomial()
    degree = ZZ(commutant.dimension() // len(center)).isqrt()
    return (
        f"the centralizer of {commutant.isometries} is the unitary group over O_K, K = Q[x]/({polynomial}), "
        f"of a rank-{degree} form with an involution of the {kind} kind on M_{degree}(K) or a K-division "
        f"algebra; its generation is not implemented (#54)"
    )


def _quaternion(commutant):
    r"""Return the quaternion model when `A` is a central division quaternion algebra over `Q`, else ``None``."""
    if commutant.dimension() != 4 or commutant.radical_basis() or len(commutant._center_basis()) != 1:
        return None
    model = QuaternionCommutant(commutant)
    return model if model.is_division() else None


def _central_field(commutant):
    r"""Decide whether `A` is semisimple, noncommutative, and has a field of degree `> 1` as centre."""
    if commutant.radical_basis():
        return False
    center = commutant._center_basis()
    if len(center) in (1, commutant.dimension()):
        return False
    for values in product(range(-3, 4), repeat=len(center)):
        element = sum(c * m for c, m in zip(values, center))
        polynomial = element.minimal_polynomial()
        if polynomial.degree() == len(center):
            return polynomial.is_irreducible()
    raise ArithmeticError(f"no generator of the centre of the commutant of {commutant.isometries} in the box")


def centralizer_generators(commutant):
    r"""Return generators of `C(Z)`, dispatching on the structure of `A`."""
    if commutant.is_compact():
        return commutant.compact_centralizer_elements()
    model = _quaternion(commutant)
    if model is not None:
        return model.centralizer_generators()
    if _central_field(commutant):
        raise NotImplementedError(_missing_group(commutant))
    return commutant.centralizer_generators()


def equivariant_transporter(commutant, target_gram, targets):
    r"""Return `t \in GL_r(Z)` with `t^t B' t = B` and `t F_i = F'_i t`, or ``None`` as a proof of absence."""
    target_gram = matrix(ZZ, target_gram)
    targets = tuple(matrix(ZZ, commutant.rank, commutant.rank, T) for T in targets)
    assert len(targets) == len(commutant.isometries)
    if commutant.is_compact():
        return commutant.compact_equivariant_transporter(target_gram, targets)
    radical = commutant.radical_basis()
    if radical and commutant.dimension() - len(radical) == 1:
        return unipotent_transporter(commutant, target_gram, targets)
    model = _quaternion(commutant)
    if model is not None:
        return model.equivariant_transporter(target_gram, targets)
    commutative = all(x * y == y * x for x in commutant.order_basis for y in commutant.order_basis)
    if not radical and commutative and commutant.dimension() == 2 and len(targets) == 1:
        return commutant.equivariant_transporter(target_gram, targets[0])
    raise NotImplementedError(_missing_group(commutant))
