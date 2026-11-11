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
- `A / rad = M_n(Q)`, `n \ge 2`: the transporter is reduced to an isometry of integral forms on
  the Levi module and a coset search in the unipotent radical (:func:`levi_transporter`).
- `A` central simple over a number field `K \ne Q`: `A = End_K(W)` for a simple module `W`, the
  adjoint is the adjoint of a Hermitian, symmetric or alternating form `h` on `W`, and `C(Z)` is
  the subgroup of `U(\Lambda, h)` lying in `O`, `\Lambda = O_K O f` (:class:`CentralFieldCommutant`).

Other noncompact commutants raise :class:`NotImplementedError` naming the missing group.
"""

from itertools import combinations, product

from sage.arith.functions import lcm
from sage.arith.misc import binomial, gcd
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.libs.pari import pari
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.quadratic_forms.binary_qf import BinaryQF
from sage.quadratic_forms.qfsolve import qfsolve
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage_extensions.algebras.isometry_commutant import _form_generators, _integer_solution, _short_vectors
from sage_extensions.categories.lattice_centralizers import finite_preimage_generators
from sage_extensions.categories.lattices import Lattices
from sage_indefinite_port.indefinite import (
    act_eisenstein,
    binary_hyperbolic_unitary,
    eisenstein_hyperbolic_unitary,
    euclidean_picard,
    gaussian_picard,
)
from sage_indefinite_port.indefinite.quaternion_voronoi import QuaternionVoronoi
from sage_indefinite_port.indefinite.witt_one_hilbert import HilbertWittOneLattice


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


def _scaled_integral(*grams):
    r"""Return the Gram matrices scaled by one common positive integer to integral matrices."""
    scale = lcm([matrix(QQ, g).denominator() for g in grams])
    return tuple(matrix(ZZ, scale * matrix(QQ, g)) for g in grams)


def _form_isometry(source, target):
    r"""Return `P \in GL_n(Z)` with `P^t T P = S`, or ``None`` when the forms are not isometric.

    Symmetric forms are decided by the provider's isometry Mor of :class:`Lattices`. Alternating
    forms are isometric exactly when their symplectic (Frobenius) normal forms agree
    (Newman, *Integral Matrices*, Thm. IV.1), and ``symplectic_form`` returns `C` with `C S C^t = F`.
    """
    source, target = _scaled_integral(source, target)
    if source == -source.transpose():
        normal, change = source.symplectic_form()
        other, other_change = target.symplectic_form()
        if normal != other:
            return None
        result = other_change.transpose() * change.transpose().inverse()
    else:
        assert source == source.transpose() and target == target.transpose()
        left = Lattices(ZZ)([[int(c) for c in row] for row in source.rows()])
        right = Lattices(ZZ)([[int(c) for c in row] for row in target.rows()])
        isometry = left.isometry_to(right)
        if isometry is None:
            return None
        result = matrix(QQ, [[QQ(c) for c in isometry(v)] for v in left.module_generators()]).transpose()
    assert result.transpose() * target * result == source and abs(result.det()) == 1
    return matrix(ZZ, result)


def _sublattices(size, index):
    r"""Yield every sublattice of `Z^n` of the given index, as the columns of its Hermite normal form."""
    def diagonals(remaining, length):
        if length == 0:
            if remaining == 1:
                yield ()
            return
        for d in remaining.divisors():
            for rest in diagonals(remaining // d, length - 1):
                yield (d, *rest)

    for diagonal in diagonals(ZZ(index), size):
        slots = [(i, j) for i in range(size) for j in range(i + 1, size)]
        for values in product(*(range(diagonal[i]) for i, _j in slots)):
            basis = matrix(ZZ, size, size)
            for k in range(size):
                basis[k, k] = diagonal[k]
            for (i, j), value in zip(slots, values):
                basis[i, j] = value
            yield basis


class LeviFrame:
    r"""The Levi module of a commutant with `A / rad A \cong M_n(Q)` and `O / (O \cap rad) = M_n(Z)`.

    This is the frame of :meth:`IsometryCommutant.radical_centralizer_generators`: an idempotent `e`
    of rank `r / n`, the module `\bar M = O e / (rad\, e) \cong Z^n` on which `A / rad` acts faithfully,
    and the form `b` with `\bar x^\dagger \bar y = b(x, y)\, t` modulo the radical. For a self-adjoint
    `h \in A` the twisted form `b_h(x, y)` is read from `x^\dagger h y` in the same way.
    """

    def __init__(self, commutant):
        self.commutant = commutant
        radical = commutant.radical_basis()
        assert radical
        size = commutant.dimension()
        radical_rows = matrix(QQ, [commutant.coordinates(x) for x in radical])
        radical_lattice = matrix(ZZ, radical_rows.denominator() * radical_rows).saturation()
        _diagonal, _left, change = radical_lattice.smith_form()
        kernel_rank = radical_lattice.nrows()
        self._project = lambda element: (commutant.coordinates(element) * change.change_ring(QQ))[kernel_rank:]
        inverse_change = change.inverse().change_ring(QQ)
        order = [m.change_ring(QQ) for m in commutant.order_basis]
        self._section = tuple(
            sum(c * m for c, m in zip(inverse_change.row(kernel_rank + i), order)) for i in range(size - kernel_rank)
        )
        quotient = size - kernel_rank
        degree = ZZ(quotient).isqrt()
        if degree**2 != quotient or degree < 2:
            raise NotImplementedError(
                f"the semisimple quotient of dimension {quotient} of the commutant of {commutant.isometries} is not "
                f"M_n(Q) with n >= 2; its transporter needs the factor-wise Levi decision (#54)"
            )
        self.degree = degree
        try:
            e = commutant._primitive_idempotent(degree)
        except ArithmeticError as error:
            raise NotImplementedError(
                f"the Levi quotient of the commutant of {commutant.isometries} has no idempotent of rank "
                f"{commutant.rank // degree} in the sparse box; a division Levi quotient has no Levi-module "
                f"transporter route (#54)"
            ) from error
        images = matrix(QQ, [self._project(m * e) for m in order])
        hermite = matrix(ZZ, images.denominator() * images).hermite_form(include_zero_rows=False)
        assert hermite.nrows() == degree
        self._basis = matrix(QQ, hermite) / images.denominator()
        self._lifts = [
            sum(c * m for c, m in zip(images.solve_left(row), order)) * e for row in self._basis.rows()
        ]
        action_rows = matrix(QQ, [self.on_module(m).list() for m in self._section])
        if action_rows.denominator() != 1 or abs(action_rows.det()) != 1:
            raise NotImplementedError(
                f"the semisimple quotient of the order of the commutant of {commutant.isometries} is not End(M); "
                f"its Levi transporter needs the nonmaximal Levi order (#54)"
            )
        self._action_inverse = action_rows.inverse()
        adjoint_e = commutant.adjoint(e)
        line = matrix(QQ, [self._project(adjoint_e * m * e) for m in order])
        generator = matrix(QQ, matrix(ZZ, line.denominator() * line).hermite_form(include_zero_rows=False))
        self._line = generator / line.denominator()
        assert self._line.nrows() == 1
        self.gram = matrix(ZZ, self.form(identity_matrix(QQ, commutant.rank)))

    def on_module(self, element):
        return matrix(QQ, [self._basis.solve_left(self._project(element * m)) for m in self._lifts]).transpose()

    def form(self, middle):
        r"""Return the Gram matrix of `b_h` on `\bar M`, `h =` ``middle``."""
        adjoint = self.commutant.adjoint
        return matrix(
            QQ, [[self._line.solve_left(self._project(adjoint(x) * middle * y))[0] for y in self._lifts] for x in self._lifts]
        )

    def lift(self, start, middle, isometry):
        r"""Return an integral `t_0 a` with `\bar a =` ``isometry`` and `a^\dagger h a = 1`, or ``None``.

        A preimage `p` has `d = p^\dagger h p - 1` self-adjoint in the radical, so `a = p (1 + d)^{-1/2}`
        solves the equation. All solutions over `\bar a` are `a N(Q)`, decided by `_integral_point`.
        """
        commutant = self.commutant
        combination = vector(QQ, matrix(QQ, isometry).list()) * self._action_inverse
        preimage = sum(c * m for c, m in zip(combination, self._section))
        defect = commutant.adjoint(preimage) * middle * preimage - identity_matrix(QQ, commutant.rank)
        nilpotency = commutant._radical_nilpotency()[0]
        assert (defect ** (nilpotency + 1)).is_zero() and commutant.adjoint(defect) == defect
        root = sum(binomial(QQ((-1, 2)), k) * defect**k for k in range(nilpotency + 1))
        return commutant._integral_point(start * preimage * root, 1)

    def coset_representatives(self):
        r"""Return representatives of the left cosets of the image `\bar\Gamma` of `C(Z)` in `O(\bar M, b)`.

        `\bar\Gamma` is the set of liftable isometries, of finite index (Platonov-Rapinchuk, Thm. 4.1).
        `x \bar\Gamma = y \bar\Gamma` exactly when `y^{-1} x` lifts; generators permute the finite coset set.
        """
        if hasattr(self, "_cosets"):
            return self._cosets
        one = identity_matrix(QQ, self.commutant.rank)
        generators = [matrix(QQ, g) for g in _form_generators(self.gram)]
        representatives = [identity_matrix(QQ, self.degree)]
        position = 0
        while position < len(representatives):
            current = representatives[position]
            position += 1
            for generator in generators:
                moved = generator * current
                if not any(self.lift(one, one, known.inverse() * moved) is not None for known in representatives):
                    representatives.append(moved)
        self._cosets = tuple(representatives)
        return self._cosets


def levi_transporter(commutant, target_gram, targets):
    r"""Return an integral transporter when `A / rad A = M_n(Q)`, or ``None`` as a proof of absence.

    Every transporter is `t_0 a` with `a \in X = \{a : t_0 a \in M_r(Z)\}` and `a^\dagger h a = 1`. The image
    `\bar X` is a right `M_n(Z)`-module, so `\bar X = \{\phi : \phi \bar M \subseteq Y\}` for `Y = \bar X \bar M`
    (Reiner, *Maximal Orders*, Thm. 21.6). The equation makes `\bar a` an isometry of `(\bar M, b)` onto
    a sublattice `Z \subseteq Y` for `b_h`, of index `k` with `k^2 \det b_h|_Y = \det b`. For each such
    `Z` and one isometry `\bar a_0`, the integral solutions over `\bar a_0 O(\bar M)` are right
    `C(Z)`-invariant, so one representative of each left `\bar\Gamma`-coset is lifted by :meth:`LeviFrame.lift`.
    """
    start = rational_intertwiner(commutant.isometries, targets)
    if start is None:
        return None
    square = commutant.gram.change_ring(QQ).inverse() * start.transpose() * target_gram * start
    assert commutant.contains(square)
    frame = LeviFrame(commutant)
    twisted = frame.form(square)
    columns = [column for x in solution_lattice(commutant, start) for column in frame.on_module(x).columns()]
    generators = matrix(QQ, columns)
    span = matrix(ZZ, generators.denominator() * generators).hermite_form(include_zero_rows=False)
    image = (matrix(QQ, span) / generators.denominator()).transpose()
    assert image.nrows() == image.ncols() == frame.degree
    ratio = QQ(frame.gram.det()) / (image.transpose() * twisted * image).det()
    if ratio <= 0 or not ratio.is_square() or ratio.sqrt().denominator() != 1:
        return None
    for hermite in _sublattices(frame.degree, ratio.sqrt()):
        sublattice = image * hermite
        isometry = _form_isometry(frame.gram, sublattice.transpose() * twisted * sublattice)
        if isometry is None:
            continue
        initial = sublattice * isometry
        assert initial.transpose() * twisted * initial == frame.gram
        for representative in frame.coset_representatives():
            found = frame.lift(start, square, initial * representative)
            if found is not None:
                return _verified(commutant, matrix(QQ, found), target_gram, targets)
    return None


# Real quadratic fields with norm-Euclidean maximal order (Chatland-Davenport, Canad. J. Math. 2 (1950);
# Barnes-Swinnerton-Dyer, Acta Math. 87 (1952)): the Euclidean reduction of witt_one_hilbert applies.
_NORM_EUCLIDEAN_REAL = (2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41, 57, 73)

# Imaginary quadratic discriminants and the `d` of `Q(\sqrt{-d})` of the port's unitary models.
_IMAGINARY_MODELS = {-4: 1, -8: 2, -3: 3, -7: 7, -11: 11}


class CentralFieldCommutant:
    r"""An explicit identification of a central simple commutant over `K \ne Q` with a form over `O_K`.

    Let `z` generate the centre, `K = Q[z]` and `\sigma` the adjoint on `K`. A rank-one element `f`
    (rank `r / n` on `Q^r` for `A \cong M_n(K)`) gives the simple module `W = A f`, a `K`-space of
    dimension `n`, and `A = End_K(W)`. The space `f^\dagger A f` is a `K`-line `K m_0`, and
    `y^\dagger x = h(x, y)\, m_0` defines a `\sigma`-sesquilinear form with `(a y)^\dagger (a x) = y^\dagger a^\dagger a x`,
    so `a^\dagger a = 1` exactly when `a` preserves `h` (`A` is simple, so the pairing is nondegenerate).
    `h` is `\lambda`-Hermitian; a scalar makes it Hermitian (second kind) or leaves it symmetric or
    alternating (first kind) (Knus-Merkurjev-Rost-Tignol, *The Book of Involutions*, Thm. 4.2).
    `\Lambda = O_K O f` is `O`-stable, so `C(Z) = U(\Lambda, h) \cap O`, of finite index in `U(\Lambda, h)`.
    """

    def __init__(self, commutant):
        self.commutant = commutant
        r = commutant.rank
        center = commutant._center_basis()
        generic = next(
            z
            for values in product(range(-3, 4), repeat=len(center))
            for z in (sum(c * m for c, m in zip(values, center)),)
            if z.minimal_polynomial().degree() == len(center)
        )
        generic = generic * generic.denominator()
        self.degree = len(center)
        polynomial = generic.minimal_polynomial()
        assert polynomial.is_irreducible()
        self.field = NumberField(polynomial, "z")
        self._powers = tuple(generic**k for k in range(self.degree))
        powers = matrix(QQ, [p.list() for p in self._powers])
        image = powers.solve_left(vector(QQ, commutant.adjoint(generic).list()))
        self.sigma = self.field.hom([self.field(list(image))])
        self.kind = "first" if self.sigma(self.field.gen()) == self.field.gen() else "second"
        n = ZZ(commutant.dimension() // self.degree).isqrt()
        assert n**2 * self.degree == commutant.dimension()
        self.n = n
        order = [m.change_ring(QQ) for m in commutant.order_basis]
        self._order = order
        f = self._rank_one_element()
        basis, span = [], matrix(QQ, 0, r * r)
        for x in (m * f for m in order):
            stacked = span.stack(matrix(QQ, [(p * x).list() for p in self._powers]))
            if stacked.rank() > span.rank():
                basis.append(x)
                span = stacked
        assert len(basis) == n and span.rank() == n * self.degree
        self._basis = tuple(basis)
        self._basis_rows = span
        flat = matrix(QQ, [self._flatten(self.on_module(m)) for m in order])
        assert flat.is_square() and flat.det() != 0
        self._flat_inverse = flat.inverse()
        adjoint_f = commutant.adjoint(f)
        line = next(v for v in (adjoint_f * m * f for m in [identity_matrix(QQ, r), *order]) if not v.is_zero())
        self._line_rows = matrix(QQ, [(p * line).list() for p in self._powers])
        gram = matrix(self.field, [[self._pairing(b, c) for b in basis] for c in basis])
        self.form_kind, self.gram = self._normalized(gram)
        integers = self.field.ring_of_integers().basis()
        rows = matrix(QQ, [self._q_coordinates(self.embed(k) * m * f) for k in integers for m in order])
        hermite = matrix(ZZ, rows.denominator() * rows).hermite_form(include_zero_rows=False)
        self._lattice = matrix(QQ, hermite) / rows.denominator()
        self._lattice_vectors = tuple(self._k_vector(row) for row in self._lattice.rows())

    def embed(self, k):
        return sum(QQ(c) * p for c, p in zip(self.field(k).list(), self._powers))

    def _evaluate(self, polynomial, y):
        return sum(self.embed(c) * y**k for k, c in enumerate(polynomial.list()))

    def _rank_one_element(self):
        r"""Return an element of rank `r / n`, from `q(y)` for `K`-factors `q` of minimal polynomials of `y`."""
        target = self.commutant.rank // self.n
        candidates = [*self._order, *(x * y - y * x for x, y in combinations(self._order, 2))]
        found = []
        for y in candidates:
            factors = y.minimal_polynomial().change_ring(self.field).factor()
            if len(factors) == 1 and factors[0][1] == 1:
                continue
            for q, _multiplicity in factors:
                value = self._evaluate(q, y)
                if not value.is_zero():
                    found.append(value)
        found.sort(key=lambda x: x.rank())
        for x in found:
            if x.rank() == target:
                return x
        for x, y in product(found, repeat=2):
            for m in self._order:
                value = x * m * y
                if not value.is_zero() and value.rank() == target:
                    return value
        raise NotImplementedError(
            f"no element of rank {target} in the commutant of {self.commutant.isometries}: it is M_{self.n} over a "
            f"division algebra over K = Q[x]/({self.field.polynomial()}), or the box missed its zero divisors; "
            f"the unitary group of a K-division algebra has no implemented route (#54)"
        )

    def _q_coordinates(self, w):
        return self._basis_rows.solve_left(vector(QQ, w.list()))

    def _k_vector(self, coordinates):
        d = self.degree
        return vector(self.field, [self.field(list(coordinates[i * d : (i + 1) * d])) for i in range(self.n)])

    def _q_vector(self, kvector):
        return vector(QQ, [c for k in kvector for c in self.field(k).list()])

    def coordinates(self, w):
        return self._k_vector(self._q_coordinates(w))

    def on_module(self, element):
        return matrix(self.field, [self.coordinates(element * b) for b in self._basis]).transpose()

    def _flatten(self, block):
        return [c for k in block.list() for c in self.field(k).list()]

    def element_of(self, block):
        r"""Return the element of `A` acting on `W` by ``block`` in the basis `b_i`."""
        combination = vector(QQ, self._flatten(block)) * self._flat_inverse
        return sum(c * m for c, m in zip(combination, self._order))

    def _pairing(self, x, y):
        coordinates = self._line_rows.solve_left(vector(QQ, (self.commutant.adjoint(y) * x).list()))
        return self.field(list(coordinates))

    def _star(self, block):
        return block.apply_map(self.sigma).transpose()

    def _normalized(self, gram):
        r"""Return the kind and Gram `G` (`G_{ij} = h(b_j, b_i)`) after the scalar making `h` (skew-)Hermitian."""
        i, j = next((i, j) for i in range(self.n) for j in range(self.n) if gram[i, j] != 0)
        ratio = self._star(gram)[i, j] / gram[i, j]
        assert self._star(gram) == ratio * gram
        if self.kind == "first":
            assert ratio in (1, -1)
            return ("symmetric" if ratio == 1 else "alternating"), gram
        z = self.field.gen()
        scalar = 1 + ratio if ratio != -1 else z - self.sigma(z)
        gram = scalar * gram
        assert self._star(gram) == gram
        return "hermitian", gram

    def h(self, x, y, gram):
        return (self._star(matrix(self.field, [y])) * gram * matrix(self.field, [x]).transpose())[0, 0]

    def _span_is_lattice(self, vectors):
        integers = self.field.ring_of_integers().basis()
        rows = matrix(QQ, [self._q_vector(k * v) for k in integers for v in vectors])
        hermite = matrix(ZZ, rows.denominator() * rows).hermite_form(include_zero_rows=False)
        return rows.rank() == self._lattice.nrows() and matrix(QQ, hermite) / rows.denominator() == self._lattice

    def _integral_combination(self, values, target):
        r"""Return `x \in Z^N` with `\sum x_a` ``values[a]`` `=` ``target`` in `K`, or ``None``."""
        rows = matrix(QQ, [self.field(v).list() for v in values] + [self.field(target).list()])
        scale = rows.denominator()
        equations = matrix(ZZ, scale * rows[:-1]).transpose()
        return _integer_solution(equations, vector(ZZ, scale * rows[-1]))

    def _primitive(self, v):
        r"""Return `\gamma v` with `I = \gamma O_K` the fractional ideal `\{k : k v \in \Lambda\}`."""
        rows = matrix(QQ, [self._lattice.solve_left(self._q_vector(p * v)) for p in self.field.power_basis()])
        span = matrix(ZZ, rows.denominator() * rows.transpose()).hermite_form(include_zero_rows=False)
        dual = (matrix(QQ, span) / rows.denominator()).inverse().transpose()
        ideal = self.field.ideal([self.field(list(row)) for row in dual.rows()])
        generator = ideal.gens_reduced()
        assert len(generator) == 1
        return generator[0] * v

    def _hyperbolic_pair(self, gram):
        r"""Return `(v, f, G')` with `v, f` null, `h'(f, v) = 1`, `h'(\Lambda, v) = O_K`, `h' = h / q`, or ``None``.

        Isotropy of `h` is isotropy of the trace form `Tr\, h(x, x) = 2 h(x, x)` (PARI ``qfsolve``).
        `f_1` with `h'(f_1, v) = 1` is corrected to `f = f_1 - t v`, `t + \bar t = h'(f_1, f_1)`.
        """
        vectors = self._lattice_vectors
        trace = matrix(QQ, [[self.h(x, y, gram).trace() for y in vectors] for x in vectors])
        found = qfsolve(trace)
        if found in ZZ:
            return None
        found = vector(QQ, found)
        found = found * found.denominator()
        found = found / gcd(list(found))
        v = self._primitive(sum(c * x for c, x in zip(found, vectors)))
        assert self.h(v, v, gram) == 0
        ideal = self.field.ideal([self.h(x, v, gram) for x in vectors])
        norm = QQ(ideal.norm())
        if not norm.is_square() or ideal != self.field.ideal(norm.sqrt()):
            return None
        gram = gram / norm.sqrt()
        solution = self._integral_combination([self.h(x, v, gram) for x in vectors], 1)
        assert solution is not None
        first = sum(c * x for c, x in zip(solution, vectors))
        height = self.h(first, first, gram)
        if height not in ZZ:
            return None
        integers = self.field.ring_of_integers().basis()
        shift = _integer_solution(matrix(ZZ, [[ZZ(k.trace()) for k in integers]]), vector(ZZ, [ZZ(height)]))
        if shift is None:
            return None
        f = first - sum(c * k for c, k in zip(shift, integers)) * v
        assert self.h(f, f, gram) == 0 and self.h(f, v, gram) == 1
        return v, f, gram

    def _orthogonal(self, vectors, against, gram):
        r"""Return a `Z`-basis, as `K`-vectors, of the elements of span(``vectors``) orthogonal to ``against``."""
        rows = matrix(QQ, [[c for w in against for c in self.h(x, w, gram).list()] for x in vectors])
        kernel = matrix(ZZ, rows.denominator() * rows).left_kernel().basis()
        return [sum(c * x for c, x in zip(k, vectors)) for k in kernel]

    def _orthonormal_basis(self, vectors, gram):
        r"""Return an `O_K`-orthonormal basis of a definite span, found among trace-norm vectors, or ``None``."""
        result = []
        current = list(vectors)
        while current:
            trace = matrix(QQ, [[self.h(x, y, gram).trace() for y in current] for x in current])
            if not trace.is_positive_definite():
                return None
            unit = next(
                (
                    w
                    for c in _short_vectors(trace, self.degree)
                    for w in (sum(a * x for a, x in zip(c, current)),)
                    if self.h(w, w, gram) == 1
                ),
                None,
            )
            if unit is None:
                return None
            result.append(unit)
            current = self._orthogonal(current, [unit], gram)
        return result

    def _unitary_model(self):
        r"""Return ``(field, generators, form, basis, gram)`` of a recognized unitary model of `(\Lambda, h)`."""
        d = _IMAGINARY_MODELS.get(self.field.discriminant())
        if self.form_kind != "hermitian" or d is None:
            return None
        pair = self._hyperbolic_pair(self.gram)
        if pair is None:
            return None
        v, f, gram = pair
        if self.n == 2:
            if d not in (1, 2, 3):
                return None
            table = binary_hyperbolic_unitary.generators(d)
            return binary_hyperbolic_unitary.field(d), tuple(table.values()), binary_hyperbolic_unitary.form(d), [v, f], gram
        complement = self._orthogonal(self._lattice_vectors, [v, f], gram)
        trace = matrix(QQ, [[self.h(x, y, gram).trace() for y in complement] for x in complement])
        if trace.is_negative_definite():
            gram, f = -gram, -f
        units = self._orthonormal_basis(complement, gram)
        if units is None:
            return None
        if self.n == 3 and d == 1:
            generators = (*gaussian_picard.GENERATORS.values(), gaussian_picard.SCALAR)
            return gaussian_picard.FIELD, generators, gaussian_picard.FORM, [v, *units, -f], gram
        if self.n == 3 and d in (2, 7, 11):
            generators = euclidean_picard.full_unitary_generators(d)
            return euclidean_picard.field(d), generators, euclidean_picard.form(d), [v, *units, f], gram
        if d == 3 and self.n in (3, 4):
            m = self.n - 2
            generators = tuple(eisenstein_hyperbolic_unitary.generators(m).values())
            return eisenstein_hyperbolic_unitary.FIELD, generators, eisenstein_hyperbolic_unitary.form(m), [v, *units, f], gram
        if d == 3 and self.n == 5:
            field = act_eisenstein.FIELD
            scalar = -act_eisenstein.OMEGA * identity_matrix(field, 5)
            generators = (*act_eisenstein.FULL_GENERATORS, scalar)
            return field, generators, act_eisenstein.FORM, [*units, v, f], gram
        return None

    def _symplectic_model(self):
        r"""Return the `SL_2(O_K)` model of an alternating `h` on `\Lambda` of rank two over a real quadratic `K`.

        For a primitive `v` and `\phi = h(\cdot, v)`, `\phi(\Lambda) = j O_K` and an `f` with `\phi(f) = j` give
        `\Lambda = O_K v \oplus O_K f`; there `U(\Lambda, h) = \{g : g^t h g = h\} = SL_2(O_K)`. With `O_K`
        norm-Euclidean, the Euclidean reduction of :class:`HilbertWittOneLattice` writes `g \in SL_2(O_K)` as
        swaps and shears times `diag(u, u^{-1})`, `u \in \pm \varepsilon^Z`.
        """
        if self.form_kind != "alternating" or self.n != 2 or self.field.degree() != 2:
            return None
        D = ZZ(self.field.discriminant()).squarefree_part()
        if D < 0 or D not in _NORM_EUCLIDEAN_REAL:
            return None
        model = HilbertWittOneLattice(D)
        E = model.field
        one, zero = E.one(), E.zero()
        unit = E(E.unit_group().fundamental_units()[0])
        generators = (
            matrix(E, 2, 2, [zero, one, -one, zero]),
            *(matrix(E, 2, 2, [one, b, zero, one]) for b in model._beta_basis),
            matrix(E, 2, 2, [unit, zero, zero, ~unit]),
            -identity_matrix(E, 2),
        )
        gram = self.gram
        vectors = self._lattice_vectors
        v = self._primitive(vectors[0])
        values = [self.h(x, v, gram) for x in vectors]
        ideal = self.field.ideal(values).gens_reduced()
        assert len(ideal) == 1
        solution = self._integral_combination(values, ideal[0])
        f = sum(c * x for c, x in zip(solution, vectors))
        form = matrix(E, 2, 2, [zero, one, -one, zero])
        return E, generators, form, [v, f], gram / self.h(f, v, gram)

    def description(self):
        return (
            f"the centralizer of {self.commutant.isometries} is the unitary group of the {self.form_kind} form "
            f"{list(self.gram)} of rank {self.n} over O_K, K = Q[x]/({self.field.polynomial()}), with an involution "
            f"of the {self.kind} kind; no unitary model of the port matches it (#54)"
        )

    def centralizer_generators(self):
        r"""Return generators of `C(Z)`: model generators conjugated to `\Lambda`, then Schreier descent to `O`."""
        model = self._unitary_model() if self.kind == "second" else self._symplectic_model()
        if model is None:
            raise NotImplementedError(self.description())
        E, generators, form, basis, gram = model
        if not self._span_is_lattice(basis):
            raise NotImplementedError(f"{self.description()}: the model basis does not span Lambda")
        embedding = E.embeddings(self.field)[0]
        change = matrix(self.field, basis).transpose()
        assert self._star(change) * gram * change == form.apply_map(embedding)
        inverse = change.inverse()
        moves = [self.element_of(change * g.apply_map(embedding) * inverse) for g in generators]
        one = identity_matrix(QQ, self.commutant.rank)
        assert all(self.commutant.adjoint(m) * m == one for m in moves)
        return self.commutant._order_subgroup_generators(moves)


def _missing_group(commutant):
    r"""Return the description of the arithmetic group of a commutant with no implemented route."""
    center = commutant._center_basis()
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
        f"the transporter of {commutant.isometries} is an isometry of rank-{degree} forms over O_K, "
        f"K = Q[x]/({polynomial}), with an involution of the {kind} kind on M_{degree}(K) or a K-division "
        f"algebra; the unitary-model isometry decision is not implemented (#54)"
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
        return CentralFieldCommutant(commutant).centralizer_generators()
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
    if radical:
        return levi_transporter(commutant, target_gram, targets)
    model = _quaternion(commutant)
    if model is not None:
        return model.equivariant_transporter(target_gram, targets)
    commutative = all(x * y == y * x for x in commutant.order_basis for y in commutant.order_basis)
    if not radical and commutative and commutant.dimension() == 2 and len(targets) == 1:
        return commutant.equivariant_transporter(target_gram, targets[0])
    raise NotImplementedError(_missing_group(commutant))
