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
- `A / rad = \prod_k M_{n_k}(Q)` otherwise: the transporter is reduced, factor by factor along the adjoint
  orbits, to isometries of integral forms (or `GL` pairs) on the Levi modules and a coset search in the
  unipotent radical (:func:`levi_transporter`).
- `A` central simple over a number field `K \ne Q`: `A = End_K(W)` for a simple module `W`, the
  adjoint is the adjoint of a Hermitian, symmetric or alternating form `h` on `W`, and `C(Z)` is
  the subgroup of `U(\Lambda, h)` lying in `O`, `\Lambda = O_K O f` (:class:`CentralFieldCommutant`).
  Hyperbolic forms use the port's unitary models; an anisotropic binary Hermitian form over `Z[i]`
  uses the Bestvina-Savin ocean (:mod:`hermitian_ocean`), and over another imaginary quadratic `O_K` of
  class number one the quaternion order of `SU(h)` (:class:`_BinaryHermitianOrder`). The transporter reduces both forms to one
  model, composes the model isometry, and searches the coset representatives of `C(Z)` in `U(\Lambda, h)`.

  A symmetric form of rank two over `O_K` (first kind) is decided on the torus `K[c]^1` of the skew line
  `K c` and one reflection coset, through relative norm equations (:meth:`CentralFieldCommutant._binary_orthogonal_generators`);
  its transporter composes a `K`-isometry from a conic point with the same norm-equation search
  (:meth:`CentralFieldCommutant._binary_orthogonal_transporter`).

:class:`VectorOrbit` decides the orbit of a vector under `C(Z)` and returns its stabilizer, through the
marked transporter of the cyclic `F`-span of the vector.

Other noncompact commutants raise :class:`NotImplementedError` naming the missing group.
"""

from itertools import combinations, product

from sage.arith.functions import lcm
from sage.arith.misc import binomial, gcd
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.libs.pari import pari
from sage.matrix.constructor import block_matrix, matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.quadratic_forms.binary_qf import BinaryQF
from sage.quadratic_forms.qfsolve import qfsolve
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.schemes.plane_conics.constructor import Conic
from sage_extensions.algebras.isometry_commutant import IsometryCommutant, LeviQuotient, _integer_solution, _short_vectors
from sage_extensions.categories.lattice_centralizers import finite_preimage_generators
from sage_extensions.categories.lattices import Lattices
from sage_indefinite_port.indefinite import (
    act_eisenstein,
    binary_hyperbolic_unitary,
    eisenstein_hyperbolic_unitary,
    euclidean_picard,
    gaussian_picard,
    hermitian_ocean,
)
from sage_indefinite_port.indefinite.equivariant_markings import EquivariantMarking, _glue_action, _schreier
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
        right `O`-module, so `\bar M` is a left `O`-module, decided by :func:`_reduced_norm_point`.
        """
        eta = square[0, 0]
        assert square == eta * identity_matrix(QQ, self.commutant.rank)
        conjugates = [self.to_model(a).conjugate() for a in lattice]
        element = _reduced_norm_point(self.voronoi(), conjugates, 1 / eta)
        if element is None:
            return None
        return _verified(self.commutant, start * self.from_model(element.conjugate()), target_gram, targets)

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
        result = matrix(QQ, [[int(c) for c in isometry(v)] for v in left.module_generators()]).transpose()
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


def _factor_initial_blocks(frame, j, image, twisted):
    r"""Return representatives `\bar a_0` of the solutions on the orbit of factor `j`, modulo the right action of `H`.

    ``image`` lists the bases of `Y_k = (\bar X \bar M)_k`, ``twisted`` the Grams `G^h_k` on `\bar M_k`. For `\pi(j) = j`,
    `\bar a_j` is an isometry of `(\bar M_j, G_j)` onto a sublattice `Z \subseteq Y_j` of index `i` with
    `i^2 \det G^h|_{Y_j} = \det G_j`; one isometry per `Z`. For a pair `k = \pi(j) \ne j`, `\bar a_j = Y_j X` and
    `\bar a_k = Y_k W` with `X^t (Y_j^t G^h_j Y_k) W = G_j`; modulo `GL(\bar M_j)` the matrix `X` is a Hermite sublattice basis
    of index dividing `\det G_j / \det(Y_j^t G^h_j Y_k)`, and `W` is forced and must be integral.
    Returns pairs (block of `j`, block of `\pi(j)` or ``None``).
    """
    k = frame.partner[j]
    gram = frame.grams[j]
    restricted = image[j].transpose() * twisted[j] * image[k]
    ratio = QQ(gram.det()) / restricted.det()
    result = []
    if k == j:
        if ratio <= 0 or not ratio.is_square() or ratio.sqrt().denominator() != 1:
            return result
        for hermite in _sublattices(frame.degrees[j], ratio.sqrt()):
            sublattice = image[j] * hermite
            isometry = _form_isometry(gram, sublattice.transpose() * twisted[j] * sublattice)
            if isometry is not None:
                result.append((sublattice * isometry, None))
        return result
    if ratio.denominator() != 1:
        return result
    for index in ZZ(abs(ratio)).divisors():
        for hermite in _sublattices(frame.degrees[j], index):
            partner = restricted.inverse() * matrix(QQ, hermite).transpose().inverse() * gram
            if partner.denominator() == 1:
                result.append((image[j] * hermite, image[k] * partner))
    return result


def levi_transporter(commutant, target_gram, targets):
    r"""Return an integral transporter when `A / rad A = \prod_k M_{n_k}(Q)`, or ``None`` as a proof of absence.

    Every transporter is `t_0 a` with `a \in X = \{a : t_0 a \in M_r(Z)\}` and `a^\dagger h a = 1`, `h = B^{-1} t_0^t B' t_0`.
    With `Y = \bar X \bar M = \bigoplus_k Y_k` (:class:`LeviQuotient` of the provider), `\bar a = (\bar a_k)` maps each
    `\bar M_k` into `Y_k` and satisfies `\bar a_j^t G^h_j \bar a_{\pi(j)} = G_j`. These solutions are `\bar a_0 H` for the
    finitely many `\bar a_0` of :func:`_factor_initial_blocks`, one product over the adjoint orbits of factors. The integral
    solutions over `\bar a_0 g` are right `C(Z)`-invariant, so one representative of each left `\bar\Gamma`-coset of `H` is
    lifted by :meth:`LeviQuotient.lift`. A nonmaximal Levi order `\bar O` only enlarges the coset set.
    """
    start = rational_intertwiner(commutant.isometries, targets)
    if start is None:
        return None
    square = commutant.gram.change_ring(QQ).inverse() * start.transpose() * target_gram * start
    assert commutant.contains(square)
    frame = LeviQuotient(commutant)
    factors = range(len(frame.degrees))
    twisted = [frame.form(square, j) for j in factors]
    actions = [frame.on_module(x) for x in solution_lattice(commutant, start)]
    image = []
    for k in factors:
        o, n = frame.offsets[k], frame.degrees[k]
        generators = matrix(QQ, [column for action in actions for column in action[o : o + n, o : o + n].columns()])
        span = matrix(ZZ, generators.denominator() * generators).hermite_form(include_zero_rows=False)
        image.append((matrix(QQ, span) / generators.denominator()).transpose())
        assert image[-1].nrows() == image[-1].ncols() == n
    orbits = [j for j in factors if frame.partner[j] >= j]
    choices = [_factor_initial_blocks(frame, j, image, twisted) for j in orbits]
    for selection in product(*choices):
        blocks = [None] * len(frame.degrees)
        for j, (block, partner) in zip(orbits, selection):
            blocks[j] = block
            if partner is not None:
                blocks[frame.partner[j]] = partner
        initial = frame.block_diagonal(blocks)
        for representative in frame.coset_representatives():
            found = frame.lift(start, square, initial * representative)
            if found is not None:
                return _verified(commutant, matrix(QQ, found), target_gram, targets)
    return None


# Imaginary quadratic discriminants and the `d` of `Q(\sqrt{-d})` of the port's unitary models.
_IMAGINARY_MODELS = {-4: 1, -8: 2, -3: 3, -7: 7, -11: 11}


def _ideal_divisors(ideal):
    r"""Return the integral ideals dividing ``ideal``, the products of the prime powers of its factorization."""
    from itertools import product as cartesian

    factors = list(ideal.factor())
    one = ideal.number_field().ideal(1)
    result = []
    for exponents in cartesian(*[range(e + 1) for _p, e in factors]):
        divisor = one
        for (prime, _e), k in zip(factors, exponents):
            divisor = divisor * prime**k
        result.append(divisor)
    return result


def _reduced_norm_point(voronoi, elements, value):
    r"""Return an element of reduced norm ``value`` in the left `\Lambda`-module spanned by ``elements``, or ``None``.

    `\Lambda` is the order of ``voronoi``. By :meth:`QuaternionVoronoi.cell_trace_bound`, every `b` of the module
    with `n(b) = \nu` has a unit translate `ub` with `\mathrm{Tr}(\beta\beta^t)` at most `|\nu|` times the bound, `\beta`
    its real image, and `n(ub) = \pm\nu`. A unit of norm `-1`, if any, is among the unit generators. So the finite
    enumeration decides existence.
    """
    value = QQ(value)
    images = [voronoi.real_matrix(q) for q in elements]
    gram = matrix(RDF, [[float((x * y.transpose()).trace()) for y in images] for x in images])
    bound = float(abs(value) * voronoi.cell_trace_bound()) * (1 + 1e-9)
    found = pari(gram).qfminim(bound, None, 2)[2]
    columns = int(pari(found).matsize()[1])
    reversing = next((u for u in voronoi.unit_generators() if u.reduced_norm() == -1), None)
    for column in range(columns):
        coordinates = [int(round(float(found[row, column]))) for row in range(len(elements))]
        element = sum(c * q for c, q in zip(coordinates, elements))
        norm = element.reduced_norm()
        if norm == -value and reversing is not None:
            element, norm = reversing * element, value
        if norm == value:
            return element
    return None


def _hermitian_binary_isometry(source, target, sigma):
    r"""Return `P \in GL_2(K)` with `P^* T P = S` for anisotropic binary Hermitian ``source`` `S`, ``target`` `T`, or ``None``.

    `S = Q^* \mathrm{diag}(\alpha_1, \alpha_2) Q` with `\alpha_1 = S_{11}`, `\alpha_2 = \det S / S_{11}`. A vector `x` with
    `T(x) = \alpha_1` is a nonzero rational point of the quinary `Q`-form `\mathrm{Re}\, T(x) - \alpha_1 t^2` (PARI ``qfsolve``),
    with `t \ne 0` because `T` is anisotropic. The `T`-orthogonal line of `x` is spanned by `y`, and `\lambda y` has value
    `\alpha_2` exactly when `\alpha_2 / T(y)` is a norm from `K` (Witt cancellation; Scharlau, *Quadratic and Hermitian
    Forms*, Ch. 7, Thm. 6.3).
    """
    K = source.base_ring()
    star = lambda X: X.apply_map(sigma).transpose()
    value = lambda x, y: (star(matrix(K, [y]).transpose()) * target * matrix(K, [x]).transpose())[0, 0]
    a1 = QQ(source[0, 0])
    a2 = QQ(source.det() / source[0, 0])
    frame = matrix(K, 2, 2, [1, source[0, 1] / source[0, 0], 0, 1])
    assert star(frame) * matrix(K, 2, 2, [a1, 0, 0, a2]) * frame == source
    z = K.gen()
    basis = [vector(K, [1, 0]), vector(K, [z, 0]), vector(K, [0, 1]), vector(K, [0, z])]
    quinary = matrix(QQ, 5, 5, 0)
    for p, x in enumerate(basis):
        for q, y in enumerate(basis):
            quinary[p, q] = QQ(value(x, y).trace()) / 2
    quinary[4, 4] = -a1
    found = qfsolve(quinary)
    if found in ZZ:
        return None
    found = vector(QQ, found)
    assert found[4] != 0
    x = sum(found[p] / found[4] * basis[p] for p in range(4))
    assert value(x, x) == a1
    row = star(matrix(K, [x]).transpose()) * target
    y = vector(K, [-row[0, 1], row[0, 0]])
    assert value(y, x) == 0
    is_norm, scale = QQ(a2 / QQ(value(y, y))).is_norm(K, element=True)
    if not is_norm:
        return None
    result = matrix(K, [x, K(scale) * y]).transpose() * frame
    assert star(result) * target * result == source
    return result


class _BinaryHermitianOrder:
    r"""The quaternion order of `SU(\Lambda, S)` for an anisotropic binary Hermitian `S` over `O_K`, `K` imaginary quadratic.

    `\Lambda = O_K^2` and `X^\dagger = S^{-1} X^* S`, a unitary involution of `M_2(K)`. There is a unique quaternion
    `Q`-subalgebra `B` with `\dagger = \mathrm{adj} \otimes \sigma`, namely `B = \{X : X^\dagger = \mathrm{adj}\, X\}`
    (Knus-Merkurjev-Rost-Tignol, *The Book of Involutions*, Prop. 2.22). On `B`, `X^\dagger X = \det X` is the reduced
    norm, so `SU(S) = B^1` and `SU(\Lambda, S) = (B \cap M_2(O_K))^1`, generated by the Voronoi algorithm of
    :class:`QuaternionVoronoi` for the order `B \cap M_2(O_K)` (Braun-Coulangeon-Nebe-Schoennenbeck).

    `\det` maps `U(S)` onto a subgroup of `K^1`. For `\zeta \in K^1`, `\zeta = \mu / \bar\mu` with `\mu = 1 + \zeta`
    (or `\mu = z - \bar z` when `\zeta = -1`), and `g` with `\det g = \zeta` lies in `U(S)` exactly when `g = \mu Z`,
    `Z \in B`, `n(Z) = 1 / N(\mu)`: `(\mu Z)^* S (\mu Z) = N(\mu) n(Z) S` and `\det(\mu Z) = \mu^2 n(Z)`. The integral `Z`
    form the right `B \cap M_2(O_K)`-module `\{Z : \mu Z \in M_2(O_K)\}`, decided by :func:`_reduced_norm_point`.
    `U(\Lambda, S)` has `\det` in the roots of unity, so one such `g` per root of unity of the image completes `SU`.
    """

    def __init__(self, field, sigma, form):
        self.field, self.sigma, self.form = field, sigma, form
        identity = identity_matrix(field, 2)
        star = lambda X: X.apply_map(sigma).transpose()
        self._star = star
        inverse = form.inverse()
        self._integral_rows = matrix(QQ, [k.list() for k in field.integral_basis()])
        cells = [
            matrix(field, 2, 2, {(i, j): k}) for i in range(2) for j in range(2) for k in field.integral_basis()
        ]
        equation = lambda X: inverse * star(X) * form - X.adjugate()
        rows = matrix(QQ, [self._flat(equation(X)) for X in cells])
        kernel = matrix(ZZ, rows.denominator() * rows).left_kernel().basis()
        order = [sum(c * X for c, X in zip(k, cells)) for k in kernel]
        assert len(order) == 4
        pure = lambda X: X - X.trace() / 2 * identity
        first = next(pure(X) for X in order if not pure(X).is_zero())
        first = first * QQ((first * first)[0, 0]).denominator()
        second = next(c for c in (first * pure(X) - pure(X) * first for X in order) if not c.is_zero())
        second = second * QQ((second * second)[0, 0]).denominator()
        a, b = QQ((first * first)[0, 0]), QQ((second * second)[0, 0])
        assert first * first == a * identity and second * second == b * identity and first * second == -second * first
        if a > 0 and b > 0:
            first, a = first * second, -a * b
        elif a > 0 > b:
            first, second, a, b = second, first, b, a
        assert a < 0 < b, f"SU of the Hermitian form {list(form)} is compact: B = ({a}, {b})_Q"
        self.algebra = QuaternionAlgebra(QQ, a, b)
        self.frame = (identity, first, second, first * second)
        self._frame_rows = matrix(QQ, [self._flat(X) for X in self.frame])
        order_rows = matrix(QQ, [list(self.to_model(X)) for X in order])
        self.voronoi = QuaternionVoronoi(int(ZZ(a)), int(ZZ(b)), order_rows)

    def _flat(self, X):
        return [c for e in X.list() for c in self.field(e).list()]

    def to_model(self, X):
        return self.algebra(list(self._frame_rows.solve_left(vector(QQ, self._flat(X)))))

    def from_model(self, quaternion):
        return sum(QQ(c) * X for c, X in zip(list(quaternion), self.frame))

    def _integral_coordinates(self, X):
        return [c for e in X.list() for c in self._integral_rows.solve_left(vector(QQ, self.field(e).list()))]

    def module(self, prefix):
        r"""Return a `Z`-basis, as quaternions, of the right module `\{Z \in B : \text{prefix}\, Z \in M_2(O_K)\}`."""
        rows = matrix(QQ, [self._integral_coordinates(prefix * X) for X in self.frame])
        saturated = matrix(ZZ, rows.denominator() * rows).saturation()
        return [self.algebra(list(rows.solve_left(vector(QQ, row)))) for row in saturated.rows()]

    def _hilbert_ninety(self, zeta):
        z = self.field.gen()
        mu = 1 + zeta if zeta != -1 else z - self.sigma(z)
        assert mu / self.sigma(mu) == zeta
        return mu

    def _point(self, prefix, mu):
        r"""Return `\mu\, \text{prefix}\, Z` integral with `n(Z) = 1 / N(\mu)`, `Z \in B`, or ``None``."""
        conjugates = [q.conjugate() for q in self.module(mu * prefix)]
        found = _reduced_norm_point(self.voronoi, conjugates, ~QQ(mu.norm()))
        return None if found is None else mu * prefix * self.from_model(found.conjugate())

    def unitary_generators(self):
        r"""Return generators of `U(\Lambda, S)`: `SU(\Lambda, S)` and one element of each determinant in the image."""
        if hasattr(self, "_generators"):
            return self._generators
        norm_one, central = self.voronoi.norm_one_generators()
        result = [self.from_model(u) for u in (*norm_one, *central)]
        assert all(X.det() == 1 for X in result)
        identity = identity_matrix(self.field, 2)
        for zeta in self.field.roots_of_unity():
            if zeta == 1:
                continue
            found = self._point(identity, self._hilbert_ninety(zeta))
            if found is not None:
                assert found.det() == zeta
                result.append(found)
        assert all(self._star(X) * self.form * X == self.form for X in result)
        self._generators = tuple(result)
        return self._generators

    def isometry_from(self, target):
        r"""Return `D \in GL_2(O_K)` with `D^* T D = S` for ``target`` `T`, or ``None`` as a proof of absence.

        Every `K`-isometry is `P_0 g`, `g \in U(S)`, with `P_0` from :func:`_hermitian_binary_isometry`. An integral one
        has `\det D` a root of unity `\zeta'` (the units of `O_K`), so `\det g = \zeta = \zeta' / \det P_0` has norm one,
        and `g = \mu Z` as in the class description, with `\mu P_0 Z` integral; this right module is searched.
        """
        start = _hermitian_binary_isometry(self.form, target, self.sigma)
        if start is None:
            return None
        delta = start.det()
        for root in self.field.roots_of_unity():
            zeta = root / delta
            if zeta * self.sigma(zeta) != 1:
                return None
            found = self._point(start, self._hilbert_ninety(zeta))
            if found is not None:
                assert self._star(found) * target * found == self.form and found.det() == root
                return found
        return None


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
    For a self-adjoint `m \in A` the form `h_m` is read from `y^\dagger m x` in the same way.
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
        self._f = f
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
        self.form_kind, self._scale = self._normalization(self._raw_gram(identity_matrix(QQ, r)))
        self.gram = self.gram_of(identity_matrix(QQ, r))
        self._lattice = self._lattice_rows(order)

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
            f"division algebra D over K = Q[x]/({self.field.polynomial()}), or the box missed its zero divisors. "
            f"Its centralizer is the unitary group of a form over an order of D, not SL_n of that order, so elementary "
            f"generation (Bass-Milnor-Serre; Vaserstein) does not apply; the missing statement is a Voronoi generation "
            f"of the unitary group of a D-form (Braun-Coulangeon-Nebe-Schoennenbeck, J. Algebra 435 (2015), for "
            f"M_n(D)^x), not implemented (#54)"
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

    def _lattice_rows(self, elements):
        r"""Return the Hermite basis, in the `Q`-coordinates of `W`, of `O_K X f` for ``elements`` `= X`."""
        integers = self.field.ring_of_integers().basis()
        rows = matrix(QQ, [self._q_coordinates(self.embed(k) * m * self._f) for k in integers for m in elements])
        hermite = matrix(ZZ, rows.denominator() * rows).hermite_form(include_zero_rows=False)
        return matrix(QQ, hermite) / rows.denominator()

    def _vectors(self, lattice):
        return tuple(self._k_vector(row) for row in lattice.rows())

    def _pairing(self, x, y, middle):
        coordinates = self._line_rows.solve_left(vector(QQ, (self.commutant.adjoint(y) * middle * x).list()))
        return self.field(list(coordinates))

    def _raw_gram(self, middle):
        return matrix(self.field, [[self._pairing(b, c, middle) for b in self._basis] for c in self._basis])

    def gram_of(self, middle):
        r"""Return the Gram `G_{ij} = h_m(b_j, b_i)` of the normalized form `h_m`, `m =` ``middle``."""
        return self._scale * self._raw_gram(middle)

    def _star(self, block):
        return block.apply_map(self.sigma).transpose()

    def _normalization(self, gram):
        r"""Return the kind of `h` and the scalar making it (skew-)Hermitian; `G_{ij} = h(b_j, b_i)`."""
        i, j = next((i, j) for i in range(self.n) for j in range(self.n) if gram[i, j] != 0)
        ratio = self._star(gram)[i, j] / gram[i, j]
        assert self._star(gram) == ratio * gram
        if self.kind == "first":
            assert ratio in (1, -1)
            return ("symmetric" if ratio == 1 else "alternating"), self.field.one()
        z = self.field.gen()
        scalar = 1 + ratio if ratio != -1 else z - self.sigma(z)
        assert self._star(scalar * gram) == scalar * gram
        return "hermitian", scalar

    def h(self, x, y, gram):
        return (self._star(matrix(self.field, [y]).transpose()) * gram * matrix(self.field, [x]).transpose())[0, 0]

    def _span_is_lattice(self, vectors, lattice):
        integers = self.field.ring_of_integers().basis()
        rows = matrix(QQ, [self._q_vector(k * v) for k in integers for v in vectors])
        hermite = matrix(ZZ, rows.denominator() * rows).hermite_form(include_zero_rows=False)
        return rows.rank() == lattice.nrows() and matrix(QQ, hermite) / rows.denominator() == lattice

    def _integral_combination(self, values, target):
        r"""Return `x \in Z^N` with `\sum x_a` ``values[a]`` `=` ``target`` in `K`, or ``None``."""
        rows = matrix(QQ, [self.field(v).list() for v in values] + [self.field(target).list()])
        scale = rows.denominator()
        equations = matrix(ZZ, scale * rows[:-1]).transpose()
        return _integer_solution(equations, vector(ZZ, scale * rows[-1]))

    def _primitive(self, v, lattice):
        r"""Return `\gamma v` with `I = \gamma O_K` the fractional ideal `\{k : k v \in \Lambda\}`."""
        rows = matrix(QQ, [lattice.solve_left(self._q_vector(p * v)) for p in self.field.power_basis()])
        span = matrix(ZZ, rows.denominator() * rows.transpose()).hermite_form(include_zero_rows=False)
        dual = (matrix(QQ, span) / rows.denominator()).inverse().transpose()
        ideal = self.field.ideal([self.field(list(row)) for row in dual.rows()])
        generator = ideal.gens_reduced()
        assert len(generator) == 1
        return generator[0] * v

    def _isotropic_vector(self, lattice, gram):
        r"""Return a primitive null vector of `\Lambda`, or ``None`` when `h` is anisotropic.

        Isotropy of `h` is isotropy of the trace form `Tr\, h(x, x) = 2 h(x, x)` (PARI ``qfsolve``).
        """
        vectors = self._vectors(lattice)
        trace = matrix(QQ, [[self.h(x, y, gram).trace() for y in vectors] for x in vectors])
        found = qfsolve(trace)
        if found in ZZ:
            return None
        found = vector(QQ, found)
        found = found * found.denominator()
        found = found / gcd(list(found))
        v = self._primitive(sum(c * x for c, x in zip(found, vectors)), lattice)
        assert self.h(v, v, gram) == 0
        return v

    def _hyperbolic_pair(self, lattice, gram, v):
        r"""Return `(f, G')` with `f` null, `h'(f, v) = 1`, `h'(\Lambda, v) = O_K`, `h' = h / q`, or ``None``.

        `f_1` with `h'(f_1, v) = 1` is corrected to `f = f_1 - t v`, `t + \bar t = h'(f_1, f_1)`.
        """
        vectors = self._vectors(lattice)
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
        return f, gram

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

    def _unitary_model(self, lattice, gram):
        r"""Return ``(kind, field, generators, form, basis, gram)`` of a recognized unitary model of `(\Lambda, h)`."""
        if self.form_kind != "hermitian":
            return None
        imaginary = self.degree == 2 and self.field.discriminant() < 0
        if self.n == 2 and imaginary and self._isotropic_vector(lattice, gram) is None:
            if self.field.discriminant() == -4:
                return self._ocean_model(lattice, gram)
            return self._quaternion_model(lattice, gram)
        d = _IMAGINARY_MODELS.get(self.field.discriminant())
        if d is None:
            return None
        v = self._isotropic_vector(lattice, gram)
        if v is None:
            return None
        pair = self._hyperbolic_pair(lattice, gram, v)
        if pair is None:
            return None
        f, gram = pair
        if self.n == 2:
            if d not in (1, 2, 3):
                return None
            table = binary_hyperbolic_unitary.generators(d)
            E = binary_hyperbolic_unitary.field(d)
            return "hyperbolic", E, tuple(table.values()), binary_hyperbolic_unitary.form(d), [v, f], gram
        complement = self._orthogonal(self._vectors(lattice), [v, f], gram)
        trace = matrix(QQ, [[self.h(x, y, gram).trace() for y in complement] for x in complement])
        if trace.is_negative_definite():
            gram, f = -gram, -f
        units = self._orthonormal_basis(complement, gram)
        if units is None:
            return None
        if self.n == 3 and d == 1:
            generators = (*gaussian_picard.GENERATORS.values(), gaussian_picard.SCALAR)
            return "hyperbolic", gaussian_picard.FIELD, generators, gaussian_picard.FORM, [v, *units, -f], gram
        if self.n == 3 and d in (2, 7, 11):
            generators = euclidean_picard.full_unitary_generators(d)
            return "hyperbolic", euclidean_picard.field(d), generators, euclidean_picard.form(d), [v, *units, f], gram
        if d == 3 and self.n in (3, 4):
            m = self.n - 2
            generators = tuple(eisenstein_hyperbolic_unitary.generators(m).values())
            E = eisenstein_hyperbolic_unitary.FIELD
            return "hyperbolic", E, generators, eisenstein_hyperbolic_unitary.form(m), [v, *units, f], gram
        if d == 3 and self.n == 5:
            field = act_eisenstein.FIELD
            scalar = -act_eisenstein.OMEGA * identity_matrix(field, 5)
            generators = (*act_eisenstein.FULL_GENERATORS, scalar)
            return "hyperbolic", field, generators, act_eisenstein.FORM, [*units, v, f], gram
        return None

    def _binary_basis(self, lattice):
        r"""Return an `O_K`-basis `(\gamma_1 e_1, w)` of a rank-two `\Lambda` over a principal `O_K`.

        The second coordinate maps `\Lambda` onto an ideal `\gamma_2 O_K`, and `w` is a preimage of `\gamma_2`;
        the kernel is `\Lambda \cap K e_1 = \gamma_1 O_K e_1`. So `\Lambda = O_K \gamma_1 e_1 \oplus O_K w`.
        """
        vectors = self._vectors(lattice)
        values = [x[1] for x in vectors]
        second = self.field.ideal(values).gens_reduced()
        assert len(second) == 1
        w = sum(c * x for c, x in zip(self._integral_combination(values, second[0]), vectors))
        rows = matrix(QQ, [self.field(value).list() for value in values])
        kernel = matrix(ZZ, rows.denominator() * rows).left_kernel().basis()
        first = self.field.ideal([sum(c * x for c, x in zip(k, vectors))[0] for k in kernel]).gens_reduced()
        assert len(first) == 1
        return [vector(self.field, [first[0], 0]), w]

    def _pullback(self, block, embedding):
        r"""Return the preimage of ``block`` under the isomorphism ``embedding`` of a quadratic field onto `K`."""
        E = embedding.domain()
        frame = matrix(QQ, [self.field.one().list(), embedding(E.gen()).list()])
        return block.apply_map(lambda x: E(list(frame.solve_left(vector(QQ, self.field(x).list())))), E)

    def _ocean_model(self, lattice, gram):
        r"""Return the Bestvina-Savin model of an anisotropic `h` of rank two over `Z[i]`.

        In an `O_K`-basis `C` of `\Lambda` with `S = C^* G C`, `g` preserves `S` exactly when `\bar g` preserves
        `\bar S`, that is when `g^T H \bar g = H` for `H = q \bar S` and any `q \in Q^\times`; `q > 0` makes `H`
        integral in the sense of :mod:`hermitian_ocean`, whose generators of `U(f)` are then those of `U(\Lambda, h)`.
        """
        basis = self._binary_basis(lattice)
        E = hermitian_ocean.FIELD
        embedding = E.embeddings(self.field)[0]
        change = matrix(self.field, basis).transpose()
        form = self._pullback(self._star(change) * gram * change, embedding)
        ocean = hermitian_ocean.HermitianOcean(_ocean_forms(form)[0])
        return "ocean", E, ocean.unitary_generators(), form, basis, gram

    def _binary_hermitian_order(self, form):
        key = tuple(form.list())
        orders = self.__dict__.setdefault("_hermitian_orders", {})
        if key not in orders:
            orders[key] = _BinaryHermitianOrder(self.field, self.sigma, form)
        return orders[key]

    def _quaternion_model(self, lattice, gram):
        r"""Return the quaternion model of an anisotropic binary Hermitian `h` over `O_K`, `K` imaginary quadratic.

        In an `O_K`-basis `C` of `\Lambda` with `S = C^* G C`, `U(\Lambda, h) = U(O_K^2, S)`, generated through the
        quaternion order of :class:`_BinaryHermitianOrder`.
        """
        if self.field.class_number() != 1:
            raise NotImplementedError(
                f"{self.description()}: O_K has class number {self.field.class_number()}, so Lambda need not be free, "
                f"and the quaternion order of SU(h) is built on an O_K-basis of Lambda (#54)"
            )
        basis = self._binary_basis(lattice)
        change = matrix(self.field, basis).transpose()
        form = self._star(change) * gram * change
        generators = self._binary_hermitian_order(form).unitary_generators()
        return "quaternion", self.field, generators, form, basis, gram

    def _model_embedding(self, E):
        return self.field.hom([self.field.gen()]) if E is self.field else E.embeddings(self.field)[0]

    def _symplectic_model(self, lattice, gram):
        r"""Return the `SL_2(O_K)` model of an alternating `h` on `\Lambda` of rank two over a real quadratic `K`.

        For a primitive `v` and `\phi = h(\cdot, v)`, `\phi(\Lambda) = j O_K` and an `f` with `\phi(f) = j` give
        `\Lambda = O_K v \oplus O_K f` when `j O_K` is principal, which class number one guarantees; there
        `U(\Lambda, h) = \{g : g^t h g = h\} = SL_2(O_K)`. `O_K^*` is infinite, so `SL_2(O_K)` is generated by the
        swap and the shears `E_{12}(b)`, `b` in a `Z`-basis of `O_K` (Vaserstein, Math. USSR-Sb. 18 (1972)
        321-332; see :mod:`witt_one_hilbert`).
        """
        if self.form_kind != "alternating" or self.n != 2 or self.field.degree() != 2:
            return None
        D = ZZ(self.field.discriminant()).squarefree_part()
        if D < 0 or self.field.class_number() != 1:
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
        vectors = self._vectors(lattice)
        v = self._primitive(vectors[0], lattice)
        values = [self.h(x, v, gram) for x in vectors]
        ideal = self.field.ideal(values).gens_reduced()
        assert len(ideal) == 1
        solution = self._integral_combination(values, ideal[0])
        f = sum(c * x for c, x in zip(solution, vectors))
        form = matrix(E, 2, 2, [zero, one, -one, zero])
        return "alternating", E, generators, form, [v, f], gram / self.h(f, v, gram)

    def description(self):
        return (
            f"the centralizer of {self.commutant.isometries} is the unitary group of the {self.form_kind} form "
            f"{list(self.gram)} of rank {self.n} over O_K, K = Q[x]/({self.field.polynomial()}), with an involution "
            f"of the {self.kind} kind; no unitary model of the port matches it (#54)"
        )

    def model(self, lattice, gram):
        r"""Return ``(kind, field, generators, form, basis, gram')`` of `(\Lambda, h)`, `\Lambda =` ``lattice``.

        The basis spans `\Lambda`, and in it `G' = s G` is the model ``form``, `s \in K^\times`.
        """
        model = self._unitary_model(lattice, gram) if self.kind == "second" else self._symplectic_model(lattice, gram)
        if model is None and self.kind == "first" and self.form_kind == "symmetric":
            raise NotImplementedError(
                f"{self.description()}: O(Lambda, h) of a symmetric O_K-lattice of rank {self.n} >= 3 needs (i) a decision "
                f"whether Lambda splits a hyperbolic O_K-plane H, (ii) generation of O(Lambda) by the Eichler transformations "
                f"E(e, u), O of the complement and spinor-norm and determinant coset representatives, which is proved when "
                f"Lambda = H + H + L' (Hahn-O'Meara, The Classical Groups and K-Theory, Ch. 5-6), not for one hyperbolic "
                f"plane, and (iii) the Eichler criterion for orbits of primitive vectors (same hypothesis) for the "
                f"transporter recursion; none of these is implemented (#54)"
            )
        if model is None:
            raise NotImplementedError(self.description())
        _kind, E, _generators, form, basis, model_gram = model
        if not self._span_is_lattice(basis, lattice):
            raise NotImplementedError(f"{self.description()}: the model basis does not span Lambda")
        change = matrix(self.field, basis).transpose()
        assert self._star(change) * model_gram * change == form.apply_map(self._model_embedding(E))
        return model

    def unitary_generators(self):
        r"""Return the model generators of `U(\Lambda, h)`, conjugated to elements of `A`."""
        if not hasattr(self, "_unitary"):
            _kind, E, generators, _form, basis, _gram = self.model(self._lattice, self.gram)
            embedding = self._model_embedding(E)
            change = matrix(self.field, basis).transpose()
            inverse = change.inverse()
            moves = [self.element_of(change * g.apply_map(embedding) * inverse) for g in generators]
            one = identity_matrix(QQ, self.commutant.rank)
            assert all(self.commutant.adjoint(m) * m == one for m in moves)
            self._unitary = tuple(moves)
        return self._unitary

    def centralizer_generators(self):
        r"""Return generators of `C(Z)`: model generators conjugated to `\Lambda`, then Schreier descent to `O`."""
        if self.kind == "first" and self.form_kind == "symmetric" and self.n == 2:
            return self._binary_orthogonal_generators()
        return self.commutant._order_subgroup_generators(list(self.unitary_generators()))

    def _to_field(self, element):
        r"""Return the element of `K` acting as the central ``element``."""
        rows = matrix(QQ, [p.list() for p in self._powers])
        return self.field(list(rows.solve_left(vector(QQ, matrix(QQ, element).list()))))

    def _left_transversal(self, moves):
        r"""Return representatives of the left cosets of `\langle` ``moves`` `\rangle \cap O` in `\langle` ``moves`` `\rangle`."""
        representatives = [identity_matrix(QQ, self.commutant.rank)]
        position = 0
        while position < len(representatives):
            current = representatives[position]
            position += 1
            for move in moves:
                moved = move * current
                if not any(self.commutant.in_order(known.inverse() * moved) for known in representatives):
                    representatives.append(moved)
        return representatives

    def _binary_orthogonal_generators(self):
        r"""Return generators of `C(Z) = O(\Lambda, h) \cap O` for a symmetric `h` of rank two over `O_K`.

        The skew elements of `A = M_2(K)` form a `K`-line `K c` with `c^2 = \delta \in K^\times`, and for a binary
        form `g c g^{-1} = \det(g) c`. So `SO(h) = L^1` for the commutative algebra `L = K[c]` with the involution
        `c \mapsto -c`, and `O(h) \setminus SO(h)` is the set of `w` with `wc = -cw`, `w^\dagger w = 1` (Knus-Merkurjev-Rost-
        Tignol, *The Book of Involutions*, Prop. 7.4 for `n = 2`). `C(Z) \cap SO` is the subgroup in `O` of the norm-one
        units of the maximal order of `L` (Schreier). Every anticommuting `w` is `x w_1` for one fixed `w_1` and `x \in L`,
        and `w^\dagger w = N_{L/K}(x)\, q`, `q = w_1^\dagger w_1 \in K`. With `M = \{x : x w_1 \in O\}` and `m M \subseteq O_L`,
        the solutions `x \in M` of `N(x) = 1/q` are `y u / m` with `y` among the finitely many elements of `O_L` of absolute
        norm `N_{K/Q}(m^2/q)` modulo units, and `u` a unit of relative norm `m^2 / (q N(y))`, unique up to the norm-one
        units `E`. `M` is stable under `E \cap O`, so a transversal of `E / (E \cap O)` finishes the decision.
        """
        torus = self._binary_torus()
        one = identity_matrix(QQ, self.commutant.rank)
        rotations = self.commutant._order_subgroup_generators(torus["moves"])
        reflection = next(
            (
                w
                for x in self._norm_points(one, torus["w1"], ~torus["q"])
                for t in torus["transversal"]
                for w in (x * t * torus["w1"],)
                if self.commutant.is_centralizer_element(w)
            ),
            None,
        )
        result = [matrix(ZZ, g) for g in rotations]
        if reflection is not None:
            result.append(matrix(ZZ, reflection))
        assert all(self.commutant.is_centralizer_element(g) for g in result)
        return tuple(result)

    def _binary_torus(self):
        r"""Return the torus data of :meth:`_binary_orthogonal_generators`: `c`, the `Q`-basis of `L = K[c]`, the
        generators ``moves`` of `\pm` the norm-one units of the maximal order of `L` as elements of `A`, a left
        transversal of `\langle` ``moves`` `\rangle \cap O`, an anticommuting `w_1 \in O` and `q = w_1^\dagger w_1`."""
        if hasattr(self, "_torus"):
            return self._torus
        commutant = self.commutant
        one = identity_matrix(QQ, commutant.rank)
        c = _integral_kernel(self._order, lambda x: x + commutant.adjoint(x))[0]
        delta = self._to_field(c * c)
        K = self.field
        assert delta != 0
        torus = {"c": c, "delta": delta, "basis": [*self._powers, *(p * c for p in self._powers)]}
        if delta.is_square():
            s = delta.sqrt()
            e = (one + self.embed(~s) * c) / 2
            split = lambda x1, x2: self.embed(x1) * e + self.embed(x2) * (one - e)
            torus.update(s=s, split=split)
            moves = [split(K(u), ~K(u)) for u in K.unit_group().gens()]
        else:
            L = K.extension(PolynomialRing(K, "t")([-delta, 0, 1]), "c")
            M = L.absolute_field("b")
            from_abs, to_abs = M.structure()
            to_matrix = lambda x: sum(self.embed(k) * c**j for j, k in enumerate(from_abs(x).list()))
            relative = from_abs(M.gen()).list()
            theta = to_abs(L(relative[0]) - L(relative[1]) * L.gen())
            norm_one = IsometryCommutant._norm_one_generators(M, theta, M.polynomial())
            torus.update(L=L, M=M, from_abs=from_abs, to_abs=to_abs, to_matrix=to_matrix)
            moves = [to_matrix(u) for u in norm_one]
        torus["moves"] = [-one, *moves]
        anticommuting = _integral_kernel(self._order, lambda x: x * c + c * x)
        rank = commutant.rank
        w1 = next(
            w
            for values in product(range(rank + 1), repeat=len(anticommuting))
            for w in (sum(a * x for a, x in zip(values, anticommuting)),)
            if w.det() != 0
        )
        torus["w1"] = w1
        torus["q"] = self._to_field(commutant.adjoint(w1) * w1)
        torus["transversal"] = self._left_transversal(torus["moves"])
        self._torus = torus
        return torus

    def _norm_points(self, prefix, suffix, nu):
        r"""Return `x_1, \ldots \in L` with `N_{L/K}(x_i) = \nu` such that every `x \in L` with `N(x) = \nu` and
        `\text{prefix}\, x\, \text{suffix} \in M_r(Z)` is some `x_i \varepsilon`, `\varepsilon \in \langle` ``moves`` `\rangle`.

        `M = \{x : \text{prefix}\, x\, \text{suffix} \in M_r(Z)\}` is a lattice with `m M \subseteq O_L`, and the solutions are
        `y u / m` for `y \in O_L` of absolute norm `N_{K/Q}(m^2 \nu)` modulo units and units `u` of relative norm
        `m^2 \nu / N(y)` (see :meth:`_binary_orthogonal_generators`).
        """
        torus = self._binary_torus()
        K = self.field
        rows = matrix(QQ, [(prefix * x * suffix).list() for x in torus["basis"]])
        scale = rows.denominator()
        saturated = matrix(ZZ, scale * rows).saturation()
        span = [rows.solve_left(vector(QQ, row)) for row in saturated.rows()]
        d = self.degree
        pairs = [(K(list(x[:d])), K(list(x[d:]))) for x in span]
        integral_rows = matrix(QQ, [k.list() for k in K.integral_basis()])
        denominator = lambda k: integral_rows.solve_left(vector(QQ, K(k).list())).denominator()
        candidates = []
        if "split" in torus:
            s, split = torus["s"], torus["split"]
            m = lcm([denominator(a + b * s) for a, b in pairs] + [denominator(a - b * s) for a, b in pairs])
            target = m**2 * nu
            if target.is_integral():
                for divisor in _ideal_divisors(K.ideal(target)):
                    if divisor.is_principal():
                        y1 = divisor.gens_reduced()[0]
                        candidates.append(split(y1 / m, target / y1 / m))
        else:
            L, M, from_abs, to_abs = torus["L"], torus["M"], torus["from_abs"], torus["to_abs"]
            abs_rows = matrix(QQ, [b.list() for b in M.integral_basis()])
            abs_denominator = lambda x: abs_rows.solve_left(vector(QQ, x.list())).denominator()
            m = lcm([abs_denominator(to_abs(L(a) + L(b) * L.gen())) for a, b in pairs])
            target = m**2 * nu
            absolute = target.norm()
            relative_units = K.unit_group()
            units = M.unit_group()
            values = [M(u) for u in units.gens()]
            torsion = relative_units.torsion_generator().order()
            logs = [list(relative_units.log(from_abs(u).relative_norm())) for u in values]
            equations = matrix(ZZ, logs).transpose().augment(
                matrix(ZZ, [[torsion]] + [[0]] * (len(logs[0]) - 1))
            )
            if absolute in ZZ:
                for y in M.elements_of_norm(ZZ(absolute)):
                    tau = target / from_abs(y).relative_norm()
                    if not (tau.is_integral() and (~tau).is_integral()):
                        continue
                    exponents = _integer_solution(equations, vector(ZZ, list(relative_units.log(tau))))
                    if exponents is None:
                        continue
                    u = M.one()
                    for value, exponent in zip(values, exponents):
                        u = u * value ** int(exponent)
                    candidates.append(torus["to_matrix"](y * u / m))
        assert all(self._to_field(self.commutant.adjoint(x) * x) == nu for x in candidates)
        return candidates

    def _binary_orthogonal_transporter(self, target_gram, targets):
        r"""Return an integral transporter for a symmetric `h` of rank two over `O_K`, or ``None`` as a proof of absence.

        Every rational transporter is `t_0 a` with `a^\dagger h_0 a = 1`, `h_0 = B^{-1} t_0^t B' t_0`. On `W` this says
        `P^t G' P = G` for the Grams `G, G'` of `h, h_{h_0}`; :func:`_binary_isometry` decides it over `K` and returns one
        `a_0`. Then `a = a_0 g` with `g \in O(h)(K) = L^1 \cup L^1 w_1` (:meth:`_binary_orthogonal_generators`): `g = x` with
        `N(x) = 1` or `g = x w_1` with `N(x) = 1 / q`. The integral ones are found by :meth:`_norm_points` with prefix
        `t_0 a_0`; the set `M` there is stable under right multiplication by `\langle` ``moves`` `\rangle \cap O \subset C(Z)`, so the
        left transversal of that subgroup completes the search.
        """
        commutant = self.commutant
        start = rational_intertwiner(commutant.isometries, targets)
        if start is None:
            return None
        square = commutant.gram.change_ring(QQ).inverse() * start.transpose() * target_gram * start
        assert commutant.contains(square)
        block = _binary_isometry(self.gram, self.gram_of(square))
        if block is None:
            return None
        one = identity_matrix(QQ, commutant.rank)
        initial = self.element_of(block)
        assert commutant.adjoint(initial) * square * initial == one
        prefix = start * initial
        torus = self._binary_torus()
        for suffix, nu in ((one, self.field.one()), (torus["w1"], ~torus["q"])):
            for x in self._norm_points(prefix, suffix, nu):
                for t in torus["transversal"]:
                    found = _verified(commutant, prefix * x * t * suffix, target_gram, targets)
                    if found is not None:
                        return found
        return None

    def coset_representatives(self):
        r"""Return representatives `u` of the left cosets `u C(Z)` in `U(\Lambda, h)`.

        `C(Z) = U(\Lambda, h) \cap O` has finite index, and `x C(Z) = y C(Z)` exactly when `y^{-1} x \in O`. Left
        multiplication by the generators permutes the finite coset set transitively, so the orbit of `C(Z)`
        is a Schreier transversal.
        """
        if not hasattr(self, "_cosets"):
            representatives = [identity_matrix(QQ, self.commutant.rank)]
            position = 0
            while position < len(representatives):
                current = representatives[position]
                position += 1
                for generator in self.unitary_generators():
                    moved = generator * current
                    if not any(self.commutant.in_order(known.inverse() * moved) for known in representatives):
                        representatives.append(moved)
            self._cosets = tuple(representatives)
        return self._cosets

    def _model_isometry(self, source, target, source_gram, target_gram):
        r"""Return `D \in GL_n(O_K)` with `D^* T D = S` for the Grams `S, T` of the two model bases, or ``None``.

        For the ocean, `D` is the image of :func:`hermitian_ocean.isometry`. Otherwise `S = \rho T` with `\rho \in K`
        since both are scalings of one model form. Then `\rho = 1` gives `D = 1`; `\rho = -1` on the hyperbolic
        plane gives `D = diag(1, -1)`; on `\rho J` alternating, `D = diag(\rho, 1)` for a unit `\rho`. Absence: for
        Hermitian models `\det S = N(\det D) \det T = \det T` gives `\rho^n = 1`, `\rho` rational, and for `n \ge 3`
        the model form `H \oplus I_{n-2}` and its negative have different signatures. For alternating ones
        `\det(\rho J) = \det(D)^2` makes `\rho` a unit.
        """
        kind = source[0]
        S = self._star(matrix(self.field, source[4]).transpose()) * source_gram * matrix(self.field, source[4]).transpose()
        T = self._star(matrix(self.field, target[4]).transpose()) * target_gram * matrix(self.field, target[4]).transpose()
        if kind == "quaternion":
            return self._binary_hermitian_order(S).isometry_from(T) if target[0] == "quaternion" else None
        if kind == "ocean":
            embedding = hermitian_ocean.FIELD.embeddings(self.field)[0]
            found = hermitian_ocean.isometry(*_ocean_forms(self._pullback(S, embedding), self._pullback(T, embedding)))
            return None if found is None else found.apply_map(embedding)
        assert source[3] == target[3]
        i, j = next((i, j) for i in range(self.n) for j in range(self.n) if T[i, j] != 0)
        ratio = S[i, j] / T[i, j]
        assert S == ratio * T
        one = identity_matrix(self.field, self.n)
        if ratio == 1:
            return one
        if kind == "hyperbolic" and self.n == 2 and ratio == -1:
            return matrix(self.field, 2, 2, [1, 0, 0, -1])
        integers = self.field.ring_of_integers()
        if kind == "alternating" and ratio in integers and ~ratio in integers:
            return matrix(self.field, 2, 2, [ratio, 0, 0, 1])
        return None

    def equivariant_transporter(self, target_gram, targets):
        r"""Return an integral transporter to ``(target_gram, targets)``, or ``None`` as a proof of absence.

        With `t_0`, `h = B^{-1} t_0^t B' t_0` and `X = \{a : t_0 a \in M_r(Z)\}`, every transporter is `t_0 a_1`
        with `X = a_1 O`, so `\Lambda' = O_K X f = a_1 \Lambda` and `a_1` is an isometry of `(\Lambda, h)` onto
        `(\Lambda', h_h)`. Conversely, with one such isometry `b_0` (decided in the models), the isometries are
        `b_0 U(\Lambda, h)`, and `t_0 b_0 u` is integral exactly when `a_1^{-1} b_0 u \in U \cap O = C(Z)`. So the
        transporters are found on the coset representatives of `U(\Lambda, h) / C(Z)`.
        """
        commutant = self.commutant
        if abs(target_gram.det()) != abs(commutant.gram.det()):
            return None
        if self.kind == "first" and self.form_kind == "symmetric" and self.n == 2:
            return self._binary_orthogonal_transporter(target_gram, targets)
        start = rational_intertwiner(commutant.isometries, targets)
        if start is None:
            return None
        square = commutant.gram.change_ring(QQ).inverse() * start.transpose() * target_gram * start
        assert commutant.contains(square)
        target_lattice = self._lattice_rows(solution_lattice(commutant, start))
        twisted = self.gram_of(square)
        source = self.model(self._lattice, self.gram)
        target = self.model(target_lattice, twisted)
        block = self._model_isometry(source, target, self.gram, twisted)
        if block is None:
            return None
        source_change = matrix(self.field, source[4]).transpose()
        target_change = matrix(self.field, target[4]).transpose()
        isometry = self.element_of(target_change * block * source_change.inverse())
        one = identity_matrix(QQ, commutant.rank)
        assert commutant.adjoint(isometry) * square * isometry == one
        for representative in self.coset_representatives():
            found = _verified(commutant, start * isometry * representative, target_gram, targets)
            if found is not None:
                return found
        return None


def _diagonal_frame(gram):
    r"""Return `(Q, (d_1, d_2))` with `G = Q^t \mathrm{diag}(d_1, d_2) Q` for a nonsingular symmetric binary `G` over a field."""
    F = gram.base_ring()
    change = identity_matrix(F, 2)
    if gram[0, 0] == 0:
        change = matrix(F, 2, 2, [0, 1, 1, 0]) if gram[1, 1] != 0 else matrix(F, 2, 2, [1, 0, 1, 1])
    moved = change.transpose() * gram * change
    assert moved[0, 0] != 0
    shear = matrix(F, 2, 2, [1, moved[0, 1] / moved[0, 0], 0, 1])
    return shear * change.inverse(), (moved[0, 0], moved.det() / moved[0, 0])


def _binary_isometry(source, target):
    r"""Return `P \in GL_2(K)` with `P^t T P = S` for nonsingular symmetric binary ``source`` `S` and ``target`` `T`, or ``None``.

    Diagonalize `S \sim \langle \alpha_1, \alpha_2 \rangle` and `T \sim \langle \beta_1, \beta_2 \rangle`. An isometry sends the first
    basis vector to `x` with `T(x) = \alpha_1`, a point of the conic `\beta_1 X^2 + \beta_2 Y^2 = \alpha_1 Z^2` with `Z \ne 0`,
    or, when `T` is isotropic (a point with `Z = 0`), `x = u + s v` for the two null lines `u, v`. The `T`-orthogonal line of
    `x` is spanned by `y = (-\beta_2 x_2, \beta_1 x_1)` with `T(y) = \beta_1 \beta_2 \alpha_1`, so `\lambda y` has value `\alpha_2` exactly
    when `\alpha_2 / (\alpha_1 \beta_1 \beta_2) = \lambda^2` (Witt cancellation; Lam, *Introduction to Quadratic Forms over
    Fields*, Ch. I, Thm. 4.2). The conic is decided by Sage's number-field conic (PARI ``rnfisnorm``).
    """
    K = source.base_ring()
    source_change, (a1, a2) = _diagonal_frame(source)
    target_change, (b1, b2) = _diagonal_frame(target)
    ratio = a2 / (a1 * b1 * b2)
    if not ratio.is_square():
        return None
    found, point = Conic(K, [b1, b2, -a1]).has_rational_point(point=True)
    if not found:
        return None
    X, Y, Z = list(point)
    if Z == 0:
        u, v = vector(K, [X, Y]), vector(K, [X, -Y])
        pairing = b1 * X**2 - b2 * Y**2
        assert pairing != 0
        x = u + (a1 / (2 * pairing)) * v
    else:
        x = vector(K, [X / Z, Y / Z])
    y = ratio.sqrt() * vector(K, [-b2 * x[1], b1 * x[0]])
    result = target_change.inverse() * matrix(K, [x, y]).transpose() * source_change
    assert result.transpose() * target * result == source
    return result


def _ocean_forms(*forms):
    r"""Return `q \bar S` for the Hermitian ``forms`` `S` over `Q(i)`, with one `q > 0` making all integral.

    :mod:`hermitian_ocean` requires `a, c \in Z` and `2 \nu \in Z[i]` for `H = [[a, \nu], [\bar\nu, c]]`.
    """
    q = lcm([QQ(c).denominator() for S in forms for c in (S[0, 0], S[1, 1], *(2 * S[0, 1]).list())])
    return tuple(q * S.apply_map(lambda c: c.conjugate()) for S in forms)


def _missing_group(commutant):
    r"""Return the description of the arithmetic group of a commutant with no implemented route."""
    center = commutant._center_basis()
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
    return (
        f"the semisimple commutant of dimension {commutant.dimension()} of {commutant.isometries}, with a centre of "
        f"dimension {len(center)} that is not a field, needs the factor-wise unitary-group transporter (#54)"
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
    if _central_field(commutant):
        return CentralFieldCommutant(commutant).equivariant_transporter(target_gram, targets)
    commutative = all(x * y == y * x for x in commutant.order_basis for y in commutant.order_basis)
    if not radical and commutative and commutant.dimension() == 2 and len(targets) == 1:
        return commutant.equivariant_transporter(target_gram, targets[0])
    raise NotImplementedError(_missing_group(commutant))


def _cyclic_words(isometries, v):
    r"""Return words `w_1 = 1, w_2, \ldots` in the `F_i` with `w_j v` a basis of the `Q\langle F \rangle`-module of `v`."""
    one = identity_matrix(QQ, len(v))
    words, span = [one], matrix(QQ, [v])
    position = 0
    while position < len(words):
        current = words[position]
        position += 1
        for F in isometries:
            moved = F.change_ring(QQ) * current
            stacked = span.stack(matrix(QQ, [moved * v]))
            if stacked.rank() > span.rank():
                words.append(moved)
                span = stacked
    return words


class VectorOrbit:
    r"""The orbit of a vector `v \in L` under the centralizer `C(Z)` of isometries `F_i`, through the commutant.

    `g \in C(Z)` commutes with every word `w` in the `F_i`, so `g v = v'` exactly when `g X = X'` for the markings
    `X = [w_j v]`, `X' = [w_j v']` on the words of :func:`_cyclic_words` (the first word is `1`). With a primitive
    basis `S = X C^{-1}` of the saturated span, `g S = S' = X' C^{-1}`, and `F_i S = S E_i`. This is the marked
    transporter of :mod:`equivariant_markings`. When `S^t B S` is nonsingular, `L_Q = S_Q \oplus K_Q` with `K = S^\perp`,
    and `g = 1 \oplus h` on the frames, `h \in GL(K)` an equivariant isometry `(K, F|_K) \to (K', F|_{K'})`: `h = h_0 c`
    for one transporter `h_0` (:func:`equivariant_transporter`) and `c \in C_K(Z)` (:func:`centralizer_generators`).
    `g` preserves `L` exactly when `S' p_S + K' h p_K` is integral, which depends only on the point `c p_K` of the
    finite orbit of the glue `p_K` modulo `K`; Schreier's lemma on that orbit gives the stabilizer and the
    transporter decision (Holt-Eick-O'Brien, *Handbook of Computational Group Theory*, Thm. 4.13). A degenerate span
    is decided by :class:`EquivariantMarking`.
    """

    def __init__(self, gram, isometries, v):
        self.gram = matrix(ZZ, gram)
        self.isometries = tuple(matrix(ZZ, F) for F in isometries)
        self.vector = vector(ZZ, v)
        assert not self.vector.is_zero()
        self.words = _cyclic_words(self.isometries, vector(QQ, self.vector))
        columns = matrix(QQ, [w * self.vector for w in self.words]).transpose()
        self.marking = matrix(ZZ, columns.transpose().change_ring(ZZ).saturation().transpose())
        self._change = self.marking.change_ring(QQ).solve_right(columns)
        self.source_actions = tuple(
            matrix(ZZ, self.marking.change_ring(QQ).solve_right(F * self.marking)) for F in self.isometries
        )
        self.source_gram = self.marking.transpose() * self.gram * self.marking

    def _target_marking(self, w):
        columns = matrix(QQ, [word * vector(QQ, w) for word in self.words]).transpose()
        target = columns * self._change.inverse()
        if target.denominator() != 1:
            return None
        target = matrix(ZZ, target)
        if target.transpose() * self.gram * target != self.source_gram:
            return None
        if any(F * target != target * E for F, E in zip(self.isometries, self.source_actions, strict=True)):
            return None
        return target

    def _frame(self, marking):
        basis = (marking.transpose() * self.gram).right_kernel().basis_matrix().transpose()
        gram = basis.transpose() * self.gram * basis
        actions = tuple(matrix(ZZ, basis.change_ring(QQ).solve_right(F * basis)) for F in self.isometries)
        frame = block_matrix(QQ, [[marking, basis]], subdivide=False)
        return basis, gram, actions, frame

    def _glue(self):
        basis, gram, actions, frame = self._frame(self.marking)
        k = self.marking.ncols()
        coordinates = frame.inverse()
        generators = centralizer_generators(IsometryCommutant(gram, actions)) if basis.ncols() else ()
        return basis, gram, actions, coordinates[:k], coordinates[k:], generators

    def stabilizer_generators(self):
        r"""Return generators of `\{g \in C(Z) : g v = v\}`."""
        if self.source_gram.det() == 0:
            marking = EquivariantMarking(self.gram, self.marking, self.isometries, self.source_actions)
            return marking.pointwise_stabilizer_generators()
        basis, _gram, _actions, to_source, to_complement, generators = self._glue()
        if basis.ncols() == 0:
            return ()
        _transversal, stabilizer = _schreier(generators, to_complement, _glue_action, identity_matrix(ZZ, basis.ncols()))
        result = tuple(
            matrix(ZZ, self.marking * to_source + basis * h.change_ring(QQ) * to_complement) for h in stabilizer
        )
        assert all(g * self.vector == self.vector for g in result)
        return result

    def transporter(self, w):
        r"""Return `g \in C(Z)` with `g v = w`, or ``None`` as a proof that `w` is not in the orbit of `v`."""
        target = self._target_marking(w)
        if target is None:
            return None
        if self.source_gram.det() == 0:
            marking = EquivariantMarking(self.gram, self.marking, self.isometries, self.source_actions)
            return marking.transporter(target)
        basis, gram, actions, to_source, to_complement, generators = self._glue()
        target_basis, target_gram, target_actions, _frame = self._frame(target)
        found = None
        if basis.ncols() == 0:
            candidates = [target * to_source]
        else:
            start = equivariant_transporter(IsometryCommutant(gram, actions), target_gram, target_actions)
            if start is None:
                return None
            transversal, _stabilizer = _schreier(generators, to_complement, _glue_action, identity_matrix(ZZ, basis.ncols()))
            candidates = (
                target * to_source + target_basis * start.change_ring(QQ) * c.change_ring(QQ) * to_complement
                for c, _point in transversal.values()
            )
        found = next((matrix(ZZ, g) for g in candidates if g.denominator() == 1), None)
        if found is not None:
            assert IsometryCommutant(self.gram, self.isometries).is_centralizer_element(found)
            assert found * self.vector == vector(ZZ, w)
        return found
