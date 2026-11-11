r"""Cusp classes, genus neighbours, lattice isometries and cusp links for `\Lambda = U \perp M` over `O_K` (#54).

`K` is totally real of degree `d`, `M` is anisotropic over `K`, and `\Lambda = \mathfrak a_e e \oplus \mathfrak a_e^{-1} f \oplus \bigoplus \mathfrak a_i m_i`
is given on a pseudo-basis (:class:`witt_one_number_field.WittOneLattice`). Every lattice below is a `Z`-basis of an `O_K`-module in
`K^n` with an ambient Gram matrix. An isotropic line is `\Lambda \cap K v = \mathfrak b v` (`\mathfrak b` the intersection of the
`\mathfrak a_i v_i^{-1}`); its divisor is `\mathfrak b\, h(v, \Lambda)`, and the statements below for `v \in \Lambda` hold for `\mathfrak b v`
with the partner in `\mathfrak b \Lambda`, parities read in `\mathfrak b^2`, and the neighbour taken in `\mathfrak b^{-1} C` with values in
`\mathfrak b^{-2}`. The orbit of `\mathfrak a_e e` consists of the lines with `\mathfrak b \mathfrak a_e^{-1}` principal whose rescaled
vector of content `\mathfrak a_e` passes the test below; the heights are `N(y \mathfrak a_e^2)`.

Isometry test. For totally definite `M` of sign `s`, an `O_K`-isometry is a `Z`-isometry of the pairs of `Z`-lattices with the forms
`s \mathrm{Tr}(h), \mathrm{Tr}(\theta h), \ldots, \mathrm{Tr}(\theta^{d-1} h)`, `\theta` a generator of `K`: these forms determine `h`, and a
`Z`-map preserving `h` is `O_K`-linear by nondegeneracy (`h(g(\lambda x), g y) = \lambda h(x, y) = h(\lambda g x, g y)`). PARI ``qfisom``
(Plesken-Souvignier) decides it. Rank one: `\mathfrak a t \cong \mathfrak b t'` iff `\mathfrak b \mathfrak a^{-1} = (\gamma)` is principal and
`h(t, t) / (\gamma^2 h(t', t'))` is the square of a unit.

Cusps. For a primitive isotropic `v \in \Lambda` of divisor `h(v, \Lambda) = O_K` there is `v'` with `h(v, v') = 1` (an integral solution
of a `Z`-linear system), `\langle v, v' \rangle` is unimodular and `\Lambda = \langle v, v' \rangle \perp C_v` with `C_v \cong v^\perp / v`.
Changing `v'` by `t \in v^\perp` changes `h(v', v') \bmod 2` by `h(t, t)` (additive mod 2), so an even partner (`h(v', v') = 0` after
`v' \mapsto v' - \tfrac12 h(v', v') v`) exists iff a `GF(2)`-linear system is solvable; even partners with the same `v` give isometric
complements (they differ by an Eichler transvection fixing `v`). Hence `v = g e` for some `g \in O(\Lambda)` iff the divisor is `O_K`, an even
partner exists and `C_v \cong M`; then `g` is `e \mapsto v`, `f \mapsto v'`, `M \to C_v` by the isometry test (Witt extension).

Kneser neighbours (Kneser 1957; Kirschmer, *Definite quadratic and hermitian forms with small class number*, Ch. 2). For a principal prime
`\mathfrak p = \pi O_K \nmid 2 \det M` and `x \in C` with `h(x, x) \in 2 \pi^2 O_K`, `x \notin \pi C`, the vector
`u = \pi v + \lambda v' + x`, `\lambda = -h(x, x)/(2\pi) \in \mathfrak p`, is isotropic and `\psi(a v + b v' + c) = c - (a/\pi) x` is an
isometry `u^\perp / u \to C_x + O_K x/\pi` with `C_x = \{c : h(x, c) \in \mathfrak p\}` (`b = -(a\lambda + h(x, c))/\pi` is integral
iff `h(x, c) \in \mathfrak p`), the `\mathfrak p`-neighbour of `C` at `x`. The lift of `x \bmod \mathfrak p` with `h(x, x) \in \mathfrak p` to
`h(x, x) \in 2\pi^2` is `x_0 + \pi(A t_1 + B t_2)`, `A \in 2 O_K`, `A \equiv 1 \bmod \mathfrak p`, `B = 1 - A`: `t_1` solves the linear
condition mod `\mathfrak p`, `t_2` the `GF(2)`-linear condition mod 2. So every class reached is realised by an explicit cusp of `\Lambda`.
For rank `\ge 3` the iterated `\mathfrak p`-neighbours exhaust the spinor genus of `M` (strong approximation for the spin group at
`\mathfrak p`), which by Eichler's theorem is the set of `M'` with `U \perp M' \cong \Lambda`.

Cusp links (Swan, *Generators and relations for certain special linear groups*, 1971, in the norm form of the Euclidean covering of
:mod:`witt_one_number_field`). For `v = y(-q(z) e + f + z)` and a link `g_i e = v_i = y_i(-q(z_i) e + f + z_i)`, `q = h/2`,
`h(v_i, v) = -y y_i q(z - z_i)`, so `g_i^{-1} E(e, a)` sends `v` to a vector with `|N(y')| = |N(y)| |N(y_i)| |N(q(z + a - z_i))|`. Let
`G = \langle \mathrm{Stab}(e), g_i \rangle`. If every `z \in M \otimes R / M_{ev}` either has some `(i, a)` with
`|N(y_i)| |N(q(z + a - z_i))| < 1` or is the point of an isotropic line containing no vector of `O(\Lambda) e` (a singular cusp), then a
vector of `O(\Lambda) e` of least `|N(y)|` in its `G`-orbit has `y = 0`, i.e. lies in `O_K^\times e`, and `G = O(\Lambda)`. The links are the
`O(\Lambda) e`-vectors of bounded `|N(y)|`, found by the classification above, and `w` is the link `v_i = f`. Boxes of `[0,1]^{kd}` are
closed by interval arithmetic; a singular point `t_s` is closed by a punctured cube `0 < |t - t_s|_\infty \le r`: with `t = t_s + \varepsilon u`,
`|u|_\infty = 1`, and the links through `t_s` (`|N(y_i)| |N(q)| = 1` there, of sign `\sigma`), `\sigma |N(y_i)| N(q) - 1 = \varepsilon \psi_i(\varepsilon, u)`
for polynomials `\psi_i`, and the cube is covered once on `[0, r] \times \partial[-1, 1]^{kd}` some `\psi_i < 0` with `\sigma |N(y_i)| N(q) > -1`.
The link height grows until the certificate holds.
"""

from itertools import product

from sage.libs.pari import pari
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import RealIntervalField

from sage_indefinite_port.indefinite.witt_one_number_field import BOX_BUDGET, _pairing

_RIF = RealIntervalField(64)
_RADII = (QQ(1) / 4, QQ(1) / 16, QQ(1) / 64, QQ(1) / 256)


class _Coordinates:
    r"""`Z`-coordinates of `K` and `K^n` in the `Z`-basis of `O_K` returned by ``field.ideal(1).basis()``."""

    def __init__(self, field):
        self.field = field
        self.basis = [field(w) for w in field.ideal(1).basis()]
        self.d = len(self.basis)
        self.inverse = matrix(QQ, [w.list() for w in self.basis]).inverse()

    def scalar(self, x):
        return vector(QQ, self.field(x).list()) * self.inverse

    def flat(self, x):
        return vector(QQ, [c for entry in x for c in self.scalar(entry)])

    def unflat(self, row, n):
        return vector(self.field, [sum((row[i * self.d + j] * self.basis[j] for j in range(self.d)), self.field(0)) for i in range(n)])

    def mod2(self, x):
        return vector(GF(2), [ZZ(c) % 2 for c in self.scalar(x)])


def standard_z_basis(field, n):
    r"""Return the `Z`-basis `\{w \epsilon_i\}` of `O_K^n`."""
    coordinates = _Coordinates(field)
    return [vector(field, [w if i == j else 0 for i in range(n)]) for j in range(n) for w in coordinates.basis]


def _z_span(coordinates, vectors, rank, n):
    r"""Return a `Z`-basis of the `Z`-span of ``vectors`` (Hermite form of the coordinates)."""
    rows = matrix(QQ, [coordinates.flat(x) for x in vectors])
    denominator = rows.denominator()
    hermite = matrix(ZZ, denominator * rows).hermite_form(include_zero_rows=False)
    assert hermite.nrows() == rank
    return [coordinates.unflat(row / denominator, n) for row in hermite.rows()]


def _integer_solution(A, target):
    r"""Return `c \in Z^N` with `A c = \text{target}`, or ``None``; `A` integral (Smith form `D = U A V`)."""
    D, U, V = A.smith_form()
    right = U * vector(ZZ, target)
    y = [ZZ(0)] * A.ncols()
    for i in range(A.nrows()):
        pivot = D[i, i] if i < A.ncols() else 0
        if pivot == 0:
            if right[i] != 0:
                return None
            continue
        if right[i] % pivot != 0:
            return None
        y[i] = right[i] // pivot
    return V * vector(ZZ, y)


def _gf2_solution(rows, target):
    r"""Return a `0/1` integer vector `\beta` with `\sum \beta_i \text{rows}_i = \text{target}` over `GF(2)`, or ``None``."""
    A = matrix(GF(2), rows)
    if A.rank() != A.stack(matrix(GF(2), [target])).rank():
        return None
    return [ZZ(c) for c in A.solve_left(vector(GF(2), target))]


def _trace_forms(field, gram, basis, sign):
    theta = field.gen()
    forms = []
    for r in range(field.degree()):
        c = sign if r == 0 else theta**r
        forms.append(matrix(QQ, [[QQ(field(c * _pairing(gram, x, y)).trace()) for y in basis] for x in basis]))
    return forms


def _rank_one_module(basis):
    r"""Return `(\mathfrak a, t, j)` with `\mathfrak a t` the rank-one module spanned by ``basis`` and `t_j \ne 0`."""
    t = basis[0]
    j = next(i for i in range(len(t)) if t[i] != 0)
    field = t.base_ring()
    return field.ideal([x[j] / t[j] for x in basis]), t, j


def lattice_isometry(field, sign, gram_a, basis_a, gram_b, basis_b):
    r"""Return the images in `K^n` of ``basis_a`` under an `O_K`-isometry onto the span of ``basis_b``, or ``None`` (module docstring).

    ``sign`` is `\pm 1` when the lattices are totally definite of that sign and `0` otherwise."""
    if len(basis_a) != len(basis_b):
        return None
    d = field.degree()
    rank = len(basis_a) // d
    if rank == 1:
        # phi(t_a) = s t_b is an isometry of a t_a onto b t_b iff s a = b and s^2 h(t_b) = h(t_a): b a^{-1} = (gamma), s = gamma u with
        # u^2 = h(t_a) / (gamma^2 h(t_b)) a unit square.
        (ia, ta, ja), (ib, tb, _jb) = _rank_one_module(basis_a), _rank_one_module(basis_b)
        quotient = ib / ia
        if not quotient.is_principal():
            return None
        gamma = quotient.gens_reduced()[0]
        ratio = field(_pairing(gram_a, ta, ta) / (gamma**2 * _pairing(gram_b, tb, tb)))
        if not (ratio.is_integral() and ratio.is_unit() and ratio.is_square()):
            return None
        s = gamma * ratio.sqrt()
        return [(x[ja] / ta[ja]) * s * tb for x in basis_a]
    if sign == 0:
        raise NotImplementedError(
            f"isometry test of rank-{rank} O_K-lattices over {field} that are not totally definite: the commutant transporter route of "
            f"commutant_groups.equivariant_transporter for the anisotropic indefinite case is needed (#54)"
        )
    forms_a = _trace_forms(field, gram_a, basis_a, sign)
    forms_b = _trace_forms(field, gram_b, basis_b, sign)
    scale = ZZ(1)
    for F in (*forms_a, *forms_b):
        scale = scale.lcm(F.denominator())
    found = pari([matrix(ZZ, scale * F) for F in forms_a]).qfisom([matrix(ZZ, scale * F) for F in forms_b])
    if found.type() == "t_INT":
        return None
    S = matrix(ZZ, found.sage())
    n = len(basis_b[0])
    images = [sum((S[i, j] * basis_b[i] for i in range(len(basis_b))), vector(field, n)) for j in range(len(basis_a))]
    assert all(
        _pairing(gram_b, images[i], images[j]) == _pairing(gram_a, basis_a[i], basis_a[j]) for i in range(len(images)) for j in range(len(images))
    )
    return images


class Cusps:
    r"""Frames, classification, neighbours and links of isotropic lines of `U \perp M` on a pseudo-basis (module docstring)."""

    def __init__(self, lattice):
        self.lattice = lattice
        self.field = lattice.field
        self.n = lattice.k + 2
        self.coordinates = _Coordinates(self.field)
        self.d = self.coordinates.d
        self.ae = lattice.ideals[0]
        self.e, self.f = lattice.e, lattice.f
        self.lam = self.scaled(self.field.ideal(1))
        signs = set(lattice.signs)
        self.sign = lattice.signs[0] if len(signs) == 1 and 0 not in signs else 0

    def h(self, x, y):
        return self.lattice.h(x, y)

    def scaled(self, ideal):
        r"""Return a `Z`-basis of `\mathfrak c \Lambda = \mathfrak c \mathfrak a_e e \oplus \mathfrak c \mathfrak a_e^{-1} f \oplus \mathfrak c M`."""
        return (
            [w * self.e for w in (ideal * self.ae).basis()]
            + [w * self.f for w in (ideal / self.ae).basis()]
            + self.lattice._scaled_z_basis(ideal)
        )

    def content(self, v):
        r"""Return `\mathfrak b` with `\Lambda \cap K v = \mathfrak b v`: the intersection of the `\mathfrak a_i v_i^{-1}`, `v_i \ne 0`."""
        result = None
        for a, c in zip(self.lattice.ideals, v):
            if c != 0:
                term = a * self.field.ideal(c) ** -1
                result = term if result is None else result.intersection(term)
        return result

    def divisor(self, v):
        r"""Return the integral ideal `\mathfrak b\, h(v, \Lambda)` of the line `\mathfrak b v = \Lambda \cap K v`."""
        return self.field.ideal([self.h(v, b) for b in self.lam]) * self.content(v)

    def _parity(self, ideal, x):
        r"""Return the coordinates of `x \in` ``ideal`` modulo 2 in the `Z`-basis of ``ideal``."""
        basis = matrix(QQ, [self.field(w).list() for w in ideal.basis()])
        coordinates = basis.solve_left(vector(QQ, self.field(x).list()))
        assert all(c.is_integral() for c in coordinates)
        return vector(GF(2), [ZZ(c) % 2 for c in coordinates])

    def _in(self, x, ideal):
        return x == 0 or (self.field.ideal(x) / ideal).is_integral()

    def frame(self, v):
        r"""Return ``(content, partner, complement, even)`` for isotropic ``v`` of divisor `O_K`, or ``None`` otherwise.

        With `\Lambda \cap K v = \mathfrak b v`, the partner `v' \in \mathfrak b \Lambda` has `h(v, v') = 1`, so `\mathfrak b v \oplus \mathfrak b^{-1} v'`
        is a unimodular hyperbolic summand; it is even when `h(v', v') \in 2 \mathfrak b^2` can be reached by `v' \mapsto v' + t`, `t \in v^\perp`."""
        assert self.h(v, v) == 0
        b = self.content(v)
        basis = self.scaled(b)
        A = matrix(ZZ, [[ZZ(c) for c in self.coordinates.scalar(self.h(v, x))] for x in basis]).transpose()
        c = _integer_solution(A, [ZZ(x) for x in self.coordinates.scalar(1)])
        if c is None:
            return None
        partner = sum((ci * x for ci, x in zip(c, basis)), vector(self.field, self.n))
        kernel = A.right_kernel_matrix()
        perp = [sum((ZZ(r) * x for r, x in zip(row, basis)), vector(self.field, self.n)) for row in kernel.rows()]
        square = b**2
        beta = _gf2_solution([self._parity(square, self.h(t, t)) for t in perp], self._parity(square, self.h(partner, partner)))
        even = beta is not None
        if even:
            partner = partner + sum((ZZ(r) * t for r, t in zip(beta, perp)), vector(self.field, self.n))
            partner = partner - (self.h(partner, partner) / 2) * v
            assert self.h(partner, partner) == 0
        assert self.h(v, partner) == 1
        c0 = self.h(partner, partner)
        projected = [x - (self.h(x, partner) - c0 * self.h(x, v)) * v - self.h(x, v) * partner for x in self.lam]
        complement = _z_span(self.coordinates, projected, self.lattice.k * self.d, self.n)
        return b, partner, complement, even

    def link(self, v):
        r"""Return `g \in O(\Lambda)` with `g e \in K^\times v`, or ``None`` when the line of `v` is not in `O(\Lambda) \cdot \mathfrak a_e e`.

        The line `\mathfrak b v` is in the orbit iff `\mathfrak b \mathfrak a_e^{-1} = (c)` is principal and `c v` (content `\mathfrak a_e`) has an
        even partner with complement isometric to `M` (module docstring)."""
        quotient = self.content(v) / self.ae
        if not quotient.is_principal():
            return None
        v = quotient.gens_reduced()[0] * v
        framed = self.frame(v)
        if framed is None or not framed[3]:
            return None
        _b, partner, complement, _even = framed
        images = lattice_isometry(self.field, self.sign, self.lattice.gram, self.lattice.z_basis, self.lattice.gram, complement)
        if images is None:
            return None
        g = matrix(self.field, self.n, self.n)
        g.set_column(0, v)
        g.set_column(1, partner)
        for j in range(self.lattice.k):
            w = self.lattice.z_basis[j * self.d][2 + j]
            g.set_column(2 + j, images[j * self.d] / w)
        assert g.transpose() * self.lattice.gram * g == self.lattice.gram and self.lattice.preserves(g)
        return g

    # Kneser neighbours.

    def good_prime(self):
        r"""Return `(\mathfrak p, \pi)`: a principal prime of least norm prime to `2 \det M \prod \mathfrak a_i^2`."""
        bad = self.field.ideal(2 * self.lattice.inner.det())
        for a in self.lattice.ideals[2:]:
            bad = bad * a**2
        bound = 2
        while True:
            for _norm, ideals in sorted(self.field.ideals_of_bdd_norm(bound).items()):
                for ideal in ideals:
                    if ideal.is_prime() and ideal.is_principal() and bad.valuation(ideal) == 0:
                        return ideal, ideal.gens_reduced()[0]
            bound *= 2

    def _residues(self, basis, sub):
        r"""Return representatives of `\text{span}(basis) / \text{span}(sub)` (Smith form)."""
        flat = matrix(QQ, [self.coordinates.flat(b) for b in basis])
        A = matrix(ZZ, [flat.solve_left(self.coordinates.flat(s)) for s in sub])
        D, U, V = A.smith_form()
        Vinv = V.inverse()
        rows = [Vinv.row(i) for i in range(len(basis))]
        for r in product(*[range(D[i, i]) for i in range(len(basis))]):
            yield sum((ri * sum((c * b for c, b in zip(row, basis)), vector(self.field, self.n)) for ri, row in zip(r, rows)), vector(self.field, self.n))

    def scaled_complement(self, b, complement):
        r"""Return a `Z`-basis of `\mathfrak b^{-1} C`."""
        return _z_span(self.coordinates, [w * c for w in (b**-1).basis() for c in complement], len(complement), self.n)

    def neighbour_vector(self, v, b, partner, tilde, x0, prime, pi):
        r"""Return the cusp `u = \pi v + \lambda v' + x` of the `\mathfrak p`-neighbour at the line of ``x0`` `\in \tilde C = \mathfrak b^{-1} C`
        (``tilde``), or ``None`` when the mod-2 lift does not exist (module docstring). ``partner`` must be even.

        `\tilde C` carries values in `J = \mathfrak b^{-2}`; the conditions `h(x, x) \in 2\pi^2 J` and `\lambda = -h(x, x)/(2\pi) \in \pi J` make
        `u \in \mathfrak b^{-1} \Lambda`."""
        h = self.h
        J = b**-2
        assert h(partner, partner) == 0 and self._in(h(x0, x0), prime * J)
        pivot = next(c for c in tilde if not self._in(h(x0, c), prime * J))
        residue = prime.residue_field()
        s = residue.lift(residue(-h(x0, x0) / (pi * 2 * h(x0, pivot))))
        A = 2 * residue.lift(1 / residue(2))
        B = 1 - A
        x1 = x0 + pi * A * s * pivot
        square = pi**2
        rho = _gf2_solution([self.coordinates.mod2(square * w) for w in self.coordinates.basis], self.coordinates.mod2(1))
        assert rho is not None
        rho = sum((r * w for r, w in zip(rho, self.coordinates.basis)), self.field(0))
        beta = _gf2_solution([self._parity(J, h(c, c)) for c in tilde], self._parity(J, rho * h(x1, x1)))
        if beta is None:
            return None
        x = x1 + pi * B * sum((ZZ(r) * c for r, c in zip(beta, tilde)), vector(self.field, self.n))
        assert self._in(h(x, x) / (2 * square), J)
        u = pi * v + (-h(x, x) / (2 * pi)) * partner + x
        assert h(u, u) == 0
        return u

    def cusp_classes(self, prime=None):
        r"""Return the cusp classes reached by iterated `\mathfrak p`-neighbours from `\mathfrak a_e e`, as dicts with the cusp ``vector``, its
        ``content``, ``partner``, ``complement``, ``even`` and, when the line lies in `O(\Lambda) \mathfrak a_e e`, its ``link``. Classes are
        distinct by the invariant (parity of the partner, isometry class of `v^\perp/v`); the first is `e`."""
        if prime is None:
            prime, pi = self.good_prime()
        else:
            pi = prime.gens_reduced()[0]
        classes = [
            dict(vector=self.e, content=self.ae, partner=self.f, complement=list(self.lattice.z_basis), even=True, link=identity_matrix(self.field, self.n))
        ]
        queue = [0]
        while queue:
            current = classes[queue.pop()]
            if not current["even"]:
                continue
            tilde = self.scaled_complement(current["content"], current["complement"])
            for x0 in self._residues(tilde, [pi * c for c in tilde]):
                if x0.is_zero() or not self._in(self.h(x0, x0), prime * current["content"] ** -2):
                    continue
                u = self.neighbour_vector(current["vector"], current["content"], current["partner"], tilde, x0, prime, pi)
                if u is None:
                    continue
                framed = self.frame(u)
                if framed is None:
                    continue
                content, partner, new, even = framed
                if any(
                    c["even"] == even and lattice_isometry(self.field, self.sign, self.lattice.gram, c["complement"], self.lattice.gram, new) is not None
                    for c in classes
                ):
                    continue
                classes.append(dict(vector=u, content=content, partner=partner, complement=new, even=even, link=self.link(u)))
                queue.append(len(classes) - 1)
        return classes

    # Links and the generalised Euclidean certificate.

    def _t_coordinates(self, z):
        even = matrix(QQ, [self.coordinates.flat(b[2:]) for b in self.lattice.even_e])
        return even.solve_left(self.coordinates.flat(z))

    def _cusps_of_height(self, height):
        r"""Yield ``(kind, t, g, norm)`` for the isotropic lines `\mathfrak b v`, `v = x e + y f + m \in \mathfrak a_e^{-1} \Lambda`, with
        `N(y \mathfrak a_e^2) = ` ``height``, modulo `y P(\mathfrak a_e)` and units: ``kind`` is ``"link"`` (``g e`` spans the line,
        ``norm`` `= |N(h(e, g e))|`) or ``"singular"`` (the line is not in `O(\Lambda) \mathfrak a_e e`)."""
        lattice = self.lattice
        module = lattice._scaled_z_basis(self.ae**-1)
        for ideal in self.field.ideals_of_bdd_norm(height).get(height, []):
            scaled = ideal / self.ae**2
            if not scaled.is_principal():
                continue
            y = scaled.gens_reduced()[0]
            for m in self._residues(module, [y * b for b in lattice.even_e]):
                x = -lattice.h(m, m) / (2 * y)
                if not x.is_integral():
                    continue
                v = x * self.e + y * self.f + m
                t = self._t_coordinates(m[2:] / y)
                content = self.content(v)
                if content != self.ae:
                    if (content / self.ae).is_principal():
                        continue
                    yield "singular", t, None, None
                    continue
                g = self.link(v)
                yield ("link", t, g, abs(g[1, 0].norm())) if g is not None else ("singular", t, None, None)

    def _norm_polynomial(self):
        r"""Return `Q(\tau) = N(q(\sum \tau_j b_j))` over `Q`, `b_j` the `Z`-basis of `M_{ev}`."""
        basis = self.lattice.even_e
        n = len(basis)
        R = PolynomialRing(QQ, "tau", n)
        tau = R.gens()
        parts = [R(0)] * self.d
        for j in range(n):
            for l in range(n):
                for r, a in enumerate(self.field(self.lattice.h(basis[j], basis[l]) / 2).list()):
                    parts[r] += a * tau[j] * tau[l]
        theta = self.field.gen()
        total = matrix(R, self.d, self.d)
        for r in range(self.d):
            total += parts[r] * (theta**r).matrix().change_ring(R)
        return total.det()

    def _local_radius(self, polynomial, t_s, links):
        r"""Return `r` such that the punctured cube of radius `r` about ``t_s`` is covered (module docstring), or ``None``."""
        n = len(t_s)
        S = PolynomialRing(QQ, ["eps", *(f"u{j}" for j in range(n))])
        eps, us = S.gens()[0], S.gens()[1:]
        pairs = []
        for _kind, t_i, _g, norm in links:
            centre = [ZZ(round(t_i[j] - t_s[j])) for j in range(n)]
            for shift in product((-1, 0, 1), repeat=n):
                tau0 = [t_s[j] + centre[j] + shift[j] - t_i[j] for j in range(n)]
                value = norm * polynomial(*tau0)
                if abs(value) != 1:
                    continue
                sigma = 1 if value > 0 else -1
                full = sigma * norm * polynomial(*[tau0[j] + eps * us[j] for j in range(n)]) - 1
                psi, remainder = full.quo_rem(eps)
                assert remainder == 0
                pairs.append((psi, full))
        if not pairs:
            return None
        for radius in _RADII:
            if all(self._face_covered(pairs, n, j, sign, radius) for j in range(n) for sign in (-1, 1)):
                return radius
        return None

    def _face_covered(self, pairs, n, j, sign, radius):
        stack = [((QQ(0), radius), *(((QQ(-1), QQ(1)) for _ in range(n - 1))))]
        visited = 0
        while stack:
            box = stack.pop()
            visited += 1
            if visited > BOX_BUDGET:
                return False
            free = list(box[1:])
            u = [_RIF(sign) if l == j else None for l in range(n)]
            for l in range(n):
                if l != j:
                    a, b = free.pop(0)
                    u[l] = _RIF(a, b)
            args = [_RIF(box[0][0], box[0][1]), *u]
            if any(psi(*args).upper() < 0 and full(*args).lower() > -2 for psi, full in pairs):
                continue
            k = max(range(len(box)), key=lambda i: box[i][1] - box[i][0])
            middle = (box[k][0] + box[k][1]) / 2
            stack.append(tuple((box[i][0], middle) if i == k else box[i] for i in range(len(box))))
            stack.append(tuple((middle, box[i][1]) if i == k else box[i] for i in range(len(box))))
        return True

    def _certified(self, grams, links, cubes):
        field = grams[0].base_ring()
        n = grams[0].nrows()
        limit = field(2) ** self.d
        stack = [tuple((QQ(0), QQ(1)) for _ in range(n))]
        visited = 0
        while stack:
            box = stack.pop()
            visited += 1
            if visited > BOX_BUDGET:
                return False
            centre = [(a + b) / 2 for a, b in box]
            if any(
                all(
                    t_s[j] + ZZ(round(centre[j] - t_s[j])) - r <= box[j][0] and box[j][1] <= t_s[j] + ZZ(round(centre[j] - t_s[j])) + r
                    for j in range(n)
                )
                for t_s, r in cubes
            ):
                continue
            if self._box_covered(box, centre, grams, field, limit, links):
                continue
            k = max(range(n), key=lambda i: box[i][1] - box[i][0])
            middle = (box[k][0] + box[k][1]) / 2
            stack.append(tuple((box[i][0], middle) if i == k else box[i] for i in range(n)))
            stack.append(tuple((middle, box[i][1]) if i == k else box[i] for i in range(n)))
        return True

    def _box_covered(self, box, centre, grams, field, limit, links):
        n = len(box)
        for _kind, t_i, _g, norm in links:
            bound = limit / norm
            base = [ZZ(round(t_i[j] - centre[j])) for j in range(n)]
            for shift in product((-1, 0, 1), repeat=n):
                offsets = [base[j] + shift[j] - t_i[j] for j in range(n)]
                tau = [field(box[j][0] + offsets[j], box[j][1] + offsets[j]) for j in range(n)]
                value = field(1)
                for gram in grams:
                    value = value * sum(gram[i, l] * tau[i] * tau[l] for i in range(n) for l in range(n))
                if abs(value).upper() < bound.lower():
                    return True
        return False

    def links(self):
        r"""Return links `g_i` with `\langle \mathrm{Stab}(e), g_i \rangle = O(\Lambda)`, raising the link height until the certificate holds."""
        grams = self.lattice._quadratic_intervals(self.lattice.even_e)
        polynomial = self._norm_polynomial()
        links, singular, cubes = [], [], []
        height = 0
        while True:
            height += 1
            for kind, t, g, norm in self._cusps_of_height(height):
                reduced = vector(QQ, [c - c.floor() for c in t])
                (links if kind == "link" else singular).append((kind, reduced, g, norm))
            cubes = []
            for _kind, t_s, _g, _norm in singular:
                radius = self._local_radius(polynomial, t_s, links)
                if radius is not None:
                    cubes.append((t_s, radius))
            if links and self._certified(grams, links, cubes):
                return tuple(g for _kind, _t, g, _norm in links)


def cusp_links(lattice):
    r"""Return the cusp links completing `\mathrm{Stab}(e)` to generators of `O(\Lambda)` when the Euclidean certificate fails."""
    return Cusps(lattice).links()
