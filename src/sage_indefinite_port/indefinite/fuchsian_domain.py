r"""Fundamental domains and generators of arithmetic Fuchsian groups from quaternion orders (#54).

Algorithm. Voight, *Computing fundamental domains for Fuchsian groups*, J. Théor. Nombres Bordeaux 21 (2009), 469-491 (cited as [V]),
and Voight, *Quaternion Algebras*, Ch. 38-39 (cited as [QA]). Let `F` be a totally real field of degree `n`, `B = (a, b)_F` a division
quaternion algebra split at exactly one real place `v` and ramified at the other real places `w`, and `O` an order. Then
`\Gamma = \iota_v(O^1) / \{\pm 1\} \subset PSL_2(R)` is a cocompact Fuchsian group ([QA] Sec. 38.4).

Splitting and the disc ([V] Sec. 1-2). At `v` one of `a, b` is positive; for `\sigma_v(a) > 0`,
`i \mapsto \mathrm{diag}(\sqrt a, -\sqrt a)`, `j \mapsto [[0, b], [1, 0]]` is an embedding `\iota_v : B \to M_2(R)` (symmetrically for `b`). A point `p`
of the upper half-plane `H` is moved to `0` of the unit disc `D` by `\phi(z) = (z - p)/(z - \bar p)`, and `\gamma \mapsto \gamma_D = M \gamma M^{-1}`,
`M = [[1, -p], [1, -\bar p]]`, carries `SL_2(R)` to `SU(1, 1)`.

Isometric circles and the exterior domain ([V] Sec. 1, Prop. 1.3-1.5). For `\gamma_D = [[a, b], [c, d]]` with `c \ne 0`, the isometric circle
`I(\gamma) = \{z : |c z + d| = 1\}` has centre `-d/c` and radius `1/|c|`, and it is the perpendicular bisector of `0` and `\gamma^{-1}(0)`. For a
finite set `G` the exterior domain `E(G) = \bigcap_{g \in G} \{|c_g z + d_g| > 1\}` is an intersection of hyperbolic half-planes containing
`0`. In the Klein model `w = 2z/(1 + |z|^2)` the circle with centre `C` becomes the chord `\mathrm{Re}(w \bar C) = 1`, so `E(G)` is a convex
Euclidean polygon cut from the disc; it is computed by clipping half-planes in the Klein model. Its hyperbolic area is
`(k - 2)\pi - \sum \theta_i` (Gauss-Bonnet), the angles `\theta_i` measured in the conformal Poincaré model. When `p` has trivial stabilizer
in `\Gamma` the Dirichlet domain of `\Gamma` centred at `p` is `E(\Gamma \setminus \{\pm 1\})`, a fundamental domain ([V] Prop. 1.5).

Enumeration ([V] Sec. 4, the absolute reduced norm). `N(x) = \|\iota_v(x)_D\|^2 + \sum_w \sigma_w(\mathrm{nrd}\, x)` (Frobenius norm of the
disc matrix plus the reduced norms at the ramified real places, positive definite there) is a positive definite quadratic form on
`O \otimes_Z R`. For `\gamma \in O^1`, `\|\gamma_D\|^2 = 2 (|c|^2 + |d|^2) = 2 \cosh d(0, \gamma 0)` and the ramified terms sum to `n - 1`, so the
Fincke-Pohst ellipsoids of `N` enumerate `O^1` by the hyperbolic distance of `\gamma(0)` from `0`. The real Gram matrix is rounded to a
rational one and the bound enlarged by the rounding error divided by the least eigenvalue, so no vector is lost; every candidate is
then tested exactly in `B`.

Reduction ([V] Alg. 4.3). For `z \in D`, while `z` lies inside some `I(g)`, `g \in G`, replace `z` by `g z` for the `g` minimizing `|c_g z + d_g|`;
each step decreases `d(z, 0)`. Reducing `\gamma(0)` gives `\delta \in \langle G \rangle` with `\delta \gamma (0) \in \overline{E(G)}`; if
`\bar\gamma = \delta \gamma \ne \pm 1`, its isometric circles cut `E(G)` and `\bar\gamma^{\pm 1}` are added. An enumerated `\gamma` whose circle cuts
`E(G)` is also added (with its inverse), so after the enumeration passes the radius of the Dirichlet domain, `E(G)` equals it.

Stop condition ([V] Sec. 5, (5.1); [QA] Main Thm. 39.1.13). For a maximal order the area of `\Gamma \backslash H` is Shimizu's
`\mathrm{vol} = 8\pi (4\pi^2)^{-n} d_F^{3/2} \zeta_F(2) \Phi(\mathfrak D)`, `\Phi(\mathfrak D) = \prod_{\mathfrak p \mid \mathfrak D} (N\mathfrak p - 1)`.
`\zeta_F(2)` is computed by PARI's ``lfun``. `\mathrm{vol}/2\pi = 2g - 2 + \sum (1 - 1/e_i)` is rational with denominator dividing the lcm
of the possible elliptic orders `e` (those with `\varphi(2e) \le 2n`, as `\mathbb Q(\cos \pi/e) \subseteq F`); both `\mathrm{vol}/2\pi` and the
area of `E(G)` over `2\pi` are rounded to that denominator and compared exactly. Always `\mathrm{area}\, E(G) \ge [\Gamma : \langle G \rangle]
\mathrm{vol}`, since `E(G)` contains the Dirichlet domain of `\langle G \rangle`; so equality proves `\langle G \rangle = \Gamma` and `E(G)` is the
Dirichlet domain.

Side pairing ([V] Prop. 1.6, Poincaré's theorem). Each side lies on `I(g)` and `g` carries it onto the side on `I(g^{-1})`; these side
pairings generate `\Gamma`, and with `-1` they generate `O^1`. The vertex cycles give the elliptic generators: the cycle transformation of
angle sum `\theta` has order `2\pi/\theta`; a side paired with itself contributes an elliptic element of order 2. The genus follows from
the area.

Units and the normalizer. `\mathrm{nrd}(O^\times) = O_{F, \Omega}^\times`, the units positive at the ramified real places, for a maximal (more
generally Eichler) order, and a two-sided ideal of a maximal order is principal exactly when its reduced norm has a generator positive at
those places ([QA] Sec. 28.5, Eichler's norm theorem and strong approximation; the Eichler condition holds as `v` splits). So
`O^\times` is generated by `O^1`, the central units, and elements of each norm in a basis of `O_{F, \Omega}^\times / (O_F^\times)^2`, found by
the same enumeration, which ends because they exist; generators of two-sided ideals are found by the enumeration in the ideal. A
non-maximal `O` is treated through a maximal `O' \supseteq O`: `O^1` and `O^\times` are the elements of `O'^1` and `O'^\times` lying in `O`, a
subgroup of finite index (`O \supseteq m O'`) with Schreier generators.
"""

from sage.arith.functions import lcm
from sage.arith.misc import euler_phi
from sage.libs.pari import pari
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.complex_mpfr import ComplexField
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RealField
from sage_extensions.categories.lattice_centralizers import finite_preimage_generators

from sage_indefinite_port.indefinite.quaternion_orders_number_field import OrderLattices, maximal_orders_containing
from sage_indefinite_port.indefinite.quaternion_voronoi_number_field import QuaternionOrderVoronoi


def _real_places(field, real_field):
    r"""Return the real embeddings of the totally real ``field`` (`\mathbb Q` included) as functions into ``real_field``."""
    if field is QQ:
        return [lambda a: real_field(QQ(a))]
    return [lambda a, s=s: s(field(a)) for s in field.embeddings(real_field)]


def fuchsian_places(algebra, real_field):
    r"""Return ``(split, ramified)`` real places when ``algebra`` is split at exactly one real place of its totally real centre, else ``None``.

    `B = (a, b)_F` is ramified at a real `\sigma` exactly when `\sigma(a) < 0` and `\sigma(b) < 0` (the Hilbert symbol at `\sigma`)."""
    field = algebra.base_ring()
    if field is not QQ and not field.is_totally_real():
        return None
    a, b = algebra.invariants()
    places = _real_places(field, real_field)
    ramified = [s for s in places if s(a) < 0 and s(b) < 0]
    split = [s for s in places if not (s(a) < 0 and s(b) < 0)]
    if len(split) != 1:
        return None
    return split[0], ramified


def ramified_real_places(algebra, real_field):
    r"""Return the real places of the totally real centre at which ``algebra`` ramifies, as functions into ``real_field``."""
    a, b = algebra.invariants()
    return [s for s in _real_places(algebra.base_ring(), real_field) if s(a) < 0 and s(b) < 0]


def _central_units(field):
    if field is QQ:
        return [QQ(-1)]
    return [field(u) for u in field.unit_group().gens_values()]


def _sign_matrix(units, places):
    return matrix(GF(2), len(units), len(places), [1 if s(u) < 0 else 0 for u in units for s in places])


def _product(elements, exponents, one):
    result = one
    for x, e in zip(elements, exponents):
        result *= x ** int(e)
    return result


def positive_unit_classes(field, places):
    r"""Return units whose classes form an `F_2`-basis of the units positive at ``places`` modulo squares."""
    units = _central_units(field)
    signs = _sign_matrix(units, places)
    return [_product(units, k, field(1)) for k in signs.left_kernel().basis()]


def positive_generator(field, ideal, places):
    r"""Return a generator of the principal ``ideal`` positive at the real ``places``, or ``None`` when there is none."""
    if not ideal.is_principal():
        return None
    g = field(ideal.gens_reduced()[0])
    if not places:
        return g
    units = _central_units(field)
    target = vector(GF(2), [1 if s(g) < 0 else 0 for s in places])
    try:
        exponents = _sign_matrix(units, places).solve_left(target)
    except ValueError:
        return None
    result = g * _product(units, exponents, field(1))
    assert all(s(result) > 0 for s in places)
    return result


def quaternion_unit_engine(algebra, order_basis):
    r"""Return the unit-group engine of an order: :class:`FuchsianDomain` for an algebra with real ramification split at exactly one
    real place, :class:`QuaternionOrderVoronoi` otherwise (totally definite or totally indefinite, Albert's classification)."""
    places = fuchsian_places(algebra, RealField(64))
    if places is not None and places[1]:
        return FuchsianDomain(algebra, order_basis)
    return QuaternionOrderVoronoi(algebra, order_basis)


class FuchsianFundamentalDomain:
    r"""The Dirichlet domain of `\iota_v(O^1)/\{\pm 1\}` centred at `p`: vertices in the disc, side elements, pairing, elliptic cycles."""

    def __init__(self, vertices, sides, pairing, elliptic, area, genus):
        self.vertices = tuple(vertices)
        self.sides = tuple(sides)
        self.pairing = tuple(pairing)
        self.elliptic = tuple(elliptic)
        self.area_over_2pi = area
        self.genus = genus

    def signature(self):
        r"""Return `(g; e_1, \ldots, e_r)` with the elliptic orders sorted."""
        return self.genus, tuple(sorted(e for _g, e in self.elliptic))


class FuchsianDomain:
    r"""Voight's fundamental-domain algorithm for the order ``order_basis`` of ``algebra`` (module docstring)."""

    def __init__(self, algebra, order_basis, precision=256):
        self.algebra = algebra
        self.center = algebra.base_ring()
        self.degree = 1 if self.center is QQ else self.center.degree()
        self.precision = precision
        self.R = RealField(precision)
        self.C = ComplexField(precision)
        self.epsilon = self.R(2) ** (-(precision // 2))
        self.tolerance = self.R(2) ** (-(precision // 4))
        places = fuchsian_places(algebra, self.R)
        assert places is not None, f"{algebra} is not split at exactly one real place of a totally real centre"
        self.split, self.ramified = places
        assert self._finite_ramified() or self.ramified, f"{algebra} is M_2(F): its unit groups are not cocompact"
        self.lattices = OrderLattices(algebra)
        self.order = tuple(algebra(x) for x in order_basis)
        assert len(self.order) == 4 * self.degree
        self._order_rows = matrix(QQ, [self.lattices.coordinates(x) for x in self.order])
        assert self._order_rows.det() != 0
        self._frame = self._split_frame()
        self._set_center(0)

    # Exact arithmetic.
    def _finite_ramified(self):
        if self.center is QQ:
            return [ZZ(p) for p in self.algebra.ramified_primes()]
        return list(self.algebra.ramified_places(inf=False))

    def _prime_norm(self, prime):
        return abs(ZZ(prime)) if self.center is QQ else prime.norm()

    def order_coordinates(self, x):
        return self._order_rows.solve_left(vector(QQ, self.lattices.coordinates(x)))

    def in_order(self, x):
        return all(c in ZZ for c in self.order_coordinates(x))

    def is_unit(self, x):
        return self.in_order(x) and self.in_order(~x)

    def trace(self, x):
        t = self.algebra(x).reduced_trace()
        return QQ(t) if self.center is QQ else QQ(self.center(t).trace())

    def coordinates(self, x):
        return self.lattices.coordinates(x)

    def element(self, row):
        return self.lattices.element(row)

    def expected_discriminant(self):
        reduced = ZZ(1)
        for prime in self._finite_ramified():
            reduced *= self._prime_norm(prime)
        field_discriminant = 1 if self.center is QQ else abs(self.center.discriminant())
        return field_discriminant**4 * reduced**2

    def is_maximal(self):
        gram = matrix(QQ, [[self.trace(x * y) for y in self.order] for x in self.order])
        return abs(gram.det()) == self.expected_discriminant()

    # The split embedding and the disc model.
    def _split_frame(self):
        r"""Return `\iota_v(1), \iota_v(i), \iota_v(j), \iota_v(k)` as complex `2 \times 2` tuples ``(a, b, c, d)`` (module docstring)."""
        C = self.C
        a, b = self.algebra.invariants()
        A, B = self.split(a), self.split(b)
        if A > 0:
            r = A.sqrt()
            i, j = (C(r), C(0), C(0), C(-r)), (C(0), C(B), C(1), C(0))
        else:
            assert B > 0
            r = B.sqrt()
            i, j = (C(0), C(A), C(1), C(0)), (C(r), C(0), C(0), C(-r))
        one = (C(1), C(0), C(0), C(1))
        return one, i, j, _multiply(i, j)

    def embed(self, x):
        r"""Return `\iota_v(x)`."""
        coefficients = [self.C(self.split(c)) for c in self.algebra(x).coefficient_tuple()]
        return tuple(sum(c * f[e] for c, f in zip(coefficients, self._frame)) for e in range(4))

    def _set_center(self, attempt):
        R, C = self.R, self.C
        p = C(R(2).sqrt() / 10 * (attempt + 1), 1 + R(3).sqrt() / 7 + R(attempt) / 11)
        self.p = p
        conjugate = p.conjugate()
        det = p - conjugate
        self._M = (C(1), -p, C(1), -conjugate)
        self._M_inverse = (-conjugate / det, p / det, C(-1) / det, C(1) / det)
        self._disc_cache = {}

    def disc(self, x):
        r"""Return `\gamma_D = M \iota_v(x) M^{-1}`."""
        key = tuple(self.lattices.coordinates(x))
        if key not in self._disc_cache:
            self._disc_cache[key] = _multiply(_multiply(self._M, self.embed(x)), self._M_inverse)
        return self._disc_cache[key]

    # The absolute reduced norm and its enumeration.
    def absolute_gram(self, basis):
        r"""Return the real Gram matrix of `N` on ``basis`` (module docstring, Enumeration)."""
        discs = [self.disc(x) for x in basis]
        m = len(basis)
        gram = matrix(self.R, m, m)
        for s in range(m):
            for t in range(s, m):
                split = sum((discs[s][e] * discs[t][e].conjugate()).real() for e in range(4))
                pairing = (basis[s] * basis[t].conjugate()).reduced_trace()
                ramified = sum((w(pairing) for w in self.ramified), self.R(0)) / 2
                gram[s, t] = gram[t, s] = split + ramified
        return gram

    def _short(self, basis, gram, bound):
        r"""Return ``(value, x)`` for every `x = \sum v_m b_m \ne 0` (one of `\pm x`) with `N(x) \le` ``bound``, sorted by value."""
        m = gram.nrows()
        least = min(e.real() for e in gram.change_ring(RDF).eigenvalues())
        assert least > 0, "the absolute reduced norm is positive definite"
        scale = ZZ(2) ** 64
        integral = matrix(ZZ, m, m, lambda s, t: (gram[s, t] * scale).round())
        # |v^t (G - G') v| <= (m / 2 scale) |v|^2 <= (m / 2 scale) v^t G v / least; double the margin for the eigenvalue estimate.
        slack = 1 + self.R(2 * m) / (scale * self.R(least))
        found = pari(integral).qfminim((bound * slack * scale).ceil())[2].sage()
        result = []
        for column in found.columns():
            v = vector(ZZ, column)
            x = sum((ZZ(c) * b for c, b in zip(v, basis)), self.algebra(0))
            result.append(((v.change_ring(self.R) * gram * v.change_ring(self.R)), x))
        result.sort(key=lambda pair: pair[0])
        return result

    # Isometric circles and the exterior domain.
    def _circle_center(self, g):
        _a, _b, c, d = self.disc(g)
        if c.abs() < self.epsilon:
            return None
        return -d / c

    def _clip(self, polygon, center, label):
        r"""Clip the convex Klein polygon by `\mathrm{Re}(w \bar C) \le 1`; each vertex carries the label of the edge starting there."""

        def h(w):
            return (w * center.conjugate()).real() - 1

        result = []
        k = len(polygon)
        for index in range(k):
            P, edge = polygon[index]
            Q, _next = polygon[(index + 1) % k]
            hp, hq = h(P), h(Q)
            if hp <= self.epsilon:
                result.append((P, edge))
                if hq > self.epsilon and hp < -self.epsilon:
                    result.append((P + (Q - P) * (hp / (hp - hq)), label))
                elif hq > self.epsilon:
                    result[-1] = (P, label)
            elif hq < -self.epsilon:
                result.append((P + (Q - P) * (hp / (hp - hq)), edge))
        cleaned = []
        for index, (w, edge) in enumerate(result):
            following = result[(index + 1) % len(result)][0]
            if len(result) > 1 and (following - w).abs() < self.epsilon:
                continue
            cleaned.append((w, edge))
        return cleaned

    def _cuts(self, polygon, center):
        return any((w * center.conjugate()).real() - 1 > self.epsilon for w, _edge in polygon)

    def _is_finite(self, polygon):
        return all(edge is not None and w.abs() < 1 - self.epsilon for w, edge in polygon)

    def _poincare(self, w):
        return w / (1 + (1 - w.abs() ** 2).sqrt())

    def _mobius(self, g, z):
        a, b, c, d = self.disc(g)
        return (a * z + b) / (c * z + d)

    def _angle(self, z, following, center_following, previous, center_previous):
        I = self.C(0, 1)
        t1 = I * (z - center_following)
        if (t1 * (following - z).conjugate()).real() < 0:
            t1 = -t1
        t2 = I * (z - center_previous)
        if (t2 * (previous - z).conjugate()).real() < 0:
            t2 = -t2
        return (t1 / t2).arg().abs()

    def _angles(self, polygon, sides):
        k = len(polygon)
        points = [self._poincare(w) for w, _edge in polygon]
        centers = [sides[edge][1] for _w, edge in polygon]
        return points, [self._angle(points[s], points[(s + 1) % k], centers[s], points[s - 1], centers[s - 1]) for s in range(k)]

    def _area_over_2pi(self, polygon, sides):
        pi = self.R.pi()
        _points, angles = self._angles(polygon, sides)
        return ((len(polygon) - 2) * pi - sum(angles)) / (2 * pi)

    # Reduction.
    def reduce_point(self, z, sides):
        r"""Return `(z', \delta)` with `z' = \delta z \in \overline{E(G)}` ([V] Alg. 4.3)."""
        delta = self.algebra(1)
        while True:
            best, best_value = None, 1 - self.epsilon
            for g, _center in sides:
                _a, _b, c, d = self.disc(g)
                value = (c * z + d).abs()
                if value < best_value:
                    best, best_value = g, value
            if best is None:
                return z, delta
            z = self._mobius(best, z)
            delta = best * delta

    def reduce_element(self, gamma, sides):
        r"""Return `\bar\gamma = \delta \gamma` with `\bar\gamma(0) \in \overline{E(G)}` ([V] Alg. 4.3)."""
        _z, delta = self.reduce_point(self._mobius(gamma, self.C(0)), sides)
        return delta * gamma

    # Shimizu's formula.
    def covolume_over_2pi(self):
        r"""Return `\mathrm{vol}(\Gamma \backslash H) / 2\pi` for the maximal order, exactly (module docstring, Stop condition)."""
        assert self.is_maximal(), "Shimizu's formula is applied to a maximal order"
        R, n = self.R, self.degree
        pi = R.pi()
        if self.center is QQ:
            L, field_discriminant = pari.lfuncreate(1), 1
        else:
            L, field_discriminant = pari.lfuncreate(pari(self.center.absolute_polynomial())), abs(self.center.discriminant())
        zeta = R(pari.lfun(L, 2, precision=self.precision).sage())
        phi = ZZ(1)
        for prime in self._finite_ramified():
            phi *= self._prime_norm(prime) - 1
        volume = 8 * pi / (4 * pi**2) ** n * R(field_discriminant) ** (R(3) / 2) * zeta * phi
        return self._rationalize(volume / (2 * pi))

    def _denominator(self):
        n = self.degree
        return lcm([ZZ(e) for e in range(2, 4 * n * n + 3) if euler_phi(2 * e) <= 2 * n])

    def _rationalize(self, value):
        denominator = self._denominator()
        exact = QQ((value * denominator).round()) / denominator
        assert (value - exact).abs() < self.tolerance, f"{value} is not a rational with denominator dividing {denominator}"
        return exact

    # The algorithm.
    def fundamental_domain(self):
        r"""Return the Dirichlet domain of `O^1` for a maximal order (module docstring)."""
        if not hasattr(self, "_domain"):
            exact = self.covolume_over_2pi()
            attempt = 0
            while True:
                self._set_center(attempt)
                domain = self._search(exact)
                if domain is not None:
                    self._domain = domain
                    break
                attempt += 1
        return self._domain

    def _search(self, exact):
        r"""Run the enumeration with reduction until the area equals ``exact``; ``None`` when `p` has a nontrivial stabilizer."""
        R, C = self.R, self.C
        one, minus = self.algebra(1), self.algebra(-1)
        gram = self.absolute_gram(self.order)
        sides = []
        polygon = [(C(-2, -2), None), (C(2, -2), None), (C(2, 2), None), (C(-2, 2), None)]
        seen = set()

        def add(g):
            nonlocal polygon
            if any(g == s or g == -s for s, _center in sides):
                return True
            center = self._circle_center(g)
            if center is None:
                return False
            sides.append((g, center))
            polygon = self._clip(polygon, center, len(sides) - 1)
            return True

        radius = R(1)
        while True:
            bound = (self.degree - 1) + 2 * radius.cosh() + self.tolerance
            for _value, x in self._short(self.order, gram, bound):
                if x.reduced_norm() != 1 or x == one or x == minus:
                    continue
                key = tuple(self.lattices.coordinates(x))
                negative = tuple(-c for c in key)
                if key in seen or negative in seen:
                    continue
                seen.add(key)
                center = self._circle_center(x)
                if center is None:
                    return None
                if self._cuts(polygon, center):
                    if not (add(x) and add(~x)):
                        return None
                reduced = self.reduce_element(x, sides)
                if reduced != one and reduced != minus:
                    center = self._circle_center(reduced)
                    if center is None:
                        return None
                    if self._cuts(polygon, center) and not (add(reduced) and add(~reduced)):
                        return None
            if self._is_finite(polygon):
                area = self._area_over_2pi(polygon, sides)
                assert area > exact - self.tolerance, f"area(E(G))/2pi = {area} < covolume/2pi = {exact}"
                if (area - exact).abs() < self.tolerance:
                    assert self._rationalize(area) == exact
                    return self._pairing(polygon, sides, exact)
            radius += 1

    def _pairing(self, polygon, sides, exact):
        r"""Return the side pairing, the elliptic cycles and the genus of the Dirichlet domain (module docstring, Side pairing)."""
        k = len(polygon)
        points, angles = self._angles(polygon, sides)
        elements = [sides[edge][0] for _w, edge in polygon]
        pairing = []
        for s in range(k):
            g = elements[s]
            inverse = ~g
            partner = [t for t in range(k) if elements[t] == inverse or elements[t] == -inverse]
            assert len(partner) == 1, f"the side on I({g}) has no paired side on I(g^-1)"
            t = partner[0]
            assert (self._mobius(g, points[s]) - points[(t + 1) % k]).abs() < self.tolerance
            assert (self._mobius(g, points[(s + 1) % k]) - points[t]).abs() < self.tolerance
            pairing.append(t)
        one, minus = self.algebra(1), self.algebra(-1)
        elliptic = [(elements[s], 2) for s in range(k) if pairing[s] == s]
        visited = set()
        two_pi = 2 * self.R.pi()
        for start in range(k):
            if start in visited:
                continue
            vertex, transformation, angle = start, one, self.R(0)
            while True:
                visited.add(vertex)
                angle += angles[vertex]
                transformation = elements[vertex] * transformation
                vertex = (pairing[vertex] + 1) % k
                if vertex == start:
                    break
            order = ZZ((two_pi / angle).round())
            assert (order * angle - two_pi).abs() < self.tolerance, f"vertex cycle angle sum {angle} is not 2 pi / e"
            if order == 1:
                continue
            power = one
            for exponent in range(1, order + 1):
                power *= transformation
                is_scalar = power == one or power == minus
                assert is_scalar == (exponent == order), f"the cycle transformation {transformation} does not have order {order}"
            elliptic.append((transformation, order))
        genus = (exact - sum(1 - QQ((1, e)) for _g, e in elliptic) + 2) / 2
        assert genus in ZZ and genus >= 0, f"the genus {genus} from area {exact} and elliptic orders is not a natural number"
        return FuchsianFundamentalDomain(points, elements, pairing, elliptic, exact, ZZ(genus))

    # Groups.
    def _maximal(self):
        if not hasattr(self, "_maximal_domain"):
            if self.is_maximal():
                self._maximal_domain = self
            else:
                maximal = maximal_orders_containing(self.algebra, self.order, self.expected_discriminant())
                self._maximal_domain = FuchsianDomain(self.algebra, maximal[0], self.precision)
        return self._maximal_domain

    def sl1_generators(self):
        r"""Return generators of `O^1`: `-1` and the side pairings of the Dirichlet domain; for a non-maximal order, Schreier generators of
        the elements of `O'^1` lying in `O`."""
        maximal = self._maximal()
        generators = (self.algebra(-1), *maximal.fundamental_domain().sides)
        if maximal is not self:
            subgroup, _index = finite_preimage_generators(generators, lambda u: u, self.in_order)
            generators = (self.algebra(-1), *subgroup)
        assert all(self.in_order(u) and u.reduced_norm() == 1 for u in generators)
        return tuple(generators)

    def _element_of_norm(self, norm):
        r"""Return `x \in O` with `\mathrm{nrd}(x) =` ``norm``, by the enumeration of `N` in growing ellipsoids."""
        gram = self.absolute_gram(self.order)
        bound = sum((w(norm) for w in self.ramified), self.R(0)) + 2 * self.split(norm).abs() + 1
        while True:
            for _value, x in self._short(self.order, gram, bound):
                if x.reduced_norm() == norm:
                    return x
            bound *= 2

    def unit_generators(self):
        r"""Return generators of `O^\times` (module docstring, Units and the normalizer)."""
        maximal = self._maximal()
        if maximal is not self:
            subgroup, _index = finite_preimage_generators(maximal.unit_generators(), lambda u: u, self.in_order)
            result = (self.algebra(-1), *subgroup)
            assert all(self.is_unit(u) for u in result)
            return result
        result = [*self.sl1_generators(), *(self.algebra(u) for u in _central_units(self.center))]
        for norm in positive_unit_classes(self.center, self.ramified):
            if norm != 1:
                result.append(self._element_of_norm(norm))
        assert all(self.is_unit(u) for u in result)
        return tuple(result)

    def norm_generator(self, elements, target):
        r"""Return `x` in the `Z`-span of ``elements`` with `\mathrm{nrd}(x) O_F =` ``target``; for an invertible one-sided ideal of a maximal
        order such `x` generates it, and it exists exactly when ``target`` has a generator positive at the ramified real places."""
        field = self.center
        assert positive_generator(field, target, self.ramified) is not None, (
            f"{target} has no generator positive at the ramified real places, so no ideal of that reduced norm is principal (Eichler)"
        )
        basis = self.lattices.basis(elements)
        assert len(basis) == 4 * self.degree
        gram = self.absolute_gram(basis)
        bound = min(gram[s, s] for s in range(gram.nrows()))
        while True:
            for _value, x in self._short(basis, gram, bound):
                if field.ideal(x.reduced_norm()) == target:
                    return x
            bound *= 2

    def normalizer_generators(self):
        r"""Return generators of `N(O)` through :func:`quaternion_orders_number_field.normalizer_generators`."""
        from sage_indefinite_port.indefinite.quaternion_orders_number_field import normalizer_generators

        maximal = maximal_orders_containing(self.algebra, self.order, self.expected_discriminant())
        return normalizer_generators(self.algebra, self.order, maximal)


def _multiply(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return (a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h)
