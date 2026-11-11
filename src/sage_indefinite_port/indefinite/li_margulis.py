r"""Bounded searches of Li-Margulis, with exact verification (#55, #56).

Theorem 4 and Lemma 20 of Li-Margulis bound one member of every orbit of primitive embeddings
of a Gram matrix `A` into `L`. Given a box containing such members, this module lists the
primitive tuples in the box with Gram matrix `A` and keeps one per `O(L)`-orbit, by the exact
marked transporter. The result is a statement about the supplied box: the bound is the input.
"""

from dataclasses import dataclass
from itertools import product
from math import comb as binomial, factorial, isqrt

from sage.matrix.constructor import matrix
from sage.functions.other import ceil
from sage.rings.integer_ring import ZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ
from sage.quadratic_forms.quadratic_form import QuadraticForm

from sage_indefinite_port.indefinite.pairing_slices import marked_transporter


def generator_bound(gram, constant):
    r"""Return the integer box size of Theorem 2, `C_r |\det B|^{r^6} \|B\|^{3 r^4}`, rounded up and made strict.

    ``constant`` is a certified majorant of the dimension constant `C_r`. Section 1.7 of the
    source computes it from the proof and appendix, and gives no numerical value.
    """
    gram = matrix(ZZ, gram)
    rank = gram.nrows()
    if rank < 3 or gram.det() == 0:
        raise ValueError("Theorem 2 needs a nonsingular form of rank at least three")
    norm = max(abs(entry) for entry in gram.list())
    return ceil(constant * abs(gram.det()) ** (rank**6) * norm ** (3 * rank**4)) + 1


def isometry_box_bound(left, right, constant, anisotropic=False):
    r"""Return an integer strictly above the bound of Theorem 1 for `\gamma^t B \gamma = A`.

    With `A` = ``left``, `B` = ``right`` and `n` the rank, equation (3) is
    `C_n |\det B|^{-(13n^3 - 39n^2 + 36n + 52)/(20n)} \|A\|^{(13n^3 - 13n^2 + 52n + 20)/40}
    \|B\|^{(13n^3 - 13n^2 + 72n - 20)/40}`, and equation (4), for anisotropic `Q_B`, is
    `C_n |\det B|^{(13n^2 + 8n - 26)/(10n)} \|A\|^{1/2} \|B\|^{(n-1)/2}`. Each integer base
    is at least one, so raising it to the ceiling of its exponent is an upper bound, and the
    negative determinant power in (3) is at most one.
    """
    left, right = matrix(ZZ, left), matrix(ZZ, right)
    n = right.nrows()
    if n < 3 or left.det() == 0 or right.det() == 0:
        raise ValueError("Theorem 1 needs nonsingular forms of rank at least three")
    norm_left = max(abs(entry) for entry in left.list())
    norm_right = max(abs(entry) for entry in right.list())
    determinant = abs(right.det())
    if anisotropic:
        exponents = QQ(13 * n**2 + 8 * n - 26) / (10 * n), QQ(1) / 2, QQ(n - 1) / 2
    else:
        exponents = QQ(0), QQ(13 * n**3 - 13 * n**2 + 52 * n + 20) / 40, QQ(13 * n**3 - 13 * n**2 + 72 * n - 20) / 40
    return ceil(
        constant
        * determinant ** ceil(exponents[0])
        * norm_left ** ceil(exponents[1])
        * norm_right ** ceil(exponents[2])
    ) + 1


def isotropic_pair_grams(lattice):
    r"""Return the Gram matrices `[[0, d], [d, a]]` of primitive pairs `(v, w)` with `v` primitive isotropic.

    `d = div(v)` divides `exp(A_L)`, and `w -> w + k v` reduces `a` modulo `2d`. An even lattice has
    even `a`. Every primitive isotropic `v` occurs in such a pair with `b(v, w) = d`.
    """
    frame = tuple(lattice.module_generators())
    gram = matrix(ZZ, [[int(lattice.b(x, y)) for y in frame] for x in frame])
    exponent = max(gram.elementary_divisors())
    even = all(gram[i, i] % 2 == 0 for i in range(gram.nrows()))
    return tuple(
        ((0, int(d)), (int(d), a)) for d in exponent.divisors() for a in range(0, 2 * int(d)) if not even or a % 2 == 0
    )


def _vectors_in_box(lattice, bound):
    generators = tuple(lattice.module_generators())
    for coordinates in product(range(-bound, bound + 1), repeat=len(generators)):
        yield sum(
            (lattice.scalar_multiple(c, g) for c, g in zip(coordinates, generators, strict=True)), lattice.zero()
        )


def primitive_tuples_in_box(lattice, gram, bound):
    r"""Yield the primitive tuples with Gram matrix ``gram`` and coordinates of absolute value at most ``bound``."""
    gram = matrix(ZZ, gram)
    size = gram.nrows()
    shell = {}
    for index in range(size):
        target = gram[index, index]
        shell[index] = [v for v in _vectors_in_box(lattice, bound) if int(lattice.b(v, v)) == int(target)]

    def extend(prefix):
        index = len(prefix)
        if index == size:
            rows = matrix(ZZ, [[int(c) for c in v] for v in prefix])
            if rows.rank() == size and all(d == 1 for d in rows.elementary_divisors()):
                yield tuple(prefix)
            return
        for candidate in shell[index]:
            if all(int(lattice.b(prefix[i], candidate)) == int(gram[i, index]) for i in range(index)):
                yield from extend([*prefix, candidate])

    yield from extend([])


def embedding_orbits_in_box(lattice, gram, bound):
    r"""Return one tuple per ``O(L)``-orbit among the primitive tuples in the box."""
    representatives = []
    for candidate in primitive_tuples_in_box(lattice, gram, bound):
        if all(marked_transporter(lattice, known, candidate) is None for known in representatives):
            representatives.append(candidate)
    return tuple(representatives)


def is_isometric_in_box(source, target, bound):
    r"""Return an isometry `source -> target` with images in the box, or ``None`` if the box has none."""
    frame = tuple(source.module_generators())
    gram = [[int(source.b(x, y)) for y in frame] for x in frame]
    labels = tuple(source.module_generating_set())
    for images in primitive_tuples_in_box(target, gram, bound):
        rows = matrix(ZZ, [[int(c) for c in v] for v in images])
        if abs(rows.det()) == 1:
            table = dict(zip(labels, images, strict=True))
            return source.Hom(target)(lambda label: table[label])
    return None


Majorant = Rational  # an exact rational majorant, at least one


@dataclass(frozen=True)
class LiMargulisConstants:
    r"""Majorants of the multiplicative constants that Li-Margulis leave as `\ll_n` (rank ``rank``).

    Section 1.7 states that each `\ll_n` constant is computable from its proof and gives no value.
    The proof of Theorem 1 (section 6) and Theorem 2 (section 9.5) compose explicit factors, which
    this module evaluates, with the following implicit factors, which :func:`default_constants`
    derives. Each field is a rational number at least one that majorizes the constant in the quoted
    statement.

    - ``theorem5``: equation (13), `\|h\| \le K \operatorname{vol}_H(H.x)^{13/5}
      (\|g_1\| \|g_2\|)^{(13n/5)(1 + \dim H/2)}` for closed orbits, by :func:`theorem5_constant` from
      the corrected (12) and the `\delta \asymp_n \|g\|^{-n}` of Lemma 11.
    - ``theorem5_compact``: the compact case of Theorem 5, `\|h\| \le K \operatorname{vol}_H(H.x)^{13/5}`,
      which also uses `\eta_0` of Lemma 10.
    - ``theorem6``: equation (18) for `N = 1`, `\operatorname{vol}_H(H.[g]) \le K
      |\det Q|^{(\dim G - \dim H)/n}`. Its proof (section 5.6) uses `c_n`, `b_n`, `\eta_0` of
      Lemma 15 and `R_n` of Lemma 14, whose compact set is `\mathfrak S(\epsilon_n)`,
      `\epsilon_n = (4 n^3 6^n (n^3 + n)^{1/n^2})^{-n^2}` (appendix).
    - ``lemma18``: equation (30), `\|\gamma\| \le K |\det B|^{-(13n^3-39n^2+26n+52)/(20n)}
      \|B\|^{(13n^3-13n^2+62n)/20}`.
    - ``corollary2``: Corollary 2 at `(\alpha, \beta) = (2, 1)`, from the corrected Lemma 24, `\|\gamma\| \le K` for
      `\gamma \in (S_{2,1}^{-1} S_{2,1}) \cap GL_n(\mathbb Z)`.
    - ``corollary3``: corrected Corollary 3, `\|B\| \le K |\det B|^{n+1}` for `(2,1)`-reduced indefinite `B`.
    - ``lemma20``: equation (31), `\|\gamma^t B \gamma\| \le K (|\det A|^{n-m} |\det B| + \|A\|)`.
    """

    rank: int
    theorem5: Majorant
    theorem5_compact: Majorant
    theorem6: Majorant
    lemma18: Majorant
    corollary2: Majorant
    corollary3: Majorant
    lemma20: Majorant

    def __post_init__(self):
        if self.rank < 3:
            raise ValueError("the Li-Margulis bounds need rank at least three")
        for name in ("theorem5", "theorem5_compact", "theorem6", "lemma18", "corollary2", "corollary3", "lemma20"):
            if QQ(getattr(self, name)) < 1:
                raise ValueError(f"the majorant {name} must be at least one")


def _power_majorant(base, exponent):
    r"""Return `base^{\lceil exponent \rceil}`, a majorant of `base^{exponent}` for `base \ge 1`, `exponent \ge 0`."""
    return QQ(base) ** ceil(QQ(exponent))


def lemma1_constant(n):
    r"""Return an integer `K \ge n^{1 + n/2}`, the constant of (24): `\|g_A\| \le K |\det A|^{-1/(2n)} \|A\|^{1/2}`.

    Proof of Lemma 1 made explicit. Write `X = O_1^t D_1 O_1` and `Y = O_2^t D_2 O_2` with `O_i`
    orthogonal and `P` a permutation matching signs; `g = O_2^t |D_2|^{-1/2} P |D_1|^{1/2} O_1`
    satisfies `g^t Y g = X`. Every eigenvalue of `X` has absolute value at most `n \|X\|`, and the
    least eigenvalue of `Y` in absolute value is at least `|\det Y| / (n \|Y\|)^{n-1}`. With
    `\|MN\| \le n \|M\| \|N\|` for the entrywise maximum norm, `\|g\| \le n^{1+n/2}
    \|X\|^{1/2} \|Y\|^{(n-1)/2} |\det Y|^{-1/2}`. Equation (24) takes `X = |\det A|^{-1/n} A`
    and `Y = X_0` of (5), with `\|X_0\| = |\det X_0| = 1`.
    """
    square = n ** (n + 2)
    root = isqrt(square)
    return root if root * root == square else root + 1


def _theorem5_exponent(n):
    return QQ(13 * n) / 5 * (1 + QQ(n * (n - 1)) / 4)


def _codimension(n):
    r"""Return `\dim G - \dim H = (n^2 - 1) - n(n-1)/2 = \dim \mathfrak h'`."""
    return (n * n + n - 2) // 2


def lemma10_radius(n):
    r"""Return `\eta_0 = (16 n^2 \binom{n}{\lfloor n/2 \rfloor})^{-1}`, an injectivity radius of Lemma 10.

    Proof of Lemma 10 made explicit. For `u = 1 + E`, `\|E\| < \eta`, the coefficient `e_k(u)` of
    the characteristic polynomial is a sum of `\binom{n}{k}` principal minors `\det(1_k + E_S)`;
    a `j \times j` minor of `E` is at most `j! \eta^j`, so `|\det(1_k + E_S) - 1| \le \sum_j
    \binom{k}{j} j! \eta^j \le k\eta/(1 - k\eta)`. With `\eta = (4n\binom{n}{\lfloor n/2 \rfloor})^{-1}`
    each `|e_k(u) - \binom{n}{k}| < 1`, so an integral characteristic polynomial is `(x-1)^n` and `u`
    is unipotent. Lemma 9: for `h_i \in B_G(\eta_0)`, `\eta_0 \le 1/(2n)`, `\|h_2^{-1}\| \le 2`
    (Neumann series) and `\|h_2^{-1} h_1 - 1\| \le n \|h_2^{-1}\| \|h_1 - h_2\| < 4 n \eta_0 = \eta`.
    """
    return QQ(1) / (16 * n * n * binomial(n, n // 2))


def lemma11_constant(n):
    r"""Return `\kappa = 8 n^3 (n-1)!`: Lemma 11 holds with `\delta = (\kappa \|g\|^n)^{-1}`.

    Proof of Lemma 11 made explicit. `\gamma - 1 = g^{-1} M g` with `M = g\gamma g^{-1} - 1`, so
    `1 \le \|\gamma - 1\| \le n^2 \|g^{-1}\| \|M\| \|g\| \le n^2 (n-1)! \|g\|^n \|M\|` (Cramer,
    `\|g^{-1}\| \le (n-1)! \|g\|^{n-1}` for `\det g = 1`). For `g_1, g_2 \in B_G(\delta)` with
    `g_1.[g] = g_2.[g]`, `M = g_2^{-1} g_1 - 1` has `\|M\| < 4 n \delta` (as in :func:`lemma10_radius`).
    """
    return 8 * n**3 * factorial(n - 1)


def lemma15_inverse_radius():
    r"""Return 3, a majorant of `1/\eta_0` for `\eta_0 = \log(2)/2` in the proof of Lemma 15.

    `\mathfrak h' = \{X \in \mathfrak{sl}_n : X_0 X^t X_0 = X\}` is the Frobenius complement of
    `\mathfrak h`, and `\sigma(g) = X_0 (g^t)^{-1} X_0` fixes `O_{Q_0}(\mathbb R)` with
    `\sigma(\exp w) = \exp(-w)` on `\mathfrak h'`. If `\exp(-w_2)\exp(w_1) = s \in SO_{Q_0}(\mathbb R)`,
    then `g \sigma(g)^{-1}` for `g = \exp(w_1) = \exp(w_2) s` gives `\exp(2w_1) = \exp(2w_2)`, and
    `\log \exp X = X` for `\|X\|_{op} < \log 2`. So `\eta_0 = \log(2)/2` in Frobenius norm, and
    `2/\log 2 < 2.89`.
    """
    return QQ(3)


def zeta_product_majorant(n):
    r"""Return `n \ge \zeta(2)\cdots\zeta(n) = m_{G/\Gamma}(G/\Gamma)` of (8): `\zeta(k) \le k/(k-1)`."""
    return QQ(n)


def lemma15_constant_power(n):
    r"""Return a majorant of `c_n^d`, `d = \dim G - \dim H`, in the proof of Lemma 15: `4 n^2 d^{\lceil d/2 \rceil}`.

    From (22), `c_n^d = 4 \zeta(2)\cdots\zeta(n) / b_n`. The tube map `(w, h) \mapsto \exp(w) h`
    has, for `m_H` and Lebesgue measure on `\mathfrak h'`, Jacobian `\det(P_{\mathfrak h'}
    \phi(\operatorname{ad} w)|_{\mathfrak h'})`, `\phi(z) = (1 - e^{-z})/z`; for `\|w\| < \eta_0`,
    `\|\operatorname{ad} w\|_{op} < \log 2` and `\|\phi(\operatorname{ad} w) - 1\| \le 1/\log 2 - 1 < 1/2`,
    so the Jacobian is at least `2^{-d}`. At the identity `m_G` of (7) is `n^{-1/2}` times
    Lebesgue measure on `\mathfrak{sl}_n` (the factor `t^{1/n} 1_n` has Frobenius speed
    `n^{-1/2}`). The unit `d`-ball has volume at least `(2/\sqrt d)^d`. So
    `b_n = n^{-1/2} d^{-d/2}`, `c_n^d \le 4 n \cdot n^{1/2} d^{d/2}`.
    """
    d = _codimension(n)
    return QQ(4 * n * n * d ** ((d + 1) // 2))


def appendix_inverse_epsilon(n):
    r"""Return `\epsilon_n^{-1} = (2 a_n)^{n^2} = (4 n^3 6^n)^{n^2} (n^3 + n)`, appendix, after (32)."""
    return ZZ(4 * n**3 * 6**n) ** (n * n) * (n**3 + n)


def lemma14_constant(n):
    r"""Return an integer `R_n` of Lemma 14 for `\Omega_{cpt} = \mathfrak S(\epsilon_n)` of the appendix.

    A unimodular lattice with `\lambda_1 \ge \epsilon_n` has `\lambda_n \le \lambda_1 \cdots
    \lambda_n / \lambda_1^{n-1} \le n^{n/2} \epsilon_n^{1-n}` (Minkowski's second theorem,
    `2^n / V_n \le n^{n/2}` since the ball contains the cube of side `2/\sqrt n`), and a basis
    with `|b_i| \le \max(1, i/2) \lambda_i \le (n/2) \lambda_n` [Cassels, Geometry of Numbers, V,
    Lemma 8]; one sign change puts it in `SL_n`. Entries are at most `|b_i|`.
    """
    return ZZ(n ** ((n + 1) // 2 + 1) * appendix_inverse_epsilon(n) ** (n - 1)) // 2 + 1


def theorem6_constant(n):
    r"""Return a majorant of the constant of (18) for `N = 1`: `(9 n^2 R_n^2)^d \cdot 4 n^2 d^{\lceil d/2 \rceil}`.

    Section 5.6 made explicit, `d = \dim G - \dim H`, `\eta_1 = c_n \operatorname{vol}^{-1/d}`.
    If `\eta_1 < \eta_0`, Lemma 15 gives `u = \exp(-w_2)\exp(w_1)` with `\|w_i\| < \eta_1`, so
    `\|u - 1\| \le e^{2\eta_1} - 1 \le 2\eta_1/\log 2 < 3 \eta_1`, and (23) gives
    `\operatorname{vol} < (9 n^2 R_n^2 c_n)^d |\det Q|^{d/n}`. Otherwise `\operatorname{vol} \le
    (c_n/\eta_0)^d \le (3 c_n)^d`, a smaller bound since `|\det Q| \ge 1`; this replaces the
    finiteness argument of Lemma 15.
    """
    d = _codimension(n)
    return QQ(9 * n * n * lemma14_constant(n) ** 2) ** d * lemma15_constant_power(n)


def corollary1_constant(k):
    r"""Return `C` with `\|\gamma^t A \gamma\| \le C |\det A|` for some `\gamma`, every nonsingular integral `A` of rank `k`.

    Corollary 1, made explicit per signature. Definite: an LLL-reduced basis has
    `\prod a_{ii} \le 2^{k(k-1)/2} |\det A|` [Lenstra-Lenstra-Lovasz 1982, (1.8)], integral
    `a_{ii} \ge 1` and `|a_{ij}| \le \sqrt{a_{ii} a_{jj}}`. Indefinite binary: a Gauss-reduced form
    has entries below `2 |\det A|^{1/2}`. Indefinite, `k \ge 3`: the proof of Theorem 3 gives
    `\|A_1\| < k R_k^2 |\det A|^{1/k}`.
    """
    definite = ZZ(2) ** (k * (k - 1) // 2)
    if k == 1:
        return QQ(1)
    if k == 2:
        return QQ(max(definite, 2))
    return QQ(max(definite, k * lemma14_constant(k) ** 2))


def lemma20_constant(n):
    r"""Return a majorant of the constant of (31), uniform in `1 \le m < n`.

    Proof of Lemma 20 made explicit, `k = n - m`, `\lambda = \det A`. `T^t E T = \lambda^{-1}
    T^t (\lambda E) T` has norm at most `C_1(k) |\lambda|^k |\det B|` by Corollary 1, since
    `|\det \lambda E| \le |\lambda|^{k+1} |\det B|`. Rounding `r` makes `\|s\| \le 1/2`, so
    `\|s A\| \le (m/2) \|A\|` and `\|s A s^t\| \le (m^2/4) \|A\|`.
    """
    return max(max(corollary1_constant(n - m), QQ(m * m) / 4, QQ(1)) for m in range(1, n))


def section6_constant(n, theorem5, theorem6):
    r"""Return `n^2 (n-1)! K_{24}^{n + 2e} K_{(13)} K_{(18)}^{13/5}`, the composite of section 6.

    It majorizes `C_n` of (3) and the constant of (30): Lemma 18 bounds `g^{-1} h s g` with
    `\|s\| = 1` and `\|h\|` from (13), (18) and (28) = (24), the same composite.
    """
    exponent = n + 2 * _theorem5_exponent(n)
    return (
        n**2
        * factorial(n - 1)
        * _power_majorant(lemma1_constant(n), exponent)
        * QQ(theorem5)
        * _power_majorant(theorem6, QQ(13) / 5)
    )


def matrix_coefficient_constant():
    r"""Return `16/5 > \pi \coth \pi`, an explicit `K` for the corrected (12).

    Corrected (12): for a unitary representation `\pi` of `SL_2(\mathbb R)` that is strongly
    `L^{64/25+\epsilon}` for every `\epsilon > 0`, with `W` the generator of `K = SO(2)`
    (`\exp(\theta W)` is rotation by `\theta`) and `S(v)^2 = \|v\|^2 + \|W v\|^2`,
    `|\langle \pi(g) v, w \rangle| \le \pi \coth(\pi)\, \Xi(g)^{25/32} S(v) S(w)`.
    The exponent `25/32` of the printed statement is unchanged.

    Proof. Decompose `\pi` as a direct integral of irreducibles; a.e. component is tempered or a
    complementary series `I(\theta\rho)` with `\phi_{\theta\rho}(a_t) \asymp e^{-(1-\theta)t}`
    in `L^q` iff `q(1-\theta) > 2`, so `1 - \theta \ge 25/32`. Each `K`-type `m \in \mathbb Z` is
    one-dimensional. Tempered: `|\langle \pi(g) v_m, w_n \rangle| \le \Xi(g) \|v_m\| \|w_n\|`
    [Cowling-Haagerup-Howe 1988, Thm 2, with `\dim \langle K v \rangle = 1`]. Complementary: in the
    compact picture the unitary product is `\langle A f_1, f_2 \rangle_{L^2(K)}` with `A e_m = c_m e_m`,
    `c_m > 0`, and `A \pi_{\lambda}(g) = \pi_{-\lambda}(g) A`. The operator `\pi_{-\lambda}(g)`
    multiplies by a positive function after a change of variable, so `|\langle \pi_{-\lambda}(g) e_m,
    e_n \rangle_{L^2}| \le \langle \pi_{-\lambda}(g) 1, 1 \rangle = \phi_\lambda(g)` and the normalized
    coefficient is at most `\sqrt{c_m / c_n}\, \phi_\lambda(g)`; applied to `g^{-1}` it is at most
    `\sqrt{c_n / c_m}\, \phi_\lambda(g)`, hence at most `\phi_\lambda(g)`. Hoelder in `\lambda` for
    `\phi_\lambda = \int_K e^{(\lambda - \rho) H}` with `\phi_0 = \Xi`, `\phi_\rho = 1` gives
    `\phi_{\theta\rho} \le \Xi^{1-\theta} \le \Xi^{25/32}`. Integrating over components and summing
    over `K`-types, `|\langle \pi(g) v, w \rangle| \le \Xi(g)^{25/32} \sum_m \|v_m\| \sum_n \|w_n\|`,
    and Cauchy-Schwarz with `\sum_{m \in \mathbb Z} (1 + m^2)^{-1} = \pi \coth \pi < 3.1534` bounds each
    sum by `(\pi \coth \pi)^{1/2} S`.
    """
    return QQ(16) / 5


def harish_chandra_decay_constant():
    r"""Return 21: `\Xi(g_t)^{25/32} < 21 e^{-5|t|/13}` for `g_t` of (11), the explicit form of (16).

    `\Xi(g_t) = (2\pi)^{-1} \int_0^{2\pi} (e^t \cos^2\theta + e^{-t} \sin^2\theta)^{-1/2} d\theta`. For
    `t \ge 0` the integrand is at most `\min(e^{t/2}, e^{-t/2} / \cos\theta)` and `\cos\theta \ge
    2s/\pi` with `s = \pi/2 - \theta`, so `\Xi(g_t) \le (1 + |t|) e^{-|t|/2}`. Then `\Xi^{25/32}
    e^{5|t|/13} \le ((1 + |t|) e^{-|t|/130})^{25/32} \le (130 e^{-129/130})^{25/32}`, and
    `130 e^{-129/130} < 130 (10/27)(66/65) = 1320/27 < 21^{32/25}` (`e > 2.7`, `e^{1/130} < 66/65`).
    """
    return QQ(21)


def _dim_h(n):
    r"""Return `\dim H = \dim SO(p, q) = n(n-1)/2`."""
    return n * (n - 1) // 2


def haar_ball_inverse_constant(n):
    r"""Return `c^{-1} = 2^{D} D^{\lceil D/2 \rceil}`, `D = \dim H`: `m_H(B_H(\delta)) \ge c \delta^D` for `\delta \le 1/(2n)`.

    With `m_H` of (9): `B_H(\delta)` contains `\exp` of the Frobenius ball of radius `r = \delta/2` in
    `\mathfrak h`, since `\|\exp X - 1\| \le e^{\|X\|_{op}} - 1 \le e^{\delta/2} - 1 < \delta`; `\exp` is
    injective there (`r < \log 2`). The eigenvalues of `\operatorname{ad}_{\mathfrak h} X` have modulus at
    most `2r \le 1/(2n)`, and `|\phi(z) - 1| \le |z|` for `\phi(z) = (1 - e^{-z})/z`, so the Jacobian of
    (9) is at least `2^{-D}`. The ball contains a cube of side `2r/\sqrt D`. Hence
    `m_H(B_H(\delta)) \ge (\delta / (2\sqrt D))^D`.
    """
    d = _dim_h(n)
    return QQ(2**d * d ** ((d + 1) // 2))


def _theorem5_base(n, inverse_radius):
    return harish_chandra_decay_constant() * matrix_coefficient_constant() * haar_ball_inverse_constant(n) * QQ(
        inverse_radius
    ) ** _dim_h(n)


def theorem5_constant(n):
    r"""Return a majorant of the constant of (13): `4 n^4 n^{\lceil 13n/5 \rceil} (21 K c^{-1} \kappa^{D})^3`.

    Section 4.4 made explicit, with `K` of :func:`matrix_coefficient_constant`, `c` of
    :func:`haar_ball_inverse_constant`, `\kappa` of :func:`lemma11_constant`, `D = \dim H`,
    `V = \operatorname{vol}_H(H.x)`. Lemma 12 is replaced by `K`-invariant bumps, so no derivative enters.
    Put `\delta_i = (\kappa \|g_i\|^n)^{-1} \le 1/(2n)` (Lemma 11; `\|g\| \ge n^{-1/2}` by Hadamard),
    `\psi_i = \int_{\iota(K)} k.\mathbf 1_{B_H(\delta_i).x_i} dk` and `b_i = m_H(B_H(\delta_i)) \le V`. By (10),
    `\mu(\psi_i) = b_i / V` and `\|\psi_i\| \le (b_i/V)^{1/2}`; `\psi_i - \mu(\psi_i)` is `K`-fixed in
    `L^2_0`, strongly `L^{64/25+\epsilon}` by Lemmas 2 and 4, so the corrected (12) gives
    `\langle a_t \psi_1, \psi_2 \rangle \ge b_1 b_2 / V^2 - K \Xi(g_t)^{25/32} (b_1 b_2)^{1/2} / V`.
    With `X = 21 K V (b_1 b_2)^{-1/2} \ge 21 K` and `e^t = X^{13/5}` this is positive by
    :func:`harish_chandra_decay_constant`, so `h = h_2^{-1} k_2^{-1} a_t k_1 h_1` with `h_i \in B_H(\delta_i)`,
    `k_i \in \iota(K)` (rotations, entries at most one) carries `x_1` to `x_2`, and
    `\|h\| \le n^4 \cdot 2 \cdot \cosh t \cdot 2 \le 4 n^4 X^{13/5}` (`\|h_2^{-1}\| \le 2` by Neumann).
    `X \le 21 K c^{-1} \kappa^D V (\|g_1\| \|g_2\|)^{nD/2}`, and `(\|g_1\| \|g_2\|)^{-13n/5} \le n^{13n/5}`
    converts the exponent `(13n/5)(D/2)` to the `(13n/5)(1 + D/2)` of (13).
    """
    return 4 * n**4 * _power_majorant(n, QQ(13 * n) / 5) * _power_majorant(_theorem5_base(n, lemma11_constant(n)), QQ(13) / 5)


def theorem5_compact_constant(n):
    r"""Return a majorant of the compact case of Theorem 5: `4 n^4 (21 K c^{-1} \eta_0^{-D})^3`.

    As :func:`theorem5_constant` with `\delta_i = \eta_0` of :func:`lemma10_radius` for both points.
    """
    return 4 * n**4 * _power_majorant(_theorem5_base(n, 1 / lemma10_radius(n)), QQ(13) / 5)


def lemma24_constant(n, alpha, beta, index):
    r"""Return `n^n \beta'^n (1+\beta)^{n-1} d^{n+1} \alpha^{3n^2-2n+1}`, `\beta' = \max(1, \beta)`, `d` = ``index``.

    Corrected Lemma 24 (with Lemma 25 as the case `d > 1`). Let `Q = (au)^t (au)` be
    `(\alpha, \beta)`-reduced, `\alpha \ge 1`, `\Lambda \subseteq \mathbb Z^n` of index `d`, and
    `v_1, \dots, v_n \in \Lambda` successive-minima vectors chosen greedily: `v_k` minimizes `Q` on
    `\Lambda \setminus \operatorname{span}(v_1, \dots, v_{k-1})`. Then `\|v_k\| \le` this constant.

    Proof. Write `y = u v`, `\lambda_k = Q(v_k)`. (i) Any `i` independent integral vectors include one with
    top nonzero coordinate `j \ge i`, and `Q \ge a_j^2 \ge \alpha^{-2(n-1)} a_i^2` on it; Minkowski's
    second theorem `\prod \lambda_i \le n^n d^2 \prod a_i^2` (`2^n/V_n \le n^{n/2}`) then gives
    `\lambda_k \le C a_k^2`, `C = n^n d^2 \alpha^{2(n-1)^2}`. (ii) Call `l` a gap if `a_l < \theta a_{l+1}`,
    `\theta^{-1} = d \sqrt n \beta' \alpha^{2n}`. Then `Q(d e_i) \le d^2 n \beta'^2 \alpha^{2n} a_l^2 < Q(w)`
    for `i \le l` and every integral `w \notin E_l = \langle e_1, \dots, e_l \rangle`, so
    `v_1, \dots, v_l` span `E_l`, and for `k > l` greedy choice makes `v_k` minimal on
    `v_k + d\mathbb Z^l`. (iii) Bounds on `|y_m|`: for `m \ge k`, `\sqrt{\lambda_k}/a_m \le \sqrt C
    \alpha^n`; for `m < k` with no gap in `[m, k)`, `\le \sqrt C \theta^{-n}`; otherwise let `l \ge m` be
    the least gap, `l < k`; size reduction by `d e_i`, `i \le l`, gives a point of the coset with
    `|y_i| \le d/2`, so `a_m^2 y_m^2 \le (d^2/4) \sum_{i \le l} a_i^2` and `|y_m| \le (d\sqrt n/2)
    \alpha^n \theta^{-n}`. All three are at most `Y = \sqrt C \theta^{-n} = n^n d^{n+1} \beta'^n
    \alpha^{3n^2-2n+1}`, and `v = u^{-1} y` gives `|v_m| \le Y (1+\beta)^{n-m}`.
    """
    alpha, beta, index = QQ(alpha), QQ(beta), ZZ(index)
    return (
        QQ(n) ** n
        * max(QQ(1), beta) ** n
        * (1 + beta) ** (n - 1)
        * index ** (n + 1)
        * alpha ** (3 * n * n - 2 * n + 1)
    )


def corollary2_constant(n):
    r"""Return `n! P^n`, `P` = :func:`lemma24_constant` at `(2, 1)`, index one: Corollary 2 at `(2, 1)`.

    `\gamma = \tau_2 \tau_1^{-1}` with `\tau_1 = (w_i)` greedy successive-minima vectors of `Q_1` and
    `\tau_2 = (\gamma w_i)`, greedy for `Q_2`; `\|\tau_1^{-1}\| \le (n-1)! \|\tau_1\|^{n-1}` (Cramer,
    `\det \tau_1 \in \mathbb Z \setminus 0`) and `\|\gamma\| \le n \|\tau_2\| \|\tau_1^{-1}\|`.
    """
    return factorial(n) * lemma24_constant(n, 2, 1, 1) ** n


def corollary3_constant(n):
    r"""Return `n! R P^{n-1}` with `\|B\| \le` it `\cdot |\det B|^{n+1}`: corrected Corollary 3.

    `P` = :func:`lemma24_constant` at `(2, 1)`, index one, and `R d^{n+1}` = :func:`lemma24_constant`
    at `(2, (n-1)!)`, index `d = |\det B|`. With `W = (w_i)` greedy for the `(2,1)`-reduced Hermite
    majorant `A` and `T = JB`, `C = JA^{-1}J` (Lemmas 22, 23), `T^t C T = A`, so `T w_i` are greedy
    successive-minima vectors of `Q_C` on `T \mathbb Z^n`, of index `d`. `\|B\| = \|T\| \le n \|TW\|
    \|W^{-1}\|`, `\|W^{-1}\| \le (n-1)! P^{n-1}`. The exponent `n+1` replaces the printed `3n-1`.
    """
    unit = lemma24_constant(n, 2, 1, 1)
    scaled = lemma24_constant(n, 2, factorial(n - 1), 1)
    return factorial(n) * scaled * unit ** (n - 1)


def default_constants(n):
    r"""Return the constants of rank ``n``, every majorant computed by this module."""
    theorem6 = theorem6_constant(n)
    theorem5 = theorem5_constant(n)
    return LiMargulisConstants(
        n,
        theorem5,
        theorem5_compact_constant(n),
        theorem6,
        section6_constant(n, theorem5, theorem6),
        corollary2_constant(n),
        corollary3_constant(n),
        lemma20_constant(n),
    )


def theorem_one_constant(constants):
    r"""Return a majorant of `C_n` in equation (3) of Theorem 1, for indefinite forms.

    Section 6: `\gamma = g_B^{-1} h g_A` (after `g_B \mapsto s_0 g_B` with `s_0` a signed diagonal
    matrix, which keeps `\|g_B\|`; the case `\det \gamma_0 = -1` multiplies by a signed
    permutation `\xi^{-1}`, which keeps the norm). Hence
    `\|\gamma\| \le n^2 \|g_B^{-1}\| \|h\| \|g_A\|` with `\|g_B^{-1}\| \le (n-1)! \|g_B\|^{n-1}`
    (Cramer, section 1.7), `\|h\|` from (13), the volume from (18) and `\|g_A\|, \|g_B\|` from (24).
    The exponents of `|\det B|, \|A\|, \|B\|` then sum to those of (3), and the constant is
    `n^2 (n-1)! K_{24}^{n + 2e} K_{(13)} K_{(18)}^{13/5}`, `e = (13n/5)(1 + \dim H / 2)`.
    """
    return section6_constant(constants.rank, constants.theorem5, constants.theorem6)


def theorem_one_anisotropic_constant(constants):
    r"""Return a majorant of `C_n` in equation (4) of Theorem 1, for anisotropic indefinite `Q_B`.

    As :func:`theorem_one_constant`, with the compact case of Theorem 5 (no `\|g_i\|` factor):
    `n^2 (n-1)! K_{24}^n K_{(13), cpt} K_{(18)}^{13/5}`.
    """
    n = constants.rank
    return (
        n**2
        * factorial(n - 1)
        * QQ(lemma1_constant(n)) ** n
        * QQ(constants.theorem5_compact)
        * _power_majorant(constants.theorem6, QQ(13) / 5)
    )


def theorem_two_constant(constants):
    r"""Return a majorant of `C_n` in Theorem 2, from the proof in section 9.5.

    With `a, b` the exponents of `\|A\|, \|B\|` in (3) and `e_{30} = (13n^3-13n^2+62n)/20`:
    `\|\tau_i\| \le C_1 \|B\|^a K_{C3}^{b} |\det B|^{(n+1)b}` (Theorem 1 for `\tau_i^t A_i \tau_i = B`
    and Corollary 3 for `A_i`; the negative power of `|\det A_i| = |\det B|` is at most one),
    `\|\sigma_j\| \le K_{30} \|B\|^{e_{30}}` by (30), `G = n \|\tau_i\| \|\sigma_j\|` bounds `\gamma_l`,
    `\|\xi\| \le n^2 (n-1)! K_{C2} G^n` by Corollary 2 and Cramer, and each generator
    `\sigma_j \xi` has norm at most `n K_{30} \|B\|^{e_{30}} \|\xi\|`. The resulting exponents
    `n(n+1)b \le n^6` and `e_{30} + n(a + e_{30}) \le 3n^4` are checked, so the bound has the form
    of Theorem 2 with constant `n^3 (n-1)! K_{C2} K_{30} (n C_1 K_{C3}^{b} K_{30})^n`.
    """
    n = constants.rank
    a = QQ(13 * n**3 - 13 * n**2 + 52 * n + 20) / 40
    b = QQ(13 * n**3 - 13 * n**2 + 72 * n - 20) / 40
    e30 = QQ(13 * n**3 - 13 * n**2 + 62 * n) / 20
    assert n * (n + 1) * b <= n**6
    assert e30 + n * (a + e30) <= 3 * n**4
    inner = n * theorem_one_constant(constants) * _power_majorant(constants.corollary3, b) * QQ(constants.lemma18)
    return n**3 * factorial(n - 1) * QQ(constants.corollary2) * QQ(constants.lemma18) * inner**n


def _gram(lattice):
    frame = tuple(lattice.module_generators())
    return matrix(ZZ, [[int(lattice.b(x, y)) for y in frame] for x in frame])


def _require_indefinite(gram):
    positive, negative, zero = QuadraticForm(QQ, 2 * gram.change_ring(QQ)).signature_vector()
    if zero or not positive or not negative:
        raise ValueError("the Li-Margulis route of section 6 needs a nonsingular indefinite form")


def _check_rank(gram, constants):
    if gram.nrows() != constants.rank:
        raise ValueError("the constants are for another rank")


def isometries_in_box(lattice, bound):
    r"""Return every isometry of ``lattice`` whose basis images have coordinates at most ``bound``."""
    gram = _gram(lattice)
    labels = tuple(lattice.module_generating_set())
    found = []
    for images in primitive_tuples_in_box(lattice, gram, bound):
        if abs(matrix(ZZ, [[int(c) for c in v] for v in images]).det()) == 1:
            table = dict(zip(labels, images, strict=True))
            found.append(lattice.Hom(lattice)(lambda label, table=table: table[label]))
    return tuple(found)


def orthogonal_group_generators(lattice, constants):
    r"""Return generators of the full `O(L)` by Theorem 2: every isometry in the theorem's box."""
    gram = _gram(lattice)
    _check_rank(gram, constants)
    _require_indefinite(gram)
    return isometries_in_box(lattice, generator_bound(gram, theorem_two_constant(constants)))


def isometry_decision(source, target, constants):
    r"""Return an isometry ``source -> target``, or ``None`` when none exists, by Theorem 1, (3)."""
    left, right = _gram(source), _gram(target)
    _check_rank(right, constants)
    _require_indefinite(right)
    if left.nrows() != right.nrows():
        return None
    return is_isometric_in_box(source, target, isometry_box_bound(left, right, theorem_one_constant(constants)))


def embedding_orbit_bound(gram_b, gram_a, constants):
    r"""Return a box meeting every `O(B)`-orbit of primitive embeddings of `A` into `B` (issue 56, B).

    For `m < n`, Lemma 20 completes a primitive `X` to `T = (X | Y)` with `\det T = 1` and
    `\|T^t B T\| \le K_{20} (|\det A|^{n-m} |\det B| + \|A\|)`; Theorem 1 gives `T'` with
    `T'^t B T' = T^t B T` and `\|T'\| < C_1 \|T^t B T\|^a \|B\|^b`; and `T' T^{-1} \in O(B)` carries
    `X` to the first `m` columns of `T'`. For `m = n` this is Theorem 1 itself.
    """
    gram_b, gram_a = matrix(ZZ, gram_b), matrix(ZZ, gram_a)
    _check_rank(gram_b, constants)
    _require_indefinite(gram_b)
    n, m = gram_b.nrows(), gram_a.nrows()
    if gram_a.det() == 0 or m > n:
        raise ValueError("the source Gram matrix must be nonsingular of rank at most n")
    if m == n:
        return isometry_box_bound(gram_a, gram_b, theorem_one_constant(constants))
    norm_a = max(abs(entry) for entry in gram_a.list())
    norm_b = max(abs(entry) for entry in gram_b.list())
    completed = QQ(constants.lemma20) * (abs(gram_a.det()) ** (n - m) * abs(gram_b.det()) + norm_a)
    a = QQ(13 * n**3 - 13 * n**2 + 52 * n + 20) / 40
    b = QQ(13 * n**3 - 13 * n**2 + 72 * n - 20) / 40
    return ceil(theorem_one_constant(constants) * _power_majorant(completed, a) * _power_majorant(norm_b, b)) + 1


def embedding_orbits(lattice, gram, constants):
    r"""Return one primitive embedding of ``gram`` per `O(L)`-orbit, complete by :func:`embedding_orbit_bound`."""
    return embedding_orbits_in_box(lattice, gram, embedding_orbit_bound(_gram(lattice), gram, constants))


def isotropic_vector_orbits_in_box(lattice, bound):
    r"""Return one primitive isotropic vector per `O(L)`-orbit among first vectors of pairs in the box.

    Issue 56, D: every primitive isotropic `v` is the first vector of a primitive pair with Gram
    matrix in :func:`isotropic_pair_grams`. The pair orbits are reduced to vector orbits by the
    marked transporter of the single vector.
    """
    representatives = []
    for gram in isotropic_pair_grams(lattice):
        for pair in embedding_orbits_in_box(lattice, gram, bound):
            vector = pair[0]
            if all(marked_transporter(lattice, (known,), (vector,)) is None for known in representatives):
                representatives.append(vector)
    return tuple(representatives)


def isotropic_vector_orbits(lattice, constants):
    r"""Return one primitive isotropic vector per `O(L)`-orbit, with the box of every pair Gram matrix."""
    gram_b = _gram(lattice)
    bound = max(embedding_orbit_bound(gram_b, gram, constants) for gram in isotropic_pair_grams(lattice))
    return isotropic_vector_orbits_in_box(lattice, bound)
