r"""Bounded searches of Li-Margulis, with exact verification (#55, #56).

Theorem 4 and Lemma 20 of Li-Margulis bound one member of every orbit of primitive embeddings
of a Gram matrix `A` into `L`. Given a box containing such members, this module lists the
primitive tuples in the box with Gram matrix `A` and keeps one per `O(L)`-orbit, by the exact
marked transporter. The result is a statement about the supplied box: the bound is the input.
"""

from dataclasses import dataclass
from itertools import product
from math import factorial, isqrt

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
    this module evaluates, with the following implicit factors, which the caller supplies. Each
    field is a rational number at least one that majorizes the constant in the quoted statement.

    - ``theorem5``: equation (13), `\|h\| \le K \operatorname{vol}_H(H.x)^{13/5}
      (\|g_1\| \|g_2\|)^{(13n/5)(1 + \dim H/2)}` for closed orbits. Its proof (section 4.4) uses
      the constant of (12) [Shalom, Thm 2.1; Kontorovich-Oh, Prop 3.2], the `\delta \asymp_n
      \|g\|^{-n}` of Lemma 11 and the bump-function constant of Lemma 12.
    - ``theorem5_compact``: the compact case of Theorem 5, `\|h\| \le K \operatorname{vol}_H(H.x)^{13/5}`,
      which also uses `\eta_0` of Lemma 10.
    - ``theorem6``: equation (18) for `N = 1`, `\operatorname{vol}_H(H.[g]) \le K
      |\det Q|^{(\dim G - \dim H)/n}`. Its proof (section 5.6) uses `c_n`, `b_n`, `\eta_0` of
      Lemma 15 and `R_n` of Lemma 14, whose compact set is `\mathfrak S(\epsilon_n)`,
      `\epsilon_n = (4 n^3 6^n (n^3 + n)^{1/n^2})^{-n^2}` (appendix).
    - ``lemma18``: equation (30), `\|\gamma\| \le K |\det B|^{-(13n^3-39n^2+26n+52)/(20n)}
      \|B\|^{(13n^3-13n^2+62n)/20}`.
    - ``corollary2``: Corollary 2 at `(\alpha, \beta) = (2, 1)`, `\|\gamma\| \le K` for
      `\gamma \in (S_{2,1}^{-1} S_{2,1}) \cap GL_n(\mathbb Z)`.
    - ``corollary3``: Corollary 3, `\|B\| \le K |\det B|^{3n-1}` for `(2,1)`-reduced indefinite `B`.
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
    n = constants.rank
    k24 = lemma1_constant(n)
    exponent = n + 2 * _theorem5_exponent(n)
    return (
        n**2
        * factorial(n - 1)
        * _power_majorant(k24, exponent)
        * QQ(constants.theorem5)
        * _power_majorant(constants.theorem6, QQ(13) / 5)
    )


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
    `\|\tau_i\| \le C_1 \|B\|^a K_{C3}^{b} |\det B|^{(3n-1)b}` (Theorem 1 for `\tau_i^t A_i \tau_i = B`
    and Corollary 3 for `A_i`; the negative power of `|\det A_i| = |\det B|` is at most one),
    `\|\sigma_j\| \le K_{30} \|B\|^{e_{30}}` by (30), `G = n \|\tau_i\| \|\sigma_j\|` bounds `\gamma_l`,
    `\|\xi\| \le n^2 (n-1)! K_{C2} G^n` by Corollary 2 and Cramer, and each generator
    `\sigma_j \xi` has norm at most `n K_{30} \|B\|^{e_{30}} \|\xi\|`. The resulting exponents
    `n(3n-1)b \le n^6` and `e_{30} + n(a + e_{30}) \le 3n^4` are checked, so the bound has the form
    of Theorem 2 with constant `n^3 (n-1)! K_{C2} K_{30} (n C_1 K_{C3}^{b} K_{30})^n`.
    """
    n = constants.rank
    a = QQ(13 * n**3 - 13 * n**2 + 52 * n + 20) / 40
    b = QQ(13 * n**3 - 13 * n**2 + 72 * n - 20) / 40
    e30 = QQ(13 * n**3 - 13 * n**2 + 62 * n) / 20
    assert n * (3 * n - 1) * b <= n**6
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
