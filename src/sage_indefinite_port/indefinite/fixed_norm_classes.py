r"""``O(L)``-classes of fixed-norm vectors of a rationally anisotropic quaternion norm lattice (#40).

``L = (Lambda, s Nrd)`` with ``Lambda = ZZ<1, i, j, k>`` in ``A = (a, b)_QQ``, ``a < 0 < b``, as built by
:class:`QuaternionNormLattice`; ``diag(2, 2, -6, -6)`` is ``(-1, 3)`` with ``s = 2``.

**Reduction bound.** Let ``S0`` be the perfect form of :class:`QuaternionVoronoi`, normalized to minimum one,
and ``phi(z) = <S0, Z Z^t>`` its value on ``z``. By Voronoi theory for ``Lambda^x`` (Koecher; Braun, Coulangeon,
Nebe and Schoennenbeck, J. Algebra 435 (2015), Sec. 2-3) the Voronoi cells
``V(u^t S0 u) = cone{u^{-1} x x^t u^{-t} : x in Min(S0)}``, ``u in Lambda^x``, cover the positive cone, since the
perfect forms form one ``Lambda^x``-class. Hence ``P = Z Z^t`` equals ``sum lambda_k r_k`` with ``lambda_k >= 0``
and rays ``r_k`` of one cell, and ``phi(u z) = <u^t S0 u, P> = sum lambda_k``. The function ``sqrt det`` is
concave and homogeneous on positive ``2 x 2`` matrices, so
``|Nrd z| = sqrt det P >= sum lambda_k sqrt det r_k = sum lambda_k |Nrd x_k| >= m sum lambda_k`` with
``m = min |Nrd x_k| >= 1``. So every ``Lambda^x``-orbit meets ``phi(z) <= |Nrd z| / m``. A left factor ``u`` of
norm ``-1`` is corrected by a fixed right factor ``e`` of norm ``-1``; ``z -> u z e`` lies in ``O(L)``.

**Exact equivalence.** A proper isometry is ``z -> a z a^{-1} w`` with ``a in N(Lambda)`` and ``w in Lambda^1``
(see :func:`orthogonal_group_generators`); improper ones precompose conjugation. Writing ``a = q nu m`` with
``nu`` a class of ``N(Lambda) / QQ^x Lambda^x`` and ``m in Lambda^x``, ``y ~ x`` iff for some ``nu`` and
``x' in {x, conj x}`` there are units ``t, w'' in Lambda^x`` with ``nu^{-1} y nu = t^{-1} x' w''``, i.e.
``x'^{-1} t y'`` is integral for ``y' = nu^{-1} y nu``. Since ``x'^{-1} n Lambda y' = conj(x') Lambda y'`` is
integral for ``n = Nrd x``, this depends only on ``t`` modulo ``n Lambda``. The image of ``Lambda^x`` in
``Lambda / n Lambda`` is the finite closure of the images of the Voronoi unit generators, so the test is finite.

The normalizer classes are supplied by the caller (``normalizer_classes_in_box``); their completeness is
the caller's.
"""

from math import gcd

from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ


def _coordinates(element):
    return tuple(QQ(c) for c in element)


def _integral(element):
    return all(c.denominator() == 1 for c in element)


def _unit_residues(model, voronoi, modulus):
    r"""Return ``{residue: unit}``: one unit of ``Lambda`` per class of the image of ``Lambda^x`` mod ``modulus``."""
    algebra = model.algebra()
    generators = [algebra(_coordinates(g)) for g in voronoi.unit_generators()]
    key = lambda unit: tuple(int(c) % modulus for c in unit)
    found = {key(algebra.one()): algebra.one()}
    frontier = [algebra.one()]
    while frontier:
        following = []
        for unit in frontier:
            for generator in generators:
                product = unit * generator
                residue = key(product)
                if residue not in found:
                    found[residue] = product
                    following.append(product)
        frontier = following
    return found


class FixedNormClasses:
    r"""Classes of ``O(L)`` on ``{z in Lambda : Nrd z = n}`` with exact transporters."""

    def __init__(self, model, voronoi, normalizer_classes) -> None:
        self._model = model
        self._voronoi = voronoi
        algebra = model.algebra()
        self._normalizers = (algebra.one(), *(algebra(_coordinates(nu)) for nu in normalizer_classes))
        self._form = voronoi.perfect_forms()[0]
        _, minimal = voronoi._minimal_vectors(self._form)
        self._m = min(abs(model.element(*c).reduced_norm()) for c in minimal)
        units = [algebra(_coordinates(u)) for u in voronoi.unit_generators()]
        self._flip = next(u for u in units if u.reduced_norm() == -1)

    def candidates(self, nrd: int):
        r"""Return a finite set of ``z`` with ``Nrd z = nrd`` meeting every ``O(L)``-orbit."""
        if nrd == 0:
            raise ValueError("the lattice is anisotropic: no nonzero vector has norm zero")
        bound = QQ(abs(nrd)) / self._m
        found = []
        for c in self._voronoi._short(self._form, bound):
            z = self._model.element(*c)
            if z.reduced_norm() == nrd:
                found.append(z)
            elif z.reduced_norm() == -nrd:
                found.append(z * self._flip)
        return found

    def transporter(self, x, y):
        r"""Return the matrix on ``1, i, j, k`` of an isometry carrying ``x`` to ``y``, or ``None`` if none exists."""
        n = x.reduced_norm()
        if y.reduced_norm() != n:
            return None
        model = self._model
        residues = _unit_residues(model, self._voronoi, int(abs(n)))
        for nu in self._normalizers:
            y_prime = ~nu * y * nu
            for conjugate in (False, True):
                source = x.conjugate() if conjugate else x
                for t in residues.values():
                    w = ~source * t * y_prime
                    if not _integral(w):
                        continue
                    left, right = nu * ~t, w * ~nu
                    columns = model.pair_matrix(left, right.conjugate())
                    if conjugate:
                        columns = columns * model.conjugation_matrix()
                    if columns * matrix(QQ, [list(x)]).transpose() != matrix(QQ, [list(y)]).transpose():
                        raise ArithmeticError("the constructed transporter does not carry x to y")
                    return columns
        return None

    def classes(self, nrd: int, primitive: bool):
        r"""Return one quaternion per ``O(L)``-class of ``z`` with ``Nrd z = nrd`` (primitive ones if asked)."""
        representatives = []
        for z in self.candidates(nrd):
            if primitive and gcd(*(int(c) for c in z)) != 1:
                continue
            if not any(self.transporter(known, z) is not None for known in representatives):
                representatives.append(z)
        return tuple(representatives)


def fixed_norm_vector_classes(model, voronoi, normalizer_classes, square: int, primitive: bool):
    r"""Return one element of ``L`` per ``O(L)``-class of vectors with ``b(v, v) = square``."""
    lattice = model.lattice()
    generators = tuple(lattice.module_generators())
    scale = int(lattice.b(generators[0], generators[0]))
    if square % scale:
        return ()
    quaternions = FixedNormClasses(model, voronoi, normalizer_classes).classes(square // scale, primitive)
    return tuple(
        sum((int(c) * g for c, g in zip(z, generators, strict=True)), lattice.zero()) for z in quaternions
    )


__all__ = ["FixedNormClasses", "fixed_norm_vector_classes"]
