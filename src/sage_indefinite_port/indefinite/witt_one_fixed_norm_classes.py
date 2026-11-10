r"""``O(L)``-classes of fixed-norm vectors of the Witt-index-one Hilbert model (#40).

``L = U + (O, Tr)`` is :class:`HilbertWittOneLattice` for ``E = Q(sqrt D)`` with norm-Euclidean ``O``. A vector
is a Hermitian matrix ``X = [[a, beta], [conj beta, d]]`` with ``b(X, X) = -2 det X``. Write ``delta = det X``.

**Assumptions.** ``O`` is norm-Euclidean, as :meth:`HilbertWittOneLattice._reduce_pair` requires; hence ``O``
is a principal ideal domain and every coprime row of ``O^2`` is the first row of an element of ``SL_2(O)``.
The group is ``O(L) = {+-1} x G ⋊ <c>``: ``G`` is the image of ``{g in GL_2(O) : det g in U^1}``
(``U^1`` the norm-one units) acting by ``X -> g X conj(g)^t``, and ``c`` is Galois conjugation of ``beta``.
This holds because ``O(L) = <SL_2(O), -1, Stab(e)>`` (module :mod:`witt_one_hilbert`) and the unit rotation
``beta -> eta beta``, ``eta in U^1``, is the action of ``diag(eta, 1)``. Only the generators of
:meth:`HilbertWittOneLattice.generators` and matrices built from ``SL_2(O)`` and ``diag(eta, 1)`` are used.

**Reduction bound.** For a coprime row ``w`` the entry ``(g X conj(g)^t)_{11}`` of a ``g`` with first row ``w``
is ``h_X(w) = w X conj(w)^t``, an integral quadratic form on ``O^2 = Z^4``. With ``M = sigma_1(X)``,
``h_X(w) = sigma_1(w) M sigma_2(w)^t``, whose Gram matrix ``[[0, M/2], [M^t/2, 0]]`` has determinant
``delta^2 / 16`` in the embedding coordinates; ``O^2`` has covolume ``Disc(E)`` there, so the Gram matrix of
``h_X`` on ``Z^4`` has ``|det| = delta^2 Disc^2 / 16``. A positive majorant ``F >= |h_X|`` of the same
``|det|`` exists (diagonalize the Gram matrix orthogonally and take absolute eigenvalues). Hermite's
inequality with ``gamma_4 = sqrt 2`` (Korkine and Zolotarev, Math. Ann. 6 (1873); Cassels, *Rational
Quadratic Forms*, Ch. 12) gives a nonzero ``w`` with ``|h_X(w)| <= F(w) <= C``, ``C^2 = |delta| Disc / 2``.
Dividing ``w`` by a generator ``lambda`` of the ideal ``(w_1, w_2)`` divides ``h_X`` by ``|N lambda| >= 1``.
So every orbit meets ``|a| <= C``. This replaces the Siegel-set reduction of ``SL_2(O)`` on ``H^2``: the cusp
and the compact core do not enter, because the bound concerns one value of ``h_X``.

**Normal forms.** After ``-1``, ``a >= 0``. If ``a > 0``, the lower shears ``[[1, 0], [t, 1]]`` send ``beta`` to
``beta + a conj(t)``, so ``beta`` is reduced to coordinates in ``[0, a)``. If ``a = 0``, then ``N beta = -delta``;
``diag(eta, 1)`` moves ``beta`` into the domain ``sigma_1(beta) > 0``, ``|delta| <= sigma_1(beta)^2 <
|delta| sigma_1(epsilon)^2`` for the norm-one generator ``epsilon > 1``; then the lower shears send ``d`` to
``d + Tr(t beta)``, so ``d`` is reduced modulo ``gcd(Tr beta, Tr(beta omega))``. The finite set of normal forms
meets every orbit.

**Exact decision.** Two normal forms are joined when a group element maps one to the other: the images of
the generators and of every coprime row ``w`` in a box with ``|h(w)| <= C``. Joined vectors are equivalent,
with a verified transporter. Different classes are proved by invariants preserved by ``O(L)``: the content,
the divisor ``div(v)`` and the ``O(L)``-orbit of ``v / div(v)`` in ``L^v / L``. When two components carry equal
invariants the procedure has no proof either way and raises :class:`ArithmeticError`; it never returns an
unproved answer.
"""

from math import ceil, gcd, sqrt

from sage.arith.misc import xgcd
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.qqbar import AA
from sage.rings.rational_field import QQ

from sage_indefinite_port.indefinite.witt_one_hilbert import recognized_model


class HilbertFixedNormClasses:
    r"""Classes of ``O(L)`` on ``{v in L : b(v, v) = square}``, ``square != 0``, in model coordinates."""

    def __init__(self, model) -> None:
        self._model = model
        E = model.field
        self._disc = ZZ(E.discriminant())
        epsilon = model.norm_one_unit()
        if AA(epsilon) < 0:
            epsilon = -epsilon
        if AA(epsilon) < 1:
            epsilon = ~epsilon
        self._epsilon = epsilon
        self._generators = model.generators()
        self._gram = model.gram()
        self._components = {}

    # -- arithmetic of the model --

    def _delta(self, square):
        if square == 0:
            raise ValueError("norm zero: use isotropic_line_representative and transporter_to_e")
        return -ZZ(square) / 2

    def _within_bound(self, value, delta):
        return 2 * QQ(value) ** 2 <= abs(delta) * self._disc

    def _h(self, v, w1, w2):
        X = self._model._matrix_of(v)
        conj = self._model.conjugate
        value = (
            w1 * conj(w1) * X[0, 0] + w1 * X[0, 1] * conj(w2) + w2 * X[1, 0] * conj(w1) + w2 * conj(w2) * X[1, 1]
        )
        return QQ(value)

    def _rows(self, radius):
        beta = self._model._beta
        for c in _box(radius):
            if any(c):
                yield beta(c[0], c[1]), beta(c[2], c[3])

    def _first_row(self, w1, w2):
        r"""Return ``g in SL_2(O)`` with first row ``(w1, w2)``, for a coprime row."""
        E = self._model.field
        reduction, end = self._model._reduce_pair(vector(E, [w1, w2]))
        nu = end[0]
        g = matrix(E, 2, 2, [nu, 0, 0, ~nu]) * (~reduction).transpose()
        if g[0, 0] != w1 or g[0, 1] != w2 or g.det() != 1:
            raise ArithmeticError(f"the completion of ({w1}, {w2}) is not in SL_2(O)")
        return g

    def _in_unit_domain(self, beta, delta):
        first = AA(beta)
        return first > 0 and abs(delta) <= first**2 < abs(delta) * AA(self._epsilon) ** 2

    def _normalize(self, v):
        r"""Return ``(m, c)`` with ``m in O(L)``, ``m v = c`` and ``c`` a normal form."""
        model = self._model
        E = model.field
        m = identity_matrix(ZZ, 4)
        v = vector(ZZ, v)
        if v[0] < 0:
            m, v = -m, -v
        a = v[0]
        if a > 0:
            x, y = v[2], v[3]
            t_bar = model._beta(-(x // a), -(y // a))
            shear = model.action(matrix(E, 2, 2, [1, 0, model.conjugate(t_bar), 1]))
            return shear * m, shear * v
        beta = model._beta(v[2], v[3])
        delta = -beta.norm()
        ceiling = abs(delta) * AA(self._epsilon) ** 2
        eta = E.one()
        while AA(beta) ** 2 >= ceiling:
            beta, eta = beta / self._epsilon, eta / self._epsilon
        while AA(beta) ** 2 < abs(delta):
            beta, eta = beta * self._epsilon, eta * self._epsilon
        if AA(beta) < 0:
            eta = -eta
        rotation = model._unit_rotation(eta)
        m, v = rotation * m, rotation * v
        beta = model._beta(v[2], v[3])
        if not self._in_unit_domain(beta, delta):
            raise ArithmeticError(f"the unit normalization of {beta} left the domain")
        g, s1, s2 = xgcd(ZZ(beta.trace()), ZZ((beta * model.omega).trace()))
        k = -(v[1] // g)
        shear = model.action(matrix(E, 2, 2, [1, 0, k * model._beta(s1, s2), 1]))
        return shear * m, shear * v

    def _move(self, v, w1, w2):
        r"""Return ``(m, c)``: ``m`` puts the coprime row ``(w1, w2)`` first, then normalizes."""
        first = self._model.action(self._first_row(w1, w2))
        m, c = self._normalize(first * v)
        return m * first, c

    def _reduce(self, v):
        r"""Return ``(m, c)`` with ``m v = c`` a normal form; the row search ends by Hermite's bound."""
        delta = self._delta(vector(ZZ, v) * self._gram * vector(ZZ, v))
        radius = 1
        while True:
            for w1, w2 in self._rows(radius):
                if self._within_bound(self._h(v, w1, w2), delta):
                    generator = self._model.field.ideal(w1, w2).gens_reduced()[0]
                    return self._move(v, w1 / generator, w2 / generator)
            radius += 1

    # -- invariants --

    def _invariant(self, v):
        r"""Return ``(content, div, orbit of v/div in L^v/L)``, each preserved by ``O(L)``."""
        v = vector(ZZ, v)
        content = gcd(*(int(c) for c in v))
        divisor = gcd(*(int(c) for c in self._gram * v))
        key = lambda u: tuple(c - c.floor() for c in u)
        start = vector(QQ, v) / divisor
        orbit = {key(start)}
        frontier = [start]
        while frontier:
            following = []
            for u in frontier:
                for g in self._generators:
                    image = g * u
                    if key(image) not in orbit:
                        orbit.add(key(image))
                        following.append(vector(QQ, key(image)))
            frontier = following
        return content, divisor, frozenset(orbit)

    # -- public interface --

    def candidates(self, square: int):
        r"""Return the normal forms with ``b(v, v) = square``: a finite set meeting every ``O(L)``-orbit."""
        if square % 2:
            return ()
        delta = self._delta(square)
        model = self._model
        found = []
        a = 1
        while self._within_bound(a, delta):
            for x in range(a):
                for y in range(a):
                    numerator = delta + model._beta(x, y).norm()
                    if numerator % a == 0:
                        found.append(vector(ZZ, [a, numerator / a, x, y]))
            a += 1
        root, unit = sqrt(abs(float(delta))), float(AA(self._epsilon))
        omega = (float(AA(model.omega)), float(AA(model.conjugate(model.omega))))
        y_bound = ceil(root * (unit + 1) / abs(omega[0] - omega[1])) + 1
        x_bound = ceil(root * unit + y_bound * max(abs(omega[0]), abs(omega[1]))) + 1
        for x in range(-x_bound, x_bound + 1):
            for y in range(-y_bound, y_bound + 1):
                beta = model._beta(x, y)
                if beta.norm() != -delta or not self._in_unit_domain(beta, delta):
                    continue
                period = gcd(int(beta.trace()), int((beta * model.omega).trace()))
                found.extend(vector(ZZ, [0, d, x, y]) for d in range(period))
        for v in found:
            v.set_immutable()
        return tuple(found)

    def _neighbours(self, c, delta):
        for g in self._generators:
            m, image = self._reduce(g * c)
            yield m * g, image
        for w1, w2 in self._rows(1):
            if self._model.field.ideal(w1, w2).is_one() and self._within_bound(self._h(c, w1, w2), delta):
                yield self._move(c, w1, w2)

    def _component_data(self, square):
        r"""Return ``(root_of, to_root)`` with ``to_root[c] * c = root_of[c]`` for the candidates."""
        if square in self._components:
            return self._components[square]
        delta = self._delta(square)
        vertices = set(self.candidates(square))
        root_of, to_root = {}, {}
        for start in self.candidates(square):
            if start in root_of:
                continue
            root_of[start], to_root[start] = start, identity_matrix(ZZ, 4)
            queue = [start]
            while queue:
                c = queue.pop()
                for m, image in self._neighbours(c, delta):
                    image.set_immutable()
                    if image not in vertices:
                        raise ArithmeticError(f"the normal form {image} is not a candidate")
                    if image not in root_of:
                        root_of[image] = start
                        to_root[image] = (to_root[c] * m.inverse()).change_ring(ZZ)
                        queue.append(image)
        self._components[square] = (root_of, to_root)
        return root_of, to_root

    def transporter(self, x, y):
        r"""Return the matrix of an element of ``O(L)`` with ``x -> y``, ``None`` if none exists.

        Raises :class:`ArithmeticError` when neither a transporter nor a separating invariant is found.
        """
        x, y = vector(ZZ, x), vector(ZZ, y)
        square = x * self._gram * x
        if y * self._gram * y != square or self._invariant(x) != self._invariant(y):
            return None
        root_of, to_root = self._component_data(square)
        mx, cx = self._reduce(x)
        my, cy = self._reduce(y)
        cx.set_immutable()
        cy.set_immutable()
        if root_of[cx] != root_of[cy]:
            raise ArithmeticError(f"{x} and {y} have equal invariants and lie in different components")
        h = (my.inverse() * to_root[cy].inverse() * to_root[cx] * mx).change_ring(ZZ)
        if h * x != y or h.transpose() * self._gram * h != self._gram:
            raise ArithmeticError("the constructed transporter is not an isometry carrying x to y")
        return h

    def classes(self, square: int, primitive: bool):
        r"""Return one vector per ``O(L)``-class with ``b(v, v) = square`` (primitive ones if asked)."""
        if square % 2:
            return ()
        root_of, _to_root = self._component_data(square)
        roots = sorted(set(root_of.values()))
        if primitive:
            roots = [r for r in roots if gcd(*(int(c) for c in r)) == 1]
        invariants = [self._invariant(r) for r in roots]
        if len(set(invariants)) != len(invariants):
            raise ArithmeticError(f"two components of norm {square} have equal invariants")
        return tuple(roots)


def _box(radius):
    span = range(-radius, radius + 1)
    return ((p, q, r, s) for p in span for q in span for r in span for s in span)


def witt_one_fixed_norm_vector_classes(lattice, square: int, primitive: bool):
    r"""Return one element of ``L`` per ``O(L)``-class with ``b(v, v) = square``, for a recognized Hilbert model."""
    recognized = recognized_model(lattice)
    if recognized is None:
        raise ValueError("the lattice was not recognized as a Hilbert model U + (O_E, Tr)")
    model, basis = recognized
    frame = tuple(lattice.module_generators())
    vectors = HilbertFixedNormClasses(model).classes(square, primitive)
    return tuple(
        sum((lattice.scalar_multiple(int(c), g) for c, g in zip(basis * v, frame, strict=True)), lattice.zero())
        for v in vectors
    )


__all__ = ["HilbertFixedNormClasses", "witt_one_fixed_norm_vector_classes"]
