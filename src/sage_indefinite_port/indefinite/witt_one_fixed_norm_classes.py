r"""``O(L)``-classes of fixed-norm and isotropic vectors of the Witt-index-one Hilbert model (#40).

``L = U + (O, Tr)`` is :class:`HilbertWittOneLattice` for a real quadratic field ``E = Q(sqrt D)`` with
maximal order ``O``; ``O`` need not be norm-Euclidean or principal. A vector is a Hermitian matrix
``X = [[a, beta], [conj beta, d]]`` with ``b(X, X) = -2 det X``. Write ``delta = det X``.

**Assumptions.** ``O`` is the maximal order of ``E`` and ``L`` is the model itself (or a lattice that
:func:`recognized_model` identifies with it). Only isometries of the following kinds are used, and each is an
element of ``O(L)``: ``-1``; ``X -> g X conj(g)^t`` for ``g in GL_2(O)`` with ``N(det g) = 1``; and the
generators of :meth:`HilbertWittOneLattice.generators`, which generate ``<SL_2(O), -1, Stab(e)>`` by
Vaserstein's theorem (see :mod:`witt_one_hilbert`). With class number above one they need not generate
``O(L)``: every orbit of the group they generate lies in one ``O(L)``-orbit, so a finite set meeting every such
orbit meets every ``O(L)``-orbit, and equal ``O(L)``-orbits are recognized by the exact decision below.

**Cusps.** For a nonzero ``w in O^2`` put ``h_X(w) = w X conj(w)^t``; for ``g in SL_2(O)`` the entry
``(g X conj(g)^t)_{11}`` is ``h_X(first row of g)``. Fix one integral ideal ``c`` in each ideal class, with
``c = O`` for the principal class, and generators ``c = c_1 O + c_2 O``. The rows ``u in c^2`` with
``u_1 O + u_2 O = c`` form one orbit under right multiplication by ``SL_2(O)``: complete ``u`` to
``M_u = [[u_1, u_2], [s_1, s_2]]`` with ``s_i in c^{-1}`` and ``det M_u = 1`` (solvable because
``c c^{-1} = O``); then ``M_c^{-1} M_u in SL_2(O)`` carries the row ``c`` to ``u``. These are the cusps of the
Hilbert modular group, one per ideal class (van der Geer, *Hilbert Modular Surfaces*, Springer 1988,
Prop. I.1.1). The frame ``Y = M_c X conj(M_c)^t`` has ``Y_{11} = h_X(c) in N(c) Z``,
``Y_{12} in c conj(c)^{-1}`` and ``Y_{22} in N(c)^{-1} Z``.

**Reduction bound.** ``h_X`` is an integral quadratic form on ``O^2 = Z^4``. With ``M = sigma_1(X)``,
``h_X(w) = sigma_1(w) M sigma_2(w)^t``, whose Gram matrix ``[[0, M/2], [M^t/2, 0]]`` has determinant
``delta^2 / 16`` in the embedding coordinates; ``O^2`` has covolume ``Disc(E)`` there, so the Gram matrix of
``h_X`` on ``Z^4`` has ``|det| = delta^2 Disc^2 / 16``. A positive majorant ``F >= |h_X|`` of the same
``|det|`` exists (diagonalize the Gram matrix orthogonally and take absolute eigenvalues). Hermite's
inequality with ``gamma_4 = sqrt 2`` (Korkine and Zolotarev, Math. Ann. 6 (1873); Cassels, *Rational
Quadratic Forms*, Ch. 12) gives a nonzero ``w`` with ``|h_X(w)| <= C``, ``C^2 = |delta| Disc / 2``. This uses
no Euclidean algorithm. Let ``a = w_1 O + w_2 O = lambda c`` with ``c`` the representative of its class. Then
``u = w / lambda`` generates ``c``, ``|N lambda| = N(a) / N(c)`` and ``h_X(u) = h_X(w) / N lambda``; moving
``u`` to the row ``c`` gives ``Y_{11} = N(c) a'`` with ``|a'| = |h_X(w)| / N(a) <= C``. So every orbit meets a
frame with ``|Y_{11}| <= N(c) C`` at some cusp. This is the Siegel-set reduction of ``SL_2(O)`` on
``H^2`` restricted to one coordinate (Strömberg, *A reduction algorithm for Hilbert modular groups*,
arXiv:2111.13274, Sec. 3, uses the same cusp decomposition).

**Normal forms at a cusp ``c``.** The row stabilizer of ``c`` in ``SL_2(O)`` is ``M_c^{-1} [[1, 0], [t, 1]]
M_c``, ``t in c^{-2}``; it sends ``Y_{12}`` to ``Y_{12} + Y_{11} conj(t)``. After ``-1``, ``Y_{11} >= 0``. If
``Y_{11} > 0``, ``Y_{12}`` is reduced to coordinates in ``[0, 1)`` with respect to ``Y_{11}`` times a basis of
``conj(c)^{-2}``, and ``Y_{22}`` is fixed by ``det``. If ``Y_{11} = 0``, then ``N(Y_{12}) = -delta``;
``M_c^{-1} diag(eta, 1) M_c``, ``eta in O^*`` of norm one, moves ``Y_{12}`` into the domain
``sigma_1 > 0``, ``|delta| <= sigma_1^2 < |delta| sigma_1(epsilon)^2`` for the norm-one generator ``epsilon > 1``;
then the shears send ``Y_{22}`` to ``Y_{22} + Tr(t Y_{12})``, so ``Y_{22}`` is reduced modulo the generator of
``Tr(c^{-2} Y_{12})``. The finite set of normal forms meets every orbit. For ``c = O`` these are the normal
forms ``0 <= x, y < a`` and ``0 <= d < gcd(Tr beta, Tr(beta omega))`` of the principal cusp.

**Isotropic vectors.** A primitive null ``X`` is ``mu u conj(u)^t`` with ``u`` a nonzero column of ``X``.
Moving ``u`` to the column ``c^t`` of its class carries ``X`` to ``+-`` the primitive vector ``e_c`` on the line
of ``c^t conj(c)``. So the ``e_c`` meet every orbit of primitive null vectors; ``-1`` identifies lines with
vectors.

**Exact decision.** Two normal forms are joined when a group element maps one to the other: the images of
the generators and of every row ``w`` in a box with ``|h(w)| <= C``. Joined vectors are equivalent, with a
verified transporter. The content and the divisor ``div(v)`` separate classes. Components with equal
invariants, and the cusp vectors ``e_c``, are decided by :func:`marked_transporter` on the coslice
``<b(x, x)> -> L``: for ``b(x, x) != 0`` this is ``EichlerRecursion.pointed_isometry``
(IndefiniteForms.g, INDEF_FORM_EquivalenceVector) on the rank-three complements ``x^perp``, which uses
``O(x^perp)`` and never ``O(L)``; for null ``x`` it is the isotropic pointed isometry through the binary
anisotropic reduction ``x^perp / x``. Every answer is proved: a transporter is verified, and ``None`` is a
negative decision.
"""

from math import ceil, gcd, lcm, sqrt

from sage.arith.misc import xgcd
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.qqbar import AA
from sage.rings.rational_field import QQ
from sage_extensions.categories.lattices import Lattices

from sage_indefinite_port.indefinite.pairing_slices import marked_transporter

from sage_indefinite_port.indefinite.witt_one_hilbert import recognized_model


def _integral_solution(system, target):
    r"""Return ``x in Z^n`` with ``system * x = target``, by the Smith form ``D = U system V``."""
    diagonal, left, right = system.smith_form()
    image = left * target
    solution = [ZZ(0)] * system.ncols()
    for i in range(system.nrows()):
        pivot = diagonal[i, i] if i < system.ncols() else 0
        if pivot == 0:
            if image[i] != 0:
                raise ArithmeticError(f"{system} x = {target} has no rational solution")
            continue
        if image[i] % pivot:
            raise ArithmeticError(f"{system} x = {target} has no integral solution")
        solution[i] = image[i] // pivot
    return right * vector(ZZ, solution)


def _complete_row(model, row, inverse_ideal):
    r"""Return ``(s_1, s_2) in (c^{-1})^2`` with ``u_1 s_2 - u_2 s_1 = 1`` for a row ``u`` generating ``c``."""
    E = model.field
    basis = tuple(E(b) for b in inverse_ideal.basis())
    u1, u2 = row
    columns = [-u2 * b for b in basis] + [u1 * b for b in basis]
    system = matrix(ZZ, [[ZZ(c) for c in model._beta_coordinates(column)] for column in columns]).transpose()
    x = _integral_solution(system, vector(ZZ, [1, 0]))
    first = x[0] * basis[0] + x[1] * basis[1]
    second = x[2] * basis[0] + x[3] * basis[1]
    if u1 * second - u2 * first != 1:
        raise ArithmeticError(f"the completion of {row} does not have determinant one")
    return first, second


class _Cusp:
    r"""The cusp of ``SL_2(O)`` attached to an integral ideal ``c``, with its frame ``M_c`` and lattices."""

    def __init__(self, model, ideal):
        E = model.field
        conj = model.conjugate
        self.ideal = ideal
        self.norm = ZZ(ideal.norm())
        principal = ideal == E.ideal(1)
        if principal:
            self.frame = identity_matrix(E, 2)
        else:
            row = tuple(E(g) for g in ideal.gens_two())
            second = _complete_row(model, row, ~ideal)
            self.frame = matrix(E, 2, 2, [row[0], row[1], second[0], second[1]])
        self.inverse = self.frame.inverse()
        self._conjugate = conj
        conjugate_ideal = E.ideal([conj(g) for g in ideal.gens()])
        if principal:
            self.shear_basis = model._beta_basis
            self.trace_basis = model._beta_basis
            self.pairing_basis = model._beta_basis
        else:
            self.shear_basis = tuple(E(b) for b in (~(conjugate_ideal * conjugate_ideal)).basis())
            self.trace_basis = tuple(E(b) for b in (~(ideal * ideal)).basis())
            self.pairing_basis = tuple(E(b) for b in (ideal * ~conjugate_ideal).basis())
        self._shear_coordinates = matrix(QQ, [list(b.vector()) for b in self.shear_basis]).inverse()

    def _star(self, g):
        return g.apply_map(self._conjugate).transpose()

    def to_frame(self, X):
        return self.frame * X * self._star(self.frame)

    def from_frame(self, Y):
        return self.inverse * Y * self._star(self.inverse)

    def move(self, g):
        r"""Return ``M_c^{-1} g M_c``: the action ``g`` on frames, as an action on ``L``."""
        return self.inverse * g * self.frame

    def shear_coordinates(self, z):
        return vector(QQ, list(z.vector())) * self._shear_coordinates

    def residues(self, A):
        r"""Yield ``B in c conj(c)^{-1}`` with coordinates of ``B / A`` in ``[0, 1)^2`` over ``conj(c)^{-2}``."""
        rows = matrix(QQ, [list(self.shear_coordinates(q / A)) for q in self.pairing_basis])
        inverse = rows.inverse()
        if inverse.denominator() != 1:
            raise ArithmeticError(f"{A} conj(c)^-2 is not contained in c conj(c)^-1")
        denominator = rows.denominator()
        for i in range(denominator):
            for k in range(denominator):
                point = vector(QQ, [QQ(i) / denominator, QQ(k) / denominator])
                if (point * inverse).denominator() == 1:
                    yield A * (point[0] * self.shear_basis[0] + point[1] * self.shear_basis[1])

    def trace_period(self, B):
        r"""Return ``(g, t)``: ``g > 0`` generates ``Tr(c^{-2} B)`` and ``t in c^{-2}`` has ``Tr(t B) = g``."""
        traces = [QQ((t * B).trace()) for t in self.trace_basis]
        denominator = lcm(*(int(t.denominator()) for t in traces))
        g, s1, s2 = xgcd(ZZ(traces[0] * denominator), ZZ(traces[1] * denominator))
        if g == 0:
            raise ArithmeticError(f"Tr(c^-2 {B}) is zero")
        return QQ(g) / denominator, s1 * self.trace_basis[0] + s2 * self.trace_basis[1]


class HilbertFixedNormClasses:
    r"""Classes of ``O(L)`` on ``{v in L : b(v, v) = square}`` in model coordinates."""

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
        self._lattice = Lattices(ZZ)([[int(c) for c in row] for row in self._gram.rows()])
        self._class_group = E.class_group()
        ideals = [E.ideal(1)]
        for element in self._class_group:
            if element.is_one():
                continue
            ideal = element.ideal()
            if not ideal.is_integral():
                raise ArithmeticError(f"the class representative {ideal} is not integral")
            ideals.append(ideal)
        self._cusps = tuple(_Cusp(model, ideal) for ideal in ideals)
        self._components = {}
        self._isotropic = None

    # -- arithmetic of the model --

    def _delta(self, square):
        if square == 0:
            raise ValueError("norm zero: the null vectors are classified by the cusps")
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

    def _class_of(self, ideal):
        r"""Return ``(j, lambda)`` with ``ideal = lambda c_j``."""
        target = self._class_group(ideal)
        for j, cusp in enumerate(self._cusps):
            if self._class_group(cusp.ideal) == target:
                generators = (ideal / cusp.ideal).gens_reduced()
                if len(generators) != 1:
                    raise ArithmeticError(f"{ideal} / {cusp.ideal} is not principal")
                return j, generators[0]
        raise ArithmeticError(f"the class of {ideal} has no representative")

    def _row_move(self, j, row):
        r"""Return ``g in SL_2(O)`` with ``c_j g = row``, for a row generating ``c_j``."""
        cusp = self._cusps[j]
        E = self._model.field
        second = _complete_row(self._model, row, ~cusp.ideal)
        g = cusp.inverse * matrix(E, 2, 2, [row[0], row[1], second[0], second[1]])
        if g.det() != 1 or not all(entry.is_integral() for entry in g.list()):
            raise ArithmeticError(f"the move of the cusp {cusp.ideal} to {row} is not in SL_2(O)")
        if tuple(cusp.frame.row(0) * g) != tuple(row):
            raise ArithmeticError(f"the move of the cusp {cusp.ideal} does not reach {row}")
        return g

    def _in_unit_domain(self, beta, delta):
        first = AA(beta)
        return first > 0 and abs(delta) <= first**2 < abs(delta) * AA(self._epsilon) ** 2

    def _normalize(self, v, j):
        r"""Return ``(m, c)`` with ``m in O(L)``, ``m v = c`` and ``c`` the normal form at the cusp ``c_j``."""
        model, cusp = self._model, self._cusps[j]
        E = model.field
        m = identity_matrix(ZZ, 4)
        v = vector(ZZ, v)
        Y = cusp.to_frame(model._matrix_of(v))
        if QQ(Y[0, 0]) < 0:
            m, v, Y = -m, -v, -Y
        A = QQ(Y[0, 0])
        if A > 0:
            x, y = cusp.shear_coordinates(Y[0, 1] / A)
            t_bar = -(x.floor() * cusp.shear_basis[0] + y.floor() * cusp.shear_basis[1])
            shear = model.action(cusp.move(matrix(E, 2, 2, [1, 0, model.conjugate(t_bar), 1])))
            return shear * m, shear * v
        beta = Y[0, 1]
        delta = -beta.norm()
        ceiling = abs(delta) * AA(self._epsilon) ** 2
        eta = E.one()
        while AA(beta) ** 2 >= ceiling:
            beta, eta = beta / self._epsilon, eta / self._epsilon
        while AA(beta) ** 2 < abs(delta):
            beta, eta = beta * self._epsilon, eta * self._epsilon
        if AA(beta) < 0:
            eta = -eta
        rotation = model.action(cusp.move(matrix(E, 2, 2, [eta, 0, 0, 1])))
        m, v = rotation * m, rotation * v
        Y = cusp.to_frame(model._matrix_of(v))
        beta = Y[0, 1]
        if not self._in_unit_domain(beta, delta):
            raise ArithmeticError(f"the unit normalization of {beta} left the domain")
        period, unit_shift = cusp.trace_period(beta)
        k = -(QQ(Y[1, 1]) / period).floor()
        shear = model.action(cusp.move(matrix(E, 2, 2, [1, 0, k * unit_shift, 1])))
        return shear * m, shear * v

    def _move(self, v, w1, w2):
        r"""Return ``(m, c)``: ``m`` moves the row ``(w1, w2) / lambda`` to its cusp, then normalizes."""
        j, scalar = self._class_of(self._model.field.ideal(w1, w2))
        first = self._model.action(self._row_move(j, (w1 / scalar, w2 / scalar)))
        m, c = self._normalize(first * v, j)
        return m * first, c

    def _reduce(self, v):
        r"""Return ``(m, c)`` with ``m v = c`` a normal form; the row search ends by Hermite's bound."""
        delta = self._delta(vector(ZZ, v) * self._gram * vector(ZZ, v))
        radius = 1
        while True:
            for w1, w2 in self._rows(radius):
                if self._within_bound(self._h(v, w1, w2), delta):
                    return self._move(v, w1, w2)
            radius += 1

    # -- invariants --

    def _invariant(self, v):
        r"""Return ``(content, div)``, both preserved by ``O(L)``."""
        v = vector(ZZ, v)
        return gcd(*(int(c) for c in v)), gcd(*(int(c) for c in self._gram * v))

    def _decided_transporter(self, x, y):
        r"""Return the matrix of ``g in O(L)`` with ``g x = y``, or ``None`` when none exists (exact)."""
        lattice = self._lattice
        frame = tuple(lattice.module_generators())

        def element(v):
            return sum((lattice.scalar_multiple(int(c), g) for c, g in zip(v, frame, strict=True)), lattice.zero())

        witness = marked_transporter(lattice, (element(x),), (element(y),))
        if witness is None:
            return None
        return matrix(ZZ, [[int(c) for c in witness(g)] for g in frame]).transpose()

    # -- fixed nonzero norm --

    def _keep(self, found, X):
        coordinates = self._model._vector_of(X)
        if all(c in ZZ for c in coordinates):
            v = vector(ZZ, coordinates)
            v.set_immutable()
            found.setdefault(v, None)

    def _unit_domain_elements(self, cusp, delta):
        r"""Yield ``B in c conj(c)^{-1}`` with ``N(B) = -delta`` in the unit domain."""
        conj = self._model.conjugate
        q1, q2 = cusp.pairing_basis
        a, b, c, d = (float(AA(q1)), float(AA(q2)), float(AA(conj(q1))), float(AA(conj(q2))))
        determinant = a * d - b * c
        radius = sqrt(abs(float(delta))) * float(AA(self._epsilon))
        x_bound = ceil(radius * (abs(d) + abs(b)) / abs(determinant)) + 1
        y_bound = ceil(radius * (abs(c) + abs(a)) / abs(determinant)) + 1
        for x in range(-x_bound, x_bound + 1):
            for y in range(-y_bound, y_bound + 1):
                B = x * q1 + y * q2
                if B.norm() == -delta and self._in_unit_domain(B, delta):
                    yield B

    def candidates(self, square: int):
        r"""Return the normal forms with ``b(v, v) = square``: a finite set meeting every ``O(L)``-orbit."""
        if square % 2:
            return ()
        delta = self._delta(square)
        model = self._model
        E = model.field
        conj = model.conjugate
        found = {}
        for cusp in self._cusps:
            a = 1
            while self._within_bound(a, delta):
                A = cusp.norm * a
                for B in cusp.residues(A):
                    self._keep(found, cusp.from_frame(matrix(E, 2, 2, [A, B, conj(B), (delta + B.norm()) / A])))
                a += 1
            for B in self._unit_domain_elements(cusp, delta):
                period, _shift = cusp.trace_period(B)
                k = 0
                while QQ(k) / cusp.norm < period:
                    self._keep(found, cusp.from_frame(matrix(E, 2, 2, [0, B, conj(B), QQ(k) / cusp.norm])))
                    k += 1
        return tuple(found)

    def _neighbours(self, c, delta):
        for g in self._generators:
            m, image = self._reduce(g * c)
            yield m * g, image
        for w1, w2 in self._rows(1):
            if self._within_bound(self._h(c, w1, w2), delta):
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
        kept = {}
        for root in sorted(set(root_of.values())):
            invariant = self._invariant(root)
            link, target = None, None
            for known in kept.get(invariant, ()):
                link = self._decided_transporter(root, known)
                if link is not None:
                    target = known
                    break
            if link is None:
                kept.setdefault(invariant, []).append(root)
                continue
            for c in [c for c, r in root_of.items() if r == root]:
                root_of[c] = target
                to_root[c] = (link * to_root[c]).change_ring(ZZ)
        self._components[square] = (root_of, to_root)
        return root_of, to_root

    # -- norm zero --

    def _cusp_vector(self, j):
        r"""Return the primitive vector ``e_c`` on the null line of ``c^t conj(c)``."""
        model = self._model
        column = matrix(model.field, 2, 1, list(self._cusps[j].frame.row(0)))
        coordinates = vector(ZZ, model._vector_of(column * column.apply_map(model.conjugate).transpose()))
        e = coordinates / gcd(*(int(c) for c in coordinates))
        e = vector(ZZ, e)
        e.set_immutable()
        return e

    def _isotropic_move(self, v):
        r"""Return ``(m, j)`` with ``m in O(L)`` and ``m v = e_{c_j}``, for a primitive null ``v``."""
        model = self._model
        X = model._matrix_of(v)
        column = (X[0, 0], X[1, 0]) if X[0, 0] != 0 else (X[0, 1], X[1, 1])
        j, scalar = self._class_of(model.field.ideal(list(column)))
        row_move = self._row_move(j, (column[0] / scalar, column[1] / scalar))
        moved = model.action(row_move.inverse().transpose())
        image, target = moved * vector(ZZ, v), self._cusp_vector(j)
        if image == -target:
            moved = -moved
        elif image != target:
            raise ArithmeticError(f"the null vector {v} did not reach the cusp vector {target}")
        return moved, j

    def _isotropic_data(self):
        r"""Return ``(root_of, to_root)`` on cusp indices: ``to_root[j] e_j = e_{root_of[j]}``."""
        if self._isotropic is not None:
            return self._isotropic
        roots, root_of, to_root = [], {}, {}
        for j in range(len(self._cusps)):
            e = self._cusp_vector(j)
            for r in roots:
                known = self._cusp_vector(r)
                if self._invariant(e) != self._invariant(known):
                    continue
                link = self._decided_transporter(e, known)
                if link is not None:
                    root_of[j], to_root[j] = r, link
                    break
            else:
                roots.append(j)
                root_of[j], to_root[j] = j, identity_matrix(ZZ, 4)
        self._isotropic = (root_of, to_root)
        return self._isotropic

    def cusp_links(self):
        r"""Return, for each cusp vector ``e_c != e`` in the ``O(L)``-orbit of ``e = e_O``, an isometry ``e_c -> e``.

        With the model generators these generate ``O(L)`` (:mod:`witt_one_hilbert`).
        """
        root_of, to_root = self._isotropic_data()
        return tuple(to_root[j] for j in sorted(root_of) if j != 0 and root_of[j] == 0)

    def _isotropic_transporter(self, x, y):
        content = gcd(*(int(c) for c in x))
        if content == 0 or content != gcd(*(int(c) for c in y)):
            return None
        mx, jx = self._isotropic_move(vector(ZZ, x / content))
        my, jy = self._isotropic_move(vector(ZZ, y / content))
        root_of, to_root = self._isotropic_data()
        if root_of[jx] != root_of[jy]:
            return None
        return my.inverse() * to_root[jy].inverse() * to_root[jx] * mx

    # -- public interface --

    def transporter(self, x, y):
        r"""Return the matrix of an element of ``O(L)`` with ``x -> y``, ``None`` if none exists."""
        x, y = vector(ZZ, x), vector(ZZ, y)
        square = x * self._gram * x
        if y * self._gram * y != square or self._invariant(x) != self._invariant(y):
            return None
        if square == 0:
            h = self._isotropic_transporter(x, y)
        else:
            root_of, to_root = self._component_data(square)
            mx, cx = self._reduce(x)
            my, cy = self._reduce(y)
            cx.set_immutable()
            cy.set_immutable()
            if root_of[cx] != root_of[cy]:
                return None
            h = my.inverse() * to_root[cy].inverse() * to_root[cx] * mx
        if h is None:
            return None
        h = h.change_ring(ZZ)
        if h * x != y or h.transpose() * self._gram * h != self._gram:
            raise ArithmeticError("the constructed transporter is not an isometry carrying x to y")
        return h

    def classes(self, square: int, primitive: bool):
        r"""Return one vector per ``O(L)``-class with ``b(v, v) = square`` (primitive ones if asked).

        For ``square = 0`` only the primitive classes form a finite set: the multiples ``k e``, ``k >= 2``, give
        one class for every ``k``.
        """
        if square % 2:
            return ()
        if square == 0:
            if not primitive:
                raise ValueError("the null vectors k e, k >= 2, form infinitely many classes; ask for primitive ones")
            root_of, _to_root = self._isotropic_data()
            return tuple(self._cusp_vector(j) for j in sorted(set(root_of.values())))
        root_of, _to_root = self._component_data(square)
        roots = sorted(set(root_of.values()))
        if primitive:
            roots = [r for r in roots if gcd(*(int(c) for c in r)) == 1]
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
