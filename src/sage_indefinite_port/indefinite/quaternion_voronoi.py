r"""Voronoi generators of the unit group of an indefinite quaternion order.

For ``A = (a, b)_QQ`` with ``a < 0 < b`` and ``b`` not a square, and ``Lambda = ZZ<1,i,j,k>``,
``A`` embeds in ``M_2(K)`` for ``K = QQ(sqrt b)``. Positive ``S`` in ``Sym^2(K)`` give forms
``phi_S(x) = Tr(x^t S x)`` on ``Lambda``. Following Braun, Coulangeon, Nebe and Schoennenbeck,
J. Algebra 435 (2015), the perfect forms modulo ``Lambda^x`` are enumerated by flips across
the facets of the three-dimensional cone on the rays ``x x^t``. The stabilizers and the
neighbour transports generate ``Lambda^x`` (issue #57).
"""

from itertools import combinations
from math import acos, atan2, pi, sqrt

from sage.libs.pari import pari
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RR
from sage_extensions.algebras.norm_one_units import norm_one_generators


class QuaternionVoronoi:
    r"""Perfect forms and unit generators of ``ZZ<1, i, j, k>`` in ``(a, b)_QQ``."""

    def __init__(self, a: int, b: int) -> None:
        if a >= 0 or b <= 0 or QQ(b).is_square():
            raise ValueError("the model requires a < 0 < b with b not a square")
        self._field = QuadraticField(b, "s")
        root = self._field.gen()
        self._embedding = next(e for e in self._field.embeddings(RR) if e(root) > 0)
        self._algebra = QuaternionAlgebra(QQ, a, b)
        self._basis = (self._algebra.one(), *self._algebra.gens())
        i, j = matrix(self._field, [[0, a], [1, 0]]), matrix(self._field, [[root, 0], [0, -root]])
        self._matrices = (matrix.identity(self._field, 2), i, j, i * j)

    def _matrix(self, coordinates):
        return sum(int(c) * m for c, m in zip(coordinates, self._matrices, strict=True))

    def _element(self, coordinates):
        return sum(int(c) * q for c, q in zip(coordinates, self._basis, strict=True))

    def form_value(self, form, element):
        r"""Return ``phi_S`` of a quaternion with integer coordinates."""
        x = self._matrix([int(c) for c in element])
        return (x.transpose() * form * x).trace()

    def _value(self, form, x):
        return (x.transpose() * form * x).trace()

    def _gram(self, form):
        return matrix(self._field, 4, 4, [self._value_pair(form, a, b) for a in self._matrices for b in self._matrices])

    def _value_pair(self, form, x, y):
        return (x.transpose() * form * y).trace()

    def _short(self, form, bound):
        approximate = self._gram(form).apply_map(self._embedding).change_ring(RDF)
        vectors = pari(approximate).qfminim(float(bound) * 1.0000001, None, 2)[2]
        count = int(pari(vectors).matsize()[1])
        return [tuple(int(round(float(vectors[row, column]))) for row in range(4)) for column in range(count)]

    def _minimum(self, form):
        gram = self._gram(form)
        diagonal = min(gram[n, n] for n in range(4))
        values = [self._value(form, self._matrix(c)) for c in self._short(form, diagonal)]
        return min([diagonal, *(v for v in values if v > 0)])

    def _minimal_vectors(self, form):
        minimum = self._minimum(form)
        return minimum, [c for c in self._short(form, minimum) if self._value(form, self._matrix(c)) == minimum]

    def _ray(self, coordinates):
        x = self._matrix(coordinates)
        product = x * x.transpose()
        return vector(self._field, [product[0, 0], 2 * product[0, 1], product[1, 1]])

    def _from_coordinates(self, direction):
        return matrix(self._field, 2, 2, [direction[0], direction[1], direction[1], direction[2]])

    @staticmethod
    def _positive(form):
        return form[0, 0] > 0 and form.det() > 0

    def _flip(self, form, direction):
        minimum, _ = self._minimal_vectors(form)
        change = self._from_coordinates(direction)
        bound = 2 * minimum
        while True:
            steps = []
            for coordinates in self._short(form, bound):
                x = self._matrix(coordinates)
                slope = self._value(change, x)
                if slope < 0:
                    step = (self._value(form, x) - minimum) / (-slope)
                    if step > 0:
                        steps.append(step)
            if steps:
                candidate = form + min(steps) * change
                if self._positive(candidate) and self._minimum(candidate) == minimum:
                    return candidate
            bound *= 2

    def _perfect(self, form):
        while True:
            minimum, vectors = self._minimal_vectors(form)
            rays = matrix(self._field, [self._ray(c) for c in vectors])
            if rays.rank() == 3:
                return form / minimum
            form = self._flip(form, vector(self._field, rays.right_kernel().basis_matrix()[0]))

    def _facets(self, rays):
        distinct = []
        for ray in rays:
            if not any(matrix(self._field, [ray, other]).rank() == 1 for other in distinct):
                distinct.append(ray)
        normals = []
        for first, second in combinations(distinct, 2):
            normal = vector(self._field, [first[1] * second[2] - first[2] * second[1], first[2] * second[0] - first[0] * second[2], first[0] * second[1] - first[1] * second[0]])
            if normal == 0:
                continue
            pairings = [normal.dot_product(ray) for ray in distinct]
            if all(value <= 0 for value in pairings):
                normal = -normal
            elif not all(value >= 0 for value in pairings):
                continue
            if not any(matrix(self._field, [normal, other]).rank() == 1 for other in normals):
                normals.append(normal)
        return normals

    def _equivalences(self, target, source):
        r"""Return the units ``g`` of the order with ``g^t target g = source``."""
        _, target_vectors = self._minimal_vectors(target)
        _, source_vectors = self._minimal_vectors(source)
        base = self._element(source_vectors[0])
        found = []
        for candidate in target_vectors:
            for sign in (1, -1):
                unit = self._element([sign * c for c in candidate]) * ~base
                coordinates = list(unit)
                if all(c in ZZ for c in coordinates) and abs(unit.reduced_norm()) == 1:
                    image = self._matrix(coordinates)
                    if image.transpose() * target * image == source:
                        found.append(unit)
        return found

    def perfect_forms(self):
        r"""Return representatives of the perfect forms modulo the unit group."""
        return self._enumerate()[0]

    def unit_generators(self):
        r"""Return generators of ``Lambda^x``."""
        return self._enumerate()[1]

    def _enumerate(self):
        representatives = [self._perfect(matrix.identity(self._field, 2))]
        generators = [-self._algebra.one()]
        position = 0
        while position < len(representatives):
            form = representatives[position]
            position += 1
            generators.extend(self._equivalences(form, form))
            _, vectors = self._minimal_vectors(form)
            for normal in self._facets([self._ray(c) for c in vectors]):
                neighbour = self._flip(form, normal)
                neighbour = neighbour / self._minimum(neighbour)
                for known in representatives:
                    transports = self._equivalences(known, neighbour)
                    if transports:
                        generators.append(transports[0])
                        break
                else:
                    representatives.append(neighbour)
        reduced = []
        for unit in generators:
            if unit != 1 and unit not in reduced:
                reduced.append(unit)
        return tuple(representatives), tuple(reduced)

    def reduce_unit(self, unit, generators):
        r"""Return a word in ``generators`` and inverses ending in an element fixing ``phi_S``, by descent."""
        form = self.perfect_forms()[0]
        moves = [*generators, *(~g for g in generators)]
        value = lambda element: self.form_value(form, list(element))
        current = unit
        while True:
            best = min((move * current for move in moves), key=value)
            if value(best) >= value(current):
                return current
            current = best

    def _cell_area(self, form):
        r"""Return the hyperbolic area of the polygon of the minimal vectors of ``form``."""
        _, vectors = self._minimal_vectors(form)
        rays = [self._ray(c) for c in vectors]
        normals = self._facets(rays)
        distinct = []
        for ray in rays:
            if not any(matrix(self._field, [ray, other]).rank() == 1 for other in distinct):
                distinct.append(ray)
        vertices = [ray for ray in distinct if sum(1 for normal in normals if normal.dot_product(ray) == 0) >= 2]
        points = []
        for ray in vertices:
            x, y, z = (float(self._embedding(ray[0])), float(self._embedding(ray[1])) / 2, float(self._embedding(ray[2])))
            root = sqrt(x * z - y * y)
            points.append((x / root, y / root, z / root))
        pairing = lambda p, q: (p[0] * q[2] + p[2] * q[0]) / 2 - p[1] * q[1]
        centre = (sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points), sum(p[2] for p in points) / len(points))
        klein = lambda p: ((p[0] - p[2]) / (p[0] + p[2]), 2 * p[1] / (p[0] + p[2]))
        origin = klein(centre)
        points.sort(key=lambda p: atan2(klein(p)[1] - origin[1], klein(p)[0] - origin[0]))
        total = 0.0
        for index, point in enumerate(points):
            before, after = points[index - 1], points[(index + 1) % len(points)]
            tangent = lambda q: tuple(a - pairing(point, q) * b for a, b in zip(q, point, strict=True))
            first, second = tangent(before), tangent(after)
            total += acos(-pairing(first, second) / sqrt(pairing(first, first) * pairing(second, second)))
        return (len(points) - 2) * pi - total

    def covolume(self):
        r"""Return the area of ``Lambda^x / {1, -1}`` acting on the hyperbolic plane."""
        forms, _ = self._enumerate()
        return sum(self._cell_area(form) / (len(self._equivalences(form, form)) / 2) for form in forms)

    def norm_one_generators(self):
        r"""Return generators of the reduced-norm-one units, with the central kernel."""
        return norm_one_generators(list(self.unit_generators()), [QQ(-1)], 2, lambda x: x.reduced_norm(), QQ)


def orthogonal_group_generators(model, voronoi, normalizer_classes):
    r"""Return generators of ``O(L)`` for the norm lattice ``L`` of the order.

    Every proper isometry is ``x -> u x v^{-1}``. From ``u Lambda v^{-1} = Lambda``, the
    left orders give ``u`` in ``N(Lambda)``, and then ``w = u v^{-1}`` lies in ``Lambda^1``,
    so it equals ``x -> u x u^{-1} w``. ``N(Lambda)`` is generated by ``QQ^x``, ``Lambda^x`` and
    ``normalizer_classes``. Quaternion conjugation gives the other determinant.
    """
    norm_one, central = voronoi.norm_one_generators()
    right = [model.isometry(model.algebra()(unit), "right") for unit in (*norm_one, *central)]
    conjugations = [model.conjugation_by_isometry(unit) for unit in (*voronoi.unit_generators(), *normalizer_classes)]
    return (*right, *conjugations, model.conjugation_isometry())

