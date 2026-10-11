"""Allcock's edgewalk on a hyperbolic lattice, over preamble operations.

A translation of the walk in the pinned polyhedral_common edgewalk source. The walk along the
edges of the fundamental polyhedron of W(L), the cusp handling and the termination certificate
are this module's; Vinberg's chamber, cones, binary-form representations, affine lattice lines,
orders of induced actions and orbit closures are the preamble's. Edge planes are computed by
private rational linear algebra on coordinate rows, and every result is rebuilt as a lattice
element before it reaches a preamble operation.

Conventions: L has signature (n, 1); roots have positive norm; a vertex is a primitive k with
q(k) <= 0, and the walls through it are the roots r with b(r, k) = 0, all others pairing negatively.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key

from dzack_research.preamble.all import ZZ, HyperbolicLattices, Lattices
from dzack_research.preamble.categories.hyperbolic_lattices import AllcockEdgewalkReport
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones
from dzack_research.preamble.categories.vector_configurations import VectorConfigurations
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.edgewalk_extensions import NormedDynkinExtension, finite_type_extensions, parabolic_extensions
from sage_indefinite_port.indefinite.edgewalk_rank2 import anisotropic_cycle
from sage_indefinite_port.readiness import unfinished

type _Vector = Lattices.ElementMethods
type _Roots = tuple[_Vector, ...]
type _Key = tuple[int, ...]
type _RationalRow = FreeModuleElement[Rational]

@unfinished(15)
def edgewalk_fundamental_domain(lattice: HyperbolicLattices.ParentMethods) -> Parent:
    r"""Allcock's edgewalk on a hyperbolic lattice.

    Returns the preamble's ``AllcockEdgewalkReport`` of ``lattice``: the simple roots of a
    fundamental polyhedron of the reflection group W(L), as elements of ``lattice``, and whether
    ``lattice`` is reflective.
    """
    return _Edgewalk(lattice).walk()


def _key(vector_: _Vector) -> _Key:
    coordinates = vector_.to_vector()
    return tuple(int(coordinates(label)) for label in vector_.parent().module_generating_set())


def _oriented_determinant(left: _Key, right: _Key) -> int:
    return left[0] * right[1] - left[1] * right[0]


def _orientation_compare(left: _Key, right: _Key) -> int:
    determinant = _oriented_determinant(left, right)
    if determinant > 0:
        return -1
    if determinant < 0:
        return 1
    return 0


def _signed_sqrt_compare(sign_left: int, square_left: Rational, sign_right: int, square_right: Rational) -> int:
    if sign_left != sign_right:
        return -1 if sign_left < sign_right else 1
    if sign_left == 0:
        return 0
    if square_left == square_right:
        return 0
    if sign_left > 0:
        return -1 if square_left < square_right else 1
    return -1 if square_left > square_right else 1


@dataclass(frozen=True)
class _Vertex:
    generator: _Vector
    roots: _Roots


@dataclass(frozen=True)
class _Candidate:
    alpha: _Vector
    residual_norm: Rational
    root_norm: Rational
    vertex: _Vertex


@dataclass(frozen=True)
class _EdgePlane:
    r"""The rank-two lattice of an edge: the projection Pi of the root sublattice C onto the edge plane (span of the edge's roots)^perp.

    ``plane`` is Pi with its form scaled by ``scale`` to be integral, in the oriented basis ``basis``
    (rational rows of L tensor Q); ``r0`` is the edge's orientation vector in Pi; ``quotient`` is
    Pi / (C cap edge plane).
    """

    constraint: Parent
    basis: Matrix_rational_dense
    plane: Lattices.ParentMethods
    scale: Integer
    r0: _Vector
    quotient: Parent


class _Edgewalk:
    def __init__(self, lattice: HyperbolicLattices.ParentMethods) -> None:
        positive, negative = lattice.signature_pair()
        if int(negative) != 1 or int(positive) + 1 != int(lattice.rank()):
            raise ValueError(f"Allcock's edgewalk needs a nondegenerate lattice of signature (n, 1), but {lattice} has signature {lattice.signature_pair()}")
        self.lattice = lattice
        self.rank = int(lattice.rank())
        self.labels = tuple(lattice.module_generating_set())
        generators = tuple(lattice.module_generators())
        self.gram = matrix(SageQQ, [[SageQQ(int(lattice.b(left, right))) for right in generators] for left in generators])
        self.rational = lattice.vector_space()
        self.field = lattice.base_ring().fraction_field()
        self.norms = tuple(int(norm) for norm in lattice.possible_root_lengths())
        self.identity_key = tuple(_key(generator) for generator in generators)

    def element(self, row: _RationalRow | _Key) -> _Vector:
        ring = self.lattice.base_ring()
        return self.lattice.linear_combination({label: ring(int(entry)) for label, entry in zip(self.labels, row, strict=True) if entry})

    def rational_element(self, row: _RationalRow) -> Lattices.ElementMethods:
        return self.rational.linear_combination(
            {label: self.field(int(SageQQ(entry).numerator())) / self.field(int(SageQQ(entry).denominator())) for label, entry in zip(self.labels, row, strict=True) if entry}
        )

    def row(self, vector_: _Vector) -> _RationalRow:
        return vector(SageQQ, _key(vector_))

    def rational_row(self, vector_: Lattices.ElementMethods) -> _RationalRow:
        coordinates = vector_.to_vector()
        return vector(SageQQ, [SageQQ(int(coordinates(label).numerator())) / SageQQ(int(coordinates(label).denominator())) for label in self.labels])

    def pair(self, left: _RationalRow, right: _RationalRow) -> Rational:
        r"""b on rational rows, for points of L tensor Q that are not lattice vectors."""
        return (left * self.gram * right.column())[0]

    def primitive(self, row: _RationalRow) -> _Vector:
        r"""The primitive vector of L on the ray of a nonzero rational row."""
        denominator = SageZZ.one()
        for entry in row:
            denominator = denominator.lcm(SageQQ(entry).denominator())
        scaled = self.element(tuple(SageZZ(denominator * entry) for entry in row))
        content = int(scaled.content())
        return self.element(tuple(entry // content for entry in _key(scaled)))

    def span_rank(self, vectors: _Roots) -> int:
        return matrix(SageQQ, [self.row(item) for item in vectors]).rank() if vectors else 0

    def perpendicular_rows(self, vectors: tuple[_RationalRow, ...]) -> Matrix_rational_dense:
        return (matrix(SageQQ, list(vectors)) * self.gram).right_kernel_matrix()

    def initial_vertex(self) -> _Vertex:
        r"""A vertex of the chamber cut out by Vinberg's walls: the negative of a chamber ray, whose walls span rank n - 1."""
        for bound in (8, 16, 32, 64, 128):
            chamber = self.lattice.fundamental_chamber(max_roots=bound, max_decompositions=max(256, 8 * bound))
            roots = tuple(chamber.wall_roots())
            candidates: list[_Vertex] = []
            for ray in chamber.primitive_rays():
                generator = -ray
                if generator.q() > 0:
                    continue
                incident = tuple(root for root in roots if self.lattice.b(root, generator) == 0)
                if self.span_rank(incident) == self.rank - 1:
                    candidates.append(_Vertex(generator, incident))
            if candidates:
                return min(candidates, key=lambda vertex: (vertex.generator.q() != 0, tuple(abs(entry) for entry in _key(vertex.generator)), _key(vertex.generator)))
        raise RuntimeError("bounded Vinberg initialization did not expose one chamber corner")

    def constraint(self, norm: Rational) -> Parent:
        r"""C = {v in L : 2 b(v, L) / norm in Z}, the vectors whose reflection of that norm is integral."""
        integral_norm = int(norm)
        return self.lattice.divisible_sublattice(integral_norm // 2 if integral_norm % 2 == 0 else integral_norm)

    def edge_plane(self, roots: _Roots, discarded: _Vector, k: _Vector, norm: Rational) -> _EdgePlane:
        root_rows = tuple(self.row(root) for root in roots)
        plane_basis = self.perpendicular_rows(root_rows)
        if plane_basis.nrows() != 2:
            raise ArithmeticError("an edge must have a two-dimensional perpendicular plane")
        r0_space = self.perpendicular_rows((*root_rows, self.row(k)))
        if r0_space.nrows() != 1:
            raise ArithmeticError("the edge orientation vector must be one-dimensional")
        r0 = r0_space.row(0)
        if self.pair(r0, self.row(discarded)) > 0:
            r0 = -r0
        if k.q() < 0:
            oriented_basis = matrix(SageQQ, [r0, self.row(k)])
        else:
            plane = self.rational.vector_configuration([self.rational_element(row) for row in plane_basis.rows()])
            positive = self.rational_row(plane.inclusion()(plane.vector_of_sign(1)))
            if self.pair(positive, self.row(k)) > 0:
                positive = -positive
            oriented_basis = matrix(SageQQ, [positive, self.row(k)])
            r0 = -self.row(k)

        constraint = self.constraint(norm)
        constraint_rows = matrix(SageQQ, [self.row(constraint.inclusion()(generator)) for generator in constraint.inclusion().domain().module_generators()])
        root_matrix = matrix(SageQQ, list(root_rows))
        root_gram = root_matrix * self.gram * root_matrix.transpose()
        projected = [row - (row * self.gram * root_matrix.transpose()) * root_gram.inverse() * root_matrix for row in constraint_rows.rows()]
        basis = matrix(SageQQ, list((SageQQ**self.rank).span(projected, SageZZ).basis_matrix().rows()))
        if basis.nrows() != 2:
            raise ArithmeticError("the projected constraint lattice must have rank two")
        if oriented_basis.transpose().solve_right(basis.transpose()).transpose().det() < 0:
            basis.rescale_row(0, -1)
        work_gram = basis * self.gram * basis.transpose()
        scale = SageZZ.one()
        for entry in work_gram.list():
            scale = scale.lcm(entry.denominator())
        plane_lattice = Lattices(ZZ)([[ZZ(int(scale * work_gram[row, column])) for column in range(2)] for row in range(2)])
        r0_coordinates = basis.transpose().solve_right(r0.column()).column(0)
        r0_vector = self.plane_primitive(plane_lattice, r0_coordinates)

        perpendicular = self.lattice.vector_configuration(list(roots)).orthogonal_complement()
        intersection = constraint.intersection(perpendicular)
        intersection_rows = matrix(SageQQ, [self.row(intersection.inclusion()(generator)) for generator in intersection.inclusion().domain().module_generators()])
        intersection_coordinates = basis.transpose().solve_right(intersection_rows.transpose()).transpose()
        if any(entry.denominator() != 1 for entry in intersection_coordinates.list()):
            raise ArithmeticError("the perpendicular intersection must be integral in the projected lattice")
        sublattice = plane_lattice.vector_configuration([self.plane_element(plane_lattice, row) for row in intersection_coordinates.rows()])
        return _EdgePlane(constraint, basis, plane_lattice, scale, r0_vector, sublattice.inclusion().cokernel())

    @staticmethod
    def plane_element(plane: Lattices.ParentMethods, coordinates: _RationalRow) -> _Vector:
        ring = plane.base_ring()
        return plane.linear_combination({label: ring(int(entry)) for label, entry in zip(plane.module_generating_set(), coordinates, strict=True) if entry})

    @staticmethod
    def plane_primitive(plane: Lattices.ParentMethods, coordinates: _RationalRow) -> _Vector:
        r"""The primitive vector of the plane lattice on the ray of nonzero rational coordinates."""
        denominator = SageZZ.one()
        for entry in coordinates:
            denominator = denominator.lcm(SageQQ(entry).denominator())
        scaled = _Edgewalk.plane_element(plane, vector(SageQQ, [denominator * entry for entry in coordinates]))
        content = int(scaled.content())
        return _Edgewalk.plane_element(plane, vector(SageQQ, [entry // content for entry in _key(scaled)]))

    def extension_roots(self, roots: _Roots, discarded: _Vector, k: _Vector, extension: NormedDynkinExtension) -> _Roots:
        residual = SageQQ(extension.residual_norm)
        if residual <= 0:
            return ()
        edge = self.edge_plane(roots, discarded, k, extension.norm)
        plane, scale, r0 = edge.plane, edge.scale, edge.r0
        work_vectors: list[_Vector] = []
        if plane.is_isotropic():
            square = scale * residual
            if square.denominator() != 1:
                return ()
            vectors = [item for item in plane.vectors_of_square(ZZ(int(square))) if _oriented_determinant(_key(r0), _key(item)) > 0]
            if k.q() < 0:
                vectors = [item for item in vectors if plane.b(r0, item) > 0]
            vectors.sort(key=lambda item: cmp_to_key(_orientation_compare)(_key(item)))
            work_vectors.extend(vectors)
        else:
            cycle = anisotropic_cycle(plane, scale * SageQQ(extension.norm), r0)
            if cycle is None:
                return ()
            automorph, representatives = cycle
            multiples: list[_Vector] = []
            for item in representatives:
                square = SageQQ(int(item.q()))
                multiplier = 1
                while multiplier * multiplier * square <= scale * extension.norm:
                    candidate = plane.scalar_multiple(plane.base_ring()(multiplier), item)
                    if SageQQ(int(candidate.q())) == scale * residual:
                        multiples.append(candidate)
                    multiplier += 1
            power = plane.Aut().identity()
            for _ in range(int(automorph.order_on(edge.quotient))):
                work_vectors.extend(power(item) for item in multiples)
                power = automorph * power

        component = vector(SageQQ, extension.u_component)
        answer: dict[_Key, _Vector] = {}
        for work_vector in work_vectors:
            ambient = component + edge.basis.transpose() * vector(SageQQ, _key(work_vector))
            if any(entry.denominator() != 1 for entry in ambient):
                continue
            root = self.element(tuple(ambient))
            if not edge.constraint.inclusion().is_in_image(root):
                continue
            if SageQQ(int(root.q())) != extension.norm:
                raise ArithmeticError("rank-two reconstruction changed the root norm")
            answer.setdefault(_key(root), root)
        return tuple(answer.values())

    def candidate_compare(self, left: _Candidate, right: _Candidate, k: _Vector) -> int:
        left_scal = -SageQQ(int(self.lattice.b(left.alpha, k)))
        right_scal = -SageQQ(int(self.lattice.b(right.alpha, k)))
        left_sign = (left_scal > 0) - (left_scal < 0)
        right_sign = (right_scal > 0) - (right_scal < 0)
        first = _signed_sqrt_compare(left_sign, left_scal * left_scal / left.residual_norm, right_sign, right_scal * right_scal / right.residual_norm)
        if first:
            return first
        second = _signed_sqrt_compare(left_sign, left_scal * left_scal / left.root_norm, right_sign, right_scal * right_scal / right.root_norm)
        if second:
            return second
        if left.root_norm == right.root_norm:
            return 0
        return -1 if left.root_norm < right.root_norm else 1

    def vertex_from_root(self, k: _Vector, roots: _Roots, discarded: _Vector, alpha: _Vector) -> _Vertex | None:
        extended = (*roots, alpha)
        kernel = self.perpendicular_rows(tuple(self.row(root) for root in extended))
        if kernel.nrows() != 1:
            return None
        generator = self.primitive(kernel.row(0))
        old_pairing = self.lattice.b(generator, k)
        if old_pairing > 0:
            generator = -generator
            old_pairing = -old_pairing
        if old_pairing == 0 or self.lattice.b(discarded, generator) >= 0:
            return None
        return _Vertex(generator, extended)

    def cusp_root(self, constraint: Parent, component: tuple[Rational, ...], k: _Vector) -> _Vector | None:
        r"""The point u + c k of C on the line through the component u in direction k with the least c > 0."""
        points = constraint.affine_line_points(self.rational_element(vector(SageQQ, component)), k)
        if points is None:
            return None
        base, step = points
        position = next(index for index, entry in enumerate(_key(k)) if entry)

        def coefficient(row: _RationalRow) -> Rational:
            return row[position] / _key(k)[position]

        offset = coefficient(self.row(base) - vector(SageQQ, component))
        increment = coefficient(self.row(step))
        threshold = -offset / increment
        shift = threshold.floor() + 1 if increment > 0 else threshold.ceil() - 1
        root = base + self.lattice.scalar_multiple(self.lattice.base_ring()(int(shift)), step)
        if offset + shift * increment <= 0:
            raise ArithmeticError("the cusp equation must return c > 0")
        return root

    def cusp_roots(self, roots: _Roots, k: _Vector, previous: _Vector) -> _Roots:
        candidates: list[tuple[int, Rational, Rational, _Vector]] = []
        for extension in parabolic_extensions(self.lattice, roots, self.norms):
            root = self.cusp_root(self.constraint(extension.norm), extension.u_component, k)
            if root is None:
                continue
            scalar = -SageQQ(int(self.lattice.b(previous, root)))
            sign = (scalar > 0) - (scalar < 0)
            candidates.append((sign, scalar * scalar / extension.norm, extension.norm, root))

        def compare(left: tuple[int, Rational, Rational, _Vector], right: tuple[int, Rational, Rational, _Vector]) -> int:
            first = _signed_sqrt_compare(left[0], left[1], right[0], right[1])
            if first:
                return first
            if left[2] == right[2]:
                return 0
            return -1 if left[2] < right[2] else 1

        candidates.sort(key=cmp_to_key(compare))
        accepted = list(roots)
        for _sign, _quant, _norm, root in candidates:
            if all(self.lattice.b(old, root) <= 0 for old in accepted):
                accepted.append(root)
        return tuple(accepted)

    def isotropic_next_vertex(self, roots: _Roots, discarded: _Vector, k: _Vector) -> _Vertex:
        perpendicular = self.lattice.vector_configuration(list(roots)).orthogonal_complement()
        plane = perpendicular.inclusion().domain()
        choices: dict[_Key, _Vector] = {}
        for isotropic in plane.primitive_isotropic_vectors():
            candidate = perpendicular.inclusion()(isotropic)
            if self.lattice.b(candidate, k) < 0 and self.lattice.b(discarded, candidate) < 0:
                choices.setdefault(_key(candidate), candidate)
        if len(choices) != 1:
            raise ArithmeticError("the edge must determine one forward isotropic endpoint")
        (generator,) = choices.values()
        return _Vertex(generator, self.cusp_roots(roots, generator, k))

    def edge_step(self, vertex: _Vertex, roots: _Roots, discarded: _Vector) -> _Vertex:
        candidates: list[_Candidate] = []
        for extension in finite_type_extensions(self.lattice, roots, self.norms):
            for alpha in self.extension_roots(roots, discarded, vertex.generator, extension):
                next_vertex = self.vertex_from_root(vertex.generator, roots, discarded, alpha)
                if next_vertex is not None:
                    candidates.append(_Candidate(alpha, extension.residual_norm, extension.norm, next_vertex))
        if candidates:
            candidates.sort(key=cmp_to_key(lambda left, right: self.candidate_compare(left, right, vertex.generator)))
            return candidates[0].vertex
        return self.isotropic_next_vertex(roots, discarded, vertex.generator)

    def edge_directions(self, vertex: _Vertex) -> tuple[tuple[_Roots, _Vector], ...]:
        r"""For each facet of the cone spanned by the vertex's roots, the roots on it and the one root it omits."""
        if len(vertex.roots) == 1:
            return ()
        cone = RationalPolyhedralCones(self.lattice).from_rays(vertex.roots)
        directions: list[tuple[_Roots, _Vector]] = []
        for covector in cone.facet_covectors():
            included = tuple(root for root in vertex.roots if cone.evaluate_covector(covector, root) == 0)
            omitted = tuple(root for root in vertex.roots if cone.evaluate_covector(covector, root) != 0)
            if len(omitted) != 1:
                raise ArithmeticError("Allcock adjacency expects a root-cone facet to omit one wall")
            directions.append((included, omitted[0]))
        return tuple(directions)

    def vertex_configuration(self, vertex: _Vertex) -> tuple[_Roots, tuple[int, ...]]:
        r"""The coloured configuration of the vertex: its roots, the vertex, and at an ideal vertex its neighbours."""
        vectors = [*vertex.roots, vertex.generator]
        colours = [1] * len(vertex.roots) + [2]
        if vertex.generator.q() == 0:
            adjacent: dict[_Key, _Vector] = {}
            for roots, discarded in self.edge_directions(vertex):
                neighbour = self.edge_step(vertex, roots, discarded).generator
                adjacent.setdefault(_key(neighbour), neighbour)
            vectors.extend(adjacent[key] for key in sorted(adjacent))
            colours.extend([3] * len(adjacent))
        return tuple(vectors), tuple(colours)

    def isometry_key(self, isometry: LatticeIsometryMethods) -> tuple[_Key, ...]:
        return tuple(_key(isometry(generator)) for generator in self.lattice.module_generators())

    @staticmethod
    def common_fixed_space(invariant: Parent, isometry: LatticeIsometryMethods) -> Parent | None:
        r"""The sublattice fixed by the group generated so far, or None once that group is infinite.

        The group generated by the symmetries of the polyhedron found so far is finite exactly when
        every generator has finite order and they fix a common vector of negative norm.
        """
        if isometry.multiplicative_order() == Infinity:
            return None
        fixed = invariant.intersection(isometry.invariant_lattice())
        _positive, negative = fixed.inclusion().domain().signature_pair()
        return None if int(negative) == 0 else fixed

    def report(self, vertices: list[_Vertex], reflective: bool, simple_roots: _Roots, isometries: list[LatticeIsometryMethods]) -> Parent:
        cones = RationalPolyhedralCones(self.lattice)
        return AllcockEdgewalkReport(
            self.lattice,
            reflective=reflective,
            simple_roots=simple_roots,
            vertices=tuple(cones.from_rays([vertex.generator], wall_roots=vertex.roots) for vertex in vertices),
            isometry_generators=tuple(isometries),
        )

    def walk(self) -> Parent:
        initial = self.initial_vertex()
        if self.rank == 2:
            if len(initial.roots) != 1:
                raise ArithmeticError("a rank-two chamber corner must have one wall")
            return self.report([initial], True, initial.roots, [])

        vertices = [initial]
        configurations = [self.vertex_configuration(initial)]
        seen: set[tuple[_Key, ...]] = {self.identity_key}
        isometries: list[LatticeIsometryMethods] = []
        invariant = self.lattice.vector_configuration(list(self.lattice.module_generators()))

        def admit(isometry: LatticeIsometryMethods) -> bool:
            r"""Record a new symmetry of the polyhedron; False once the group they generate is infinite."""
            nonlocal invariant
            key = self.isometry_key(isometry)
            if key in seen:
                return True
            seen.add(key)
            isometries.append(isometry)
            fixed = self.common_fixed_space(invariant, isometry)
            if fixed is None:
                return False
            invariant = fixed
            return True

        position = 0
        configuration_category = VectorConfigurations(self.lattice.base_ring())
        while position < len(vertices):
            vertex = vertices[position]
            configuration = configurations[position]
            for automorphism in configuration_category.spanning_automorphism_generators(self.lattice, *configuration):
                if not admit(automorphism):
                    return self.report(vertices, False, (), isometries)
            for roots, discarded in self.edge_directions(vertex):
                neighbour = self.edge_step(vertex, roots, discarded)
                neighbour_configuration = self.vertex_configuration(neighbour)
                isometries_to_neighbour = (configuration_category.spanning_isometry(self.lattice, *known, self.lattice, *neighbour_configuration) for known in configurations)
                transporter = next((isometry for isometry in isometries_to_neighbour if isometry is not None), None)
                if transporter is None:
                    vertices.append(neighbour)
                    configurations.append(neighbour_configuration)
                elif not admit(transporter):
                    return self.report(vertices, False, (), isometries)
            position += 1

        group = self.lattice.O().subgroup(tuple(isometries))
        seeds = list({_key(root): root for vertex in vertices for root in vertex.roots}.values())
        closure = group.orbit_closure(seeds, lambda isometry, root: isometry(root))
        return self.report(vertices, True, tuple(sorted(closure, key=_key)), isometries)


__all__ = ["edgewalk_fundamental_domain"]
