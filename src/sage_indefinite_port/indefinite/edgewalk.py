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
from dzack_research.preamble.categories.definite_lattices import _ExactCVPEngine
from dzack_research.preamble.categories.group.g_sets import FiniteGSets
from dzack_research.preamble.categories.hyperbolic_lattices import AllcockEdgewalkReport
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones
from dzack_research.preamble.categories.sets.set_categories import Sets
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.element import Element
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.edgewalk_extensions import (
    NormedDynkinExtension,
    finite_type_extensions,
    parabolic_extensions,
)
from sage_indefinite_port.indefinite.edgewalk_rank2 import anisotropic_cycle
from sage_indefinite_port.indefinite.edgewalk_scalars import _integer, _rational
from sage_indefinite_port.indefinite.lorentzian_cells import (
    binary_split_vectors_of_square,
)

type _Vector = Lattices.ElementMethods
type _Roots = tuple[_Vector, ...]
type _Key = tuple[int, ...]
type _RationalRow = FreeModuleElement[Rational]


def edgewalk_fundamental_domain(
    lattice: HyperbolicLattices.ParentMethods,
) -> AllcockEdgewalkReport:
    r"""Allcock's edgewalk on a hyperbolic lattice.

    Returns the preamble's ``AllcockEdgewalkReport`` of ``lattice``: the simple roots of a
    fundamental polyhedron of the reflection group W(L), as elements of ``lattice``, and whether
    ``lattice`` is reflective.
    """
    if int(lattice.rank()) > 2 and not lattice.is_isotropic():
        from sage_indefinite_port.indefinite.indefinite_reduction import (
            indefinite_reduced_presentation,
        )

        reduction = indefinite_reduced_presentation(lattice)
        result = _Edgewalk(reduction.reduced).walk()
        inverse = reduction.isometry
        cones = RationalPolyhedralCones(lattice)
        return AllcockEdgewalkReport(
            lattice,
            reflective=result.is_reflective(),
            simple_roots=tuple(inverse(root) for root in result.simple_roots()),
            vertices=tuple(
                cones.from_rays(
                    tuple(inverse(ray) for ray in cone.primitive_rays()),
                    wall_roots=tuple(inverse(root) for root in cone.wall_roots()),
                )
                for cone in result.vertices()
            ),
            isometry_generators=tuple(inverse * generator * inverse.inverse() for generator in result.polyhedron_isometry_generators()),
        )
    return _Edgewalk(lattice).walk()


def _key(vector_: _Vector) -> _Key:
    coordinates = vector_.to_vector()
    return tuple(_integer(coordinates(label)) for label in vector_.parent().module_generating_set())


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


def _primitive_positive_direction(entries) -> tuple[int, ...]:
    """Primitive integral direction, retaining its orientation."""
    from functools import reduce
    from math import gcd

    coordinates = tuple(_integer(entry) for entry in entries)
    divisor = reduce(gcd, coordinates, 0)
    if not divisor:
        return coordinates
    return tuple(entry // divisor for entry in coordinates)


def _same_positive_ray(left, right) -> bool:
    """Exact positive collinearity of nonzero integral vectors."""
    a = tuple(_integer(x) for x in left)
    b = tuple(_integer(x) for x in right)
    pivot = next((i for i, x in enumerate(a) if x), None)
    if pivot is None:
        return False
    return a[pivot] * b[pivot] > 0 and all(a[pivot] * b[i] == b[pivot] * a[i] for i in range(len(a)))


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
    (rational rows of L tensor Q); ``r0`` is the edge's orientation vector in Pi.
    ``constraint_inverse`` inverts the rows of a basis of C, so a row v of L tensor Q lies in C
    exactly when ``v * constraint_inverse`` is integral (edgewalk.h:905, SolutionIntMat(Latt, v_T)).
    ``intersection`` holds the rows, in the coordinates of ``basis``, of a basis of C cap P, the
    sublattice of Pi whose quotient controls the period (edgewalk.h:803, Expr_t).
    """

    constraint_inverse: Matrix_rational_dense
    basis: Matrix_rational_dense
    plane: Lattices.ParentMethods
    scale: Integer
    r0: _Vector
    intersection: Matrix_rational_dense


class _Edgewalk:
    def __init__(self, lattice: HyperbolicLattices.ParentMethods) -> None:
        positive, negative = lattice.signature_pair()
        if int(negative) != 1 or int(positive) + 1 != int(lattice.rank()):
            raise ValueError(f"Allcock's edgewalk needs a nondegenerate lattice of signature (n, 1), but {lattice} has signature {lattice.signature_pair()}")
        self.lattice = lattice
        self.rank = int(lattice.rank())
        self.labels = tuple(lattice.module_generating_set())
        generators = tuple(lattice.module_generators())
        self.gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageQQ)
        self.rational = lattice.vector_space()
        self.field = lattice.base_ring().fraction_field()
        self.norms = tuple(int(norm) for norm in lattice.possible_root_lengths())
        self.identity_key = tuple(_key(generator) for generator in generators)
        # The congruence matrix and its homogeneous integral kernel in
        # Get_DataMapping depend on the root norm, not the affine coset.
        # They are reused for every height of one Vinberg search.
        self._vinberg_mapping_engine = None
        self._vinberg_mapping_cache = {}
        self._vinberg_polytope_cache = {}

    def element(self, row: _RationalRow | _Key) -> _Vector:
        ring = self.lattice.base_ring()
        return self.lattice.linear_combination({label: ring(_integer(entry)) for label, entry in zip(self.labels, row, strict=True) if entry})

    def rational_element(self, row: _RationalRow) -> Lattices.ElementMethods:
        return self.rational.linear_combination(
            {
                label: self.field(int(_rational(entry).numerator())) / self.field(int(_rational(entry).denominator()))
                for label, entry in zip(self.labels, row, strict=True)
                if entry
            }
        )

    def row(self, vector_: _Vector) -> _RationalRow:
        return vector(SageQQ, _key(vector_))

    def rational_row(self, vector_: Lattices.ElementMethods) -> _RationalRow:
        coordinates = vector_.to_vector()
        return vector(
            SageQQ,
            [_rational(coordinates(label)) for label in self.labels],
        )

    def pair(self, left: _RationalRow, right: _RationalRow) -> Rational:
        r"""b on rational rows, for points of L tensor Q that are not lattice vectors."""
        return (left * self.gram * right.column())[0]

    def pairing(self, left: _Vector, right: _Vector) -> Rational:
        return self.pair(self.row(left), self.row(right))

    def norm(self, point: _Vector) -> Rational:
        row = self.row(point)
        return self.pair(row, row)

    def primitive(self, row: _RationalRow) -> _Vector:
        r"""The primitive vector of L on the ray of a nonzero rational row."""
        denominator = SageZZ.one()
        for entry in row:
            denominator = denominator.lcm(_rational(entry).denominator())
        scaled = self.element(tuple(SageZZ(_integer(denominator * entry)) for entry in row))
        content = int(scaled.content())
        return self.element(tuple(entry // content for entry in _key(scaled)))

    def span_rank(self, vectors: _Roots) -> int:
        return matrix(SageQQ, [self.row(item) for item in vectors]).rank() if vectors else 0

    def perpendicular_rows(self, vectors: tuple[_RationalRow, ...]) -> Matrix_rational_dense:
        return (
            matrix(
                SageQQ,
                len(vectors),
                self.rank,
                [entry for row in vectors for entry in row],
            )
            * self.gram
        ).right_kernel_matrix()

    def initial_isotropic_vertex(self, k: _Vector | None = None) -> _Vertex | None:
        r"""The ``isotropic_vinberg`` starting corner of pinned edgewalk.h:1880–2055.

        The source finds a primitive null vector, enumerates each admissible
        root-norm shell in the projection of C_N ∩ k-perp to k-perp/k,
        selects one simple cone, lifts its walls and completes the affine
        cusp roots. The preamble supplies the isotropic reduction and exact
        definite-lattice fixed-square enumeration.
        """
        if k is None:
            k = self.lattice.isotropic_vector()
        reduction = k.isotropic_reduction()
        perpendicular = reduction.orthogonal_complement()
        perpendicular_inclusion = perpendicular.inclusion()
        quotient_projection = reduction.projection()
        rank = int(reduction.rank())
        projected_roots: list[tuple[_Vector, _Vector]] = []

        # Coordinate rows, as for the edge planes: the rows of a basis of
        # k^perp in L, and the rows of their images in the framing of the
        # quotient k^perp/k (the projection pi of edgewalk.h:1955-1968,
        # LatticeProjectionFramework onto the complement of k in k^perp).
        perpendicular_rows = matrix(SageQQ, [self.row(perpendicular_inclusion(generator)) for generator in perpendicular.module_generators()])
        quotient_rows = matrix(SageQQ, [_key(quotient_projection(generator)) for generator in perpendicular.module_generators()])
        k_column = self.gram * self.row(k).column()

        for norm in self.norms:
            # edgewalk.h:1951-1954: Latt_inter_NSP = IntersectionLattice(Latt, NSP),
            # the integral left kernel of B G k on the rows B of a basis of C_N.
            constraint = self.constraint(_rational(norm))
            constraint_rows = matrix(SageQQ, [self.row(constraint.inclusion()(generator)) for generator in constraint.inclusion().domain().module_generators()])
            if constraint_rows.nrows() != self.rank or constraint_rows.rank() != self.rank:
                raise ArithmeticError("the module generators of C_N must be a basis of a full-rank sublattice")
            equations = constraint_rows * k_column
            denominator = SageZZ.one()
            for entry in equations.list():
                denominator = denominator.lcm(entry.denominator())
            intersection_rows = (denominator * equations).change_ring(SageZZ).left_kernel().basis_matrix().change_ring(SageQQ) * constraint_rows
            # fr.BasisProj (edgewalk.h:1955-1960): the projected generators,
            # read through the perpendicular basis and the quotient framing.
            projected_rows = perpendicular_rows.solve_left(intersection_rows) * quotient_rows
            if any(entry.denominator() != 1 for entry in projected_rows.list()):
                raise ArithmeticError("the projection of C_N cap k^perp must be integral in k^perp/k")
            hermite, transformation = projected_rows.change_ring(SageZZ).echelon_form(transformation=True)
            image_rank = hermite.rank()
            image_basis = hermite.matrix_from_rows(range(image_rank)).change_ring(SageQQ)
            lift_rows = transformation.matrix_from_rows(range(image_rank)).change_ring(SageQQ) * intersection_rows
            image = reduction.subobject_on(tuple(self.plane_element(reduction, row) for row in image_basis.rows()))
            for root in image.vectors_of_square(reduction.base_ring()(_integer(norm))):
                projected = image.inclusion()(root)
                # fr.GetOnePreimage (edgewalk.h:2010): one preimage in C_N cap k^perp.
                coefficients = image_basis.solve_left(vector(SageQQ, _key(projected)))
                lifted = self.element(tuple(coefficients * lift_rows))
                if self.norm(lifted) != norm:
                    raise ArithmeticError("the isotropic-reduction lift changed a prescribed root norm")
                projected_roots.append((projected, lifted))

        if not projected_roots:
            # get_simple_cone: an empty projected root system certifies
            # non-reflectivity, rather than an arithmetic failure.
            return None
        labels = tuple(reduction.module_generating_set())

        def quotient_row(point: _Vector) -> _RationalRow:
            coordinates = point.to_vector()
            return vector(SageQQ, [_rational(coordinates(label)) for label in labels])

        projected_rows = tuple(quotient_row(projected) for projected, _root in projected_roots)
        if matrix(SageQQ, projected_rows).rank() != rank:
            # edgewalk.h: the projected roots must span the rank-(n-2)
            # quotient for a reflective cusp.
            return None

        # GetFacetOneDomain_ListIdx (lorentzian_linalg.h:267–339):
        # regular splitting vector, positive half and its extreme rays.
        # A regular splitting is outside the finitely many hyperplanes
        # annihilating the nonzero projected roots. On the moment curve
        # (1,t,...,t^(rank-1)), each forbidden hyperplane has at most
        # rank-1 integer parameters t. Exhausting t=0,1,... therefore
        # constructs a regular splitting without probabilistic retries.
        if any(not any(row) for row in projected_rows):
            raise ArithmeticError("the projected root list contains a zero vector")
        parameter = 0
        while True:
            splitting = vector(SageQQ, [SageQQ(int(parameter)) ** degree for degree in range(rank)])
            if all(splitting.dot_product(row) != 0 for row in projected_rows):
                break
            parameter += 1
        positive = [(projected, root) for (projected, root), row in zip(projected_roots, projected_rows, strict=True) if splitting.dot_product(row) > 0]
        cone = RationalPolyhedralCones(reduction).from_rays([projected for projected, _root in positive])
        simple_roots: list[_Vector] = []
        simple_projected_rows: list[_RationalRow] = []
        for ray in cone.primitive_rays():
            ray_row = quotient_row(ray)
            for projected, root in positive:
                root_row = quotient_row(projected)
                if matrix(SageQQ, [ray_row, root_row]).rank() == 1 and reduction.b(ray, projected) > 0:
                    simple_roots.append(root)
                    simple_projected_rows.append(root_row)
                    break
            else:
                raise ArithmeticError("a simple cone ray has no original admissible root lift")
        # Source get_simple_cone accepts non-simplicial faces. Test rank in
        # the isotropic quotient, not in the ambient lattice of their lifts.
        if matrix(SageQQ, simple_projected_rows).rank() != rank:
            raise ArithmeticError("the isotropic starting root cone does not span the projected lattice")

        # edgewalk.h:2020–2049: choose a negative kP orthogonal to the
        # simple walls, orient by b(kP,k)<0 and complete the cusp by the
        # source's affine-line equation fibres.
        plane_rows = self.perpendicular_rows(tuple(self.row(root) for root in simple_roots))
        plane = self.rational.subobject_on([self.rational_element(row) for row in plane_rows.rows()])
        negative = self.rational_row(plane.inclusion()(plane.vector_of_sign(-1)))
        if self.pair(negative, self.row(k)) > 0:
            negative = -negative
        previous = self.primitive(negative)
        return _Vertex(k, self.cusp_roots(tuple(simple_roots), k, previous))

    def _vinberg_polytope_roots(self, engine, shift, directions, positive_gram, center, radius, root_norm, walls) -> tuple:
        r"""Exact full-rank ``FindRoot_filter/fct_CVP_Poly`` branch.

        This is the Bareiss spherical recursion and LP facet pruning of
        Shvec_exact.h:334-599, called by vinberg.h:611-697. Sage's PPL
        backend provides the source's exact linear-programming primitive.
        """
        dimension = positive_gram.nrows()
        prepared = self._vinberg_polytope_cache.get(root_norm)
        if prepared is None:
            # LLL_gram returns U with U^t G U reduced. Source Pmat = U^t.
            # Both this change of coordinates and FullGramInfo's exact
            # Bareiss data depend on the root-norm kernel, not on the
            # affine shift, height, or current system of facet constraints.
            change = positive_gram.LLL_gram()
            if abs(change.det()) != 1:
                raise ArithmeticError("Vinberg's spherical-coordinate change is not unimodular")
            transformed_directions = directions * change
            bareiss = (change.transpose() * positive_gram * change).change_ring(SageQQ)
            minors = []
            previous = SageQQ.one()
            for i in range(dimension):
                pivot = bareiss[i, i]
                if pivot <= 0:
                    raise ArithmeticError("Vinberg's perpendicular Gram is not positive definite")
                minors.append(pivot)
                for row in range(i + 1, dimension):
                    for column in range(i + 1, dimension):
                        bareiss[row, column] = (bareiss[row, column] * pivot - bareiss[row, i] * bareiss[i, column]) / previous
                previous = pivot
            prepared = (positive_gram, change, transformed_directions, bareiss, tuple(minors))
            self._vinberg_polytope_cache[root_norm] = prepared
        original_gram, change, transformed_directions, bareiss, minors = prepared
        if original_gram != positive_gram or transformed_directions != directions * change:
            raise ArithmeticError("Vinberg's cached spherical kernel changed inside one root norm")
        coset = change.inverse() * (-center)
        ambient_gram = matrix(SageZZ, self.gram)
        inequalities = []
        for wall in walls:
            covector = -(ambient_gram * vector(SageZZ, wall))
            inequalities.append((SageQQ(covector.dot_product(shift)),) + tuple(SageQQ(covector.dot_product(column)) for column in transformed_directions.columns()))

        def last_true(predicate):
            # Shvec_LastTrue: exact doubling and bisection.
            if predicate(0):
                lower, upper = 0, 1
                while predicate(upper):
                    lower, upper = upper, 2 * upper
            else:
                lower, upper = -1, 0
                while not predicate(lower):
                    upper, lower = lower, 2 * lower
            while upper - lower > 1:
                middle = (lower + upper) // 2
                if predicate(middle):
                    lower = middle
                else:
                    upper = middle
            return lower

        def spherical_bounds(square, offset, denominator):
            # Exact Bound_Ceil / Bound_Floor, with no inexact seed.
            def at_most_upper(index):
                value = denominator * index + offset
                return value <= 0 or value * value <= square

            def at_least_lower(negative_index):
                value = offset - denominator * negative_index
                return value >= 0 or value * value <= square

            return -last_true(at_least_lower), last_true(at_most_upper)

        coordinates = [SageZZ.zero()] * dimension
        # The section at level i substitutes precisely the coordinates
        # j > i already chosen by the spherical recursion.  Carry these
        # constant terms down the search tree, just as the fraction-free
        # Shvec recursion carries its affine offset.  Constructing them
        # afresh from every suffix for each LP node redundantly recomputes
        # all earlier coordinate assignments.
        section_constants = [row[0] for row in inequalities]
        roots = []

        def enumerate_level(level, remaining):
            offset = minors[level] * coset[level] + sum(
                (bareiss[level, j] * (coordinates[j] + coset[j]) for j in range(level + 1, dimension)),
                SageQQ.zero(),
            )
            square = remaining * (minors[level - 1] if level else 1)
            if square < 0:
                return
            lower, upper = spherical_bounds(square, offset, minors[level])
            section = Polyhedron(
                ieqs=[(constant,) + row[1 : level + 2] for constant, row in zip(section_constants, inequalities, strict=True)],
                base_ring=SageQQ,
                backend="ppl",
            )
            if section.is_empty():
                return
            vertices = tuple(vertex.vector()[level] for vertex in section.vertices())
            rays = tuple(ray.vector()[level] for ray in section.rays())
            lines = tuple(line.vector()[level] for line in section.lines())
            # A nonempty section need not be pointed.  For example, a
            # bounded strip times a line has no extreme points, although
            # its intersection with the positive-definite sphere is finite.
            # Such a section cannot supply a vertex-derived coordinate
            # bound.  Retain the exhaustive spherical bounds instead of
            # discarding its valid integral solutions.
            if not vertices and not lines:
                raise ArithmeticError("a nonempty pointed Vinberg facet section has no vertex")
            if vertices and not any(ray < 0 for ray in rays) and not any(line for line in lines):
                lower = max(lower, int(min(vertices).ceil()))
            if vertices and not any(ray > 0 for ray in rays) and not any(line for line in lines):
                upper = min(upper, int(max(vertices).floor()))
            for value in range(lower, upper + 1):
                coordinates[level] = SageZZ(value)
                term = minors[level] * coordinates[level] + offset
                if level:
                    next_remaining = (minors[level - 1] * remaining - term * term) / minors[level]
                    if next_remaining >= 0:
                        for index, row in enumerate(inequalities):
                            section_constants[index] += row[level + 1] * value
                        enumerate_level(level - 1, next_remaining)
                        for index, row in enumerate(inequalities):
                            section_constants[index] -= row[level + 1] * value
                    continue
                enumerated_square = radius - (remaining - term * term) / minors[0]
                if enumerated_square != radius:
                    continue
                candidate = shift + transformed_directions * vector(SageZZ, coordinates)
                if candidate.dot_product(ambient_gram * candidate) != root_norm:
                    raise ArithmeticError("Vinberg's constrained spherical norm witness is invalid")
                if any((2 * pairing) % root_norm for pairing in ambient_gram * candidate):
                    raise ArithmeticError("Vinberg's constrained root-integrality witness is invalid")
                if any(candidate.dot_product(ambient_gram * vector(SageZZ, wall)) > 0 for wall in walls):
                    raise ArithmeticError("Vinberg's constrained facet witness is invalid")
                roots.append(engine.V(candidate))

        enumerate_level(dimension - 1, minors[-1] * radius)
        return tuple(roots)

    def _vinberg_decomposed_roots(self, engine, component, root_norm: int, walls=None) -> tuple:
        r"""Exact ``Get_DataMapping`` and spherical solutions for one Vinberg coset.

        The source solves ``2 G(a+Mw) = k h`` over the integers before
        enumerating positive-definite solutions.  Its integer solution and
        homogeneous kernel give the *entire* admissible affine coset.  The
        preamble's exact PARI CVP backend enumerates that coset's spherical
        fibre, with an exact norm and reflection-divisibility certificate.

        See polyhedral_common src_lorentzian/vinberg.h:407-566.  The later
        full-rank facet-constrained CVP optimization is separate from this
        congruence reduction.
        """
        gram = matrix(SageZZ, self.gram)
        n = gram.nrows()
        a = vector(SageZZ, component)
        if self._vinberg_mapping_engine is not engine:
            self._vinberg_mapping_engine = engine
            self._vinberg_mapping_cache.clear()
            self._vinberg_polytope_cache.clear()
        prepared = self._vinberg_mapping_cache.get(root_norm)
        if prepared is None:
            orthogonal = matrix(SageZZ, tuple(engine.V1.gens())).transpose()
            if orthogonal.nrows() != n or orthogonal.ncols() != n - 1:
                raise ArithmeticError("Vinberg's perpendicular lattice has incorrect rank")
            if root_norm in (1, 2):
                diagonal = left = right = None
                congruence = None
                directions = orthogonal
            else:
                # A = [2GM | -kI] has full row rank, and its Smith
                # transformations are independent of the affine offset a.
                # The last n-1 columns of the right transformation give an
                # integral kernel basis for every congruence fibre.
                congruence = (2 * gram * orthogonal).augment(-root_norm * identity_matrix(SageZZ, n))
                diagonal, left, right = congruence.smith_form()
                kernel = right.matrix_from_rows_and_columns(range(n - 1), range(n, 2 * n - 1))
                directions = orthogonal * kernel
                if directions.rank() != n - 1 or not (congruence * right.matrix_from_columns(range(n, 2 * n - 1))).is_zero():
                    raise ArithmeticError("Vinberg's integral congruence kernel is invalid")
            positive_gram = directions.transpose() * gram * directions
            if not positive_gram.is_positive_definite():
                raise ArithmeticError("Vinberg's perpendicular form is not positive definite")
            backend = _ExactCVPEngine._from_positive_engine_gram(self.lattice.base_ring(), positive_gram)
            prepared = (orthogonal, congruence, diagonal, left, right, directions, positive_gram, backend)
            self._vinberg_mapping_cache[root_norm] = prepared
        orthogonal, congruence, diagonal, left, right, directions, positive_gram, backend = prepared
        if root_norm in (1, 2):
            shift = a
        else:
            # Sage smith_form returns (D,P,Q) with P*A*Q=D.  Since
            # -k*I is a block of A, its row rank is n.  The final n-1
            # columns of Q project to the *integral* kernel in w, not
            # merely the rational nullspace of the congruence.
            residual = left * (-2 * gram * a)
            if any(residual[i] % diagonal[i, i] for i in range(n)):
                return ()
            solution = vector(
                SageZZ,
                [residual[i] // diagonal[i, i] for i in range(n)] + [0] * (n - 1),
            )
            integer_solution = right * solution
            if congruence * integer_solution != -2 * gram * a:
                raise ArithmeticError("Vinberg's integral congruence witness is invalid")
            w0 = vector(SageZZ, integer_solution[: n - 1])
            shift = a + orthogonal * w0
            if any((2 * pairing) % root_norm for pairing in gram * shift):
                raise ArithmeticError("Vinberg's affine congruence solution is not integral")
        linear = directions.transpose() * gram * shift
        center = -positive_gram.change_ring(SageQQ).solve_right(linear)
        radius = SageQQ(root_norm - shift.dot_product(gram * shift)) + center.dot_product(positive_gram * center)
        if radius < 0:
            return ()
        if walls and matrix(SageZZ, tuple(walls)).rank() == n:
            return self._vinberg_polytope_roots(engine, shift, directions, positive_gram, center, radius, root_norm, walls)
        roots = []
        for coordinates, square in backend.close_vector_coordinates(center, radius):
            if square != radius:
                continue
            row = shift + directions * vector(SageZZ, coordinates)
            if row.dot_product(gram * row) != root_norm:
                raise ArithmeticError("Vinberg's spherical solution has incorrect norm")
            if any((2 * pairing) % root_norm for pairing in gram * row):
                raise ArithmeticError("Vinberg's spherical solution is not a root")
            roots.append(engine.V(row))
        return tuple(roots)

    def initial_vertex(self) -> _Vertex | None:
        r"""Find a corner, first checking the finite root system at a timelike point.

        A full-rank perpendicular root system makes the point a corner of
        some chamber. Otherwise Vinberg's increasing prefixes locate a corner;
        no finite search bound is used as a nonexistence certificate.
        """
        if self.lattice.is_isotropic():
            return self.initial_isotropic_vertex()
        timelike = self.lattice.negative_vector()
        roots = self.ordinary_vertex_roots(timelike)
        if self.span_rank(roots) == self.rank - 1:
            # ordinary_vertex_roots has already enumerated the complete
            # perpendicular root family and selected its cone rays.
            # The require_full path would enumerate all fixed-norm vectors
            # again, although its sole remaining condition is this count.
            if len(roots) != self.rank - 1:
                raise ArithmeticError("ordinary initial cone has incorrect simple-root count")
            return _Vertex(timelike, roots)
        # vinberg.h FindRoots_Kernel / FindOneInitialRay: examine the accepted
        # root system after its initial finite cone and after each new root,
        # stopping on the first feasible nonpositive ray. The preamble's
        # VinAl backend supplies the actual root-decomposition enumeration.
        from vinal import IterateIntegerPoints, VinAl

        engine = VinAl(self.gram, vector(SageZZ, _key(timelike)), False, None, initialize_cosets=False)
        # vinberg.h GetVinbergAux: enumerate translation classes of the
        # saturated perpendicular lattice together with v0, and reduce
        # each representative by the quotient of -<w,v0> and -v0^2.
        # VinAl sorts the unnormalized representatives, which changes the
        # source's decomposition sequence; reconstruct W before searching.
        coset_basis = matrix(SageZZ, tuple(engine.V1.gens()) + (engine.v0,))
        norm_shift = -engine.v0.inner_product(engine.v0)
        engine.W = []
        for coordinates in IterateIntegerPoints(coset_basis):
            representative = engine.V(coordinates)
            quotient = (-representative.inner_product(engine.v0)) // norm_shift
            normalized = representative - quotient * engine.v0
            # GetVinbergAux uses QuoInt/ResInt to place each class in
            # the fundamental timelike strip, 0 <= -<w,v0> < -v0^2.
            height = -normalized.inner_product(engine.v0)
            if not 0 <= height < norm_shift:
                raise ArithmeticError("Vinberg coset representative leaves timelike strip")
            engine.W.append(normalized)
        # vinberg.h GetVinbergFromG supplies the caller's root-length set.
        # VinAl otherwise derives its own list from the form determinant,
        # which need not equal the admissible lengths of this lattice.
        engine.root_lengths = list(self.norms)
        # vinberg.h FundCone -> GetFacetOneDomain: select one member of
        # each perpendicular root pair by a regular linear functional,
        # and retain the original normed roots on the extreme cone rays.
        # Lexicographic sign is a deterministic regular functional for this
        # finite integral root family, unlike VinAl.FundCone's primitive
        # ray representatives.
        initial_roots = [root for norm in engine.root_lengths for root in self._vinberg_decomposed_roots(engine, engine.V([0] * engine.n), norm)]
        # GetFacetOneDomain first reduces the root coordinate columns to
        # their actual span, and GetFacetOneDomain_ListIdx requires at least
        # twice that reduced dimension in the nonempty root family.
        if initial_roots:
            reduced_dimension = matrix(SageQQ, initial_roots).rank()
            if len(initial_roots) < 2 * reduced_dimension:
                raise ArithmeticError("initial Vinberg root family has fewer than twice its rank in roots")
        positive = [root for root in initial_roots if next((entry for entry in root if entry), 0) > 0]
        # GetFacetOneDomain_ListIdx allocates exactly half the roots on
        # the assumption that the family is centrally symmetric. A count
        # alone cannot establish the required opposite-root correspondence.
        if initial_roots:
            root_rows = [tuple(_integer(entry) for entry in root) for root in initial_roots]
            multiplicities = {}
            for row in root_rows:
                multiplicities[row] = multiplicities.get(row, 0) + 1
            if any(multiplicities.get(row, 0) != multiplicities.get(tuple(-entry for entry in row), 0) for row in multiplicities):
                raise ArithmeticError("Vinberg perpendicular roots do not form opposite pairs")
        if positive:
            positive_cone = RationalPolyhedralCones(self.lattice).from_rays(self.element(tuple(_integer(entry) for entry in root)) for root in positive)
            rays = tuple(self.row(ray) for ray in positive_cone.primitive_rays())
            # GetFacetOneDomain_ListIdx selects nonredundant *indices*;
            # parallel root vectors cannot both survive its reduction.
            chosen = []
            used_rays = set()
            ray_by_direction = {_primitive_positive_direction(ray): index for index, ray in enumerate(rays)}
            for root in positive:
                index = ray_by_direction.get(_primitive_positive_direction(root))
                if index is not None and index not in used_rays:
                    chosen.append(root)
                    used_rays.add(index)
            if len(used_rays) != len(rays):
                raise ArithmeticError("initial Vinberg cone ray has no normed root")
            engine.roots = chosen
        else:
            engine.roots = []

        def accepted_vertex() -> _Vertex | None:
            roots = tuple(self.element(tuple(_integer(entry) for entry in root)) for root in engine.roots)
            if not roots:
                return None
            equations = matrix(SageQQ, [self.row(root) * self.gram for root in roots])
            rank = equations.rank()
            if rank < self.rank - 1:
                return None
            if rank == self.rank - 1:
                generator = self.primitive(equations.right_kernel_matrix().row(0))
                if self.pairing(generator, timelike) > 0:
                    generator = -generator
                if self.norm(generator) > 0 or self.pairing(generator, timelike) >= 0:
                    return None
            else:
                # GetOneInteriorVertex: exact dual description of -R*G,
                # tested on each extreme ray in the selected time component.
                # GetOneInteriorVertex uses the covectors -G*r. The
                # lattice's algebraic correlation produces these directly,
                # without assuming that primal and dual generator labels agree.
                inequalities = tuple(-root.to_covector() for root in roots)
                cone = RationalPolyhedralCones(self.lattice).from_halfspaces(inequalities)
                # GetOneInteriorVertex (vinberg.h:833-909) accepts the first
                # nonspacelike extreme ray of the feasible cone, regardless
                # of its time orientation.  The perpendicular-root chamber
                # was selected by an independent regular functional, so it
                # need not lie in the component containing the initial v0.
                generator = next((ray for ray in cone.primitive_rays() if self.norm(ray) <= 0), None)
                if generator is None:
                    return None
            if self.norm(generator) == 0:
                return self.initial_isotropic_vertex(generator)
            return _Vertex(generator, self.ordinary_vertex_roots(generator, require_full=True))

        initial = accepted_vertex()
        if initial is not None:
            return initial
        # vinberg.h FindRoots_Kernel tests its exit only after every root
        # belonging to one decomposition request has been inserted. VinAl's
        # NextRoot yields those roots individually; use the backend's called
        # decomposition helpers so an intermediate root cannot end the search.
        # vinberg.h IterateRootDecompositions::get_k (pinned 1592b246)
        # prioritizes each length's next coset by val/sqrt(val), with
        # val = -<v0, a>. VinAl instead uses val/sqrt(k), so its iterator
        # does not reproduce the C++ candidate sequence.
        candidates = {norm: 1 for norm in engine.root_lengths}
        coset_heights = tuple(-engine.v0.inner_product(w) for w in engine.W)
        period = len(coset_heights)
        if not period:
            raise ArithmeticError("Vinberg decomposition has no coset representatives")
        while True:
            norm = min(
                engine.root_lengths,
                key=lambda length: coset_heights[candidates[length] % period] + (candidates[length] // period) * norm_shift,
            )
            index = candidates[norm]
            component = engine.W[index % period] + (index // period) * engine.v0
            candidates[norm] = index + 1
            # vinberg.h FindRoot_filter/fct_CVP retains every norm-k
            # candidate satisfying all accepted wall inequalities. VinAl's
            # IsNewRoot additionally rejects linearly dependent candidates;
            # that is not the upstream filter and can discard actual walls.
            new_roots = [
                root for root in self._vinberg_decomposed_roots(engine, component, norm, walls=engine.roots) if all(root.inner_product(wall) <= 0 for wall in engine.roots)
            ]
            if not new_roots:
                continue
            engine.roots.extend(new_roots)
            # vinberg.h ReduceListRoot removes generators redundant in the
            # accepted root cone. Keep the original normed root vectors,
            # using primitive rays only to identify their directions.
            ambient_roots = tuple(self.element(tuple(_integer(entry) for entry in root)) for root in engine.roots)
            extreme = RationalPolyhedralCones(self.lattice).from_rays(ambient_roots).primitive_rays()
            # ReduceListRoot/get_non_redundant_indices preserves the order
            # of the original family. Keep its first normed representative
            # of each extremal direction rather than reordering by rays.
            ray_rows = tuple(self.row(ray) for ray in extreme)
            reduced = []
            retained_directions = set()
            # Compare primitive directions once per accepted root rather
            # than testing every root against every extremal ray.
            ray_by_direction = {_primitive_positive_direction(row): index for index, row in enumerate(ray_rows)}
            for root in engine.roots:
                direction = _primitive_positive_direction(root)
                index = ray_by_direction.get(direction)
                if index is not None and index not in retained_directions:
                    reduced.append(root)
                    retained_directions.add(index)
            if len(retained_directions) != len(ray_rows):
                raise ArithmeticError("Vinberg cone ray has no accepted normed root")
            engine.roots = reduced
            initial = accepted_vertex()
            if initial is not None:
                return initial

    def ordinary_vertex_roots(self, timelike: _Vector, *, require_full: bool = False) -> _Roots:
        """A simple system of all ambient reflective roots perpendicular to k.

        Lexicographic positivity is a regular ordering on this finite root
        system. Its positive-root cone has precisely the simple-root rays.
        """
        inclusion = timelike.orthogonal_complement().inclusion()
        complement = inclusion.domain()
        positive: list[_Vector] = []
        for norm in self.norms:
            for point in complement.vectors_of_square(complement.base_ring()(_integer(norm))):
                root = inclusion(point)
                if next(entry for entry in _key(root) if entry) < 0:
                    continue
                if all((2 * int(self.pairing(root, basis))) % norm == 0 for basis in self.lattice.module_generators()):
                    positive.append(root)
        if not positive:
            if require_full:
                # The upstream ordinary-cone constructor rejects a root
                # family whose facet count is not dimSpace (including zero).
                raise ArithmeticError("ordinary initial cone has incorrect simple-root count")
            return ()
        # edgewalk.h get_simple_cone_from_lattice selects actual normed
        # roots from FindFixedDistVectors, not primitive multiples of them.
        # Use geometric rays only to identify which original roots span
        # the selected cone, retaining each root's prescribed norm.
        rays = RationalPolyhedralCones(self.lattice).from_rays(positive).primitive_rays()
        if require_full and len(rays) != self.rank - 1:
            # edgewalk.h get_simple_cone_from_lattice requires exactly
            # dimSpace selected walls for an ordinary starting vertex.
            raise ArithmeticError("ordinary initial cone has incorrect simple-root count")
        roots_by_direction = {}
        for root in positive:
            roots_by_direction.setdefault(_primitive_positive_direction(_key(root)), root)
        return tuple(roots_by_direction[_primitive_positive_direction(_key(ray))] for ray in rays)

    def constraint(self, norm: Rational) -> Parent:
        r"""C = {v in L : 2 b(v, L) / norm in Z}, the vectors whose reflection of that norm is integral."""
        integral_norm = int(norm)
        return self.constraint_for_divisor(integral_norm // 2 if integral_norm % 2 == 0 else integral_norm)

    @cached_method
    def constraint_for_divisor(self, divisor: int) -> Parent:
        return self.lattice.divisible_sublattice(divisor)

    def edge_plane(self, roots: _Roots, discarded: _Vector, k: _Vector, norm: Rational) -> _EdgePlane:
        integral_norm = int(norm)
        return self.edge_plane_for_divisor(
            roots,
            discarded,
            k,
            integral_norm // 2 if integral_norm % 2 == 0 else integral_norm,
        )

    @cached_method
    def edge_perpendicular(self, roots: _Roots) -> Parent:
        r"""The sublattice L cap P of the edge plane P = (span of the edge's roots)^perp.

        edgewalk.h:598 computes Pplane once per edge; only Latt = C_N varies
        with the root norm (edgewalk.h:761-764, IntersectionLattice_VectorSpace).
        """
        return self.lattice.subobject_on(list(roots)).orthogonal_complement()

    @cached_method
    def edge_plane_for_divisor(self, roots: _Roots, discarded: _Vector, k: _Vector, divisor: int) -> _EdgePlane:
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
        if self.norm(k) < 0:
            oriented_basis = matrix(SageQQ, [r0, self.row(k)])
        else:
            plane = self.rational.subobject_on([self.rational_element(row) for row in plane_basis.rows()])
            positive = self.rational_row(plane.inclusion()(plane.vector_of_sign(1)))
            if self.pair(positive, self.row(k)) > 0:
                positive = -positive
            oriented_basis = matrix(SageQQ, [positive, self.row(k)])
            r0 = -self.row(k)

        constraint = self.constraint_for_divisor(divisor)
        constraint_rows = matrix(
            SageQQ,
            [self.row(constraint.inclusion()(generator)) for generator in constraint.inclusion().domain().module_generators()],
        )
        wall_matrix = matrix(
            SageQQ,
            len(root_rows),
            self.rank,
            [entry for row in root_rows for entry in row],
        )
        # The walls through an edge need not be linearly independent.  The
        # projection only depends on their rational span; invert the pairing
        # on an independent basis of that span, not the (possibly singular)
        # Gram matrix of every wall.  The edge plane is Lorentzian, so its
        # orthogonal complement must be positive definite and nondegenerate.
        root_matrix = wall_matrix.row_space().basis_matrix()
        root_gram = root_matrix * self.gram * root_matrix.transpose()
        if not root_gram.is_invertible():
            raise ArithmeticError("the edge-wall span must be nondegenerate")
        projected = [row - (row * self.gram * root_matrix.transpose()) * root_gram.inverse() * root_matrix for row in constraint_rows.rows()]
        basis = matrix(
            SageQQ,
            list((SageQQ**self.rank).span(projected, SageZZ).basis_matrix().rows()),
        )
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

        # edgewalk.h:763-764, Basis_P_inter_LN = IntersectionLattice_VectorSpace(Latt, Pplane):
        # C cap P on the same coordinate rows as the projection above. A
        # vector c * B of C (B the rows of a basis of C) lies in P exactly
        # when it pairs to zero with every edge wall, so the coefficient
        # lattice of C cap P is the integral left kernel of B G R^T.
        if constraint_rows.nrows() != self.rank or constraint_rows.rank() != self.rank:
            raise ArithmeticError("the module generators of C_N must be a basis of a full-rank sublattice")
        equations = constraint_rows * self.gram * wall_matrix.transpose()
        denominator = SageZZ.one()
        for entry in equations.list():
            denominator = denominator.lcm(entry.denominator())
        coefficients = (denominator * equations).change_ring(SageZZ).left_kernel().basis_matrix()
        intersection_rows = coefficients.change_ring(SageQQ) * constraint_rows
        if intersection_rows.nrows() != 2:
            raise ArithmeticError("C cap P must have rank two")
        intersection_coordinates = basis.transpose().solve_right(intersection_rows.transpose()).transpose()
        if any(entry.denominator() != 1 for entry in intersection_coordinates.list()):
            raise ArithmeticError("the perpendicular intersection must be integral in the projected lattice")
        return _EdgePlane(
            constraint_rows.inverse(),
            basis,
            plane_lattice,
            scale,
            r0_vector,
            intersection_coordinates,
        )

    @staticmethod
    def plane_element(plane: Lattices.ParentMethods, coordinates: _RationalRow) -> _Vector:
        ring = plane.base_ring()
        return plane.linear_combination({label: ring(_integer(entry)) for label, entry in zip(plane.module_generating_set(), coordinates, strict=True) if entry})

    @staticmethod
    def plane_primitive(plane: Lattices.ParentMethods, coordinates: _RationalRow) -> _Vector:
        r"""The primitive vector of the plane lattice on the ray of nonzero rational coordinates."""
        denominator = SageZZ.one()
        for entry in coordinates:
            denominator = denominator.lcm(_rational(entry).denominator())
        scaled = _Edgewalk.plane_element(plane, vector(SageQQ, [denominator * entry for entry in coordinates]))
        entries = _key(scaled)
        content = SageZZ(int(scaled.content()))
        if content == 0:
            raise ValueError("the zero vector does not determine a primitive ray")
        primitive = vector(SageZZ, (SageZZ(int(entry)).quo_rem(content)[0] for entry in entries))
        if any(SageZZ(int(entry)) % content for entry in entries):
            raise ArithmeticError("content does not divide the integral ray coordinates")
        return _Edgewalk.plane_element(plane, primitive)

    @cached_method
    def _anisotropic_edge_candidate_family(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        norm: Rational,
    ) -> _Roots:
        """The complete anisotropic edge-plane family for one root length.

        polyhedral_common edgewalk.h:753-845 computes Anisotropic once per
        e_norm, scales *each* returned primitive vector by all positive
        integers remaining under e_norm, and unfolds the finite quotient
        period.  Only get_next_anisotropic subsequently filters res_norm.
        Reusing the full family across extensions also preserves that
        original enumeration order.
        """
        edge = self.edge_plane(roots, discarded, k, norm)
        plane, scale, r0 = edge.plane, edge.scale, edge.r0
        cycle = anisotropic_cycle(plane, scale * _rational(norm), r0)
        if cycle is None:
            return ()
        automorph, representatives = cycle
        bound = _rational(scale * norm)
        multiples: list[_Vector] = []
        for item in representatives:
            square = _rational(item.q())
            if square <= 0:
                raise ArithmeticError("Anisotropic must return primitive positive-square vectors")
            multiplier = 1
            while multiplier * multiplier * square <= bound:
                candidate = plane.scalar_multiple(plane.base_ring()(multiplier), item)
                if _rational(candidate.q()) != multiplier * multiplier * square:
                    raise ArithmeticError("anisotropic scalar multiplication changed the square")
                multiples.append(candidate)
                multiplier += 1

        # lorentzian_linalg.h:70-143,
        # GetMatrixExponentSublattice_TrivClass: first find the least
        # positive power preserving C cap Pi, then the order of its
        # action on Pi/(C cap Pi). The source tests integral membership
        # in the sublattice and then the induced finite quotient action,
        # not the generally larger order in GL(2, Z/exp(Pi/S)).
        # On coordinate rows of Pi the period acts by x -> x * A, A the rows
        # of the images of the basis (edgewalk.h:833, TheMat.transpose() * v).
        # With M the rows of a basis of S = C cap Pi, h preserves S exactly
        # when M h M^{-1} is integral, and acts trivially on Pi/S exactly
        # when (h - 1) M^{-1} is integral.
        action = matrix(SageQQ, [_key(automorph(point)) for point in plane.module_generators()])
        sublattice = edge.intersection
        sublattice_inverse = sublattice.inverse()
        identity = identity_matrix(SageQQ, 2)

        def integral(product: Matrix_rational_dense) -> bool:
            return all(entry.denominator() == 1 for entry in product.list())

        period = action
        preserving_order = 1
        while not integral(sublattice * period * sublattice_inverse):
            period = period * action
            preserving_order += 1
        quotient_order = 1
        quotient_power = period
        while not integral((quotient_power - identity) * sublattice_inverse):
            quotient_power = quotient_power * period
            quotient_order += 1
        period_order = preserving_order * quotient_order

        family = list(multiples)
        for _ in range(1, period_order):
            multiples = [automorph(item) for item in multiples]
            family.extend(multiples)
        return tuple(family)

    @cached_method
    def binary_candidates(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        norm: Rational,
        residual: Rational,
    ) -> _Roots:
        # The pinned fixed-norm factorization admits negative residual
        # squares as well as positive ones. Only the zero-square case lacks
        # a finite signed factor-pair enumeration.
        if residual == 0:
            return ()
        edge = self.edge_plane(roots, discarded, k, norm)
        plane, scale, r0 = edge.plane, edge.scale, edge.r0
        if plane.is_isotropic():
            square = scale * residual
            if square.denominator() != 1:
                return ()
            # Pinned edgewalk.h:921-961 uses the *entire* signed factor-pair
            # enumeration, including imprimitive vectors, before ordering
            # roots and checking their ambient integral reconstructions.
            # A vector-orbit representative family is not that candidate set.
            vectors = [item for item in binary_split_vectors_of_square(plane, int(square)) if _oriented_determinant(_key(r0), _key(item)) > 0]
            if self.norm(k) < 0:
                vectors = [item for item in vectors if plane.b(r0, item) > 0]
            vectors.sort(key=lambda item: cmp_to_key(_orientation_compare)(_key(item)))
            return tuple(vectors)
        target_square = _rational(scale * residual)
        if target_square <= 0 or target_square > _rational(scale * norm) or target_square.denominator() != 1:
            return ()
        return tuple(item for item in self._anisotropic_edge_candidate_family(roots, discarded, k, norm) if _rational(item.q()) == target_square)

    def extension_roots(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        extension: NormedDynkinExtension,
    ) -> _Roots:
        edge = self.edge_plane(roots, discarded, k, extension.norm)
        work_vectors = self.binary_candidates(roots, discarded, k, extension.norm, extension.residual_norm)
        component = vector(SageQQ, extension.u_component)
        for work_vector in work_vectors:
            ambient = component + edge.basis.transpose() * vector(SageQQ, _key(work_vector))
            if any(entry.denominator() != 1 for entry in ambient):
                continue
            # edgewalk.h:905,988 get_next_anisotropic/get_next_isotropic:
            # SolutionIntMat(Latt, v_T), every integral candidate must
            # belong to the norm-specific lattice.
            if any(entry.denominator() != 1 for entry in ambient * edge.constraint_inverse):
                continue
            root = self.element(tuple(ambient))
            if self.pair(ambient, ambient) != extension.norm:
                raise ArithmeticError("rank-two reconstruction changed the root norm")
            # edgewalk.h get_next_anisotropic/get_next_isotropic returns the
            # FIRST integral member of this ordered family. Later members
            # are not inspected, even if vertex_from_root subsequently
            # rejects the resulting endpoint. Preserve that stopping rule
            # rather than reconstructing the entire family for each Dynkin
            # extension (the family itself is cached by binary_candidates).
            return (root,)
        return ()

    def candidate_compare(self, left: _Candidate, right: _Candidate, k: _Vector) -> int:
        left_scal = -_rational(self.pairing(left.alpha, k))
        right_scal = -_rational(self.pairing(right.alpha, k))
        left_sign = (left_scal > 0) - (left_scal < 0)
        right_sign = (right_scal > 0) - (right_scal < 0)
        first = _signed_sqrt_compare(
            left_sign,
            left_scal * left_scal / left.residual_norm,
            right_sign,
            right_scal * right_scal / right.residual_norm,
        )
        if first:
            return first
        second = _signed_sqrt_compare(
            left_sign,
            left_scal * left_scal / left.root_norm,
            right_sign,
            right_scal * right_scal / right.root_norm,
        )
        if second:
            return second
        if left.root_norm == right.root_norm:
            return 0
        return -1 if left.root_norm < right.root_norm else 1

    def vertex_from_root(self, k: _Vector, roots: _Roots, discarded: _Vector, alpha: _Vector) -> _Vertex | None:
        extended = (*roots, alpha)
        kernel = self.perpendicular_rows(tuple(self.row(root) for root in extended))
        if kernel.nrows() != 1:
            raise ArithmeticError("an edge-root extension must determine a one-dimensional perpendicular space")
        generator = self.primitive(kernel.row(0))
        old_pairing = self.pairing(generator, k)
        if old_pairing > 0:
            generator = -generator
            old_pairing = -old_pairing
        if old_pairing == 0 or self.pairing(discarded, generator) >= 0:
            return None
        # edgewalk.h:1032-1058 stores alpha itself in MatRoot. A scalar
        # multiple defines the same hyperplane, but not the same normed root
        # configuration; replacing alpha by its primitive part here changes
        # the reference algorithm's subsequent extensions and stabilizers.
        return _Vertex(generator, (*roots, alpha))

    def cusp_root(self, constraint: Parent, component: tuple[Rational, ...], k: _Vector) -> _Vector | None:
        r"""The point u + c k of C on the line through the component u in direction k with the least c > 0."""
        points = constraint.affine_line_points(self.rational_element(vector(SageQQ, component)), k)
        if points is None:
            return None
        base, step = (constraint.inclusion()(point) for point in points)
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
            scalar = -_rational(self.pairing(previous, root))
            sign = (scalar > 0) - (scalar < 0)
            candidates.append((sign, scalar * scalar / extension.norm, extension.norm, root))

        def compare(
            left: tuple[int, Rational, Rational, _Vector],
            right: tuple[int, Rational, Rational, _Vector],
        ) -> int:
            first = _signed_sqrt_compare(left[0], left[1], right[0], right[1])
            if first:
                return first
            if left[2] == right[2]:
                return 0
            return -1 if left[2] < right[2] else 1

        candidates.sort(key=cmp_to_key(compare))
        accepted = list(roots)
        for _sign, _quant, _norm, root in candidates:
            if all(self.pairing(old, root) <= 0 for old in accepted):
                # edgewalk.h:393-408, DetermineRootsCuspidalCase, appends the
                # actual chosen vector, without primitive normalization.
                accepted.append(root)
        return tuple(accepted)

    def isotropic_next_vertex(self, roots: _Roots, discarded: _Vector, k: _Vector) -> _Vertex:
        perpendicular = self.edge_perpendicular(roots)
        plane = perpendicular.inclusion().domain()
        choices: dict[_Key, _Vector] = {}
        for isotropic in plane.primitive_isotropic_vectors():
            candidate = perpendicular.inclusion()(isotropic)
            # edgewalk.h:1094-1120, get_can_gen orients each isotropic
            # line by b(k, v) < 0 *before* testing the discarded wall.
            # Primitive-isotropic-line enumerators need not return the
            # representative in this time component; rejecting the wrong
            # sign loses the unique forward cusp altogether.
            direction = self.pairing(candidate, k)
            if direction == 0:
                # In the split edge plane this is the old isotropic line.
                continue
            if direction > 0:
                candidate = -candidate
            if self.pairing(discarded, candidate) < 0:
                choices.setdefault(_key(candidate), candidate)
        if len(choices) != 1:
            raise ArithmeticError("the edge must determine one forward isotropic endpoint")
        (generator,) = choices.values()
        return _Vertex(generator, self.cusp_roots(roots, generator, k))

    @cached_method
    def finite_extensions(self, roots: _Roots) -> tuple[NormedDynkinExtension, ...]:
        """Complete pinned Dynkin extensions of a fixed edge-wall system."""
        return finite_type_extensions(self.lattice, roots, self.norms)

    @cached_method
    def edge_step(self, vertex: _Vertex, roots: _Roots, discarded: _Vector) -> _Vertex:
        candidates: list[_Candidate] = []
        for extension in self.finite_extensions(roots):
            # edgewalk.h:966-1061 calls get_next(extension) once: it takes
            # the first integral candidate in the binary search ordering.
            # Only that candidate is tested as a possible next vertex.
            chosen = self.extension_roots(roots, discarded, vertex.generator, extension)
            if chosen:
                alpha = chosen[0]
                next_vertex = self.vertex_from_root(vertex.generator, roots, discarded, alpha)
                if next_vertex is not None:
                    candidates.append(_Candidate(alpha, extension.residual_norm, extension.norm, next_vertex))
        if candidates:
            candidates.sort(key=cmp_to_key(lambda left, right: self.candidate_compare(left, right, vertex.generator)))
            return candidates[0].vertex
        return self.isotropic_next_vertex(roots, discarded, vertex.generator)

    @cached_method
    def edge_directions(self, vertex: _Vertex) -> tuple[tuple[_Roots, _Vector], ...]:
        r"""The walls along each outgoing edge and a wall orienting that edge."""
        if len(vertex.roots) == 1:
            return ()
        cone = RationalPolyhedralCones(self.lattice).from_rays(vertex.roots)
        directions: list[tuple[_Roots, _Vector]] = []
        for covector in cone.facet_covectors():
            included = tuple(root for root in vertex.roots if cone.evaluate_covector(covector, root) == 0)
            omitted = tuple(root for root in vertex.roots if cone.evaluate_covector(covector, root) != 0)
            if not omitted:
                raise ArithmeticError("a root-cone facet must omit a wall")
            # An ideal corner with several affine components can omit more
            # than one wall. Any omitted wall fixes the same forward side.
            directions.append((included, omitted[0]))
        return tuple(directions)

    def coloured_family(self, pairs: tuple[tuple[_Vector, int], ...]) -> Element:
        r"""The coloured family i -> (v_i, c_i) in Sets().Mor(I, L x {0, 1})."""
        colours = Sets.Δ[1]
        points = Sets().product((self.lattice, colours))
        return Sets().Mor(Sets.Δ[len(pairs) - 1], points)(lambda index: points((pairs[int(index)][0], colours[pairs[int(index)][1]])))

    @cached_method
    def vertex_configuration(self, vertex: _Vertex) -> Element:
        r"""The coloured configuration of the vertex: its roots (colour 0) and the vertex (colour 1).

        # source: fund_domain_vertices.h:128-162 GetInitialComputation (map_v: roots 1, gen 2)

        At an ideal vertex the family spans the degenerate hyperplane k^perp,
        and its stabilizer and transporters are the ``isotropstabequiv``
        computations of fund_domain_vertices.h:466-838 (selected in
        edgewalk.h:1170-1199 gen_fund_domain_fund_info), which the
        coloured-family stabilizer of ``O(L)`` performs for such families.
        """
        vectors = [*vertex.roots, vertex.generator]
        colours = [0] * len(vertex.roots) + [1]
        return self.coloured_family(tuple(zip(vectors, colours, strict=True)))

    def isometry_key(self, isometry: LatticeIsometryMethods) -> tuple[_Key, ...]:
        return tuple(_key(isometry(generator)) for generator in self.lattice.module_generators())

    @staticmethod
    def common_fixed_space(invariant: Parent, isometry: LatticeIsometryMethods) -> Parent | None:
        r"""The sublattice fixed by the group generated so far, or None once that group is infinite.

        A finite cone-preserving group fixes a timelike vector by averaging.
        Conversely its integral stabilizer is discrete in a compact real
        orthogonal group, hence finite. The common fixed lattice decides both
        directions, without separate order computations for the generators.
        """
        fixed = invariant.intersection(isometry.invariant_lattice())
        _positive, negative = fixed.inclusion().domain().signature_pair()
        return None if int(negative) == 0 else fixed

    def root_chamber_certificate(self, roots: _Roots) -> bool:
        r"""Upstream LORENTZ_RunEdgewalkAlgorithm dual-description ray test.

        The inequalities are -b(root, x) >= 0.  A pointed rational
        chamber is tested on its exact extreme rays, as in the pinned
        upstream dual-description termination condition.
        """
        if not roots:
            return False
        inequalities = tuple(-root.to_covector() for root in roots)
        chamber = RationalPolyhedralCones(self.lattice).from_halfspaces(inequalities)
        if not chamber.is_pointed():
            return False
        rays = tuple(chamber.primitive_rays())
        # The dual-description norm check alone is vacuous on the zero
        # cone.  A hyperbolic chamber must have full span and contain a
        # timelike vector.  Moreover, different nonspacelike rays must
        # belong to the same closed time cone.  Otherwise their positive
        # span contains spacelike vectors despite each ray having q <= 0.
        # Pairwise nonpositive inner products make the entire positive span
        # nonspacelike by expanding the square of any positive combination.
        if not rays or self.span_rank(rays) != self.rank:
            return False
        if any(self.norm(ray) > 0 for ray in rays):
            return False
        if any(self.pairing(left, right) > 0 for index, left in enumerate(rays) for right in rays[index + 1 :]):
            return False
        return any(self.norm(ray) < 0 for ray in rays) or any(self.pairing(left, right) < 0 for index, left in enumerate(rays) for right in rays[index + 1 :])

    def report(
        self,
        vertices: list[_Vertex],
        reflective: bool,
        simple_roots: _Roots,
        isometries: list[LatticeIsometryMethods],
    ) -> AllcockEdgewalkReport:
        cones = RationalPolyhedralCones(self.lattice)
        return AllcockEdgewalkReport(
            self.lattice,
            reflective=reflective,
            simple_roots=simple_roots,
            vertices=tuple(cones.from_rays([vertex.generator], wall_roots=vertex.roots) for vertex in vertices),
            isometry_generators=tuple(isometries),
        )

    def edge_orbits(self, vertex: _Vertex, local_generators: tuple[LatticeIsometryMethods, ...]) -> tuple[tuple[_Roots, _Vector], ...]:
        """Outgoing edges modulo the supplied full vertex stabilizer."""
        directions = self.edge_directions(vertex)
        direction_positions = {frozenset(_key(root) for root in roots): index for index, (roots, _discarded) in enumerate(directions)}

        def act_on_direction(isometry: LatticeIsometryMethods, index: int) -> int:
            return direction_positions[frozenset(_key(isometry(root)) for root in directions[index][0])]

        local_group = self.lattice.O().subgroup(local_generators)
        direction_orbits = FiniteGSets(local_group)(tuple(range(len(directions))), act_on_direction).orbits()
        return tuple(directions[int(orbit.representative())] for orbit in direction_orbits)

    def binary_walk(self) -> AllcockEdgewalkReport:
        """The edgewalk of a binary lattice: a compact segment between two walls, or no finite chamber.

        # source: edgewalk.h:2147-2188 StandardEdgewalkAnalysis, with rank two.

        ``isotropic_vinberg`` starts at an isotropic vector when there is one
        (edgewalk.h:2119-2126). Its quotient k^perp/k is zero, so the
        projected root list of get_simple_cone is empty and the source raises
        NonReflectivityException "list of roots is empty" (edgewalk.h:1984-1989).
        Otherwise Vinberg supplies a vertex k on one wall r (GetOneVertex,
        get_simple_cone), its single facet is l_ui = [], v_disc = r, and
        EdgewalkProcedure finds the next wall. That vertex is inserted unless
        it is equivalent to k (edgewalk.h:1449-1537 func_insert_vertex), and
        its own edge returns to k.
        """
        if self.lattice.is_isotropic():
            return self.report([], False, (), [])
        initial = self.initial_vertex()
        (root,) = initial.roots
        neighbor = self.edge_step(initial, (), root)
        (neighbor_root,) = neighbor.roots
        symmetry = self.lattice.O().transporter(self.vertex_configuration(initial), self.vertex_configuration(neighbor))
        vertices = [initial] if symmetry is not None else [initial, neighbor]
        return self.report(vertices, True, (root, neighbor_root), [] if symmetry is None else [symmetry])

    def walk(self) -> AllcockEdgewalkReport:
        if self.rank == 2:
            return self.binary_walk()
        initial = self.initial_vertex()
        if initial is None:
            return self.report([], False, (), [])

        vertices = [initial]
        configurations = [self.vertex_configuration(initial)]
        seen: set[tuple[_Key, ...]] = {self.identity_key}
        isometries: list[LatticeIsometryMethods] = []
        invariant = self.lattice.subobject_on(list(self.lattice.module_generators()))

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

        def certified_roots() -> _Roots:
            # LORENTZ_RunEdgewalkAlgorithm::f_try_terminate: close the
            # accumulated walls under every discovered chamber symmetry.
            group = self.lattice.O().subgroup(tuple(isometries))
            seeds = tuple({_key(root): root for vertex in vertices for root in vertex.roots}.values())
            closure = FiniteGSets(group).from_finite_orbits(seeds, lambda isometry, root: isometry(root))
            return tuple(sorted(closure, key=_key))

        certified_root_count = 0

        def try_reflectivity_certificate() -> AllcockEdgewalkReport | None:
            # LORENTZ_RunEdgewalkAlgorithm::f_try_terminate tests the dual
            # description only when the closed root orbit has grown.
            nonlocal certified_root_count
            roots = certified_roots()
            if len(roots) <= certified_root_count:
                return None
            certified_root_count = len(roots)
            if self.root_chamber_certificate(roots):
                return self.report(vertices, True, roots, isometries)
            return None

        # The source checks after insertion of the initial vertex and after
        # every newly discovered vertex, not only once the quotient closes.
        certificate = try_reflectivity_certificate()
        if certificate is not None:
            return certificate

        position = 0
        while position < len(vertices):
            vertex = vertices[position]
            configuration = configurations[position]
            local_generators = tuple(self.lattice.O().stabilizer(configuration).group_generators())
            for automorphism in local_generators:
                previous_generators = len(isometries)
                if not admit(automorphism):
                    return self.report(vertices, False, (), isometries)
                if len(isometries) > previous_generators:
                    certificate = try_reflectivity_certificate()
                    if certificate is not None:
                        return certificate
            for roots, discarded in self.edge_orbits(vertex, local_generators):
                neighbour = self.edge_step(vertex, roots, discarded)
                neighbour_configuration = self.vertex_configuration(neighbour)
                isometries_to_neighbour = (self.lattice.O().transporter(known, neighbour_configuration) for known in configurations)
                transporter = next(
                    (isometry for isometry in isometries_to_neighbour if isometry is not None),
                    None,
                )
                if transporter is None:
                    vertices.append(neighbour)
                    configurations.append(neighbour_configuration)
                    certificate = try_reflectivity_certificate()
                    if certificate is not None:
                        return certificate
                else:
                    previous_generators = len(isometries)
                    if not admit(transporter):
                        return self.report(vertices, False, (), isometries)
                    if len(isometries) > previous_generators:
                        certificate = try_reflectivity_certificate()
                        if certificate is not None:
                            return certificate
            position += 1

        simple_roots = certified_roots()
        return self.report(vertices, self.root_chamber_certificate(simple_roots), simple_roots, isometries)


__all__ = ["edgewalk_fundamental_domain"]
