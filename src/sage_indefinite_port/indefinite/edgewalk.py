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

from collections.abc import Callable

from dataclasses import dataclass
from functools import cache, cmp_to_key

from sage_extensions.all import ZZ, HyperbolicLattices, Lattices
from sage_extensions.categories.definite_lattices import _ExactCVPEngine
from sage_extensions.categories.hyperbolic_lattices import AllcockEdgewalkReport
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage_extensions.categories.polyhedral_cones import RationalPolyhedralCones, _primitive_integral_coordinates
from sage_extensions.categories.sets.set_categories import Sets
from sage_extensions.tensors.tensor import _engine_component_matrix
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.libs.pari import pari
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
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
    binary_split_form_vectors_of_square,
    binary_split_null_basis,
)

type _Vector = Lattices.ElementMethods
type _Roots = tuple[_Vector, ...]
type _Key = tuple[int, ...]
type _RationalRow = FreeModuleElement[Rational]


def edgewalk_fundamental_domain(
    lattice: HyperbolicLattices.ParentMethods,
    *,
    option_norms: str = "all",
    early_termination: bool = True,
    initial_vertex: _Vector | tuple[_Vector, _Roots] | None = None,
) -> AllcockEdgewalkReport:
    r"""Allcock's edgewalk on a hyperbolic lattice.

    Returns the preamble's ``AllcockEdgewalkReport`` of ``lattice``: the simple roots of a
    fundamental polyhedron of the reflection group W(L), as elements of ``lattice``, and whether
    ``lattice`` is reflective.

    edgewalk.h:55-67, :2100-2126. ``option_norms`` is ``"all"`` or ``"K3"``
    (norm-two roots only). With ``early_termination`` false the walk continues
    past an infinite polyhedron symmetry group and lists every vertex orbit.
    ``initial_vertex`` is a vertex, or a vertex with its simple roots.
    """
    if initial_vertex is not None:
        # A prescribed corner belongs to the original presentation and must
        # be used as the start, including in the special rank-two cycle.
        return _Edgewalk(lattice, option_norms=option_norms,
                         early_termination=early_termination,
                         initial_vertex=initial_vertex).walk()
    if int(lattice.rank()) > 2 and not lattice.is_isotropic():
        from sage_indefinite_port.indefinite.indefinite_reduction import (
            indefinite_reduced_presentation,
        )

        reduction = indefinite_reduced_presentation(lattice)
        result = _Edgewalk(reduction.reduced, option_norms=option_norms,
                           early_termination=early_termination, initial_vertex=initial_vertex).walk()
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
            reason_non_reflective=result.reason_non_reflective(),
        )
    return _Edgewalk(lattice, option_norms=option_norms,
                     early_termination=early_termination, initial_vertex=initial_vertex).walk()


def edgewalk_isometry(
    source: HyperbolicLattices.ParentMethods,
    target: HyperbolicLattices.ParentMethods,
) -> LatticeIsometryMethods | None:
    r"""Exact vertex-configuration equivalence along Allcock's edgewalk.

    The source of ``LORENTZ_RunEdgewalkAlgorithm_Isomorphism`` traverses
    vertex orbits of the first form and compares each coloured vertex
    configuration with the initial vertex of the second form. Its positive
    answer is a certified integral ambient isometry. Use the same finite
    vertex-orbit traversal and the owned configuration transporter.
    """
    if (int(source.rank()) != int(target.rank())
            or source.signature_pair() != target.signature_pair()):
        return None
    if source is target:
        return source.O().identity()
    if (_engine_component_matrix(source.gram_tensor()).det()
            != _engine_component_matrix(target.gram_tensor()).det()):
        return None
    if source.is_even() != target.is_even():
        return None
    if int(source.rank()) == 2:
        # The rank-two edgewalk has a binary cycle rather than the
        # higher-dimensional finite vertex-orbit graph. A chosen initial
        # wall need not match one from a different integral presentation.
        # The exact binary Lorentzian isometry owner handles both the
        # isotropic and anisotropic cases, including the cycle period.
        from sage_indefinite_port.indefinite.lorentzian_cells import lorentzian_isometry

        return lorentzian_isometry(source, target)
    first, second = _Edgewalk(source), _Edgewalk(target)
    first_vertex, second_vertex = first.initial_vertex(), second.initial_vertex()
    if first_vertex is None or second_vertex is None:
        # Absence of an edgewalk starting corner is not an isometry
        # obstruction. The complete Lorentzian perfect-cell equivalence
        # algorithm remains available for this non-vertex branch.
        from sage_indefinite_port.indefinite.lorentzian_cells import lorentzian_isometry

        return lorentzian_isometry(source, target)
    goal = second.vertex_configuration(second_vertex)
    seen = []
    frontier = [first_vertex]
    incomplete_vertex_graph = False
    while frontier:
        vertex = frontier.pop(0)
        configuration = first.vertex_configuration(vertex)
        witness = source.Isom(target).transporter(configuration, goal)
        if witness is not None:
            if not all(target.b(witness(x), witness(y)) == source.b(x, y)
                       for x in _coordinate_basis(source) for y in _coordinate_basis(source)):
                raise ArithmeticError("the edgewalk vertex transporter is not an ambient isometry")
            return witness
        seen.append(configuration)
        local_generators = tuple(source.O().stabilizer(configuration).group_generators())
        for roots, discarded in first.edge_orbits(vertex, local_generators):
            neighbour = first.edge_step(vertex, roots, discarded)
            if neighbour is None:
                # A nonreflective edge has no next chamber vertex. The
                # vertex-orbit graph is then not a complete isometry
                # obstruction: use the complete perfect-cell equivalence
                # algorithm rather than returning a false negative.
                incomplete_vertex_graph = True
                continue
            next_configuration = first.vertex_configuration(neighbour)
            if any(source.O().transporter(previous, next_configuration) is not None
                   for previous in (*seen, *(first.vertex_configuration(item) for item in frontier))):
                continue
            frontier.append(neighbour)
    if incomplete_vertex_graph:
        from sage_indefinite_port.indefinite.lorentzian_cells import lorentzian_isometry

        return lorentzian_isometry(source, target)
    return None


def _key(vector_: _Vector) -> _Key:
    coordinates = vector_.to_vector()
    return tuple(_integer(coordinates(label)) for label in vector_.parent().module_generating_set())


def _coordinate_basis(lattice: Lattices.ParentMethods) -> tuple[_Vector, ...]:
    """Basis elements in the exact order used by coordinate-label matrices."""
    return tuple(lattice.module_generator(label) for label in lattice.module_generating_set())


def _isometry_row_matrix_in_coordinate_order(lattice: Lattices.ParentMethods, isometry: LatticeIsometryMethods) -> Matrix_rational_dense:
    """Rows of image coordinates in the lattice's ordered coordinate labels.

    The binary edge-period code compares this matrix with sublattice
    matrices in coordinate-label order. Do not use ``module_generators()``
    to choose the rows: its enumeration is a separate API from the labels.
    """
    labels = tuple(lattice.module_generating_set())
    basis = _coordinate_basis(lattice)
    return matrix(
        SageQQ,
        [tuple(_rational(isometry(generator).to_vector()(coordinate)) for coordinate in labels) for generator in basis],
    )


def _hermite_rows(rows: list[list[int]]) -> tuple[list[list[int]], list[list[int]]]:
    r"""The nonzero rows of the row Hermite normal form H = U A of integer rows A, and the rows of U spanning the left kernel of A.

    Cohen, *A Course in Computational Algebraic Number Theory*, Algorithm 2.4.4
    (Hermite normal form by Euclidean row operations, transposed to rows) and
    Proposition 2.4.9: with U unimodular and U A = [H; 0], the rows of U beside
    the zero rows of U A are a basis of the integral left kernel of A. Pivots are
    positive and the entries above a pivot lie in [0, pivot), Sage's convention
    for ``hermite_form``; the form of a lattice is unique, so this is the basis
    Sage computes.
    """
    a = [list(row) for row in rows]
    u = [[int(i == j) for j in range(len(a))] for i in range(len(a))]
    width = len(a[0]) if a else 0
    pivot_row = 0
    for column in range(width):
        while True:
            live = [i for i in range(pivot_row, len(a)) if a[i][column]]
            if not live:
                break
            least = min(live, key=lambda i: abs(a[i][column]))
            a[pivot_row], a[least] = a[least], a[pivot_row]
            u[pivot_row], u[least] = u[least], u[pivot_row]
            if len(live) == 1:
                break
            pivot = a[pivot_row][column]
            for i in range(pivot_row + 1, len(a)):
                if a[i][column]:
                    q = a[i][column] // pivot
                    a[i] = [x - q * y for x, y in zip(a[i], a[pivot_row])]
                    u[i] = [x - q * y for x, y in zip(u[i], u[pivot_row])]
        if pivot_row == len(a) or not a[pivot_row][column]:
            continue
        if a[pivot_row][column] < 0:
            a[pivot_row] = [-x for x in a[pivot_row]]
            u[pivot_row] = [-x for x in u[pivot_row]]
        pivot = a[pivot_row][column]
        for i in range(pivot_row):
            q = a[i][column] // pivot
            if q:
                a[i] = [x - q * y for x, y in zip(a[i], a[pivot_row])]
                u[i] = [x - q * y for x, y in zip(u[i], u[pivot_row])]
        pivot_row += 1
    return a[:pivot_row], u[pivot_row:]


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


def _cone_of_rays(rank: int, rows) -> Polyhedron:
    r"""The cone spanned by coordinate rows, as the engine polyhedron of ``RationalPolyhedralCones.from_rays``.

    The owned cone rebuilds the rays' polyhedron from the primitive integral
    rows of its H-representation, which fixes the order of its facets and
    its extreme rays; this is the same engine polyhedron, without the owned
    covectors and vectors.
    """
    rays = Polyhedron(vertices=[(0,) * rank], rays=rows, base_ring=SageQQ)
    return Polyhedron(
        ambient_dim=rank,
        ieqs=[[0, *_primitive_integral_coordinates(tuple(row.A()))] for row in rays.inequalities()],
        eqns=[[0, *_primitive_integral_coordinates(tuple(row.A()))] for row in rays.equations()],
        base_ring=SageQQ,
    )


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

    ``form`` is the Gram coefficients (a, b, c) of Pi with its form scaled by ``scale`` to be integral,
    in the oriented basis ``basis`` (rational rows of L tensor Q); ``r0`` is the coordinate pair of
    the edge's primitive orientation vector in Pi (edgewalk.h:747, r0_work).
    ``basis`` holds integer rows h_0, h_1 with b_i = h_i / ``basis_denominator``.
    ``constraint_inverse`` is (M, e) with M / e the inverse of the rows of a basis of C, so an
    integral row v lies in C exactly when e divides v * M (edgewalk.h:905, SolutionIntMat(Latt, v_T)).
    ``intersection`` holds the rows, in the coordinates of ``basis``, of a basis of C cap P, the
    sublattice of Pi whose quotient controls the period (edgewalk.h:803, Expr_t).
    """

    constraint_inverse: tuple[tuple[tuple[int, ...], ...], int]
    basis: tuple[_Key, _Key]
    basis_denominator: int
    form: tuple[int, int, int]
    scale: Integer
    r0: _Key
    intersection: Matrix_rational_dense

    @property
    def plane(self) -> Parent:
        r"""The exact integral rank-two form represented by the cached coefficients."""
        a, b, c = self.form
        return Lattices(ZZ)([[a, b], [b, c]])


class _Edgewalk:
    def __init__(self, lattice: HyperbolicLattices.ParentMethods, *,
                 option_norms: str = "all", early_termination: bool = True,
                 initial_vertex: _Vector | tuple[_Vector, _Roots] | None = None) -> None:
        positive, negative = lattice.signature_pair()
        if int(negative) != 1 or int(positive) + 1 != int(lattice.rank()):
            raise ValueError(f"Allcock's edgewalk needs a nondegenerate lattice of signature (n, 1), but {lattice} has signature {lattice.signature_pair()}")
        self.lattice = lattice
        self.rank = int(lattice.rank())
        self.labels = tuple(lattice.module_generating_set())
        coordinate_basis = _coordinate_basis(lattice)
        self.gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageQQ)
        self.rational = lattice.vector_space()
        self.field = lattice.base_ring().fraction_field()
        if option_norms not in ("all", "K3"):
            raise ValueError(f"unsupported edgewalk root-norm option {option_norms!r}")
        self.norms = ((2,) if option_norms == "K3"
                      else tuple(int(norm) for norm in lattice.possible_root_lengths()))
        self.early_termination = bool(early_termination)
        self.supplied_initial_vertex = initial_vertex
        # edgewalk.h:1985-1999 raise NonReflectivityException from
        # get_simple_cone; StandardEdgewalkAnalysis (edgewalk.h:2175-2183)
        # reports the reason. A starting corner that certifies nonreflectivity
        # therefore has to carry that reason to the report.
        self.initial_nonreflectivity_reason: str | None = None
        self.identity_key = tuple(_key(generator) for generator in coordinate_basis)
        # The congruence matrix and its homogeneous integral kernel in
        # Get_DataMapping depend on the root norm, not the affine coset.
        # They are reused for every height of one Vinberg search.
        self._vinberg_mapping_perpendicular = None
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

    @cached_method
    def key(self, vector_: _Vector) -> _Key:
        r"""The integral coordinate tuple of a lattice vector, read once per vector."""
        return _key(vector_)

    @cached_method
    def row(self, vector_: _Vector) -> _RationalRow:
        r"""The immutable rational coordinate row of a lattice vector, read once per vector."""
        coordinates = vector(SageQQ, self.key(vector_))
        coordinates.set_immutable()
        return coordinates

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
        r"""The echelon basis of the rational perpendicular of ``vectors``.

        Sage's generic Gaussian elimination returns the same unique echelon
        kernel basis as its default route through FLINT, without the
        conversions that cost ten times as much on these small matrices.
        """
        return (
            matrix(
                SageQQ,
                len(vectors),
                self.rank,
                [entry for row in vectors for entry in row],
            )
            * self.gram
        ).right_kernel_matrix(algorithm="generic")

    def initial_isotropic_vertex(self, k: _Vector | None = None) -> _Vertex | None:
        r"""The ``isotropic_vinberg`` starting corner of pinned edgewalk.h:1880–2055.

        The source finds a primitive null vector, enumerates each admissible
        root-norm shell in the projection of C_N ∩ k-perp to k-perp/k,
        selects one simple cone, lifts its walls and completes the affine
        cusp roots. The quotient k^perp/k, its norm shells (PARI qfminim) and
        the simple cone (pplpy) are computed on engine data.
        """
        if k is None:
            k = self.lattice.isotropic_vector()
        # k^perp/k on engine data (edgewalk.h:1955-1968, LatticeProjectionFramework):
        # the rows K of an integral basis of k^perp (the integral left kernel
        # of G k), and the Smith transform U of the column c^T of the
        # coordinates of k in K, U c^T = (+-1, 0, ..., 0)^T, as Sage frames
        # the cokernel of k -> k^perp (the owned isotropic reduction's Smith
        # generators).  A vector x K of k^perp projects to the last
        # coordinates of U x^T; the columns of U^-1 after the first lift
        # the quotient generators.
        k_column = self.gram * self.row(k).column()
        perpendicular_rows = (k_column.denominator() * k_column).change_ring(SageZZ).left_kernel_matrix(algorithm="pari").echelon_form().change_ring(SageQQ)
        k_coordinates = perpendicular_rows.solve_left(self.row(k))
        smith, unimodular, _right = matrix(SageZZ, [k_coordinates]).transpose().smith_form()
        if smith[0, 0] != 1:
            raise ArithmeticError("an isotropic vector must be primitive in k^perp")
        rank = int(perpendicular_rows.nrows()) - 1
        quotient_rows = unimodular.matrix_from_rows(range(1, rank + 1)).transpose().change_ring(SageQQ)
        complement = unimodular.inverse().matrix_from_columns(range(1, rank + 1)).transpose().change_ring(SageQQ) * perpendicular_rows
        quotient_gram = complement * self.gram * complement.transpose()
        projected_roots: list[tuple[_RationalRow, _Vector]] = []

        for norm in self.norms:
            # edgewalk.h:1951-1954: Latt_inter_NSP = IntersectionLattice(Latt, NSP),
            # the integral left kernel of B G k on the rows B of a basis of C_N.
            constraint_rows, _inverse = self.constraint_rows(self.divisor(norm))
            equations = constraint_rows * k_column
            denominator = SageZZ.one()
            for entry in equations.list():
                denominator = denominator.lcm(entry.denominator())
            intersection_rows = (denominator * equations).change_ring(SageZZ).left_kernel_matrix(algorithm="pari").echelon_form().change_ring(SageQQ) * constraint_rows
            # fr.BasisProj (edgewalk.h:1955-1960): the projected generators,
            # read through the perpendicular basis and the quotient framing.
            projected_rows = perpendicular_rows.solve_left(intersection_rows) * quotient_rows
            if any(entry.denominator() != 1 for entry in projected_rows.list()):
                raise ArithmeticError("the projection of C_N cap k^perp must be integral in k^perp/k")
            hermite, transformation = projected_rows.change_ring(SageZZ).echelon_form(transformation=True)
            image_rank = hermite.rank()
            image_basis = hermite.matrix_from_rows(range(image_rank)).change_ring(SageQQ)
            lift_rows = transformation.matrix_from_rows(range(image_rank)).change_ring(SageQQ) * intersection_rows
            # The norm shell of the image lattice (edgewalk.h:1969-2000) on its
            # Gram matrix in the basis image_basis: PARI qfminim, the vectors
            # of square exactly N, each followed by its negative, as the
            # definite shell enumeration orders them.
            image_gram = (image_basis * quotient_gram * image_basis.transpose()).change_ring(SageZZ)
            if not image_gram.is_positive_definite():
                raise ArithmeticError("the projected root lattice k^perp/k must be positive definite")
            target = _integer(norm)
            shell = []
            for column in pari(image_gram).qfminim(target, None, 0)[2].sage().columns():
                if column * image_gram * column == target:
                    shell.extend((column, -column))
            for coefficients in shell:
                projected = coefficients * image_basis
                # fr.GetOnePreimage (edgewalk.h:2010): one preimage in C_N cap k^perp.
                lifted = self.element(tuple(coefficients * lift_rows))
                if self.norm(lifted) != norm:
                    raise ArithmeticError("the isotropic-reduction lift changed a prescribed root norm")
                projected_roots.append((projected, lifted))

        if not projected_roots:
            # get_simple_cone raises NonReflectivityException here
            # (edgewalk.h:1985-1988) and StandardEdgewalkAnalysis turns it
            # into the report's reason (edgewalk.h:2175-2183): "The list of
            # vectors is empty. Cannot be reflective."
            self.initial_nonreflectivity_reason = "list of roots is empty"
            return None
        projected_rows = tuple(projected for projected, _root in projected_roots)
        if matrix(SageQQ, projected_rows).rank() != rank:
            # edgewalk.h:1990-1999: rnk < G.rows() - 2 means the roots at the
            # cusp cannot be a reflective cusp's simple system.
            self.initial_nonreflectivity_reason = "the rank of the roots is too low"
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
        positive = [(projected, root, row) for (projected, root), row in zip(projected_roots, projected_rows, strict=True) if splitting.dot_product(row) > 0]
        # The extreme rays of the cone of the positive projected roots in
        # k^perp/k (pplpy dual description, as the owned cone's engine).
        cone = Polyhedron(rays=[row for _projected, _root, row in positive], base_ring=SageQQ, backend="ppl")
        simple_roots: list[_Vector] = []
        simple_projected_rows: list[_RationalRow] = []
        for ray in cone.rays():
            ray_row = vector(SageQQ, ray.vector())
            # get_one_root (edgewalk.h:2001-2022): the admissible root on this
            # ray, i.e. a positive multiple of it (the quotient is definite).
            for _projected, root, root_row in positive:
                if matrix(SageQQ, [ray_row, root_row]).rank() == 1 and ray_row.dot_product(root_row) > 0:
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

    def _vinberg_polytope_roots(self, shift, directions, positive_gram, center, radius, root_norm, walls) -> _Roots:
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
            covector = -(ambient_gram * vector(SageZZ, self.key(wall)))
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
                if any(candidate.dot_product(ambient_gram * vector(SageZZ, self.key(wall))) > 0 for wall in walls):
                    raise ArithmeticError("Vinberg's constrained facet witness is invalid")
                roots.append(self.element(tuple(candidate)))

        enumerate_level(dimension - 1, minors[-1] * radius)
        return tuple(roots)

    def _vinberg_decomposed_roots(self, perpendicular, component: _Vector, root_norm: int, walls=None) -> _Roots:
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
        a = vector(SageZZ, self.key(component))
        if self._vinberg_mapping_perpendicular is not perpendicular:
            self._vinberg_mapping_perpendicular = perpendicular
            self._vinberg_mapping_cache.clear()
            self._vinberg_polytope_cache.clear()
        prepared = self._vinberg_mapping_cache.get(root_norm)
        if prepared is None:
            inclusion = perpendicular.inclusion()
            perpendicular = inclusion.domain()
            perpendicular_basis = _coordinate_basis(perpendicular)
            orthogonal = matrix(SageZZ, [self.key(inclusion(generator)) for generator in perpendicular_basis]).transpose()
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
        if walls and matrix(SageZZ, [self.key(wall) for wall in walls]).rank() == n:
            return self._vinberg_polytope_roots(shift, directions, positive_gram, center, radius, root_norm, walls)
        roots = []
        for coordinates, square in backend.close_vector_coordinates(center, radius):
            if square != radius:
                continue
            row = shift + directions * vector(SageZZ, coordinates)
            if row.dot_product(gram * row) != root_norm:
                raise ArithmeticError("Vinberg's spherical solution has incorrect norm")
            if any((2 * pairing) % root_norm for pairing in gram * row):
                raise ArithmeticError("Vinberg's spherical solution is not a root")
            roots.append(self.element(tuple(row)))
        return tuple(roots)

    def prescribed_vertex(self, initial: _Vector | tuple[_Vector, _Roots]) -> _Vertex | None:
        r"""The vertex given by edgewalk.h:2100-2115 FileVertexRoots or FileVertex.

        A bare vertex has its simple cone computed, as get_simple_cone does
        (edgewalk.h:1848-2049): a null ray goes through the cusp construction
        and a timelike one through the corner root system. When that
        construction certifies nonreflectivity (edgewalk.h:1985-1999) the
        caller must report the source's reason instead of raising
        (StandardEdgewalkAnalysis catches it at edgewalk.h:2175-2183), so None
        is returned and ``initial_nonreflectivity_reason`` carries the text.
        A prescribed ray is a ray: its primitive part is the same vertex, and
        RemoveFractionVector at edgewalk.h:2107/2143 is the source's own sign
        that the caller's scaling is not part of the vertex.
        """
        if isinstance(initial, tuple):
            generator, roots = initial
        else:
            generator, roots = initial, None
        if generator.parent() is not self.lattice or self.norm(generator) > 0:
            raise ValueError("the starting vertex must be a nonspacelike vector of the walking lattice")
        if not any(self.key(generator)):
            raise ValueError("the zero vector does not prescribe a vertex ray")
        # Both FileVertexRoots (:2107) and the FileVertex/common return
        # (:2143) apply RemoveFractionVector. The roots of FileVertexRoots are
        # data and are not primitive-normalised with the generator.
        ray = self.primitive(self.row(generator))
        if roots is not None:
            return _Vertex(ray, tuple(roots))
        if self.norm(ray) == 0:
            vertex = self.initial_isotropic_vertex(ray)
            if vertex is None:
                return None
            return vertex
        simple_roots = self.ordinary_vertex_roots(ray)
        if self.span_rank(simple_roots) != self.rank - 1:
            raise ValueError(f"{generator} is not a vertex: its perpendicular roots have rank {self.span_rank(simple_roots)}")
        return _Vertex(ray, simple_roots)

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
        # stopping on the first feasible nonpositive ray.
        #
        # GetVinbergAux (vinberg.h:213-307) takes the saturated integral
        # perpendicular and ComputeTranslationClasses for
        # L / (v0^perp + Z*v0).  Its pinned basic_common_cpp@995fd625
        # MAT_MatrixInt.h:2881-2958 traverses the classes breadth first,
        # starting at zero and adding the ambient basis vectors in order.
        # Preserve that order through the actual module cokernel: the root
        # iterator starts at index 1 and consumes W sequentially.  Neither
        # ComputeTranslationClasses nor GetVinbergAux sorts by height.
        perpendicular = timelike.orthogonal_complement()
        perpendicular_inclusion = perpendicular.inclusion()
        perpendicular_generators = tuple(perpendicular_inclusion(generator) for generator in _coordinate_basis(perpendicular_inclusion.domain()))
        split = self.lattice.subobject_on((*perpendicular_generators, timelike))
        projection = split.inclusion().cokernel_projection()
        quotient = projection.codomain()
        representatives = [self.lattice.zero()]
        residues = [quotient.zero()]
        generators = _coordinate_basis(self.lattice)
        next_representative = 0
        while next_representative < len(representatives):
            representative = representatives[next_representative]
            for generator in generators:
                candidate = representative + generator
                residue = projection(candidate)
                if residue not in residues:
                    residues.append(residue)
                    representatives.append(candidate)
            next_representative += 1
        norm_shift = -self.norm(timelike)
        cosets = []
        for representative in representatives:
            # GetVinbergAux uses integer QuoInt (vinberg.h:258-260).
            # Division in QQ computes the field quotient, including for //.
            multiple = _integer(-self.pairing(representative, timelike)) // _integer(norm_shift)
            normalized = representative - multiple * timelike
            # GetVinbergAux uses QuoInt/ResInt to place each class in
            # the fundamental timelike strip, 0 <= -<w,v0> < -v0^2.
            height = -self.pairing(normalized, timelike)
            if not 0 <= height < norm_shift:
                raise ArithmeticError("Vinberg coset representative leaves timelike strip")
            cosets.append(normalized)
        # vinberg.h FundCone -> GetFacetOneDomain: select one member of
        # each perpendicular root pair by a regular linear functional,
        # and retain the original normed roots on the extreme cone rays.
        # Lexicographic sign is a deterministic regular functional for this
        # finite integral root family.
        initial_roots = [root for norm in self.norms for root in self._vinberg_decomposed_roots(perpendicular, self.lattice.zero(), norm)]
        # GetFacetOneDomain first reduces the root coordinate columns to
        # their actual span, and GetFacetOneDomain_ListIdx requires at least
        # twice that reduced dimension in the nonempty root family.
        if initial_roots:
            reduced_dimension = self.span_rank(tuple(initial_roots))
            if len(initial_roots) < 2 * reduced_dimension:
                raise ArithmeticError("initial Vinberg root family has fewer than twice its rank in roots")
        positive = [root for root in initial_roots if next((entry for entry in self.key(root) if entry), 0) > 0]
        # GetFacetOneDomain_ListIdx allocates exactly half the roots on
        # the assumption that the family is centrally symmetric. A count
        # alone cannot establish the required opposite-root correspondence.
        if initial_roots:
            root_rows = [self.key(root) for root in initial_roots]
            multiplicities = {}
            for row in root_rows:
                multiplicities[row] = multiplicities.get(row, 0) + 1
            if any(multiplicities.get(row, 0) != multiplicities.get(tuple(-entry for entry in row), 0) for row in multiplicities):
                raise ArithmeticError("Vinberg perpendicular roots do not form opposite pairs")
        if positive:
            positive_cone = RationalPolyhedralCones(self.lattice).from_rays(positive)
            rays = tuple(self.row(ray) for ray in positive_cone.primitive_rays())
            # GetFacetOneDomain_ListIdx selects nonredundant *indices*;
            # parallel root vectors cannot both survive its reduction.
            chosen = []
            used_rays = set()
            ray_by_direction = {_primitive_positive_direction(ray): index for index, ray in enumerate(rays)}
            for root in positive:
                index = ray_by_direction.get(_primitive_positive_direction(self.key(root)))
                if index is not None and index not in used_rays:
                    chosen.append(root)
                    used_rays.add(index)
            if len(used_rays) != len(rays):
                raise ArithmeticError("initial Vinberg cone ray has no normed root")
            accepted_roots = chosen
        else:
            accepted_roots = []

        def accepted_vertex() -> _Vertex | None:
            roots = tuple(accepted_roots)
            if not roots:
                return None
            equations = matrix(SageQQ, [self.row(root) * self.gram for root in roots])
            rank = equations.rank()
            if rank < self.rank - 1:
                return None
            if rank == self.rank - 1:
                generator = self.primitive(equations.right_kernel_matrix(algorithm="generic").row(0))
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
                inequalities = tuple(self.wall_covector(root) for root in roots)
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
        # belonging to one decomposition request has been inserted, so an
        # intermediate root cannot end the search.
        # vinberg.h IterateRootDecompositions::get_k (pinned 1592b246)
        # prioritizes each length's next coset by val/sqrt(val), with
        # val = -<v0, a>.
        candidates = {norm: 1 for norm in self.norms}
        coset_heights = tuple(-self.pairing(timelike, representative) for representative in cosets)
        period = len(coset_heights)
        if not period:
            raise ArithmeticError("Vinberg decomposition has no coset representatives")
        while True:
            norm = min(
                self.norms,
                key=lambda length: coset_heights[candidates[length] % period] + (candidates[length] // period) * norm_shift,
            )
            index = candidates[norm]
            component = cosets[index % period] + (index // period) * timelike
            candidates[norm] = index + 1
            # vinberg.h FindRoot_filter/fct_CVP retains every norm-k
            # candidate satisfying all accepted wall inequalities.
            new_roots = [
                root
                for root in self._vinberg_decomposed_roots(perpendicular, component, norm, walls=accepted_roots)
                if all(self.pairing(root, wall) <= 0 for wall in accepted_roots)
            ]
            if not new_roots:
                continue
            accepted_roots.extend(new_roots)
            # vinberg.h ReduceListRoot removes generators redundant in the
            # accepted root cone. Keep the original normed root vectors,
            # using primitive rays only to identify their directions.
            ambient_roots = tuple(accepted_roots)
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
            for root in accepted_roots:
                direction = _primitive_positive_direction(self.key(root))
                index = ray_by_direction.get(direction)
                if index is not None and index not in retained_directions:
                    reduced.append(root)
                    retained_directions.add(index)
            if len(retained_directions) != len(ray_rows):
                raise ArithmeticError("Vinberg cone ray has no accepted normed root")
            accepted_roots = reduced
            initial = accepted_vertex()
            if initial is not None:
                return initial

    def ordinary_vertex_roots(self, timelike: _Vector, *, require_full: bool = False) -> _Roots:
        """A simple system of all ambient reflective roots perpendicular to k.

        Lexicographic positivity is a regular ordering on this finite root
        system. Its positive-root cone has precisely the simple-root rays.
        """
        # edgewalk.h:1845-1866: one FindFixedDistVectors per norm on the
        # Gram matrix of the definite complement k^perp.  The complement's
        # basis rows and Gram matrix are its engine data; PARI qfminim (the
        # engine of the complement's vectors_of_square) returns one vector
        # of each sign pair, in the order that shell lists them, and only
        # the accepted roots are wrapped as lattice vectors.
        # The rows of k^perp: the echelon basis of the integral left kernel
        # of G k, the basis of the owned orthogonal complement.
        gram = self.gram.change_ring(SageZZ)
        basis = (gram * vector(SageZZ, self.key(timelike))).column().left_kernel_matrix(algorithm="pari").echelon_form()
        complement_gram = pari(basis * gram * basis.transpose())
        positive: list[_Vector] = []
        for norm in self.norms:
            for point in complement_gram.qfminim(norm, None, 0)[2].sage().columns():
                root = point * basis
                if root * gram * root != norm:
                    continue
                if next(entry for entry in root if entry) < 0:
                    root = -root
                if any((2 * entry) % norm for entry in gram * root):
                    continue
                positive.append(self.element(tuple(root)))
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
        rays = tuple(_primitive_integral_coordinates(tuple(ray.vector())) for ray in _cone_of_rays(self.rank, [self.row(root) for root in positive]).rays())
        if require_full and len(rays) != self.rank - 1:
            # edgewalk.h get_simple_cone_from_lattice requires exactly
            # dimSpace selected walls for an ordinary starting vertex.
            raise ArithmeticError("ordinary initial cone has incorrect simple-root count")
        roots_by_direction = {}
        for root in positive:
            roots_by_direction.setdefault(_primitive_positive_direction(_key(root)), root)
        return tuple(roots_by_direction[_primitive_positive_direction(ray)] for ray in rays)

    @staticmethod
    def divisor(norm: Rational) -> int:
        r"""The divisor d with C_N = {v in L : b(v, L) in d Z}: N/2 for even N, N for odd N."""
        integral_norm = _integer(norm)
        return integral_norm // 2 if integral_norm % 2 == 0 else integral_norm

    @cached_method
    def constraint_rows(self, divisor: int) -> tuple[Matrix_rational_dense, Matrix_rational_dense]:
        r"""The rows of a basis of C = {v in L : b(v, L) in d Z} and their inverse, computed once per norm.

        edgewalk.h:110-131 ComputeLattice_LN and ComputeSublatticeInfos:
        C = IntersectionLattice(I, d G^-1), the integral vectors of the
        lattice d L^# written in the framing of L.  Generic callees: the
        intersection is Sage's exact ZZ-module intersection (echelon basis)
        in place of the pinned IntersectionLattice, and the pinned
        SublatticeBasisReduction (LLL) is not applied.  Every consumer
        (the edge plane projection, C cap k^perp and the cusp root) reads
        only the integral span of these rows, so the basis does not matter.
        """
        # C = {v in Z^n : G v in d Z^n}.  The integral kernel of [G | -d I]
        # is {(v, w) : G v = d w}; G is nondegenerate, so w is determined
        # by v and the projection to v is an isomorphism onto C.  Its
        # Hermite form is the echelon basis of the ZZ-module intersection.
        gram = self.gram.change_ring(SageZZ)
        kernel = gram.augment(-divisor * identity_matrix(SageZZ, self.rank)).right_kernel_matrix(algorithm="pari")
        rows = kernel.matrix_from_columns(range(self.rank)).hermite_form(include_zero_rows=False).change_ring(SageQQ)
        if rows.nrows() != self.rank or rows.rank() != self.rank:
            raise ArithmeticError("C_N must be a full-rank sublattice")
        rows.set_immutable()
        inverse = rows.inverse()
        inverse.set_immutable()
        return rows, inverse

    @cached_method
    def constraint_engine(self, divisor: int) -> tuple[tuple[_Key, ...], tuple[tuple[tuple[int, ...], ...], int]]:
        r"""The integer rows of the basis of C and (M, e) with M / e their inverse, read once per norm from ``constraint_rows``."""
        rows, inverse = self.constraint_rows(divisor)
        denominator = inverse.denominator()
        return (
            tuple(tuple(int(entry) for entry in row) for row in rows.rows()),
            (tuple(tuple(int(entry) for entry in row) for row in (denominator * inverse).rows()), int(denominator)),
        )

    @cached_method
    def gram_entries(self) -> tuple[tuple[Rational, ...], ...]:
        r"""The rows of the Gram matrix as tuples of rationals, read once per walk."""
        return tuple(tuple(row) for row in self.gram.rows())

    def edge_plane(self, roots: _Roots, discarded: _Vector, k: _Vector, norm: Rational) -> _EdgePlane:
        return self.edge_plane_for_divisor(roots, discarded, k, self.divisor(norm))

    @cached_method
    def edge_orientation(self, roots: _Roots, discarded: _Vector, k: _Vector) -> tuple[Matrix_rational_dense, _RationalRow]:
        r"""The oriented basis of the edge plane and the orientation vector r0, computed once per edge (edgewalk.h:598-668)."""
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
            # A vector of positive square in the Lorentzian plane, on its
            # basis rows b_0, b_1 with Gram g: b_0 if g_00 > 0; else, for
            # g_00 < 0, b_1 - (g_01/g_00) b_0 of square det(g)/g_00 > 0; and for
            # g_00 = 0, b_1 + t b_0 of square g_11 + 2 t g_01 = 1.  Every such
            # vector with b(., k) < 0 lies in one component of the positive
            # cone, which k bounds, so the orientation of (positive, k) and
            # the isotropy of the plane do not depend on the choice.
            first, second = plane_basis.rows()
            g00, g01, g11 = self.pair(first, first), self.pair(first, second), self.pair(second, second)
            if g00 > 0:
                positive = first
            elif g00 < 0:
                positive = second - (g01 / g00) * first
            else:
                positive = second + ((1 - g11) / (2 * g01)) * first
            if self.pair(positive, positive) <= 0:
                raise ArithmeticError("the edge plane must contain a vector of positive square")
            if self.pair(positive, self.row(k)) > 0:
                positive = -positive
            oriented_basis = matrix(SageQQ, [positive, self.row(k)])
            r0 = -self.row(k)
        return oriented_basis, r0

    @cached_method
    def edge_is_isotropic(self, roots: _Roots, discarded: _Vector, k: _Vector) -> bool:
        r"""Whether the edge plane P is isotropic, decided once per edge (edgewalk.h:871-889).

        The source factors Gtest = Pplane G Pplane^T (GetIsotropicFactorization).
        P is Lorentzian of rank two, so it is isotropic exactly when -det(Gtest)
        is a square in QQ.
        """
        oriented_basis, _r0 = self.edge_orientation(roots, discarded, k)
        return (-(oriented_basis * self.gram * oriented_basis.transpose()).det()).is_square()

    @cached_method
    def edge_wall_data(self, roots: _Roots) -> tuple[tuple[_Key, ...], int, tuple[_Key, ...]]:
        r"""The projection of L (x) QQ onto the edge plane and the pairings G W^T with the edge walls, computed once per edge (edgewalk.h:598, Pplane).

        The walls through an edge need not be linearly independent.  The
        projection only depends on their rational span; invert the pairing
        on an independent basis of that span (the reduced row echelon rows),
        not the possibly singular Gram matrix of every wall.  The edge plane
        is Lorentzian, so its orthogonal complement must be positive definite
        and nondegenerate.  Row vectors x map to x * K / e for the integer
        rows K and denominator e returned; the wall pairings are returned as
        integer rows of a positive multiple of G W^T (only their integral
        left kernel is read).
        """
        root_rows = tuple(self.row(root) for root in roots)
        wall_matrix = matrix(SageQQ, len(root_rows), self.rank, [entry for row in root_rows for entry in row])
        echelon = wall_matrix.echelon_form()
        root_matrix = echelon.matrix_from_rows(range(echelon.rank()))
        root_gram = root_matrix * self.gram * root_matrix.transpose()
        if not root_gram.is_invertible():
            raise ArithmeticError("the edge-wall span must be nondegenerate")
        complement = identity_matrix(SageQQ, self.rank) - self.gram * root_matrix.transpose() * root_gram.inverse() * root_matrix
        wall_pairings = self.gram * wall_matrix.transpose()
        complement_denominator = complement.denominator()
        wall_pairings = wall_pairings.denominator() * wall_pairings
        return (
            tuple(tuple(int(entry) for entry in row) for row in (complement_denominator * complement).rows()),
            int(complement_denominator),
            tuple(tuple(int(entry) for entry in row) for row in wall_pairings.rows()),
        )

    @cached_method
    def edge_plane_for_divisor(self, roots: _Roots, discarded: _Vector, k: _Vector, divisor: int) -> _EdgePlane:
        r"""The rank-two lattice Pi of an edge for one norm, on integer rows (edgewalk.h:747-803)."""
        oriented_basis, r0 = self.edge_orientation(roots, discarded, k)
        constraint, constraint_inverse = self.constraint_engine(divisor)
        complement, denominator, wall_pairings = self.edge_wall_data(roots)
        columns = range(self.rank)
        # The ZZ span of the projected rows C K / e is (1/e) times the span of
        # the integer rows C K; the Hermite normal form is unique and scales
        # with a positive factor, so this is the echelon basis of Sage's ZZ span.
        hermite, _kernel = _hermite_rows([[sum(c * row[j] for c, row in zip(row_c, complement)) for j in columns] for row_c in constraint])
        if len(hermite) != 2:
            raise ArithmeticError("the projected constraint lattice must have rank two")
        first, second = (tuple(row) for row in hermite)
        # The basis (h_0, h_1) / e is positively oriented against (o_0, o_1)
        # exactly when the 2 x 2 minors of both on one pair of columns i, j
        # (with nonzero minor of o) have one sign (Cramer's rule).
        top, bottom = oriented_basis.row(0), oriented_basis.row(1)
        i, j, oriented_minor = next((i, j, top[i] * bottom[j] - top[j] * bottom[i]) for i in columns for j in range(i + 1, self.rank) if top[i] * bottom[j] - top[j] * bottom[i])
        if (first[i] * second[j] - first[j] * second[i] > 0) != (oriented_minor > 0):
            first = tuple(-entry for entry in first)
        i, j, minor = next((i, j, first[i] * second[j] - first[j] * second[i]) for i in columns for j in range(i + 1, self.rank) if first[i] * second[j] - first[j] * second[i])
        gram = self.gram_entries()
        gram_first = [sum(g * h for g, h in zip(row, first)) for row in gram]
        gram_second = [sum(g * h for g, h in zip(row, second)) for row in gram]
        square = denominator * denominator
        work = [
            SageQQ(sum(h * g for h, g in zip(first, gram_first))) / square,
            SageQQ(sum(h * g for h, g in zip(first, gram_second))) / square,
            SageQQ(sum(h * g for h, g in zip(second, gram_second))) / square,
        ]
        scale = work[0].denominator().lcm(work[1].denominator()).lcm(work[2].denominator())
        form = (int(scale * work[0]), int(scale * work[1]), int(scale * work[2]))

        def coordinates(row) -> tuple[Rational, Rational]:
            # e * row = x h_0 + y h_1 on columns i, j.
            x, y = denominator * row[i], denominator * row[j]
            return (SageQQ(x * second[j] - y * second[i]) / minor, SageQQ(first[i] * y - first[j] * x) / minor)

        r0_work = self.primitive_coordinates(coordinates(r0))

        # edgewalk.h:763-764, Basis_P_inter_LN = IntersectionLattice_VectorSpace(Latt, Pplane):
        # C cap P on the same coordinate rows as the projection above. A
        # vector c * B of C (B the rows of a basis of C) lies in P exactly
        # when it pairs to zero with every edge wall, so the coefficient
        # lattice of C cap P is the integral left kernel of B G R^T.
        _image, kernel = _hermite_rows([[sum(c * row[w] for c, row in zip(row_c, wall_pairings)) for w in range(len(wall_pairings[0]))] for row_c in constraint])
        # The Hermite normal form of a basis of the integral left kernel is the
        # echelon basis Sage's left_kernel() carries (the form is unique for
        # the kernel lattice).
        coefficients, _zero = _hermite_rows(kernel)
        if len(coefficients) != 2:
            raise ArithmeticError("C cap P must have rank two")
        intersection_rows = [[sum(c * row[j] for c, row in zip(row_c, constraint)) for j in columns] for row_c in coefficients]
        intersection_coordinates = [entry for row in intersection_rows for entry in coordinates(row)]
        if any(entry.denominator() != 1 for entry in intersection_coordinates):
            raise ArithmeticError("the perpendicular intersection must be integral in the projected lattice")
        return _EdgePlane(
            constraint_inverse,
            (first, second),
            denominator,
            form,
            scale,
            r0_work,
            matrix(SageQQ, 2, 2, intersection_coordinates),
        )

    @staticmethod
    def plane_element(plane: Lattices.ParentMethods, coordinates: _RationalRow) -> _Vector:
        ring = plane.base_ring()
        return plane.linear_combination({label: ring(_integer(entry)) for label, entry in zip(plane.module_generating_set(), coordinates, strict=True) if entry})

    @staticmethod
    def primitive_coordinates(coordinates: _RationalRow) -> _Key:
        r"""The coordinates of the primitive integral vector on the ray of nonzero rational coordinates (RemoveFractionVector)."""
        denominator = SageZZ.one()
        for entry in coordinates:
            denominator = denominator.lcm(_rational(entry).denominator())
        integral = [SageZZ(denominator * entry) for entry in coordinates]
        content = SageZZ(0)
        for entry in integral:
            content = content.gcd(entry)
        if content == 0:
            raise ValueError("the zero vector does not determine a primitive ray")
        return tuple(int(entry // content) for entry in integral)

    @cached_method
    def _anisotropic_edge_candidate_family(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        norm: Rational,
    ) -> tuple[_Key, ...]:
        """The complete anisotropic edge-plane family for one root length.

        polyhedral_common edgewalk.h:753-845 computes Anisotropic once per
        e_norm, scales *each* returned primitive vector by all positive
        integers remaining under e_norm, and unfolds the finite quotient
        period.  Only get_next_anisotropic subsequently filters res_norm.
        Reusing the full family across extensions also preserves that
        original enumeration order.
        """
        edge = self.edge_plane(roots, discarded, k, norm)
        form, scale = edge.form, edge.scale
        # edge.form is scale * Gwork and the bound is scale * N, while
        # Anisotropic receives Gwork and res_norm itself (edgewalk.h:765-777);
        # pass the factor so Canonical's floor is taken at the source's scale.
        cycle = anisotropic_cycle(form, scale * _rational(norm), edge.r0, scale=scale)
        if cycle is None:
            return ()
        action, representatives = cycle
        a, b, c = form
        bound = _rational(scale * norm)
        multiples: list[_Key] = []
        for x, y in representatives:
            square = a * x * x + 2 * b * x * y + c * y * y
            if square <= 0:
                raise ArithmeticError("Anisotropic must return primitive positive-square vectors")
            multiplier = 1
            while multiplier * multiplier * square <= bound:
                multiples.append((multiplier * x, multiplier * y))
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

        (p00, p01), (p10, p11) = ((int(entry) for entry in row) for row in action.rows())
        family = list(multiples)
        for _ in range(1, period_order):
            multiples = [(x * p00 + y * p10, x * p01 + y * p11) for x, y in multiples]
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
    ) -> tuple[_Key, ...]:
        # The pinned fixed-norm factorization admits negative residual
        # squares as well as positive ones. Only the zero-square case lacks
        # a finite signed factor-pair enumeration.
        if residual == 0:
            return ()
        edge = self.edge_plane(roots, discarded, k, norm)
        scale = edge.scale
        if self.edge_is_isotropic(roots, discarded, k):
            square = scale * residual
            if square.denominator() != 1:
                return ()
            a, b, c = edge.form
            x0, y0 = edge.r0
            # Pinned edgewalk.h:921-961 uses the *entire* signed factor-pair
            # enumeration, including imprimitive vectors, before ordering
            # roots and checking their ambient integral reconstructions.
            # A vector-orbit representative family is not that candidate set.
            # The source factors the form on integers (Factor_GP_LN); so does
            # binary_split_form_vectors_of_square, with no plane lattice.
            vectors = [item for item in binary_split_form_vectors_of_square(edge.form, int(square)) if _oriented_determinant(edge.r0, item) > 0]
            if self.norm(k) < 0:
                vectors = [(x, y) for x, y in vectors if a * x0 * x + b * (x0 * y + y0 * x) + c * y0 * y > 0]
            vectors.sort(key=cmp_to_key(_orientation_compare))
            return tuple(vectors)
        target_square = _rational(scale * residual)
        if target_square <= 0 or target_square > _rational(scale * norm) or target_square.denominator() != 1:
            return ()
        a, b, c = edge.form
        return tuple((x, y) for x, y in self._anisotropic_edge_candidate_family(roots, discarded, k, norm) if a * x * x + 2 * b * x * y + c * y * y == target_square)

    def extension_roots(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        extension: NormedDynkinExtension,
    ) -> _Roots:
        edge = self.edge_plane(roots, discarded, k, extension.norm)
        work_vectors = self.binary_candidates(roots, discarded, k, extension.norm, extension.residual_norm)
        first, second = edge.basis
        denominator = edge.basis_denominator
        inverse, inverse_denominator = edge.constraint_inverse
        component = extension.u_component
        for x, y in work_vectors:
            ambient = [u + SageQQ(x * a + y * b) / denominator for u, a, b in zip(component, first, second, strict=True)]
            if any(entry.denominator() != 1 for entry in ambient):
                continue
            ambient = [int(entry) for entry in ambient]
            # edgewalk.h:905,988 get_next_anisotropic/get_next_isotropic:
            # SolutionIntMat(Latt, v_T), every integral candidate must
            # belong to the norm-specific lattice.
            if any(sum(a * row[j] for a, row in zip(ambient, inverse)) % inverse_denominator for j in range(self.rank)):
                continue
            root = self.element(tuple(ambient))
            if sum(a * g * b for a, row in zip(ambient, self.gram_entries()) for g, b in zip(row, ambient)) != extension.norm:
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
        # A Lorentzian chamber corner is a nonspacelike ray.  The root
        # congruence and side-pairing tests alone do not imply this: a
        # positive-norm intersection can satisfy both and is not a vertex
        # of the hyperbolic polyhedron.
        if self.norm(generator) > 0:
            return None
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

    def cusp_root(self, divisor: int, component: tuple[Rational, ...], k: _Vector) -> _Vector | None:
        r"""The point u + c k of C on the line through the component u in direction k with the least c > 0.

        lorentzian_linalg.h:160-215 ResolveLattEquation: express u and k in a
        basis of C cap span(u, k), make the k coordinate a single gcd entry by
        a unimodular change, and read the solvability condition and the least
        positive parameter from the remaining entries.
        """
        rows, inverse = self.constraint_rows(divisor)
        u = vector(SageQQ, component) * inverse
        d = self.row(k) * inverse
        span = matrix(SageQQ, [u, d])
        if span.rank() != 2:
            raise ArithmeticError("the parabolic component and the cusp must be independent")
        denominator = SageZZ.one()
        for entry in span.list():
            denominator = denominator.lcm(entry.denominator())
        # IntersectionLattice_VectorSpace(Latt, {u, k}): Z^n cap span in the
        # coordinates of the basis of C.
        intersection = (denominator * span).change_ring(SageZZ).saturation().change_ring(SageQQ)
        u1, u2 = intersection.solve_left(u)
        k1, k2 = intersection.solve_left(d)
        # ComputePairGcd(k1, k2): a unimodular P with (k1, k2) P = (g, 0).
        scale = k1.denominator().lcm(k2.denominator())
        gcd, s, t = SageZZ(scale * k1).xgcd(SageZZ(scale * k2))
        left, right = SageZZ(scale * k1) // gcd, SageZZ(scale * k2) // gcd
        # P = [[s, -right], [t, left]] has determinant s*left + t*right = 1.
        u1_norm = s * u1 + t * u2
        u2_norm = -right * u1 + left * u2
        k1_norm = gcd / scale
        if u2_norm.denominator() != 1:
            return None
        c0 = -u1_norm / k1_norm
        step = 1 / k1_norm
        threshold = -c0 / step
        h = threshold.floor() + 1 if step > 0 else threshold.ceil() - 1
        c = c0 + h * step
        if c <= 0:
            raise ArithmeticError("the cusp equation must return c > 0")
        return self.element(tuple(vector(SageQQ, component) + c * self.row(k)))

    def cusp_roots(self, roots: _Roots, k: _Vector, previous: _Vector) -> _Roots:
        candidates: list[tuple[int, Rational, Rational, _Vector]] = []
        for extension in parabolic_extensions(self.lattice, roots, self.norms):
            root = self.cusp_root(self.divisor(extension.norm), extension.u_component, k)
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
        # The two null lines of the split edge plane L cap P on engine rows:
        # P is the integral left kernel of G R^T (R the rows of the walls),
        # and each line is primitive in L cap P, which is saturated in L.
        gram = self.gram.change_ring(SageZZ)
        # A binary edge has no retained wall roots (l_ui = []). Preserve the
        # ambient column count explicitly: matrix([]) is 0 x 0, not the
        # 0 x rank wall matrix whose transpose multiplies G.
        walls = matrix(SageZZ, len(roots), self.rank,
                       [entry for root in roots for entry in self.key(root)])
        plane_rows = (gram * walls.transpose()).left_kernel_matrix(algorithm="pari").echelon_form()
        plane_gram = plane_rows * gram * plane_rows.transpose()
        choices: dict[_Key, _Vector] = {}
        for line in binary_split_null_basis((plane_gram[0, 0], plane_gram[0, 1], plane_gram[1, 1])):
            candidate = self.element(tuple(vector(SageZZ, line) * plane_rows))
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
    def edge_step(self, vertex: _Vertex, roots: _Roots, discarded: _Vector) -> _Vertex | None:
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
        # With no finite root extension, an anisotropic rational edge plane
        # has no rational ideal endpoint either.  The chamber edge therefore
        # has no vertex in the projective closure, so the reflection chamber
        # cannot be a finite-volume fundamental polyhedron.  This is a
        # nonreflective terminal edge, not a failed cusp computation.
        #
        # Substitution, recorded deliberately: the source has no output for
        # this configuration.  Its isotropic fallback returns an empty
        # generator list on an anisotropic edge plane and it aborts with
        # TerminalException{"Rethink needed"} (edgewalk.h:1117-1122), after
        # which no verdict exists.  The negative verdict comes from Allcock's
        # dichotomy instead (main.tex:1971-1985, 2082-2089): a finite-volume
        # chamber has every edge shared with a neighbouring chamber, so a
        # terminal anisotropic edge rules reflectivity out.
        if not self.edge_is_isotropic(roots, discarded, vertex.generator):
            return None
        return self.isotropic_next_vertex(roots, discarded, vertex.generator)

    @cached_method
    def edge_directions(self, vertex: _Vertex) -> tuple[tuple[_Roots, _Vector], ...]:
        r"""The walls along each outgoing edge and a wall orienting that edge."""
        if len(vertex.roots) == 1:
            return ()
        rows = [self.row(root) for root in vertex.roots]
        cone = _cone_of_rays(self.rank, rows)
        if cone.dim() <= 0:
            return ()
        directions: list[tuple[_Roots, _Vector]] = []
        for inequality in cone.inequalities():
            covector = vector(SageQQ, inequality.A())
            included = tuple(root for root, row in zip(vertex.roots, rows, strict=True) if covector.dot_product(row) == 0)
            omitted = tuple(root for root, row in zip(vertex.roots, rows, strict=True) if covector.dot_product(row) != 0)
            if not omitted:
                raise ArithmeticError("a root-cone facet must omit a wall")
            # edgewalk.h:525-539 GetAdjacencyDirection: v_disc is the last
            # root outside the facet. An ideal vertex with several affine
            # components has facets that omit several walls (U+<2>+<2>
            # under K3 norms). Every omitted wall gives the same edge: the
            # roots vanish at k, so a point x = a k + c d of the edge plane
            # P = (included)^perp pairs with an omitted root r as
            # b(r, x) = c b(r, d), where d spans the edge ray and
            # b(r, d) < 0 because r is not a wall of that edge. The tests
            # in edge_plane, vertex_from_root and isotropic_next_vertex
            # read only this sign, which is -sign(c) for every omitted r.
            directions.append((included, omitted[-1]))
        return tuple(directions)

    @cached_method
    def coloured_points(self) -> tuple[Parent, tuple[Element, Element]]:
        r"""The product L x {0, 1} and its two colours, built once per walk."""
        colours = Sets.Δ[1]
        return Sets().product((self.lattice, colours)), (colours[0], colours[1])

    def coloured_family(self, pairs: tuple[tuple[_Vector, int], ...]) -> Element:
        r"""The coloured family i -> (v_i, c_i) in Sets().Mor(I, L x {0, 1})."""
        points, colours = self.coloured_points()

        # The stabilizer and transporter read each value many times; each
        # point of L x {0, 1} is built on its first reading only.
        @cache
        def point(index: int) -> Element:
            return points((pairs[index][0], colours[pairs[index][1]]))

        return Sets().Mor(Sets.Δ[len(pairs) - 1], points)(lambda index: point(int(index)))

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

    @cached_method
    def configuration_invariant(self, vertex: _Vertex) -> tuple:
        r"""An isometry invariant of the coloured configuration of ``vertex``.

        # source: edgewalk.h:1457 func_insert_vertex tests LORENTZ_TestEquivalence
        # only for known vertices of equal hash, the hash of the canonical
        # weighted matrix (fund_domain_vertices.h:365, ComputeHashWeightMatrix_raw).

        Each vector contributes its colour, its norm and the sorted coloured
        pairings with the other vectors. An isometry carrying one
        configuration to the other up to a colour-preserving reindexing
        preserves every entry, so configurations with different invariants
        have no transporter, and only equal invariants are tested.
        """
        vectors = (*vertex.roots, vertex.generator)
        colours = (0,) * len(vertex.roots) + (1,)
        return tuple(
            sorted(
                (
                    colour,
                    self.norm(vector_),
                    tuple(sorted((colours[other], self.pairing(vector_, vectors[other])) for other in range(len(vectors)) if other != index)),
                )
                for index, (vector_, colour) in enumerate(zip(vectors, colours, strict=True))
            )
        )

    def isometry_key(self, isometry: LatticeIsometryMethods) -> tuple[_Key, ...]:
        return tuple(self.key(isometry(generator)) for generator in _coordinate_basis(self.lattice))

    def common_fixed_space(self, invariant: Matrix_rational_dense, isometry_matrix: Matrix_integer_dense) -> Matrix_rational_dense | None:
        r"""The rows of the rational space fixed by the group generated so far, or None once that group is infinite.

        A finite cone-preserving group fixes a timelike vector by averaging.
        Conversely its integral stabilizer is discrete in a compact real
        orthogonal group, hence finite. The common fixed space decides both
        directions, without separate order computations for the generators.
        An isometry acts on rows by x -> x M, so its fixed space is the left
        kernel of M - 1; the negative index of the form on the intersection
        is PARI qfsign of its Gram matrix.
        """
        fixed_rows = (isometry_matrix - 1).change_ring(SageQQ).left_kernel_matrix()
        fixed = invariant.row_space().intersection(fixed_rows.row_space()).basis_matrix()
        if not fixed.nrows():
            return None
        _positive, negative = pari(fixed * self.gram * fixed.transpose()).qfsign()
        return None if int(negative) == 0 else fixed

    @cached_method
    def wall_covector(self, root: _Vector) -> Element:
        r"""The covector -b(root, -) of a wall inequality, computed once per root (-G * root, vinberg.h:1366 and edgewalk.h:1721)."""
        return -root.to_covector()

    def root_chamber_certificate(self, roots: _Roots) -> bool:
        r"""Upstream LORENTZ_RunEdgewalkAlgorithm dual-description ray test.

        The inequalities are -b(root, x) >= 0.  A pointed rational
        chamber is tested on its exact extreme rays, as in the pinned
        upstream dual-description termination condition.
        """
        if not roots:
            return False
        # DualDescriptionStandard of ListIneq = -LVect * G (edgewalk.h:1721-1726)
        # on the Gram matrix and the root coordinates: the cone
        # {x : -b(root, x) >= 0} and its extreme rays.
        chamber = Polyhedron(
            ieqs=[(0, *(-(self.gram * self.row(root)))) for root in roots],
            base_ring=SageQQ,
            backend="ppl",
        )
        if chamber.lines():
            return False
        rays = tuple(vector(SageQQ, ray.vector()) for ray in chamber.rays())
        # The dual-description norm check alone is vacuous on the zero
        # cone.  A hyperbolic chamber must have full span and contain a
        # timelike vector.  Moreover, different nonspacelike rays must
        # belong to the same closed time cone.  Otherwise their positive
        # span contains spacelike vectors despite each ray having q <= 0.
        # Pairwise nonpositive inner products make the entire positive span
        # nonspacelike by expanding the square of any positive combination.
        if not rays or matrix(SageQQ, rays).rank() != self.rank:
            return False
        if any(self.pair(ray, ray) > 0 for ray in rays):
            return False
        if any(self.pair(left, right) > 0 for index, left in enumerate(rays) for right in rays[index + 1 :]):
            return False
        return any(self.pair(ray, ray) < 0 for ray in rays) or any(self.pair(left, right) < 0 for index, left in enumerate(rays) for right in rays[index + 1 :])

    def report(
        self,
        vertices: list[_Vertex],
        reflective: bool,
        simple_roots: _Roots,
        isometries: list[LatticeIsometryMethods],
        reason_non_reflective: str | None = None,
    ) -> AllcockEdgewalkReport:
        cones = RationalPolyhedralCones(self.lattice)
        return AllcockEdgewalkReport(
            self.lattice,
            reflective=reflective,
            simple_roots=simple_roots,
            vertices=tuple(cones.from_rays([vertex.generator], wall_roots=vertex.roots) for vertex in vertices),
            isometry_generators=tuple(isometries),
            reason_non_reflective=reason_non_reflective,
        )

    def edge_orbits(self, vertex: _Vertex, local_generators: tuple[LatticeIsometryMethods, ...]) -> tuple[tuple[_Roots, _Vector], ...]:
        r"""Outgoing edges modulo the supplied full vertex stabilizer, the first direction of each orbit in facet order.

        edgewalk.h:1555 OrbitSplittingSet(vf, GRP1_integral): the facets of
        the root cone split into orbits of the stabilizer acting on the roots.
        A generator acts on root coordinates by x -> x * M, M the rows of the
        images of the basis (``isometry_key``), so the orbits are swept on the
        root coordinates and the integral matrices of the generators.
        """
        directions = self.edge_directions(vertex)
        root_rows = [tuple(vector(SageZZ, self.key(root)) for root in roots) for roots, _discarded in directions]
        position = {frozenset(tuple(row) for row in rows): index for index, rows in enumerate(root_rows)}
        actions = tuple(matrix(SageZZ, self.isometry_key(generator)) for generator in local_generators)
        reached = [False] * len(directions)
        representatives = []
        for start in range(len(directions)):
            if reached[start]:
                continue
            reached[start] = True
            representatives.append(start)
            frontier = [start]
            while frontier:
                images = []
                for index in frontier:
                    for action in actions:
                        image = position[frozenset(tuple(row * action) for row in root_rows[index])]
                        if not reached[image]:
                            reached[image] = True
                            images.append(image)
                frontier = images
        return tuple(directions[index] for index in representatives)

    def binary_walk(self, initial: _Vertex | None = None) -> AllcockEdgewalkReport:
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
        supplied = self.supplied_initial_vertex
        if supplied is not None and initial is None:
            point = supplied[0] if isinstance(supplied, tuple) else supplied
            if point.parent() is not self.lattice or self.norm(point) > 0:
                raise ValueError("the starting vertex must be a nonspacelike vector of the walking lattice")
        if supplied is None and initial is None and self.lattice.is_isotropic():
            # The DEFAULT isotropic_vinberg start is a null ray. In rank two
            # its perpendicular quotient is zero-dimensional, so the source's
            # get_simple_cone raises NonReflectivityException("list of roots
            # is empty"). A FileVertex/FileVertexRoots start takes the other
            # branch of get_initial_vertex and must not be intercepted here.
            return self.report([], False, (), [], "list of roots is empty")
        if initial is None:
            if supplied is None:
                initial = self.initial_vertex()
            else:
                # As in the higher-rank walk, a bare vertex has its simple
                # cone computed while a (vertex, roots) pair is taken as
                # given and checked below (edgewalk.h:2100-2115).
                initial = self.prescribed_vertex(supplied)
        if initial is None:
            return self.report([], False, (), [], self.initial_nonreflectivity_reason)
        if initial is not None and (initial.generator.parent() is not self.lattice
                                   or self.norm(initial.generator) > 0
                                   or any(root.parent() is not self.lattice
                                          or self.pairing(initial.generator, root) != 0
                                          or self.norm(root) not in self.norms
                                          for root in initial.roots)):
            raise ValueError("the binary starting vertex has invalid wall roots")
        if initial is None or len(initial.roots) != 1:
            raise ValueError("a binary starting vertex must have exactly one wall root")
        (root,) = initial.roots
        neighbor = self.edge_step(initial, (), root)
        if neighbor is None:
            return self.report(
                [initial], False, (), [],
                "An anisotropic binary edge has no finite root extension or rational ideal endpoint",
            )
        symmetry = self.lattice.O().transporter(self.vertex_configuration(initial), self.vertex_configuration(neighbor))
        vertices = [initial] if symmetry is not None else [initial, neighbor]
        if self.lattice.is_isotropic():
            # In a hyperbolic basis an integral isometry maps each of the two
            # rational null lines to itself or the other. Its image on a null
            # basis consequently has diagonal or antidiagonal entries a,d
            # with a*d = 1, hence a,d are units in either case. O(L) is
            # therefore FINITE. A finite group acting on H^1 = R has quotient
            # R or a half-line, both of infinite
            # length, so no reflection subgroup can be reflective. This exact
            # arithmetic argument covers prescribed starts too; unlike the
            # default null start it does not pretend that their wall list was
            # empty. In particular an ideal neighbour has no wall root and
            # must not be unpacked as another reflective endpoint.
            return self.report(
                vertices, False, (), [] if symmetry is None else [symmetry],
                "An isotropic binary lattice has finite integral orthogonal group, so its reflection group has infinite covolume",
            )
        (neighbor_root,) = neighbor.roots
        return self.report(vertices, True, (root, neighbor_root), [] if symmetry is None else [symmetry])

    def walk(self, initial: _Vertex | None = None, visit: Callable[[_Vertex], bool] | None = None) -> AllcockEdgewalkReport:
        r"""edgewalk.h:1697-1792. A visit callback returning True stops at that vertex."""
        if self.rank == 2 and visit is None:
            return self.binary_walk(initial)
        if initial is None:
            supplied = self.supplied_initial_vertex
            if supplied is None:
                initial = self.initial_vertex()
            else:
                # edgewalk.h:2100-2143 FileVertex / FileVertexRoots. A bare
                # vertex has its simple cone computed by get_simple_cone,
                # including the cusp route for a null ray and the check that
                # the perpendicular roots are a full-rank corner (:1848-1873);
                # a (vertex, roots) pair is taken as given, exactly as the
                # source trusts FileVertexRoots (:2100-2110).
                if not isinstance(supplied, tuple) and (
                    supplied.parent() is not self.lattice or self.norm(supplied) > 0
                ):
                    raise ValueError("the starting vertex must be a nonspacelike lattice vector")
                initial = self.prescribed_vertex(supplied)
        if initial is not None and (
            initial.generator.parent() is not self.lattice
            or self.norm(initial.generator) > 0
            or any(root.parent() is not self.lattice
                   or self.pairing(initial.generator, root) != 0
                   or self.norm(root) not in self.norms for root in initial.roots)
        ):
            raise ValueError("the supplied starting vertex has invalid wall roots")
        if initial is None:
            # The corner search can end in the source's own nonreflectivity
            # decision (edgewalk.h:1985-1999); report its reason rather than a
            # negative verdict with no mathematical justification.
            return self.report([], False, (), [], self.initial_nonreflectivity_reason)
        if visit is not None and visit(initial):
            return self.report([initial], False, (), [])

        vertices = [initial]
        configurations = [self.vertex_configuration(initial)]
        seen: set[tuple[_Key, ...]] = {self.identity_key}
        isometries: list[LatticeIsometryMethods] = []
        # A rational row x is fixed by a chamber symmetry with row matrix M
        # precisely when x(M-I)=0. A finite cone-preserving group fixes a
        # timelike vector by averaging; conversely the integral stabilizer of
        # a timelike vector is discrete in a compact orthogonal group and is
        # therefore finite. Store the complete common fixed space exactly.
        isometry_matrices: list[Matrix_integer_dense] = []
        invariant = identity_matrix(SageQQ, self.rank)
        nonreflective = False
        nonreflective_reason = None

        def admit(isometry: LatticeIsometryMethods) -> bool:
            r"""Adjoin the isometry, preserving an exact infinite-group certificate."""
            nonlocal invariant, nonreflective, nonreflective_reason
            key = self.isometry_key(isometry)
            if key in seen:
                return True
            seen.add(key)
            isometries.append(isometry)
            action = matrix(SageZZ, key)
            isometry_matrices.append(action)
            if not self.early_termination:
                # edgewalk.h:1764-1774 calls group_tester.GeneratorUpdate and
                # asks for its finiteness status ONLY inside this branch. The
                # non-terminating source mode continues collecting every
                # vertex orbit and does not manufacture a nonreflectivity
                # reason from the generators it happens to have found.
                return True
            if nonreflective:
                return True
            fixed = self.common_fixed_space(invariant, action)
            if fixed is None:
                nonreflective = True
                nonreflective_reason = "The generated chamber-isometry group has no nonzero common fixed timelike vector"
                return False
            invariant = fixed
            return True

        def certified_roots() -> _Roots:
            # LORENTZ_RunEdgewalkAlgorithm::f_try_terminate: close the
            # accumulated walls under every discovered chamber symmetry.
            # compute_full_root_orbit_iter (edgewalk.h:1216-1258): the
            # orbits of the simple roots under the generators, on the
            # coordinates of the roots and the integral matrices of the
            # isometries; each root is wrapped once per walk.
            for vertex in vertices:
                for root in vertex.roots:
                    root_elements.setdefault(self.key(root), root)
            frontier = [vector(SageZZ, key) for key in root_elements]
            closed = set(root_elements)
            while frontier:
                images = []
                for row in frontier:
                    for generator in isometry_matrices:
                        image = row * generator
                        key = tuple(image)
                        if key not in closed:
                            closed.add(key)
                            images.append(image)
                frontier = images
            for key in closed:
                if key not in root_elements:
                    root_elements[key] = self.element(key)
            return tuple(root_elements[key] for key in sorted(closed))

        certified_root_count = 0
        root_elements: dict[_Key, _Vector] = {}

        def try_reflectivity_certificate() -> AllcockEdgewalkReport | None:
            # LORENTZ_RunEdgewalkAlgorithm::f_try_terminate tests the dual
            # description only when the closed root orbit has grown.
            nonlocal certified_root_count
            if nonreflective:
                return None
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
                    return self.report(vertices, False, (), isometries, nonreflective_reason)
                if len(isometries) > previous_generators:
                    certificate = try_reflectivity_certificate()
                    if certificate is not None:
                        return certificate
            for roots, discarded in self.edge_orbits(vertex, local_generators):
                neighbour = self.edge_step(vertex, roots, discarded)
                if neighbour is None:
                    return self.report(
                        vertices, False, (), isometries,
                        "An anisotropic chamber edge has no finite root extension or rational ideal endpoint",
                    )
                neighbour_configuration = self.vertex_configuration(neighbour)
                invariant_of_neighbour = self.configuration_invariant(neighbour)
                isometries_to_neighbour = (
                    self.lattice.O().transporter(known, neighbour_configuration)
                    for known, known_vertex in zip(configurations, vertices, strict=True)
                    if self.configuration_invariant(known_vertex) == invariant_of_neighbour
                )
                transporter = next(
                    (isometry for isometry in isometries_to_neighbour if isometry is not None),
                    None,
                )
                if transporter is None:
                    vertices.append(neighbour)
                    configurations.append(neighbour_configuration)
                    if visit is not None and visit(neighbour):
                        return self.report(vertices, False, (), isometries)
                    certificate = try_reflectivity_certificate()
                    if certificate is not None:
                        return certificate
                else:
                    previous_generators = len(isometries)
                    if not admit(transporter):
                        return self.report(vertices, False, (), isometries, nonreflective_reason)
                    if len(isometries) > previous_generators:
                        certificate = try_reflectivity_certificate()
                        if certificate is not None:
                            return certificate
            position += 1

        # Early detection would already have returned above. In the source's
        # non-terminating mode no group-finiteness verdict is made; closing
        # the roots (or discovering that this closure does not terminate) is
        # intentionally left to the same dual-description computation as the
        # pinned implementation (edgewalk.h:1311-1316,1700-1747).
        simple_roots = () if nonreflective else certified_roots()
        return self.report(vertices, not nonreflective and self.root_chamber_certificate(simple_roots),
                           simple_roots, isometries, nonreflective_reason)


__all__ = ["edgewalk_fundamental_domain", "edgewalk_isometry"]
