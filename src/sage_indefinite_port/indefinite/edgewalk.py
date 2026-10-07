"""Native Sage implementation of Allcock's Lorentzian edgewalk.

This is a direct translation of the mathematical path in the pinned
polyhedral_common edgewalk source. Heavy primitive operations are delegated
to Sage exact linear algebra, integer lattices, cones and finite groups; no
polyhedral_common binary is invoked.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import cmp_to_key
from itertools import combinations
from math import gcd
from typing import TypedDict

from dzack_research.preamble.all import ZZ, HyperbolicLattices, Lattices
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from sage.arith.misc import xgcd
from sage.geometry.cone import Cone
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.modules.free_module import FreeModule_submodule_with_basis_pid
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ

from sage_indefinite_port.backends.canonization import (
    CellConfiguration,
    cell_stabilizer,
    cell_transporter,
)
from sage_indefinite_port.indefinite.edgewalk_extensions import (
    NormedDynkinExtension,
    compute_possible_extensions,
)
from sage_indefinite_port.indefinite.edgewalk_rank2 import (
    anisotropic_cycle,
    fixed_norm_vectors_isotropic,
    isotropic_factorization,
    oriented_determinant,
    primitive_isotropic_vectors,
    quadratic_eval,
    scalar_eval,
)
from sage_indefinite_port.readiness import unfinished

type _Row = tuple[int, ...]
type _Roots = tuple[_Row, ...]
type _Gram = Matrix_integer_dense | Matrix_rational_dense
type _RationalRow = Sequence[int | Integer | Rational] | FreeModuleElement[Rational] | FreeModuleElement[Integer]
type _ConstraintLattice = FreeModule_submodule_with_basis_pid[Integer]


class EdgewalkRecord(TypedDict):
    """Allcock's fundamental-domain record, in coordinate rows of the input Gram basis."""

    simple_root_rows: _Roots
    vertices: tuple[tuple[_Row, _Roots], ...]
    is_reflective: bool
    isometry_generator_rows: tuple[_Roots, ...]


def _signature_pair(gram: _Gram) -> tuple[int, int]:
    diagonal = QuadraticForm(SageQQ, 2 * matrix(SageQQ, gram)).rational_diagonal_form().matrix()
    positive = sum(diagonal[index, index] > 0 for index in range(diagonal.nrows()))
    negative = sum(diagonal[index, index] < 0 for index in range(diagonal.nrows()))
    return positive, negative


@dataclass(frozen=True)
class _Vertex:
    generator: _Row
    roots: _Roots


@dataclass(frozen=True)
class _Candidate:
    alpha: _Row
    residual_norm: Rational
    root_norm: Rational
    vertex: _Vertex


def _integer_gram(gram: Matrix_integer_dense) -> Matrix_integer_dense:
    result = matrix(SageZZ, [[int(entry) for entry in row] for row in gram])
    if not result.is_square() or result != result.transpose() or result.det() == 0:
        raise ValueError("Allcock edgewalk needs a nondegenerate symmetric Gram matrix")
    positive, negative = _signature_pair(result)
    if negative != 1:
        raise ValueError(f"Allcock edgewalk uses signature (n,1), but the Gram matrix has signature ({positive},{negative})")
    return result


def _q(gram: _Gram, row: _RationalRow) -> Rational:
    value = vector(SageQQ, row)
    return (value * gram * value.column())[0]


def _pair(gram: _Gram, left: _RationalRow, right: _RationalRow) -> Rational:
    return (vector(SageQQ, left) * gram * vector(SageQQ, right).column())[0]


def _primitive_row(row: _RationalRow) -> _Row:
    rational = vector(SageQQ, row)
    denominator = SageZZ.one()
    for entry in rational:
        denominator = denominator.lcm(entry.denominator())
    integers = [SageZZ(denominator * entry) for entry in rational]
    content = 0
    for value in integers:
        content = gcd(content, abs(int(value)))
    if content == 0:
        raise ValueError("the zero ray has no primitive integral generator")
    integers = [entry // content for entry in integers]
    return tuple(int(entry) for entry in integers)


def _row_of_element(lattice: Lattices.ParentMethods, element: Lattices.ElementMethods) -> _Row:
    coordinates = element.to_vector()
    return tuple(int(coordinates(label)) for label in lattice.module_generating_set())


def _corner_from_roots(gram: Matrix_integer_dense, roots: _Roots) -> _Vertex | None:
    rank = gram.nrows()
    if len(roots) < rank - 1:
        return None
    candidates: list[_Vertex] = []
    root_matrix = matrix(SageQQ, roots)
    for positions in combinations(range(len(roots)), rank - 1):
        selected = root_matrix.matrix_from_rows(positions)
        null = (selected * gram).right_kernel_matrix()
        if null.nrows() != 1:
            continue
        generator = vector(SageQQ, null.row(0))
        pairings = [(vector(SageQQ, root) * gram * generator.column())[0] for root in roots]
        if all(value >= 0 for value in pairings):
            generator = -generator
            pairings = [-value for value in pairings]
        if any(value > 0 for value in pairings):
            continue
        if (generator * gram * generator.column())[0] > 0:
            continue
        primitive = _primitive_row(generator)
        incident = tuple(root for root in roots if _pair(gram, root, primitive) == 0)
        if matrix(SageQQ, incident).rank() != rank - 1:
            continue
        candidates.append(_Vertex(primitive, incident))
    if not candidates:
        return None
    candidates.sort(
        key=lambda vertex: (
            _q(gram, vertex.generator) != 0,
            tuple(abs(entry) for entry in vertex.generator),
            vertex.generator,
        )
    )
    return candidates[0]


def _initial_vertex(lattice: HyperbolicLattices.ParentMethods, gram: Matrix_integer_dense) -> _Vertex:
    for bound in (8, 16, 32, 64, 128):
        _complete, roots = lattice._vinberg_search(None, bound, max(256, 8 * bound))
        rows = tuple(_row_of_element(lattice, root) for root in roots)
        corner = _corner_from_roots(gram, rows)
        if corner is not None:
            return corner
    raise RuntimeError("bounded Vinberg initialization did not expose one chamber corner")


def _constraint_lattice(gram: Matrix_integer_dense, norm: Rational) -> _ConstraintLattice:
    rank = gram.nrows()
    ambient_integer = SageZZ**rank
    ambient_rational = SageQQ**rank
    scaled_dual = (SageQQ(norm) / 2) * gram.change_ring(SageQQ).inverse()
    dual_lattice = ambient_rational.span(scaled_dual.rows(), SageZZ)
    return ambient_integer.intersection(dual_lattice)


def _positive_vector_in_plane(gram: Matrix_integer_dense, plane_basis: Matrix_rational_dense, k: _Row) -> FreeModuleElement[Rational]:
    half = SageQQ(1)
    k_row = vector(SageQQ, k)
    while True:
        for basis_row in plane_basis.rows():
            for sign in (-1, 1):
                candidate = k_row + sign * half * vector(SageQQ, basis_row)
                if _q(gram, candidate) > 0:
                    if _pair(gram, candidate, k) > 0:
                        candidate = -candidate
                    return candidate
        half /= 2


def _projection_data(
    gram: Matrix_integer_dense,
    roots: _Roots,
    discarded: _Row,
    k: _Row,
    norm: Rational,
) -> tuple[_ConstraintLattice, Matrix_rational_dense, Matrix_rational_dense, FreeModuleElement[Integer], Matrix_integer_dense]:
    rank = gram.nrows()
    root_matrix = matrix(SageQQ, roots)
    if roots:
        root_gram = root_matrix * gram * root_matrix.transpose()
        plane_basis = (root_matrix * gram).right_kernel().basis_matrix()
    else:
        root_gram = matrix(SageQQ, 0, 0)
        plane_basis = matrix.identity(SageQQ, rank)
    if plane_basis.nrows() != 2:
        raise ArithmeticError("an edge must have a two-dimensional perpendicular plane")

    equations = list((root_matrix * gram).rows()) if roots else []
    equations.append(tuple(vector(SageQQ, k) * gram))
    r0_space = matrix(SageQQ, equations).right_kernel_matrix()
    if r0_space.nrows() != 1:
        raise ArithmeticError("the edge orientation vector must be one-dimensional")
    r0 = vector(SageQQ, r0_space.row(0))
    if _pair(gram, r0, discarded) > 0:
        r0 = -r0

    if _q(gram, k) < 0:
        oriented_basis = matrix(SageQQ, [r0, vector(SageQQ, k)])
    else:
        positive = _positive_vector_in_plane(gram, plane_basis, k)
        oriented_basis = matrix(SageQQ, [positive, vector(SageQQ, k)])
        r0 = -vector(SageQQ, k)

    constraint = _constraint_lattice(gram, norm)
    constraint_rows = matrix(SageQQ, constraint.basis_matrix())
    projected: list[FreeModuleElement[Rational]] = []
    for row in constraint_rows.rows():
        if roots:
            coefficients = (row * gram * root_matrix.transpose()) * root_gram.inverse()
            root_component = coefficients * root_matrix
        else:
            root_component = vector(SageQQ, [0] * rank)
        projected.append(vector(SageQQ, row) - vector(SageQQ, root_component))
    projected_module = (SageQQ**rank).span(projected, SageZZ)
    basis = matrix(SageQQ, list(projected_module.basis_matrix().rows()))
    if basis.nrows() != 2:
        raise ArithmeticError("the projected constraint lattice must have rank two")

    expression = oriented_basis.transpose().solve_right(basis.transpose()).transpose()
    if expression.det() < 0:
        basis.rescale_row(0, -1)
    work_gram = basis * gram * basis.transpose()
    r0_coordinates = basis.transpose().solve_right(r0.column()).column(0)
    r0_integral = _primitive_row(r0_coordinates)

    plane = (SageQQ**rank).subspace(basis.rows())
    intersection = constraint.intersection(plane)
    intersection_rows = matrix(SageQQ, intersection.basis_matrix())
    intersection_coordinates = basis.transpose().solve_right(intersection_rows.transpose()).transpose()
    if any(entry.denominator() != 1 for entry in intersection_coordinates.list()):
        raise ArithmeticError("the perpendicular intersection must be integral in the projected lattice")
    intersection_coordinates = matrix(SageZZ, intersection_coordinates)
    return constraint, basis, work_gram, vector(SageZZ, r0_integral), intersection_coordinates


def _class_action_order(transform: Matrix_integer_dense, intersection_coordinates: Matrix_integer_dense) -> int:
    transform = matrix(SageZZ, transform).transpose()
    if intersection_coordinates.det() == 0:
        raise ArithmeticError("the perpendicular sublattice must have finite index")
    sublattice = (SageZZ**2).span(intersection_coordinates.rows(), SageZZ)
    standard = (vector(SageZZ, [1, 0]), vector(SageZZ, [0, 1]))
    power = matrix.identity(SageZZ, 2)
    quotient_size = abs(int(intersection_coordinates.det()))
    guard = max(12, 12 * quotient_size * quotient_size)
    for order in range(1, guard + 1):
        power *= transform
        preserves = all(vector(SageZZ, power * vector(SageZZ, row)) in sublattice for row in intersection_coordinates.rows())
        trivial = all(vector(SageZZ, power * basis - basis) in sublattice for basis in standard)
        if preserves and trivial:
            return order
    raise ArithmeticError("failed to close the finite projected-lattice class action")


def _orientation_compare(left: FreeModuleElement[Integer], right: FreeModuleElement[Integer]) -> int:
    determinant = oriented_determinant(left, right)
    if determinant > 0:
        return -1
    if determinant < 0:
        return 1
    return 0


def _extension_roots(
    gram: Matrix_integer_dense,
    roots: _Roots,
    discarded: _Row,
    k: _Row,
    extension: NormedDynkinExtension,
) -> _Roots:
    constraint, basis, work_gram, r0, intersection = _projection_data(gram, roots, discarded, k, extension.norm)
    residual = SageQQ(extension.residual_norm)
    if residual <= 0:
        return ()
    work_vectors: list[FreeModuleElement[Integer]] = []
    if isotropic_factorization(work_gram) is not None:
        vectors = list(fixed_norm_vectors_isotropic(work_gram, residual))
        if _q(gram, k) < 0:
            vectors = [item for item in vectors if scalar_eval(work_gram, r0, item) > 0]
        vectors = [item for item in vectors if oriented_determinant(r0, item) > 0]
        vectors.sort(key=cmp_to_key(_orientation_compare))
        work_vectors.extend(vectors)
    else:
        cycle = anisotropic_cycle(work_gram, extension.norm, r0)
        if cycle is None:
            return ()
        transform, representatives = cycle
        primitive: list[FreeModuleElement[Integer]] = []
        for item in representatives:
            square = quadratic_eval(work_gram, item)
            multiplier = 1
            while multiplier * multiplier * square <= extension.norm:
                candidate = multiplier * item
                if quadratic_eval(work_gram, candidate) == residual:
                    primitive.append(candidate)
                multiplier += 1
        order = _class_action_order(transform, intersection)
        action = matrix.identity(SageZZ, 2)
        for _ in range(order):
            work_vectors.extend(action * item for item in primitive)
            action = transform.transpose() * action

    answer: list[_Row] = []
    component = vector(SageQQ, extension.u_component)
    for work_vector in work_vectors:
        ambient = component + basis.transpose() * vector(SageQQ, work_vector)
        if any(entry.denominator() != 1 for entry in ambient):
            continue
        integral = vector(SageZZ, ambient)
        if integral not in constraint:
            continue
        row = tuple(int(entry) for entry in integral)
        if _q(gram, row) != extension.norm:
            raise ArithmeticError("rank-two reconstruction changed the root norm")
        answer.append(row)
    return tuple(dict.fromkeys(answer))


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


def _candidate_compare(left: _Candidate, right: _Candidate, gram: Matrix_integer_dense, k: _Row) -> int:
    left_scal = -_pair(gram, left.alpha, k)
    right_scal = -_pair(gram, right.alpha, k)
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


def _vertex_from_root(gram: Matrix_integer_dense, k: _Row, roots: _Roots, discarded: _Row, alpha: _Row) -> _Vertex | None:
    extended = tuple(roots) + (alpha,)
    kernel = (matrix(SageQQ, extended) * gram).right_kernel_matrix()
    if kernel.nrows() != 1:
        return None
    generator = vector(SageQQ, kernel.row(0))
    old_pairing = _pair(gram, generator, k)
    if old_pairing > 0:
        generator = -generator
        old_pairing = -old_pairing
    if old_pairing == 0 or _pair(gram, discarded, generator) >= 0:
        return None
    return _Vertex(_primitive_row(generator), extended)


def _rational_gcd_pair(left: Rational, right: Rational) -> tuple[Rational, Matrix_rational_dense]:
    left, right = SageQQ(left), SageQQ(right)
    denominator = left.denominator().lcm(right.denominator())
    first = SageZZ(left * denominator)
    second = SageZZ(right * denominator)
    common, s, t = xgcd(first, second)
    if common < 0:
        common, s, t = -common, -s, -t
    if common == 0:
        raise ArithmeticError("the cusp direction vanished in lattice coordinates")
    transform = matrix(
        SageQQ,
        [[s, -second // common], [t, first // common]],
    )
    return SageQQ(common) / denominator, transform


def _resolve_cusp_lattice_equation(constraint: _ConstraintLattice, component: tuple[Rational, ...], k: _Row) -> _Row | None:
    space = (SageQQ ** len(k)).span((vector(SageQQ, component), vector(SageQQ, k)))
    intersection = constraint.intersection(space)
    basis = matrix(SageQQ, intersection.basis_matrix())
    if basis.nrows() != 2:
        return None
    u_coordinates = basis.transpose().solve_right(vector(SageQQ, component).column()).column(0)
    k_coordinates = basis.transpose().solve_right(vector(SageQQ, k).column()).column(0)
    gcd_value, transform = _rational_gcd_pair(k_coordinates[0], k_coordinates[1])
    normalized_u = transform.transpose() * u_coordinates
    normalized_k = transform.transpose() * k_coordinates
    if normalized_k[1] != 0 or normalized_u[1].denominator() != 1:
        return None
    c0 = -normalized_u[0] / gcd_value
    step = 1 / gcd_value
    threshold = -c0 / step
    if step > 0:
        h = threshold.ceil()
        if h == threshold:
            h += 1
    else:
        h = threshold.floor()
        if h == threshold:
            h -= 1
    c = c0 + h * step
    if c <= 0:
        raise ArithmeticError("the source cusp equation must return c>0")
    result = vector(SageQQ, component) + c * vector(SageQQ, k)
    if any(entry.denominator() != 1 for entry in result):
        raise ArithmeticError("the solved cusp root is not integral")
    integral = vector(SageZZ, result)
    if integral not in constraint:
        raise ArithmeticError("the solved cusp root missed the constraint lattice")
    return tuple(int(entry) for entry in integral)


type _CuspCandidate = tuple[int, Rational, Rational, _Row]


def _cusp_roots(gram: Matrix_integer_dense, roots: _Roots, k: _Row, previous: _Row, norms: tuple[int, ...]) -> _Roots:
    candidates: list[_CuspCandidate] = []
    for extension in compute_possible_extensions(gram, roots, norms, only_spherical=False):
        if extension.residual_norm != 0:
            continue
        constraint = _constraint_lattice(gram, extension.norm)
        root = _resolve_cusp_lattice_equation(constraint, extension.u_component, k)
        if root is None:
            continue
        scalar = -_pair(gram, previous, root)
        sign = (scalar > 0) - (scalar < 0)
        candidates.append((sign, scalar * scalar / extension.norm, extension.norm, root))

    def compare(left: _CuspCandidate, right: _CuspCandidate) -> int:
        first = _signed_sqrt_compare(left[0], left[1], right[0], right[1])
        if first:
            return first
        if left[2] == right[2]:
            return 0
        return -1 if left[2] < right[2] else 1

    candidates.sort(key=cmp_to_key(compare))
    accepted = list(roots)
    for _sign, _quant, _norm, root in candidates:
        if all(_pair(gram, old, root) <= 0 for old in accepted):
            accepted.append(root)
    return tuple(accepted)


def _isotropic_next_vertex(gram: Matrix_integer_dense, roots: _Roots, discarded: _Row, k: _Row, norms: tuple[int, ...]) -> _Vertex:
    root_matrix = matrix(SageQQ, roots)
    plane_basis = (root_matrix * gram).right_kernel().basis_matrix() if roots else matrix.identity(SageQQ, gram.nrows())
    reduced = plane_basis * gram * plane_basis.transpose()
    choices: list[_Row] = []
    for coordinates in primitive_isotropic_vectors(reduced):
        generator = plane_basis.transpose() * vector(SageQQ, coordinates)
        for sign in (1, -1):
            candidate = sign * generator
            if _pair(gram, candidate, k) < 0 and _pair(gram, discarded, candidate) < 0:
                choices.append(_primitive_row(candidate))
    if len(set(choices)) != 1:
        raise ArithmeticError("the edge must determine one forward isotropic endpoint")
    generator = choices[0]
    cusp_roots = _cusp_roots(gram, roots, generator, k, norms)
    return _Vertex(generator, cusp_roots)


def _edge_step(gram: Matrix_integer_dense, vertex: _Vertex, roots: _Roots, discarded: _Row, norms: tuple[int, ...]) -> _Vertex:
    candidates: list[_Candidate] = []
    for extension in compute_possible_extensions(gram, roots, norms, only_spherical=True):
        for alpha in _extension_roots(gram, roots, discarded, vertex.generator, extension):
            next_vertex = _vertex_from_root(gram, vertex.generator, roots, discarded, alpha)
            if next_vertex is not None:
                candidates.append(_Candidate(alpha, extension.residual_norm, extension.norm, next_vertex))
    if candidates:
        candidates.sort(key=cmp_to_key(lambda left, right: _candidate_compare(left, right, gram, vertex.generator)))
        return candidates[0].vertex
    return _isotropic_next_vertex(gram, roots, discarded, vertex.generator, norms)


def _edge_directions(vertex: _Vertex) -> tuple[tuple[_Roots, _Row], ...]:
    if len(vertex.roots) == 1:
        return ()
    cone = Cone(rays=vertex.roots)
    directions: list[tuple[_Roots, _Row]] = []
    for facet in cone.facets():
        included = tuple(root for root in vertex.roots if facet.contains(vector(SageQQ, root)))
        omitted = tuple(root for root in vertex.roots if root not in included)
        if len(omitted) != 1:
            raise ArithmeticError("Allcock adjacency expects a root-cone facet to omit one wall")
        directions.append((included, omitted[0]))
    return tuple(directions)


def _vertex_configuration(
    lattice: HyperbolicLattices.ParentMethods,
    gram: Matrix_integer_dense,
    vertex: _Vertex,
    norms: tuple[int, ...],
) -> CellConfiguration:
    rows = list(vertex.roots) + [vertex.generator]
    roles = [1] * len(vertex.roots) + [2]
    if _q(gram, vertex.generator) == 0:
        adjacent_generators: list[_Row] = []
        for roots, discarded in _edge_directions(vertex):
            adjacent = _edge_step(gram, vertex, roots, discarded, norms)
            adjacent_generators.append(adjacent.generator)
        auxiliary = sorted(set(adjacent_generators))
        rows.extend(auxiliary)
        roles.extend([3] * len(auxiliary))
    return CellConfiguration.from_coordinate_rows(
        lattice,
        tuple(rows),
        tuple(roles),
    )


def _matrix_rows(lattice: HyperbolicLattices.ParentMethods, isometry: LatticeIsometryMethods) -> _Roots:
    action = lattice.Aut()._row_action_matrix(isometry)
    return tuple(tuple(int(entry) for entry in row) for row in action.rows())


def _update_finiteness(gram: Matrix_integer_dense, invariant_basis: Matrix_rational_dense, generator_rows: _Roots) -> Matrix_rational_dense | None:
    generator = matrix(SageZZ, generator_rows)
    if generator.multiplicative_order() == Infinity:
        return None
    difference = invariant_basis * generator - invariant_basis
    kernel = difference.left_kernel_matrix()
    if kernel.nrows() == 0:
        return None
    invariant = kernel * invariant_basis
    restricted = invariant * gram * invariant.transpose()
    _positive, negative = _signature_pair(restricted)
    if negative == 0:
        return None
    return invariant


def _full_root_orbit(root_rows: set[_Row], isometry_rows: set[_Roots]) -> _Roots:
    roots = set(root_rows)
    pending = list(roots)
    matrices = tuple(matrix(SageZZ, rows) for rows in isometry_rows)
    while pending:
        root = pending.pop()
        row = vector(SageZZ, root)
        for generator in matrices:
            image = tuple(int(entry) for entry in row * generator)
            if image not in roots:
                roots.add(image)
                pending.append(image)
    return tuple(sorted(roots))


@unfinished(31)
def edgewalk_fundamental_domain(gram: Matrix_integer_dense) -> EdgewalkRecord:
    """Return Allcock's fundamental-domain record for a Lorentzian Gram matrix."""
    gram = _integer_gram(gram)
    lattice = HyperbolicLattices(ZZ)(Lattices(ZZ)([[int(entry) for entry in row] for row in gram.rows()]))
    norms = tuple(int(norm) for norm in lattice.possible_root_lengths())
    initial = _initial_vertex(lattice, gram)

    if gram.nrows() == 2:
        if len(initial.roots) != 1:
            raise ArithmeticError("a rank-two chamber corner must have one wall")
        return {
            "simple_root_rows": initial.roots,
            "vertices": ((initial.generator, initial.roots),),
            "is_reflective": True,
            "isometry_generator_rows": (),
        }

    vertices = [initial]
    configurations = [_vertex_configuration(lattice, gram, initial, norms)]
    isometries: set[_Roots] = set()
    invariant_basis: Matrix_rational_dense = matrix.identity(SageQQ, gram.nrows())
    identity = tuple(tuple(int(i == j) for j in range(gram.nrows())) for i in range(gram.nrows()))
    position = 0
    while position < len(vertices):
        vertex = vertices[position]
        stabilizer = cell_stabilizer(configurations[position])
        for generator in stabilizer.generators():
            rows = _matrix_rows(lattice, generator)
            if rows != identity:
                isometries.add(rows)
                updated_basis = _update_finiteness(gram, invariant_basis, rows)
                if updated_basis is None:
                    return {
                        "simple_root_rows": (),
                        "vertices": tuple((item.generator, item.roots) for item in vertices),
                        "is_reflective": False,
                        "isometry_generator_rows": tuple(sorted(isometries)),
                    }
                invariant_basis = updated_basis

        for roots, discarded in _edge_directions(vertex):
            neighbor = _edge_step(gram, vertex, roots, discarded, norms)
            neighbor_configuration = _vertex_configuration(lattice, gram, neighbor, norms)
            found = None
            transporter = None
            for index, configuration in enumerate(configurations):
                candidate = cell_transporter(configuration, neighbor_configuration)
                if candidate is not None:
                    found = index
                    transporter = candidate
                    break
            if found is None:
                vertices.append(neighbor)
                configurations.append(neighbor_configuration)
            elif transporter is not None:
                rows = _matrix_rows(lattice, transporter)
                if rows != identity:
                    isometries.add(rows)
                    updated_basis = _update_finiteness(gram, invariant_basis, rows)
                    if updated_basis is None:
                        return {
                            "simple_root_rows": (),
                            "vertices": tuple((item.generator, item.roots) for item in vertices),
                            "is_reflective": False,
                            "isometry_generator_rows": tuple(sorted(isometries)),
                        }
                    invariant_basis = updated_basis
        position += 1

    root_rows = {root for vertex in vertices for root in vertex.roots}
    simple_roots = _full_root_orbit(root_rows, isometries)
    return {
        "simple_root_rows": simple_roots,
        "vertices": tuple((item.generator, item.roots) for item in vertices),
        "is_reflective": True,
        "isometry_generator_rows": tuple(sorted(isometries)),
    }
