"""Exact Lorentzian perfect-cell local backend."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import gcd
from typing import Literal, TypedDict

from dzack_research.preamble.all import ZZ, Lattices
from dzack_research.preamble.categories.lattice_engines import _rational_positive_vector
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.rings.ring_foundation import _engine_element
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ

from sage_indefinite_port.backends.canonization import (
    CellConfiguration,
    _gram_matrix,
)
from sage_indefinite_port.backends.canonization import (
    cell_stabilizer as configuration_stabilizer,
)
from sage_indefinite_port.backends.canonization import (
    cell_transporter as configuration_transporter,
)
from sage_indefinite_port.backends.polyhedral import (
    FacetIncidence,
    _configuration_permutation,
)
from sage_indefinite_port.backends.polyhedral import (
    facet_orbits as configuration_facet_orbits,
)
from sage_indefinite_port.groups.integral_structures import RationalMatrixGroup

type PerfectMode = Literal["total", "isotropic"]


class TraversalObjectRecord(TypedDict):
    EXT: list[list[int]]
    GRP: list[list[int]]


class TraversalAdjacencyData(TypedDict):
    eInc: list[int]
    eBigMat: list[list[int]]


class TraversalAdjacencyRecord(TypedDict):
    iOrb: int
    x: TraversalAdjacencyData


class TraversalRecord(TypedDict):
    x: TraversalObjectRecord
    ListAdj: list[TraversalAdjacencyRecord]


@dataclass(frozen=True)
class LorentzianPerfectCell:
    lattice: Lattices.ParentMethods
    vector_configuration: tuple[Lattices.ElementMethods, ...]
    mode: PerfectMode

    def __post_init__(self) -> None:
        if self.mode not in ("total", "isotropic"):
            raise ValueError(f"unknown Lorentzian perfect-cell mode {self.mode!r}")
        positive, negative = self.lattice.signature_pair()
        if int(positive) != 1 and int(negative) != 1:
            raise ValueError(
                f"a Lorentzian perfect cell needs signature (1,n) or (n,1), but {self.lattice} has signature {self.lattice.signature_pair()}"
            )
        configuration = self.configuration()
        zero = self.lattice.base_ring().zero()
        if self.mode == "isotropic" and any(
            vector.q() != zero for vector in configuration.vectors
        ):
            raise ValueError(
                "an isotropic perfect cell may contain only isotropic vectors"
            )
        sign = _normalizing_sign(self.lattice)
        if self.mode == "total" and any(
            sign * vector.q() < zero for vector in configuration.vectors
        ):
            raise ValueError(
                "a total perfect cell may contain only nonnegative normalized-norm vectors"
            )

    def configuration(self) -> CellConfiguration:
        return CellConfiguration(self.lattice, self.vector_configuration)


@dataclass(frozen=True)
class LorentzianCellAdjacency:
    source: LorentzianPerfectCell
    facet: FacetIncidence
    target: LorentzianPerfectCell
    transporter: LatticeIsometryMethods | None


class LorentzianPerfectLocalBackend:
    def initial_cell(
        self,
        lattice: Lattices.ParentMethods,
        mode: PerfectMode = "total",
    ) -> LorentzianPerfectCell:
        normalized = _normalized_lattice(lattice)
        direction = _positive_direction(normalized)
        normalized_vectors = _search_initial_vectors(normalized, direction, mode)
        gram = _gram_matrix(normalized).change_ring(SageQQ)
        central = vector(SageQQ, direction)
        first = _coordinate_row(normalized, normalized_vectors[0], SageQQ)
        scalar = (central * gram * first.column())[0]
        base_normal = vector(SageQQ, [-scalar, *(gram * central.column()).column(0)])
        while matrix(
            SageQQ,
            [_coordinate_row(normalized, item, SageQQ) for item in normalized_vectors],
        ).rank() < int(lattice.module_rank()):
            expanded = matrix(
                SageQQ,
                [
                    [1, *_coordinate_row(normalized, item, SageQQ)]
                    for item in normalized_vectors
                ],
            )
            nullspace = expanded.right_kernel_matrix()
            directions = (
                nullspace.matrix_from_columns(range(1, nullspace.ncols()))
                * gram.inverse()
            )
            normal_direction = _negative_direction_in_span(gram, directions)
            normalized_vectors, base_normal, _test_direction, _max_scal = (
                _kernel_flipping(
                    normalized, normalized_vectors, base_normal, normal_direction, mode
                )
            )
        vectors = tuple(
            _same_coordinates(normalized, lattice, item) for item in normalized_vectors
        )
        return LorentzianPerfectCell(lattice, vectors, mode)

    def cell_stabilizer(self, cell: LorentzianPerfectCell) -> RationalMatrixGroup:
        return configuration_stabilizer(cell.configuration())

    def facet_orbits(
        self, cell: LorentzianPerfectCell
    ) -> tuple[tuple[FacetIncidence, ...], ...]:
        return configuration_facet_orbits(
            cell.configuration(), self.cell_stabilizer(cell)
        )

    def flip_across(
        self, cell: LorentzianPerfectCell, facet: FacetIncidence
    ) -> LorentzianPerfectCell:
        normalized = _normalized_lattice(cell.lattice)
        vectors = tuple(
            _same_coordinates(cell.lattice, normalized, item)
            for item in cell.vector_configuration
        )
        if not facet or any(
            position < 0 or position >= len(vectors) for position in facet
        ):
            raise ValueError("a perfect-cell facet must be a nonempty incidence subset")
        nonincident = next(
            (position for position in range(len(vectors)) if position not in facet),
            None,
        )
        if nonincident is None:
            raise ValueError("a facet must omit at least one configuration vector")
        selected = matrix(
            SageQQ,
            [
                _coordinate_row(normalized, vectors[position], SageQQ)
                for position in sorted(facet)
            ],
        )
        facet_nullspace = selected.right_kernel_matrix()
        if facet_nullspace.nrows() != 1:
            raise ArithmeticError("a perfect-cell facet must have one normal direction")
        direction = facet_nullspace.row(0)
        outside = _coordinate_row(normalized, vectors[nonincident], SageQQ)
        if direction.dot_product(outside) < 0:
            direction = -direction
        normal_direction = vector(SageQQ, [0, *direction])
        expanded = matrix(
            SageQQ,
            [[1, *_coordinate_row(normalized, item, SageQQ)] for item in vectors],
        )
        base_nullspace = expanded.right_kernel_matrix()
        if base_nullspace.nrows() != 1:
            raise ArithmeticError(
                "a perfect cell must have a unique affine supporting functional"
            )
        base_normal = base_nullspace.row(0)
        if vector(SageQQ, base_normal[1:]).dot_product(outside) <= 0:
            base_normal = -base_normal
        critical = tuple(vectors[position] for position in sorted(facet))
        flipped, _normal, _direction, _max_scal = _kernel_flipping(
            normalized, critical, base_normal, normal_direction, cell.mode
        )
        returned = tuple(
            _same_coordinates(normalized, cell.lattice, item) for item in flipped
        )
        return LorentzianPerfectCell(cell.lattice, returned, cell.mode)

    def cell_transporter(
        self, source: LorentzianPerfectCell, target: LorentzianPerfectCell
    ) -> LatticeIsometryMethods | None:
        if source.mode != target.mode:
            return None
        return configuration_transporter(source.configuration(), target.configuration())


class LorentzianPerfectComplex:
    r"""Complete quotient traversal of one Lorentzian perfect-domain complex.

    A quotient adjacency stores a selected target representative and the exact
    isometry carrying it onto the actual neighbor obtained by flipping the
    retained source facet.
    """

    def __init__(
        self,
        lattice: Lattices.ParentMethods,
        mode: PerfectMode = "total",
        *,
        local_backend: LorentzianPerfectLocalBackend | None = None,
    ) -> None:
        self._lattice = lattice
        self._mode = mode
        self._local_backend = local_backend or LorentzianPerfectLocalBackend()

    def lattice(self) -> Lattices.ParentMethods:
        return self._lattice

    def mode(self) -> PerfectMode:
        return self._mode

    def local_backend(self) -> LorentzianPerfectLocalBackend:
        return self._local_backend

    @cached_property
    def _traversal(
        self,
    ) -> tuple[tuple[LorentzianPerfectCell, ...], tuple[LorentzianCellAdjacency, ...]]:
        backend = self.local_backend()
        representatives: list[LorentzianPerfectCell] = [
            backend.initial_cell(self.lattice(), self.mode())
        ]
        adjacencies: list[LorentzianCellAdjacency] = []
        position = 0
        while position < len(representatives):
            source = representatives[position]
            for orbit in backend.facet_orbits(source):
                facet = orbit[0]
                neighbor = backend.flip_across(source, facet)
                target = None
                transporter = None
                for representative in representatives:
                    candidate = backend.cell_transporter(representative, neighbor)
                    if candidate is not None:
                        target = representative
                        transporter = candidate
                        break
                if target is None:
                    target = neighbor
                    representatives.append(target)
                    transporter = self.lattice().Aut().identity()
                if transporter is None:
                    raise ArithmeticError(
                        "a quotient-cell adjacency has no transporter to its selected target representative"
                    )
                adjacencies.append(
                    LorentzianCellAdjacency(source, facet, target, transporter)
                )
            position += 1
        return tuple(representatives), tuple(adjacencies)

    def quotient_cells(self) -> tuple[LorentzianPerfectCell, ...]:
        return self._traversal[0]

    def adjacencies(self) -> tuple[LorentzianCellAdjacency, ...]:
        return self._traversal[1]

    def component_preserving_group(self) -> RationalMatrixGroup:
        generators: list[LatticeIsometryMethods] = []
        for cell in self.quotient_cells():
            generators.extend(self.local_backend().cell_stabilizer(cell).generators())
        generators.extend(
            adjacency.transporter
            for adjacency in self.adjacencies()
            if adjacency.transporter is not None
        )
        return RationalMatrixGroup(self.lattice(), tuple(generators))

    def full_orthogonal_group(self) -> RationalMatrixGroup:
        component = self.component_preserving_group()
        negation = self.lattice().Aut()(
            {
                label: -self.lattice().module_generator(label)
                for label in self.lattice().module_generating_set()
            }
        )
        return RationalMatrixGroup(self.lattice(), component.generators() + (negation,))

    def records(self) -> tuple[TraversalRecord, ...]:
        cells = self.quotient_cells()
        positions = {id(cell): position for position, cell in enumerate(cells)}
        by_source: dict[int, list[TraversalAdjacencyRecord]] = {
            position: [] for position in range(len(cells))
        }
        for adjacency in self.adjacencies():
            source_position = positions[id(adjacency.source)]
            target_position = positions[id(adjacency.target)]
            transporter = adjacency.transporter
            if transporter is None:
                raise ArithmeticError(
                    "a completed quotient adjacency has no transporter"
                )
            row_action = self.lattice().Aut()._row_action_matrix(transporter)
            incidence = [
                int(position in adjacency.facet)
                for position in range(len(adjacency.source.vector_configuration))
            ]
            by_source[source_position].append(
                {
                    "iOrb": target_position,
                    "x": {
                        "eInc": incidence,
                        "eBigMat": [
                            [int(entry) for entry in row] for row in row_action.rows()
                        ],
                    },
                }
            )
        records: list[TraversalRecord] = []
        for position, cell in enumerate(cells):
            configuration = cell.configuration()
            stabilizer = self.local_backend().cell_stabilizer(cell)
            permutations = [
                list(_configuration_permutation(configuration, generator))
                for generator in stabilizer.generators()
            ]
            records.append(
                {
                    "x": {
                        "EXT": [
                            [
                                int(entry)
                                for entry in _coordinate_row(
                                    self.lattice(), vector, SageZZ
                                )
                            ]
                            for vector in cell.vector_configuration
                        ],
                        "GRP": permutations,
                    },
                    "ListAdj": by_source[position],
                }
            )
        return tuple(records)


def perfect_domain_traversal(
    gram_rows, option: PerfectMode = "total"
) -> tuple[TraversalRecord, ...]:
    r"""Return complete native perfect-domain traversal records for the Gram rows."""
    lattice = Lattices(ZZ)(gram_rows)
    return LorentzianPerfectComplex(lattice, option).records()


class NoLocalMarkTheoremError(RuntimeError):
    """The local perfect-cell mark theorem is unavailable for this norm."""


class MarkedCellOrbitAlgorithm:
    """Global isotropic-vertex orbits from local cell-stabilizer orbits."""

    def __init__(self, complex_: LorentzianPerfectComplex, norm=0) -> None:
        if norm != 0:
            raise NoLocalMarkTheoremError(
                "the Lorentzian marked-cell source proves finite local marks only for norm zero"
            )
        self._complex = complex_

    def complex(self) -> LorentzianPerfectComplex:
        return self._complex

    def local_marks(
        self, cell: LorentzianPerfectCell
    ) -> tuple[Lattices.ElementMethods, ...]:
        zero = cell.lattice.base_ring().zero()
        return tuple(
            vector for vector in cell.vector_configuration if vector.q() == zero
        )

    def local_stabilizer_orbits(
        self, cell: LorentzianPerfectCell
    ) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        marks = self.local_marks(cell)
        stabilizer = self.complex().local_backend().cell_stabilizer(cell)
        generators = stabilizer.generators_and_inverses()
        unseen = set(range(len(marks)))
        orbits: list[tuple[Lattices.ElementMethods, ...]] = []
        while unseen:
            start = min(unseen)
            pending = [start]
            orbit_positions: set[int] = set()
            while pending:
                position = pending.pop()
                if position in orbit_positions:
                    continue
                orbit_positions.add(position)
                unseen.discard(position)
                mark = marks[position]
                for generator in generators:
                    moved = generator(mark)
                    matches = [
                        index
                        for index, candidate in enumerate(marks)
                        if candidate == moved
                    ]
                    if len(matches) != 1:
                        raise ArithmeticError(
                            "a cell stabilizer does not permute the isotropic local marks"
                        )
                    if matches[0] not in orbit_positions:
                        pending.append(matches[0])
            orbits.append(tuple(marks[index] for index in sorted(orbit_positions)))
        return tuple(orbits)

    def transport_marks(
        self, adjacency: LorentzianCellAdjacency
    ) -> tuple[tuple[Lattices.ElementMethods, Lattices.ElementMethods], ...]:
        transporter = adjacency.transporter
        if transporter is None:
            raise ArithmeticError("a marked quotient adjacency needs a transporter")
        inverse = ~transporter
        target_marks = self.local_marks(adjacency.target)
        transported = []
        for position in sorted(adjacency.facet):
            source_mark = adjacency.source.vector_configuration[position]
            if source_mark.q() != source_mark.parent().base_ring().zero():
                continue
            target_mark = inverse(source_mark)
            matches = [
                candidate for candidate in target_marks if candidate == target_mark
            ]
            if len(matches) != 1:
                raise ArithmeticError(
                    "an isotropic facet mark does not transport to a unique mark of the target representative"
                )
            transported.append((source_mark, matches[0]))
        return tuple(transported)

    def global_orbits(self) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        cells = self.complex().quotient_cells()
        local = tuple(self.local_stabilizer_orbits(cell) for cell in cells)
        nodes = [
            (cell_position, orbit_position)
            for cell_position, cell_orbits in enumerate(local)
            for orbit_position in range(len(cell_orbits))
        ]
        graph = {node: set() for node in nodes}

        def local_orbit_position(cell_position, mark):
            matches = [
                orbit_position
                for orbit_position, orbit in enumerate(local[cell_position])
                if any(candidate == mark for candidate in orbit)
            ]
            if len(matches) != 1:
                raise ArithmeticError(
                    "a transported isotropic mark does not belong to a unique local stabilizer orbit"
                )
            return matches[0]

        positions = {id(cell): position for position, cell in enumerate(cells)}
        for adjacency in self.complex().adjacencies():
            source_position = positions[id(adjacency.source)]
            target_position = positions[id(adjacency.target)]
            for source_mark, target_mark in self.transport_marks(adjacency):
                source_node = (
                    source_position,
                    local_orbit_position(source_position, source_mark),
                )
                target_node = (
                    target_position,
                    local_orbit_position(target_position, target_mark),
                )
                graph[source_node].add(target_node)
                graph[target_node].add(source_node)

        unseen = set(nodes)
        global_orbits: list[tuple[Lattices.ElementMethods, ...]] = []
        while unseen:
            start = min(unseen)
            pending = [start]
            component = []
            while pending:
                node = pending.pop()
                if node not in unseen:
                    continue
                unseen.remove(node)
                component.append(node)
                pending.extend(graph[node])
            marks: list[Lattices.ElementMethods] = []
            for cell_position, orbit_position in sorted(component):
                for mark in local[cell_position][orbit_position]:
                    if not any(existing == mark for existing in marks):
                        marks.append(mark)
            global_orbits.append(tuple(marks))
        return tuple(global_orbits)


def _coordinate_row(lattice, element, ring):
    coordinates = element.to_vector()
    base_ring = lattice.base_ring()
    return vector(
        ring,
        tuple(
            ring(_engine_element(base_ring, coordinates(label)))
            for label in lattice.module_generating_set()
        ),
    )


def _normalizing_sign(lattice):
    positive, negative = lattice.signature_pair()
    match int(positive), int(negative):
        case 1, _:
            return lattice.base_ring().one()
        case _, 1:
            return -lattice.base_ring().one()
        case _:
            raise ValueError(
                f"{lattice} is not Lorentzian: its signature is {lattice.signature_pair()}"
            )


def _normalized_lattice(lattice):
    sign = _normalizing_sign(lattice)
    return lattice if sign == lattice.base_ring().one() else lattice.twist(sign)


def _same_coordinates(source, target, element):
    coordinates = element.to_vector()
    return target(
        tuple(
            target.base_ring()(coordinates(label))
            for label in source.module_generating_set()
        )
    )


def _ambient_element(lattice, coordinates):
    ring = lattice.base_ring()
    return lattice(tuple(ring(int(SageZZ(entry))) for entry in coordinates))


def _positive_direction(lattice) -> tuple[SageQQ, ...]:
    raw = _rational_positive_vector(lattice.gram_tensor())
    rationals = raw.base_ring()
    coordinates = tuple(SageQQ(_engine_element(rationals, entry)) for entry in raw)
    primitive, _scale = _primitive_integral_direction(coordinates)
    return tuple(SageQQ(entry) for entry in primitive)


def _primitive_integral_direction(coordinates):
    denominator = SageZZ.one()
    for entry in coordinates:
        denominator = denominator.lcm(SageZZ(SageQQ(entry).denominator()))
    integers = [SageZZ(SageQQ(entry) * denominator) for entry in coordinates]
    content = SageZZ.zero()
    for value in integers:
        content = SageZZ(gcd(int(content), abs(int(value))))
    if content == 0:
        raise ValueError("a Lorentzian direction cannot be zero")
    return tuple(value // content for value in integers), SageQQ(denominator) / SageQQ(
        content
    )


def _bezout_partner(lattice, timelike):
    ring = lattice.base_ring()
    labels = tuple(lattice.module_generating_set())
    pairings = tuple(
        lattice.b(timelike, lattice.module_generator(label)) for label in labels
    )
    gcd_value = ring.zero()
    coefficients = []
    for pairing in pairings:
        new_gcd, old_coefficient, new_coefficient = gcd_value.xgcd(pairing)
        coefficients = [old_coefficient * coefficient for coefficient in coefficients]
        coefficients.append(new_coefficient)
        gcd_value = new_gcd
    divisibility = timelike.div()
    if gcd_value != divisibility:
        coefficients = [-coefficient for coefficient in coefficients]
        gcd_value = -gcd_value
    if gcd_value != divisibility:
        raise ArithmeticError(
            "Bezout coefficients do not realize the timelike divisibility"
        )
    partner = lattice.linear_combination(
        {
            label: coefficient
            for label, coefficient in zip(labels, coefficients, strict=True)
            if coefficient
        }
    )
    if lattice.b(timelike, partner) != divisibility:
        raise ArithmeticError("the Bezout partner has the wrong pairing")
    return partner


def _find_positive_vectors(
    lattice, rational_direction, max_scal, mode: PerfectMode, *, only_shortest: bool
):
    direction_q = tuple(SageQQ(entry) for entry in rational_direction)
    primitive, scale = _primitive_integral_direction(direction_q)
    timelike = _ambient_element(lattice, primitive)
    if timelike.q() <= lattice.base_ring().zero():
        raise ValueError("positive-vector enumeration needs a timelike direction")
    scaled_max = SageQQ(max_scal) * scale
    divisibility = timelike.div()
    partner = _bezout_partner(lattice, timelike)
    complement = timelike.orthogonal_complement()
    if not complement.is_negative_definite():
        raise ArithmeticError(
            "a timelike vector must have negative-definite orthogonal complement"
        )
    inclusion = complement.inclusion()
    ambient_rows = matrix(
        SageQQ,
        [
            _coordinate_row(lattice, inclusion(generator), SageQQ)
            for generator in complement.module_generators()
        ],
    )
    partner_row = _coordinate_row(lattice, partner, SageQQ)
    timelike_row = _coordinate_row(lattice, timelike, SageQQ)
    square = SageQQ(_engine_element(lattice.base_ring(), timelike.q()))
    d = SageQQ(_engine_element(lattice.base_ring(), divisibility))
    perpendicular_part = partner_row - (d / square) * timelike_row
    perpendicular_coordinates = ambient_rows.transpose().solve_right(
        perpendicular_part.column()
    )
    result = []
    multiplier = 1
    while True:
        level = SageQQ(multiplier) * d
        if scaled_max > 0 and level > scaled_max:
            break
        target = tuple(
            -SageQQ(multiplier) * perpendicular_coordinates[index, 0]
            for index in range(perpendicular_coordinates.nrows())
        )
        bound = -(SageQQ(multiplier) ** 2) * d * d / square
        field = complement.base_ring().fraction_field()
        owned_target = tuple(
            field(int(entry.numerator())) / field(int(entry.denominator()))
            for entry in target
        )
        owned_bound = field(int(bound.numerator())) / field(int(bound.denominator()))
        close = complement.close_vectors(owned_target, owned_bound)
        for perpendicular in close.index_set():
            if mode == "isotropic" and close[perpendicular] != owned_bound:
                continue
            candidate = lattice.scalar_multiple(
                lattice.base_ring()(multiplier), partner
            ) + inclusion(perpendicular)
            if mode == "isotropic" and candidate.q() != lattice.base_ring().zero():
                raise ArithmeticError(
                    "an isotropic shell returned a nonisotropic vector"
                )
            if mode == "total" and candidate.q() < lattice.base_ring().zero():
                raise ArithmeticError("a total shell returned a negative-norm vector")
            result.append(candidate)
        if only_shortest and result:
            break
        multiplier += 1
    return tuple(result)


def _search_initial_vectors(lattice, direction, mode: PerfectMode):
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    direction_row = vector(SageQQ, direction)
    max_scal = (direction_row * gram * direction_row.column())[0]
    while True:
        vectors = _find_positive_vectors(
            lattice, direction, max_scal, mode, only_shortest=True
        )
        if vectors:
            return vectors
        max_scal *= 2


def _upper_bound(gram, base_normal, direction_normal):
    inverse = gram.inverse()
    base_constant = SageQQ(base_normal[0])
    direction_constant = SageQQ(direction_normal[0])
    base_covector = vector(SageQQ, base_normal[1:])
    direction_covector = vector(SageQQ, direction_normal[1:])
    bounds = []
    constant_bound = None
    if direction_constant > 0:
        constant_bound = -base_constant / direction_constant
        if constant_bound <= 0:
            raise ArithmeticError(
                "the affine-normal constant gives a nonpositive flip bound"
            )
        bounds.append(constant_bound)
    shift = SageQQ.one()
    while True:
        covector = base_covector + shift * direction_covector
        value = inverse * covector.column()
        square = (value.transpose() * gram * value)[0, 0]
        if square < 0:
            bounds.append(shift)
            break
        shift *= 2
    base_vector = inverse * base_covector.column()
    direction_vector = inverse * direction_covector.column()
    basis = matrix(SageQQ, [base_vector.column(0), direction_vector.column(0)])
    restricted = basis * gram * basis.transpose()
    a, b, c = restricted[0, 0], restricted[0, 1], restricted[1, 1]
    discriminant = b * b - a * c
    isotropic_bound = None
    if discriminant > 0 and SageQQ(discriminant).is_square():
        root = SageQQ(discriminant).sqrt()
        candidates = (
            ((-b + root) / a, (-b - root) / a)
            if a != 0
            else ((SageQQ(-c) / (2 * b),) if b != 0 else ())
        )
        positive = [value for value in candidates if value > 0]
        isotropic_bound = min(positive) if positive else None
    if isotropic_bound is not None:
        bounds.append(isotropic_bound)
    if (
        constant_bound is not None
        and isotropic_bound is not None
        and constant_bound == isotropic_bound
    ):
        return None
    return min(bounds)


def _negative_direction_in_span(gram, spanning_rows):
    restricted = spanning_rows * gram * spanning_rows.transpose()
    diagonal, change = QuadraticForm(SageQQ, 2 * restricted).rational_diagonal_form(
        return_matrix=True
    )
    negative = [
        index
        for index in range(diagonal.matrix().nrows())
        if diagonal.matrix()[index, index] < 0
    ]
    if not negative:
        raise ArithmeticError("the perfect-cell nullspace has no negative direction")
    ambient = spanning_rows.transpose() * change.column(negative[0])
    return vector(SageQQ, [0, *ambient])


def _vector_rows(vectors):
    return {
        tuple(int(entry) for entry in _coordinate_row(item.parent(), item, SageZZ))
        for item in vectors
    }


def _kernel_flipping(
    lattice, critical, base_normal, direction_normal, mode: PerfectMode
):
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    upper = _upper_bound(gram, base_normal, direction_normal)
    if upper is None:
        raise ArithmeticError(
            "the flip direction ends at a forbidden isotropic boundary"
        )
    lower = SageQQ.zero()
    inverse = gram.inverse()
    critical_first = _coordinate_row(lattice, critical[0], SageQQ)
    critical_rows = _vector_rows(critical)
    total = ()
    while True:
        middle = (lower + upper) / 2
        normal = base_normal + middle * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = test_direction.column(0)
        square = (test_direction.transpose() * gram * test_direction)[0, 0]
        max_scal = (critical_first * gram * test_direction)[0]
        if square <= 0 or max_scal <= 0:
            upper = middle
            continue
        total = _find_positive_vectors(
            lattice, tuple(test_direction_row), max_scal, mode, only_shortest=True
        )
        total_rows = _vector_rows(total)
        if total_rows == critical_rows:
            lower = middle
            continue
        if critical_rows.issubset(total_rows):
            return total, normal, test_direction_row, max_scal
        break
    while True:
        if not total:
            raise ArithmeticError(
                "perfect-cell flipping reached an empty positive-vector shell"
            )
        expanded = vector(SageQQ, [1, *_coordinate_row(lattice, total[0], SageQQ)])
        denominator = direction_normal.dot_product(expanded)
        if denominator == 0:
            raise ArithmeticError(
                "the flip direction is parallel to a newly found wall"
            )
        shift = -base_normal.dot_product(expanded) / denominator
        normal = base_normal + shift * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = test_direction.column(0)
        max_scal = (critical_first * gram * test_direction)[0]
        total = _find_positive_vectors(
            lattice, tuple(test_direction_row), max_scal, mode, only_shortest=True
        )
        if critical_rows.issubset(_vector_rows(total)):
            return total, normal, test_direction_row, max_scal
