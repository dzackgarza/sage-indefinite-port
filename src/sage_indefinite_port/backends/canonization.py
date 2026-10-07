"""Exact configuration canonization through Sage's Bliss backend.

Bliss only proposes configuration permutations. Every returned map is
reconstructed over QQ and verified over ZZ against the full configuration and
the ambient bilinear forms.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache

from dzack_research.preamble.all import QQ, ZZ, Lattices, Modules, RestrictedScalarsModules
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.module_morphisms.module_morphisms import ModuleEmbeddingMethods
from dzack_research.preamble.categories.rings.ring_foundation import _engine_element
from sage.graphs.graph import Graph
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ

from sage_indefinite_port.groups.integral_structures import (
    FiniteIntegralRepresentation,
    IntegralStructureAction,
    RationalMatrixGroup,
)

type CoordinateRows = tuple[tuple[int | Integer | Rational, ...], ...]


@dataclass(frozen=True)
class CellConfiguration:
    """A finite full-rank vector configuration in one integral lattice."""

    lattice: Lattices.ParentMethods
    vectors: tuple[Lattices.ElementMethods, ...]
    roles: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if not self.vectors:
            raise ValueError("a cell configuration must contain at least one vector")
        if any(vector.parent() is not self.lattice for vector in self.vectors):
            raise ValueError("every configuration vector must lie in the configured lattice")
        if self.roles and len(self.roles) != len(self.vectors):
            raise ValueError("configuration roles must label every configuration vector")
        if _coordinate_matrix(self).rank() != int(self.lattice.module_rank()):
            raise ValueError("a cell configuration must span its ambient lattice over QQ")

    @classmethod
    def from_coordinate_rows(
        cls,
        lattice: Lattices.ParentMethods,
        rows: tuple[tuple[int, ...], ...],
        roles: tuple[int, ...] = (),
    ) -> CellConfiguration:
        """Construct a configuration from integral coordinate rows."""
        labels = tuple(lattice.module_generating_set())
        if any(len(row) != len(labels) for row in rows):
            raise ValueError("configuration rows must have the ambient lattice rank")
        vectors = tuple(lattice.linear_combination({label: lattice.base_ring()(int(entry)) for label, entry in zip(labels, row, strict=True) if entry}) for row in rows)
        return cls(lattice, vectors, roles)


@dataclass(frozen=True)
class _BlissCanonization:
    canonical_signature: tuple[
        tuple[tuple[int, int, int], ...],
        tuple[tuple[int, int], ...],
    ]
    labeling: tuple[int, ...]


def presentation_bucket_key(
    lattice: Lattices.ParentMethods,
) -> tuple[int, tuple[int, int], bool, int]:
    """Return cheap necessary isometry invariants for presentation bucketing."""
    positive, negative = lattice.signature_pair()
    return (
        int(lattice.module_rank()),
        (int(positive), int(negative)),
        bool(lattice.is_even()),
        int(_engine_element(lattice.base_ring(), lattice.discriminant())),
    )


@cache
def _configuration_integral_structure_action(
    configuration: CellConfiguration,
) -> tuple[IntegralStructureAction, Lattices.ParentMethods, RestrictedScalarsModules.ParentMethods, ModuleEmbeddingMethods]:
    canon = _bliss_canonization(configuration)
    rational_group = _configuration_rational_automorphism_group(configuration, canon)
    rational_lattice = configuration.lattice.base_change(ZZ.fraction_field_map())
    restriction = Modules(QQ).restriction_of_scalars(ZZ.Mor(QQ)(lambda element: QQ(element)))
    integral_structure_space = restriction(rational_lattice)
    inclusion = _rational_lattice_with_integral_structure(
        configuration.lattice,
        rational_lattice,
        integral_structure_space,
        _identity_rows(int(configuration.lattice.module_rank())),
    )
    if rational_group.rational_lattice() is not rational_lattice:
        raise ArithmeticError("configuration automorphisms and integral structure use different rational lattices")
    return IntegralStructureAction(rational_group, inclusion), rational_lattice, integral_structure_space, inclusion


@cache
def cell_transporter(
    source: CellConfiguration,
    target: CellConfiguration,
) -> LatticeIsometryMethods | None:
    """Return a verified integral isometry carrying source to target."""
    if presentation_bucket_key(source.lattice) != presentation_bucket_key(target.lattice):
        return None
    if len(source.vectors) != len(target.vectors):
        return None

    source_canon = _bliss_canonization(source)
    target_canon = _bliss_canonization(target)
    if source_canon.canonical_signature != target_canon.canonical_signature:
        return None

    target_by_canonical = {target_canon.labeling[position]: position for position in range(len(target.vectors))}
    try:
        matching = tuple(target_by_canonical[source_canon.labeling[position]] for position in range(len(source.vectors)))
    except KeyError as exc:
        raise ArithmeticError("Bliss mapped a vector vertex outside the vector color classes") from exc

    base_action = _configuration_row_action(source, target, matching)
    verified = _verified_integral_isometry(source, target, matching, base_action)
    if verified is not None:
        return verified

    target_action, target_rational, integral_structure_space, target_inclusion = _configuration_integral_structure_action(target)

    source_image = _rational_lattice_with_integral_structure(
        target.lattice,
        target_rational,
        integral_structure_space,
        _rational_rows(base_action),
    )
    correction = target_action.transporter(
        source_image,
        target_inclusion,
    )
    if correction is None:
        return None

    correction_action = target_rational.Aut()._row_action_matrix(correction)
    corrected = base_action * correction_action
    target_coordinates = _coordinate_matrix(target)
    target_position = {tuple(row): position for position, row in enumerate(target_coordinates.rows())}
    corrected_permutation: list[int] = []
    for row in (matrix(SageQQ, target_coordinates) * correction_action).rows():
        if any(entry.denominator() != 1 for entry in row):
            raise ArithmeticError("an integral-structure correction moves a configuration vector to nonintegral coordinates")
        key = tuple(SageZZ(entry) for entry in row)
        position = target_position.get(key)
        if position is None:
            raise ArithmeticError("an integral-structure correction does not preserve the target configuration")
        corrected_permutation.append(position)
    corrected_matching = tuple(corrected_permutation[position] for position in matching)
    return _verified_integral_isometry(
        source,
        target,
        corrected_matching,
        corrected,
    )


def cell_stabilizer(configuration: CellConfiguration) -> RationalMatrixGroup:
    """Return the integral isometries preserving the configuration setwise."""
    graph, partition, _colors = _configuration_graph(configuration)
    permutation_group = graph.automorphism_group(
        partition=partition,
        algorithm="bliss",
        edge_labels=True,
    )
    permutation_generators = tuple(permutation_group.gens())
    rational_actions = tuple(
        _configuration_row_action(
            configuration,
            configuration,
            tuple(int(generator(position)) for position in range(graph.order())),
        )
        for generator in permutation_generators
    )
    rank = int(configuration.lattice.module_rank())
    invariant_basis = matrix(SageQQ, _identity_rows(rank))
    action_matrices = rational_actions + tuple(action.inverse() for action in rational_actions)
    while True:
        inverse_basis = invariant_basis.inverse()
        spanning_rows = list(invariant_basis.rows())
        stable = True
        for action_matrix in action_matrices:
            moved = invariant_basis * action_matrix
            spanning_rows.extend(moved.rows())
            if any(entry.denominator() != 1 for row in moved.rows() for entry in row * inverse_basis):
                stable = False
        if stable:
            break
        denominator = SageZZ.one()
        for row in spanning_rows:
            for entry in row:
                denominator = denominator.lcm(entry.denominator())
        integral_rows = tuple(tuple(SageZZ(denominator * entry) for entry in row) for row in spanning_rows)
        integral_basis = matrix(SageZZ, integral_rows).row_module().basis_matrix()
        invariant_basis = matrix(
            SageQQ,
            [[SageQQ(entry) / SageQQ(denominator) for entry in row] for row in integral_basis.rows()],
        )

    restricted_actions: list[Matrix_integer_dense] = []
    for action_matrix in rational_actions:
        restricted = invariant_basis * action_matrix * invariant_basis.inverse()
        if any(entry.denominator() != 1 for entry in restricted.list()):
            raise ArithmeticError("a configuration automorphism does not preserve its invariant over-lattice")
        restricted_actions.append(
            matrix(
                SageZZ,
                restricted.nrows(),
                restricted.ncols(),
                tuple(SageZZ(entry) for entry in restricted.list()),
            )
        )

    selected_coordinates = (
        matrix(
            SageQQ,
            _identity_rows(rank),
        )
        * invariant_basis.inverse()
    )
    if any(entry.denominator() != 1 for entry in selected_coordinates.list()):
        raise ArithmeticError("the configured lattice is not integral in its invariant over-lattice")
    selected_key = FiniteIntegralRepresentation._row_lattice_key(tuple(tuple(SageZZ(entry) for entry in row) for row in selected_coordinates.rows()))
    orbit_keys = [selected_key]
    orbit_position = {selected_key: 0}
    frontier = [0]
    while frontier:
        source_position = frontier.pop()
        source_key = orbit_keys[source_position]
        for restricted_action in restricted_actions:
            candidate_key = FiniteIntegralRepresentation._row_lattice_key((matrix(SageZZ, source_key) * restricted_action).rows())
            if candidate_key not in orbit_position:
                orbit_position[candidate_key] = len(orbit_keys)
                orbit_keys.append(candidate_key)
                frontier.append(len(orbit_keys) - 1)

    finite_generators: list[GapElement] = []
    for restricted_action in restricted_actions:
        images: list[int] = []
        for key in orbit_keys:
            image_key = FiniteIntegralRepresentation._row_lattice_key((matrix(SageZZ, key) * restricted_action).rows())
            images.append(orbit_position[image_key] + 1)
        finite_generators.append(libgap.PermList(images))
    graph_gap_generators = tuple(libgap.PermList([int(generator(position)) + 1 for position in range(graph.order())]) for generator in permutation_generators)
    graph_gap_group = libgap.Group(graph_gap_generators)
    finite_group = libgap.Group(finite_generators)
    orbit_homomorphism = libgap.GroupHomomorphismByImages(
        graph_gap_group,
        finite_group,
        graph_gap_generators,
        finite_generators,
    )
    point_stabilizer = libgap.Stabilizer(finite_group, 1)
    stabilizer_group = libgap.PreImage(orbit_homomorphism, point_stabilizer)
    stabilizer_permutations = libgap.SmallGeneratingSet(stabilizer_group)

    integral_generators: list[LatticeIsometryMethods] = []
    for permutation in stabilizer_permutations:
        matching = tuple(int(position) - 1 for position in libgap.ListPerm(permutation, graph.order()).sage())
        row_action = _configuration_row_action(
            configuration,
            configuration,
            matching,
        )
        verified = _verified_integral_isometry(
            configuration,
            configuration,
            matching,
            row_action,
        )
        if verified is None:
            raise ArithmeticError("a finite configuration-stabilizer lift is not an integral isometry")
        integral_generators.append(verified)
    return RationalMatrixGroup(
        configuration.lattice,
        tuple(integral_generators),
    )


@cache
def _bliss_canonization(configuration: CellConfiguration) -> _BlissCanonization:
    graph, partition, colors = _configuration_graph(configuration)
    canonical_graph, certificate = graph.canonical_label(
        partition=partition,
        algorithm="bliss",
        certificate=True,
        edge_labels=True,
    )
    labeling = tuple(int(certificate[position]) for position in range(graph.order()))
    canonical_colors = [(0, 0)] * graph.order()
    for position, color in enumerate(colors):
        canonical_colors[labeling[position]] = color
    canonical_edges = tuple(
        sorted(
            (
                min(int(left), int(right)),
                max(int(left), int(right)),
                int(label),
            )
            for left, right, label in canonical_graph.edges(labels=True, sort=False)
        )
    )
    return _BlissCanonization(
        (canonical_edges, tuple(canonical_colors)),
        labeling,
    )


@cache
def _configuration_graph(
    configuration: CellConfiguration,
) -> tuple[Graph, list[list[int]], tuple[tuple[int, int], ...]]:
    coordinates = _coordinate_matrix(configuration)
    pairing = coordinates * _gram_matrix(configuration.lattice) * coordinates.transpose()
    count = len(configuration.vectors)
    roles = configuration.roles or tuple(0 for _ in range(count))
    colors = [(int(roles[position]), int(pairing[position, position])) for position in range(count)]
    edges: list[tuple[int, int, int]] = []

    for left in range(count):
        for right in range(left + 1, count):
            edges.append((left, right, int(pairing[left, right])))

    graph = Graph(count)
    graph.add_edges(edges)
    color_classes: dict[tuple[int, int], list[int]] = {}
    for vertex, color in enumerate(colors):
        color_classes.setdefault(color, []).append(vertex)
    partition = [color_classes[color] for color in sorted(color_classes)]
    return graph, partition, tuple(colors)


@cache
def _coordinate_matrix(configuration: CellConfiguration) -> Matrix_integer_dense:
    labels = tuple(configuration.lattice.module_generating_set())
    base_ring = configuration.lattice.base_ring()
    return matrix(
        SageZZ,
        tuple(tuple(Integer(_engine_element(base_ring, vector.to_vector()(label))) for label in labels) for vector in configuration.vectors),
    )


def _configuration_row_action(
    source: CellConfiguration,
    target: CellConfiguration,
    matching: tuple[int, ...],
) -> Matrix_rational_dense:
    source_coordinates = matrix(SageQQ, _coordinate_matrix(source))
    target_coordinates = matrix(SageQQ, _coordinate_matrix(target))
    mapped_target = matrix(
        SageQQ,
        tuple(target_coordinates.row(position) for position in matching),
    )
    rank = int(source.lattice.module_rank())
    pivots = tuple(int(position) for position in source_coordinates.transpose().pivots())
    if len(pivots) < rank:
        raise ArithmeticError("the source configuration lost full rank during transporter reconstruction")
    basis_rows = pivots[:rank]
    source_basis = matrix(
        SageQQ,
        tuple(source_coordinates.row(position) for position in basis_rows),
    )
    target_basis = matrix(
        SageQQ,
        tuple(mapped_target.row(position) for position in basis_rows),
    )
    action = source_basis.inverse() * target_basis
    if source_coordinates * action != mapped_target:
        raise ArithmeticError("the Bliss configuration bijection is not induced by one linear map")
    return action


def _verified_integral_isometry(
    source: CellConfiguration,
    target: CellConfiguration,
    matching: tuple[int, ...],
    row_action: Matrix_rational_dense,
) -> LatticeIsometryMethods | None:
    if any(entry.denominator() != 1 for entry in row_action.list()):
        return None
    integral = matrix(
        SageZZ,
        row_action.nrows(),
        row_action.ncols(),
        tuple(Integer(entry) for entry in row_action.list()),
    )
    if abs(Integer(integral.determinant())) != 1:
        return None

    source_coordinates = _coordinate_matrix(source)
    target_coordinates = _coordinate_matrix(target)
    mapped_target = matrix(
        SageZZ,
        tuple(target_coordinates.row(position) for position in matching),
    )
    if source_coordinates * integral != mapped_target:
        raise ArithmeticError("an alleged integral transporter does not move the full configuration")
    if integral * _gram_matrix(target.lattice) * integral.transpose() != _gram_matrix(source.lattice):
        raise ArithmeticError("an alleged cell transporter does not preserve the ambient form")

    isometry = source.lattice.Isom(target.lattice)._isometry_from_column_matrix(integral.transpose())
    if not isinstance(isometry, LatticeIsometryMethods):
        raise TypeError("the verified transporter was not realized as a lattice isometry")
    return isometry


@cache
def _gram_matrix(lattice: Lattices.ParentMethods) -> Matrix_integer_dense:
    generators = tuple(lattice.module_generators())
    base_ring = lattice.base_ring()
    return matrix(
        SageZZ,
        tuple(tuple(Integer(_engine_element(base_ring, lattice.b(left, right))) for right in generators) for left in generators),
    )


def _configuration_rational_automorphism_group(
    configuration: CellConfiguration,
    canon: _BlissCanonization,
) -> RationalMatrixGroup:
    rational_lattice = configuration.lattice.base_change(ZZ.fraction_field_map())
    automorphisms: list[LatticeIsometryMethods] = []
    graph, partition, _colors = _configuration_graph(configuration)
    permutation_group = graph.automorphism_group(
        partition=partition,
        algorithm="bliss",
        edge_labels=True,
    )
    for generator in permutation_group.gens():
        permutation = tuple(int(generator(position)) for position in range(graph.order()))
        matching = tuple(permutation[position] for position in range(len(configuration.vectors)))
        if any(position >= len(configuration.vectors) for position in matching):
            raise ArithmeticError("a Bliss generator moved a vector vertex to an incidence vertex")
        action = _configuration_row_action(configuration, configuration, matching)
        automorphisms.append(rational_lattice.Aut()._isometry_from_column_matrix(action.transpose()))
    return RationalMatrixGroup(rational_lattice, tuple(automorphisms))


def _rational_lattice_with_integral_structure(
    lattice: Lattices.ParentMethods,
    rational_lattice: Lattices.ParentMethods,
    space: RestrictedScalarsModules.ParentMethods,
    rows: CoordinateRows,
) -> ModuleEmbeddingMethods:
    rank = int(lattice.module_rank())
    domain = ZZ.free_module(rank)
    domain_labels = tuple(domain.module_generating_set())
    rational_labels = tuple(rational_lattice.module_generating_set())
    inclusion = domain.Mono(space)(
        {
            domain_labels[position]: space.wrap(
                rational_lattice.linear_combination(
                    {
                        label: QQ(int(SageQQ(entry).numerator())) / QQ(int(SageQQ(entry).denominator()))
                        for label, entry in zip(
                            rational_labels,
                            rows[position],
                            strict=True,
                        )
                        if entry
                    }
                )
            )
            for position in range(rank)
        }
    )
    return inclusion


def _identity_rows(rank: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(1 if row == column else 0 for column in range(rank)) for row in range(rank))


def _rational_rows(
    matrix_value: Matrix_rational_dense,
) -> tuple[tuple[Rational, ...], ...]:
    return tuple(tuple(matrix_value[row, column] for column in range(matrix_value.ncols())) for row in range(matrix_value.nrows()))
