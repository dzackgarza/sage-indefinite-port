"""Exact Lorentzian perfect-cell local backend."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cache, cached_property
from math import gcd
from typing import Literal, TypedDict

from dzack_research.preamble.all import ZZ, Lattices
from dzack_research.preamble.categories.definite_lattices import (
    _element_from_coordinates,
    _ExactCVPEngine,
)
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
from sage_indefinite_port.groups.integral_structures import GeneratedSubgroup, IntegralStructureAction, RationalMatrixGroup
from sage_indefinite_port.indefinite.eichler import EichlerOrbitCover, build_eichler_envelope
from sage_indefinite_port.indefinite.vector_sections import NonIsotropicVectorSection, orthogonal_section
from sage_indefinite_port.invariants import AttackProfile, VectorPrefilter

type PerfectMode = Literal["total", "isotropic"]


class IndefiniteOrthogonalAlgorithm:
    def __init__(self) -> None:
        self._lorentzian_complex_cache: dict[
            Lattices.ParentMethods,
            LorentzianPerfectComplex,
        ] = {}

    def attack_profile(self, lattice: Lattices.ParentMethods) -> AttackProfile:
        return AttackProfile.from_lattice(lattice)

    def orthogonal_group(self, lattice: Lattices.ParentMethods):
        profile = self.attack_profile(lattice)
        match profile.positive_index:
            case 0:
                from sage_indefinite_port.backends.definite import definite_orthogonal_group

                return definite_orthogonal_group(profile.signed_view)
            case 1:
                if int(profile.signed_view.module_rank()) == 2:
                    labels = tuple(profile.signed_view.module_generating_set())
                    if len(labels) == 2:
                        left = profile.signed_view.module_generator(labels[0])
                        right = profile.signed_view.module_generator(labels[1])
                        zero = profile.signed_view.base_ring().zero()
                        if profile.signed_view.q(left) == zero and profile.signed_view.q(right) == zero and profile.signed_view.b(left, right) != zero:
                            automorphisms = profile.signed_view.Aut()
                            swap = automorphisms((right, left))
                            negation = automorphisms((-left, -right))
                            return RationalMatrixGroup(profile.signed_view, (swap, negation))
                return LorentzianPerfectComplex(profile.signed_view, "total").full_orthogonal_group()
            case _:
                return self._higher_witt_orthogonal_group(profile.signed_view)

    def _two_u_cover_model(self, lattice: Lattices.ParentMethods) -> EichlerOrbitCover:
        if lattice.splits_two_hyperbolic_planes():
            return EichlerOrbitCover(lattice.two_u_eichler_model_from_represented_biproduct())
        envelope = build_eichler_envelope(lattice)
        return EichlerOrbitCover(envelope.two_u_decomposition.lattice.two_u_eichler_model_from_represented_biproduct())

    def _higher_witt_orthogonal_group(self, lattice: Lattices.ParentMethods) -> GeneratedSubgroup:
        if lattice.splits_two_hyperbolic_planes():
            model = self._two_u_cover_model(lattice)
            vector = model.choose_splitting_vector()
            approximate_family = model.subgroup()
            approximate = tuple(approximate_family[label] for label in approximate_family.index_set())
            stabilizer = self.vector_stabilizer(vector).generators()
            transporters = []
            for candidate in model.covering_representatives(vector.q(), primitive=vector.is_primitive()):
                witness = self.vector_transporter(vector, candidate)
                if witness is not None:
                    transporters.append(witness)
            generators = approximate + stabilizer + tuple(transporters)
            return GeneratedSubgroup(RationalMatrixGroup(lattice, generators), generators)

        envelope = build_eichler_envelope(lattice)
        envelope_group = self._higher_witt_orthogonal_group(envelope.envelope)
        rational_envelope = envelope.integral_action.rational_group().rational_lattice()
        fraction_map = lattice.base_ring().fraction_field_map()
        rational_automorphisms = rational_envelope.Aut()
        rational_generators = tuple(rational_automorphisms(generator.base_change(fraction_map)) for generator in envelope_group.generators())
        rational_group = RationalMatrixGroup(rational_envelope, rational_generators)
        action = IntegralStructureAction(rational_group, envelope.inclusion)
        integral = action.lattice_stabilizer()
        lifted = []
        for generator in integral.generators():
            images = tuple(envelope.lattice_to_envelope.lift(generator(envelope.lattice_to_envelope(source_generator))) for source_generator in lattice.module_generators())
            lifted.append(lattice.Aut()(images))
        generators = tuple(lifted)
        return GeneratedSubgroup(RationalMatrixGroup(lattice, generators), generators)

    def isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods):
        if self.attack_profile(source).positive_index != self.attack_profile(target).positive_index:
            return None
        if source.is_definite() and target.is_definite():
            from sage_indefinite_port.backends.definite import definite_isometry

            return definite_isometry(source, target)
        source_profile = self.attack_profile(source)
        target_profile = self.attack_profile(target)
        if source_profile.positive_index == 1 and target_profile.positive_index == 1:
            signed_source = source_profile.signed_view
            signed_target = target_profile.signed_view
            source_complex = self._lorentzian_complex_cache.get(signed_source)
            if source_complex is None:
                source_complex = LorentzianPerfectComplex(signed_source, "total")
                self._lorentzian_complex_cache[signed_source] = source_complex
            target_complex = self._lorentzian_complex_cache.get(signed_target)
            if target_complex is None:
                target_complex = LorentzianPerfectComplex(signed_target, "total")
                self._lorentzian_complex_cache[signed_target] = target_complex
            return source_complex.isometry_to(target_complex)
        source_model = self._two_u_cover_model(source)
        target_model = self._two_u_cover_model(target)
        witness = source_model.model.isometry_to(target_model.model)
        if witness is None:
            return None
        if witness.domain().gram_tensor() != source.gram_tensor() or witness.codomain().gram_tensor() != target.gram_tensor():
            raise ArithmeticError("the 2U model isometry has the wrong represented endpoints")
        return witness

    def vector_stabilizer(self, vector: Lattices.ElementMethods) -> GeneratedSubgroup:
        section = orthogonal_section(vector)
        if not isinstance(section, NonIsotropicVectorSection):
            raise NotImplementedError("the isotropic stabilizer branch lands with PHASE-T5")
        reduced_group = self.orthogonal_group(section.reduced_object())
        reduced_generators = reduced_group.generators() if isinstance(reduced_group, RationalMatrixGroup) else tuple(reduced_group.group_generators())
        ambient = vector.parent()
        lifted = []
        for generator in reduced_generators:
            torsor = section.rational_lift(generator, target=section)
            if torsor.integral_parameters(ambient, ambient) is None:
                raise ArithmeticError("a reduced stabilizer generator has no integral ambient lift")
            lifted.append(torsor.one_integral_extension())
        generators = tuple(lifted)
        if any(generator(vector) != vector for generator in generators):
            raise ArithmeticError("a lifted vector-stabilizer generator does not fix the selected vector")
        supergroup = RationalMatrixGroup(ambient, generators)
        return GeneratedSubgroup(supergroup, generators)

    def vector_transporter(self, source_vector: Lattices.ElementMethods, target_vector: Lattices.ElementMethods):
        if VectorPrefilter.from_vector(source_vector) != VectorPrefilter.from_vector(target_vector):
            return None
        source_section = orthogonal_section(source_vector)
        target_section = orthogonal_section(target_vector)
        if isinstance(source_section, NonIsotropicVectorSection) != isinstance(target_section, NonIsotropicVectorSection):
            return None
        reduced_source = source_section.reduced_object()
        reduced_target = target_section.reduced_object()
        reduced_isometry = reduced_source.isometry_to(reduced_target)
        if reduced_isometry is None:
            return None
        source = source_vector.parent()
        target = target_vector.parent()
        candidates = [reduced_isometry]
        if reduced_target.is_definite():
            for automorphism in reduced_target.O():
                candidates.append(automorphism * reduced_isometry)
        for candidate in candidates:
            torsor = source_section.rational_lift(candidate, target=target_section)
            if torsor.integral_parameters(source, target) is None:
                continue
            witness = torsor.one_integral_extension()
            if witness(source_vector) != target_vector:
                raise ArithmeticError("an integral vector transporter does not carry the selected source vector to the target")
            return witness
        return None


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
            raise ValueError(f"a Lorentzian perfect cell needs signature (1,n) or (n,1), but {self.lattice} has signature {self.lattice.signature_pair()}")
        gram = _gram_matrix(self.lattice)
        squares = tuple((row * gram * row.column())[0] for row in (_coordinate_row(self.lattice, vector, SageZZ) for vector in self.vector_configuration))
        if self.mode == "isotropic" and any(square != 0 for square in squares):
            raise ValueError("an isotropic perfect cell may contain only isotropic vectors")
        normalized_sign = 1 if int(positive) == 1 else -1
        if self.mode == "total" and any(normalized_sign * square < 0 for square in squares):
            raise ValueError("a total perfect cell may contain only nonnegative normalized-norm vectors")

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
                [[1, *_coordinate_row(normalized, item, SageQQ)] for item in normalized_vectors],
            )
            nullspace = expanded.right_kernel_matrix()
            directions = nullspace.matrix_from_columns(range(1, nullspace.ncols())) * gram.inverse()
            outside_direction = _negative_direction_in_span(gram, directions)
            direction_covector = gram * outside_direction.column()
            critical_row = _coordinate_row(normalized, normalized_vectors[0], SageQQ)
            direction_scalar = (critical_row * direction_covector)[0]
            normal_direction = vector(
                SageQQ,
                [-direction_scalar, *direction_covector.column(0)],
            )
            normalized_vectors, base_normal, _test_direction, _max_scal = _kernel_flipping(normalized, normalized_vectors, base_normal, normal_direction, mode)
        vectors = tuple(_same_coordinates(normalized, lattice, item) for item in normalized_vectors)
        return LorentzianPerfectCell(lattice, vectors, mode)

    @cache
    def cell_stabilizer(self, cell: LorentzianPerfectCell) -> RationalMatrixGroup:
        return configuration_stabilizer(cell.configuration())

    def facet_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[FacetIncidence, ...], ...]:
        return configuration_facet_orbits(cell.configuration(), self.cell_stabilizer(cell))

    def flip_across(self, cell: LorentzianPerfectCell, facet: FacetIncidence) -> LorentzianPerfectCell:
        normalized = _normalized_lattice(cell.lattice)
        vectors = cell.vector_configuration
        if not facet or any(position < 0 or position >= len(vectors) for position in facet):
            raise ValueError("a perfect-cell facet must be a nonempty incidence subset")
        nonincident = next(
            (position for position in range(len(vectors)) if position not in facet),
            None,
        )
        if nonincident is None:
            raise ValueError("a facet must omit at least one configuration vector")
        selected = matrix(
            SageQQ,
            [_coordinate_row(normalized, vectors[position], SageQQ) for position in sorted(facet)],
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
            raise ArithmeticError("a perfect cell must have a unique affine supporting functional")
        base_normal = base_nullspace.row(0)
        if vector(SageQQ, base_normal[1:]).dot_product(outside) <= 0:
            base_normal = -base_normal
        critical = tuple(vectors[position] for position in sorted(facet))
        flipped, _normal, _direction, _max_scal = _kernel_flipping(normalized, critical, base_normal, normal_direction, cell.mode)
        returned = tuple(_same_coordinates(normalized, cell.lattice, item) for item in flipped)
        return LorentzianPerfectCell(cell.lattice, returned, cell.mode)

    def cell_transporter(self, source: LorentzianPerfectCell, target: LorentzianPerfectCell) -> LatticeIsometryMethods | None:
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
        self._traversal_representatives: list[LorentzianPerfectCell] | None = None
        self._traversal_representative_buckets: dict[
            tuple[tuple[int, ...], tuple[int, ...]],
            list[LorentzianPerfectCell],
        ] = {}
        self._traversal_adjacencies: list[LorentzianCellAdjacency] = []
        self._traversal_position = 0

    def lattice(self) -> Lattices.ParentMethods:
        return self._lattice

    def mode(self) -> PerfectMode:
        return self._mode

    def local_backend(self) -> LorentzianPerfectLocalBackend:
        return self._local_backend

    def _ensure_traversal_started(self) -> None:
        if self._traversal_representatives is not None:
            return
        initial = self.local_backend().initial_cell(self.lattice(), self.mode())
        self._traversal_representatives = [initial]
        self._traversal_representative_buckets = {
            _perfect_form_hash_key(initial): [initial],
        }

    def _expand_traversal_once(self) -> tuple[LorentzianPerfectCell, ...]:
        self._ensure_traversal_started()
        representatives = self._traversal_representatives
        if representatives is None or self._traversal_position >= len(representatives):
            return ()
        backend = self.local_backend()
        source = representatives[self._traversal_position]
        new_representatives: list[LorentzianPerfectCell] = []
        for orbit in backend.facet_orbits(source):
            facet = orbit[0]
            neighbor = backend.flip_across(source, facet)
            target = None
            transporter = None
            neighbor_key = _perfect_form_hash_key(neighbor)
            for representative in self._traversal_representative_buckets.get(neighbor_key, ()):
                candidate = backend.cell_transporter(representative, neighbor)
                if candidate is not None:
                    target = representative
                    transporter = candidate
                    break
            if target is None:
                target = neighbor
                representatives.append(target)
                new_representatives.append(target)
                self._traversal_representative_buckets.setdefault(neighbor_key, []).append(target)
                transporter = self.lattice().Aut().identity()
            if transporter is None:
                raise ArithmeticError("a quotient-cell adjacency has no transporter to its selected target representative")
            self._traversal_adjacencies.append(LorentzianCellAdjacency(source, facet, target, transporter))
        self._traversal_position += 1
        return tuple(new_representatives)

    def isometry_to(self, target: LorentzianPerfectComplex):
        if self.mode() != target.mode():
            return None
        self._ensure_traversal_started()
        target._ensure_traversal_started()
        representatives = self._traversal_representatives
        target_representatives = target._traversal_representatives
        if representatives is None or target_representatives is None:
            raise ArithmeticError("a started Lorentzian traversal has no initial cell")
        target_cell = target_representatives[0]
        target_key = _perfect_form_hash_key(target_cell)
        backend = self.local_backend()
        for source_cell in self._traversal_representative_buckets.get(target_key, ()):
            witness = backend.cell_transporter(source_cell, target_cell)
            if witness is not None:
                return witness
        while self._traversal_position < len(representatives):
            for source_cell in self._expand_traversal_once():
                if _perfect_form_hash_key(source_cell) == target_key:
                    witness = backend.cell_transporter(source_cell, target_cell)
                    if witness is not None:
                        return witness
        return None

    @cached_property
    def _traversal(
        self,
    ) -> tuple[tuple[LorentzianPerfectCell, ...], tuple[LorentzianCellAdjacency, ...]]:
        self._ensure_traversal_started()
        representatives = self._traversal_representatives
        if representatives is None:
            raise ArithmeticError("a started Lorentzian traversal has no initial cell")
        while self._traversal_position < len(representatives):
            self._expand_traversal_once()
        return tuple(representatives), tuple(self._traversal_adjacencies)

    def quotient_cells(self) -> tuple[LorentzianPerfectCell, ...]:
        return self._traversal[0]

    def adjacencies(self) -> tuple[LorentzianCellAdjacency, ...]:
        return self._traversal[1]

    def component_preserving_group(self) -> RationalMatrixGroup:
        generators: list[LatticeIsometryMethods] = []
        for cell in self.quotient_cells():
            generators.extend(self.local_backend().cell_stabilizer(cell).generators())
        generators.extend(adjacency.transporter for adjacency in self.adjacencies() if adjacency.transporter is not None)
        return RationalMatrixGroup(self.lattice(), tuple(generators))

    def full_orthogonal_group(self) -> RationalMatrixGroup:
        component = self.component_preserving_group()
        negation = self.lattice().Aut()({label: -self.lattice().module_generator(label) for label in self.lattice().module_generating_set()})
        return RationalMatrixGroup(self.lattice(), component.generators() + (negation,))

    def records(self) -> tuple[TraversalRecord, ...]:
        cells = self.quotient_cells()
        positions = {id(cell): position for position, cell in enumerate(cells)}
        by_source: dict[int, list[TraversalAdjacencyRecord]] = {position: [] for position in range(len(cells))}
        for adjacency in self.adjacencies():
            source_position = positions[id(adjacency.source)]
            target_position = positions[id(adjacency.target)]
            transporter = adjacency.transporter
            if transporter is None:
                raise ArithmeticError("a completed quotient adjacency has no transporter")
            row_action = self.lattice().Aut()._row_action_matrix(transporter)
            incidence = [int(position in adjacency.facet) for position in range(len(adjacency.source.vector_configuration))]
            by_source[source_position].append(
                {
                    "iOrb": target_position,
                    "x": {
                        "eInc": incidence,
                        "eBigMat": [[int(entry) for entry in row] for row in row_action.rows()],
                    },
                }
            )
        records: list[TraversalRecord] = []
        for position, cell in enumerate(cells):
            configuration = cell.configuration()
            stabilizer = self.local_backend().cell_stabilizer(cell)
            permutations = [list(_configuration_permutation(configuration, generator)) for generator in stabilizer.generators()]
            records.append(
                {
                    "x": {
                        "EXT": [[int(entry) for entry in _coordinate_row(self.lattice(), vector, SageZZ)] for vector in cell.vector_configuration],
                        "GRP": permutations,
                    },
                    "ListAdj": by_source[position],
                }
            )
        return tuple(records)


def _perfect_form_hash_key(
    cell: LorentzianPerfectCell,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    r"""Return the exact pairing multisets used by upstream f_hash.

    ComputeInvariantPerfectForm hashes the multiplicities of diagonal and
    off-diagonal Gram pairings of the configuration. Keeping the exact
    multisets as the Python dictionary key avoids hash collisions while
    retaining precisely the same necessary isometry invariant.
    """
    gram = _gram_matrix(cell.lattice)
    rows = tuple(_coordinate_row(cell.lattice, item, SageZZ) for item in cell.vector_configuration)
    diagonal = tuple(sorted(int((row * gram * row.column())[0]) for row in rows))
    off_diagonal = tuple(sorted(int((rows[left] * gram * rows[right].column())[0]) for left in range(len(rows)) for right in range(left + 1, len(rows))))
    return diagonal, off_diagonal


def perfect_domain_traversal(gram_rows, option: PerfectMode = "total") -> tuple[TraversalRecord, ...]:
    r"""Return complete native perfect-domain traversal records for the Gram rows."""
    lattice = Lattices(ZZ)(gram_rows)
    return LorentzianPerfectComplex(lattice, option).records()


class NoLocalMarkTheoremError(RuntimeError):
    """The local perfect-cell mark theorem is unavailable for this norm."""


class MarkedCellOrbitAlgorithm:
    """Global isotropic-vertex orbits from local cell-stabilizer orbits."""

    def __init__(self, complex_: LorentzianPerfectComplex, norm=0) -> None:
        if norm != 0:
            raise NoLocalMarkTheoremError("the Lorentzian marked-cell source proves finite local marks only for norm zero")
        self._complex = complex_

    def complex(self) -> LorentzianPerfectComplex:
        return self._complex

    def local_marks(self, cell: LorentzianPerfectCell) -> tuple[Lattices.ElementMethods, ...]:
        zero = cell.lattice.base_ring().zero()
        return tuple(vector for vector in cell.vector_configuration if vector.q() == zero)

    def local_stabilizer_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
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
                    matches = [index for index, candidate in enumerate(marks) if candidate == moved]
                    if len(matches) != 1:
                        raise ArithmeticError("a cell stabilizer does not permute the isotropic local marks")
                    if matches[0] not in orbit_positions:
                        pending.append(matches[0])
            orbits.append(tuple(marks[index] for index in sorted(orbit_positions)))
        return tuple(orbits)

    def transport_marks(self, adjacency: LorentzianCellAdjacency) -> tuple[tuple[Lattices.ElementMethods, Lattices.ElementMethods], ...]:
        transporter = adjacency.transporter
        if transporter is None:
            raise ArithmeticError("a marked quotient adjacency needs a transporter")
        source_marks = tuple(
            adjacency.source.vector_configuration[position]
            for position in sorted(adjacency.facet)
            if adjacency.source.vector_configuration[position].q() == adjacency.source.lattice.base_ring().zero()
        )
        target_marks = self.local_marks(adjacency.target)
        transported = []
        for target_mark in target_marks:
            source_mark = transporter(target_mark)
            matches = [candidate for candidate in source_marks if candidate == source_mark]
            if len(matches) > 1:
                raise ArithmeticError("an isotropic facet mark does not transport to a unique mark of the target representative")
            if matches:
                transported.append((matches[0], target_mark))
        if len(transported) != len(source_marks):
            raise ArithmeticError("an isotropic facet mark has no unique mark in the target representative")
        return tuple(transported)

    def global_orbits(self) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        cells = self.complex().quotient_cells()
        local = tuple(self.local_stabilizer_orbits(cell) for cell in cells)
        nodes = [(cell_position, orbit_position) for cell_position, cell_orbits in enumerate(local) for orbit_position in range(len(cell_orbits))]
        graph = {node: set() for node in nodes}

        def local_orbit_position(cell_position, mark):
            matches = [orbit_position for orbit_position, orbit in enumerate(local[cell_position]) if any(candidate == mark for candidate in orbit)]
            if len(matches) != 1:
                raise ArithmeticError("a transported isotropic mark does not belong to a unique local stabilizer orbit")
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
        tuple(ring(_engine_element(base_ring, coordinates(label))) for label in lattice.module_generating_set()),
    )


def _normalizing_sign(lattice):
    positive, negative = lattice.signature_pair()
    match int(positive), int(negative):
        case 1, _:
            return lattice.base_ring().one()
        case _, 1:
            return -lattice.base_ring().one()
        case _:
            raise ValueError(f"{lattice} is not Lorentzian: its signature is {lattice.signature_pair()}")


@cache
def _normalized_lattice(lattice):
    sign = _normalizing_sign(lattice)
    return lattice if sign == lattice.base_ring().one() else lattice.twist(sign)


def _same_coordinates(source, target, element):
    if source is target:
        return element
    coordinates = element.to_vector()
    return _element_from_coordinates(
        target,
        tuple(coordinates(label) for label in source.module_generating_set()),
    )


def _ambient_element(lattice, coordinates):
    return _element_from_coordinates(lattice, tuple(SageZZ(entry) for entry in coordinates))


def _positive_direction(lattice) -> tuple[SageQQ, ...]:
    raw = _rational_positive_vector(lattice.gram_tensor())
    rationals = raw.base_ring()
    coordinates = tuple(SageQQ(_engine_element(rationals, entry)) for entry in raw)
    primitive, _scale = _primitive_integral_direction(coordinates)
    gram = _gram_matrix(lattice)
    rank = len(primitive)
    directions = []
    for first in range(rank):
        for sign in (-1, 1):
            row = [SageZZ.zero()] * rank
            row[first] = SageZZ(sign)
            directions.append(vector(SageZZ, row))
    for first in range(rank):
        for second in range(first + 1, rank):
            for first_sign in (-1, 1):
                for second_sign in (-1, 1):
                    row = [SageZZ.zero()] * rank
                    row[first] = SageZZ(first_sign)
                    row[second] = SageZZ(second_sign)
                    directions.append(vector(SageZZ, row))

    current = vector(SageZZ, primitive)
    while True:
        changes = 0
        for direction in directions:
            current_norm = (current * gram * current.column())[0]
            alpha = 0
            while True:
                candidate = current + (alpha + 1) * direction
                candidate_norm = (candidate * gram * candidate.column())[0]
                if 0 < candidate_norm < current_norm:
                    current_norm = candidate_norm
                    alpha += 1
                else:
                    break
            if alpha:
                current += alpha * direction
                changes += alpha
        if changes == 0:
            return tuple(SageQQ(entry) for entry in current)


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
    return tuple(value // content for value in integers), SageQQ(denominator) / SageQQ(content)


def _bezout_partner(lattice, timelike):
    ring = lattice.base_ring()
    labels = tuple(lattice.module_generating_set())
    pairings = tuple(lattice.b(timelike, lattice.module_generator(label)) for label in labels)
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
        raise ArithmeticError("Bezout coefficients do not realize the timelike divisibility")
    partner = lattice.linear_combination({label: coefficient for label, coefficient in zip(labels, coefficients, strict=True) if coefficient})
    if lattice.b(timelike, partner) != divisibility:
        raise ArithmeticError("the Bezout partner has the wrong pairing")
    return partner


class _PositiveVectorEnumerator:
    r"""Prepared exact positive-vector enumeration for one timelike direction."""

    def __init__(self, lattice, rational_direction, max_scal, mode: PerfectMode):
        self.lattice = lattice
        self.mode = mode
        direction_q = tuple(SageQQ(entry) for entry in rational_direction)
        primitive, scale = _primitive_integral_direction(direction_q)
        self.scaled_max = SageQQ(max_scal) * scale

        gram = _gram_matrix(lattice)
        self.gram = gram
        timelike_row = vector(SageZZ, primitive)
        square = SageZZ((timelike_row * gram * timelike_row.column())[0])
        if square <= 0:
            raise ValueError("positive-vector enumeration needs a timelike direction")

        pairing_row = timelike_row * gram
        kernel_rows = matrix(SageZZ, [pairing_row]).right_kernel_matrix()
        kernel_gram = -(kernel_rows * gram * kernel_rows.transpose())
        transform = kernel_gram.LLL_gram()
        self.kernel_rows = transform.transpose() * kernel_rows
        reduced_gram = -(self.kernel_rows * gram * self.kernel_rows.transpose())
        self.kernel_cvp_engine = _ExactCVPEngine._from_positive_engine_gram(
            lattice.base_ring(),
            reduced_gram,
        )

        gcd_value = SageZZ.zero()
        coefficients: list[SageZZ] = []
        for pairing in pairing_row:
            new_gcd, old_coefficient, new_coefficient = gcd_value.xgcd(SageZZ(pairing))
            coefficients = [old_coefficient * coefficient for coefficient in coefficients]
            coefficients.append(new_coefficient)
            gcd_value = new_gcd
        if gcd_value < 0:
            gcd_value = -gcd_value
            coefficients = [-coefficient for coefficient in coefficients]
        self.d = SageQQ(gcd_value)
        self.bezout_row = vector(SageZZ, coefficients)
        if self.bezout_row.dot_product(pairing_row) != gcd_value:
            raise ArithmeticError("Bezout coefficients do not realize timelike divisibility")

        alpha = self.d / SageQQ(square)
        translation = alpha * vector(SageQQ, timelike_row) - vector(SageQQ, self.bezout_row)
        self.base_target = tuple(self.kernel_rows.change_ring(SageQQ).transpose().solve_right(translation.column()).column(0))
        self.base_bound = (self.d * self.d) / SageQQ(square)

    def max_multiplier(self):
        if self.scaled_max <= 0:
            return None
        return int(self.scaled_max / self.d)

    def first_shell(self):
        maximum = self.max_multiplier()
        if maximum is None:
            return None
        shell = self.kernel_cvp_engine.first_close_vector_scale_coordinates(
            self.base_target,
            self.base_bound,
            maximum,
            exact_distance=self.mode == "isotropic",
        )
        if shell is None:
            return None
        multiplier, close_coordinates = shell
        return multiplier, self._ambient_vectors_from_shell(multiplier, close_coordinates)

    def vectors_at_multiplier(self, multiplier):
        target = tuple(SageQQ(multiplier) * entry for entry in self.base_target)
        bound = (SageQQ(multiplier) ** 2) * self.base_bound
        field = self.lattice.base_ring().fraction_field()
        owned_target = tuple(field(int(entry.numerator())) / field(int(entry.denominator())) for entry in target)
        owned_bound = field(int(bound.numerator())) / field(int(bound.denominator()))
        close_coordinates = self.kernel_cvp_engine.close_vector_coordinates(
            owned_target,
            owned_bound,
        )
        return self._ambient_vectors_from_shell(multiplier, close_coordinates)

    def _ambient_vectors_from_shell(self, multiplier, close_coordinates):
        bound = (SageQQ(multiplier) ** 2) * self.base_bound
        result = []
        for kernel_coordinates, kernel_square in close_coordinates:
            if self.mode == "isotropic" and kernel_square != bound:
                continue
            kernel_row = vector(
                SageZZ,
                kernel_coordinates,
            )
            ambient_row = SageZZ(multiplier) * self.bezout_row + self.kernel_rows.transpose() * kernel_row
            candidate_square = (ambient_row * self.gram * ambient_row.column())[0]
            if self.mode == "isotropic" and candidate_square != 0:
                raise ArithmeticError("an isotropic shell returned a nonisotropic vector")
            if self.mode == "total" and candidate_square < 0:
                raise ArithmeticError("a total shell returned a negative-norm vector")
            candidate = _ambient_element(self.lattice, tuple(ambient_row))
            result.append(candidate)
        return tuple(result)


def _find_positive_vectors(lattice, rational_direction, max_scal, mode: PerfectMode, *, only_shortest: bool):
    enumerator = _PositiveVectorEnumerator(lattice, rational_direction, max_scal, mode)
    if only_shortest and enumerator.scaled_max > 0:
        first = enumerator.first_shell()
        if first is None:
            return ()
        _multiplier, vectors = first
        return vectors

    result = []
    multiplier = 1
    while True:
        level = SageQQ(multiplier) * enumerator.d
        if enumerator.scaled_max > 0 and level > enumerator.scaled_max:
            break
        shell = enumerator.vectors_at_multiplier(multiplier)
        result.extend(shell)
        if only_shortest and shell:
            break
        multiplier += 1
    return tuple(result)


def _search_initial_vectors(lattice, direction, mode: PerfectMode):
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    direction_row = vector(SageQQ, direction)
    max_scal = (direction_row * gram * direction_row.column())[0]
    while True:
        vectors = _find_positive_vectors(lattice, direction, max_scal, mode, only_shortest=True)
        if vectors:
            return vectors
        max_scal *= 2


def _source_mid_value(lower, upper):
    r"""Port the upstream continued-fraction middle-value rule."""
    low = Fraction(int(lower.numerator()), int(lower.denominator()))
    upp = Fraction(int(upper.numerator()), int(upper.denominator()))
    midpoint = (low + upp) / 2
    target_low = (2 * low + upp) / 3
    target_upp = (low + 2 * upp) / 3

    terms = []
    work = midpoint
    while True:
        floor_value = work.numerator // work.denominator
        terms.append(floor_value)
        if work == floor_value:
            break
        work = 1 / (work - floor_value)

    for end in range(len(terms)):
        approximant = Fraction(terms[end], 1)
        for position in range(end - 1, -1, -1):
            approximant = terms[position] + 1 / approximant
        if target_low <= approximant <= target_upp:
            return SageQQ(approximant.numerator) / SageQQ(approximant.denominator)
    raise ArithmeticError("continued-fraction middle value was not found")


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
            raise ArithmeticError("the affine-normal constant gives a nonpositive flip bound")
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
        candidates = ((-b + root) / c, (-b - root) / c) if c != 0 else ((SageQQ(-a) / (2 * b),) if b != 0 else ())
        positive = [value for value in candidates if value > 0]
        isotropic_bound = min(positive) if positive else None
    if isotropic_bound is not None:
        bounds.append(isotropic_bound)
    if constant_bound is not None and isotropic_bound is not None and constant_bound == isotropic_bound:
        return None
    return min(bounds)


def _negative_direction_in_span(gram, spanning_rows):
    restricted = spanning_rows * gram * spanning_rows.transpose()
    diagonal, change = QuadraticForm(SageQQ, 2 * restricted).rational_diagonal_form(return_matrix=True)
    negative = [index for index in range(diagonal.matrix().nrows()) if diagonal.matrix()[index, index] < 0]
    if not negative:
        raise ArithmeticError("the perfect-cell nullspace has no negative direction")
    ambient = spanning_rows.transpose() * change.column(negative[0])
    return vector(SageQQ, ambient)


def _vector_rows(vectors):
    return {tuple(int(entry) for entry in _coordinate_row(item.parent(), item, SageZZ)) for item in vectors}


def _kernel_flipping(lattice, critical, base_normal, direction_normal, mode: PerfectMode):
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    upper = _upper_bound(gram, base_normal, direction_normal)
    if upper is None:
        raise ArithmeticError("the flip direction ends at a forbidden isotropic boundary")
    lower = SageQQ.zero()
    inverse = gram.inverse()
    critical_first = _coordinate_row(lattice, critical[0], SageQQ)
    critical_rows = _vector_rows(critical)
    total = ()
    while True:
        middle = _source_mid_value(lower, upper)
        normal = base_normal + middle * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = test_direction.column(0)
        square = (test_direction.transpose() * gram * test_direction)[0, 0]
        max_scal = (critical_first * gram * test_direction)[0]
        if square <= 0 or max_scal <= 0:
            upper = middle
            continue
        enumerator = _PositiveVectorEnumerator(lattice, tuple(test_direction_row), max_scal, mode)
        first = enumerator.first_shell()
        if first is None:
            raise ArithmeticError("a critical perfect-cell shell disappeared during flipping")
        _first_multiplier, total = first
        total_rows = _vector_rows(total)
        if total_rows == critical_rows:
            lower = middle
            continue
        if critical_rows.issubset(total_rows):
            return total, normal, test_direction_row, max_scal
        break
    while True:
        if not total:
            raise ArithmeticError("perfect-cell flipping reached an empty positive-vector shell")
        expanded = vector(SageQQ, [1, *_coordinate_row(lattice, total[0], SageQQ)])
        denominator = direction_normal.dot_product(expanded)
        if denominator == 0:
            raise ArithmeticError("the flip direction is parallel to a newly found wall")
        shift = -base_normal.dot_product(expanded) / denominator
        normal = base_normal + shift * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = test_direction.column(0)
        max_scal = (critical_first * gram * test_direction)[0]
        total = _find_positive_vectors(lattice, tuple(test_direction_row), max_scal, mode, only_shortest=True)
        if critical_rows.issubset(_vector_rows(total)):
            return total, normal, test_direction_row, max_scal
