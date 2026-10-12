"""Exact Lorentzian perfect-cell local backend."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property
from random import choice, randrange, shuffle
from typing import Literal

from dzack_research.preamble.all import Groups, Lattices
from dzack_research.preamble.categories.group.g_sets import FiniteGSets
from dzack_research.preamble.categories.group.groups import GroupMorphism
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones
from dzack_research.preamble.categories.reduction_complexes import LorentzianPerfectDomainTraversal, PerfectDomainOrbitAdjacency
from dzack_research.preamble.categories.sets.finite_ordered_sets import finite_ordered_set
from dzack_research.preamble.categories.sets.indexed_families import finite_indexed_family
from dzack_research.preamble.categories.vector_configurations import VectorConfigurations
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.matrix.special import identity_matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ
from sage.sets.real_set import RealSet
from sage.structure.element import RingElement
from sage.structure.parent import Parent

from sage_indefinite_port.invariants import AttackProfile

type PerfectMode = Literal["total", "isotropic"]
type FacetIncidence = frozenset[int]

# The admissible squares in the normalized signature-(1,n) form are actual
# subsets of R, not a source-level string selector.
NONNEGATIVE_SQUARES = RealSet.unbounded_above_closed(0)
ZERO_SQUARES = RealSet.point(0)


def lorentzian_orthogonal_group_generators(lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
    r"""Generators of O(L) for L of signature (1, n) or (n, 1): O^Omega(L) from the perfect-domain traversal, and -1."""
    return LorentzianPerfectComplex(lattice, NONNEGATIVE_SQUARES).full_orthogonal_group_generators()


def lorentzian_isometry(source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
    r"""An isometry between two lattices of signature (1, n) or (n, 1), found by matching perfect cells, or None."""
    return LorentzianPerfectComplex(source, NONNEGATIVE_SQUARES).isometry_to(LorentzianPerfectComplex(target, NONNEGATIVE_SQUARES))


def _gram_matrix(lattice: Lattices.ParentMethods) -> Matrix_integer_dense:
    return _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)


@dataclass(frozen=True)
class LorentzianPerfectCell:
    lattice: Lattices.ParentMethods
    vector_configuration: tuple[Lattices.ElementMethods, ...]
    allowed_squares: RealSet

    def __post_init__(self) -> None:
        if self.allowed_squares not in (NONNEGATIVE_SQUARES, ZERO_SQUARES):
            raise ValueError(f"unsupported normalized-norm locus {self.allowed_squares!r}")
        positive, negative = self.lattice.signature_pair()
        if int(positive) != 1 and int(negative) != 1:
            raise ValueError(f"a Lorentzian perfect cell needs signature (1,n) or (n,1), but {self.lattice} has signature {self.lattice.signature_pair()}")
        gram = _gram_matrix(self.lattice)
        rows = tuple(_coordinate_row(self.lattice, point, SageZZ) for point in self.vector_configuration)
        squares = tuple(int(row * gram * row) for row in rows)
        normalized_sign = 1 if int(positive) == 1 else -1
        if any(normalized_sign * square not in self.allowed_squares for square in squares):
            raise ValueError(f"a perfect cell may contain only vectors with normalized square in {self.allowed_squares}")

    def configuration(self) -> Parent:
        r"""The marked lattice (L, Z^S -> L) of the cell's vectors, an object of Ar(Lattices)."""
        return self.lattice.vector_configuration(list(self.vector_configuration))


@dataclass(frozen=True)
class LorentzianCellAdjacency:
    source: LorentzianPerfectCell
    facet: FacetIncidence
    target: LorentzianPerfectCell
    transporter: LatticeIsometryMethods | None


class LorentzianPerfectLocalBackend:
    def initial_cell(self, lattice: Lattices.ParentMethods, allowed_squares: RealSet) -> LorentzianPerfectCell:
        normalized = _lorentzian_view(lattice)
        direction = _positive_direction(normalized)
        normalized_vectors = _search_initial_vectors(normalized, direction, allowed_squares)
        gram = _gram_matrix(normalized).change_ring(SageQQ)
        central = vector(SageQQ, direction)
        first = _coordinate_row(normalized, normalized_vectors[0], SageQQ)
        scalar = (central * gram * first.column())[0]
        base_normal = vector(SageQQ, [-scalar, *(gram * central.column()).column(0)])
        while matrix(SageQQ, [_coordinate_row(normalized, item, SageQQ) for item in normalized_vectors]).rank() < int(lattice.module_rank()):
            expanded = matrix(SageQQ, [[1, *_coordinate_row(normalized, item, SageQQ)] for item in normalized_vectors])
            nullspace = expanded.right_kernel_matrix()
            directions = nullspace.matrix_from_columns(range(1, nullspace.ncols())) * gram.inverse()
            outside_direction = _negative_direction_in_span(normalized, directions, normalized_vectors, base_normal)
            direction_covector = gram * outside_direction.column()
            critical_row = _coordinate_row(normalized, normalized_vectors[0], SageQQ)
            direction_scalar = (critical_row * direction_covector)[0]
            normal_direction = vector(SageQQ, [-direction_scalar, *direction_covector.column(0)])
            normalized_vectors, base_normal, _test_direction, _max_scal = _kernel_flipping(normalized, normalized_vectors, base_normal, normal_direction, allowed_squares)
        vectors = tuple(_transfer(lattice, item) for item in normalized_vectors)
        return LorentzianPerfectCell(lattice, vectors, allowed_squares)

    @cached_method
    def cell_stabilizer(self, cell: LorentzianPerfectCell) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, S)) for the cell's vector configuration S."""
        return tuple(
            VectorConfigurations(cell.lattice.base_ring()).spanning_automorphism_generators(cell.lattice, cell.vector_configuration, (0,) * len(cell.vector_configuration))
        )

    @cached_method
    def facet_lift(self, cell: LorentzianPerfectCell) -> GroupMorphism:
        r"""Evaluate words in the given integral cell-stabilizer generators.

        The source's finite facet-orbit calculation requires the action of
        these generators, not construction of the infinite arithmetic matrix
        subgroup of O(L) in its general group backend.
        """
        generators = self.cell_stabilizer(cell)
        free = Groups.Free(index_set=finite_ordered_set(range(len(generators))))
        return free.Mor(cell.lattice.O())({position: generator for position, generator in enumerate(generators)})

    @cached_method
    def facet_action(self, cell: LorentzianPerfectCell) -> Parent:
        r"""The finite action of Aut((L, S)) on the facets spanned by S."""
        vectors = cell.vector_configuration
        cone = RationalPolyhedralCones(cell.lattice).from_rays(vectors)
        facets = tuple(frozenset(position for position, item in enumerate(vectors) if cone.evaluate_covector(covector, item) == 0) for covector in cone.facet_covectors())

        lift = self.facet_lift(cell)

        def act(word, facet: FacetIncidence) -> FacetIncidence:
            isometry = lift(word)
            return frozenset(vectors.index(isometry(vectors[position])) for position in facet)

        return FiniteGSets(lift.domain())(facets, act)

    def facet_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[FacetIncidence, ...], ...]:
        return tuple(tuple(orbit.points()) for orbit in self.facet_action(cell).orbits())

    def flip_across(self, cell: LorentzianPerfectCell, facet: FacetIncidence) -> LorentzianPerfectCell | None:
        normalized = _lorentzian_view(cell.lattice)
        vectors = cell.vector_configuration
        if not facet or any(position < 0 or position >= len(vectors) for position in facet):
            raise ValueError("a perfect-cell facet must be a nonempty incidence subset")
        nonincident = next((position for position in range(len(vectors)) if position not in facet), None)
        if nonincident is None:
            raise ValueError("a facet must omit at least one configuration vector")
        # In a Lorentzian plane, either null ray is an actual boundary of the
        # chosen positive cone, not a side shared with a second perfect cell.
        # With null basis (u,v), the cone is {a u + b v : a,b >= 0}; crossing
        # the face b=0 yields ab < 0 and leaves the admissible cone. Do not
        # try the source's interior-facet flip across a boundary face.
        if int(cell.lattice.module_rank()) == 2 and all(cell.lattice.q(vectors[position]) == 0 for position in facet):
            return None
        selected = matrix(SageQQ, [_coordinate_row(normalized, vectors[position], SageQQ) for position in sorted(facet)])
        facet_nullspace = selected.right_kernel_matrix()
        if facet_nullspace.nrows() != 1:
            raise ArithmeticError("a perfect-cell facet must have one normal direction")
        direction = facet_nullspace.row(0)
        outside = _coordinate_row(normalized, vectors[nonincident], SageQQ)
        if direction.dot_product(outside) < 0:
            direction = -direction
        normal_direction = vector(SageQQ, [0, *direction])
        expanded = matrix(SageQQ, [[1, *_coordinate_row(normalized, item, SageQQ)] for item in vectors])
        base_nullspace = expanded.right_kernel_matrix()
        if base_nullspace.nrows() != 1:
            raise ArithmeticError("a perfect cell must have a unique affine supporting functional")
        base_normal = base_nullspace.row(0)
        if vector(SageQQ, base_normal[1:]).dot_product(outside) <= 0:
            base_normal = -base_normal
        critical = tuple(vectors[position] for position in sorted(facet))
        flipped, _normal, _direction, _max_scal = _kernel_flipping(normalized, critical, base_normal, normal_direction, cell.allowed_squares)
        returned = tuple(_transfer(cell.lattice, item) for item in flipped)
        return LorentzianPerfectCell(cell.lattice, returned, cell.allowed_squares)

    def cell_transporter(self, source: LorentzianPerfectCell, target: LorentzianPerfectCell) -> LatticeIsometryMethods | None:
        r"""An element of Isom((L, S), (L, T)) for the cells' configurations, or None."""
        if source.allowed_squares != target.allowed_squares:
            raise ValueError("an isometry of perfect cells must preserve their prescribed normalized-norm locus")
        return VectorConfigurations(source.lattice.base_ring()).spanning_isometry(
            source.lattice,
            source.vector_configuration,
            (0,) * len(source.vector_configuration),
            target.lattice,
            target.vector_configuration,
            (0,) * len(target.vector_configuration),
        )


class LorentzianPerfectComplex:
    r"""Complete quotient traversal of one Lorentzian perfect-domain complex.

    A quotient adjacency stores a selected target representative and the exact
    isometry carrying it onto the actual neighbor obtained by flipping the
    retained source facet.
    """

    def __init__(self, lattice: Lattices.ParentMethods, allowed_squares: RealSet) -> None:
        if allowed_squares not in (NONNEGATIVE_SQUARES, ZERO_SQUARES):
            raise ValueError(f"unsupported normalized-norm locus {allowed_squares!r}")
        self._lattice = lattice
        self._allowed_squares = allowed_squares
        self._local_backend = LorentzianPerfectLocalBackend()
        self._traversal_representatives: list[LorentzianPerfectCell] | None = None
        self._traversal_representative_buckets: dict[tuple[tuple[int, ...], tuple[int, ...]], list[LorentzianPerfectCell]] = {}
        self._traversal_adjacencies: list[LorentzianCellAdjacency] = []
        self._traversal_position = 0

    def lattice(self) -> Lattices.ParentMethods:
        return self._lattice

    def allowed_squares(self) -> RealSet:
        return self._allowed_squares

    def local_backend(self) -> LorentzianPerfectLocalBackend:
        return self._local_backend

    def _ensure_traversal_started(self) -> None:
        if self._traversal_representatives is not None:
            return
        initial = self.local_backend().initial_cell(self.lattice(), self.allowed_squares())
        self._traversal_representatives = [initial]
        self._traversal_representative_buckets = {_perfect_form_hash_key(initial): [initial]}

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
            if neighbor is None:
                # A boundary face has no adjacent perfect cell or side pairing.
                continue
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
            action = backend.facet_action(source)
            for incident in orbit:
                symmetry = backend.facet_lift(source)(action.transporter_witness(facet, incident))
                self._traversal_adjacencies.append(LorentzianCellAdjacency(source, incident, target, symmetry * transporter))
        self._traversal_position += 1
        return tuple(new_representatives)

    def isometry_to(self, target: LorentzianPerfectComplex) -> LatticeIsometryMethods | None:
        if self.allowed_squares() != target.allowed_squares():
            raise ValueError("an isometry of perfect-domain complexes must preserve the prescribed normalized-norm locus")
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
    def _traversal(self) -> tuple[tuple[LorentzianPerfectCell, ...], tuple[LorentzianCellAdjacency, ...]]:
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

    def component_preserving_generators(self) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of O^Omega(L): the cells' automorphism groups and the side pairings."""
        generators: list[LatticeIsometryMethods] = []
        for cell in self.quotient_cells():
            generators.extend(self.local_backend().cell_stabilizer(cell))
        generators.extend(adjacency.transporter for adjacency in self.adjacencies() if adjacency.transporter is not None)
        return tuple(generators)

    def full_orthogonal_group_generators(self) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of O(L) = O^Omega(L) x {±1}."""
        negation = self.lattice().Aut()({label: -self.lattice().module_generator(label) for label in self.lattice().module_generating_set()})
        return self.component_preserving_generators() + (negation,)


def _perfect_form_hash_key(cell: LorentzianPerfectCell) -> tuple[tuple[int, ...], tuple[int, ...]]:
    r"""Return the exact pairing multisets used by upstream f_hash.

    ComputeInvariantPerfectForm hashes the multiplicities of diagonal and
    off-diagonal Gram pairings of the configuration. Keeping the exact
    multisets as the Python dictionary key avoids hash collisions while
    retaining precisely the same necessary isometry invariant.
    """
    rows = matrix(SageZZ, [_coordinate_row(cell.lattice, point, SageZZ) for point in cell.vector_configuration])
    pairings = rows * _gram_matrix(cell.lattice) * rows.transpose()
    diagonal = tuple(sorted(int(pairings[index, index]) for index in range(rows.nrows())))
    off_diagonal = tuple(sorted(int(pairings[left, right]) for left in range(rows.nrows()) for right in range(left + 1, rows.nrows())))
    return diagonal, off_diagonal


def perfect_domain_traversal(lattice: Lattices.ParentMethods, mode: PerfectMode) -> LorentzianPerfectDomainTraversal:
    r"""The perfect-domain decomposition of a Lorentzian lattice in the given mode.

    Returns the preamble's ``LorentzianPerfectDomainTraversal`` of ``lattice``: the quotient
    cells, their automorphism groups and the side pairings, from which O^Omega(L) is generated.
    Each cell is the cone on its vectors; an adjacency carries the source cell, the
    representative of the neighbour's orbit, the neighbour itself as the image of that
    representative, the facet they share, and the isometry carrying the representative onto it.
    """
    # Decode the frozen #33 historical selector at its sole entry point;
    # the actual algorithm receives a mathematical norm locus.
    if mode == "total":
        allowed_squares = NONNEGATIVE_SQUARES
    elif mode == "isotropic":
        allowed_squares = ZERO_SQUARES
    else:
        raise ValueError(f"unsupported perfect-domain selection {mode!r}")
    complex_ = LorentzianPerfectComplex(lattice, allowed_squares)
    backend = complex_.local_backend()
    cones = RationalPolyhedralCones(lattice)
    cell_cones = {id(cell): cones.from_rays(cell.vector_configuration) for cell in complex_.quotient_cells()}
    adjacencies = []
    for adjacency in complex_.adjacencies():
        source = adjacency.source.vector_configuration
        target_to_neighbour = adjacency.transporter
        adjacencies.append(
            PerfectDomainOrbitAdjacency(
                cell_cones[id(adjacency.source)],
                cell_cones[id(adjacency.target)],
                cones.from_rays(tuple(target_to_neighbour(vector_) for vector_ in adjacency.target.vector_configuration)),
                cones.from_rays(tuple(source[position] for position in sorted(adjacency.facet))),
                target_to_neighbour,
            )
        )
    cells = complex_.quotient_cells()
    cones_by_position = tuple(cell_cones[id(cell)] for cell in cells)
    stabilizers = finite_indexed_family(finite_ordered_set(cones_by_position), lambda cone: backend.cell_stabilizer(cells[cones_by_position.index(cone)]))
    return LorentzianPerfectDomainTraversal(lattice, tuple(cell_cones.values()), adjacencies, stabilizers)


class NoLocalMarkTheoremError(RuntimeError):
    """The local perfect-cell mark theorem does not cover this norm: vectors outside the closed positive cone."""


class MarkedCellOrbitAlgorithm:
    r"""O(L)-orbits of primitive vectors of square n in the closed positive cone, from the quotient cells.

    The marks of a cell are the primitive lattice vectors of square n in its cone: its isotropic
    rays for n = 0, and for timelike n the finite set the preamble enumerates. It is finite: a
    lattice vector of the cone is a combination of the rays with coefficients in (1/D)Z_{>=0} for
    a D fixed by the cell, distinct rays of the closed positive cone pair positively, and every
    term of q is nonnegative, so q = n bounds each coefficient. Marks equivalent under a cell's
    automorphism group are one local orbit, a mark on a shared facet is identified with its
    preimage under the side pairing, and the cells around any face are connected through facet
    adjacencies, so the connected components are the orbits. An element of O(L) carrying one
    vector of the cone to another preserves the cone's sheet, so these are the O(L)-orbits.
    """

    def __init__(self, complex_: LorentzianPerfectComplex, norm: int) -> None:
        sign = AttackProfile.from_lattice(complex_.lattice()).sign
        if sign * norm < 0:
            raise NoLocalMarkTheoremError(f"vectors of square {norm} in {complex_.lattice()} are spacelike: they lie outside the closed positive cone that the cells cover")
        self._complex = complex_
        self._norm = norm

    def complex(self) -> LorentzianPerfectComplex:
        return self._complex

    def local_marks(self, cell: LorentzianPerfectCell) -> tuple[Lattices.ElementMethods, ...]:
        if self._norm == 0:
            zero = cell.lattice.base_ring().zero()
            return tuple(vector for vector in cell.vector_configuration if vector.q() == zero)
        cone = RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
        return tuple(vector for vector in cone.vectors_of_square(self._norm) if vector.content() == 1)

    def local_stabilizer_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        marks = self.local_marks(cell)
        lift = self.complex().local_backend().facet_lift(cell)
        positions = tuple(range(len(marks)))

        def act(word, position: int) -> int:
            isometry = lift(word)
            return marks.index(isometry(marks[position]))

        return tuple(tuple(marks[position] for position in sorted(orbit.points())) for orbit in FiniteGSets(lift.domain())(positions, act).orbits())

    def transport_marks(self, adjacency: LorentzianCellAdjacency) -> tuple[tuple[Lattices.ElementMethods, Lattices.ElementMethods], ...]:
        transporter = adjacency.transporter
        if transporter is None:
            raise ArithmeticError("a marked quotient adjacency needs a transporter")
        source = adjacency.source.vector_configuration
        facet = RationalPolyhedralCones(adjacency.source.lattice).from_rays(tuple(source[position] for position in sorted(adjacency.facet)))
        source_marks = tuple(mark for mark in self.local_marks(adjacency.source) if facet.contains(mark))
        target_marks = self.local_marks(adjacency.target)
        transported = []
        for target_mark in target_marks:
            source_mark = transporter(target_mark)
            matches = [candidate for candidate in source_marks if candidate == source_mark]
            if len(matches) > 1:
                raise ArithmeticError("a facet mark does not transport to a unique mark of the target representative")
            if matches:
                transported.append((matches[0], target_mark))
        if len(transported) != len(source_marks):
            raise ArithmeticError("a facet mark has no unique mark in the target representative")
        return tuple(transported)

    def global_orbits(self) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        cells = self.complex().quotient_cells()
        local = tuple(self.local_stabilizer_orbits(cell) for cell in cells)
        nodes = [(cell_position, orbit_position) for cell_position, cell_orbits in enumerate(local) for orbit_position in range(len(cell_orbits))]
        graph: dict[tuple[int, int], set[tuple[int, int]]] = {node: set() for node in nodes}

        def local_orbit_position(cell_position: int, mark: Lattices.ElementMethods) -> int:
            matches = [orbit_position for orbit_position, orbit in enumerate(local[cell_position]) if any(candidate == mark for candidate in orbit)]
            if len(matches) != 1:
                raise ArithmeticError("a transported isotropic mark does not belong to a unique local stabilizer orbit")
            return matches[0]

        positions = {id(cell): position for position, cell in enumerate(cells)}
        for adjacency in self.complex().adjacencies():
            source_position = positions[id(adjacency.source)]
            target_position = positions[id(adjacency.target)]
            for source_mark, target_mark in self.transport_marks(adjacency):
                source_node = (source_position, local_orbit_position(source_position, source_mark))
                target_node = (target_position, local_orbit_position(target_position, target_mark))
                graph[source_node].add(target_node)
                graph[target_node].add(source_node)

        unseen = set(nodes)
        global_orbits: list[tuple[Lattices.ElementMethods, ...]] = []
        while unseen:
            start = min(unseen)
            pending = [start]
            component: list[tuple[int, int]] = []
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


def _coordinate_row[S: RingElement](lattice: Lattices.ParentMethods, element: Lattices.ElementMethods, ring: Parent[S]) -> FreeModuleElement[S]:
    coordinates = element.to_vector()
    return vector(ring, tuple(ring(int(coordinates(label))) for label in lattice.module_generating_set()))


def _rational(entry: RingElement) -> Rational:
    return SageQQ(int(entry.numerator())) / SageQQ(int(entry.denominator()))


def _lorentzian_view(lattice: Lattices.ParentMethods) -> Lattices.ParentMethods:
    r"""L or L(-1), whichever has signature (1, n)."""
    profile = AttackProfile.from_lattice(lattice)
    if profile.positive_index != 1:
        raise ValueError(f"{lattice} is not Lorentzian: its signature is {lattice.signature_pair()}")
    return profile.signed_view


def _transfer(target: Lattices.ParentMethods, element: Lattices.ElementMethods) -> Lattices.ElementMethods:
    r"""The element with the same coordinates in a lattice on the same module (L and its twist)."""
    coordinates = element.to_vector()
    return target.linear_combination({label: coordinates(label) for label in target.module_generating_set()})


def _ambient_element(lattice: Lattices.ParentMethods, coordinates: tuple[Integer, ...]) -> Lattices.ElementMethods:
    ring = lattice.base_ring()
    return lattice.linear_combination({label: ring(int(entry)) for label, entry in zip(lattice.module_generating_set(), coordinates, strict=True) if entry})


def _owned_rational(lattice: Lattices.ParentMethods, value: Rational) -> RingElement:
    field = lattice.base_ring().fraction_field()
    return field(int(value.numerator())) / field(int(value.denominator()))


def _positive_direction(lattice: Lattices.ParentMethods) -> tuple[Rational, ...]:
    r"""A positive vector of L, primitive, then greedily shortened along e_i and e_i +- e_j."""
    positive = lattice.positive_vector()
    content = int(positive.content())
    primitive = tuple(int(entry) // content for entry in _coordinate_row(lattice, positive, SageZZ))
    gram = _gram_matrix(lattice)
    rank = len(primitive)
    directions: list[FreeModuleElement[Integer]] = []
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


def _primitive_integral_direction(lattice: Lattices.ParentMethods, coordinates: tuple[Rational, ...]) -> tuple[Lattices.ElementMethods, Rational]:
    r"""The primitive vector t of L on the ray of a rational direction u, and the scale s with u = t / s."""
    denominator = SageZZ.one()
    for entry in coordinates:
        denominator = denominator.lcm(SageZZ(SageQQ(entry).denominator()))
    scaled = _ambient_element(lattice, tuple(SageZZ(SageQQ(entry) * denominator) for entry in coordinates))
    content = int(scaled.content())
    primitive = _ambient_element(lattice, tuple(int(entry) // content for entry in _coordinate_row(lattice, scaled, SageZZ)))
    return primitive, SageQQ(denominator) / SageQQ(content)


class _PositiveVectorEnumerator:
    r"""Vectors v of L with q(v) >= 0 (total) or q(v) = 0 (isotropic) and 0 < b(v, t) <= max, by multiplier.

    With t primitive timelike, d = div(t) and h a Bezout partner (b(t, h) = d), a vector with
    b(v, t) = m d is v = m h + x for x in K = t^perp, and q(v) = m^2 d^2 / q(t) + q(x - m w) for
    w = (d / q(t)) t - h in K tensor Q. K is negative definite, so the admissible x are the close
    vectors of m w in K within m^2 d^2 / q(t) (on the sphere, in isotropic mode).
    """

    def __init__(self, lattice: Lattices.ParentMethods, rational_direction: tuple[Rational, ...], max_scal: Rational, allowed_squares: RealSet) -> None:
        self.lattice = lattice
        self.gram = _gram_matrix(lattice)
        self.allowed_squares = allowed_squares
        timelike, scale = _primitive_integral_direction(lattice, tuple(SageQQ(entry) for entry in rational_direction))
        self.scaled_max = SageQQ(max_scal) * scale
        square = SageQQ(int(timelike.q()))
        if square <= 0:
            raise ValueError("positive-vector enumeration needs a timelike direction")
        self.d = SageQQ(int(timelike.div()))
        self.partner = timelike.bezout_partner()
        self.inclusion = timelike.orthogonal_complement().inclusion()
        self.kernel = self.inclusion.domain()
        fraction_map = lattice.base_ring().fraction_field_map()
        rational = lattice.base_change(fraction_map)
        alpha = self.d / square
        timelike_row = _coordinate_row(lattice, timelike, SageZZ)
        partner_row = _coordinate_row(lattice, self.partner, SageZZ)
        translation = rational.linear_combination(
            {label: _owned_rational(lattice, alpha * timelike_row[position] - partner_row[position]) for position, label in enumerate(lattice.module_generating_set())}
        )
        self.base_target = self.inclusion.base_change(fraction_map).lift(translation)
        self.base_bound = (self.d * self.d) / square

    def max_multiplier(self) -> int | None:
        if self.scaled_max <= 0:
            return None
        return int(self.scaled_max / self.d)

    def first_shell(self) -> tuple[int, tuple[Lattices.ElementMethods, ...]] | None:
        maximum = self.max_multiplier()
        if maximum is None:
            return None
        for multiplier in range(1, maximum + 1):
            vectors = self.vectors_at_multiplier(multiplier)
            if vectors:
                return multiplier, vectors
        return None

    def vectors_at_multiplier(self, multiplier: int) -> tuple[Lattices.ElementMethods, ...]:
        target = self.base_target.parent().scalar_multiple(_owned_rational(self.lattice, SageQQ(multiplier)), self.base_target)
        bound = _owned_rational(self.lattice, -(SageQQ(multiplier) ** 2) * self.base_bound)
        # Lorentzian.g:369-373 delegates this finite exact distance fibre to
        # the general closest-vector computation; use its semantic API.
        distances = self.kernel.close_vectors(target, bound)
        return self._ambient_vectors(multiplier, tuple(distances.index_set()))

    def _ambient_vectors(self, multiplier: int, close_vectors: tuple[Lattices.ElementMethods, ...]) -> tuple[Lattices.ElementMethods, ...]:
        shift = self.lattice.scalar_multiple(self.lattice.base_ring()(multiplier), self.partner)
        result: list[Lattices.ElementMethods] = []
        for close in close_vectors:
            candidate = shift + self.inclusion(close)
            row = _coordinate_row(self.lattice, candidate, SageZZ)
            candidate_square = int(row * self.gram * row)
            if candidate_square < 0:
                raise ArithmeticError("a close-vector shell returned a negative-norm vector")
            if candidate_square not in self.allowed_squares:
                continue
            result.append(candidate)
        return tuple(result)


def _first_positive_shell(
    lattice: Lattices.ParentMethods, rational_direction: tuple[Rational, ...], max_scal: Rational, allowed_squares: RealSet
) -> tuple[Lattices.ElementMethods, ...]:
    r"""The vectors of the first nonempty shell with 0 < b(v, t) <= max_scal, or ()."""
    enumerator = _PositiveVectorEnumerator(lattice, rational_direction, max_scal, allowed_squares)
    if enumerator.scaled_max <= 0:
        raise ArithmeticError(f"the positive-vector shells below {max_scal} along {rational_direction} are empty by definition: the bound must be positive")
    first = enumerator.first_shell()
    return () if first is None else first[1]


def _search_initial_vectors(lattice: Lattices.ParentMethods, direction: tuple[Rational, ...], allowed_squares: RealSet) -> tuple[Lattices.ElementMethods, ...]:
    r"""The first shell along the direction, doubling the bound until it is nonempty; it terminates because L has vectors of the mode's kind."""
    if allowed_squares == ZERO_SQUARES and not lattice.is_isotropic():
        raise ValueError(f"{lattice} is anisotropic: it has no isotropic perfect cells")
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    direction_row = vector(SageQQ, direction)
    max_scal = (direction_row * gram * direction_row.column())[0]
    while True:
        vectors = _first_positive_shell(lattice, direction, max_scal, allowed_squares)
        if vectors:
            return vectors
        max_scal *= 2


def _source_mid_value(lower: Rational, upper: Rational) -> Rational:
    r"""Port the upstream continued-fraction middle-value rule."""
    low = Fraction(int(lower.numerator()), int(lower.denominator()))
    upp = Fraction(int(upper.numerator()), int(upper.denominator()))
    midpoint = (low + upp) / 2
    target_low = (2 * low + upp) / 3
    target_upp = (low + 2 * upp) / 3

    terms: list[int] = []
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


def _upper_bound(gram: Matrix_rational_dense, base_normal: FreeModuleElement[Rational], direction_normal: FreeModuleElement[Rational]) -> Rational | None:
    inverse = gram.inverse()
    base_constant = SageQQ(base_normal[0])
    direction_constant = SageQQ(direction_normal[0])
    base_covector = vector(SageQQ, base_normal[1:])
    direction_covector = vector(SageQQ, direction_normal[1:])
    bounds: list[Rational] = []
    constant_bound: Rational | None = None
    if direction_constant > 0:
        constant_bound = -base_constant / direction_constant
        if constant_bound <= 0:
            raise ArithmeticError("the affine-normal constant gives a nonpositive flip bound")
        bounds.append(constant_bound)
    base_vector = inverse * base_covector.column()
    direction_vector = inverse * direction_covector.column()
    base_square = (base_vector.transpose() * gram * base_vector)[0, 0]
    cross = (base_vector.transpose() * gram * direction_vector)[0, 0]
    direction_square = (direction_vector.transpose() * gram * direction_vector)[0, 0]
    if base_square <= 0:
        raise ArithmeticError("LORENTZ_Kernel_Flipping requires a positive-norm base normal")
    if direction_square >= 0:
        raise ArithmeticError("LORENTZ_Kernel_Flipping requires a negative-norm direction normal")
    # Lorentzian.g:778-790. The required negative-norm direction guarantees
    # eventual negativity, so the reference's doubling search terminates.
    shift = SageQQ.one()
    while base_square + 2 * cross * shift + direction_square * shift * shift >= 0:
        shift *= 2
    bounds.append(shift)
    basis = matrix(SageQQ, [base_vector.column(0), direction_vector.column(0)])
    restricted = basis * gram * basis.transpose()
    a, b, c = restricted[0, 0], restricted[0, 1], restricted[1, 1]
    discriminant = b * b - a * c
    isotropic_bound: Rational | None = None
    if discriminant > 0 and SageQQ(discriminant).is_square():
        root = SageQQ(discriminant).sqrt()
        candidates = ((-b + root) / c, (-b - root) / c) if c != 0 else ((SageQQ(-a) / (2 * b),) if b != 0 else ())
        positive = [value for value in candidates if value > 0]
        isotropic_bound = min(positive) if positive else None
    if isotropic_bound is not None:
        bounds.append(isotropic_bound)
    if constant_bound is not None and isotropic_bound is not None and constant_bound == isotropic_bound:
        return None
    upper = min(bounds)
    # Lorentzian.g:823-831: a returned upper bound must lie on or outside
    # the positive cone. A smaller affine-constant bound is not silently
    # accepted as an admissible cone exit.
    if base_square + 2 * cross * upper + direction_square * upper * upper > 0:
        raise ArithmeticError("the computed Lorentzian upper bound remains inside the positive cone")
    return upper


def _negative_direction_in_span(
    lattice: Lattices.ParentMethods,
    spanning_rows: Matrix_rational_dense,
    critical: tuple[Lattices.ElementMethods, ...],
    base_normal: FreeModuleElement[Rational],
) -> FreeModuleElement[Rational]:
    r"""The negative direction accepted by Lorentzian.g's GetOneOutsideRay.

    The source perturbs the basis until GetUpperBound does not report its
    inadmissible isotropic endpoint. A single arbitrary negative vector is
    insufficient for that algorithm.
    """
    rational = lattice.base_change(lattice.base_ring().fraction_field_map())
    labels = tuple(lattice.module_generating_set())
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    perturbation = identity_matrix(SageZZ, spanning_rows.nrows())
    while True:
        rows = perturbation * spanning_rows
        span = rational.vector_configuration(
            [rational.linear_combination({label: _owned_rational(lattice, entry) for label, entry in zip(labels, row, strict=True) if entry}) for row in rows.rows()]
        )
        negative = span.inclusion()(span.vector_of_sign(-1)).to_vector()
        direction = vector(SageQQ, tuple(_rational(negative(label)) for label in labels))
        covector = gram * direction.column()
        scalars = {(_coordinate_row(lattice, point, SageQQ) * covector)[0] for point in critical}
        if len(scalars) != 1:
            raise ArithmeticError("the GetOneOutsideRay direction must pair equally with every critical vector")
        (scalar,) = scalars
        normal = vector(SageQQ, [-scalar, *covector.column(0)])
        if _upper_bound(gram, base_normal, normal) is not None:
            return direction
        # IndefiniteFormsFundamental.g:11-30, GetRandomMatrixPerturbation.
        # These changes of basis preserve the rational span of the rows.
        dimension = spanning_rows.nrows()
        modification = identity_matrix(SageZZ, dimension)
        if choice((1, 2)) == 1:
            permutation = list(range(dimension))
            shuffle(permutation)
            modification = matrix(SageZZ, dimension, dimension)
            for index, image in enumerate(permutation):
                modification[index, image] = choice((-1, 1))
        else:
            index = randrange(dimension)
            image = randrange(dimension)
            modification[index, image] = choice((-1, 1))
        perturbation *= modification


def _vector_rows(vectors: Iterable[Lattices.ElementMethods]) -> set[tuple[int, ...]]:
    return {tuple(int(entry) for entry in _coordinate_row(item.parent(), item, SageZZ)) for item in vectors}


def _kernel_flipping(
    lattice: Lattices.ParentMethods,
    critical: tuple[Lattices.ElementMethods, ...],
    base_normal: FreeModuleElement[Rational],
    direction_normal: FreeModuleElement[Rational],
    allowed_squares: RealSet,
) -> tuple[tuple[Lattices.ElementMethods, ...], FreeModuleElement[Rational], FreeModuleElement[Rational], Rational]:
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    # Lorentzian.g:853-860. Both affine normals annihilate the critical
    # configuration and are independent. Those are hypotheses of the
    # one-parameter wall crossing, not optional diagnostics.
    if matrix(SageQQ, [base_normal, direction_normal]).rank() != 2:
        raise ArithmeticError("Lorentzian flip normals must be independent")
    for point in critical:
        augmented = vector(SageQQ, [1, *_coordinate_row(lattice, point, SageQQ)])
        if base_normal.dot_product(augmented) != 0 or direction_normal.dot_product(augmented) != 0:
            raise ArithmeticError("both Lorentzian flip normals must vanish on the critical vectors")
        if allowed_squares == ZERO_SQUARES and point.q() != lattice.base_ring().zero():
            raise ArithmeticError("an isotropic flip requires isotropic critical vectors")
    upper = _upper_bound(gram, base_normal, direction_normal)
    if upper is None:
        raise ArithmeticError("the flip direction ends at a forbidden isotropic boundary")
    lower = SageQQ.zero()
    inverse = gram.inverse()
    critical_first = _coordinate_row(lattice, critical[0], SageQQ)
    critical_rows = _vector_rows(critical)
    total: tuple[Lattices.ElementMethods, ...] = ()
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
        enumerator = _PositiveVectorEnumerator(lattice, tuple(test_direction_row), max_scal, allowed_squares)
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
        # Lorentzian.g:949-965. The next supporting hyperplane is defined by
        # the new vector from the previous shell; the new shell must contain
        # that vector as well as every critical vector. Preserving only the
        # critical set is insufficient to certify the source's second exit.
        new_vector = total[0]
        expanded = vector(SageQQ, [1, *_coordinate_row(lattice, new_vector, SageQQ)])
        denominator = direction_normal.dot_product(expanded)
        if denominator == 0:
            raise ArithmeticError("the flip direction is parallel to a newly found wall")
        shift = -base_normal.dot_product(expanded) / denominator
        normal = base_normal + shift * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = test_direction.column(0)
        max_scal = (critical_first * gram * test_direction)[0]
        total = _first_positive_shell(lattice, tuple(test_direction_row), max_scal, allowed_squares)
        required = critical_rows | _vector_rows((new_vector,))
        if required.issubset(_vector_rows(total)):
            return total, normal, test_direction_row, max_scal
