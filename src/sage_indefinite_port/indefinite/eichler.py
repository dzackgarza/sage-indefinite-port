"""Eichler-cover primitives delegated to the preamble lattice owner."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import gcd, isqrt

from dzack_research.preamble.all import QQ, ZZ, Modules
from dzack_research.preamble.categories.eichler_criterion import TwoUEichlerModel
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.categories.modules.module_morphisms.module_morphisms import ModuleEmbeddingMethods
from dzack_research.preamble.categories.rings.ring_foundation import _engine_element, _engine_ring
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.modules.free_module_element import FreeModuleElement
from sage.quadratic_forms.qfsolve import qfsolve
from sage.rings.integer import Integer as SageInteger
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ

from sage_indefinite_port.backends.canonization import _rational_lattice_with_integral_structure
from sage_indefinite_port.groups.integral_structures import ArithmeticSubgroup, GeneratedSubgroup, IntegralStructureAction, RationalMatrixGroup


class InfiniteLocusError(ValueError):
    """Raised when a requested orbit cover is not finite."""


@dataclass(frozen=True)
class OrbitCover:
    """A finite covering family, not necessarily a full-orbit decomposition."""

    representatives: tuple[Lattices.ElementMethods, ...]

    def __iter__(self):
        return iter(self.representatives)

    def __len__(self) -> int:
        return len(self.representatives)


@dataclass(frozen=True)
class OrbitCoverModel:
    """Immutable finite-cover model with an optional arithmetic refinement."""

    base: EichlerOrbitCover
    refinement: ArithmeticSubgroup | None = None
    envelope: EichlerEnvelope | None = None

    @classmethod
    def from_lattice(cls, lattice: Lattices.ParentMethods) -> OrbitCoverModel:
        match lattice.is_decomposable():
            case True:
                factors = tuple(lattice.biproduct_factors())
                hyperbolic_gram = Lattices(lattice.base_ring())("U").gram_tensor()
                match (
                    len(factors) >= 2,
                    factors[0].gram_tensor() == hyperbolic_gram if factors else False,
                    factors[1].gram_tensor() == hyperbolic_gram if len(factors) >= 2 else False,
                ):
                    case (True, True, True):
                        return cls(EichlerOrbitCover(lattice.two_u_eichler_model_from_represented_biproduct()))
                    case _:
                        pass
            case False:
                pass
        envelope = build_eichler_envelope(lattice)
        base = EichlerOrbitCover(envelope.two_u_decomposition.lattice.two_u_eichler_model_from_represented_biproduct())
        return cls(base, envelope=envelope)

    def lattice(self):
        return self.base.lattice() if self.envelope is None else self.envelope.lattice

    def subgroup(self):
        if self.refinement is not None:
            return self.refinement
        if self.envelope is None:
            return self.base.subgroup()
        embedding = self.envelope._embedding_matrix
        embedding_inverse = self.envelope._embedding_inverse
        lattice_automorphisms = self.envelope.lattice.Aut()
        generators = []
        for column_matrix in self.envelope.integral_action.finite_representation()._lattice_stabilizer_column_matrices:
            rational_action = column_matrix.transpose()
            pulled_action = embedding * rational_action * embedding_inverse
            if any(entry.denominator() != 1 for entry in pulled_action.list()):
                raise ArithmeticError("a selected-lattice stabilizer does not preserve the Eichler embedding")
            if pulled_action * embedding != embedding * rational_action:
                raise ArithmeticError("the conjugated Eichler stabilizer does not commute with the lattice embedding")
            generators.append(lattice_automorphisms._isometry_from_column_matrix(pulled_action.transpose()))
        generators = tuple(generators)
        return GeneratedSubgroup(RationalMatrixGroup(self.envelope.lattice, generators), generators)

    def covering_representatives(self, norm, *, primitive: bool) -> OrbitCover:
        if self.envelope is None:
            return self.base.covering_representatives(norm, primitive=primitive)
        source = self.envelope.lattice
        owned_norm = norm if getattr(norm, "parent", lambda: None)() is source.base_ring() else source.base_ring()(int(norm))
        scaled_norm = owned_norm * self.envelope.similarity_scale()
        if scaled_norm == source.base_ring().zero():
            if not primitive:
                raise InfiniteLocusError("nonprimitive isotropic vectors form an infinite locus")
            base_cover = self.base.covering_representatives(scaled_norm, primitive=True)
        else:
            base_cover = self.base.covering_representatives(scaled_norm, primitive=False)

        rational_envelope = self.envelope.integral_action.rational_group().rational_lattice()
        fraction_map = source.base_ring().fraction_field_map()
        envelope_labels = tuple(self.envelope.envelope.module_generating_set())
        rational_labels = tuple(rational_envelope.module_generating_set())
        finite_representation = self.envelope.integral_action.finite_representation()
        representatives = []
        for representative in base_cover:
            coordinates = representative.to_vector()
            rational_representative = rational_envelope.linear_combination(
                {
                    rational_label: fraction_map(coordinates(envelope_label))
                    for envelope_label, rational_label in zip(envelope_labels, rational_labels, strict=True)
                    if coordinates(envelope_label)
                }
            )
            for moved_representative in finite_representation.orbit_images(rational_representative):
                candidate = self.envelope.pullback_element(moved_representative)
                if candidate is None or candidate.q() != owned_norm:
                    continue
                if primitive and not candidate.is_primitive():
                    continue
                if all(candidate != known for known in representatives):
                    representatives.append(candidate)
        return OrbitCover(tuple(representatives))

    def one_representative(self, norm, *, primitive: bool):
        cover = self.covering_representatives(norm, primitive=primitive)
        if not cover.representatives:
            raise ValueError(f"no covering representative exists for norm {norm}")
        return cover.representatives[0]

    def refined_by(self, subgroup: ArithmeticSubgroup) -> OrbitCoverModel:
        return OrbitCoverModel(self.base, subgroup, self.envelope)

    def choose_splitting_vector(self, *, objective: str = "minimize_recursive_complexity"):
        match self.envelope:
            case None:
                return self.base.choose_splitting_vector(objective=objective)
            case _:
                match objective:
                    case "minimize_recursive_complexity":
                        pass
                    case _:
                        raise ValueError("the implemented splitting-vector objective is 'minimize_recursive_complexity'")
                lattice = self.envelope.lattice
                basis = tuple(lattice.module_generators())
                candidates = list(basis)
                for left_position, left in enumerate(basis):
                    for right in basis[left_position + 1 :]:
                        candidates.append(left + right)
                        candidates.append(left - right)
                positive = tuple(vector for vector in candidates if vector.q() > lattice.base_ring().zero())
                match positive:
                    case ():
                        raise ValueError("the Eichler-envelope lattice has no positive vector in its framing span")
                    case _:
                        pass
                labels = tuple(lattice.module_generating_set())
                return min(
                    positive,
                    key=lambda vector: (
                        abs(int(vector.q())),
                        tuple(int(vector.to_vector()(label)) for label in labels),
                    ),
                )


@dataclass(frozen=True)
class TwoHyperbolicPlaneDecomposition:
    """The explicit ``U + U + K`` decomposition retained by an Eichler model."""

    lattice: Lattices.ParentMethods
    first_hyperbolic_plane: Lattices.ParentMethods
    second_hyperbolic_plane: Lattices.ParentMethods
    complement: Lattices.ParentMethods
    sum_isometry: tuple[ModuleEmbeddingMethods, ModuleEmbeddingMethods, ModuleEmbeddingMethods]

    @classmethod
    def from_model(cls, model: TwoUEichlerModel) -> TwoHyperbolicPlaneDecomposition:
        lattice = model.lattice()
        return cls(
            lattice,
            model.first_hyperbolic_plane(),
            model.second_hyperbolic_plane(),
            model.orthogonal_complement(),
            (lattice.injection(0), lattice.injection(1), lattice.injection(2)),
        )


@dataclass(frozen=True)
class EichlerEnvelope:
    """An ambient Eichler lattice together with the integral lattice embedded in it."""

    lattice: Lattices.ParentMethods
    envelope: Lattices.ParentMethods
    inclusion: ModuleEmbeddingMethods
    lattice_to_envelope: ModuleEmbeddingMethods
    embedding_rows: tuple[tuple[int, ...], ...]
    two_u_decomposition: TwoHyperbolicPlaneDecomposition
    integral_action: IntegralStructureAction

    @cached_property
    def _embedding_matrix(self) -> Matrix_rational_dense:
        return matrix(SageQQ, self.embedding_rows)

    @cached_property
    def _embedding_inverse(self) -> Matrix_rational_dense:
        return self._embedding_matrix.inverse()

    def similarity_scale(self):
        r"""Return the scale of ``lattice_to_envelope`` from its defining pairings."""
        source = self.lattice
        target = self.envelope
        generators = tuple(source.module_generators())
        scale = None
        for left in generators:
            for right in generators:
                source_pairing = source.b(left, right)
                if source_pairing == source.base_ring().zero():
                    continue
                target_pairing = target.b(
                    self.lattice_to_envelope(left),
                    self.lattice_to_envelope(right),
                )
                scale = target_pairing / source_pairing
                break
            if scale is not None:
                break
        if scale is None:
            raise ArithmeticError("a nondegenerate Eichler envelope has no nonzero pairing in its frame")
        if any(target.b(self.lattice_to_envelope(left), self.lattice_to_envelope(right)) != scale * source.b(left, right) for left in generators for right in generators):
            raise ArithmeticError("the Eichler embedding does not have one similarity scale")
        return scale

    def pullback_element(self, element):
        r"""Pull a rational-envelope element back to the selected lattice when it lies there."""
        rational_envelope = self.integral_action.rational_group().rational_lattice()
        rational_labels = tuple(rational_envelope.module_generating_set())
        coordinate_vector = rational_envelope(element).to_vector()
        envelope_coordinates = matrix(
            SageQQ,
            [
                [
                    SageQQ(
                        _engine_element(
                            rational_envelope.base_ring(),
                            coordinate_vector(label),
                        )
                    )
                    for label in rational_labels
                ]
            ],
        )
        source_coordinates = envelope_coordinates * self._embedding_inverse
        if any(entry.denominator() != 1 for entry in source_coordinates.list()):
            return None
        source_labels = tuple(self.lattice.module_generating_set())
        row = source_coordinates.row(0)
        return self.lattice.linear_combination({source_label: self.lattice.base_ring()(int(entry)) for source_label, entry in zip(source_labels, row, strict=True) if entry})

    def pullback_isometry(self, isometry: LatticeIsometryMethods) -> LatticeIsometryMethods:
        r"""Restrict a rational-envelope isometry preserving the selected lattice."""
        rational_envelope = self.integral_action.rational_group().rational_lattice()
        rational_automorphisms = rational_envelope.Aut()
        embedding = matrix(SageQQ, self.embedding_rows)
        rational_action = rational_automorphisms._row_action_matrix(isometry)
        pulled_action = embedding * rational_action * embedding.inverse()
        if any(entry.denominator() != 1 for entry in pulled_action.list()):
            raise ArithmeticError("a selected-lattice stabilizer does not preserve the Eichler embedding")
        witness = self.lattice.Aut()._isometry_from_column_matrix(pulled_action.transpose())
        if pulled_action * embedding != embedding * rational_action:
            raise ArithmeticError("the conjugated Eichler stabilizer does not commute with the lattice embedding")
        return witness

    def stabilizer_of_original_lattice(self, group=None) -> ArithmeticSubgroup:
        if group is not None and group is not self.integral_action.rational_group():
            raise ValueError("the requested group is not the envelope's represented rational group")
        return self.integral_action.finite_representation().lattice_stabilizer()

    def right_cosets(self):
        return self.integral_action.right_cosets()


@dataclass(frozen=True)
class EichlerOrbitCover:
    """Immutable covering model for the preamble's represented ``2U + K`` lattice."""

    model: TwoUEichlerModel

    def lattice(self):
        return self.model.lattice()

    def subgroup(self):
        return self.model.approximate_generating_family()

    def covering_representatives(self, norm, *, primitive: bool) -> OrbitCover:
        ring = self.lattice().base_ring()
        owned_norm = norm if getattr(norm, "parent", lambda: None)() is ring else ring(int(norm))
        if primitive:
            family = self.model.covering_vector_representatives(owned_norm)
            return OrbitCover(tuple(family))
        if int(owned_norm) == 0:
            raise InfiniteLocusError("nonprimitive isotropic vectors form an infinite locus")
        lattice = self.lattice()
        representatives = []
        for divisor in square_divisors(owned_norm):
            primitive_norm = ring(int(owned_norm) // (divisor * divisor))
            family = self.model.covering_vector_representatives(primitive_norm)
            scalar = lattice.base_ring()(divisor)
            representatives.extend(lattice.scalar_multiple(scalar, representative) for representative in family)
        return OrbitCover(tuple(representatives))

    def one_representative(self, norm, *, primitive: bool):
        cover = self.covering_representatives(norm, primitive=primitive)
        if not cover.representatives:
            raise ValueError(f"no covering representative exists for norm {norm}")
        return cover.representatives[0]

    def choose_splitting_vector(self, *, objective: str = "minimize_recursive_complexity"):
        if objective != "minimize_recursive_complexity":
            raise ValueError("the implemented splitting-vector objective is 'minimize_recursive_complexity'")
        lattice = self.lattice()
        basis = tuple(lattice.module_generators())
        candidates = list(basis)
        for left_position, left in enumerate(basis):
            for right in basis[left_position + 1 :]:
                candidates.append(left + right)
                candidates.append(left - right)
        positive = tuple(vector for vector in candidates if vector.q() > lattice.base_ring().zero())
        if not positive:
            raise ValueError("the Eichler model lattice has no positive vector in its framing span")
        labels = tuple(lattice.module_generating_set())
        return min(
            positive,
            key=lambda vector: (
                abs(int(vector.q())),
                tuple(int(vector.to_vector()(label)) for label in labels),
            ),
        )


def build_eichler_envelope(lattice: Lattices.ParentMethods) -> EichlerEnvelope:
    """Construct a literal-2U Eichler envelope containing ``lattice``."""
    ring = lattice.base_ring()
    source_gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
    working_lattice = Lattices(ring)([[ring(int(source_gram[row, column])) for column in range(source_gram.ncols())] for row in range(source_gram.nrows())])
    working_labels = tuple(working_lattice.module_generating_set())
    first_v, first_w = find_hyperbolic_pair(working_lattice)
    first_pair_rows = matrix(
        SageZZ,
        [
            [
                SageZZ(int(_engine_element(ring, vector.to_vector()(label))))
                for label in working_labels
            ]
            for vector in (first_v, first_w)
        ],
    )
    first_complement_rows = (first_pair_rows * source_gram).right_kernel_matrix()
    first_complement_gram = first_complement_rows * source_gram * first_complement_rows.transpose()
    first_complement = Lattices(ring)(
        [
            [
                ring(int(first_complement_gram[row, column]))
                for column in range(first_complement_gram.ncols())
            ]
            for row in range(first_complement_gram.nrows())
        ]
    )
    second_v, second_w = find_hyperbolic_pair(first_complement)
    first_complement_labels = tuple(first_complement.module_generating_set())
    second_pair_local_rows = matrix(
        SageZZ,
        [
            [
                SageZZ(int(_engine_element(ring, vector.to_vector()(label))))
                for label in first_complement_labels
            ]
            for vector in (second_v, second_w)
        ],
    )
    second_pair_rows = second_pair_local_rows * first_complement_rows
    two_u_rows = first_pair_rows.stack(second_pair_rows)
    complement_rows = (two_u_rows * source_gram).right_kernel_matrix()
    complement_gram = complement_rows * source_gram * complement_rows.transpose()
    match complement_rows.nrows():
        case 0:
            complement = Lattices(ring)(ring.free_module(0))
        case _:
            complement = Lattices(ring)(
                [
                    [
                        ring(int(complement_gram[row, column]))
                        for column in range(complement_gram.ncols())
                    ]
                    for row in range(complement_gram.nrows())
                ]
            )
    model = complement.two_u_eichler_model()
    envelope = model.lattice()

    envelope_labels = tuple(envelope.module_generating_set())
    source_gram = source_gram.change_ring(SageQQ)
    full_basis = two_u_rows.stack(complement_rows).change_ring(SageQQ)
    first_scale = (first_pair_rows.row(0) * source_gram * first_pair_rows.row(1).column())[0]
    second_scale = (second_pair_rows.row(0) * source_gram * second_pair_rows.row(1).column())[0]
    if first_scale <= 0 or second_scale <= 0:
        raise ArithmeticError("the selected hyperbolic pairs must have positive pairing")
    normalization = matrix.identity(SageQQ, int(lattice.module_rank()))
    normalization[0, 0] = SageQQ.one() / first_scale
    normalization[2, 2] = SageQQ.one() / second_scale
    envelope_basis = normalization * full_basis
    envelope_gram = envelope_basis * source_gram * envelope_basis.transpose()
    if any(entry.denominator() != 1 for entry in envelope_gram.list()):
        raise ArithmeticError("the Eichler reduction did not produce an integral envelope form")
    represented_envelope_gram = _engine_component_matrix(envelope.gram_tensor()).change_ring(SageQQ)
    if envelope_gram != represented_envelope_gram:
        raise ArithmeticError("the Eichler reduction does not match the represented 2U envelope")

    pre_embedding = envelope_basis.inverse()
    denominator = 1
    for entry in pre_embedding.list():
        entry_denominator = int(entry.denominator())
        denominator = denominator * entry_denominator // gcd(denominator, entry_denominator)
    embedding_matrix = SageQQ(denominator) * pre_embedding
    if any(entry.denominator() != 1 for entry in embedding_matrix.list()):
        raise ArithmeticError("clearing denominators did not produce an integral Eichler embedding")
    embedding_rows = tuple(tuple(int(entry) for entry in embedding_matrix.row(position)) for position in range(embedding_matrix.nrows()))

    def image(label):
        position = int(lattice.module_generating_set().ranking_map()(label))
        row = embedding_rows[position]
        return envelope.linear_combination({envelope_label: ring(coefficient) for envelope_label, coefficient in zip(envelope_labels, row, strict=True) if coefficient})

    inclusion_to_envelope = lattice.module_category().Mor(lattice, envelope)(image)
    rational_envelope, rational_matrix_family = model._approximate_generator_column_matrices_after_base_change(ring.fraction_field_map())
    rational_generator_matrices = tuple(rational_matrix_family[label] for label in rational_matrix_family.index_set())
    if not rational_generator_matrices:
        raise ArithmeticError("the Eichler approximate family has no generators")
    restriction = Modules(QQ).restriction_of_scalars(ZZ.Mor(QQ)(lambda element: QQ(element)))
    integral_structure_space = restriction(rational_envelope)
    inclusion = _rational_lattice_with_integral_structure(
        lattice,
        rational_envelope,
        integral_structure_space,
        embedding_rows,
    )
    rational_group = RationalMatrixGroup._from_column_matrices(
        rational_envelope,
        rational_generator_matrices,
    )
    action = IntegralStructureAction(rational_group, inclusion)
    return EichlerEnvelope(
        lattice,
        envelope,
        inclusion,
        inclusion_to_envelope,
        embedding_rows,
        TwoHyperbolicPlaneDecomposition.from_model(model),
        action,
    )


def eichler_transvection(
    isotropic: Lattices.ElementMethods,
    orthogonal: Lattices.ElementMethods,
) -> LatticeIsometryMethods:
    r"""Return the exact Eichler transvection ``E_(f,x)`` in the ambient lattice."""
    lattice = isotropic.parent()
    if orthogonal.parent() is not lattice:
        raise ValueError("an Eichler transvection needs two vectors in one lattice")
    if not lattice.is_even():
        raise ValueError("an Eichler transvection in this port requires an even lattice")
    if not isotropic.is_isotropic():
        raise ValueError("the first Eichler-transvection vector must be isotropic")
    if lattice.b(isotropic, orthogonal) != lattice.base_ring().zero():
        raise ValueError("the second Eichler-transvection vector must lie in f^perp")
    return lattice.eichler_transvection(isotropic, orthogonal)


def square_divisors(integer) -> tuple[int, ...]:
    r"""Return positive ``c`` such that ``c^2`` divides the nonzero integer ``integer``."""
    value = abs(int(integer))
    if value == 0:
        raise ValueError("zero has infinitely many square divisors")
    return tuple(divisor for divisor in range(1, isqrt(value) + 1) if value % (divisor * divisor) == 0)


def find_hyperbolic_pair(
    lattice: Lattices.ParentMethods,
    *,
    primitive: bool = True,
    method: str = "auto",
) -> tuple[Lattices.ElementMethods, Lattices.ElementMethods]:
    r"""Return verified integral isotropic ``v,w`` with ``b(v,w)>0``.

    PARI's exact ``qfsolve`` supplies one rational isotropic direction.  Clearing
    denominators and dividing the coordinate gcd gives an integral primitive
    vector ``v``.  For any integral ``h`` with ``d=b(v,h)>0``, the vector

    ``w = 2 d h - q(h) v``

    is integral isotropic and satisfies ``b(v,w)=2d^2``.
    """
    if method != "auto":
        raise ValueError("the implemented hyperbolic-pair method is 'auto'")
    labels = tuple(lattice.module_generating_set())
    rank = len(labels)
    gram = _engine_component_matrix(lattice.gram_tensor())
    candidate_coordinates = []
    for position in range(rank):
        coordinates = tuple(1 if index == position else 0 for index in range(rank))
        candidate_coordinates.append(coordinates)
    for left_position in range(rank):
        for right_position in range(left_position + 1, rank):
            candidate_coordinates.append(tuple(1 if index in (left_position, right_position) else 0 for index in range(rank)))
            candidate_coordinates.append(tuple(1 if index == left_position else -1 if index == right_position else 0 for index in range(rank)))
    isotropic_coordinates = []
    for coordinates in candidate_coordinates:
        norm = sum(coordinates[row] * gram[row, column] * coordinates[column] for row in range(rank) for column in range(rank))
        if norm == 0:
            isotropic_coordinates.append(coordinates)
    small_pairs = []
    for left in isotropic_coordinates:
        for right in isotropic_coordinates:
            pairing = sum(left[row] * gram[row, column] * right[column] for row in range(rank) for column in range(rank))
            if pairing == 0:
                continue
            oriented_right = right
            if pairing < 0:
                oriented_right = tuple(-coordinate for coordinate in right)
                pairing = -pairing
            small_pairs.append(
                (
                    int(pairing),
                    left,
                    oriented_right,
                )
            )
    match small_pairs:
        case []:
            pass
        case _:
            _pairing, left_coordinates, right_coordinates = min(small_pairs)
            ring = lattice.base_ring()
            left = lattice.linear_combination({label: ring(coordinate) for label, coordinate in zip(labels, left_coordinates, strict=True) if coordinate})
            right = lattice.linear_combination({label: ring(coordinate) for label, coordinate in zip(labels, right_coordinates, strict=True) if coordinate})
            return left, right
    gram = gram.change_ring(SageQQ)
    solution = qfsolve(gram)
    if isinstance(solution, SageInteger):
        raise ValueError("the lattice is anisotropic over QQ")
    if not isinstance(solution, FreeModuleElement):
        raise ValueError("qfsolve returned a degenerate isotropic subspace instead of one vector")
    denominators = [entry.denominator() for entry in solution]
    denominator = 1
    for entry_denominator in denominators:
        denominator = denominator * int(entry_denominator) // gcd(denominator, int(entry_denominator))
    coordinates = [int(denominator * entry) for entry in solution]
    content = 0
    for coordinate in coordinates:
        content = gcd(content, abs(coordinate))
    if content == 0:
        raise ArithmeticError("qfsolve returned the zero vector")
    coordinates = [coordinate // content for coordinate in coordinates]
    ring = lattice.base_ring()
    v = lattice.linear_combination({label: ring(coordinate) for label, coordinate in zip(labels, coordinates, strict=True) if coordinate})
    if not v.is_isotropic():
        raise ArithmeticError("the saturated qfsolve witness is not isotropic")
    if primitive and not v.is_primitive():
        raise ArithmeticError("the saturated qfsolve witness is not primitive")

    h = None
    d = ring.zero()
    for label in labels:
        candidate = lattice.module_generator(label)
        pairing = lattice.b(v, candidate)
        if pairing != ring.zero():
            h = candidate
            d = pairing
            break
    if h is None:
        raise ArithmeticError("a nonzero isotropic vector pairs trivially with every lattice generator")
    if d < ring.zero():
        h = lattice.scalar_multiple(-ring.one(), h)
        d = -d
    two_d = ring(2) * d
    w = lattice.scalar_multiple(two_d, h) - lattice.scalar_multiple(h.q(), v)
    if not w.is_isotropic() or lattice.b(v, w) <= ring.zero():
        raise ArithmeticError("the constructed partner does not form a verified hyperbolic pair")
    return v, w


__all__ = [
    "EichlerEnvelope",
    "build_eichler_envelope",
    "EichlerOrbitCover",
    "InfiniteLocusError",
    "OrbitCover",
    "OrbitCoverModel",
    "TwoHyperbolicPlaneDecomposition",
    "eichler_transvection",
    "find_hyperbolic_pair",
    "square_divisors",
]
