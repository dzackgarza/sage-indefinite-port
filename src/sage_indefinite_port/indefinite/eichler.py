"""Eichler-cover primitives delegated to the preamble lattice owner."""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt

from dzack_research.preamble.all import QQ, ZZ, Modules
from dzack_research.preamble.categories.eichler_criterion import TwoUEichlerModel
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.categories.modules.module_morphisms.module_morphisms import ModuleEmbeddingMethods
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import FreeModuleElement
from sage.quadratic_forms.qfsolve import qfsolve
from sage.rings.integer import Integer as SageInteger
from sage.rings.rational_field import QQ as SageQQ

from sage_indefinite_port.backends.canonization import _rational_lattice_with_integral_structure
from sage_indefinite_port.groups.integral_structures import ArithmeticSubgroup, IntegralStructureAction, RationalMatrixGroup


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

    def lattice(self):
        return self.base.lattice()

    def subgroup(self):
        return self.refinement if self.refinement is not None else self.base.subgroup()

    def covering_representatives(self, norm, *, primitive: bool) -> OrbitCover:
        return self.base.covering_representatives(norm, primitive=primitive)

    def one_representative(self, norm, *, primitive: bool):
        return self.base.one_representative(norm, primitive=primitive)

    def refined_by(self, subgroup: ArithmeticSubgroup) -> OrbitCoverModel:
        return OrbitCoverModel(self.base, subgroup)

    def choose_splitting_vector(self, *, objective: str = "minimize_recursive_complexity"):
        return self.base.choose_splitting_vector(objective=objective)


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
    two_u_decomposition: TwoHyperbolicPlaneDecomposition
    integral_action: IntegralStructureAction

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
            return OrbitCover(tuple(family[label] for label in family.index_set()))
        if int(owned_norm) == 0:
            raise InfiniteLocusError("nonprimitive isotropic vectors form an infinite locus")
        lattice = self.lattice()
        representatives = []
        for divisor in square_divisors(owned_norm):
            primitive_norm = ring(int(owned_norm) // (divisor * divisor))
            family = self.model.covering_vector_representatives(primitive_norm)
            scalar = lattice.base_ring()(divisor)
            representatives.extend(lattice.scalar_multiple(scalar, family[label]) for label in family.index_set())
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
    first_v, first_w = find_hyperbolic_pair(lattice)
    first_plane = lattice.subobject_on((first_v, first_w)).saturation()
    first_complement = lattice.orthogonal_complement(first_plane)
    second_v, second_w = find_hyperbolic_pair(first_complement)
    second_v_ambient = first_complement.inclusion()(second_v)
    second_w_ambient = first_complement.inclusion()(second_w)
    two_u_span = lattice.subobject_on((first_v, first_w, second_v_ambient, second_w_ambient))
    complement = lattice.orthogonal_complement(two_u_span)
    model = complement.two_u_eichler_model()
    envelope = model.lattice()

    source_labels = tuple(lattice.module_generating_set())
    envelope_labels = tuple(envelope.module_generating_set())
    selected = (first_v, first_w, second_v_ambient, second_w_ambient) + tuple(complement.inclusion()(generator) for generator in complement.module_generators())
    from dzack_research.preamble.categories.rings.ring_foundation import _engine_element, _engine_ring

    source_gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageQQ)
    full_basis = matrix(
        _engine_ring(ring),
        [[_engine_element(ring, vector.to_vector()(label)) for label in source_labels] for vector in selected],
    ).change_ring(SageQQ)
    first_scale = SageQQ(_engine_element(ring, lattice.b(first_v, first_w)))
    second_scale = SageQQ(_engine_element(ring, lattice.b(second_v_ambient, second_w_ambient)))
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
    rational_envelope = envelope.base_change(ring.fraction_field_map())
    restriction = Modules(QQ).restriction_of_scalars(ZZ.Mor(QQ)(lambda element: QQ(element)))
    integral_structure_space = restriction(rational_envelope)
    inclusion = _rational_lattice_with_integral_structure(
        lattice,
        rational_envelope,
        integral_structure_space,
        embedding_rows,
    )
    rational_group = RationalMatrixGroup(rational_envelope, ())
    action = IntegralStructureAction(rational_group, inclusion)
    return EichlerEnvelope(
        lattice,
        envelope,
        inclusion,
        inclusion_to_envelope,
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
    basis = tuple(lattice.module_generators())
    for left in basis:
        if not left or not left.is_isotropic():
            continue
        for right in basis:
            if not right or not right.is_isotropic():
                continue
            pairing = lattice.b(left, right)
            if pairing == lattice.base_ring().zero():
                continue
            if pairing < lattice.base_ring().zero():
                right = -right
            return left, right
    if int(lattice.module_rank()) == 2 and lattice.signature_pair() == (lattice.base_ring().one(), lattice.base_ring().one()):
        candidates = basis + (basis[0] + basis[1], basis[0] - basis[1], -basis[0] + basis[1], -basis[0] - basis[1])
        isotropic = tuple(vector for vector in candidates if vector and vector.is_isotropic())
        for left in isotropic:
            for right in isotropic:
                pairing = lattice.b(left, right)
                if pairing != lattice.base_ring().zero():
                    if pairing < lattice.base_ring().zero():
                        right = -right
                    return left, right
    gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageQQ)
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
    labels = tuple(lattice.module_generating_set())
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
