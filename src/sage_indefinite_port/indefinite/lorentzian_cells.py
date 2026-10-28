"""Exact Lorentzian perfect-cell local backend."""

from __future__ import annotations

from collections import Counter, deque
from array import array
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from fractions import Fraction
from functools import cache, cached_property
from itertools import combinations, islice, product
from math import gcd, isfinite, isqrt
from sys import float_info
from random import choice, randrange, shuffle
from typing import Literal, cast, overload

import numpy as np
from cypari2.handle_error import PariError
from dzack_research.preamble.all import FramedFreeModules, Lattices, Sets
from dzack_research.preamble.categories.definite_lattices import LatticeReduction, _ExactCVPEngine
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.rings.ring_foundation import _owned_engine_element
from dzack_research.preamble.categories.modules.pure.modules import _engine_matrix
from dzack_research.preamble.categories.polyhedral_cones import (
    RationalPolyhedralCones,
    _integral_covector,
    _rational_polyhedral_cone,
)
from dzack_research.preamble.categories.reduction_complexes import (
    LorentzianPerfectDomainTraversal,
    PerfectDomainOrbitAdjacency,
)
from dzack_research.preamble.categories.sets.finite_ordered_sets import (
    finite_ordered_set,
)
from dzack_research.preamble.categories.sets.indexed_families import (
    finite_indexed_family,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.libs.gap.libgap import libgap
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import diagonal_matrix, matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.matrix.matrix_real_double_dense import Matrix_real_double_dense
from sage.matrix.special import identity_matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.complex_double import ComplexDoubleElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ
from sage.rings.real_double import RDF
from sage.sets.real_set import RealSet
from sage.structure.element import RingElement
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.indefinite_reduction import (
    _compute_indefinite_reduction,
    _indefinite_reduction,
)
from sage_indefinite_port.invariants import AttackProfile

type PerfectMode = Literal["total", "isotropic"]
type FacetIncidence = frozenset[int]
type SchreierWord = tuple[SchreierWord, int] | None
type FacetMaskKey = int | tuple[bool, tuple[int, ...] | int]


def _sparse_facet_complement_limit(ray_count: int) -> int:
    """Estimate the tuple/bitset memory crossover for this ray universe."""
    # On 64-bit CPython, a tuple slot plus a non-cached Python int costs about
    # 36 bytes, while the bitset costs one bit per ray. The conservative
    # 28-byte offset accounts for the tuple and integer object headers.
    return max(0, (ray_count - 224) // 288)


def _facet_complement_key(mask: int, ray_count: int | None = None) -> FacetMaskKey:
    """Encode a facet sparsely, or bit-pack the smaller side when both are dense."""
    if ray_count is None:
        ray_count = max(1, mask.bit_length())
    complement_count = mask.bit_count()
    incident_count = ray_count - complement_count
    limit = _sparse_facet_complement_limit(ray_count)
    if complement_count <= limit:
        return (True, tuple(_bit_positions(mask)))
    if incident_count <= limit:
        incident_mask = ((1 << ray_count) - 1) ^ mask
        return (False, tuple(_bit_positions(incident_mask)))
    if incident_count < complement_count:
        incident_mask = ((1 << ray_count) - 1) ^ mask
        return False, incident_mask
    return mask


def _facet_incidence(key: FacetMaskKey, ray_count: int) -> FacetIncidence:
    """Recover the exact incident-ray set from its compact facet key."""
    if isinstance(key, int):
        return frozenset(_bit_positions(((1 << ray_count) - 1) ^ key))
    complement, positions = key
    if isinstance(positions, int):
        mask = ((1 << ray_count) - 1) ^ positions if complement else positions
        return frozenset(_bit_positions(mask))
    if complement:
        # Avoid constructing a full frozenset(range(ray_count)) and then a
        # second full frozenset for the difference. The sparse exclusion set
        # and final incidence set are the only retained containers.
        excluded = set(positions)
        return frozenset(position for position in range(ray_count) if position not in excluded)
    return frozenset(positions)


def _permuted_facet_key(key: FacetMaskKey, permutation: array | tuple[int, ...]) -> FacetMaskKey:
    """Apply a ray permutation while preserving the sparse/dense key form."""
    if isinstance(key, int):
        return _facet_complement_key(_permuted_bitset(key, permutation), len(permutation))
    complement, positions = key
    if isinstance(positions, int):
        image = _permuted_bitset(positions, permutation)
        if complement:
            return _facet_complement_key(image, len(permutation))
        nonincident = ((1 << len(permutation)) - 1) ^ image
        return _facet_complement_key(nonincident, len(permutation))
    return complement, tuple(sorted(permutation[position] for position in positions))


def _bit_positions(mask: int) -> Iterator[int]:
    """Yield set-bit positions of a nonnegative integer, in increasing order."""
    while mask:
        least_bit = mask & -mask
        yield least_bit.bit_length() - 1
        mask ^= least_bit


def _schreier_generator_indices(word: SchreierWord) -> Iterator[int]:
    """Yield a nested Schreier word's generator indices in action order."""
    indices = []
    while word is not None:
        word, generator_index = word
        indices.append(generator_index)
    yield from reversed(indices)


def _permuted_bitset(mask: int, permutation: array | tuple[int, ...]) -> int:
    """Apply a position permutation to an exact bitset using byte storage."""
    size = len(permutation)
    source = mask.to_bytes((size + 7) // 8, "little")
    image = bytearray(len(source))
    for byte_index, value in enumerate(source):
        while value:
            least_bit = value & -value
            position = 8 * byte_index + least_bit.bit_length() - 1
            target = permutation[position]
            image[target >> 3] |= 1 << (target & 7)
            value ^= least_bit
    return int.from_bytes(image, "little")


# The admissible squares in the normalized signature-(1,n) form are actual
# subsets of R, not a source-level string selector.
NONNEGATIVE_SQUARES = RealSet.unbounded_above_closed(0)
ZERO_SQUARES = RealSet.point(0)


def _binary_split_isotropic_basis(lattice: Lattices.ParentMethods) -> Matrix_integer_dense:
    r"""The two primitive null lines of a split binary form, as row vectors.

    The preamble's finite primitive-null-line enumeration replaces
    ``GetBasisIsotropicVectors_reduced`` in pinned
    ``src_lorentzian/two_dim_lorentzian.h:443–461``. The signs of its
    two representatives are immaterial to the eight sign/permutation choices.
    """
    # GetBasisIsotropicVectors_reduced clears denominators and divides by
    # the content of each rational null vector.  Normalize explicitly so
    # the index of the null-line sublattice is independent of which primitive
    # isotropic representatives the preamble returns.
    rays_by_line = {}
    for point in lattice.primitive_isotropic_vectors():
        row = _coordinate_row(lattice, point, SageQQ)
        denominator = SageZZ(1)
        for entry in row:
            denominator = denominator.lcm(entry.denominator())
        integral = vector(SageZZ, (denominator * entry for entry in row))
        content = SageZZ(integral[0]).gcd(integral[1])
        if content == 0:
            raise ArithmeticError("a null-line generator must be nonzero")
        primitive = vector(SageZZ, (entry // content for entry in integral))
        first_nonzero = next(entry for entry in primitive if entry)
        if first_nonzero < 0:
            primitive = -primitive
        rays_by_line[tuple(primitive)] = primitive
    rays = tuple(sorted(rays_by_line.values(), key=lambda row: tuple(row)))
    if len(rays) != 2:
        raise ArithmeticError("a split nondegenerate binary lattice has exactly two rational null lines")
    # The pinned reduced null-line basis is integral and primitive.  Use
    # that exact arithmetic type for the determinant/index comparison and
    # perform the rational inverse only where the isometry is constructed.
    return matrix(SageZZ, rays)


def _binary_split_isometries(source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
    r"""Port ``TwoDimIsotropic_TestEquivalence`` and
    ``TwoDimIsotropic_AutomGenerator`` (pinned two_dim_lorentzian.h:478–612).

    Enumerate precisely the eight signed permutations of the two primitive
    null lines; keep exactly those inducing integral isometries. For different
    lattices, each returned morphism has source and target as endpoints.
    """
    source_basis = _binary_split_isotropic_basis(source)
    target_basis = _binary_split_isotropic_basis(target)
    # A null-line permutation may have determinant different from one only
    # when the source and target null-line sublattices have different index.
    # Integral unimodular isometries preserve that index.
    if abs(source_basis.det()) != abs(target_basis.det()):
        return ()
    source_gram = _gram_matrix(source)
    target_gram = _gram_matrix(target)
    source_labels = tuple(source.module_generating_set())
    target_labels = tuple(target.module_generating_set())
    result: list[LatticeIsometryMethods] = []
    inverse_source_basis = source_basis.change_ring(SageQQ).inverse()
    # The pinned generator routine enumerates the eight signed null-line
    # permutations, retaining exactly the integral Gram-preserving ones.
    for swap, sign0, sign1 in product((0, 1), repeat=3):
        first, second = (1, 0) if swap else (0, 1)
        end_basis = matrix(SageZZ, [(2 * sign0 - 1) * target_basis.row(first), (2 * sign1 - 1) * target_basis.row(second)])
        rows = inverse_source_basis * end_basis
        if any(entry.denominator() != 1 for entry in rows.list()):
            continue
        if rows * target_gram * rows.transpose() != source_gram:
            continue
        images = {
            label: target.linear_combination({t_label: target.base_ring()(int(entry)) for t_label, entry in zip(target_labels, row, strict=True) if entry})
            for label, row in zip(source_labels, rows.rows(), strict=True)
        }
        result.append(source.Isom(target)(images))
    return tuple(result)


def _binary_split_automorphisms(lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
    r"""Generators of O(L) for a split binary L: TwoDimIsotropic_AutomGenerator (two_dim_lorentzian.h:491-530).

    The source keeps the integral signed null-line permutations and requires at least two of
    them (the identity and -1).
    """
    generators = _binary_split_isometries(lattice, lattice)
    if len(generators) < 2:
        raise ArithmeticError("a split binary lattice has at least the automorphisms 1 and -1")
    return generators


def binary_split_vectors_of_square(lattice: Lattices.ParentMethods, square: int) -> tuple[Lattices.ElementMethods, ...]:
    r"""All integral vectors of nonzero square in a split Lorentzian plane.

    This is the *unreduced* finite candidate list of the pinned
    EnumerateVectorFixedNorm_Factorization (two_dim_lorentzian.h:623-685),
    used by edgewalk.h:921-961 before testing orientation and ambient
    integrality. In particular it includes imprimitive vectors and every
    signed factor pair; orbit representatives cannot replace this list.

    If u,v are primitive null generators, k=det(u,v), c=b(u,v), then
    q(z)=(2c/k^2) det(z,v) det(u,z). Thus every integral solution
    arises from one factor pair of square*k^2/(2c), followed by the
    exact integral-preimage check in the pinned routine.
    """
    square = int(square)
    if not square:
        raise ValueError("the fixed-norm factorization requires nonzero square; zero has infinitely many integral vectors")
    if int(lattice.rank()) != 2 or tuple(map(int, lattice.signature_pair())) != (1, 1) or not lattice.is_isotropic():
        raise ValueError("fixed-norm null factorization requires a split Lorentzian plane")
    null_basis = _binary_split_isotropic_basis(lattice)
    first, second = null_basis.rows()
    determinant = SageZZ(null_basis.det())
    pairing = SageZZ((first * _gram_matrix(lattice) * second.column())[0])
    if not determinant or not pairing:
        raise ArithmeticError("primitive null generators must span a nondegenerate plane")
    product_value = SageQQ(square) * determinant**2 / (2 * pairing)
    if product_value.denominator() != 1:
        return ()
    value = SageZZ(product_value)
    # Row coordinates recover z=(m1*u+m2*v)/det(u,v).
    inverse_factor_forms = null_basis.transpose().change_ring(SageQQ) / determinant
    labels = tuple(lattice.module_generating_set())
    found: list[Lattices.ElementMethods] = []
    for divisor in SageZZ(abs(value)).divisors():
        for sign in (-1, 1):
            left = SageZZ(sign * divisor)
            right = value // left
            row = inverse_factor_forms * vector(SageQQ, (left, right))
            if any(entry.denominator() != 1 for entry in row):
                continue
            point = lattice.linear_combination({label: lattice.base_ring()(int(entry)) for label, entry in zip(labels, row, strict=True) if entry})
            if point.q() != lattice.base_ring()(square):
                raise ArithmeticError("fixed-norm factorization returned an incorrect square")
            found.append(point)
    return tuple(found)


def binary_split_primitive_vector_classes(lattice: Lattices.ParentMethods, square: int) -> tuple[Lattices.ElementMethods, ...]:
    r"""Complete primitive fixed-square orbits of a split binary Lorentzian lattice.

    Implements the rank-two branch of pinned TwoDimIsotropic_OrbitRepresentative
    (two_dim_lorentzian.h:613-761), including its zero-norm branch.
    EnumerateVectorFixedNorm_Factorization enumerates factor pairs of integral
    linear forms. Construct the forms intrinsically from the primitive null
    lines instead of making an arbitrary factorization normalization.

    For primitive null generators u,v, put k=det(u,v), c=b(u,v), and
    m1(z)=det(z,v), m2(z)=det(u,z). Then
    z=(m1(z)u+m2(z)v)/k and q(z)=(2c/k^2)m1(z)m2(z).
    Enumerate every signed factor pair of square*k^2/(2c) and keep its
    integral primitive preimage. For square=0 the primitive solutions are
    +u,-u,+v,-v. The pinned signed-null-line isometries generate O(L),
    and their finite action yields exactly the orbit representatives.
    """
    if int(lattice.rank()) != 2 or tuple(map(int, lattice.signature_pair())) != (1, 1):
        raise ValueError("the split binary orbit algorithm requires signature (1,1)")
    if not lattice.is_isotropic():
        raise ValueError("a split binary orbit algorithm requires two rational null lines")
    null_basis = _binary_split_isotropic_basis(lattice)
    first, second = null_basis.rows()
    determinant = SageZZ(null_basis.det())
    gram = _gram_matrix(lattice)
    pairing = SageZZ((first * gram * second.column())[0])
    if not determinant or not pairing:
        raise ArithmeticError("the primitive null vectors must span a nondegenerate rational plane")

    if not square:
        rows = (first, -first, second, -second)
    else:
        rows = tuple(_coordinate_row(lattice, point, SageZZ) for point in binary_split_vectors_of_square(lattice, square) if point.is_primitive())

    labels = tuple(lattice.module_generating_set())

    def element(row):
        return lattice.linear_combination({label: lattice.base_ring()(int(entry)) for label, entry in zip(labels, row, strict=True) if entry})

    candidates = {tuple(int(entry) for entry in row): element(row) for row in rows}
    if not candidates:
        return ()
    generators = _binary_split_automorphisms(lattice)
    # As in TwoDimIsotropic_OrbitRepresentative, traverse the full orbit
    # under the null-line isometry generators, not just the initial list.
    # The initial list is nevertheless exhaustive: leaving it is an error.
    # The pinned algorithm acts by integral matrices on coordinate rows.
    # Compute these matrices once, rather than repeatedly constructing
    # lattice elements and evaluating morphisms at every orbit edge.
    actions = tuple(
        matrix(
            SageZZ,
            [[int(_coordinate_row(lattice, generator(lattice.linear_combination({label: lattice.base_ring().one()})), SageZZ)[j]) for j in range(2)] for label in labels],
        )
        for generator in generators
    )
    unseen = set(candidates)
    representatives: list[Lattices.ElementMethods] = []
    while unseen:
        start = min(unseen)
        representatives.append(candidates[start])
        orbit = {start}
        pending = [start]
        for key in pending:
            row = vector(SageZZ, key)
            for action in actions:
                image_key = tuple(int(entry) for entry in row * action)
                if image_key not in candidates:
                    raise ArithmeticError("signed-factor enumeration is not isometry-invariant")
                if image_key not in orbit:
                    orbit.add(image_key)
                    pending.append(image_key)
        unseen.difference_update(orbit)
    return tuple(representatives)


def perfect_complex(
    lattice: Lattices.ParentMethods,
    allowed_squares: RealSet,
    *,
    retain_facet_orbits: bool = True,
) -> LorentzianPerfectComplex:
    r"""The perfect-domain complex of a lattice, memoized as IndefiniteForms.g:1173-1299 memoizes its groups.

    The traversal of a complex is deterministic data of the lattice. Consumers
    with the same facet-retention contract share their cells, stabilizers and
    side pairings; the orthogonal-group route uses a separate cache entry that
    discards per-facet Schreier records once each orbit has been processed.
    """
    return _cached_perfect_complex(lattice, allowed_squares, retain_facet_orbits)


@cache
def _cached_perfect_complex(
    lattice: Lattices.ParentMethods,
    allowed_squares: RealSet,
    retain_facet_orbits: bool,
) -> LorentzianPerfectComplex:
    return LorentzianPerfectComplex(lattice, allowed_squares, retain_facet_orbits=retain_facet_orbits)


def lorentzian_orthogonal_group_generators(lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
    r"""Generators of O(L), with the pinned split rank-two branch."""
    if int(lattice.rank()) == 1:
        label = next(iter(lattice.module_generating_set()))
        negative = lattice.linear_combination({label: -lattice.base_ring().one()})
        return (lattice.Isom(lattice)({label: negative}),)
    if int(lattice.rank()) == 2 and lattice.is_isotropic():
        return _binary_split_automorphisms(lattice)
    # The orthogonal-group consumer needs quotient cells, cell stabilizers
    # and one side pairing per facet orbit. It does not need the per-facet
    # Schreier records retained by the richer public perfect-domain adapter.
    # Use a separate cached traversal mode so the public complex keeps its
    # complete domain data if another consumer requests it later.
    return perfect_complex(
        lattice, NONNEGATIVE_SQUARES, retain_facet_orbits=False
    ).full_orthogonal_group_generators()


def lorentzian_isometry(source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
    r"""An isometry between two Lorentzian lattices, or None."""
    # polyhedral_common@a55fcb7b, lorentzian_perfect.h:1897-1910:
    # LORENTZ_TestEquivalenceMatrices_Reduced first compares
    # INDEF_FORM_Invariant_NonDeg (rank, determinant, inertia). Use the
    # actual invariant tuple rather than its hashed representation.
    if source.rank() != target.rank() or source.signature_pair() != target.signature_pair() or _gram_matrix(source).det() != _gram_matrix(target).det():
        return None
    if int(source.rank()) == 1:
        # The rank-one integral isometry exists only for equal Gram entries;
        # the invariant check above supplies precisely that condition.
        source_label = next(iter(source.module_generating_set()))
        target_label = next(iter(target.module_generating_set()))
        image = target.linear_combination({target_label: target.base_ring().one()})
        return source.Isom(target)({source_label: image})
    if int(source.rank()) == 2 and int(target.rank()) == 2:
        if source.is_isotropic() != target.is_isotropic():
            return None
        if source.is_isotropic():
            return next(iter(_binary_split_isometries(source, target)), None)
    # LORENTZ_TestEquivalenceMatrices_Reduced (lorentzian_perfect.h:1887-1975)
    # traverses the complex of one form and tests each new cell against one
    # perfect cell of the other. Isom(L, M) is empty exactly when Isom(M, L)
    # is, so traverse whichever complex this session has already explored
    # further and invert its witness.
    source_complex = perfect_complex(source, NONNEGATIVE_SQUARES)
    target_complex = perfect_complex(target, NONNEGATIVE_SQUARES)
    if target_complex.explored_cell_count() > source_complex.explored_cell_count():
        witness = target_complex.isometry_to(source_complex)
        return None if witness is None else witness.inverse()
    return source_complex.isometry_to(target_complex)


def _gram_matrix(lattice: Lattices.ParentMethods) -> Matrix_integer_dense:
    return _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)


@dataclass(frozen=True)
class LorentzianPerfectCell:
    r"""A perfect cell, held as the coordinate rows of its vector configuration.

    The configuration is the cell's defining datum; its lattice elements are
    raised on first request.
    """

    lattice: Lattices.ParentMethods
    rows: tuple[tuple[int, ...], ...]
    allowed_squares: RealSet

    @cached_property
    def vector_configuration(self) -> tuple[Lattices.ElementMethods, ...]:
        return tuple(_ambient_element(self.lattice, row) for row in self.rows)

    def __post_init__(self) -> None:
        if self.allowed_squares not in (NONNEGATIVE_SQUARES, ZERO_SQUARES):
            raise ValueError(f"unsupported normalized-norm locus {self.allowed_squares!r}")
        positive, negative = self.lattice.signature_pair()
        if int(positive) != 1 and int(negative) != 1:
            raise ValueError(f"a Lorentzian perfect cell needs signature (1,n) or (n,1), but {self.lattice} has signature {self.lattice.signature_pair()}")
        gram = _gram_matrix(self.lattice)
        rows = tuple(vector(SageZZ, row) for row in self.rows)
        squares = tuple(int(row * gram * row) for row in rows)
        normalized_sign = 1 if int(positive) == 1 else -1
        if any(normalized_sign * square not in self.allowed_squares for square in squares):
            raise ValueError(f"a perfect cell may contain only vectors with normalized square in {self.allowed_squares}")
        # Lorentzian.g:303-312 (LORENTZ_CheckCorrectnessVectorFamily):
        # the augmented vectors (1,v) span an affine hyperplane of rank n.
        # Both GetOnePerfect and DoFlipping certify this before returning.
        augmented = matrix(SageQQ, [[1, *row] for row in rows])
        if augmented.rank() != gram.nrows():
            raise ArithmeticError("a perfect-cell vector family must have augmented rank equal to the lattice rank")


class _CertifiedPerfectDomainOrbitAdjacency(PerfectDomainOrbitAdjacency):
    """Side pairing certified by exact ray and supporting-covector data.

    The generic adjacency constructor enumerates all faces of both cells to
    check the common-face condition. On cells with 100+ facets, repeating
    that enumeration for every side pairing dominates traversal time. This
    constructor verifies the same face conditions by exact rank, membership,
    and supporting-covector checks. The actual neighbor and common-face cones
    are materialized lazily; when requested, the neighbor is checked against
    the exact primitive extreme rays transformed by the integral action.
    """

    def __init__(self, source, target, target_to_neighbor_factory, *, facet_normal_row, face_rows, target_halfspace_rows, action, neighbor_factory, common_face_factory):
        lattice = source.ambient_lattice()
        if target.ambient_lattice() is not lattice:
            raise ValueError("a certified perfect-domain adjacency must lie in one ambient lattice")
        rows = matrix(SageQQ, face_rows)
        normal = vector(SageQQ, facet_normal_row)
        # `face_rows` comes from the exact incidence key of this minimized
        # H-facet, so its codimension and vanishing on these rays are certified
        # when the incidence is built. Batch the remaining opposite-side
        # check instead of constructing a Sage vector and dot product per ray.
        source_values = matrix(SageQQ, source._engine.rays_list()) * normal.column()
        target_rays = matrix(SageZZ, target._engine.rays_list()) * action
        neighbor_values = matrix(SageQQ, target_rays) * normal.column()
        source_values = tuple(source_values.column(0))
        neighbor_values = tuple(neighbor_values.column(0))
        source_positive = bool(source_values) and all(value >= 0 for value in source_values) and any(source_values)
        source_negative = bool(source_values) and all(value <= 0 for value in source_values) and any(source_values)
        neighbor_positive = bool(neighbor_values) and all(value >= 0 for value in neighbor_values) and any(neighbor_values)
        neighbor_negative = bool(neighbor_values) and all(value <= 0 for value in neighbor_values) and any(neighbor_values)
        if not ((source_positive and neighbor_negative) or (source_negative and neighbor_positive)):
            raise ArithmeticError("the common-face covector does not support the two cells on opposite sides")
        inverse_action = action.inverse().change_ring(SageQQ)
        transformed_face_rows = rows * inverse_action
        if target.equation_covectors():
            # The usual perfect-cell cones are full-dimensional. Preserve
            # the exact equation test for a lower-dimensional target too.
            contained = all(
                target._engine.contains(tuple(row))
                for row in transformed_face_rows.rows()
            )
        else:
            face_values = transformed_face_rows * target_halfspace_rows.transpose()
            contained = all(value >= 0 for value in face_values.list())
        if not contained:
            raise ArithmeticError("an incident source ray is not contained in the proposed neighbor")
        self._source = source
        self._target = target
        self._neighbor = None
        self._neighbor_factory = neighbor_factory
        self._common_face = None
        self._common_face_factory = common_face_factory
        self._target_to_neighbor = None
        self._target_to_neighbor_factory = target_to_neighbor_factory
        self._action = action

    def neighbor(self):
        """Materialize the exact adjacent cone only when a consumer requests it."""
        if self._neighbor is None:
            neighbor = self._neighbor_factory()
            target_rays = matrix(SageZZ, self._target._engine.rays_list()) * self._action
            expected_rays = {tuple(int(entry) for entry in row) for row in target_rays.rows()}
            actual_rays = {tuple(int(entry) for entry in row) for row in neighbor._engine.rays_list()}
            if expected_rays != actual_rays:
                raise ArithmeticError("the supplied isometry does not carry the target cell onto the proposed neighbor")
            self._neighbor = neighbor
        return self._neighbor

    def common_face(self):
        """Materialize the common face cone only when a consumer requests it."""
        if self._common_face is None:
            self._common_face = self._common_face_factory()
        return self._common_face

    def target_to_neighbor(self):
        """Construct and cache the owned lattice isometry only when requested."""
        if self._target_to_neighbor is None:
            morphism = self._target_to_neighbor_factory()
            lattice = self._source.ambient_lattice()
            basis = tuple(lattice.module_generator(label) for label in lattice.module_generating_set())
            morphism_action = matrix(
                SageZZ,
                [_coordinate_row(lattice, morphism(vector), SageZZ) for vector in basis],
            )
            if morphism_action != self._action:
                raise ArithmeticError("the Schreier permutation and categorical side-pairing actions disagree")
            self._target_to_neighbor = morphism
        return self._target_to_neighbor


class _CertifiedPerfectDomainTraversal(LorentzianPerfectDomainTraversal):
    """Traversal whose exact facet-incidence completeness check is already discharged.

    The preamble's default ``unpaired_facets`` recomputes every facet cone and
    compares it with every common-face cone.  The port has already checked
    the identical incidence condition using the minimized H-presentations,
    including the totally isotropic exterior-boundary exception.  Reuse that
    certificate instead of repeating the expensive categorical construction
    whenever a consumer asks whether the traversal is complete.
    """

    def unpaired_facets(self):
        return finite_ordered_set(())


@dataclass(frozen=True)
class LorentzianCellAdjacency:
    source: LorentzianPerfectCell
    facet: FacetIncidence
    target: LorentzianPerfectCell
    transporter: LatticeIsometryMethods | None


class _ConfigurationPairings:
    r"""The Gram pairings b(v_i, v_j) of a vector configuration, read from its n x rank rows.

    Rows of the n x n pairing matrix are produced in blocks and single
    entries and columns on request, so the matrix is never stored.
    """

    def __init__(
        self,
        rows: tuple[tuple[int, ...], ...] | Matrix_integer_dense,
        gram: Matrix_integer_dense,
        *,
        cache_columns: bool = False,
    ) -> None:
        # lorentzian_perfect.h:1212-1256 uses exact scalar products without
        # a coefficient bound. Sage's integer matrices retain that contract
        # for both configuration equivalence and perfect-form invariants.
        # Configuration-invariant consumers already own a dense coordinate
        # matrix; share it instead of copying the full n-by-r array again.
        self.rows = rows if isinstance(rows, Matrix_integer_dense) else matrix(SageZZ, rows)
        self.covectors = self.rows * gram.change_ring(SageZZ)
        # Backtracking repeatedly asks for the same anchor columns while
        # comparing partial configurations. Enable a small bounded cache only
        # for that workload: streaming consumers such as perfect-form hashing
        # visit each column once and must not retain full-length Python tuples.
        self._cache_columns = cache_columns
        self._column_cache: dict[int, tuple[int, ...]] = {}
        self._column_cache_limit = max(1, min(64, 250_000 // max(1, self.rows.nrows())))

    def __len__(self) -> int:
        return self.rows.nrows()

    def __getitem__(self, index: tuple[int, int]) -> int:
        left, right = index
        return int(self.covectors.row(left).dot_product(self.rows.row(right)))

    def diagonal(self):
        return tuple(value for block in self.diagonal_blocks() for value in block)

    def diagonal_blocks(self):
        """Yield exact diagonal Gram entries in bounded row blocks."""
        # Compute the diagonal of the Gram matrix with FLINT-backed
        # elementwise products, rather than one matrix-row extraction and dot
        # product per configuration vector. Bound the temporary product by
        # entry count: large configurations must not allocate another full
        # copy of their coordinate matrix merely to obtain the diagonal.
        row_count, rank = self.rows.dimensions()
        block_size = max(1, min(1024, 100_000 // max(1, rank)))
        for start in range(0, row_count, block_size):
            stop = min(start + block_size, row_count)
            products = self.covectors.matrix_from_rows(range(start, stop)).elementwise_product(
                self.rows.matrix_from_rows(range(start, stop))
            )
            yield tuple(sum(row) for row in products.rows())

    def column(self, index: int, *, cache: bool | None = None):
        # The form is symmetric, so column j holds b(v_l, v_j) for every l.
        # A one-column integer matrix selects Sage's FLINT multiplication;
        # Matrix._matrix_times_vector_ instead builds and sums every column.
        use_cache = self._cache_columns if cache is None else cache
        cached = self._column_cache.get(index)
        if use_cache and cached is not None:
            return cached
        covector = self.covectors.matrix_from_rows([index]).transpose()
        result = tuple(int(entry) for entry in (self.rows * covector).column(0))
        if not use_cache:
            return result
        if len(self._column_cache) >= self._column_cache_limit:
            self._column_cache.pop(next(iter(self._column_cache)))
        self._column_cache[index] = result
        return result

    def blocks(self):
        # ComputeInvariantPerfectForm streams pairings (lorentzian_perfect.h:
        # 1223-1244). Multiply bounded row blocks using FLINT rather than
        # issuing one matrix-vector multiplication per row. A 16-by-n block
        # remains small for the 49,861-ray A7 configuration, while reducing
        # Python/FLINT dispatches by a factor of sixteen.
        # Bound the temporary exact matrix by entry count as well as row
        # count: configurations larger than the A7 case must not turn a fixed
        # row block into an unbounded memory allocation.
        block_size = _configuration_pairing_block_size(len(self))
        transposed_rows = self.rows.transpose()
        for start in range(0, len(self), block_size):
            stop = min(start + block_size, len(self))
            pairing_block = self.covectors.matrix_from_rows(range(start, stop)) * transposed_rows
            yield start, np.array(pairing_block.rows(), dtype=object)


def _configuration_pairing_block_size(configuration_size: int) -> int:
    r"""Target one million entries per exact pairing block, capped at 16 rows.

    A single row may itself exceed the target for configurations with more
    than one million vectors, so the row block has a minimum size of one.
    """
    if configuration_size < 0:
        raise ValueError("a configuration size cannot be negative")
    return max(1, min(16, 1_000_000 // max(1, configuration_size)))


def _iter_owned_row_blocks(rows, block_size: int) -> Iterator[tuple[tuple[int, ...], ...]]:
    """Yield immutable row blocks and release consumed Normaliz list entries."""
    if block_size < 1:
        raise ValueError("row block size must be positive")
    for start in range(0, len(rows), block_size):
        stop = min(len(rows), start + block_size)
        block = tuple(tuple(int(entry) for entry in rows[position]) for position in range(start, stop))
        if isinstance(rows, list):
            for position in range(start, stop):
                rows[position] = None
        yield block


# The full weighted pairing graph has quadratic size even though an integral
# isometry is determined by the images of a rank-sized basis. Large frozen
# configurations use only intrinsic vertex colours here and retain exact
# pairwise constraints in the basis-image search below.
_PAIRING_GRAPH_REFINEMENT_LIMIT = 4096


def _configuration_signature_multisets_match(
    source_colours: tuple[int, ...],
    source_columns: tuple[tuple[int, ...], ...],
    target_colours: tuple[int, ...],
    target_columns: tuple[tuple[int, ...], ...],
    *,
    source_signature_counts: Counter[tuple[int, tuple[int, ...]]] | None = None,
) -> bool:
    """Compare row-signature multisets with one incrementally cancelled counter.

    Retaining a nested signature tuple for every configuration row at every
    recursion depth multiplies memory by the ambient rank. Keep only the
    assigned pairing columns in the search stack. During each comparison,
    cancel matching target signatures immediately so equal configurations do
    not require two simultaneous counters of all row signatures.
    """
    if len(source_colours) != len(target_colours) or len(source_columns) != len(target_columns):
        return False
    counts = source_signature_counts
    if counts is None:
        counts = Counter(
            (colour, tuple(column[index] for column in source_columns))
            for index, colour in enumerate(source_colours)
        )
    for index, colour in enumerate(target_colours):
        signature = (colour, tuple(column[index] for column in target_columns))
        counts[signature] -= 1
        if counts[signature] < 0:
            # Do not retain a second list of every row signature just to
            # restore the reusable source counter. Recompute the already
            # visited target prefix; mismatches typically terminate early.
            for previous_index in range(index + 1):
                previous = (
                    target_colours[previous_index],
                    tuple(column[previous_index] for column in target_columns),
                )
                counts[previous] += 1
                if counts[previous] == 0:
                    del counts[previous]
            return False
    # Successful comparisons scan the target rows a second time to restore
    # counts, trading O(n * rank) tuple construction for no O(n) side list.
    for index, colour in enumerate(target_colours):
        signature = (colour, tuple(column[index] for column in target_columns))
        counts[signature] += 1
        if counts[signature] == 0:
            del counts[signature]
    return True


def _packed_pairing_fibres(
    colours: tuple[int, ...],
    pairings: tuple[int, ...],
    index_typecode: str,
) -> dict[tuple[int, int], array]:
    """Partition target positions by exact colour and pairing using packed indices."""
    if len(colours) != len(pairings):
        raise ValueError("pairing fibres require one colour per target row")
    if index_typecode not in ("I", "Q"):
        raise ValueError("pairing-fibre indices require an unsigned packed integer type")
    fibres: dict[tuple[int, int], array] = {}
    for index, (colour, pairing) in enumerate(zip(colours, pairings, strict=True)):
        fibres.setdefault((colour, int(pairing)), array(index_typecode)).append(index)
    return fibres


def _row_multiset(row, own_value: int) -> tuple[tuple[int, int], ...]:
    r"""The multiset of a pairing row's entries other than the vertex's own, as sorted (value, count)."""
    # NumPy's unique/sort on object arrays is not an integral-arithmetic
    # primitive.  Counter retains arbitrary-precision codes without casting
    # them back to a fixed-width dtype; subtract exactly one self-pairing.
    counts = Counter(int(value) for value in row)
    counts[own_value] -= 1
    if counts[own_value] < 0:
        raise ArithmeticError("the diagonal pairing is absent from its row")
    return tuple(sorted((value, count) for value, count in counts.items() if count))


def _refined_signatures(pairings: _ConfigurationPairings, colours: tuple[int, ...], modulus: int) -> tuple[tuple[int, tuple[tuple[int, int], ...]], ...]:
    r"""Each vertex's colour with the multiset of (pairing, neighbour colour), coded as pairing * modulus + colour."""
    if modulus <= 0 or len(colours) != len(pairings):
        raise ValueError("refinement requires one colour per configuration vector and a positive modulus")
    # Keep both pairings and colour encodings in arbitrary-precision arithmetic.
    colour_array = np.array(colours, dtype=object)
    signatures = []
    for start, block in pairings.blocks():
        codes = block.astype(object) * modulus + colour_array[None, :]
        for offset, row in enumerate(codes):
            index = start + offset
            signatures.append((colours[index], _row_multiset(row, int(row[index]))))
    return tuple(signatures)


class LorentzianPerfectLocalBackend:
    def _normaliz_cell_cone(self, cell: LorentzianPerfectCell):
        """Cache Normaliz's exact cone until its native facet vector is consumed.

        Unlike ``cached_method``, this small explicit cache can release a cone
        after ``NmzResult`` has copied its support hyperplanes into Python. The
        cone is still shared with an extant polyhedron wrapper, and repeated
        ordinary calls reuse it until a raw H-result is requested.
        """
        from sage.geometry.polyhedron.backend_normaliz import NmzCone

        cache = self.__dict__.setdefault("_normaliz_cell_cones", {})
        if cell in cache:
            return cache[cell]
        dimension = int(cell.lattice.module_rank())
        # PyNormaliz accepts immutable tuple rows directly. Passing the
        # canonical configuration avoids duplicating all 49,861 ray rows as
        # Python lists while Normaliz constructs its own exact cone input.
        cone = NmzCone(
            vertices=[(0,) * dimension + (1,)],
            cone=cell.rows,
            subspace=(),
        )
        cache[cell] = cone
        return cone

    def _release_normaliz_cell_cone(self, cell: LorentzianPerfectCell, cone) -> None:
        """Drop the backend's cached reference after raw H-result extraction."""
        cache = self.__dict__.get("_normaliz_cell_cones", {})
        if cache.get(cell) is cone:
            del cache[cell]

    def _normaliz_cell_polyhedron(self, cell: LorentzianPerfectCell):
        """Temporarily wrap the cached Normaliz cone in Sage's polyhedron interface.

        Do not cache this wrapper: Sage's Normaliz polyhedron retains Python
        V/H-representation objects in addition to the engine cone. The exact
        engine is already cached separately, and callers either convert the
        result into the owned cone (itself cached) or consume its ray rows
        immediately. Retaining the wrapper duplicates the large A7 H-list
        alongside the owned cone and the raw incidence iterator's result.
        """
        from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
        from sage.geometry.polyhedron.parent import Polyhedra_normaliz

        parent = Polyhedra_normaliz(SageQQ, int(cell.lattice.module_rank()), "normaliz")
        return Polyhedron_normaliz(
            parent,
            None,
            None,
            normaliz_cone=self._normaliz_cell_cone(cell),
            internal_base_ring=SageQQ,
        )

    @cached_method
    def cell_cone(self, cell: LorentzianPerfectCell):
        if len(cell.rows) < 10_000:
            return RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)

        # PPL's generic V-to-H conversion is the current A7 bottleneck: its
        # initial perfect configuration has 49,861 rays. Normaliz computes
        # the exact dual description. Build the owned cone directly from its
        # raw support covectors rather than first materializing Sage's full
        # Polyhedron_normaliz inequality/equation objects and then converting
        # those objects back into owned covectors.
        # Request H- and equation-data from the same native cone. This avoids
        # a second complete support-hyperplane query and a dense QQ copy of all
        # ray rows solely to recover the linear span's annihilator. The raw
        # coordinate rows are converted directly to the owned covectors below.
        hyperplanes, equations = self._normaliz_cell_h_presentation(
            cell, include_equations=True
        )
        return _rational_polyhedral_cone(
            cell.lattice,
            self._iter_normaliz_owned_covectors(cell, hyperplanes),
            equation_covectors=self._iter_normaliz_owned_covectors(cell, equations),
        )

    def _normaliz_cell_h_presentation(self, cell: LorentzianPerfectCell, *, include_equations: bool):
        """Extract exact H-data directly in the cone's linear coordinates."""
        from sage.geometry.polyhedron.backend_normaliz import NmzCone, NmzResult

        # Use the ray-cone input directly. The cached Normaliz cone used by
        # Polyhedron_normaliz includes a homogenized origin vertex, which adds
        # an irrelevant constant coordinate and tautological H-inequality.
        # This path needs only the native cone result, not the polyhedron
        # wrapper, and must preserve the lean rank-dimensional H rows.
        normaliz_cone = NmzCone(cone=cell.rows, subspace=())
        try:
            hyperplanes = NmzResult(normaliz_cone, "SupportHyperplanes")
            equations = NmzResult(normaliz_cone, "Equations") if include_equations else ()
        finally:
            del normaliz_cone
        return hyperplanes, equations

    def _iter_normaliz_owned_covectors(self, cell: LorentzianPerfectCell, rows):
        """Convert rank-dimensional Normaliz rows to owned covectors in place."""
        for index, row in enumerate(rows):
            if isinstance(rows, list):
                rows[index] = None
            if all(entry == 0 for entry in row):
                del row
                continue
            covector = _integral_covector(cell.lattice, row)
            del row
            yield covector

    def _iter_cell_halfspace_pairing_rows(
        self,
        cell: LorentzianPerfectCell,
        labels: tuple[int, ...],
        *,
        use_owned_cone: bool = True,
    ) -> Iterator[tuple[int, ...]]:
        """Yield exact facet pairing rows in the scale returned by Normaliz.

        The incidence computation only tests whether each pairing is zero;
        primitive normalization is therefore unnecessary here. Full-domain
        consumers reuse the owned cone's minimized H-family; group-only
        traversal requests raw Normaliz rows to avoid constructing a PPL cone
        that it will never consume.
        """
        if not use_owned_cone:
            # Orthogonal-group generation consumes exact facet incidences and
            # flips them from the original ray configuration. It does not
            # need the owned cone object, so read the raw Normaliz H-result
            # directly and avoid constructing a PPL-side representation of
            # the complete support-hyperplane family.
            hyperplanes, _equations = self._normaliz_cell_h_presentation(
                cell, include_equations=False
            )
            for index, row in enumerate(hyperplanes):
                if isinstance(hyperplanes, list):
                    hyperplanes[index] = None
                coordinates = tuple(int(entry) for entry in row)
                del row
                if any(coordinates):
                    yield coordinates
            return
        if len(cell.rows) >= 1_000:
            # The traversal's first nonboundary facet requires cell_cone for
            # flipping. Re-querying Normaliz here creates a second complete
            # support-hyperplane result while that owned cone is then built.
            # Reuse its exact minimized H-family instead: it is already cached
            # and is the same facet family, with covectors normalized only by
            # a nonzero integral scale that does not change incidence zeros.
            covectors = self.cell_cone(cell).halfspace_covectors()
            for covector in covectors:
                yield tuple(int(covector(label)) for label in labels)
            return
        # Small owned cones are reused directly rather than querying Normaliz
        # for a second H-presentation.
        for covector in self.cell_cone(cell).halfspace_covectors():
            yield tuple(int(covector(label)) for label in labels)

    def _iter_cell_halfspace_covectors(self, cell: LorentzianPerfectCell):
        """Yield the owned cone's exact primitive facet covectors."""
        yield from self.cell_cone(cell).halfspace_covectors()

    def cell_halfspace_covectors(self, cell: LorentzianPerfectCell):
        """Return the exact primitive homogeneous facet covectors of ``cell``.

        The owned cone already stores this exact tuple. Return it directly
        rather than constructing and caching a second tuple with a duplicate
        pointer array for large facet families.
        """
        return self.cell_cone(cell).halfspace_covectors()

    @cached_method
    def cell_primitive_rays(self, cell: LorentzianPerfectCell, square: int | None = None):
        """Return exact primitive extreme rays without rebuilding PPL.

        Normaliz already owns the V-representation while computing the
        minimized H-presentation. Use its exact extreme rays for large cells
        rather than constructing an owned cone from the H-presentation and
        converting that presentation back to V merely to find isotropic marks.
        If ``square`` is specified, filter in integer coordinates before
        constructing owned lattice elements.
        """
        if len(cell.rows) < 1_000:
            rays = tuple(self.cell_cone(cell).primitive_rays())
            return rays if square is None else tuple(ray for ray in rays if int(ray.q()) == square)
        from sage.geometry.polyhedron.backend_normaliz import NmzCone, NmzResult

        # Query the exact cone generators directly. Polyhedron_normaliz would
        # first allocate V/H-representation objects for every ray and facet,
        # although this caller only needs the integer extreme-ray rows.
        normaliz_cone = NmzCone(cone=cell.rows, subspace=())
        ray_rows = NmzResult(normaliz_cone, "ExtremeRays")
        del normaliz_cone
        gram = _gram_matrix(cell.lattice) if square is not None else None
        owned_rays = []
        for index, row in enumerate(ray_rows):
            if isinstance(ray_rows, list):
                ray_rows[index] = None
            coordinates = tuple(int(entry) for entry in row)
            del row
            if gram is not None:
                ray_vector = vector(SageZZ, coordinates)
                if int(ray_vector * gram * ray_vector) != square:
                    continue
            owned_rays.append(_ambient_element(cell.lattice, coordinates))
        if ray_rows:
            del coordinates
        del ray_rows
        return tuple(owned_rays)

    def initial_cell(self, lattice: Lattices.ParentMethods, allowed_squares: RealSet) -> LorentzianPerfectCell:
        normalized = _lorentzian_view(lattice)
        direction = _positive_direction(normalized)
        # The search runs on coordinate rows; L and its normalized twist share them.
        critical = _search_initial_vectors(normalized, direction, allowed_squares)
        gram = _gram_matrix(normalized).change_ring(SageQQ)
        central = vector(SageQQ, direction)
        first = vector(SageQQ, critical[0])
        scalar = (central * gram * first.column())[0]
        base_normal = vector(SageQQ, [-scalar, *(gram * central.column()).column(0)])
        while matrix(SageQQ, critical).rank() < int(lattice.module_rank()):
            expanded = matrix(SageQQ, [[1, *row] for row in critical])
            # Lorentzian.g uses NullspaceMat(TransposedMat(EXT)); retain GAP's
            # source-ordered row basis because each row feeds a separately
            # normalized candidate direction in GetOneOutsideRay.
            nullspace = _pinned_nullspace_matrix(expanded.transpose())
            # Lorentzian.g:1179-1188 applies RemoveFraction to each
            # covector-derived direction separately before GetOneOutsideRay.
            # Normalizing only the span basis globally changes the source's
            # candidate rows and therefore its negative-vector search order.
            directions = _outside_ray_directions(nullspace, gram)
            outside_direction = _negative_direction_in_span(normalized, directions, critical, base_normal)
            direction_covector = gram * outside_direction.column()
            critical_row = vector(SageQQ, critical[0])
            direction_scalar = (critical_row * direction_covector)[0]
            normal_direction = vector(SageQQ, [-direction_scalar, *direction_covector.column(0)])
            critical, flipped_normal, _test_direction, _max_scal = _kernel_flipping(
                normalized, critical, base_normal, normal_direction, allowed_squares
            )
            # Lorentzian.g:1194: normalize the returned supporting normal
            # before the next wall-crossing step.
            base_normal = _source_remove_fraction(flipped_normal)
        return LorentzianPerfectCell(lattice, critical, allowed_squares)

    @cached_method
    def cell_stabilizer(self, cell: LorentzianPerfectCell) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of the finite stabilizer of a spanning perfect cell.

        Every stabilizer element is determined by the images of any rational
        basis drawn from the configuration. The exact integral-transporter
        search grows a certified permutation subgroup as it discovers new
        isometries, reducing the remaining basis-image search modulo the
        pointwise stabilizer of the already assigned images. This returns
        generators without first materializing the full stabilizer group.
        """
        return self.cell_transporter(cell, cell, all_images=True, generators_only=True)

    def point_permutations(self, cell: LorentzianPerfectCell, points: tuple[Lattices.ElementMethods, ...]) -> tuple[tuple[int, ...], ...]:
        r"""Each cell-stabilizer generator's permutation of a finite stable family of points, as position images."""
        lattice = cell.lattice
        coordinate_rows = tuple(tuple(int(entry) for entry in _coordinate_row(lattice, point, SageZZ)) for point in points)
        return self.point_permutations_from_rows(cell, coordinate_rows)

    def point_permutations_from_rows(
        self,
        cell: LorentzianPerfectCell,
        coordinate_rows: tuple[tuple[int, ...], ...],
    ) -> tuple[tuple[int, ...], ...]:
        r"""The same permutation certificate for points already in lattice coordinates."""
        lattice = cell.lattice
        labels = tuple(lattice.module_generating_set())
        positions = {row: position for position, row in enumerate(coordinate_rows)}
        if len(positions) != len(coordinate_rows):
            raise ArithmeticError("a stable point configuration contains repeated points")
        generators = self.cell_stabilizer(cell)
        if not coordinate_rows:
            # The unique action on the empty set is the empty permutation.
            # Sage represents matrix(ZZ, ()) as 0-by-0, which cannot be
            # multiplied by the nonempty ambient action matrix.
            return tuple(() for _generator in generators)
        point_matrix = matrix(SageZZ, coordinate_rows)
        permutations = []
        for generator in generators:
            # The rows of this matrix are g(e_i). Linearity therefore gives
            # g(x) = x * action for every row of point_matrix. Evaluate the
            # categorical isometry only on the ambient basis, rather than
            # separately on every vertex of a large perfect cell.
            action = matrix(SageZZ, [_coordinate_row(lattice, generator(lattice.module_generator(label)), SageZZ) for label in labels])
            images = point_matrix * action
            try:
                permutation = tuple(positions[tuple(int(entry) for entry in row)] for row in images.rows())
            except KeyError as exc:
                raise ArithmeticError("a purported cell-stabilizer generator fails to preserve its point configuration") from exc
            permutations.append(permutation)
        return tuple(permutations)

    def facet_orbit_permutations(self, cell: LorentzianPerfectCell) -> tuple[tuple[tuple[FacetIncidence, tuple[int, ...]], ...], ...]:
        r"""Orbits of Aut((L, S)) on facets, paired with their exact ray-permutation transversals.

        Lorentzian.g:1318-1322 (PermutationStabilizer := GRP_int) and
        OnSetsGroupFormalism: the stabilizer acts on facets through its
        permutation action on the cell's rays. The Schreier transversals are
        retained as ray permutations; ``facet_orbits`` materializes the owned
        lattice isometries only for callers requiring categorical morphisms.
        """
        # Work from the compact incidence keys rather than first expanding and
        # caching a second complete collection of large frozenset incidences.
        # This public result necessarily materializes each incidence and full
        # ray permutation, but the intermediate word-orbit view is not needed.
        # Do not cache the materialized result: it can contain one full ray
        # permutation per facet, and retaining it after the caller has consumed
        # it would defeat the bounded-memory traversal used to construct it.
        mask_orbits, generator_permutations = self._facet_orbit_transversal_masks(cell)
        ray_count = len(cell.rows)
        permutation_typecode = (
            "H" if ray_count <= 2**16 else "I" if ray_count <= 2**32 else "Q"
        )
        def materialize_orbit(orbit):
            def materialize_element(item):
                key, word = item
                # Compose the Schreier word in packed storage. The public
                # return type requires a full tuple per facet, but keeping all
                # intermediate compositions as Python tuples multiplies the
                # peak allocation by the word length.
                permutation = array(permutation_typecode, range(ray_count))
                for generator_index in _schreier_generator_indices(word):
                    generator = generator_permutations[generator_index]
                    permutation = array(
                        permutation_typecode,
                        (generator[position] for position in permutation),
                    )
                return (_facet_incidence(key, ray_count), tuple(permutation))

            def consume_orbit():
                for index, item in enumerate(orbit):
                    if item is None:
                        raise ArithmeticError("a facet orbit entry was consumed more than once")
                    yield materialize_element(item)
                    # The BFS orbit is private and has no further consumer.
                    # Drop its key/Schreier-word references as each public
                    # result is committed to the output tuple.
                    orbit[index] = None

            return tuple(consume_orbit())

        return tuple(materialize_orbit(orbit) for orbit in mask_orbits)

    def _facet_mask_data(
        self,
        cell: LorentzianPerfectCell,
        *,
        use_owned_cone: bool = True,
    ) -> tuple[
        Iterator[FacetMaskKey],
        tuple[array, ...],
        dict[FacetMaskKey, FacetMaskKey | None] | None,
    ]:
        r"""Yield exact facet keys and the ray permutations of stabilizer generators.

        Facet keys are generated in H-presentation order. Keeping this as an
        iterator avoids a second dictionary containing every facet before the
        orbit traversal starts. The raw-H path necessarily retains a canonical
        lookup map, which the group-only orbit traversal also uses for its
        visited marks; the full-domain path uses a separate visited set.
        """
        ray_count = len(cell.rows)
        labels = tuple(cell.lattice.module_generating_set())
        rank = len(cell.rows[0])
        # Keep stabilizer permutation construction blockwise. Materializing
        # both R^T and A^T R^T duplicates the entire 49,861-ray configuration
        # as dense Sage matrices before any facet processing starts.
        permutation_block_size = _configuration_pairing_block_size(ray_count)
        generators = self.cell_stabilizer(cell)
        if not use_owned_cone and ray_count > _PAIRING_GRAPH_REFINEMENT_LIMIT:
            # The stabilizer computation has finished using its dense
            # row/covector matrices. They are useful again when comparing the
            # flipped neighbor, but retaining both while building the packed
            # ray permutations and Normaliz H-result needlessly raises the
            # peak. Drop this one identity-LRU entry; the first subsequent
            # transporter rebuilds it and then reuses it normally.
            invariant_cache = self.__dict__.get("_configuration_invariant_cache", {})
            cached = invariant_cache.get(id(cell))
            if cached is not None and cached[0] is cell:
                del invariant_cache[id(cell)]
                del cached
        # The canonical configuration rows are already immutable tuples, so
        # use them directly as dictionary keys rather than allocating a second
        # 49,861-tuple copy of the full A7 ray configuration.
        positions = {row: index for index, row in enumerate(cell.rows)}
        if len(positions) != ray_count:
            raise ArithmeticError("a perfect-cell configuration contains repeated rays")
        # Store ray-position permutations in packed unsigned words. A Python
        # tuple uses an eight-byte pointer per ray, and these full permutations
        # remain live throughout facet-orbit traversal. Choose a word wide
        # enough for the largest possible position (ray_count - 1).
        permutation_typecode = (
            "H" if ray_count <= 2**16 else "I" if ray_count <= 2**32 else "Q"
        )
        permutations: list[array] = []
        for generator in generators:
            action = matrix(
                SageZZ,
                [
                    _coordinate_row(
                        cell.lattice,
                        generator(cell.lattice.module_generator(label)),
                        SageZZ,
                    )
                    for label in labels
                ],
            )
            permutation = array(permutation_typecode)
            try:
                for start in range(0, ray_count, permutation_block_size):
                    stop = min(start + permutation_block_size, ray_count)
                    image_rows = (matrix(SageZZ, cell.rows[start:stop]) * action).rows()
                    permutation.extend(positions[tuple(row)] for row in image_rows)
                    del image_rows
                # Keep the packed buffer: converting it back to a tuple here
                # defeats the memory reduction while retaining both buffers
                # briefly during conversion.
                permutations.append(permutation)
            except KeyError as exc:
                raise ArithmeticError("a cell-stabilizer generator does not preserve its rays") from exc

        # Facet-key generation can dominate the live memory for a large
        # configuration. Release the row-position index and the final loop
        # temporaries before asking Normaliz for its H-result. The stabilizer
        # tuple itself is cached by cell_stabilizer and deliberately remains
        # available for the traversal's later categorical stabilizer query.
        if generators:
            del generator, action
        del positions

        def facet_keys() -> Iterator[FacetMaskKey]:
            covector_rows = iter(
                self._iter_cell_halfspace_pairing_rows(
                    cell, labels, use_owned_cone=use_owned_cone
                )
            )
            block_size = _configuration_pairing_block_size(ray_count)
            sparse_limit = _sparse_facet_complement_limit(ray_count)
            # NmzResult("SupportHyperplanes") materializes Normaliz's complete
            # exact H-result before its first row can be consumed. Delay the
            # dense transposed ray matrix until that result exists, so it does
            # not add another rank-by-ray allocation to Normaliz's peak while
            # computing the full support-hyperplane list.
            covector_block = tuple(islice(covector_rows, block_size))
            if not covector_block:
                return
            transposed_ray_rows = matrix(
                SageZZ, rank, ray_count, lambda i, j: cell.rows[j][i]
            )
            while covector_block:
                pairing_block = matrix(SageZZ, covector_block) * transposed_ray_rows
                block_keys = []
                for offset in range(len(covector_block)):
                    pairing_row = pairing_block.row(offset)
                    nonincident_positions: list[int] = []
                    incident_positions: list[int] = []
                    nonincident_count = 0
                    incident_count = 0
                    for position, value in enumerate(pairing_row):
                        if value:
                            nonincident_count += 1
                            if nonincident_count <= sparse_limit:
                                nonincident_positions.append(position)
                            elif nonincident_count == sparse_limit + 1:
                                nonincident_positions.clear()
                        else:
                            incident_count += 1
                            if incident_count <= sparse_limit:
                                incident_positions.append(position)
                            elif incident_count == sparse_limit + 1:
                                incident_positions.clear()
                    if nonincident_count == 0:
                        raise ArithmeticError("a minimized H-presentation covector contains every extreme ray")
                    if incident_count == 0:
                        raise ArithmeticError("a minimized H-presentation covector has no incident extreme rays")
                    if nonincident_count <= sparse_limit:
                        key: FacetMaskKey = (True, tuple(nonincident_positions))
                    elif incident_count <= sparse_limit:
                        key = (False, tuple(incident_positions))
                    else:
                        # Store the smaller side as the bitset. Applying a
                        # ray permutation visits each set bit, so this choice
                        # reduces the orbit cost for facets incident to only a
                        # small fraction of the 49,861 input rays.
                        store_nonincident = nonincident_count <= incident_count
                        support_bytes = bytearray((ray_count + 7) // 8)
                        for position, value in enumerate(pairing_row):
                            if bool(value) == store_nonincident:
                                support_bytes[position >> 3] |= 1 << (position & 7)
                        support_mask = int.from_bytes(support_bytes, "little")
                        key = support_mask if store_nonincident else (False, support_mask)
                    block_keys.append(key)
                # Detach the block's compact keys before yielding any of them.
                # Orbit traversal may do substantial work before asking for
                # the next key; retaining the dense pairing matrix across that
                # work overlaps its memory with the current Schreier orbit.
                del pairing_block, pairing_row, nonincident_positions, incident_positions
                # The primitive coordinate tuples have served their purpose
                # once the keys are built; do not retain them while the
                # consumer traverses an orbit for each yielded key.
                del covector_block
                yield from block_keys
                del block_keys
                covector_block = tuple(islice(covector_rows, block_size))

        if not use_owned_cone:
            # Normaliz's Python API returns its complete support-hyperplane
            # family as a list. If that list remains inside the facet-key
            # generator while the quotient traversal runs, all not-yet-read
            # H-rows overlap with the growing visited-facet set and Schreier
            # orbit queue. Drain it now into exact incidence keys, then let
            # the consumer mark values in this same dictionary as it visits
            # them. This keeps the full H-result out of the long traversal
            # phase without materializing the dense pairing graph. Canonical
            # key objects are shared with the traversal state: dense integer
            # bitsets can be several kilobytes each, so retaining a second
            # equal object for every permuted facet is prohibitive.
            canonical_keys: dict[FacetMaskKey, FacetMaskKey | None] = {}
            for key in facet_keys():
                canonical_keys.setdefault(key, key)

            def consume_keys() -> Iterator[FacetMaskKey]:
                # Dict insertion order is Normaliz H-presentation order, so
                # the canonical key table itself is the ordered key stream.
                # Do not retain a parallel deque of one pointer per facet.
                yield from canonical_keys

            return consume_keys(), tuple(permutations), canonical_keys
        return facet_keys(), tuple(permutations), None

    @staticmethod
    def _iter_facet_mask_orbits(
        facet_keys: Iterator[FacetMaskKey],
        permutations: tuple[array, ...],
        *,
        retain_words: bool = True,
        canonical_keys: dict[FacetMaskKey, FacetMaskKey | None] | None = None,
    ) -> Iterator[list[tuple[FacetMaskKey, SchreierWord] | None]]:
        """Yield exact facet orbits with one visited set and one current orbit.

        Every generator permutation comes from an exact lattice isometry
        preserving this ray configuration, hence its image of a facet is a
        facet. Full-domain traversal therefore needs only a visited set. When
        the raw-H path supplies a canonical-key dictionary, that same mapping
        doubles as visitation state rather than allocating a second table.
        """
        identity_word = None
        missing_key = object()
        # The raw-H path already owns a canonical-key dictionary. Use its
        # values as visitation marks instead of allocating a second full set
        # of references to every (potentially multi-kilobyte) facet key.
        visited: set[FacetMaskKey] | None = set() if canonical_keys is None else None
        facet_count = 0
        for start in facet_keys:
            facet_count += 1
            if canonical_keys is None:
                assert visited is not None
                if start in visited:
                    continue
                visited.add(start)
            else:
                canonical_start = canonical_keys[start]
                if canonical_start is None:
                    continue
                canonical_keys[start] = None
                start = canonical_start
            if not retain_words:
                # Group-generator assembly only needs one facet from each
                # orbit. Yield it before exploring the rest of the orbit, so
                # the caller can flip the representative immediately. Then
                # finish the exact orbit using a queue of keys alone: retaining
                # a Schreier word per facet is needed by the full-domain API,
                # but redundant for O(L) generators.
                pending = deque([start])
                yield [(start, identity_word)]
                while pending:
                    facet_mask = pending.popleft()
                    for permutation in permutations:
                        image_mask = _permuted_facet_key(facet_mask, permutation)
                        if canonical_keys is not None:
                            canonical_image = canonical_keys.get(image_mask, missing_key)
                            if canonical_image is missing_key:
                                raise ArithmeticError(
                                    "a cell-stabilizer generator sends a facet outside the exact H-presentation"
                                )
                            if canonical_image is None:
                                continue
                            image_mask = canonical_image
                            canonical_keys[image_mask] = None
                            pending.append(image_mask)
                        else:
                            assert visited is not None
                            if image_mask in visited:
                                continue
                            visited.add(image_mask)
                            pending.append(image_mask)
                continue
            orbit: list[tuple[FacetMaskKey, SchreierWord] | None] = [(start, identity_word)]
            for entry in orbit:
                if entry is None:
                    raise ArithmeticError("a facet orbit buffer was cleared before traversal completed")
                facet_mask, element_word = entry
                for generator_index, permutation in enumerate(permutations):
                    image_mask = _permuted_facet_key(facet_mask, permutation)
                    if canonical_keys is not None:
                        canonical_image = canonical_keys.get(image_mask, missing_key)
                        if canonical_image is missing_key:
                            raise ArithmeticError(
                                "a cell-stabilizer generator sends a facet outside the exact H-presentation"
                            )
                        if canonical_image is None:
                            continue
                        image_mask = canonical_image
                        canonical_keys[image_mask] = None
                        orbit.append((image_mask, (element_word, generator_index)))
                    else:
                        assert visited is not None
                        if image_mask in visited:
                            continue
                        visited.add(image_mask)
                        orbit.append((image_mask, (element_word, generator_index)))
            # Transfer the mutable BFS buffer to the caller rather than
            # copying every orbit into a tuple while the list is still live.
            yield orbit
        if canonical_keys is None:
            assert visited is not None
            visited_count = len(visited)
        else:
            visited_count = sum(value is None for value in canonical_keys.values())
        if visited_count != facet_count:
            raise ArithmeticError(
                "cell-stabilizer permutations generate incidences outside the H-presentation"
            )

    def _facet_orbit_transversal_masks(
        self, cell: LorentzianPerfectCell
    ) -> tuple[Iterator[list[tuple[FacetMaskKey, SchreierWord] | None]], tuple[array, ...]]:
        """Build compact facet-orbit buffers for one consuming operation.

        Deliberately do not cache this result: public consumers materialize
        full incidences, ray permutations, or owned isometries from these masks.
        Retaining the masks in a cached-method table would keep a second
        complete facet-orbit representation alive after that conversion, which
        is prohibitive for the 49,861-ray A7 configuration.
        """
        facet_keys, permutations, _canonical_keys = self._facet_mask_data(cell)
        # The iterator transfers ownership of each completed BFS list. Turning
        # each list into a tuple here duplicates its pointer array at the peak
        # where the full facet collection is retained, with no consumer needing
        # tuple-specific operations on these private orbit buffers.
        return self._iter_facet_mask_orbits(facet_keys, permutations), permutations

    def _facet_orbit_transversal_words(
        self, cell: LorentzianPerfectCell
    ) -> tuple[
        tuple[tuple[tuple[FacetIncidence, SchreierWord], ...], ...],
        tuple[array, ...],
    ]:
        """Materialize public-style incidences while releasing consumed mask entries."""
        mask_orbits, permutations = self._facet_orbit_transversal_masks(cell)
        def materialize_orbit(orbit):
            def consume_orbit():
                for index, entry in enumerate(orbit):
                    if entry is None:
                        raise ArithmeticError("a facet orbit entry was consumed more than once")
                    mask, word = entry
                    yield (_facet_incidence(mask, len(cell.rows)), word)
                    orbit[index] = None

            return tuple(consume_orbit())

        return tuple(materialize_orbit(orbit) for orbit in mask_orbits), permutations

    def facet_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[tuple[FacetIncidence, LatticeIsometryMethods], ...], ...]:
        """Facet orbits with owned isometries, materialized for this public view."""
        labels = tuple(cell.lattice.module_generating_set())
        ray_rows = matrix(SageQQ, cell.rows)
        basis_positions = tuple(int(index) for index in ray_rows.transpose().pivots())
        if len(basis_positions) != len(labels):
            raise ArithmeticError("a perfect-cell configuration does not span its ambient lattice")
        inverse_basis = ray_rows.matrix_from_rows(basis_positions).inverse()
        mask_orbits, generator_permutations = self._facet_orbit_transversal_masks(cell)
        def materialize_orbit(orbit):
            def materialize_element(item):
                nonincident_mask, word = item
                # Apply the Schreier word only to the rank-sized basis rows.
                # Materializing a full ray permutation for every facet costs
                # O(number_of_facets * number_of_rays) integers; this public
                # view needs only the basis images to recover each isometry.
                image_positions = list(basis_positions)
                for generator_index in _schreier_generator_indices(word):
                    permutation = generator_permutations[generator_index]
                    image_positions = [permutation[position] for position in image_positions]
                image_basis = matrix(SageQQ, [ray_rows[position] for position in image_positions])
                action = inverse_basis * image_basis
                if any(entry.denominator() != 1 for entry in action.list()):
                    raise ArithmeticError("a Schreier facet permutation does not induce an integral lattice action")
                isometry = cell.lattice.O()(tuple(cell.lattice([int(entry) for entry in row]) for row in action.rows()))
                facet = _facet_incidence(nonincident_mask, len(cell.rows))
                return (facet, isometry)

            def consume_orbit():
                for index, item in enumerate(orbit):
                    if item is None:
                        raise ArithmeticError("a facet orbit entry was consumed more than once")
                    yield materialize_element(item)
                    orbit[index] = None

            return tuple(consume_orbit())

        return tuple(materialize_orbit(orbit) for orbit in mask_orbits)

    def flip_across(self, cell: LorentzianPerfectCell, facet: FacetIncidence) -> LorentzianPerfectCell | None:
        normalized = _lorentzian_view(cell.lattice)
        vectors = cell.rows
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
        gram = _gram_matrix(normalized)
        if int(cell.lattice.module_rank()) == 2 and all(vector(SageZZ, vectors[position]) * gram * vector(SageZZ, vectors[position]) == 0 for position in facet):
            return None
        selected = matrix(SageQQ, [vectors[position] for position in sorted(facet)])
        facet_nullspace = selected.right_kernel_matrix()
        if facet_nullspace.nrows() != 1:
            raise ArithmeticError("a perfect-cell facet must have one normal direction")
        direction = facet_nullspace.row(0)
        outside = vector(SageQQ, vectors[nonincident])
        if direction.dot_product(outside) < 0:
            direction = -direction
        normal_direction = _source_remove_fraction(vector(SageQQ, [0, *direction]))
        expanded = matrix(SageQQ, [[1, *row] for row in vectors])
        base_nullspace = expanded.right_kernel_matrix()
        if base_nullspace.nrows() != 1:
            raise ArithmeticError("a perfect cell must have a unique affine supporting functional")
        base_normal = base_nullspace.row(0)
        if vector(SageQQ, base_normal[1:]).dot_product(outside) <= 0:
            base_normal = -base_normal
        base_normal = _source_remove_fraction(base_normal)
        critical = tuple(vectors[position] for position in sorted(facet))
        flipped, _normal, _direction, _max_scal = _kernel_flipping(
            normalized, critical, base_normal, normal_direction, cell.allowed_squares
        )
        return LorentzianPerfectCell(cell.lattice, flipped, cell.allowed_squares)

    def _integral_configuration_invariants(self, cell: LorentzianPerfectCell):
        r"""Exact pairing data of a configuration, computed once per cell.

        The quotient traversal compares each representative with many
        neighboring cells. The pairings and each vertex's pairing spectrum
        depend only on the configuration; the pairings are read from the
        n x rank rows on request, never stored as an n x n matrix. Cache these
        exact witnesses without replacing the integral-isometry decision by a
        hash. The two-entry identity LRU retains the source and current target
        while avoiding an unbounded cache of dense matrices for every neighbor.
        """
        cache = self.__dict__.setdefault("_configuration_invariant_cache", {})
        key = id(cell)
        cached = cache.pop(key, None)
        if cached is not None and cached[0] is cell:
            cache[key] = cached
            return cached[1]
        rows = matrix(SageZZ, cell.rows)
        pairings = _ConfigurationPairings(rows, _gram_matrix(cell.lattice), cache_columns=True)
        if len(cell.rows) > _PAIRING_GRAPH_REFINEMENT_LIMIT:
            # The complete pairing-spectrum colour requires n(n-1)/2
            # pairings. Coordinate content and Gram-covector divisibility are
            # inexpensive exact integral-isometry invariants: on the frozen
            # 49,861-row A7 configuration, their two 10-coordinate gcd passes
            # take under half a second and the covector divisibility separates
            # the nonprimitive pairings. Retain these colours rather than
            # collapsing the entire search to norms alone.
            spectra_list = []
            index = 0
            for norm_block in pairings.diagonal_blocks():
                for norm in norm_block:
                    row = cell.rows[index]
                    covector = pairings.covectors.row(index)
                    spectra_list.append(
                        (int(norm), gcd(*map(int, row)), gcd(*map(int, covector)))
                    )
                    index += 1
            if index != len(cell.rows):
                raise ArithmeticError("the diagonal pairing stream has the wrong row count")
            spectra = tuple(spectra_list)
            del spectra_list
        else:
            norms = pairings.diagonal()
            spectra = []
            for start, block in pairings.blocks():
                for offset, row in enumerate(block):
                    index = start + offset
                    spectra.append((int(norms[index]), gcd(*cell.rows[index]), _row_multiset(row, int(row[index]))))
            spectra = tuple(spectra)
        # The nonzero Smith factors are invariant under right multiplication
        # by an integral unimodular action. Taking the transpose keeps the
        # invariant list bounded by the ambient rank rather than the number
        # of configuration rows (which is 49,861 in the A7 consumer).
        smith_invariants = rows.transpose().elementary_divisors()
        row_sum = vector(SageZZ, (sum(row[j] for row in cell.rows) for j in range(rows.ncols())))
        row_moment = rows.transpose() * rows
        distinct_rows = len(set(cell.rows)) == len(cell.rows)
        data = (rows, pairings, spectra, smith_invariants, row_sum, row_moment, distinct_rows)
        cache[key] = (cell, data)
        while len(cache) > 2:
            cache.pop(next(iter(cache)))
        return data

    @overload
    def cell_transporter(
        self,
        source: LorentzianPerfectCell,
        target: LorentzianPerfectCell,
        *,
        all_images: Literal[True],
        generators_only: bool = False,
    ) -> tuple[LatticeIsometryMethods, ...]: ...

    @overload
    def cell_transporter(
        self,
        source: LorentzianPerfectCell,
        target: LorentzianPerfectCell,
        *,
        all_images: Literal[False] = False,
        generators_only: Literal[False] = False,
    ) -> LatticeIsometryMethods | None: ...

    def cell_transporter(
        self,
        source: LorentzianPerfectCell,
        target: LorentzianPerfectCell,
        *,
        all_images: bool = False,
        generators_only: bool = False,
    ) -> LatticeIsometryMethods | None | tuple[LatticeIsometryMethods, ...]:
        r"""An exact integral isometry transporting one entire perfect-cell configuration.

        Moderate configurations use the owned finite-family transporter: a
        labelled pairing-graph isomorphism followed by integral-structure
        transport, as in Lorentzian.g:660-696 and
        AutomorphismPolytope.g:288-339. Large configurations avoid that
        quadratic graph and use the exact basis-image search, stopping at the
        first fully certified witness when only one transporter is requested.

        For enumeration, a spanning subfamily of rank r determines an isometry uniquely from
        its r images. Enumerate all norm- and pairing-compatible target images
        of a fixed source basis; at each leaf verify
        integrality, unimodularity, the entire Gram form, and the exact
        multiset of configuration vectors. No partial graph isomorphism is
        treated as a certificate, and no potential integral transporter is
        excluded by these necessary conditions.
        """
        if source.allowed_squares != target.allowed_squares:
            raise ValueError("an isometry of perfect cells must preserve their prescribed normalized-norm locus")
        if generators_only and (not all_images or source is not target):
            raise ValueError("stabilizer-generator search requires one identical source and target perfect cell")
        n = len(source.rows)
        if n != len(target.rows):
            return () if all_images else None
        rank = int(source.lattice.module_rank())
        if rank != int(target.lattice.module_rank()):
            return () if all_images else None
        # An integral transporter is unimodular.  Its congruence of the
        # ambient Gram matrices preserves determinant and inertia, even
        # when the configurations have identical pairing spectra.  Reject
        # nonisometric ambient forms before constructing the potentially
        # expensive coloured-configuration search.
        source_gram = _gram_matrix(source.lattice)
        target_gram = _gram_matrix(target.lattice)
        if source_gram.det() != target_gram.det() or source.lattice.signature_pair() != target.lattice.signature_pair():
            return () if all_images else None
        if not all_images:
            if source.lattice is target.lattice and source.rows == target.rows:
                return source.lattice.Aut().identity()
            if n <= _PAIRING_GRAPH_REFINEMENT_LIMIT:
                if source.lattice is target.lattice and Counter(source.rows) == Counter(target.rows):
                    return source.lattice.Aut().identity()
                source_family = Sets().Mor(Sets.Δ[n - 1], source.lattice)(
                    lambda index: source.vector_configuration[int(index)]
                )
                target_family = Sets().Mor(Sets.Δ[n - 1], target.lattice)(
                    lambda index: target.vector_configuration[int(index)]
                )
                return source.lattice.Isom(target.lattice).transporter(source_family, target_family)
            # The preamble's labelled finite-family transporter constructs a
            # complete pairing graph. For large cells, use the same exact
            # basis-image search as stabilizer enumeration but stop at the
            # first fully certified witness instead of building that graph.
        source_data = self._integral_configuration_invariants(source)
        # Stabilizer generation is the hot path for perfect-cell automorphisms:
        # source and target are literally the same configuration. Rebuilding
        # its row matrix, Gram-covector matrix, invariant colours and Smith data
        # twice only duplicates O(n * rank) exact storage and work.
        target_data = source_data if source is target else self._integral_configuration_invariants(target)
        (
            source_rows,
            source_pairings,
            source_spectra,
            source_smith,
            source_sum,
            source_moment,
            source_distinct_rows,
        ) = source_data
        (
            target_rows,
            target_pairings,
            target_spectra,
            target_smith,
            target_sum,
            target_moment,
            target_distinct_rows,
        ) = target_data
        # For large configurations, certify equality of the row multisets by
        # exact target-set membership when both inputs have distinct rows. An
        # integral unimodular action is injective, so n distinct source rows
        # landing in the n-element target set prove equality. This avoids
        # building a second large Counter for every neighbor transporter.
        # Preserve the multiset fallback for repeated-row configurations.
        use_position_certificate = generators_only or (
            n > _PAIRING_GRAPH_REFINEMENT_LIMIT
            and source_distinct_rows
            and target_distinct_rows
        )
        target_multiset = None if use_position_certificate else Counter(target.rows)
        if source_smith != target_smith:
            return ()
        # The second moment transforms by M_target = A^T M_source A for an
        # integral unimodular transporter A, so its determinant is an exact
        # basis-independent invariant. It is only rank-by-rank and rejects
        # some incompatible configurations before the large image search.
        if source_moment.det() != target_moment.det():
            return ()
        # The row sum itself transforms by the same unimodular action. Its
        # coordinate content and its norm in the ambient Gram form are exact
        # invariants too, and cost only O(rank) arithmetic before backtracking.
        if gcd(*map(int, source_sum)) != gcd(*map(int, target_sum)):
            return ()
        if source_sum.dot_product(source_gram * source_sum) != target_sum.dot_product(
            target_gram * target_sum
        ):
            return ()
        if source is not target and Counter(source_spectra) != Counter(target_spectra):
            return ()
        signatures = set(source_spectra)
        if source is not target:
            signatures.update(target_spectra)
        colour_numbers = {signature: index for index, signature in enumerate(sorted(signatures))}
        source_colours = tuple(colour_numbers[signature] for signature in source_spectra)
        target_colours = (
            source_colours
            if source is target
            else tuple(colour_numbers[signature] for signature in target_spectra)
        )

        # For moderate configurations, refine the two coloured,
        # integral-weighted configuration graphs (the vertex-coloured Gram
        # graph of LinPolytopeIntegral_Isomorphism_GramMat) by complete
        # neighbour-colour distributions. Each refinement signature is an
        # integral-isometry invariant; the partition only splits and reaches
        # a fixed point. For larger configurations, the intrinsic colours
        # above are retained without constructing the quadratic pairing graph.
        if n <= _PAIRING_GRAPH_REFINEMENT_LIMIT:
            while True:
                modulus = len(set(source_colours) | set(target_colours))
                source_refined = _refined_signatures(source_pairings, source_colours, modulus)
                target_refined = _refined_signatures(target_pairings, target_colours, modulus)
                refined_numbers = {signature: index for index, signature in enumerate(sorted(set(source_refined) | set(target_refined)))}
                refined_source_colours = tuple(refined_numbers[signature] for signature in source_refined)
                refined_target_colours = tuple(refined_numbers[signature] for signature in target_refined)
                if Counter(refined_source_colours) != Counter(refined_target_colours):
                    return ()
                # Refinement signatures contain the old colour, so no two old
                # classes can merge. Once the number of classes stops growing,
                # the partition is equitable and further refinement cannot
                # split it. Comparing numerical colour IDs instead is wrong:
                # canonical renumbering can permute IDs even when the partition
                # has stabilized, potentially causing a renumbering cycle.
                stable = len(refined_numbers) == len(set(source_colours) | set(target_colours))
                source_colours, target_colours = refined_source_colours, refined_target_colours
                if stable:
                    break

        buckets = {}
        for j, colour in enumerate(target_colours):
            buckets.setdefault(colour, []).append(j)
        # Build each large colour fibre's hash table once. Re-hashing the
        # 22k-member A7 isotropic fibre for every unassigned basis row and
        # every recursion node is avoidable; copy the prebuilt set and remove
        # already-used images instead.
        bucket_sets = {colour: set(indices) for colour, indices in buckets.items()}
        # A partial isometry i -> j constrains every further basis image
        # k -> l by the exact Gram pairing.
        # Only images of the selected basis occur as anchors. Index their
        # possible target neighbours by the exact pairing signature, and
        # materialize a compatibility fibre only when backtracking needs it.
        # Eagerly storing all (i,j,k) fibres costs cubic space in the
        # configuration size even when the search rejects its first branch.
        compatible_images = {}
        target_neighbour_signatures = {}
        # Each target anchor has the same exact pairing signature
        # partition regardless of the source vertex assigned to it.
        # Build its immutable fibres once; backtracking only intersects them.
        # Each cached fibre/partition can contain O(n) indices. Bound both
        # caches by total entry budget rather than a fixed number of anchors,
        # otherwise the 49,861-vector A7 search can retain quadratic data.
        anchor_cache_limit = max(1, min(64, 250_000 // max(1, n)))

        if source_rows.rank() != rank or target_rows.rank() != rank:
            raise ArithmeticError("a perfect-cell vector configuration does not span its lattice")
        # A rare spectrum minimizes the initial image search. Select the
        # lexicographically first independent rows in that rarity ordering
        # from one exact elimination, rather than recomputing the rank of
        # every candidate prefix at every basis position.
        rarity_order = tuple(sorted(range(n), key=lambda i: (len(buckets[source_colours[i]]), i)))
        ordered_rows = source_rows.matrix_from_rows(rarity_order)
        basis = tuple(rarity_order[int(position)] for position in ordered_rows.transpose().pivots())
        if len(basis) != rank:
            raise ArithmeticError("failed to select a basis of a spanning perfect-cell configuration")
        source_basis = source_rows.matrix_from_rows(basis)
        inverse_source_basis = source_basis.change_ring(SageQQ).inverse()
        source_sum_in_basis = vector(SageQQ, source_sum) * inverse_source_basis
        target_sum_rational = vector(SageQQ, target_sum)

        # Any transporter carries the ordered source basis to an ordered
        # target basis with the same lattice index. This determinant
        # obstruction applies before constructing the candidate rational
        # action, and is stronger than testing the Smith invariants of
        # proper subsets alone.
        source_basis_index = abs(source_basis.det())

        source_labels = tuple(source.lattice.module_generating_set())
        target_labels = tuple(target.lattice.module_generating_set())
        # Exact aggregate invariants of the complete row multiset. Every
        # configuration transporter permutes rows, so it must carry their
        # vector sum and second-moment matrix to the corresponding target
        # aggregates. These rank-sized certificates cheaply reject a basis
        # assignment before the O(n * rank) row-membership pass.
        source_submodule_invariants = {}
        target_submodule_invariants = {}
        # Source subsets come from a fixed rank-sized basis, but target
        # subsets range over the entire 49,861-row configuration. Bound both
        # caches: the target search can visit arbitrarily many distinct
        # subsets before finding a witness, and an unbounded dictionary would
        # retain every failed branch's tuple keys and Smith factors.
        submodule_cache_limit = max(64, min(8_192, 250_000 // max(1, rank)))
        # Smith data depend on the selected row subset, not its order.
        # Cache proper partial spans; full-rank leaves instead receive
        # the complete integral transporter certificate.
        stabilizer_elements: list[LatticeIsometryMethods] = []
        first_witness: LatticeIsometryMethods | None = None
        first_only = not all_images
        if use_position_certificate:
            # The immutable coordinate tuples are already the configuration's
            # canonical rows; reuse them as keys instead of allocating another
            # set of row tuples. The stabilizer uses source positions; ordinary
            # large-cell transporters use target positions.
            position_rows = source.rows if generators_only else target.rows
            positions = {row: index for index, row in enumerate(position_rows)}
            if len(positions) != n:
                raise ArithmeticError("the exact row-position certificate requires distinct configuration vectors")
        if generators_only:
            symmetric = SymmetricGroup(n)
            subgroup = symmetric.subgroup([])
            group_permutations = []

        def forced_sum_image(images):
            """Use the exact configuration sum to force a last basis image.

            In source-basis coordinates, ``source_sum = alpha * B``. If all
            but one nonzero coefficient of alpha have assigned target images,
            the target sum determines the remaining image row uniquely. This
            is an exact necessary condition, not a heuristic; it often avoids
            trying thousands of candidates for the final basis vector.
            """
            if not use_position_certificate:
                return None
            unassigned = []
            assigned = []
            for position, coefficient in enumerate(source_sum_in_basis):
                if not coefficient:
                    continue
                source_index = basis[position]
                target_index = images.get(source_index)
                if target_index is None:
                    unassigned.append((source_index, coefficient))
                else:
                    assigned.append((coefficient, target_index))
            if len(unassigned) > 1:
                return None
            assigned_sum = vector(SageQQ, rank)
            for coefficient, target_index in assigned:
                assigned_sum += coefficient * vector(SageQQ, target_rows.row(target_index))
            if not unassigned:
                return True if assigned_sum == target_sum_rational else False
            source_index, coefficient = unassigned[0]
            required = (target_sum_rational - assigned_sum) / coefficient
            if any(entry.denominator() != 1 for entry in required):
                return False
            target_index = positions.get(tuple(int(entry) for entry in required))
            if target_index is None or target_index in images.values():
                return False
            return source_index, target_index

        def backtrack(images, source_columns, target_columns, certified_target_basis=None):
            nonlocal subgroup, first_witness
            if first_witness is not None:
                return first_witness
            if len(images) == rank:
                target_basis_integral = certified_target_basis if certified_target_basis is not None else target_rows.matrix_from_rows([images[i] for i in basis])
                if certified_target_basis is None and abs(target_basis_integral.det()) != source_basis_index:
                    return None
                target_basis = target_basis_integral.change_ring(SageQQ)
                action = inverse_source_basis * target_basis
                if any(entry.denominator() != 1 for entry in action.list()):
                    return None
                # det(action) = det(target_basis)/det(source_basis), so the
                # preceding equality of absolute determinants already proves
                # unimodularity once the rational action is integral.
                action = action.change_ring(SageZZ)
                if action * target_gram * action.transpose() != source_gram:
                    return None
                if source_sum * action != target_sum:
                    return None
                if action.transpose() * source_moment * action != target_moment:
                    return None
                if use_position_certificate:
                    # The exact row-position map certifies equality of the
                    # full row multiset without allocating another Counter or
                    # a second full n-by-r matrix at every completed basis
                    # assignment. Keep transformed rows bounded.
                    permutation_images = [] if generators_only else None
                    block_size = _configuration_pairing_block_size(n)
                    for start in range(0, n, block_size):
                        stop = min(start + block_size, n)
                        transformed_block = source_rows.matrix_from_rows(range(start, stop)) * action
                        for row in transformed_block.rows():
                            position = positions.get(tuple(int(entry) for entry in row))
                            if position is None:
                                return None
                            if permutation_images is not None:
                                permutation_images.append(position + 1)
                    if permutation_images is not None:
                        permutation = symmetric(permutation_images)
                else:
                    transformed_rows = source_rows * action
                    if target_multiset is None or Counter(tuple(row) for row in transformed_rows.rows()) != target_multiset:
                        return None
                witness = source.lattice.Isom(target.lattice)(
                    {
                        label: target.lattice.linear_combination(
                            {target_label: target.lattice.base_ring()(int(entry)) for target_label, entry in zip(target_labels, row, strict=True) if entry}
                        )
                        for label, row in zip(source_labels, action.rows(), strict=True)
                    }
                )
                if generators_only:
                    if permutation not in subgroup:
                        stabilizer_elements.append(witness)
                        group_permutations.append(permutation)
                        subgroup = symmetric.subgroup(group_permutations)
                elif first_only:
                    first_witness = witness
                    return witness
                else:
                    stabilizer_elements.append(witness)
                return None

            forced = forced_sum_image(images)
            if forced is False:
                return None
            forced_pair = forced if isinstance(forced, tuple) else None
            used = set(images.values())
            choice = None
            candidates_for_choice = None
            basis_order = basis
            forced_source = forced_target = -1
            if forced_pair is not None:
                forced_source, forced_target = forced_pair
                if target_colours[forced_target] != source_colours[forced_source]:
                    return None
                # The row-sum equation has already proved this is the only
                # possible image for forced_source. Choose it immediately;
                # constructing the other large colour fibres first would
                # discard the exact information just computed.
                basis_order = (forced_source, *(i for i in basis if i != forced_source))
            for i in basis_order:
                if i in images:
                    continue
                if forced_pair is not None and i == forced_source:
                    choice = i
                    candidates_for_choice = [forced_target]
                    break
                admissible_images = None
                for k, v in images.items():
                    key = (k, v, i)
                    if key not in compatible_images:
                        # The target neighbour partition depends on the
                        # image v alone, not on which source anchor k maps
                        # to v. Reuse it across distinct partial assignments.
                        # Packed unsigned arrays replace one Python set of
                        # boxed row indices per pairing fibre; on the frozen
                        # A7 configuration there are 49,861 indices per
                        # anchor, so frozenset-of-int storage dominates the
                        # actual exact pairing values by a large factor.
                        if v not in target_neighbour_signatures:
                            index_typecode = "I" if n <= 2**32 else "Q"
                            target_neighbour_signatures[v] = _packed_pairing_fibres(
                                target_colours,
                                target_pairings.column(v),
                                index_typecode,
                            )
                            if len(target_neighbour_signatures) > anchor_cache_limit:
                                target_neighbour_signatures.pop(next(iter(target_neighbour_signatures)))
                        signature = (source_colours[i], source_pairings[i, k])
                        compatible_images[key] = target_neighbour_signatures[v].get(signature, frozenset())
                        if len(compatible_images) > anchor_cache_limit:
                            compatible_images.pop(next(iter(compatible_images)))
                    fibre = compatible_images[key]
                    if admissible_images is None:
                        # A pairing fibre is often far smaller than the raw
                        # colour class. Start from that packed fibre rather
                        # than allocating a set for all 22k A7 isotropic rays.
                        if len(fibre) < len(bucket_sets[source_colours[i]]):
                            admissible_images = set(fibre)
                        else:
                            admissible_images = bucket_sets[source_colours[i]].copy()
                            admissible_images.intersection_update(fibre)
                    else:
                        admissible_images.intersection_update(fibre)
                    if not admissible_images:
                        break
                if admissible_images is None:
                    admissible_images = bucket_sets[source_colours[i]].copy()
                admissible_images.difference_update(used)
                candidates = sorted(admissible_images)
                if not candidates:
                    return None
                if candidates_for_choice is None or len(candidates) < len(candidates_for_choice):
                    choice = i
                    candidates_for_choice = candidates
            if choice is None or candidates_for_choice is None:
                raise ArithmeticError("basis-image search lost an unassigned basis vector")
            # Let H be the subgroup already generated by exact integral
            # witnesses. For a partial assignment a, its pointwise
            # stabilizer H_a acts on the admissible next images. Composing
            # a complete transporter with an element of H_a does not alter
            # any previously assigned basis image. Therefore just one
            # representative of each H_a-orbit of next images is needed
            # to generate the full stabilizer together with H.
            #
            # H grows when a new witness is found. Recompute the orbit
            # representatives when its generator count changes, so the
            # search can reduce previously pending sibling branches too.
            version = -1
            admissible = set(candidates_for_choice)
            next_source_columns = None
            source_signature_counts = None
            for j in candidates_for_choice:
                if generators_only:
                    # A witness found in an earlier sibling can enlarge H.
                    # The current sibling list was fixed before that
                    # discovery, so recompute the H_a-orbit partition.
                    if version != len(group_permutations):
                        fixer = subgroup
                        for assigned_target in images.values():
                            fixer = fixer.stabilizer(assigned_target + 1)
                        orbit_minima = {}
                        remaining_orbits = set(admissible)
                        while remaining_orbits:
                            candidate = min(remaining_orbits)
                            orbit = {int(position) - 1 for position in fixer.orbit(candidate + 1)}
                            if not orbit <= admissible:
                                raise ArithmeticError("an integral pointwise stabilizer does not preserve the admissible image class")
                            for member in orbit:
                                orbit_minima[member] = candidate
                            remaining_orbits.difference_update(orbit)
                        version = len(group_permutations)
                    if orbit_minima[j] != j:
                        continue
                images[choice] = j
                # This is a necessary condition for any ambient isometry,
                # including at intermediate nodes of the basis search.
                if any(target_pairings[j, images[k]] != source_pairings[choice, k] for k in images if k != choice):
                    del images[choice]
                    continue
                # The Smith invariants of every partial span are preserved
                # by a unimodular ambient transformation. This is an exact
                # integral obstruction, not a cutoff on the finite search.
                selected = tuple(sorted(images))
                # At full rank the absolute determinant is the lattice
                # index of the selected image basis. Reject a different
                # index before constructing its Smith normal form or the
                # rational transporter. This is precisely the determinant
                # obstruction checked again at the completed leaf.
                if len(selected) == rank:
                    # The leaf performs the stronger complete certificate:
                    # integrality, unimodularity, Gram preservation and
                    # equality of the entire vector multiset. Smith forms
                    # at full rank are redundant after the determinant
                    # check, and can dominate the cost of each leaf.
                    target_basis_integral = target_rows.matrix_from_rows([images[i] for i in basis])
                    if abs(target_basis_integral.det()) == source_basis_index:
                        # Extend signatures in basis-assignment order. Source
                        # and target use the same sequence of assigned pairs,
                        # so no sorting by source index is needed for this
                        # exact multiset comparison.
                        if next_source_columns is None:
                            source_column = source_pairings.column(choice)
                            next_source_columns = source_columns + (source_column,)
                            source_signature_counts = Counter(
                                (colour, tuple(column[index] for column in next_source_columns))
                                for index, colour in enumerate(source_colours)
                            )
                        assert next_source_columns is not None and source_signature_counts is not None
                        target_column = target_pairings.column(j)
                        next_target_columns = target_columns + (target_column,)
                        if _configuration_signature_multisets_match(
                            source_colours,
                            next_source_columns,
                            target_colours,
                            next_target_columns,
                            source_signature_counts=source_signature_counts,
                        ):
                            # The parent frame otherwise retains its
                            # 49,861-entry source signature counter throughout
                            # the recursive search. Rebuild it only if another
                            # sibling survives the exact filters.
                            child_source_columns = next_source_columns
                            source_signature_counts = None
                            next_source_columns = None
                            backtrack(images, child_source_columns, next_target_columns, target_basis_integral)
                            if first_witness is not None:
                                return first_witness
                    del images[choice]
                    continue
                # The chosen source rows are linearly independent. Every
                # rational isometry preserves their rank, so dependent
                # target rows cannot extend to a transporter. Check this
                # exact obstruction before computing Smith invariants.
                target_selected = tuple(sorted(images[i] for i in selected))
                target_partial = target_rows.matrix_from_rows(target_selected)
                if target_partial.rank() != len(selected):
                    del images[choice]
                    continue
                if selected not in source_submodule_invariants:
                    source_submodule_invariants[selected] = tuple(source_rows.matrix_from_rows(selected).elementary_divisors())
                    if len(source_submodule_invariants) > submodule_cache_limit:
                        source_submodule_invariants.pop(next(iter(source_submodule_invariants)))
                if target_selected not in target_submodule_invariants:
                    target_submodule_invariants[target_selected] = tuple(target_partial.elementary_divisors())
                    if len(target_submodule_invariants) > submodule_cache_limit:
                        target_submodule_invariants.pop(next(iter(target_submodule_invariants)))
                if source_submodule_invariants[selected] == target_submodule_invariants[target_selected]:
                    # Extend exact per-row signatures only after the cheaper
                    # rank and Smith obstructions have passed. The shared
                    # assignment order keeps source/target coordinates aligned.
                    if next_source_columns is None:
                        source_column = source_pairings.column(choice)
                        next_source_columns = source_columns + (source_column,)
                        source_signature_counts = Counter(
                            (colour, tuple(column[index] for column in next_source_columns))
                            for index, colour in enumerate(source_colours)
                        )
                    assert next_source_columns is not None and source_signature_counts is not None
                    target_column = target_pairings.column(j)
                    next_target_columns = target_columns + (target_column,)
                    if _configuration_signature_multisets_match(
                        source_colours,
                        next_source_columns,
                        target_colours,
                        next_target_columns,
                        source_signature_counts=source_signature_counts,
                    ):
                        child_source_columns = next_source_columns
                        source_signature_counts = None
                        next_source_columns = None
                        backtrack(images, child_source_columns, next_target_columns)
                        if first_witness is not None:
                            return first_witness
                del images[choice]
            return None

        backtrack({}, (), ())
        return first_witness if first_only else tuple(stabilizer_elements)


class LorentzianPerfectComplex:
    r"""Complete quotient traversal of one Lorentzian perfect-domain complex.

    A quotient adjacency stores a selected target representative and the exact
    isometry carrying it onto the actual neighbor obtained by flipping the
    retained source facet.
    """

    def __init__(
        self,
        lattice: Lattices.ParentMethods,
        allowed_squares: RealSet,
        *,
        retain_facet_orbits: bool = True,
    ) -> None:
        if allowed_squares not in (NONNEGATIVE_SQUARES, ZERO_SQUARES):
            raise ValueError(f"unsupported normalized-norm locus {allowed_squares!r}")
        self._lattice = lattice
        self._allowed_squares = allowed_squares
        self._retain_facet_orbits = retain_facet_orbits
        self._local_backend = LorentzianPerfectLocalBackend()
        self._traversal_representatives: list[LorentzianPerfectCell] | None = None
        self._traversal_representative_buckets: dict[PerfectFormKey, list[LorentzianPerfectCell]] = {}
        self._traversal_adjacencies: list[LorentzianCellAdjacency] = []
        self._traversal_processed: set[int] = set()
        self._traversal_facet_orbits = {}
        self._traversal_facet_generators = {}

    def lattice(self) -> Lattices.ParentMethods:
        return self._lattice

    def allowed_squares(self) -> RealSet:
        return self._allowed_squares

    def local_backend(self) -> LorentzianPerfectLocalBackend:
        return self._local_backend

    def _ensure_traversal_started(self) -> None:
        if self._traversal_representatives is not None:
            return
        backend = self.local_backend()
        initial = backend.initial_cell(self.lattice(), self.allowed_squares())
        self._traversal_representatives = [initial]
        self._traversal_representative_buckets = {_perfect_form_hash_key(initial, backend): [initial]}

    def _expand_traversal_once(self) -> tuple[LorentzianPerfectCell, ...]:
        self._ensure_traversal_started()
        representatives = self._traversal_representatives
        if representatives is None or len(self._traversal_processed) == len(representatives):
            return ()
        backend = self.local_backend()
        # LatticeDelaunays.g:125-145 selects the *smallest unprocessed*
        # Delaunay cell by its number of vertices, breaking ties by insertion
        # order. It does not expand the representatives in FIFO order.
        source_position = min(
            (index for index in range(len(representatives)) if index not in self._traversal_processed),
            key=lambda index: (len(representatives[index].rows), index),
        )
        source = representatives[source_position]
        new_representatives: list[LorentzianPerfectCell] = []
        facet_masks, generator_permutations, canonical_facet_keys = backend._facet_mask_data(
            source, use_owned_cone=self._retain_facet_orbits
        )
        source_id = id(source)
        # Side-pairing witnesses are only queried for actual adjacencies, so
        # retain the Schreier words for those facet orbits alone. Yielding one
        # orbit at a time avoids keeping a second full decomposition alive
        # while neighbors are constructed; boundary-face orbits are discarded.
        if self._retain_facet_orbits:
            self._traversal_facet_orbits[source_id] = {}
            self._traversal_facet_generators[source_id] = generator_permutations
        for orbit in backend._iter_facet_mask_orbits(
            facet_masks,
            generator_permutations,
            retain_words=self._retain_facet_orbits,
            canonical_keys=canonical_facet_keys,
        ):
            facet = _facet_incidence(orbit[0][0], len(source.rows))
            neighbor = backend.flip_across(source, facet)
            if neighbor is None:
                # A boundary face has no adjacent perfect cell or side pairing.
                continue
            target = None
            transporter = None
            neighbor_key = _perfect_form_hash_key(neighbor, backend)
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
            # LatticeDelaunays.g:145-176 records one adjacency for each
            # stabilizer orbit of facets, not an adjacency for every facet.
            # The other facets are carried to the selected one by the cell
            # stabilizer, whose generators are already retained separately.
            # Recording every orbit member would replicate the same quotient
            # edge and its side pairing, with an expensive polyhedral witness
            # check for each copy.
            if self._retain_facet_orbits:
                self._traversal_facet_orbits[source_id][facet] = orbit
            self._traversal_adjacencies.append(LorentzianCellAdjacency(source, facet, target, transporter))
        self._traversal_processed.add(source_position)
        return tuple(new_representatives)

    def explored_cell_count(self) -> int:
        r"""The number of quotient cells whose facets the traversal has already flipped."""
        return len(self._traversal_processed)

    def isometry_to(self, target: LorentzianPerfectComplex) -> LatticeIsometryMethods | None:
        if self.allowed_squares() != target.allowed_squares():
            raise ValueError("an isometry of perfect-domain complexes must preserve the prescribed normalized-norm locus")
        # An integral isometry preserves the rank, inertia and determinant.
        # Check these before either complex starts its potentially expensive
        # perfect-domain enumeration; they are necessary even when the
        # configuration hash happens to coincide.
        source_gram = _gram_matrix(self.lattice())
        target_gram = _gram_matrix(target.lattice())
        if source_gram.nrows() != target_gram.nrows() or source_gram.det() != target_gram.det() or self.lattice().signature_pair() != target.lattice().signature_pair():
            return None
        self._ensure_traversal_started()
        target._ensure_traversal_started()
        representatives = self._traversal_representatives
        target_representatives = target._traversal_representatives
        if representatives is None or target_representatives is None:
            raise ArithmeticError("a started Lorentzian traversal has no initial cell")
        target_cell = target_representatives[0]
        backend = self.local_backend()
        target_key = _perfect_form_hash_key(target_cell, backend)
        for source_cell in self._traversal_representative_buckets.get(target_key, ()):
            witness = backend.cell_transporter(source_cell, target_cell)
            if witness is not None:
                return witness
        while len(self._traversal_processed) < len(representatives):
            for source_cell in self._expand_traversal_once():
                if _perfect_form_hash_key(source_cell, backend) == target_key:
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
        while len(self._traversal_processed) < len(representatives):
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


type PerfectFormKey = tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
    tuple[int, ...],
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
]


def _perfect_form_hash_key(
    cell: LorentzianPerfectCell,
    backend: LorentzianPerfectLocalBackend | None = None,
) -> PerfectFormKey:
    r"""ComputeInvariantPerfectForm (lorentzian_perfect.h:1212-1256), with exact integral invariants.

    For moderate configurations, this uses the source's full histogram of
    diagonal and off-diagonal Gram pairings, streamed in row blocks rather
    than materializing an n-by-n matrix. For large configurations, diagonal
    norms, row contents, Gram-covector divisibilities and Smith factors remain
    in the key while only the quadratic off-diagonal histogram is omitted.
    This is only a traversal bucket key; an exact transporter check remains
    mandatory on every match.
    """
    if backend is not None and len(cell.rows) > _PAIRING_GRAPH_REFINEMENT_LIMIT:
        # The large-cell transporter needs exactly the same norm-only colour
        # and Smith invariants as this traversal bucket. Reuse that backend's
        # identity cache instead of constructing the dense row/covector
        # matrices and computing both invariants a second time. The cache then
        # also serves the subsequent stabilizer/transporter operation.
        _rows, _pairings, spectra, smith_invariants, _row_sum, _row_moment, _distinct_rows = backend._integral_configuration_invariants(cell)
        diagonal = Counter(norm for norm, _content, _divisibility in spectra)
        vector_contents = tuple(sorted(Counter(content for _norm, content, _divisibility in spectra).items()))
        divisibilities = tuple(sorted(Counter(divisibility for _norm, _content, divisibility in spectra).items()))
        return (
            tuple(sorted(diagonal.items())),
            (),
            tuple(int(entry) for entry in smith_invariants),
            vector_contents,
            divisibilities,
        )
    pairings = _ConfigurationPairings(cell.rows, _gram_matrix(cell.lattice))
    diagonal: Counter[int] = Counter()
    for block in pairings.diagonal_blocks():
        diagonal.update(int(value) for value in block)
    off_diagonal: Counter[int] = Counter()
    if len(cell.rows) <= _PAIRING_GRAPH_REFINEMENT_LIMIT:
        for start, block in pairings.blocks():
            for offset, row in enumerate(block):
                i = start + offset
                # The pinned loop counts each unordered pair once, without sorting
                # or retaining an array of Python scalar objects. Bounded row
                # blocks reduce the FLINT dispatch count while keeping the exact
                # n-by-n Gram matrix unmaterialized.
                off_diagonal.update(int(row[j]) for j in range(i + 1, len(pairings)))
    smith_invariants = tuple(int(entry) for entry in pairings.rows.transpose().elementary_divisors())
    if len(cell.rows) <= _PAIRING_GRAPH_REFINEMENT_LIMIT:
        vector_contents = tuple(sorted(Counter(gcd(*row) for row in cell.rows).items()))
        divisibilities = tuple(
            sorted(Counter(gcd(*pairings.covectors.row(index)) for index in range(len(pairings))).items())
        )
    else:
        # Row content and Gram-covector divisibility are linear-cost exact
        # invariants. Their gcd reductions are cheap for fixed ambient rank;
        # unlike the omitted off-diagonal histogram, they do not scale
        # quadratically with the number of configuration vectors.
        vector_contents = tuple(sorted(Counter(gcd(*map(int, row)) for row in cell.rows).items()))
        divisibilities = tuple(
            sorted(
                Counter(
                    gcd(*map(int, pairings.covectors.row(index)))
                    for index in range(len(pairings))
                ).items()
            )
        )
    return (
        tuple(sorted(diagonal.items())),
        tuple(sorted((value, count) for value, count in off_diagonal.items() if count)),
        smith_invariants,
        vector_contents,
        divisibilities,
    )


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
    complex_ = perfect_complex(lattice, allowed_squares)
    backend = complex_.local_backend()
    cones = RationalPolyhedralCones(lattice)
    cells = complex_.quotient_cells()
    # Construct each owned cone only when its facet incidences are processed.
    # Eager construction runs Normaliz's complete support-hyperplane algorithm
    # for every quotient cell before releasing any per-cell native temporaries.
    # Keep the resulting owned cones for the later adjacency consumers.
    cell_cones = {}
    # A supporting facet is the zero section of its defining covector.
    # Reuse the source cell's H-presentation for every member of its
    # stabilizer orbit instead of recomputing a V-to-H representation.
    # Only quotient-adjacency representatives need stored normals. The
    # traversal already carries exact Schreier words for every facet in each
    # adjacency orbit, so retaining a covector for every orbit member is
    # redundant O(number-of-facets * rank) owned data.
    facet_covectors_by_incidence = {}
    facet_keys_by_cell = {}
    # Every quotient cell is represented in this same ambient lattice; reuse
    # its canonical coordinate labels instead of rebuilding the identical
    # generator tuple once per cell.
    labels = tuple(lattice.module_generating_set())
    for cell in cells:
        cell_cone = backend.cell_cone(cell)
        cell_cones[id(cell)] = cell_cone
        ray_count = len(cell.rows)
        # The upstream Lorentzian perfect algorithm computes one facet per
        # stabilizer orbit (lorentzian_perfect.h::f_adj); the Python traversal
        # has already retained those exact orbit representatives and words.
        # Keep compact incidence keys for duplicate validation, but retain
        # covectors solely for representatives used by quotient adjacencies.
        representative_keys = {
            orbit[0][0]
            for orbit in complex_._traversal_facet_orbits.get(id(cell), {}).values()
        }
        adjacency_facet_keys = {
            orbit_key
            for orbit in complex_._traversal_facet_orbits.get(id(cell), {}).values()
            for orbit_key, _word in orbit
        }
        seen_incidence_keys = set()
        boundary_facet_keys = set()
        by_incidence = {}
        # The owned cone is already required and retained above. Re-querying
        # Normaliz for raw facets would materialize a second complete H-list;
        # reuse its exact minimized covectors while keeping incidence pairing
        # blockwise, so no full facet-by-ray matrix is constructed.
        covectors = cell_cone.halfspace_covectors()
        covector_at = covectors.ranking_map().inverse()
        # IndexedSymbolicFunctionSet.__getitem__ rebuilds this inverse on
        # every call. Retain its owned unranking map once for the whole cell;
        # Python slicing is not part of that indexed-family interface.
        # The cone's stored H-presentation already is the minimized facet
        # family. `facet_covectors()` constructs every facet cone merely to
        # recover these same normals, which is prohibitively expensive on
        # high-facet cells.
        # Evaluate bounded blocks of exact covectors against the entire ray
        # matrix. Constructing the full facet-by-rank H matrix before pairing
        # retains another dense allocation proportional to the facet count.
        # Build each coordinate block directly from its owned covectors.
        block_size = _configuration_pairing_block_size(ray_count)
        ray_transpose = None
        pairing_block = None
        covector_block = None
        all_nonincident_mask = (1 << ray_count) - 1
        covector_blocks = (
            (start // block_size, tuple(covector_at(position) for position in range(start, min(start + block_size, len(covectors)))))
            for start in range(0, len(covectors), block_size)
        )
        for block_index, covector_block in covector_blocks:
            start = block_index * block_size
            stop = min(start + block_size, len(covectors))
            # Delay the rank-by-ray transpose until a nonempty facet block is
            # available; an empty H-presentation needs no ray matrix.
            if ray_transpose is None:
                ray_transpose = matrix(
                    SageZZ, len(labels), ray_count,
                    lambda i, j: cell.rows[j][i],
                )
            # Bind only the current covector block.
            # Materialize only this bounded block through the owned unranking
            # map, reusing its inverse instead of rebuilding it per covector.
            # Every owned halfspace covector is primitive integral data in
            # the lattice dual, so exact incidence needs no QQ coercion.
            halfspace_block = matrix(
                SageZZ,
                len(covector_block),
                len(labels),
                lambda row, column: (
                    covector_block[row](labels[column])
                ),
            )
            pairing_block = halfspace_block * ray_transpose
            del halfspace_block
            for offset in range(stop - start):
                nonincident_mask = 0
                # Enumerate the row directly rather than allocating a list of
                # all nonincident positions; dense facet complements can make
                # that list large. The temporary row view is consumed inline.
                for position, pairing in enumerate(pairing_block.row(offset)):
                    if pairing:
                        nonincident_mask |= 1 << position
                if nonincident_mask == all_nonincident_mask:
                    raise ArithmeticError("a minimized H-presentation covector has no incident extreme rays")
                key = _facet_complement_key(nonincident_mask, ray_count)
                if key in seen_incidence_keys:
                    raise ArithmeticError("distinct perfect-cell facets have the same extreme-ray incidence")
                seen_incidence_keys.add(key)
                if key in representative_keys:
                    by_incidence[key] = covector_block[offset]
                elif key not in adjacency_facet_keys:
                    # Boundary facets have no side-pairing orbit retained by
                    # the quotient traversal. Keep only their compact keys;
                    # all nonboundary keys are already owned by Schreier orbits.
                    boundary_facet_keys.add(key)
            # Release the completed dense pairing block before the next one.
            # Drop the tuple as well: assignment to the next block would
            # otherwise build its replacement while this previous block is
            # still referenced, briefly doubling peak slice memory.
            del pairing_block, covector_block
        # Release the transposed configuration and covector slice after the
        # final incidence block, before adjacency construction begins.
        # ``pairing_block`` was deleted after every nonempty block above.
        # Do not delete it again here: the final iteration leaves it unbound.
        del covector_block, ray_transpose, covectors, covector_at
        if set(by_incidence) != representative_keys:
            raise ArithmeticError("a quotient-adjacency facet has no defining covector in its source cone")
        facet_covectors_by_incidence[id(cell)] = by_incidence
        # Do not retain a second full set of every facet key: the quotient
        # traversal already owns every nonboundary orbit key. Only boundary
        # facets need a separate retained set for the later completeness test.
        facet_keys_by_cell[id(cell)] = boundary_facet_keys
        del seen_incidence_keys, adjacency_facet_keys
    ambient_labels = tuple(lattice.module_generating_set())
    ambient_dual = lattice.dual_module()
    dual_labels = tuple(ambient_dual.module_generating_set())
    ambient_gram = _gram_matrix(lattice)
    neighbor_cache = {}
    # Containment checks need a target's H-coordinate matrix only while one
    # adjacency is certified. Keeping one dense matrix per quotient cell
    # duplicates every owned cone's H-presentation across the whole traversal.
    # Retain only the most recently requested target matrix.
    target_halfspace_cache_key = None
    target_halfspace_cache = None
    common_face_cache = {}
    identity_action = identity_matrix(SageZZ, len(ambient_labels))
    ambient_basis = tuple(lattice.module_generator(label) for label in ambient_labels)
    source_basis_data = {}
    adjacencies = []
    paired_facets = {id(cell): set() for cell in cells}
    for adjacency in complex_.adjacencies():
        # The GAP quotient (LatticeDelaunays.g:147-171) retains only one
        # edge per facet orbit.  The preamble's complete domain object also
        # requires a witness for each actual facet.  Recover those witnesses
        # using the Schreier transversal already computed for this orbit.
        orbit = complex_._traversal_facet_orbits[id(adjacency.source)][adjacency.facet]
        generator_permutations = complex_._traversal_facet_generators[id(adjacency.source)]
        transporter_action = matrix(
            SageZZ,
            [_coordinate_row(lattice, adjacency.transporter(basis_vector), SageZZ) for basis_vector in ambient_basis],
        )
        source_key = id(adjacency.source)
        if source_key not in source_basis_data:
            source_rows = matrix(SageQQ, adjacency.source.rows)
            # Pivot columns of the transpose select independent original
            # rows in one exact elimination, instead of recomputing rank for
            # every candidate prefix of the extreme-ray family.
            basis_positions = list(source_rows.transpose().pivots())
            if len(basis_positions) != len(ambient_labels):
                raise ArithmeticError("a source perfect-cell configuration does not span the ambient lattice")
            source_basis_data[source_key] = (tuple(basis_positions), source_rows.matrix_from_rows(basis_positions).inverse())
        basis_positions, inverse_source_basis = source_basis_data[source_key]
        source_all_ray_mask = (1 << len(adjacency.source.rows)) - 1
        for incident_mask, symmetry_word in orbit:
            incident = _facet_incidence(incident_mask, len(adjacency.source.rows))
            # In row coordinates g(v)=vM. A covector column a cutting out
            # a>=0 on the target becomes M^{-1}a on its image. Since g is
            # integral unimodular, the transformed H-data remain integral.
            # This exact transformation avoids a new V-to-H conversion for
            # every member of a facet-stabilizer orbit.
            image_positions = list(basis_positions)
            generator_indices = []
            node = symmetry_word
            while node is not None:
                node, generator_index = node
                generator_indices.append(generator_index)
            for generator_index in reversed(generator_indices):
                generator = generator_permutations[generator_index]
                image_positions = [generator[position] for position in image_positions]
            image_basis = matrix(SageQQ, [adjacency.source.rows[position] for position in image_positions])
            symmetry_action = inverse_source_basis * image_basis
            if any(entry.denominator() != 1 for entry in symmetry_action.list()):
                raise ArithmeticError("a Schreier stabilizer permutation does not induce an integral lattice action")
            symmetry_action = symmetry_action.change_ring(SageZZ)
            # In row coordinates, symmetry * transporter acts by the product
            # M_transporter M_symmetry. The compact Schreier word determines
            # the exact integral action on the selected basis; check a
            # nontrivial Schreier element in each orbit against the owned map.
            action = transporter_action * symmetry_action
            if symmetry_action * ambient_gram * symmetry_action.transpose() != ambient_gram or abs(symmetry_action.det()) != 1:
                raise ArithmeticError("a Schreier stabilizer permutation does not induce an integral lattice isometry")
            if action * ambient_gram * action.transpose() != ambient_gram or abs(action.det()) != 1:
                raise ArithmeticError("a Schreier side-pairing action is not an integral lattice isometry")
            # The Schreier word determines the exact action matrix.
            # Defer the owned categorical isometry until a caller requests the
            # side pairing, rather than constructing one for every orbit facet.
            target_to_neighbour_factory = (
                lambda symmetry_action=symmetry_action, transporter=adjacency.transporter, lattice=lattice:
                lattice.O()(tuple(lattice([int(entry) for entry in row]) for row in symmetry_action.rows())) * transporter
            )
            target_cone = cell_cones[id(adjacency.target)]
            target_key = id(adjacency.target)
            if target_halfspace_cache_key != target_key:
                target_covectors = tuple(target_cone.halfspace_covectors())
                target_halfspace_cache = matrix(
                    SageQQ,
                    len(target_covectors),
                    len(ambient_labels),
                    lambda row, column: _rational(target_covectors[row](ambient_labels[column])),
                )
                del target_covectors
                target_halfspace_cache_key = target_key
            target_halfspace_rows = target_halfspace_cache
            action_key = (id(adjacency.target), tuple(int(entry) for entry in action.list()))

            def neighbor_factory(action_key=action_key, target_cone=target_cone, action=action):
                cached = neighbor_cache.get(action_key)
                if cached is not None:
                    return cached
                if action == identity_action:
                    neighbor = target_cone
                else:
                    # Transport the target's minimized H-presentation by
                    # M^{-1} on covector columns (the row action is x -> xM).
                    # This expensive H-to-V conversion is deferred until a
                    # consumer actually asks for the concrete neighbor cone;
                    # the traversal certificate below needs only exact rays.
                    inverse_action = action.inverse().change_ring(SageQQ)

                    def transport_covector(covector):
                        coordinates = covector.to_vector()
                        column = vector(
                            SageQQ,
                            (_rational(coordinates(label)) for label in dual_labels),
                        ).column()
                        transformed = inverse_action * column
                        if any(entry.denominator() != 1 for entry in transformed.column(0)):
                            raise ArithmeticError("a unimodular side pairing produced a nonintegral facet covector")
                        return ambient_dual.linear_combination(
                            {
                                label: ambient_dual.base_ring()(int(entry))
                                for label, entry in zip(dual_labels, transformed.column(0), strict=True)
                                if entry
                            }
                        )

                    neighbor = cones.from_halfspaces(
                        tuple(transport_covector(covector) for covector in target_cone.halfspace_covectors()),
                        equation_covectors=tuple(
                            transport_covector(covector) for covector in target_cone.equation_covectors()
                        ),
                    )
                neighbor_cache[action_key] = neighbor
                return neighbor
            # The Schreier orbit records the exact source facet incidence;
            # its cached supporting covector cuts out that face.
            source_ray_count = len(adjacency.source.rows)
            incident_mask = sum(1 << position for position in incident)
            facet_key = _facet_complement_key(
                ((1 << source_ray_count) - 1) ^ incident_mask, source_ray_count
            )
            try:
                representative_normal = facet_covectors_by_incidence[id(adjacency.source)][orbit[0][0]]
            except KeyError as exc:
                raise ArithmeticError("a quotient-adjacency facet has no defining covector in its source cone") from exc
            # The stored normal belongs to the orbit representative. Transport
            # it to this Schreier facet using the same exact row action that
            # transports its ray basis: for x -> xM, covectors transform by
            # M^{-1}. This is why only representative normals need be retained.
            if symmetry_word is None:
                facet_normal = representative_normal
            else:
                representative_column = vector(
                    SageQQ,
                    (_rational(representative_normal.to_vector()(label)) for label in dual_labels),
                ).column()
                transformed_column = symmetry_action.inverse().change_ring(SageQQ) * representative_column
                if any(entry.denominator() != 1 for entry in transformed_column.column(0)):
                    raise ArithmeticError("a Schreier facet action produced a nonintegral covector")
                facet_normal = ambient_dual.linear_combination(
                    {
                        label: ambient_dual.base_ring()(int(entry))
                        for label, entry in zip(dual_labels, transformed_column.column(0), strict=True)
                        if entry
                    }
                )
            facet_normal_row = vector(
                SageQQ, (_rational(facet_normal.to_vector()(label)) for label in ambient_labels)
            )
            # Membership in this dictionary is itself the exact incidence
            # certificate: its keys were computed by evaluating every stored
            # H-covector on every extreme ray and checked to have codimension
            # one above. Repeating that full ray scan for each Schreier member
            # is quadratic in the number of facets and adds no information.
            # Retain the exact extreme-ray incidence as the face certificate.
            # The shared face cone is built below from the source's existing
            # H-presentation and this facet's supporting covector.
            face_rows = tuple(adjacency.source.rows[position] for position in sorted(incident))

            face_key = (id(adjacency.source), incident)

            def common_face_factory(face_key=face_key, source_cone=cell_cones[id(adjacency.source)], facet_normal=facet_normal):
                cached = common_face_cache.get(face_key)
                if cached is None:
                    # The facet is already given by this exact supporting
                    # covector. Intersect the existing H-presentation rather
                    # than converting its known rays through a fresh V/H
                    # polyhedron; every Schreier representative of this face
                    # shares the resulting owned cone.
                    cached = source_cone.face_on_covector(facet_normal)
                    common_face_cache[face_key] = cached
                return cached

            adjacencies.append(
                _CertifiedPerfectDomainOrbitAdjacency(
                    cell_cones[id(adjacency.source)],
                    cell_cones[id(adjacency.target)],
                    target_to_neighbour_factory,
                    facet_normal_row=facet_normal_row,
                    face_rows=face_rows,
                    target_halfspace_rows=target_halfspace_rows,
                    action=action,
                    neighbor_factory=neighbor_factory,
                    common_face_factory=common_face_factory,
                )
            )
            paired_facets[id(adjacency.source)].add(facet_key)
    # The target H-coordinate matrix is needed only while constructing and
    # certifying adjacencies. Do not keep the final dense cache alive during
    # stabilizer assembly and traversal-object construction.
    del target_halfspace_cache, target_halfspace_cache_key
    cones_by_position = tuple(cell_cones[id(cell)] for cell in cells)
    # The index set is the actual cone objects, but their source cells are
    # already known by identity. Searching the cone family with ``index``
    # invokes polyhedral equality repeatedly and can turn stabilizer assembly
    # quadratic in the number of quotient cells.
    cell_by_cone_identity = {id(cone): cell for cell, cone in zip(cells, cones_by_position, strict=True)}
    stabilizers = finite_indexed_family(
        finite_ordered_set(cones_by_position),
        lambda cone: backend.cell_stabilizer(cell_by_cone_identity[id(cone)]),
    )
    # The generic traversal constructor checks completeness by materializing
    # every facet as a categorical cone. That repeats expensive subobject
    # construction even though the minimized H-presentation above already
    # identifies every facet exactly by its extreme-ray incidence. Perform
    # the same completeness test on those incidences, retaining the precise
    # Lorentzian exception for totally isotropic exterior-boundary faces.
    # Reuse the same exact all-position mask for every candidate in a cell.
    for cell in cells:
        expected_paired_facets = {
            facet_key
            for orbit in complex_._traversal_facet_orbits.get(id(cell), {}).values()
            for facet_key, _word in orbit
        }
        if paired_facets[id(cell)] != expected_paired_facets:
            raise ArithmeticError("quotient side-pairings do not cover every retained facet orbit")
        # Boundary facets have no quotient adjacency and therefore no retained
        # normal, but their exact incidence keys must still participate in the
        # completeness check. Keep this independent of the reduced normal map.
        by_incidence = facet_keys_by_cell[id(cell)]
        ray_count = len(cell.rows)
        all_positions_mask = (1 << ray_count) - 1
        for facet_key in by_incidence:
            if facet_key in paired_facets[id(cell)]:
                continue
            # Iterate incident positions directly from the compact key.
            # The full incidence set is only needed on the error path.
            if isinstance(facet_key, int):
                incident_positions = _bit_positions(all_positions_mask ^ facet_key)
            elif isinstance(facet_key[1], int):
                incident_mask = (
                    all_positions_mask ^ facet_key[1]
                    if facet_key[0]
                    else facet_key[1]
                )
                incident_positions = _bit_positions(incident_mask)
            elif facet_key[0]:
                excluded = set(facet_key[1])
                incident_positions = (position for position in range(ray_count) if position not in excluded)
            else:
                incident_positions = iter(facet_key[1])
            # Only unmatched boundary candidates require a Gram certificate.
            # Constructing the entire ray matrix duplicates the full A7
            # configuration when a much smaller selected-row matrix suffices.
            # Total isotropy is a property of the span, so it suffices to
            # check a basis of the incident rays, not every pair of the full
            # facet. Exact minor checks find at most two independent rays,
            # avoiding a quadratic number of Gram products for large facets.
            facet_basis = []
            for position in incident_positions:
                row = cell.rows[position]
                if not any(row):
                    continue
                if not facet_basis:
                    facet_basis.append(row)
                    # A non-null ray already disproves total isotropy of the
                    # facet span; no further incidence scan is required.
                    if sum(
                        row[i] * ambient_gram[i, j] * row[j]
                        for i in range(len(row))
                        for j in range(len(row))
                    ):
                        break
                    continue
                # With a nonzero first row, rank two is equivalent to some
                # 2-by-2 minor being nonzero. Avoid a Sage rank computation
                # for every incident ray of a potentially enormous facet.
                first = facet_basis[0]
                if any(
                    first[i] * row[j] != first[j] * row[i]
                    for i in range(len(row))
                    for j in range(i + 1, len(row))
                ):
                    facet_basis.append(row)
                    break
            basis_matrix = matrix(SageQQ, facet_basis)
            basis_gram = basis_matrix * ambient_gram * basis_matrix.transpose()
            if not any(basis_gram.list()):
                continue
            incidence = _facet_incidence(facet_key, ray_count)
            raise ArithmeticError(
                f"the perfect-domain decomposition of {lattice} is incomplete: "
                f"the interior facet with extreme-ray incidence {sorted(incidence)} "
                f"of cell {cell} has no adjacent domain"
            )
    # All constructor invariants have been checked above using exact incidence
    # and Gram-matrix certificates, without constructing every face cone.
    traversal = _CertifiedPerfectDomainTraversal.__new__(_CertifiedPerfectDomainTraversal)
    traversal._lattice = lattice
    traversal._cells = finite_ordered_set(tuple(cell_cones.values()))
    traversal._adjacencies = finite_ordered_set(tuple(adjacencies))
    traversal._stabilizer_generators = stabilizers
    if traversal._cells.cardinality() == 0:
        raise ValueError(f"the perfect-domain decomposition of {lattice} has no orbit representatives")
    if any(cell.ambient_lattice() is not lattice for cell in traversal._cells):
        raise ValueError(f"the perfect-domain decomposition of {lattice} contains a cell in a different lattice")
    if any(adjacency.source() not in traversal._cells for adjacency in traversal._adjacencies):
        raise ValueError("a perfect-domain adjacency starts at a cell outside the orbit representatives")
    if any(adjacency.target_representative() not in traversal._cells for adjacency in traversal._adjacencies):
        raise ValueError("a perfect-domain adjacency ends at a cell outside the orbit representatives")
    return traversal


class NoLocalMarkTheoremError(RuntimeError):
    """The local perfect-cell mark theorem does not cover this norm: vectors outside the closed positive cone."""


class MarkedCellOrbitAlgorithm:
    r"""O(L)-orbits of primitive vectors of square n in the closed positive cone, from the quotient cells.

    For n = 0 this is LORENTZ_GetOrbitRepresentative_Kernel (lorentzian_perfect.h:1625-1835):
    the marks of a cell are its null extreme rays, marks equivalent under the cell's automorphism
    group are one local orbit, a mark on a shared facet is identified with its preimage under the
    side pairing, and the connected components of the resulting graph are the orbits.

    The source stops for n != 0 (lorentzian_perfect.h:1656-1661). For timelike n the same graph
    is built with local orbits enumerated up to the cell stabilizer G. The barycentric
    subdivision of a cell is G-invariant and the stabilizer of a maximal flag of faces fixes the
    barycenters of its faces, so the cones on the barycenters of one flag per G-orbit cover the
    cell modulo G. A point of the cell determines its chain of faces and barycentric coefficients
    intrinsically (greedy descent through minimal faces), so two marks are G-equivalent exactly
    when their coefficients agree and some element of G carries one chain onto the other. The
    marks on the representative facet of an adjacency are enumerated in the same way up to the
    facet's stabilizer, which also stabilizes the neighbouring cell. Each simplicial cone is
    finite in marks: its rays pair nonnegatively, so q = n bounds every coefficient.
    """

    def __init__(self, complex_: LorentzianPerfectComplex, norm: int) -> None:
        sign = AttackProfile.from_lattice(complex_.lattice()).sign
        if sign * norm < 0:
            raise NoLocalMarkTheoremError(f"vectors of square {norm} in {complex_.lattice()} are spacelike: they lie outside the closed positive cone that the cells cover")
        self._complex = complex_
        self._norm = norm
        self._local_mark_cache: dict[int, tuple[Lattices.ElementMethods, ...]] = {}
        self._local_mark_row_cache: dict[int, tuple[tuple[int, ...], ...]] = {}

    def complex(self) -> LorentzianPerfectComplex:
        return self._complex

    def local_marks(self, cell: LorentzianPerfectCell) -> tuple[Lattices.ElementMethods, ...]:
        key = id(cell)
        if key in self._local_mark_cache:
            return self._local_mark_cache[key]
        if self._norm != 0:
            raise ValueError("timelike marks are enumerated up to the cell stabilizer, not listed per cell")
        # Null vectors of a pointed cone in one closed Lorentzian sheet
        # lie on its null extreme rays. The configuration need not consist
        # of primitive generators; filter the actual cone's primitive rays
        # by exact Gram square before constructing owned elements.
        zero = cell.lattice.base_ring().zero()
        marks = self.complex().local_backend().cell_primitive_rays(cell, square=int(zero))
        self._local_mark_cache[key] = marks
        return marks

    def local_stabilizer_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        marks = self.local_marks(cell)
        backend = self.complex().local_backend()
        permutations = backend.point_permutations(cell, marks)
        unvisited = set(range(len(marks)))
        orbits = []
        while unvisited:
            orbit = [min(unvisited)]
            unvisited.discard(orbit[0])
            for position in orbit:
                for permutation in permutations:
                    if permutation[position] in unvisited:
                        unvisited.discard(permutation[position])
                        orbit.append(permutation[position])
            orbits.append(tuple(marks[position] for position in sorted(orbit)))
        return tuple(orbits)

    def _mark_rows(self, cell: LorentzianPerfectCell) -> tuple[tuple[int, ...], ...]:
        key = id(cell)
        if key not in self._local_mark_row_cache:
            self._local_mark_row_cache[key] = tuple(tuple(int(entry) for entry in _coordinate_row(cell.lattice, mark, SageZZ)) for mark in self.local_marks(cell))
        return self._local_mark_row_cache[key]

    def transport_marks(self, adjacency: LorentzianCellAdjacency) -> tuple[tuple[Lattices.ElementMethods, Lattices.ElementMethods], ...]:
        transporter = adjacency.transporter
        if transporter is None:
            raise ArithmeticError("a marked quotient adjacency needs a transporter")
        lattice = adjacency.source.lattice
        source = adjacency.source.rows
        # The facet is the face of the source cone cut out by the normal of
        # its rays' span, so a mark of the source cell lies on the facet
        # exactly when it is orthogonal (in coordinates) to that normal.
        facet_rows = matrix(SageQQ, [source[position] for position in sorted(adjacency.facet)])
        normals = facet_rows.right_kernel_matrix()
        if normals.nrows() != 1:
            raise ArithmeticError("a perfect-cell facet must span a hyperplane")
        normal = normals.row(0)
        source_marks = self.local_marks(adjacency.source)
        source_by_coordinates = {row: mark for row, mark in zip(self._mark_rows(adjacency.source), source_marks, strict=True) if vector(SageQQ, row).dot_product(normal) == 0}
        target_marks = self.local_marks(adjacency.target)
        labels = tuple(lattice.module_generating_set())
        # In row coordinates the transporter acts by x -> x * action.
        action = matrix(SageZZ, [_coordinate_row(lattice, transporter(lattice.module_generator(label)), SageZZ) for label in labels])
        transported = []
        for row, target_mark in zip(self._mark_rows(adjacency.target), target_marks, strict=True):
            image = tuple(int(entry) for entry in vector(SageZZ, row) * action)
            if image in source_by_coordinates:
                transported.append((source_by_coordinates[image], target_mark))
        if len(transported) != len(source_by_coordinates):
            raise ArithmeticError("a facet mark has no unique mark in the target representative")
        return tuple(transported)

    def global_orbits(self) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        r"""The orbits, each as a tuple of marks whose first entry is its representative.

        For n = 0 each tuple lists the null marks of the orbit in all quotient cells; for
        timelike n it lists one mark per local orbit of the component.
        """
        if self._norm != 0:
            return self._timelike_global_orbits()
        cells = self.complex().quotient_cells()
        local = tuple(self.local_stabilizer_orbits(cell) for cell in cells)
        nodes = [(cell_position, orbit_position) for cell_position, cell_orbits in enumerate(local) for orbit_position in range(len(cell_orbits))]
        graph: dict[tuple[int, int], set[tuple[int, int]]] = {node: set() for node in nodes}
        local_mark_positions = tuple(
            {tuple(int(entry) for entry in _coordinate_row(cell.lattice, mark, SageZZ)): orbit_position for orbit_position, orbit in enumerate(cell_orbits) for mark in orbit}
            for cell, cell_orbits in zip(cells, local, strict=True)
        )

        def local_orbit_position(cell_position: int, mark: Lattices.ElementMethods) -> int:
            key = tuple(int(entry) for entry in _coordinate_row(cells[cell_position].lattice, mark, SageZZ))
            try:
                return local_mark_positions[cell_position][key]
            except KeyError as exc:
                raise ArithmeticError("a transported mark does not belong to a local stabilizer orbit") from exc

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
            seen_marks: set[tuple[int, ...]] = set()
            for cell_position, orbit_position in sorted(component):
                for mark in local[cell_position][orbit_position]:
                    key = tuple(int(entry) for entry in _coordinate_row(cells[cell_position].lattice, mark, SageZZ))
                    if key not in seen_marks:
                        seen_marks.add(key)
                        marks.append(mark)
            global_orbits.append(tuple(marks))
        return tuple(global_orbits)

    def _timelike_global_orbits(self) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        complex_ = self.complex()
        backend = complex_.local_backend()
        cells = complex_.quotient_cells()
        structures = tuple(_CellFaces(backend, cell) for cell in cells)
        local: list[list[_Mark]] = []
        for structure in structures:
            representatives: list[_Mark] = []
            for flag in structure.flag_simplices(structure.full_face, structure.group):
                for row in _simplicial_vectors_of_square(structure.gram, [structure.barycenter(face) for face in flag], self._norm):
                    if gcd(*(int(entry) for entry in row)) == 1 and structure.classify(row, representatives) is None:
                        representatives.append(structure.mark(row))
            local.append(representatives)
        nodes = [(cell_position, orbit_position) for cell_position, representatives in enumerate(local) for orbit_position in range(len(representatives))]
        graph: dict[tuple[int, int], set[tuple[int, int]]] = {node: set() for node in nodes}
        positions = {id(cell): position for position, cell in enumerate(cells)}
        for adjacency in complex_.adjacencies():
            source_position = positions[id(adjacency.source)]
            target_position = positions[id(adjacency.target)]
            source = structures[source_position]
            target = structures[target_position]
            transporter = adjacency.transporter
            if transporter is None:
                raise ArithmeticError("a marked quotient adjacency needs a transporter")
            lattice = adjacency.source.lattice
            configuration = adjacency.source.rows
            facet_rows = matrix(SageQQ, [configuration[position] for position in sorted(adjacency.facet)])
            normals = facet_rows.right_kernel_matrix()
            if normals.nrows() != 1:
                raise ArithmeticError("a perfect-cell facet must span a hyperplane")
            normal = normals.row(0)
            facet = frozenset(index for index, row in enumerate(source.rows) if normal.dot_product(row) == 0)
            if facet not in source.incidences:
                raise ArithmeticError("an adjacency facet is not a facet of its source cone")
            facet_group = libgap.Stabilizer(source.group, sorted(index + 1 for index in facet), libgap.OnSets)
            # The transporter carries the target representative onto the
            # neighbour across the facet; in row coordinates g(x) = x * action.
            action = matrix(SageZZ, [_coordinate_row(lattice, transporter(lattice.module_generator(label)), SageZZ) for label in source.labels])
            inverse = action.change_ring(SageQQ).inverse()
            if any(entry.denominator() != 1 for entry in inverse.list()):
                raise ArithmeticError("a side pairing must be an integral automorphism")
            inverse = inverse.change_ring(SageZZ)
            for flag in source.flag_simplices(facet, facet_group):
                for row in _simplicial_vectors_of_square(source.gram, [source.barycenter(face) for face in flag], self._norm):
                    if gcd(*(int(entry) for entry in row)) != 1:
                        continue
                    source_orbit = source.classify(row, local[source_position])
                    target_orbit = target.classify(row * inverse, local[target_position])
                    if source_orbit is None or target_orbit is None:
                        raise ArithmeticError("a facet mark lies in no local orbit of its cells")
                    graph[(source_position, source_orbit)].add((target_position, target_orbit))
                    graph[(target_position, target_orbit)].add((source_position, source_orbit))
        unseen = set(nodes)
        orbits: list[tuple[Lattices.ElementMethods, ...]] = []
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
            orbits.append(tuple(structures[cell_position].element(local[cell_position][orbit_position].row) for cell_position, orbit_position in sorted(component)))
        return tuple(orbits)


@dataclass(frozen=True)
class _Mark:
    r"""A mark of a cell with its intrinsic chain of faces and barycentric coefficients."""

    row: FreeModuleElement[Integer]
    chain: tuple[frozenset[int], ...]
    coefficients: tuple[Rational, ...]


class _CellFaces:
    r"""The faces of a perfect cell as sets of extreme-ray indices, with the cell stabilizer acting on them."""

    def __init__(self, backend: LorentzianPerfectLocalBackend, cell: LorentzianPerfectCell) -> None:
        lattice = cell.lattice
        self.lattice = lattice
        self.labels = tuple(lattice.module_generating_set())
        self.gram = _gram_matrix(lattice)
        # The cell already owns the exact primitive integer ray rows. Keep
        # their order directly instead of materializing 49,861 lattice
        # elements and immediately extracting the same coordinates again.
        self.rows = tuple(vector(SageZZ, row) for row in cell.rows)
        # The traversal/complex collector already owns this cell's cone.
        # Consume its exact H-family directly into the dense coordinate matrix
        # below; retaining a tuple of owned covector elements alongside that
        # matrix duplicates one Python object and coordinate vector per facet.
        interior = sum(self.rows, vector(SageZZ, [0] * len(self.labels)))
        # Compute exact facet-to-ray incidences in bounded matrix blocks.
        # The scalar Python loop performs one exact dot product per ray and
        # facet; A7 has 49,861 rays, so batching Sage's exact matrix product
        # avoids millions of Python-level pairing calls without retaining the
        # full facet-by-ray pairing matrix.
        # Both the primitive rays and primitive integral dual covectors have
        # integer coordinates. Keep this incidence product over ZZ so Sage's
        # exact integer matrix backend can use its optimized integer kernels
        # instead of doing rational arithmetic for a zero/nonzero test.
        ray_matrix = matrix(SageZZ, self.rows)
        ray_transpose = ray_matrix.transpose()
        del ray_matrix
        covector_rows = []
        for owned_covector in backend.cell_cone(cell).halfspace_covectors():
            row = []
            for label in self.labels:
                entry = _rational(owned_covector.to_vector()(label))
                if entry.denominator() != 1:
                    raise ArithmeticError(
                        "primitive integral facet covectors must have integral lattice coordinates"
                    )
                row.append(int(entry))
            if sum(entry * coordinate for entry, coordinate in zip(row, interior)) <= 0:
                raise ArithmeticError("the facet covectors of a cell must be positive on its interior")
            covector_rows.append(row)
        covector_matrix = matrix(SageZZ, covector_rows)
        del covector_rows
        all_ray_positions = frozenset(range(len(self.rows)))
        block_size = _configuration_pairing_block_size(len(self.rows))
        incidences = []
        for start in range(0, covector_matrix.nrows(), block_size):
            stop = min(start + block_size, covector_matrix.nrows())
            pairing_block = covector_matrix.matrix_from_rows(range(start, stop)) * ray_transpose
            for offset in range(stop - start):
                incidence = all_ray_positions.difference(
                    pairing_block.row(offset).nonzero_positions()
                )
                if not incidence:
                    raise ArithmeticError("a minimized H-presentation covector has no incident extreme rays")
                incidences.append(frozenset(incidence))
        if len(set(incidences)) != len(incidences):
            raise ArithmeticError("distinct perfect-cell facets have the same extreme-ray incidence")
        self.facet_incidences = tuple(incidences)
        del incidences
        self.incidences = frozenset(self.facet_incidences)
        self.all_facet_positions = frozenset(range(len(self.facet_incidences)))
        self.full_face = frozenset(range(len(self.rows)))
        self.group = libgap.Group(
            [
                libgap.PermList([image + 1 for image in permutation])
                for permutation in backend.point_permutations_from_rows(cell, cell.rows)
            ],
            libgap.eval("()"),
        )
        self._ranks: dict[frozenset[int], int] = {}
        self._barycenter_pairings: dict[frozenset[int], FreeModuleElement[Rational]] = {}
        self._covector_matrix = covector_matrix

    def rank(self, face: frozenset[int]) -> int:
        if face not in self._ranks:
            self._ranks[face] = int(matrix(SageQQ, [self.rows[index] for index in sorted(face)]).rank())
        return self._ranks[face]

    def facets(self, face: frozenset[int]) -> frozenset[frozenset[int]]:
        r"""The facets of a face: its intersections with cell facets that drop the rank by one."""
        dimension = self.rank(face)
        return frozenset(part for incidence in self.incidences for part in (face & incidence,) if part != face and part and self.rank(part) == dimension - 1)

    def barycenter(self, face: frozenset[int]) -> FreeModuleElement[Integer]:
        indices = sorted(face)
        return sum((self.rows[index] for index in indices[1:]), self.rows[indices[0]])

    def flag_simplices(self, top: frozenset[int], group) -> tuple[tuple[frozenset[int], ...], ...]:
        r"""One maximal flag of faces below ``top`` per orbit of ``group``, which must stabilize ``top``."""
        flags: list[tuple[frozenset[int], ...]] = []

        def descend(chain: tuple[frozenset[int], ...], stabilizer) -> None:
            face = chain[-1]
            if self.rank(face) == 1:
                flags.append(chain)
                return
            facets = sorted(sorted(index + 1 for index in facet) for facet in self.facets(face))
            for orbit in libgap.Orbits(stabilizer, facets, libgap.OnSets):
                representative = [int(point) for point in orbit[0]]
                descend((*chain, frozenset(point - 1 for point in representative)), libgap.Stabilizer(stabilizer, representative, libgap.OnSets))

        descend((top,), group)
        return tuple(flags)

    def mark(self, row: FreeModuleElement[Integer]) -> _Mark:
        r"""The chain and coefficients of a point: descend through the minimal faces containing it.

        At a point p of the relative interior of its minimal face F with barycenter b, the
        coefficient of b is the largest t with p - t b in the cell; the remainder lies in a
        proper face of F.
        """
        point = vector(SageQQ, row)
        pairings = self._covector_matrix * point
        chain: list[frozenset[int]] = []
        coefficients: list[Rational] = []
        while point != 0:
            tight = self.all_facet_positions.difference(pairings.nonzero_positions())
            # Begin with the smallest tight facet incidence. Intersecting
            # large facet sets first needlessly copies many more ray indices
            # at every step of the mark-chain descent.
            face = min((self.facet_incidences[index] for index in tight), key=len) if tight else self.full_face
            for index in tight:
                incidence = self.facet_incidences[index]
                if incidence is not face:
                    face = face & incidence
            barycenter = vector(SageQQ, self.barycenter(face))
            if face not in self._barycenter_pairings:
                self._barycenter_pairings[face] = self._covector_matrix * barycenter
            barycenter_pairings = self._barycenter_pairings[face]
            coefficient = min(
                pairing / barycenter_pairing
                for pairing, barycenter_pairing in zip(pairings, barycenter_pairings, strict=True)
                if barycenter_pairing > 0
            )
            chain.append(face)
            coefficients.append(coefficient)
            point -= coefficient * barycenter
            pairings -= coefficient * barycenter_pairings
        return _Mark(row, tuple(chain), tuple(coefficients))

    def classify(self, row: FreeModuleElement[Integer], representatives: list[_Mark]) -> int | None:
        r"""The position of the representative G-equivalent to the point, or None."""
        mark = self.mark(row)
        chain = [sorted(index + 1 for index in face) for face in mark.chain]
        for position, representative in enumerate(representatives):
            if representative.coefficients != mark.coefficients or tuple(map(len, representative.chain)) != tuple(map(len, mark.chain)):
                continue
            image = [sorted(index + 1 for index in face) for face in representative.chain]
            if libgap.RepresentativeAction(self.group, chain, image, libgap.OnTuplesSets) != libgap.eval("fail"):
                return position
        return None

    def element(self, row: FreeModuleElement[Integer]) -> Lattices.ElementMethods:
        ring = self.lattice.base_ring()
        return self.lattice.linear_combination({label: ring(int(entry)) for label, entry in zip(self.labels, row, strict=True) if entry})


def _simplicial_vectors_of_square(ambient_gram: Matrix_integer_dense, ambient_rays: list[FreeModuleElement[Integer]], norm: int) -> tuple[FreeModuleElement[Integer], ...]:
    r"""The coordinate rows of the lattice vectors of square n in the simplicial cone on ambient_rays.

    The cone lies in the saturation W of its rays' span, where the search runs in coordinates of a
    basis of W. Let R be the matrix of the piece's primitive rays r_i and D = |det R|. A lattice
    vector of the piece is v = sum c_i r_i with C = D c in Z_{>=0}^k, and C mod D runs
    over the D residues of L / (sum Z r_i) (Smith form of R). The rays pair with one
    sign, so after that sign every pairing b_ij is nonnegative and
    D^2 q(v) = sum C_i C_j b_ij is monotone in each C_i. A timelike ray has
    C_i^2 b_ii <= D^2 |n|; a null ray needs a positively paired partner with C_j >= 1,
    so 2 C_i b_ij <= D^2 |n|. Each residue class is enumerated by depth-first search
    within these bounds, pruning a branch once its square exceeds |n|.
    """
    spanning = matrix(SageZZ, ambient_rays)
    saturation = spanning.saturation()
    if saturation.nrows() != spanning.nrows() or spanning.rank() != spanning.nrows():
        raise ArithmeticError("the rays of a simplicial cone must be independent")
    local = saturation.solve_left(spanning)
    if any(entry.denominator() != 1 for entry in local.list()):
        raise ArithmeticError("the rays of a cone must lie in the saturation of their span")
    gram = saturation * ambient_gram * saturation.transpose()
    ray_rows = [vector(SageZZ, row) for row in local.rows()]
    # Search the timelike rays first: a null ray's coefficient raises
    # the square only once a positively paired coefficient is nonzero.
    ray_rows.sort(key=lambda row: (row * gram * row) == 0)
    rank = len(ray_rows)
    ray_matrix = matrix(SageZZ, ray_rows)
    pairing_matrix = ray_matrix * gram * ray_matrix.transpose()
    pairings = [[int(pairing_matrix[i, j]) for j in range(rank)] for i in range(rank)]
    # A simplicial face contained in a null ray has identically zero
    # quadratic form. It contributes no positive-norm marks, regardless
    # of the ambient lattice's signature. In particular, choosing the
    # sign from the zero pairing matrix would be mathematically invalid.
    if all(entry == 0 for row in pairings for entry in row):
        return ()
    sign = 1 if sum(sum(row) for row in pairings) > 0 else -1
    signed = [[sign * entry for entry in row] for row in pairings]
    if any(entry < 0 for row in signed for entry in row):
        raise ArithmeticError(f"the rays {ray_rows} of a cell's cone do not pair with one sign, as rays of the closed positive cone do")
    target = sign * int(norm)
    if target <= 0:
        raise ArithmeticError("positive-square marks require a nonzero norm with the sign of the cone")
    determinant = abs(int(ray_matrix.det()))
    if determinant == 0:
        raise ArithmeticError("the rays of a simplicial piece must be independent")
    scaled_target = target * determinant * determinant
    bounds = []
    for position in range(rank):
        if signed[position][position] > 0:
            bounds.append(int(Integer(scaled_target // signed[position][position]).isqrt()))
        else:
            partners = [signed[position][other] for other in range(rank) if other != position and signed[position][other] > 0]
            if not partners:
                raise ArithmeticError(f"the isotropic ray {ray_rows[position]} has no positive pairing in its simplicial piece")
            # A positive-norm vector involving this null ray must have
            # a nonzero coefficient on at least one positively paired
            # ray. Zero pairings with other rays impose no constraint;
            # they must not be used as the minimum positive pairing.
            bounds.append(scaled_target // (2 * min(partners)))
    # Residues of L modulo the ray span: with S = U R V, the rows of
    # V^{-1} form a basis W of L with sum Z r_i = sum s_i Z W_i.
    smith, _left, right = ray_matrix.smith_form()
    basis = right.inverse()
    adjugate = (ray_matrix.change_ring(SageQQ).inverse() * determinant).change_ring(SageZZ)
    diagonal = [int(smith[index, index]) for index in range(rank)]
    found: list[FreeModuleElement[Integer]] = []

    last = rank - 1

    def record(scaled: list[int]) -> None:
        row = vector(SageZZ, scaled) * ray_matrix
        if any(entry % determinant for entry in row):
            raise ArithmeticError("a residue-class combination of the rays is not a lattice vector")
        found.append(vector(SageZZ, [int(entry) // determinant for entry in row]) * saturation)

    def search(scaled: list[int], paired: list[int], square: int, position: int) -> None:
        if position == last:
            # The square is a polynomial of degree at most two in the
            # last step t: D^2 q = square + 2 D t p + D^2 t^2 b, with
            # p = (bC)_last and b = b_last,last. Solve it for t >= 0.
            linear = 2 * determinant * paired[last]
            quadratic = determinant * determinant * signed[last][last]
            constant = square - scaled_target
            if quadratic == 0:
                if linear == 0:
                    if constant == 0:
                        raise ArithmeticError("a positive-square mark is orthogonal to a null ray of its cone")
                    return
                steps = (-constant,) if -constant % linear == 0 else ()
                steps = tuple(step // linear for step in steps)
            else:
                discriminant = linear * linear - 4 * quadratic * constant
                if discriminant < 0:
                    return
                root = isqrt(discriminant)
                if root * root != discriminant:
                    return
                steps = tuple({(-linear + sign * root) // (2 * quadratic) for sign in (1, -1) if (-linear + sign * root) % (2 * quadratic) == 0})
            start = scaled[last]
            for step in steps:
                if step >= 0 and start + determinant * step <= bounds[last]:
                    scaled[last] = start + determinant * step
                    record(scaled)
            scaled[last] = start
            return
        start = scaled[position]
        while scaled[position] <= bounds[position] and square <= scaled_target:
            search(scaled, paired, square, position + 1)
            # D^2 q(C + D e_i) = D^2 q(C) + 2 D (bC)_i + D^2 b_ii.
            square += 2 * determinant * paired[position] + determinant * determinant * signed[position][position]
            for other in range(rank):
                paired[other] += determinant * signed[other][position]
            scaled[position] += determinant
        steps = (scaled[position] - start) // determinant
        for other in range(rank):
            paired[other] -= steps * determinant * signed[other][position]
        scaled[position] = start

    for digits in product(*(range(entry) for entry in diagonal)):
        representative = vector(SageZZ, digits) * basis
        # D c = representative * adj(R) (up to the sign of det R); its
        # residues mod D are the scaled coefficients of the class.
        numerators = representative * adjugate
        scaled = [int(entry) % determinant for entry in numerators]
        paired = [sum(signed[row][column] * scaled[column] for column in range(rank)) for row in range(rank)]
        square = sum(scaled[row] * paired[row] for row in range(rank))
        search(scaled, paired, square, 0)
    return tuple(found)


def _coordinate_row[S: RingElement](lattice: FramedFreeModules.ParentMethods, element: FramedFreeModules.ElementMethods, ring: Parent[S]) -> FreeModuleElement[S]:
    coordinates = element.to_vector()
    # A preamble scalar belongs to its own coefficient ring, not Sage's ZZ/QQ
    # parents.  Cross the private engine boundary by exact rational components;
    # Sage's coercion discovery cannot convert those categorical scalars.
    return vector(ring, tuple(ring(_rational(coordinates(label))) for label in lattice.module_generating_set()))


def _positive_ldl_factorization(gram: Matrix_rational_dense) -> tuple[list[list[Rational]], list[Rational]]:
    """Exact ``G = L D L^T`` factorization of a positive-definite Gram matrix."""
    if gram.nrows() == 0 or gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("the fixed-distance LDL factorization requires a nonempty symmetric Gram matrix")
    if gram.base_ring() is not SageQQ:
        raise ValueError("the fixed-distance LDL factorization requires an exact QQ Gram matrix")
    rank = int(gram.nrows())
    zero = cast(Rational, SageQQ.zero())
    lower: list[list[Rational]] = [[zero for _ in range(rank)] for _ in range(rank)]
    diagonal: list[Rational] = []
    for column in range(rank):
        lower[column][column] = SageQQ.one()
        pivot = cast(
            Rational,
            gram[column, column]
            - sum(
                (lower[column][previous] ** 2 * diagonal[previous] for previous in range(column)),
                zero,
            ),
        )
        if pivot <= 0:
            raise ArithmeticError("the fixed-distance sphere requires a positive-definite Gram matrix")
        diagonal.append(pivot)
        for row in range(column + 1, rank):
            lower[row][column] = cast(
                Rational,
                (
                    gram[row, column]
                    - sum(
                        (lower[row][previous] * lower[column][previous] * diagonal[previous] for previous in range(column)),
                        zero,
                    )
                )
                / pivot,
            )
    return lower, diagonal


def _validate_fixed_distance_basis_inverse(
    basis_change: Matrix_rational_dense,
    inverse: Matrix_rational_dense,
) -> None:
    """Certify an exact two-sided inverse before transporting a sphere center."""
    rank = int(basis_change.nrows())
    if basis_change.ncols() != rank or inverse.nrows() != rank or inverse.ncols() != rank:
        raise ValueError("the fixed-distance sphere basis inverse has incorrect dimensions")
    if basis_change.base_ring() is not SageQQ or inverse.base_ring() is not SageQQ:
        raise ValueError("the fixed-distance sphere basis inverse must be over exact QQ")
    identity = identity_matrix(SageQQ, rank)
    if basis_change * inverse != identity or inverse * basis_change != identity:
        raise ArithmeticError("the fixed-distance sphere basis inverse is not a two-sided inverse")


def _validate_fixed_distance_ldl_factorization(
    gram: Matrix_rational_dense,
    factorization: tuple[list[list[Rational]], list[Rational]],
) -> None:
    """Validate a cached unit-lower LDL factorization before sphere pruning."""
    rank = int(gram.nrows())
    lower, diagonal = factorization
    if len(lower) != rank or len(diagonal) != rank or any(len(row) != rank for row in lower):
        raise ValueError("the fixed-distance LDL factorization has incorrect dimensions")
    entries = [entry for row in lower for entry in row] + list(diagonal)
    if any(not isinstance(entry, (int, Integer, Rational)) for entry in entries):
        raise ValueError("the fixed-distance LDL factorization must use exact rational entries")
    for row in range(rank):
        if lower[row][row] != 1 or any(lower[row][column] != 0 for column in range(row + 1, rank)):
            raise ArithmeticError("the fixed-distance LDL factorization is not unit lower triangular")
        if diagonal[row] <= 0:
            raise ArithmeticError("the fixed-distance LDL factorization has a nonpositive pivot")
    lower_matrix = matrix(SageQQ, lower)
    diagonal_matrix_ = diagonal_matrix(SageQQ, diagonal)
    if lower_matrix * diagonal_matrix_ * lower_matrix.transpose() != gram:
        raise ArithmeticError("the fixed-distance LDL factorization does not reconstruct its Gram matrix")


def _fixed_distance_coordinate_solutions(
    gram: Matrix_rational_dense,
    target: tuple[Rational, ...],
    squared_distance: Rational,
    ldl_factorization: tuple[list[list[Rational]], list[Rational]] | None = None,
) -> tuple[tuple[tuple[Integer, ...], Rational], ...]:
    r"""Enumerate exactly the lattice points on a rational-center sphere.

    The Gram matrix is positive definite. This is the coordinate-level
    ``fixed_dist_vectors`` operation used by isotropic Lorentzian shells; it
    avoids enumerating the entire ellipsoid and filtering it afterward.
    """
    if gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("fixed-distance coordinate enumeration requires a symmetric square Gram matrix")
    if gram.base_ring() is not SageQQ:
        raise ValueError("fixed-distance coordinate enumeration requires an exact QQ Gram matrix")
    rank = int(gram.nrows())
    if len(target) != rank:
        raise ValueError("fixed-distance target must have the Gram matrix dimension")
    if any(not isinstance(entry, (int, Integer, Rational)) for entry in target):
        raise ValueError("fixed-distance target must have exact rational coordinates")
    if not isinstance(squared_distance, (int, Integer, Rational)):
        raise ValueError("fixed squared distance must be exact rational data")
    target = tuple(SageQQ(entry) for entry in target)
    bound = SageQQ(squared_distance)
    if bound < 0:
        raise ValueError("fixed squared distance must be nonnegative")
    if rank == 0:
        if ldl_factorization is not None:
            _validate_fixed_distance_ldl_factorization(gram, ldl_factorization)
        return (((), bound),) if bound == 0 else ()
    factorization = ldl_factorization if ldl_factorization is not None else _positive_ldl_factorization(gram)
    lower, diagonal = factorization
    _validate_fixed_distance_ldl_factorization(gram, factorization)
    zero = SageQQ.zero()
    coordinates = [0] * rank
    solutions: list[tuple[tuple[Integer, ...], Rational]] = []

    def enumerate_coordinate(index: int, accumulated: Rational) -> None:
        if index < 0:
            if accumulated == bound:
                solutions.append((tuple(SageZZ(entry) for entry in coordinates), bound))
            return
        remaining = bound - accumulated
        if remaining < 0:
            return
        center = target[index] - sum(
            (lower[later][index] * (coordinates[later] - target[later]) for later in range(index + 1, rank)),
            zero,
        )
        squared_radius = remaining / diagonal[index]
        if index == 0:
            # At the innermost coordinate the remaining equation is exactly
            # D_0 (z-center)^2 = remaining. Solve it directly instead of
            # traversing every integer in the enclosing interval; this is
            # the distinction between a fixed-distance sphere and an
            # at-most-distance ellipsoid enumeration.
            numerator = int(squared_radius.numerator())
            denominator = int(squared_radius.denominator())
            root_numerator = isqrt(numerator)
            root_denominator = isqrt(denominator)
            if root_numerator * root_numerator != numerator or root_denominator * root_denominator != denominator:
                return
            radius = SageQQ(root_numerator) / SageQQ(root_denominator)
            candidates = (center,) if radius == 0 else (center - radius, center + radius)
            for candidate in candidates:
                if candidate.denominator() == 1:
                    coordinates[0] = int(candidate)
                    solutions.append((tuple(SageZZ(entry) for entry in coordinates), bound))
            return
        radius = isqrt(int(squared_radius.numerator() // squared_radius.denominator()))
        floor_center = int(center.numerator() // center.denominator())
        ceil_center = -int((-center).numerator() // (-center).denominator())
        first = ceil_center - radius - 1
        last = floor_center + radius + 1
        while first <= last and diagonal[index] * (SageQQ(first) - center) ** 2 > remaining:
            first += 1
        while first <= last and diagonal[index] * (SageQQ(last) - center) ** 2 > remaining:
            last -= 1
        for coordinate in range(first, last + 1):
            term = diagonal[index] * (SageQQ(coordinate) - center) ** 2
            if term <= remaining:
                coordinates[index] = coordinate
                enumerate_coordinate(index - 1, accumulated + term)

    enumerate_coordinate(rank - 1, zero)
    return tuple(solutions)


def _bounded_distance_coordinate_solutions(
    gram: Matrix_rational_dense,
    target: tuple[Rational, ...],
    squared_distance_bound: Rational,
    ldl_factorization: tuple[list[list[Rational]], list[Rational]] | None = None,
) -> tuple[tuple[tuple[Integer, ...], Rational], ...]:
    r"""Enumerate every integer point in an exact rational-center ellipsoid.

    This is the bounded-CVP operation needed by positive-norm Lorentzian
    shells. It uses exact LDL pruning rather than the provider's PARI-backed
    ``close_vectors`` path, whose stack can grow beyond its configured limit
    on the published E10 vector-orbit cases. Every returned distance is the
    exact value of ``(x-target)^T gram (x-target)``.
    """
    if gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("bounded-distance coordinate enumeration requires a symmetric square Gram matrix")
    if gram.base_ring() is not SageQQ:
        raise ValueError("bounded-distance coordinate enumeration requires an exact QQ Gram matrix")
    rank = int(gram.nrows())
    if len(target) != rank or any(not isinstance(entry, (int, Integer, Rational)) for entry in target):
        raise ValueError("bounded-distance target must have exact rational coordinates of the Gram dimension")
    if not isinstance(squared_distance_bound, (int, Integer, Rational)):
        raise ValueError("bounded squared distance must be exact rational data")
    target = tuple(SageQQ(entry) for entry in target)
    bound = SageQQ(squared_distance_bound)
    if bound < 0:
        raise ValueError("bounded squared distance must be nonnegative")
    if rank == 0:
        if ldl_factorization is not None:
            _validate_fixed_distance_ldl_factorization(gram, ldl_factorization)
        return (((), SageQQ.zero()),)
    factorization = ldl_factorization if ldl_factorization is not None else _positive_ldl_factorization(gram)
    lower, diagonal = factorization
    _validate_fixed_distance_ldl_factorization(gram, factorization)
    zero = SageQQ.zero()
    coordinates = [0] * rank
    solutions: list[tuple[tuple[Integer, ...], Rational]] = []

    def enumerate_coordinate(index: int, accumulated: Rational) -> None:
        if index < 0:
            solutions.append((tuple(SageZZ(entry) for entry in coordinates), accumulated))
            return
        remaining = bound - accumulated
        if remaining < 0:
            return
        center = target[index] - sum(
            (lower[later][index] * (coordinates[later] - target[later]) for later in range(index + 1, rank)),
            zero,
        )
        squared_radius = remaining / diagonal[index]
        radius = isqrt(int(squared_radius.numerator() // squared_radius.denominator()))
        floor_center = int(center.numerator() // center.denominator())
        ceil_center = -int((-center).numerator() // (-center).denominator())
        first = ceil_center - radius - 1
        last = floor_center + radius + 1
        while first <= last and diagonal[index] * (SageQQ(first) - center) ** 2 > remaining:
            first += 1
        while first <= last and diagonal[index] * (SageQQ(last) - center) ** 2 > remaining:
            last -= 1
        for coordinate in range(first, last + 1):
            term = diagonal[index] * (SageQQ(coordinate) - center) ** 2
            if term <= remaining:
                coordinates[index] = coordinate
                enumerate_coordinate(index - 1, accumulated + term)

    enumerate_coordinate(rank - 1, zero)
    return tuple(solutions)


def _fixed_distance_dual_lll_solutions(
    gram: Matrix_rational_dense,
    target: tuple[Rational, ...],
    squared_distance: Rational,
    dual_transform: Matrix_integer_dense,
    enumeration_order_transform: Matrix_integer_dense,
    ldl_factorization: tuple[list[list[Rational]], list[Rational]] | None = None,
) -> tuple[tuple[tuple[Integer, ...], Rational], ...]:
    """Enumerate an exact sphere in the pinned dual-LLL frame, then reconstruct.

    If ``Q`` is the dual transform and ``P`` the recursive coordinate map,
    the pinned identity is ``Q * P^T = I`` and the reduced Gram is ``Q G Q^T``.
    Enumerating in that frame preserves the source's reduced-basis pruning;
    returned coordinates are transported back by ``Q^T`` and ordered by
    their original pinned recursive coordinates ``P*x``.
    """
    if gram.nrows() != gram.ncols() or gram.base_ring() is not SageQQ:
        raise ValueError("fixed-distance dual-LLL enumeration requires a square exact QQ Gram matrix")
    rank = int(gram.nrows())
    if len(target) != rank or any(not isinstance(entry, (int, Integer, Rational)) for entry in target):
        raise ValueError("fixed-distance dual-LLL target must have exact rational coordinates of the Gram dimension")
    if not isinstance(squared_distance, (int, Integer, Rational)) or SageQQ(squared_distance) < 0:
        raise ValueError("fixed-distance dual-LLL squared distance must be exact and nonnegative")
    for transform in (dual_transform, enumeration_order_transform):
        if transform.base_ring() not in (SageZZ, SageQQ) or any(
            not isinstance(entry, (int, Integer, Rational)) or SageQQ(entry).denominator() != 1
            for entry in transform.list()
        ):
            raise ValueError("fixed-distance dual-LLL transforms must have exact integer entries")
    if dual_transform.nrows() != rank or dual_transform.ncols() != rank:
        raise ValueError("the fixed-distance dual-LLL transform has incorrect dimensions")
    if enumeration_order_transform.nrows() != rank or enumeration_order_transform.ncols() != rank:
        raise ValueError("the fixed-distance enumeration-order transform has incorrect dimensions")
    if dual_transform * enumeration_order_transform.transpose() != identity_matrix(SageZZ, rank):
        raise ArithmeticError("the fixed-distance dual-LLL transforms are not inverse-transposes")
    reduced_gram = dual_transform * gram * dual_transform.transpose()
    transformed_target_column = enumeration_order_transform * vector(SageQQ, target).column()
    transformed_target = tuple(SageQQ(transformed_target_column[i, 0]) for i in range(rank))
    reduced_solutions = _fixed_distance_coordinate_solutions(
        reduced_gram, transformed_target, squared_distance, ldl_factorization=ldl_factorization
    )
    inverse_basis = dual_transform.transpose()
    solutions = []
    for coordinates, distance in reduced_solutions:
        original_column = inverse_basis * vector(SageZZ, coordinates).column()
        if any(entry.denominator() != 1 for entry in original_column.list()):
            raise ArithmeticError("dual-LLL sphere reconstruction produced nonintegral coordinates")
        original = tuple(SageZZ(original_column[row, 0]) for row in range(rank))
        solutions.append((original, distance))
    return tuple(sorted(solutions, key=lambda item: _dual_lll_enumeration_key(item[0], enumeration_order_transform)))


def _bounded_distance_dual_lll_solutions(
    gram: Matrix_rational_dense,
    target: tuple[Rational, ...],
    squared_distance_bound: Rational,
    dual_transform: Matrix_integer_dense,
    enumeration_order_transform: Matrix_integer_dense,
    ldl_factorization: tuple[list[list[Rational]], list[Rational]] | None = None,
) -> tuple[tuple[tuple[Integer, ...], Rational], ...]:
    """Enumerate a bounded CVP ellipsoid in the pinned dual-LLL frame."""
    if gram.nrows() != gram.ncols() or gram.base_ring() is not SageQQ:
        raise ValueError("bounded-distance dual-LLL enumeration requires a square exact QQ Gram matrix")
    rank = int(gram.nrows())
    if len(target) != rank or any(not isinstance(entry, (int, Integer, Rational)) for entry in target):
        raise ValueError("bounded-distance dual-LLL target must have exact rational coordinates of the Gram dimension")
    if not isinstance(squared_distance_bound, (int, Integer, Rational)) or SageQQ(squared_distance_bound) < 0:
        raise ValueError("bounded-distance dual-LLL bound must be exact and nonnegative")
    for transform in (dual_transform, enumeration_order_transform):
        if transform.base_ring() not in (SageZZ, SageQQ) or any(
            not isinstance(entry, (int, Integer, Rational)) or SageQQ(entry).denominator() != 1
            for entry in transform.list()
        ):
            raise ValueError("bounded-distance dual-LLL transforms must have exact integer entries")
    if dual_transform.nrows() != rank or dual_transform.ncols() != rank:
        raise ValueError("the bounded-distance dual-LLL transform has incorrect dimensions")
    if enumeration_order_transform.nrows() != rank or enumeration_order_transform.ncols() != rank:
        raise ValueError("the bounded-distance enumeration-order transform has incorrect dimensions")
    if dual_transform * enumeration_order_transform.transpose() != identity_matrix(SageZZ, rank):
        raise ArithmeticError("the bounded-distance dual-LLL transforms are not inverse-transposes")
    reduced_gram = dual_transform * gram * dual_transform.transpose()
    transformed_target_column = enumeration_order_transform * vector(SageQQ, target).column()
    transformed_target = tuple(SageQQ(transformed_target_column[i, 0]) for i in range(rank))
    reduced_solutions = _bounded_distance_coordinate_solutions(
        reduced_gram, transformed_target, squared_distance_bound, ldl_factorization=ldl_factorization
    )
    inverse_basis = dual_transform.transpose()
    solutions = []
    for coordinates, distance in reduced_solutions:
        original_column = inverse_basis * vector(SageZZ, coordinates).column()
        if any(entry.denominator() != 1 for entry in original_column.list()):
            raise ArithmeticError("dual-LLL bounded-sphere reconstruction produced nonintegral coordinates")
        original = tuple(SageZZ(original_column[row, 0]) for row in range(rank))
        solutions.append((original, distance))
    return tuple(sorted(solutions, key=lambda item: _dual_lll_enumeration_key(item[0], enumeration_order_transform)))


def _fixed_distance_vectors(
    lattice: Lattices.ParentMethods,
    target: FramedFreeModules.ElementMethods,
    signed_distance: RingElement,
    reduction: LatticeReduction | None = None,
    ldl_factorization: tuple[list[list[Rational]], list[Rational]] | None = None,
    basis_inverse: Matrix_rational_dense | None = None,
) -> tuple[Lattices.ElementMethods, ...]:
    r"""Enumerate the exact distance sphere in a negative-definite lattice.

    This is the fixed-distance branch of ``lorentzian_perfect.h:289-299``.
    Unlike filtering ``close_vectors`` after enumerating its entire sublevel
    ellipsoid, this shares the exact coordinate-level sphere recursion with
    the positive-vector shell path. Both prune on the prescribed distance and
    solve the innermost coordinate equation directly.
    """
    bound = cast(Rational, -_rational(signed_distance))
    if bound < 0:
        raise ValueError("a fixed-distance sphere in a negative-definite lattice needs nonpositive signed distance")
    if int(lattice.module_rank()) == 0:
        return (lattice.zero(),) if bound == 0 else ()
    if reduction is None:
        reduction = lattice.lll_reduction()
    reduced = cast(Lattices.ParentMethods, reduction.reduced)
    change_of_basis = cast(Matrix_integer_dense, _engine_matrix(reduction.change_of_basis_matrix))
    basis_change = cast(Matrix_rational_dense, change_of_basis.change_ring(SageQQ))
    rank = int(basis_change.nrows())
    inverse = basis_inverse if basis_inverse is not None else cast(Matrix_rational_dense, basis_change.inverse())
    _validate_fixed_distance_basis_inverse(basis_change, inverse)
    target_row = cast("FreeModuleElement[Rational]", _coordinate_row(lattice, target, SageQQ))
    target_column = cast(Matrix_rational_dense, inverse * target_row.column())
    reduced_target = cast("FreeModuleElement[Rational]", vector(cast("Parent[Rational]", SageQQ), tuple(target_column[row, 0] for row in range(rank))))
    if basis_change * target_column != target_row.column():
        raise ArithmeticError("the reduced fixed-distance center was not transported by the basis isomorphism")
    gram = cast(Matrix_rational_dense, -_gram_matrix(reduced).change_ring(SageQQ))
    lower, diagonal = ldl_factorization if ldl_factorization is not None else _positive_ldl_factorization(gram)
    solutions = tuple(
        coordinates
        for coordinates, _distance in _fixed_distance_coordinate_solutions(
            gram,
            tuple(SageQQ(entry) for entry in reduced_target),
            bound,
            ldl_factorization=(lower, diagonal),
        )
    )

    isometry = cast(Callable[[Lattices.ElementMethods], Lattices.ElementMethods], reduction.isometry)
    results: list[Lattices.ElementMethods] = []
    original_gram = cast(Matrix_rational_dense, -_gram_matrix(lattice).change_ring(SageQQ))
    for solution in solutions:
        point = isometry(_ambient_element(reduced, tuple(SageZZ(coordinate) for coordinate in solution)))
        displacement = cast("FreeModuleElement[Rational]", _coordinate_row(lattice, point, SageQQ) - target_row)
        actual_distance = (displacement * original_gram * displacement.column())[0]
        if actual_distance != bound:
            raise ArithmeticError("the fixed-distance sphere lost its exact distance under LLL transport")
        results.append(point)
    return tuple(results)


def _rational(entry: RingElement) -> Rational:
    entry = cast(Rational, entry)
    return SageQQ(int(entry.numerator())) / SageQQ(int(entry.denominator()))


def _source_remove_fraction(row: FreeModuleElement[Rational]) -> FreeModuleElement[Rational]:
    r"""Primitive affine normal, as in Fractions.g:38-66."""
    entries = tuple(SageQQ(entry) for entry in row)
    if not any(entries):
        return vector(cast("Parent[Rational]", SageQQ), entries)
    denominator = SageZZ.one()
    for entry in entries:
        denominator = denominator.lcm(SageZZ(entry.denominator()))
    integral = tuple(SageZZ(denominator * entry) for entry in entries)
    content = SageZZ.zero()
    for integral_entry in integral:
        content = content.gcd(integral_entry)
    return vector(cast("Parent[Rational]", SageQQ), (integral_entry // content for integral_entry in integral))


def _outside_ray_directions(
    affine_nullspace: Matrix_rational_dense, gram: Matrix_rational_dense
) -> Matrix_rational_dense:
    r"""Individually primitive covector directions used by GetOneOutsideRay.

    ``Lorentzian.g:1187`` applies ``RemoveFraction`` to each affine-nullspace
    row's covector part after multiplication by the inverse Gram matrix.
    Normalization is rowwise, not a single normalization of the resulting span.
    """
    inverse = gram.inverse()
    return matrix(
        SageQQ,
        [
            _source_remove_fraction(vector(SageQQ, row[1:]) * inverse)
            for row in affine_nullspace.rows()
        ],
    )


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
    labels = lattice.module_generating_set()
    return lattice.linear_combination({label: ring(int(entry)) for label, entry in zip(labels, coordinates, strict=True) if entry})


def _owned_rational(
    lattice: Lattices.ParentMethods, value: Rational, *, field: Parent | None = None
) -> RingElement:
    r"""Raise an exact rational through the selected fraction-field realization."""
    if field is None:
        field = lattice.base_ring().fraction_field()
    return _owned_engine_element(field, SageQQ(value))


def _pinned_diagonalize_symmetric_matrix(
    qgram: Matrix_rational_dense,
) -> tuple[Matrix_rational_dense, Matrix_rational_dense]:
    r"""Port GAP ``DiagonalizeSymmetricMatrix`` for an exact rational form.

    Return the source-ordered congruence transform and its diagonal Gram
    matrix. ``NullspaceReduction`` and each Gram--Schmidt complement use the
    pinned GAP nullspace and greedy column-selection conventions.
    """
    if qgram.nrows() != qgram.ncols() or qgram != qgram.transpose():
        raise ValueError("symmetric diagonalization requires a symmetric square matrix")
    rank = int(qgram.nrows())
    # PositivitySymmetricMatrices.g:27-35 uses GAP NullspaceMat via
    # NullspaceReduction; preserve that row order for the positive-vector
    # diagonalization just as in GetSomeNegativeVector.
    radical = _pinned_nullspace_matrix(qgram)
    radical_rank = int(radical.nrows())
    zero_radical = cast(
        Matrix_rational_dense,
        matrix(cast(Parent, SageQQ), radical_rank, rank, [SageQQ.zero()] * (radical_rank * rank)),
    )
    if radical * qgram != zero_radical:
        raise ArithmeticError("symmetric diagonalization did not compute the Gram radical")
    transform = cast(Matrix_rational_dense, matrix(cast(Parent, SageQQ), rank, rank, 0))
    for row_index in range(radical_rank):
        for column_index in range(rank):
            transform[row_index, column_index] = radical[row_index, column_index]
    pivot_columns = set(_pinned_independent_column_indices(radical))
    row_index = radical_rank
    for column_index in range(rank):
        if column_index in pivot_columns:
            continue
        transform[row_index, column_index] = 1
        row_index += 1
    radical_reduced = cast(Matrix_rational_dense, transform * qgram * transform.transpose())
    if any(radical_reduced[i, j] for i in range(radical_rank) for j in range(rank)):
        raise ArithmeticError("radical complement does not split orthogonally")
    nullspace_transform = cast(Matrix_rational_dense, matrix(cast(Parent, SageQQ), transform))
    nondegenerate_rank = rank - radical_rank
    nondegenerate = cast(
        Matrix_rational_dense,
        radical_reduced.matrix_from_rows_and_columns(range(radical_rank, rank), range(radical_rank, rank)),
    )
    diagonal_basis: list[FreeModuleElement[Rational]] = []
    for basis_index in range(nondegenerate_rank):
        if basis_index == 0:
            orthogonal_basis = cast(
                Matrix_rational_dense,
                identity_matrix(cast(Parent, SageQQ), nondegenerate_rank),
            )
        else:
            previous = cast(Matrix_rational_dense, matrix(SageQQ, diagonal_basis))
            orthogonal_basis = _pinned_nullspace_matrix(
                cast(Matrix_rational_dense, nondegenerate * previous.transpose())
            )
        restricted = cast(Matrix_rational_dense, orthogonal_basis * nondegenerate * orthogonal_basis.transpose())
        orthogonal_rank = int(orthogonal_basis.nrows())
        selected = next(
            (cast("FreeModuleElement[Rational]", vector(SageQQ, orthogonal_rank, {i: 1})) for i in range(orthogonal_rank) if restricted[i, i]),
            None,
        )
        if selected is None:
            selected = next(
                (
                    cast("FreeModuleElement[Rational]", vector(SageQQ, orthogonal_rank, {i: 1, j: 1}))
                    for i in range(orthogonal_rank)
                    for j in range(i + 1, orthogonal_rank)
                    if restricted[i, j]
                ),
                None,
            )
        if selected is None:
            raise ArithmeticError("the nondegenerate restriction has no nonisotropic vector")
        diagonal_basis.append(cast("FreeModuleElement[Rational]", selected * orthogonal_basis))
    if nondegenerate_rank:
        nondegenerate_transform = cast(Matrix_rational_dense, matrix(SageQQ, diagonal_basis))
        for row_index in range(nondegenerate_rank):
            for column_index in range(rank):
                transform[radical_rank + row_index, column_index] = sum(
                    nondegenerate_transform[row_index, j] * nullspace_transform[radical_rank + j, column_index] for j in range(nondegenerate_rank)
                )
    diagonal_matrix = cast(Matrix_rational_dense, transform * qgram * transform.transpose())
    expected_diagonal = cast(Matrix_rational_dense, radical_reduced)
    if nondegenerate_rank:
        expected_nonzero_block = cast(
            Matrix_rational_dense,
            nondegenerate_transform * nondegenerate * nondegenerate_transform.transpose(),
        )
        for i in range(nondegenerate_rank):
            for j in range(nondegenerate_rank):
                expected_diagonal[radical_rank + i, radical_rank + j] = expected_nonzero_block[i, j]
    if diagonal_matrix != expected_diagonal:
        raise ArithmeticError("pinned diagonalization lost its independently assembled congruence witness")
    if any(diagonal_matrix[i, j] for i in range(rank) for j in range(rank) if i != j):
        raise ArithmeticError("pinned symmetric diagonalization left off-diagonal entries")
    if not transform.det():
        raise ArithmeticError("positive-direction diagonalization is singular")

    return transform, diagonal_matrix


def _pinned_approximation_iterator(
    direction: tuple[Rational, ...] | tuple[float, ...],
) -> tuple[
    tuple[tuple[Rational | float, tuple[tuple[int, int], ...]], ...],
    list[Rational | float],
    list[int],
]:
    r"""Construct the source ``ApproxIterator`` state for one direction."""
    groups: dict[Rational | float, list[tuple[int, int]]] = {}
    for coordinate, entry in enumerate(direction):
        if isinstance(entry, float) and not isfinite(entry):
            raise ValueError("ApproxIterator floating directions must be finite")
        if entry:
            if isinstance(entry, float):
                step: Rational | float = 1.0 / abs(entry)
                if not isfinite(step):
                    raise ValueError("ApproxIterator floating steps must be finite")
                sign = 1 if entry > 0 else -1
            else:
                rational_entry = cast(Rational, entry)
                step = SageQQ.one() / abs(rational_entry)
                sign = 1 if rational_entry > 0 else -1
            groups.setdefault(step, []).append((coordinate, sign))
    if not groups:
        raise ValueError("ApproxIterator requires a nonzero direction")
    ordered = tuple((step, tuple(group)) for step, group in sorted(groups.items()))
    nexts = [step / 2.0 if isinstance(step, float) else cast(Rational, step) / SageQQ(2) for step, _group in ordered]
    return ordered, nexts, [0] * len(direction)


def _advance_pinned_approximation_iterator(
    ordered: tuple[tuple[Rational | float, tuple[tuple[int, int], ...]], ...],
    nexts: list[Rational | float],
    approximation: list[int],
) -> list[int]:
    r"""Advance all tied minima exactly as ``ApproxIterator.increment``."""
    if len(nexts) != len(ordered) or not nexts:
        raise ValueError("ApproxIterator state must contain one next-time per step")
    least = min(nexts)
    for group_index, (step, group) in enumerate(ordered):
        if nexts[group_index] != least:
            continue
        nexts[group_index] = nexts[group_index] + step
        for coordinate, sign in group:
            approximation[coordinate] += sign
    return approximation


def _pinned_positive_diagonal_directions(
    transform: Matrix_rational_dense,
    diagonal: Matrix_rational_dense,
) -> tuple[tuple[Rational, ...], ...]:
    r"""Port ``GetPositiveDirections_diag`` from Positivity.h:139-160.

    Scan the diagonal in source order, retain strictly positive entries, and
    normalize each corresponding transformation row by its exact L1 norm.
    """
    if diagonal.nrows() != diagonal.ncols() or transform.nrows() != diagonal.nrows() or transform.ncols() != diagonal.ncols():
        raise ValueError("diagonal directions require matching square transform and Gram matrices")
    if any(diagonal[i, j] for i in range(diagonal.nrows()) for j in range(diagonal.ncols()) if i != j):
        raise ValueError("diagonal directions require a diagonal Gram matrix")
    directions: list[tuple[Rational, ...]] = []
    for index in range(int(diagonal.nrows())):
        if diagonal[index, index] > 0:
            row = tuple(SageQQ(entry) for entry in transform.row(index))
            norm = sum(abs(entry) for entry in row)
            if not norm:
                raise ArithmeticError("a positive diagonal direction has zero transformation row")
            directions.append(tuple(entry / norm for entry in row))
    return tuple(directions)


def _pinned_positive_eigen_directions(
    gram: Matrix_real_double_dense,
) -> tuple[tuple[float, ...], ...]:
    r"""Reproduce ``GetPositiveDirections_eigen``'s indexed row extraction.

    The pinned Eigen code reads ``ListVect(i,j)`` from the eigenvector matrix,
    so its direction is row ``i`` of the matrix whose columns are eigenvectors,
    not the eigenvector column itself. Keep that observable source convention
    explicit rather than substituting Sage's right-eigenvector vectors.
    """
    if gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("the Eigen positive-direction search requires a symmetric square matrix")
    if any(not isfinite(float(entry)) for entry in gram.list()):
        raise ValueError("the Eigen positive-direction search requires finite Gram entries")
    data = cast(
        list[tuple[ComplexDoubleElement, list[FreeModuleElement], int]],
        gram.eigenvectors_right(),
    )
    eigenvalues: list[float] = []
    eigenvector_columns: list[tuple[float, ...]] = []
    for eigenvalue, vectors, multiplicity in data:
        eigenvalue_real = float(eigenvalue.real())
        if not isfinite(eigenvalue_real):
            raise ArithmeticError("eigensolver returned a nonfinite eigenvalue")
        for eigenvector in vectors:
            column = tuple(float(entry.real()) for entry in eigenvector)
            if any(not isfinite(entry) for entry in column):
                raise ArithmeticError("eigensolver returned a nonfinite eigenvector")
            eigenvalues.append(eigenvalue_real)
            eigenvector_columns.append(column)
        if len(vectors) != multiplicity:
            raise ArithmeticError("eigensolver multiplicity disagrees with its returned eigenvectors")
    if len(eigenvalues) != int(gram.nrows()) or any(len(column) != len(eigenvalues) for column in eigenvector_columns):
        raise ArithmeticError("eigensolver returned a dimensionally incomplete eigenbasis")
    return tuple(
        tuple(eigenvector_columns[column][row] for column in range(len(eigenvector_columns)))
        for row, eigenvalue in enumerate(eigenvalues)
        if eigenvalue > 0
    )


def _source_require_integral_matrix(candidate, name: str) -> None:
    """Require an exact coefficient ring and integral entries."""
    if not candidate.base_ring().is_exact():
        raise ValueError(f"{name} must have entries in an exact ring")
    for entry in candidate.list():
        try:
            integer = SageZZ(entry)
        except (TypeError, ValueError, ArithmeticError) as error:
            raise ValueError(f"{name} must have integral entries") from error
        if integer != entry:
            raise ValueError(f"{name} must have integral entries")


def _source_reconstruct_integral_vector(
    gram: Matrix_integer_dense,
    basis: Matrix_integer_dense,
    reduced: Matrix_integer_dense,
    candidate: FreeModuleElement[Integer],
    critical_norm: Rational = SageQQ.zero(),
    strict_inequality: bool = True,
) -> tuple[Integer, ...] | None:
    r"""Lift one V2 candidate and certify the requested quadratic inequality.

    Return ``None`` when the reduced norm fails the requested bound; otherwise
    map row
    coordinates by ``basis.transpose()`` and verify the exact original-form
    witness before the caller returns it.
    """
    if candidate not in SageZZ**reduced.nrows():
        raise ValueError("candidate coordinates must match the reduced Gram matrix")
    reduced_norm = SageQQ((candidate * reduced * candidate.column())[0])
    critical_norm = SageQQ(critical_norm)
    accepted = reduced_norm > critical_norm if strict_inequality else reduced_norm >= critical_norm
    if not accepted:
        return None
    original_column = basis.transpose() * candidate.column()
    result = tuple(SageZZ(original_column[row, 0]) for row in range(int(original_column.nrows())))
    original_vector = vector(SageZZ, result)
    original_norm = SageQQ((original_vector * gram * original_vector.column())[0])
    if original_norm != reduced_norm:
        raise ArithmeticError("integral-vector reduction lost its exact norm witness")
    return result


def _source_integral_vector_allmeth_v2(
    gram: Matrix_integer_dense,
    basis: Matrix_integer_dense,
    reduced: Matrix_integer_dense,
    critical_norm: Rational = SageQQ.zero(),
    strict_inequality: bool = True,
    *,
    diagonalize: Callable[[Matrix_rational_dense], tuple[Matrix_rational_dense, Matrix_rational_dense]] = _pinned_diagonalize_symmetric_matrix,
    eigen_directions: Callable[[Matrix_real_double_dense], tuple[tuple[float, ...], ...]] = _pinned_positive_eigen_directions,
) -> tuple[Integer, ...]:
    r"""Port ``GetIntegralVector_allmeth_V2`` on a supplied reduced form.

    This is the ordered exact-diagonal/eigen ApproxIterator consumer. The
    caller supplies the pinned integral reduction, so this operation can be
    verified independently of the surrounding Sage category imports.
    """
    if gram.nrows() == 0 or gram.nrows() != gram.ncols():
        raise ValueError("the Gram matrix must be nonempty and square")
    _source_require_integral_matrix(gram, "the original Gram matrix")
    _source_require_integral_matrix(reduced, "the reduced Gram matrix")
    if basis.nrows() != gram.nrows() or basis.ncols() != gram.ncols():
        raise ValueError("the reduction basis must be a square matrix")
    _source_require_integral_matrix(basis, "the reduction basis")
    if abs(basis.det()) != 1:
        raise ValueError("the reduction basis must be unimodular over the integer lattice")
    if reduced != basis * gram * basis.transpose() or reduced != reduced.transpose():
        raise ArithmeticError("the supplied reduced Gram matrix is not the basis congruence")
    rank = int(reduced.nrows())
    qgram = cast(Matrix_rational_dense, reduced.change_ring(SageQQ))
    transform, diagonal_matrix = diagonalize(qgram)
    if (
        transform.nrows() != rank
        or transform.ncols() != rank
        or diagonal_matrix.nrows() != rank
        or diagonal_matrix.ncols() != rank
        or not transform.det()
    ):
        raise ArithmeticError("exact diagonalization returned a singular or dimensionally invalid transform")
    if diagonal_matrix != transform * qgram * transform.transpose():
        raise ArithmeticError("exact diagonalization failed its Gram congruence certificate")
    if any(diagonal_matrix[i, j] for i in range(rank) for j in range(rank) if i != j):
        raise ArithmeticError("exact diagonalization returned a nondiagonal Gram matrix")
    if not any(diagonal_matrix[i, i] > 0 for i in range(rank)):
        raise ValueError("the Gram matrix has no positive direction to approximate")

    iterators: list[
        tuple[
            tuple[tuple[Rational | float, tuple[tuple[int, int], ...]], ...],
            list[Rational | float],
            list[int],
        ]
    ] = [
        _pinned_approximation_iterator(direction)
        for direction in _pinned_positive_diagonal_directions(transform, diagonal_matrix)
    ]
    # Positivity.h:171-187 and :335-342. These floating eigen-directions are
    # an optimization only; the exact diagonal iterators above are still the
    # guaranteed terminating search. Test the exact integral quadratic value
    # before accepting any eigen-direction approximation.
    eigen_iterators: list[
        tuple[
            tuple[tuple[Rational | float, tuple[tuple[int, int], ...]], ...],
            list[Rational | float],
            list[int],
        ]
    ] = []
    try:
        double_gram = cast(
            Matrix_real_double_dense,
            matrix(cast(Parent, RDF), reduced),
        )
        for direction in eigen_directions(double_gram):
            eigen_norm = sum(abs(eigen_entry) for eigen_entry in direction)
            # The optional floating branch must not admit NaN or infinite
            # steps into ApproxIterator.
            if eigen_norm and all(isfinite(entry) for entry in direction) and isfinite(eigen_norm):
                normalized = tuple(eigen_entry / eigen_norm for eigen_entry in direction)
                if all(isfinite(entry) and (entry == 0 or abs(entry) >= 1.0 / float_info.max) for entry in normalized):
                    eigen_iterators.append(_pinned_approximation_iterator(normalized))
    except (ArithmeticError, OverflowError, RuntimeError, TypeError, ValueError):
        # This branch only accelerates the exact search. Overflow or failure
        # in the floating eigensolver must not disable its exact directions.
        eigen_iterators.clear()
    if not iterators:
        raise ValueError("the Gram matrix has no positive direction")
    while True:
        for ordered, nexts, approximation in iterators:
            _advance_pinned_approximation_iterator(ordered, nexts, approximation)
            candidate = vector(SageZZ, approximation)
            result = _source_reconstruct_integral_vector(
                gram, basis, reduced, candidate, critical_norm, strict_inequality
            )
            if result is not None:
                return result
        for ordered, nexts, approximation in eigen_iterators:
            _advance_pinned_approximation_iterator(ordered, nexts, approximation)
            candidate = vector(SageZZ, approximation)
            result = _source_reconstruct_integral_vector(
                gram, basis, reduced, candidate, critical_norm, strict_inequality
            )
            if result is not None:
                return result


def _source_remove_fraction_vector(vector_over_q) -> tuple[Integer, ...]:
    r"""Port ``RemoveFractionVector``: clear denominators, then divide content."""
    rational = tuple(SageQQ(entry) for entry in vector_over_q)
    denominator = SageZZ.one()
    for entry in rational:
        denominator = denominator.lcm(entry.denominator())
    integral = tuple(SageZZ(entry * denominator) for entry in rational)
    content = SageZZ.zero()
    for entry in integral:
        content = content.gcd(abs(entry))
    # Fractions.g:46-48 returns the zero vector unchanged with multiplier 1.
    if not content:
        return integral
    return tuple(entry // content for entry in integral)


def _source_validate_integral_reduction(
    gram: Matrix_integer_dense,
    basis: Matrix_integer_dense,
    reduced: Matrix_integer_dense,
) -> None:
    r"""Certify the row-basis congruence consumed by GetIntegralVector_allmeth."""
    rank = int(gram.nrows())
    if gram.ncols() != rank or basis.nrows() != rank or basis.ncols() != rank:
        raise ValueError("integral reduction must use square matrices of the Gram dimension")
    if reduced.nrows() != rank or reduced.ncols() != rank:
        raise ValueError("reduced Gram matrix has the wrong dimension")
    _source_require_integral_matrix(gram, "original Gram matrix")
    _source_require_integral_matrix(reduced, "reduced Gram matrix")
    _source_require_integral_matrix(basis, "integral reduction basis")
    if abs(basis.det()) != 1:
        raise ArithmeticError("integral reduction basis is not unimodular")
    if reduced != basis * gram * basis.transpose():
        raise ArithmeticError("integral reduction basis does not certify its reduced Gram matrix")


def _source_integral_vector_allmeth(
    gram: Matrix_integer_dense,
    critical_norm: Rational,
    strict_inequality: bool,
    *,
    reduction: Callable[[Matrix_integer_dense], tuple[Matrix_integer_dense, Matrix_integer_dense]] = _indefinite_reduction,
    isotropic_reduction: Callable[..., tuple[Matrix_integer_dense, Matrix_integer_dense, object]] = _compute_indefinite_reduction,
    find_vector: Callable[[Matrix_integer_dense, Matrix_integer_dense, Matrix_integer_dense, Rational, bool], tuple[Integer, ...]] = _source_integral_vector_allmeth_v2,
) -> tuple[Integer, ...]:
    r"""Port ``GetIntegralVector_allmeth`` from Positivity.h:410-482.

    The isotropic shortcut is enabled exactly when the requested inequality
    admits norm zero. Otherwise run IndefiniteReduction, search its reduced
    form with V2, and reconstruct the certified original-coordinate vector.
    """
    if gram.nrows() == 0 or gram.nrows() != gram.ncols():
        raise ValueError("integral-vector search requires a nonempty square Gram matrix")
    if gram != gram.transpose():
        raise ValueError("integral-vector search requires a symmetric Gram matrix")
    critical_norm = SageQQ(critical_norm)
    isotropic_fine = critical_norm < 0 if strict_inequality else critical_norm <= 0
    reduced_result = None
    if isotropic_fine:
        basis, reduced, isotropic = isotropic_reduction(gram, look_for_isotropic=True)
        if isotropic is not None:
            integral = vector(SageZZ, _source_remove_fraction_vector(isotropic))
            if not any(integral) or (integral * gram * integral.column())[0] != 0:
                raise ArithmeticError("isotropic reduction did not return a nonzero integral isotropic vector")
            return tuple(integral)
        # ComputeReductionIndefinite already returned the exact reduced pair
        # used by the source fallback; do not run IndefiniteReduction again.
        reduced_result = basis, reduced
    if reduced_result is None:
        basis, reduced = reduction(gram)
    else:
        basis, reduced = reduced_result
    _source_validate_integral_reduction(gram, basis, reduced)
    return find_vector(gram, basis, reduced, critical_norm, strict_inequality)


def _source_certify_positive_integral_vector(
    gram: Matrix_integer_dense, candidate: tuple[Integer, ...]
) -> tuple[Integer, ...]:
    r"""Enforce the public positive-vector contract at the consumer boundary."""
    if gram.nrows() != gram.ncols() or len(candidate) != gram.nrows():
        raise ValueError("positive-vector witness dimension must match the Gram matrix")
    try:
        integral_entries = tuple(SageZZ(entry) for entry in candidate)
    except (TypeError, ValueError, ArithmeticError) as error:
        raise ValueError("positive-vector witness coordinates must be integral") from error
    if any(entry != integer for entry, integer in zip(candidate, integral_entries)):
        raise ValueError("positive-vector witness coordinates must be integral")
    integral = integral_entries
    if not any(integral):
        raise ArithmeticError("positive-vector search returned the zero vector")
    witness = vector(SageZZ, integral)
    norm = SageZZ((witness * gram * witness.column())[0])
    if norm <= 0:
        raise ArithmeticError("positive-vector search returned a nonpositive-norm witness")
    return integral


def _integral_positive_vector_allmeth(gram: Matrix_integer_dense) -> tuple[Integer, ...]:
    r"""Port ``GetIntegralPositiveVector_allmeth`` from Positivity.h:499-518."""
    if gram.nrows() == 0:
        raise ValueError("the Gram matrix has no positive direction")
    return _source_certify_positive_integral_vector(
        gram, _source_integral_vector_allmeth(gram, SageQQ.zero(), True)
    )


def _source_graver_k_basis(dimension: int, degree: int) -> tuple[tuple[int, ...], ...]:
    r"""Port ``GetGraverKbasis`` from GraverBasis.h for the signed support basis.

    Enumerate supports in increasing size, lexicographic coordinate order, and
    then the source's binary sign-pattern order (−1 before +1). Each vector has
    entries in ``{-1,0,1}`` and support size at most ``degree``.
    """
    if dimension < 1 or degree < 1:
        raise ValueError("Graver K-basis needs positive dimension and degree")
    result: list[tuple[int, ...]] = []
    for support_size in range(1, min(degree, dimension) + 1):
        for support in combinations(range(dimension), support_size):
            for signs in product((-1, 1), repeat=support_size):
                row = [0] * dimension
                for coordinate, sign in zip(support, signs, strict=True):
                    row[coordinate] = sign
                result.append(tuple(row))
    return tuple(result)


def _source_random_matrix_perturbation(
    dimension: int,
    *,
    choose: Callable[[tuple[int, ...]], int] = choice,
    shuffle_values: Callable[[list[int]], None] = shuffle,
    choose_index: Callable[[int], int] = randrange,
) -> Matrix_integer_dense:
    r"""Port ``GetRandomMatrixPerturbation`` from Positivity.h:536-559."""
    if dimension < 1:
        raise ValueError("a random matrix perturbation needs positive dimension")
    if choose((0, 1)) == 0:
        permutation = list(range(dimension))
        shuffle_values(permutation)
        signs = tuple(choose((-1, 1)) for _ in range(dimension))
        return _source_permutation_sign_matrix(permutation, signs).change_ring(SageZZ)
    row = choose_index(dimension)
    column = choose_index(dimension)
    sign = choose((-1, 1))
    return _source_elementary_perturbation(dimension, row, column, sign).change_ring(SageZZ)


def _source_direct_graver_improvement(
    gram: Matrix_integer_dense,
    initial: FreeModuleElement[Integer],
    directions: tuple[FreeModuleElement[Integer], ...],
) -> FreeModuleElement[Integer]:
    r"""Port ``DirectImprovement`` from Positivity.h:588-603.

    For each source-ordered Graver direction, apply the maximal consecutive
    positive-norm descent prefix returned by ``GetAlpha``; repeat full passes
    until an entire pass makes no change.
    """
    if gram.nrows() != gram.ncols() or initial not in SageZZ**gram.nrows():
        raise ValueError("direct Graver improvement requires a matching integral Gram vector")
    if any(direction not in SageZZ**gram.nrows() for direction in directions):
        raise ValueError("every Graver direction must have the Gram matrix dimension")
    improved = vector(SageZZ, initial)
    while True:
        changes = 0
        for direction in directions:
            base_norm = int((improved * gram * improved.column())[0])
            mixed = int((improved * gram * direction.column())[0])
            direction_norm = int((direction * gram * direction.column())[0])
            alpha = _quadratic_descent_prefix(base_norm, mixed, direction_norm)
            if alpha:
                previous_norm = base_norm
                for step in range(1, alpha + 1):
                    step_norm = base_norm + 2 * mixed * step + direction_norm * step * step
                    if step_norm <= 0 or step_norm >= previous_norm:
                        raise ArithmeticError("Graver descent prefix contains a nonpositive or nondecreasing step")
                    previous_norm = step_norm
                improved += alpha * direction
                actual_norm = int((improved * gram * improved.column())[0])
                if actual_norm != previous_norm:
                    raise ArithmeticError("Graver descent prefix disagrees with the exact updated vector norm")
                changes += alpha
        if changes == 0:
            return improved


def _source_should_reset_positive_perturbation(perturbed_l1: int, original_l1: int) -> bool:
    r"""Apply Positivity.h:660's strict 1000-times L1 reset threshold."""
    if perturbed_l1 < 0 or original_l1 < 0:
        raise ValueError("entrywise Gram L1 norms must be nonnegative")
    return perturbed_l1 > 1000 * original_l1


def _source_positive_search_best(
    gram: Matrix_integer_dense,
    candidate: FreeModuleElement[Integer],
    best_vector: FreeModuleElement[Integer] | None,
    best_norm: Integer | None,
    without_improvement: int,
) -> tuple[FreeModuleElement[Integer] | None, Integer | None, int, bool]:
    r"""Apply the best-vector and termination transition from Positivity.h:641-657."""
    if gram.nrows() == 0 or gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("positive-search state requires a nonempty symmetric Gram matrix")
    if candidate not in SageZZ**gram.nrows():
        raise ValueError("positive-search candidate and Gram dimensions must agree")
    # The source increments this counter and returns immediately at 500, so
    # 500 is a terminal state, never an input to another transition.
    if not 0 <= without_improvement < 500:
        raise ValueError("the no-improvement counter must lie in the nonterminal range 0..499")
    if (best_vector is None) != (best_norm is None):
        raise ValueError("the best positive vector and its norm must be present together")
    if best_vector is not None:
        if best_vector not in SageZZ**gram.nrows():
            raise ValueError("the best positive vector has the wrong dimension")
        certified_best_norm = SageZZ((best_vector * gram * best_vector.column())[0])
        if best_norm != certified_best_norm or certified_best_norm <= 0:
            raise ArithmeticError("the saved best vector and positive norm are inconsistent")
    norm = SageZZ((candidate * gram * candidate.column())[0])
    if norm <= 0:
        raise ValueError("the specialized positive search requires a positive-norm candidate")
    if best_norm is None or norm < best_norm:
        best_vector = vector(SageZZ, candidate)
        best_norm = norm
        without_improvement = 0
    else:
        without_improvement += 1
    stop = without_improvement == 500 or best_norm < 10
    return best_vector, best_norm, without_improvement, stop


def _source_indefinite_short_positive_vector(
    gram: Matrix_integer_dense,
    find_positive: Callable[[Matrix_integer_dense], tuple[Integer, ...]],
    perturbation_factory: Callable[[int], Matrix_integer_dense] = _source_random_matrix_perturbation,
) -> tuple[Integer, ...]:
    r"""Port the complete ``INDEFINITE_GetShortPositiveVector`` loop.

    The supplied finder is precisely ``GetIntegralPositiveVector_allmeth``;
    injection exists only to permit deterministic, bounded source-operation
    tests without importing the entire research category stack.
    """
    if gram.nrows() == 0 or gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("short-positive-vector search requires a nonempty symmetric Gram matrix")
    rank = int(gram.nrows())
    directions = tuple(vector(SageZZ, row) for row in _source_graver_k_basis(rank, 2))
    gram_l1 = sum(abs(entry) for entry in gram.list())
    perturbation = identity_matrix(SageZZ, rank)
    best_vector: FreeModuleElement[Integer] | None = None
    best_norm: Integer | None = None
    without_improvement = 0
    while True:
        if perturbation.nrows() != rank or perturbation.ncols() != rank:
            raise ArithmeticError("positive-vector perturbation has the wrong dimension")
        try:
            _source_require_integral_matrix(perturbation, "positive-vector perturbation")
            integral_perturbation = matrix(SageZZ, perturbation)
        except (TypeError, ValueError, ArithmeticError) as error:
            raise ArithmeticError(
                "positive-vector perturbation is not an exact integral lattice automorphism: "
                "expected an exact integral matrix"
            ) from error
        if integral_perturbation != perturbation or abs(integral_perturbation.det()) != 1:
            raise ArithmeticError("positive-vector perturbation is not an integral lattice automorphism")
        perturbation = integral_perturbation
        perturbed_gram = perturbation * gram * perturbation.transpose()
        if perturbed_gram != perturbed_gram.transpose():
            raise ArithmeticError("positive-vector perturbation lost Gram symmetry")
        # Positivity.h:614-624 takes the *integral* positive witness of the
        # perturbed Gram matrix from GetIntegralPositiveVector_allmeth.  Its
        # V2 approximation iterators are implemented above.  Calling the
        # preamble's positive_vector instead silently substituted its
        # different exact diagonalization algorithm for this source step.
        found_raw = tuple(find_positive(perturbed_gram))
        try:
            found = vector(SageZZ, found_raw)
        except (TypeError, ValueError, ArithmeticError) as error:
            raise ValueError("positive-vector search witness coordinates must be integral") from error
        if any(entry != integer for entry, integer in zip(found_raw, found)):
            raise ArithmeticError("positive-vector search returned nonintegral coordinates")
        if len(found_raw) != rank or not any(found):
            raise ArithmeticError("positive-vector search returned a zero or dimensionally invalid witness")
        found_norm = SageZZ((found * perturbed_gram * found.column())[0])
        if found_norm <= 0:
            raise ArithmeticError("positive-vector search returned a nonpositive-norm witness")
        current = perturbation.transpose() * found.column()
        current_norm = SageZZ((current.transpose() * gram * current)[0, 0])
        if current_norm != found_norm:
            raise ArithmeticError("positive-vector perturbation changed the exact witness norm")
        improved = _source_direct_graver_improvement(gram, current.column(0), directions)
        improved_norm = SageZZ((improved * gram * improved.column())[0])
        if improved_norm <= 0 or improved_norm > current_norm:
            raise ArithmeticError("direct Graver improvement did not preserve positivity and norm descent")
        best_vector, best_norm, without_improvement, stop = _source_positive_search_best(
            gram, improved, best_vector, best_norm, without_improvement
        )
        if stop:
            if best_vector is None:
                raise ArithmeticError("the positive-vector search has no best vector")
            return tuple(SageZZ(entry) for entry in best_vector)
        # Positivity.h:660 resets at a strict factor of 1000.
        perturbed_l1 = sum(abs(entry) for entry in perturbed_gram.list())
        if _source_should_reset_positive_perturbation(perturbed_l1, gram_l1):
            perturbation = identity_matrix(SageZZ, rank)
        perturbation *= perturbation_factory(rank)



def _positive_direction(lattice: Lattices.ParentMethods) -> tuple[Rational, ...]:
    r"""Pinned polyhedral_common@a55fcb7b ``INDEFINITE_GetShortPositiveVector``."""
    gram = _gram_matrix(lattice)
    return tuple(SageQQ(entry) for entry in _source_indefinite_short_positive_vector(gram, _integral_positive_vector_allmeth))


def _primitive_integral_direction(
    coordinates: tuple[Rational, ...], expected_rank: int | None = None
) -> tuple[Rational, tuple[Integer, ...]]:
    r"""Return ``(TheMult, TheVect)`` as in ``RemoveFractionPlusCoef``.

    For nonzero rational coordinates ``x``, the returned positive multiplier
    ``m`` and primitive integral tuple ``v`` satisfy ``x = v / m`` exactly.
    This is the source contract used to scale ``MaxScal`` before enumeration.
    """
    if expected_rank is not None and len(coordinates) != expected_rank:
        raise ValueError("the rational direction must have the lattice rank")
    if not coordinates:
        raise ValueError("the rational direction must have positive dimension")
    if not any(coordinates):
        raise ValueError("the rational direction must be nonzero")
    denominator = SageZZ.one()
    for entry in coordinates:
        denominator = denominator.lcm(SageZZ(SageQQ(entry).denominator()))
    scaled_coordinates = tuple(SageZZ(SageQQ(entry) * denominator) for entry in coordinates)
    content = SageZZ(gcd(*(int(entry) for entry in scaled_coordinates)))
    if content <= 0:
        raise ArithmeticError("a nonzero rational direction has positive integral content")
    primitive_coordinates = tuple(SageZZ(entry // content) for entry in scaled_coordinates)
    return SageQQ(denominator) / SageQQ(content), primitive_coordinates


def _integral_direction_pairing_bound(max_scal: Rational, scale: Rational) -> Rational:
    r"""Convert a pairing bound on rational ``x=v/scale`` to integral ``v``.

    The source's ``MaxScal`` bounds ``b(candidate, x)``; the enumerator checks
    ``b(candidate, v)`` for the primitive integral representative. Since
    ``v=scale*x``, this bound scales linearly, not quadratically. The source
    then increments the integral-coordinate bound by one per failed search.
    """
    return SageQQ(max_scal) * SageQQ(scale)


def _next_integral_pairing_bound(bound: Rational) -> Rational:
    r"""Advance the source's integral ``MaxScal`` by exactly one."""
    return SageQQ(bound) + 1


def _pinned_positive_gcd_dot(values: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    r"""Reproduce ``ComputeGcdDot`` followed by ``PositivityNormalizeGcdDot``.

    The returned coefficients satisfy ``dot(values, coefficients) = gcd(values)``
    with positive gcd. Their exact choice fixes the affine CVP center in the
    pinned perfect-domain consumer.
    """
    if not values:
        raise ValueError("the gcd-dot input must be nonempty")
    divisor = int(values[0])
    coefficients = [1] + [0] * (len(values) - 1)
    for position in range(1, len(values)):
        other = int(values[position])
        if divisor == 0 and other == 0:
            a, b, next_divisor = 0, 0, 0
        else:
            f = abs(divisor)
            fm = 1 if divisor >= 0 else -1
            g = abs(other)
            gm = 0
            while g:
                quotient = f // g
                f, fm, g, gm = g, gm, f - quotient * g, fm - quotient * gm
            next_divisor = f
            a = fm
            b = 0 if other == 0 else (f - fm * divisor) // other
        divisor = next_divisor
        coefficients = [a * coefficient for coefficient in coefficients]
        coefficients[position] = b
    if divisor <= 0:
        # The pinned PositivityNormalizeGcdDot returns -V when gcd is zero
        # as well as when it is negative; this matters for a singleton zero.
        return -divisor, tuple(-coefficient for coefficient in coefficients)
    return divisor, tuple(coefficients)


def _dual_lll_enumeration_key(
    coordinates: tuple[int, ...], transform: Matrix_integer_dense
) -> tuple[Integer, ...]:
    """Order a CVP result by the pinned recursive coordinates ``P*x``.

    ``CVPSolver`` reconstructs the original column as ``x=P^{-1}(y-s)`` for
    the dual-LLL reduction, where ``s`` is the fixed ``ReductionMod1vector``
    translation for this target. Thus ``y=P*x+s``; the same affine shift is
    added to every candidate in one CVP call and cannot change their order.
    The recursive enumerator's last coordinate is the outer loop and zero the
    inner loop, so sorting by reversed ``P*x`` coordinates reproduces it.
    """
    reduced = transform * vector(SageZZ, coordinates).column()
    return tuple(reduced[index, 0] for index in reversed(range(transform.nrows())))


def _pinned_sublattice_basis_reduction(basis: Matrix_integer_dense) -> Matrix_integer_dense:
    r"""Reproduce ``SublatticeBasisReductionKernel`` from pinned ClassicLLL.h.

    Alternate direct and dual LLL on the row Gram matrix, accepting a
    unimodular row transform only when the ambient basis's entrywise L1 norm
    strictly decreases. This stopping rule fixes the basis used by CVP.
    """
    reduced_basis = matrix(SageZZ, basis)
    norm = sum(abs(entry) for entry in reduced_basis.list())
    while reduced_basis.nrows() > 1:
        changed = False
        for dual in (False, True):
            basis_gram = reduced_basis * reduced_basis.transpose()
            gram = basis_gram.adjugate() if dual else basis_gram
            reduction = libgap.LLLReducedGramMat(
                [[int(gram[i, j]) for j in range(gram.ncols())] for i in range(gram.nrows())],
                libgap(99) / libgap(100),
            )
            transform = matrix(
                SageZZ,
                [[SageZZ(int(entry)) for entry in row] for row in reduction["transformation"]],
            )
            reduced_gram = matrix(
                SageZZ,
                [[SageZZ(int(entry)) for entry in row] for row in reduction["remainder"]],
            )
            if abs(int(transform.det())) != 1:
                raise ArithmeticError("the pinned sublattice LLL transform is not unimodular")
            if transform * gram * transform.transpose() != reduced_gram:
                raise ArithmeticError("the pinned sublattice LLL transform misses its reduced Gram")
            if dual:
                # The GAP call above reduces adj(G) directly and returns P
                # with P*adj(G)*P^T=R. ClassicLLL.h's dual basis transform is
                # Q=P^{-T}, giving Q*G*Q^T=det(G)*R^{-1}. Preserve Q over ZZ.
                transform = matrix(SageZZ, transform.inverse().transpose())
                expected_candidate_gram = basis_gram.det() * reduced_gram.inverse()
            else:
                expected_candidate_gram = reduced_gram
            candidate = transform * reduced_basis
            if candidate * candidate.transpose() != expected_candidate_gram:
                raise ArithmeticError("the pinned sublattice transform misses its exact reduced Gram identity")
            candidate_norm = sum(abs(entry) for entry in candidate.list())
            if candidate_norm < norm:
                reduced_basis = candidate
                norm = candidate_norm
                changed = True
        if not changed:
            break
    return reduced_basis


def _pinned_integer_covector_kernel(covector: tuple[int, ...]) -> Matrix_integer_dense:
    r"""Reproduce ``basic_common_cpp::NullspaceIntVect`` for an integral row.

    The pinned implementation clears the first nonzero column, then builds
    its kernel from the recursive ``ComputeGCD_information`` unimodular
    transform. The transform matters here because its row order feeds the
    subsequent LLL and CVP enumeration.
    """
    rank = len(covector)
    if rank == 0:
        return matrix(SageZZ, 0, 0)
    pivot = next((i for i, value in enumerate(covector) if value), None)
    if pivot is None:
        return identity_matrix(SageZZ, rank)

    def pair_gcd(m: int, n: int) -> Matrix_integer_dense:
        if m == 0 and n == 0:
            return identity_matrix(SageZZ, 2)
        f = abs(m)
        fm = 1 if m >= 0 else -1
        g = abs(n)
        gm = 0
        while g:
            quotient = f // g
            f, fm, g, gm = g, gm, f - quotient * g, fm - quotient * gm
        coefficient_m = fm
        coefficient_n = 0 if n == 0 else (f - fm * m) // n
        return matrix(SageZZ, [[coefficient_m, -n // f], [coefficient_n, m // f]])

    def gcd_transform(values: tuple[int, ...]) -> Matrix_integer_dense:
        size = len(values)
        if size == 1:
            return identity_matrix(SageZZ, 1)
        if size == 2:
            return pair_gcd(values[0], values[1])
        previous = gcd_transform(values[:-1])
        pair = pair_gcd(gcd(*values[:-1]), values[-1])
        result = matrix(SageZZ, size, size)
        for i in range(size - 1):
            result[i, 0] = pair[0, 0] * previous[i, 0]
        result[size - 1, 0] = pair[1, 0]
        for i in range(size - 2):
            for j in range(size - 1):
                result[j, i + 1] = previous[j, i + 1]
        for i in range(size - 1):
            result[i, size - 1] = pair[0, 1] * previous[i, 0]
        result[size - 1, size - 1] = pair[1, 1]
        return result

    nonpivot = tuple(i for i in range(rank) if i != pivot)
    values = (int(covector[pivot]),) + tuple(int(covector[i]) for i in nonpivot)
    transform = gcd_transform(values)
    values_column = vector(SageZZ, values).column()
    transformed_values = transform.transpose() * values_column
    # ComputeGCD_information deliberately leaves a singleton unchanged,
    # including its sign; for two or more entries ComputePairGcd makes the
    # gcd nonnegative. A rank-one covector therefore has this source-specific
    # signed singleton contract even though its kernel is empty.
    expected_gcd = values[0] if rank == 1 else gcd(*values)
    if transformed_values[0, 0] != expected_gcd or any(transformed_values[i, 0] for i in range(1, rank)):
        raise ArithmeticError("pinned recursive GCD transform does not isolate the exact covector gcd")
    if abs(int(transform.det())) != 1:
        raise ArithmeticError("pinned recursive GCD transform is not unimodular")
    result = matrix(SageZZ, rank - 1, rank)
    for i in range(rank - 1):
        result[i, pivot] = transform[0, i + 1]
        for j, column in enumerate(nonpivot):
            result[i, column] = transform[j + 1, i + 1]
    if any(entry for entry in (result * vector(SageZZ, covector).column()).list()):
        raise ArithmeticError("pinned integer covector-kernel construction failed its annihilation identity")
    return result


class _PositiveVectorEnumerator:
    r"""Vectors v of L with q(v) >= 0 (total) or q(v) = 0 (isotropic) and 0 < b(v, t) <= max, by multiplier.

    With t primitive timelike, d = div(t) and h a Bezout partner (b(t, h) = d), a vector with
    b(v, t) = m d is v = m h + x for x in K = t^perp, and q(v) = m^2 d^2 / q(t) + q(x - m w) for
    w = (d / q(t)) t - h in K tensor Q. K is negative definite for q; the closest-vector
    backend uses the positive-definite Gram matrix -q|K, so admissible x are the close vectors of
    m w within m^2 d^2 / q(t) (on the sphere, in isotropic mode).
    """

    def __init__(self, lattice: Lattices.ParentMethods, rational_direction: tuple[Rational, ...], max_scal: Rational, allowed_squares: RealSet) -> None:
        self.lattice = lattice
        self.gram = _gram_matrix(lattice)
        self.allowed_squares = allowed_squares
        scale, timelike_coordinates = _primitive_integral_direction(
            tuple(SageQQ(entry) for entry in rational_direction),
            expected_rank=int(lattice.module_rank()),
        )
        self.timelike_row = vector(SageZZ, timelike_coordinates)
        self.timelike_covector = self.gram * self.timelike_row.column()
        self.scale = scale
        self.scaled_max = _integral_direction_pairing_bound(max_scal, scale)
        square = SageQQ((self.timelike_row * self.gram * self.timelike_row.column())[0])
        if square <= 0:
            raise ValueError("positive-vector enumeration needs a timelike direction")
        # div(t) is gcd(b(t,e_i)) over the integral basis.  The exact
        # Gram covector was just computed above, so use its coordinates
        # directly instead of re-entering the generic categorical pairing
        # iterator for the same finite gcd.
        pairing_coordinates = tuple(int(entry) for entry in self.timelike_covector.column(0))
        self.d = SageQQ(gcd(*pairing_coordinates))
        fraction_map = lattice.base_ring().fraction_field_map()
        self.rational_field = fraction_map.codomain()
        # Match the pinned ComputeGcdDot/PositivityNormalizeGcdDot coefficient
        # recurrence exactly: a different valid Bezout partner changes the
        # affine CVP center and therefore the first candidate used by flipping.
        divisor, bezout = _pinned_positive_gcd_dot(pairing_coordinates)
        if divisor != int(self.d) or sum(a * b for a, b in zip(pairing_coordinates, bezout, strict=True)) != int(self.d):
            raise ArithmeticError("the pinned gcd-dot partner does not realize the direction divisibility")
        partner_row = vector(SageZZ, bezout)
        # lorentzian_perfect.h calls NullspaceIntVect, not GAP's
        # NullspaceIntMat. Reproduce the pinned C++ Euclidean-gcd transform:
        # the resulting saturated basis has a source-defined row order that
        # determines subsequent LLL and recursive CVP enumeration order.
        kernel_rows = _pinned_integer_covector_kernel(pairing_coordinates)
        if kernel_rows.nrows() != int(lattice.module_rank()) - 1 or any(
            entry for entry in (kernel_rows * vector(SageZZ, pairing_coordinates).column()).list()
        ):
            raise ArithmeticError("the integral orthogonal-kernel basis does not span the timelike covector kernel")
        # A rank-correct integral kernel can still span a proper finite-index
        # sublattice. The pinned NullspaceIntMat/NullspaceIntVect contract is a
        # saturated Z-kernel; check its Smith factors before any basis reduction.
        kernel_smith = kernel_rows.smith_form()[0]
        if any(kernel_smith[i, i] != 1 for i in range(kernel_rows.nrows())):
            raise ArithmeticError("the integral orthogonal-kernel basis is not saturated")
        # The pinned C++ consumer immediately applies SublatticeBasisReduction
        # to NullspaceIntVect's ordered kernel basis. Its alternating direct/
        # dual LLL and strict ambient-L1 stopping rule control the subsequent CVP.
        kernel_rows = _pinned_sublattice_basis_reduction(kernel_rows)
        self.inclusion_rows = kernel_rows
        # Lorentzian.g:337-343 defines GramMat = -Ubasis*LorMat*Ubasis^T,
        # the positive-definite form used by ClosestAtDistanceVallentinProgram.
        kernel_gram = -(self.inclusion_rows * self.gram * self.inclusion_rows.transpose())
        kernel_gram_integer = matrix(
            SageZZ,
            kernel_gram.nrows(),
            kernel_gram.ncols(),
            [SageZZ(kernel_gram[i, j]) for i in range(kernel_gram.nrows()) for j in range(kernel_gram.ncols())],
        )
        # LATT_near's CVPSolver applies LLLreducedBasisDual. Its vendored
        # ClassicLLL.h implementation is adapted from GAP's
        # LLLReducedGramMat, which returns the exact integral row transform
        # as ``transformation``. Use that same algorithm on adj(G): Sage's
        # LLL_gram may choose a different reduced basis and therefore a
        # different first shortest vector.
        if kernel_gram_integer.nrows() > 1:
            adjugate_gram = kernel_gram_integer.adjugate()
            gap_gram = [
                [int(adjugate_gram[row, column]) for column in range(adjugate_gram.ncols())]
                for row in range(adjugate_gram.nrows())
            ]
            # ClassicLLL.h fixes the Lovasz parameter at y=99/100; GAP's
            # default is 3/4, so pass the pinned value explicitly.
            gap_reduction = libgap.LLLReducedGramMat(gap_gram, libgap(99) / libgap(100))
            gap_transform = gap_reduction["transformation"]
            direct_transform = matrix(
                SageZZ,
                [[SageZZ(int(entry)) for entry in row] for row in gap_transform],
            )
            if abs(int(direct_transform.det())) != 1:
                raise ArithmeticError("GAP's dual-LLL source transform is not unimodular")
            # CVPSolver calls LLLreducedBasisDual, not LLLreducedBasis on
            # adj(G) directly. ClassicLLL.h returns Q=P^{-T}; its reduced-Gram
            # identity uses Q G Q^T = det(G) R^{-1}, where R=P adj(G) P^T.
            # The CVP enumerator then reconstructs original columns x=P^{-1}y,
            # so its recursive coordinates are y=P*x. Keep Q for the Gram
            # identity and P separately for ordering the returned candidates.
            self.enumeration_transform = direct_transform.inverse().transpose()
            self.enumeration_order_transform = direct_transform
            reduced_adjugate = matrix(
                SageZZ,
                [[SageZZ(int(entry)) for entry in row] for row in gap_reduction["remainder"]],
            )
            if direct_transform * adjugate_gram * direct_transform.transpose() != reduced_adjugate:
                raise ArithmeticError("GAP's dual-LLL source transform does not reproduce its returned reduced adjugate")
            if self.enumeration_transform * kernel_gram_integer * self.enumeration_transform.transpose() != kernel_gram_integer.det() * reduced_adjugate.inverse():
                raise ArithmeticError("the dual-LLL enumeration transform does not reproduce the dual-reduced Gram")
        else:
            self.enumeration_transform = identity_matrix(SageZZ, kernel_gram_integer.nrows())
            self.enumeration_order_transform = self.enumeration_transform
        if self.enumeration_transform * self.enumeration_order_transform.transpose() != identity_matrix(
            SageZZ, kernel_gram_integer.nrows()
        ):
            raise ArithmeticError("dual-LLL CVP reconstruction and recursive ordering transforms are not inverse-transposes")
        self.cvp_engine = _ExactCVPEngine._from_positive_engine_gram(lattice.base_ring(), kernel_gram_integer)
        alpha = self.d / square
        timelike_row = self.timelike_row
        translation = vector(
            SageQQ,
            [alpha * timelike_row[position] - partner_row[position] for position in range(len(timelike_row))],
        ).column()
        target_coordinates = self.inclusion_rows.transpose().change_ring(SageQQ).solve_right(translation)
        if self.inclusion_rows.transpose().change_ring(SageQQ) * target_coordinates != translation:
            raise ArithmeticError("the Bezout translation is not in the rational orthogonal-kernel span")
        # Both pinned consumers have the same semantic CVP center. C++ forms
        # eTrans = alpha*eVect-eBasSol = w and passes its coordinates directly;
        # the older GAP wrapper forms -w and negates eSol at its call boundary.
        # Here target_coordinates solves Ubasis^T*w = alpha*t-h, so retain +w.
        self.base_target_coordinates = tuple(SageQQ(entry) for entry in target_coordinates.column(0))
        # The rational field's selected engine is Sage QQ for this integral
        # lattice leaf. Retain these private arithmetic coordinates locally;
        # only the resulting coefficients are raised back into the owned ring.
        self.base_target_engine_coordinates = tuple(SageQQ(coordinate) for coordinate in self.base_target_coordinates)
        self.base_bound = (self.d * self.d) / square
        self.partner_row = partner_row

    @cached_property
    def fixed_distance_ldl_factorization(self) -> tuple[list[list[Rational]], list[Rational]] | None:
        """Cache the exact LDL certificate for the fixed dual-LLL Gram."""
        gram = self.cvp_engine.engine_gram_q
        if gram.nrows() == 0:
            return None
        reduced_gram = self.enumeration_transform * gram * self.enumeration_transform.transpose()
        return _positive_ldl_factorization(reduced_gram)

    def max_multiplier(self) -> int | None:
        if self.scaled_max <= 0:
            return None
        return int(self.scaled_max / self.d)

    def first_shell(self) -> tuple[int, tuple[tuple[int, ...], ...]] | None:
        maximum = self.max_multiplier()
        if maximum is None:
            return None
        for multiplier in range(1, maximum + 1):
            vectors = self.vectors_at_multiplier(multiplier)
            if vectors:
                return multiplier, vectors
        return None

    @cached_method
    def vectors_at_multiplier(self, multiplier: int) -> tuple[tuple[int, ...], ...]:
        r"""The coordinate rows of the shell vectors at this multiplier; elements are raised only for a returned cell."""
        if multiplier <= 0:
            raise ValueError("the positive-vector shell multiplier must be positive")
        # Scale the cached exact coordinate vector directly. The generic
        # scalar-action morphism repeatedly flattens owned categories here;
        # this is the same module element as scalar_multiple, assembled once
        # through the owning module's linear-combination constructor.
        scalar = SageQQ(multiplier)
        target = tuple(
            _owned_engine_element(self.rational_field, scalar * coordinate)
            for coordinate in self.base_target_engine_coordinates
        )
        bound = _owned_rational(
            self.lattice, (SageQQ(multiplier) ** 2) * self.base_bound, field=self.rational_field
        )
        # Keep the same exact coordinate Gram, target, and squared-distance
        # bound as close_vectors. Isotropic mode uses the dedicated exact-sphere
        # recursion; total mode keeps PARI-backed bounded CVP and falls back to
        # exact LDL pruning only when PARI's stack overflows. Both stay at the
        # coordinate level without allocating owned lattice vectors or indexed
        # families for intermediate shell candidates.
        if self.allowed_squares == ZERO_SQUARES:
            # The pinned solver.fixed_dist_vectors branch enumerates the
            # exact sphere, not the full ellipsoid followed by a norm filter.
            # Keep its recursive coordinate order by applying the same dual-
            # LLL ordering key used for the at-most-distance branch.
            distances = _fixed_distance_dual_lll_solutions(
                self.cvp_engine.engine_gram_q,
                self.cvp_engine._engine_target_coordinates(target),
                self.cvp_engine._engine_scalar(bound),
                self.enumeration_transform,
                self.enumeration_order_transform,
                ldl_factorization=self.fixed_distance_ldl_factorization,
            )
        else:
            try:
                distances = self.cvp_engine.close_vector_coordinates(target, bound)
            except PariError as error:
                # Keep the pinned PARI CVP as the normal path. Its default
                # stack overflows on several published E10 vector-orbit cases;
                # only that specific failure falls back to exact bounded LDL
                # enumeration in the same dual-LLL frame and source order.
                if "stack overflows" not in str(error):
                    raise
                distances = _bounded_distance_dual_lll_solutions(
                    self.cvp_engine.engine_gram_q,
                    self.cvp_engine._engine_target_coordinates(target),
                    self.cvp_engine._engine_scalar(bound),
                    self.enumeration_transform,
                    self.enumeration_order_transform,
                    ldl_factorization=self.fixed_distance_ldl_factorization,
                )
            # CVPSolver's dual-LLL map gives recursive coordinates y=P*x;
            # the shared key follows computeIt_Gen_Kernel's outer-to-inner order.
            distances = tuple(sorted(distances, key=lambda item: _dual_lll_enumeration_key(item[0], self.enumeration_order_transform)))
        close_coordinates = tuple(coordinates for coordinates, _square in distances)
        close_squares = tuple(square for _coordinates, square in distances)
        return self._ambient_vectors(multiplier, close_coordinates, close_squares)

    def _ambient_vectors(
        self,
        multiplier: int,
        close_coordinates: tuple[tuple[Integer, ...], ...],
        close_squares: tuple[Rational, ...],
    ) -> tuple[tuple[int, ...], ...]:
        shift = multiplier * self.partner_row
        result: list[tuple[int, ...]] = []
        if len(close_squares) != len(close_coordinates):
            raise ValueError("each closest-vector coordinate must have exactly one squared distance")
        for index, close in enumerate(close_coordinates):
            distance = close_squares[index]
            if distance < 0 or distance > multiplier**2 * self.base_bound:
                raise ArithmeticError("a closest-vector result lies outside its exact squared-distance bound")
            if any(SageQQ(coordinate).denominator() != 1 for coordinate in close):
                raise ArithmeticError("a closest-vector result has nonintegral lattice coordinates")
            row = shift + vector(SageZZ, [int(coordinate) for coordinate in close]) * self.inclusion_rows
            candidate_square = int(row * self.gram * row)
            expected_square = multiplier**2 * self.base_bound - distance
            # Since the distance was checked against the full radius above,
            # expected_square is nonnegative. Equality therefore also rejects
            # negative-norm candidates without a redundant later norm branch.
            if candidate_square != expected_square:
                raise ArithmeticError(
                    "the reconstructed Lorentzian norm disagrees with its exact CVP squared distance"
                )
            # Lorentzian.g:377-385 verifies both the scalar-product shell
            # and norm after reconstructing each closest-vector solution.
            pairing = (row * self.timelike_covector)[0]
            if pairing != multiplier * self.d or pairing <= 0 or pairing > self.scaled_max:
                raise ArithmeticError("a close-vector shell returned a vector outside its prescribed pairing fibre")
            # Lorentzian.g:382-388: the isotropic option keeps fNorm = 0 among
            # the closest-vector solutions, which lie within the distance.
            if candidate_square not in self.allowed_squares:
                continue
            result.append(tuple(int(entry) for entry in row))
        # Keep the exact CVP backend's enumeration order. The pinned
        # LATT_near/FindAtMostDistVectors path orders through its LLL-reduced
        # basis and recursive coordinate bounds, not by norm; the first row
        # is consumed by LORENTZ_Kernel_Flipping, so norm sorting is not a
        # source-faithful substitute for that ordering.
        return tuple(result)


def _first_positive_shell(
    lattice: Lattices.ParentMethods, rational_direction: tuple[Rational, ...], max_scal: Rational, allowed_squares: RealSet
) -> tuple[tuple[int, ...], ...]:
    r"""The coordinate rows of the first nonempty shell with 0 < b(v, t) <= max_scal, or ()."""
    enumerator = _PositiveVectorEnumerator(lattice, rational_direction, max_scal, allowed_squares)
    if enumerator.scaled_max <= 0:
        raise ArithmeticError(f"the positive-vector shells below {max_scal} along {rational_direction} are empty by definition: the bound must be positive")
    first = enumerator.first_shell()
    return () if first is None else first[1]


def _search_initial_vectors(lattice: Lattices.ParentMethods, direction: tuple[Rational, ...], allowed_squares: RealSet) -> tuple[tuple[int, ...], ...]:
    r"""The first shell along the direction, following LORENTZ_SearchInitialVectorOpt's successive integral bounds."""
    if allowed_squares == ZERO_SQUARES and not lattice.is_isotropic():
        raise ValueError(f"{lattice} is anisotropic: it has no isotropic perfect cells")
    gram = _gram_matrix(lattice)
    _scale, primitive_direction = _primitive_integral_direction(
        tuple(SageQQ(entry) for entry in direction),
        expected_rank=int(lattice.module_rank()),
    )
    # The source initializes MaxScal to q(t) for its integral direction t.
    # Compute that norm directly on the exact primitive integral row instead
    # of recovering it from q(direction) and the fraction-clearing multiplier.
    primitive_row = vector(SageZZ, primitive_direction)
    max_scal = (primitive_row * gram * primitive_row.column())[0]
    enumerator = _PositiveVectorEnumerator(lattice, primitive_direction, max_scal, allowed_squares)
    if enumerator.scaled_max <= 0:
        raise ArithmeticError("the initial positive-vector shell bound must be positive")
    # Each failed bound has already exhausted every multiplier up to its
    # maximum.  Increasing the bound cannot change those exact shells;
    # resume at the first previously unexamined multiplier instead of
    # repeatedly scanning the same empty prefix.
    next_multiplier = 1
    while True:
        maximum = enumerator.max_multiplier()
        if maximum is None:
            raise ArithmeticError("the initial positive-vector shell bound became nonpositive")
        for multiplier in range(next_multiplier, maximum + 1):
            rows = enumerator.vectors_at_multiplier(multiplier)
            if rows:
                return rows
        next_multiplier = max(next_multiplier, maximum + 1)
        # Lorentzian.g:411–423 increases the integral MaxScal by one, rather
        # than exponentially expanding the search interval. The enumerator's
        # scaled bound is the sole state used to compute the next multiplier;
        # do not maintain a second, now-unused unscaled bound here.
        enumerator.scaled_max = _next_integral_pairing_bound(enumerator.scaled_max)


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


def _quadratic_descent_prefix(base_norm: int, mixed: int, direction_norm: int) -> int:
    r"""Return GetAlpha's exact accepted prefix for ``q(k)=A+2Bk+Ck^2``.

    The pinned loop accepts consecutive steps while ``0 < q(k) < q(k-1)``.
    Its first step can fail even if a later step would decrease, so test that
    step before bracketing. Once it succeeds, the decreasing-positive steps
    form a prefix: the difference ``q(k)-q(k-1)`` is affine in ``k`` and the
    quadratic cannot become positive again after crossing zero while still
    decreasing. Exponential bracketing and integer bisection therefore return
    the same endpoint as the pinned linear loop.
    """

    def accepted_step(step: int) -> bool:
        candidate_norm = base_norm + 2 * mixed * step + direction_norm * step * step
        return candidate_norm > 0 and 2 * mixed + (2 * step - 1) * direction_norm < 0

    if not accepted_step(1):
        return 0
    alpha = 1
    upper = 2
    while accepted_step(upper):
        alpha = upper
        upper *= 2
    lower = alpha
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if accepted_step(middle):
            lower = middle
        else:
            upper = middle
    return lower


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
    basis = matrix(
        SageQQ,
        2,
        gram.nrows(),
        [*base_vector.column(0), *direction_vector.column(0)],
    )
    restricted = basis * gram * basis.transpose()
    a, b, c = restricted[0, 0], restricted[0, 1], restricted[1, 1]
    discriminant = b * b - a * c
    # Lorentzian.g:721-728 (GetRationalIsotropyVectors): a negative
    # discriminant contradicts signature (1,1); an irrational root gives none.
    if discriminant < 0:
        raise ArithmeticError("the flip plane must have signature (1,1)")
    isotropic_bound: Rational | None = None
    if SageQQ(discriminant).is_square():
        root = SageQQ(discriminant).sqrt()
        candidates = ((-b + root) / c, (-b - root) / c) if c != 0 else ((-a / (2 * b),) if b != 0 else ())
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


def _pinned_nullspace_matrix(gram: Matrix_rational_dense) -> Matrix_rational_dense:
    r"""Return the pinned GAP ``NullspaceMat`` row basis (left nullspace).

    GAP's returned vectors have length ``gram.nrows()``, including when the
    input is rectangular. The order and normalization feed source Gram--Schmidt
    and first-negative-vector search; Sage's canonical basis can change them.
    """
    if gram.nrows() == 0:
        # NullspaceMat returns left-kernel vectors of length nrows, so the
        # left kernel of a 0-by-n matrix has the empty 0-by-0 row basis.
        return matrix(SageQQ, 0, 0)
    gap_rows = [[libgap(SageQQ(entry)) for entry in row] for row in gram.rows()]
    gap_basis = libgap.NullspaceMat(gap_rows)
    rows = [[SageQQ(entry) for entry in row] for row in gap_basis]
    return matrix(SageQQ, rows) if rows else matrix(SageQQ, 0, gram.nrows())


def _pinned_independent_column_indices(matrix_work: Matrix_rational_dense) -> tuple[int, ...]:
    r"""Greedy independent-column selection used by GAP ``ColumnReduction``."""
    selected: list[int] = []
    rank = 0
    for column in range(int(matrix_work.ncols())):
        candidate = matrix_work.matrix_from_columns([*selected, column])
        candidate_rank = int(candidate.rank())
        if candidate_rank > rank:
            selected.append(column)
            rank = candidate_rank
    return tuple(selected)


def _pinned_nearest_integer(value: Rational) -> Integer:
    r"""GAP ``NearestInteger``: round exact half-integers toward +infinity."""
    return SageZZ((SageQQ(value) + SageQQ(1) / 2).floor())


def _source_some_negative_vector(gram: Matrix_rational_dense) -> tuple[Integer, ...]:
    r"""Port ``GetSomeNegativeVector`` from PositivitySymmetricMatrices.g:173-190.

    Diagonalize the rational symmetric form by the source's radical split and
    orthogonal Gram--Schmidt procedure. Normalize each nonpositive diagonal
    direction by its largest coordinate, then apply the source's increasing
    integer-multiplier/nearest-integer search in diagonal-direction order.
    """
    if gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("negative-vector search requires a symmetric Gram matrix")
    rank = int(gram.nrows())
    if not rank:
        raise ValueError("the zero-dimensional form has no negative vector")
    # Preserve the pinned GAP NullspaceMat basis order: its row normalization
    # controls the source diagonalizer and the subsequent negative-vector scan.
    radical = _pinned_nullspace_matrix(gram)
    radical_rank = int(radical.nrows())
    pivot_columns = set(_pinned_independent_column_indices(radical))
    transform = cast(Matrix_rational_dense, matrix(cast(Parent, SageQQ), rank, rank, 0))
    for i in range(radical_rank):
        for j in range(rank):
            transform[i, j] = radical[i, j]
    row_index = radical_rank
    for column_index in range(rank):
        if column_index not in pivot_columns:
            transform[row_index, column_index] = 1
            row_index += 1
    reduced = cast(Matrix_rational_dense, transform * gram * transform.transpose())
    if any(reduced[i, j] for i in range(radical_rank) for j in range(rank)):
        raise ArithmeticError("the radical basis does not split orthogonally")
    dimension = rank - radical_rank
    nondegenerate = cast(
        Matrix_rational_dense,
        reduced.matrix_from_rows_and_columns(range(radical_rank, rank), range(radical_rank, rank)),
    )
    orthogonal_basis: list[FreeModuleElement[Rational]] = []
    for _ in range(dimension):
        if not orthogonal_basis:
            complement = identity_matrix(cast(Parent, SageQQ), dimension)
        else:
            previous = matrix(cast(Parent, SageQQ), orthogonal_basis)
            complement = _pinned_nullspace_matrix(nondegenerate * previous.transpose())
        restricted = cast(Matrix_rational_dense, complement * nondegenerate * complement.transpose())
        selected: list[Rational] | None = None
        for i in range(int(restricted.nrows())):
            if restricted[i, i]:
                selected = [SageQQ.one() if j == i else SageQQ.zero() for j in range(int(restricted.nrows()))]
                break
        if selected is None:
            for i in range(int(restricted.nrows())):
                for j in range(i + 1, int(restricted.ncols())):
                    if restricted[i, j]:
                        selected = [SageQQ.one() if k in (i, j) else SageQQ.zero() for k in range(int(restricted.nrows()))]
                        break
                if selected is not None:
                    break
        if selected is None:
            raise ArithmeticError("nondegenerate Gram--Schmidt step has no nonisotropic vector")
        orthogonal_basis.append(vector(SageQQ, selected) * complement)
    diagonalizer = matrix(cast(Parent, SageQQ), orthogonal_basis) if orthogonal_basis else identity_matrix(cast(Parent, SageQQ), 0)
    block = identity_matrix(cast(Parent, SageQQ), rank)
    for i in range(dimension):
        for j in range(dimension):
            block[radical_rank + i, radical_rank + j] = diagonalizer[i, j]
    transform = cast(Matrix_rational_dense, block * transform)
    diagonal = cast(Matrix_rational_dense, transform * gram * transform.transpose())
    if any(diagonal[i, j] for i in range(rank) for j in range(rank) if i != j):
        raise ArithmeticError("source symmetric diagonalization left off-diagonal entries")
    if not any(diagonal[i, i] < 0 for i in range(rank)):
        raise ValueError("the Gram matrix has no negative direction")
    candidates: list[tuple[Rational, ...]] = []
    for i in range(rank):
        if diagonal[i, i] > 0:
            continue
        row = tuple(SageQQ(entry) for entry in transform.row(i))
        scale = max(abs(entry) for entry in row)
        if scale:
            candidates.append(tuple(entry / scale for entry in row))
    if not candidates:
        raise ValueError("the Gram matrix has no nonpositive diagonal direction")
    multiplier = 1
    while True:
        for row in candidates:
            # Fractions.g:138-148, NearestInteger rounds ties toward +infinity.
            integral = tuple(_pinned_nearest_integer(multiplier * entry) for entry in row)
            candidate = vector(SageZZ, integral)
            if (candidate * gram * candidate.column())[0] < 0:
                return integral
        multiplier += 1


def _negative_direction_in_span(
    lattice: Lattices.ParentMethods,
    spanning_rows: Matrix_rational_dense,
    critical: tuple[tuple[int, ...], ...],
    base_normal: FreeModuleElement[Rational],
    excluded_directions: list[FreeModuleElement[Rational]] | None = None,
) -> FreeModuleElement[Rational]:
    r"""The negative direction accepted by Lorentzian.g's GetOneOutsideRay.

    The source perturbs the basis until GetUpperBound does not report its
    inadmissible isotropic endpoint. A single arbitrary negative vector is
    insufficient for that algorithm.
    """
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    dimension = spanning_rows.nrows()
    if dimension == 1 and excluded_directions:
        raise ArithmeticError("a rejected outside ray spans the entire one-dimensional negative space")

    def projective_key(row: FreeModuleElement[Rational]) -> tuple[int, ...]:
        primitive = tuple(int(entry) for entry in _source_remove_fraction(row))
        first = next((entry for entry in primitive if entry), 0)
        return tuple(-entry for entry in primitive) if first < 0 else primitive

    excluded_keys = {projective_key(direction) for direction in (excluded_directions or [])}
    # Keep the exact perturbation over QQ, matching the rational span basis.
    # Mixed ZZ/QQ matrix multiplication dispatches through the generic Sage
    # coercion graph on every random trial; embedding these unimodular entries
    # once preserves the source transformations without that repeated dispatch.
    perturbation = identity_matrix(SageQQ, dimension)
    while True:
        rows = perturbation * spanning_rows
        restricted = rows * gram * rows.transpose()
        negative_row = vector(SageQQ, _source_some_negative_vector(restricted))
        direction = negative_row * rows
        covector = gram * direction.column()
        scalars = {(vector(SageQQ, row) * covector)[0] for row in critical}
        if len(scalars) != 1:
            raise ArithmeticError("the GetOneOutsideRay direction must pair equally with every critical vector")
        (scalar,) = scalars
        normal = vector(SageQQ, [-scalar, *covector.column(0)])
        if projective_key(direction) not in excluded_keys and _upper_bound(gram, base_normal, normal) is not None:
            return direction
        # IndefiniteFormsFundamental.g:11-30, GetRandomMatrixPerturbation.
        # Match GetRandomMatrixPerturbation verbatim: a diagonal 1 is replaced
        # by +/-1, while an off-diagonal entry adds an elementary shear. The
        # diagonal -1 overwrite is singular in the pinned source; preserve it.
        if choice((1, 2)) == 1:
            permutation = list(range(dimension))
            shuffle(permutation)
            signs = tuple(choice((-1, 1)) for _ in range(dimension))
            modification = _source_permutation_sign_matrix(permutation, signs)
        else:
            index = randrange(dimension)
            image = randrange(dimension)
            sign = choice((-1, 1))
            modification = _source_elementary_perturbation(dimension, index, image, sign)
        perturbation *= modification


def _source_permutation_sign_matrix(
    permutation: list[int], signs: tuple[int, ...]
) -> Matrix_rational_dense:
    r"""The signed permutation branch of GetRandomMatrixPerturbation."""
    dimension = len(permutation)
    if sorted(permutation) != list(range(dimension)) or len(signs) != dimension or any(s not in (-1, 1) for s in signs):
        raise ValueError("a signed permutation requires a permutation and one sign per row")
    result = matrix(SageQQ, dimension, dimension)
    for index, image in enumerate(permutation):
        result[index, image] = signs[index]
    return result


def _source_elementary_perturbation(
    dimension: int, row: int, column: int, sign: int
) -> Matrix_rational_dense:
    r"""The identity-with-one-entry branch of GetRandomMatrixPerturbation."""
    if dimension < 1 or not 0 <= row < dimension or not 0 <= column < dimension or sign not in (-1, 1):
        raise ValueError("an elementary perturbation needs valid indices and sign ±1")
    result = identity_matrix(SageQQ, dimension)
    result[row, column] = sign
    return result


def _kernel_flipping(
    lattice: Lattices.ParentMethods,
    critical: tuple[tuple[int, ...], ...],
    base_normal: FreeModuleElement[Rational],
    direction_normal: FreeModuleElement[Rational],
    allowed_squares: RealSet,
) -> tuple[tuple[tuple[int, ...], ...], FreeModuleElement[Rational], FreeModuleElement[Rational], Rational]:
    r"""LORENTZ_Kernel_Flipping (Lorentzian.g:842-968) on coordinate rows of the critical set."""
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    # Lorentzian.g:853-860. Both affine normals annihilate the critical
    # configuration and are independent. Those are hypotheses of the
    # one-parameter wall crossing, not optional diagnostics.
    if matrix(SageQQ, [base_normal, direction_normal]).rank() != 2:
        raise ArithmeticError("Lorentzian flip normals must be independent")
    for row in critical:
        augmented = vector(SageQQ, [1, *row])
        if base_normal.dot_product(augmented) != 0 or direction_normal.dot_product(augmented) != 0:
            raise ArithmeticError("both Lorentzian flip normals must vanish on the critical vectors")
        point = vector(SageQQ, row)
        if allowed_squares == ZERO_SQUARES and point * gram * point != 0:
            raise ArithmeticError("an isotropic flip requires isotropic critical vectors")
    upper = _upper_bound(gram, base_normal, direction_normal)
    if upper is None:
        raise ArithmeticError("the flip direction ends at a forbidden isotropic boundary")
    lower = SageQQ.zero()
    inverse = gram.inverse()
    critical_first = vector(SageQQ, critical[0])
    critical_rows = set(critical)
    total: tuple[tuple[int, ...], ...] = ()
    while True:
        middle = _source_mid_value(lower, upper)
        normal = base_normal + middle * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = _source_remove_fraction(test_direction.column(0))
        square = (test_direction_row * gram * test_direction_row.column())[0]
        max_scal = (critical_first * gram * test_direction_row.column())[0]
        if square <= 0 or max_scal <= 0:
            upper = middle
            continue
        # Lorentzian.g:918 evaluates the shortest shell at the trial normal.
        total = _first_positive_shell(lattice, tuple(test_direction_row), max_scal, allowed_squares)
        total_rows = set(total)
        # GAP IsSubset(A, B) means B is contained in A (GAP Reference 30.5-1).
        # Lorentzian.g:922-924 rejects a proper subset of the critical shell.
        if total_rows.issubset(critical_rows) and len(critical_rows) > len(total_rows):
            raise ArithmeticError("the shortest shell cannot be a proper subset of the critical configuration")
        if total_rows == critical_rows:
            lower = middle
            continue
        # Lorentzian.g:924-930, EXIT 1: the shell contains CritSet.
        if critical_rows.issubset(total_rows):
            return total, normal, test_direction_row, max_scal
        # Lorentzian.g:922-930 enters the second scheme when the shell does
        # not contain CritSet. Do not add a transversality prefilter here:
        # the pinned algorithm selects ListTotal[1] verbatim in scheme two.
        break
    while True:
        if not total:
            raise ArithmeticError("the second flip scheme requires a nonempty shortest shell")
        # Lorentzian.g:952 selects ListTotal[1] verbatim. Preserve the
        # enumerator's first result; do not substitute a norm-based ordering.
        expanded = vector(SageQQ, [1, *total[0]])
        denominator = direction_normal.dot_product(expanded)
        if denominator == 0:
            raise ArithmeticError(
                "the first source-ordered vector does not determine a flip parameter: "
                f"row={total[0]!r}, critical={tuple(int(entry) for entry in total[0]) in critical_rows}"
            )
        shift = -base_normal.dot_product(expanded) / denominator
        normal = base_normal + shift * direction_normal
        test_direction = inverse * vector(SageQQ, normal[1:]).column()
        test_direction_row = _source_remove_fraction(test_direction.column(0))
        max_scal = (critical_first * gram * test_direction_row.column())[0]
        total = _first_positive_shell(lattice, tuple(test_direction_row), max_scal, allowed_squares)
        # Lorentzian.g:959-966, EXIT 2: ListTotal contains CritSet
        # together with the single vector selected above.
        if (critical_rows | {tuple(int(entry) for entry in expanded[1:])}).issubset(set(total)):
            return total, normal, test_direction_row, max_scal
