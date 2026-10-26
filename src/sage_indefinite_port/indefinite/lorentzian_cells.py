"""Exact Lorentzian perfect-cell local backend."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from functools import cache, cached_property
from itertools import combinations, product
from math import gcd, isfinite, isqrt
from sys import float_info
from random import choice, randrange, shuffle
from typing import Literal, cast, overload

import numpy as np
from dzack_research.preamble.all import FramedFreeModules, Lattices
from dzack_research.preamble.categories.definite_lattices import LatticeReduction, _ExactCVPEngine
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.rings.ring_foundation import _owned_engine_element
from dzack_research.preamble.categories.modules.pure.modules import _engine_matrix
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones, _cone_from_engine_polyhedron
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


@cache
def perfect_complex(lattice: Lattices.ParentMethods, allowed_squares: RealSet) -> LorentzianPerfectComplex:
    r"""The perfect-domain complex of a lattice, memoized as IndefiniteForms.g:1173-1299 memoizes its groups.

    The traversal of a complex is deterministic data of the lattice; every consumer of
    the same lattice shares its cells, stabilizers and side pairings.
    """
    return LorentzianPerfectComplex(lattice, allowed_squares)


def lorentzian_orthogonal_group_generators(lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
    r"""Generators of O(L), with the pinned split rank-two branch."""
    if int(lattice.rank()) == 1:
        label = next(iter(lattice.module_generating_set()))
        negative = lattice.linear_combination({label: -lattice.base_ring().one()})
        return (lattice.Isom(lattice)({label: negative}),)
    if int(lattice.rank()) == 2 and lattice.is_isotropic():
        return _binary_split_automorphisms(lattice)
    return perfect_complex(lattice, NONNEGATIVE_SQUARES).full_orthogonal_group_generators()


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


_PAIRING_BLOCK_ROWS = 1024


class _ConfigurationPairings:
    r"""The Gram pairings b(v_i, v_j) of a vector configuration, read from its n x rank rows.

    Rows of the n x n pairing matrix are produced in blocks and single
    entries and columns on request, so the matrix is never stored.
    """

    def __init__(self, rows: tuple[tuple[int, ...], ...], gram: Matrix_integer_dense) -> None:
        # Bound with Python integers *before* conversion: np.abs(-2**63)
        # overflows, and converting an out-of-range lattice coefficient to
        # int64 must not silently change the integral Gram pairing.
        maximum_row = max(abs(int(entry)) for row in rows for entry in row)
        maximum_gram = max(abs(int(entry)) for row in gram.rows() for entry in row)
        bound = gram.nrows() ** 2 * maximum_row ** 2 * maximum_gram
        if bound >= 2**40:
            raise ArithmeticError("perfect-cell pairings exceed exact 64-bit integer arithmetic")
        self.rows = np.array(rows, dtype=np.int64)
        form = np.array([[int(entry) for entry in row] for row in gram.rows()], dtype=np.int64)
        self.covectors = self.rows @ form

    def __len__(self) -> int:
        return self.rows.shape[0]

    def __getitem__(self, index: tuple[int, int]) -> int:
        left, right = index
        return int(self.covectors[left] @ self.rows[right])

    def diagonal(self):
        return np.einsum("ij,ij->i", self.covectors, self.rows)

    def column(self, index: int):
        # The form is symmetric, so column j holds b(v_l, v_j) for every l.
        return self.rows @ self.covectors[index]

    def blocks(self):
        for start in range(0, len(self), _PAIRING_BLOCK_ROWS):
            yield start, self.covectors[start : start + _PAIRING_BLOCK_ROWS] @ self.rows.T


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
    # Pairings are exact int64, but the refinement encoding multiplies
    # them by the number of colour classes.  That product need not fit
    # int64 even when every pairing does.  Promote before multiplication.
    colour_array = np.array(colours, dtype=object)
    signatures = []
    for start, block in pairings.blocks():
        codes = block.astype(object) * modulus + colour_array[None, :]
        for offset, row in enumerate(codes):
            index = start + offset
            signatures.append((colours[index], _row_multiset(row, int(row[index]))))
    return tuple(signatures)


def _anchor_profile(pairings: _ConfigurationPairings, colours: tuple[int, ...], anchors: tuple[int, ...]) -> Counter:
    r"""The multiset of (colour, pairings with the anchors) over the configuration."""
    columns = np.stack([pairings.column(anchor) for anchor in anchors], axis=1)
    return Counter((colours[index], tuple(int(entry) for entry in columns[index])) for index in range(len(pairings)))


class LorentzianPerfectLocalBackend:
    @cached_method
    def cell_cone(self, cell: LorentzianPerfectCell):
        return RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)

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
        labels = tuple(lattice.module_generating_set())
        coordinate_rows = tuple(tuple(int(entry) for entry in _coordinate_row(lattice, point, SageZZ)) for point in points)
        positions = {row: position for position, row in enumerate(coordinate_rows)}
        if len(positions) != len(points):
            raise ArithmeticError("a stable point configuration contains repeated points")
        point_matrix = matrix(SageZZ, len(points), len(labels), lambda row, col: coordinate_rows[row][col])
        permutations = []
        for generator in self.cell_stabilizer(cell):
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

    @cached_method
    def facet_orbit_permutations(self, cell: LorentzianPerfectCell) -> tuple[tuple[tuple[FacetIncidence, tuple[int, ...]], ...], ...]:
        r"""Orbits of Aut((L, S)) on facets, paired with their exact ray-permutation transversals.

        Lorentzian.g:1318-1322 (PermutationStabilizer := GRP_int) and
        OnSetsGroupFormalism: the stabilizer acts on facets through its
        permutation action on the cell's rays. The Schreier transversals are
        retained as ray permutations; ``facet_orbits`` materializes the owned
        lattice isometries only for callers requiring categorical morphisms.
        """
        vectors = cell.vector_configuration
        cone = self.cell_cone(cell)
        # The stored homogeneous H-presentation is already minimized. Asking
        # for facet_covectors() reconstructs each face cone only to recover
        # these same supporting normals, which dominates high-facet cells.
        covectors = tuple(cone.halfspace_covectors())
        if not covectors:
            return ()
        labels = tuple(cell.lattice.module_generating_set())
        ray_rows = matrix(SageQQ, cell.rows)
        integral_ray_rows = matrix(SageZZ, cell.rows)
        covector_rows = matrix(SageQQ, [[_rational(covector.to_vector()(label)) for label in labels] for covector in covectors])
        pairings = covector_rows * ray_rows.transpose()
        facets = [frozenset(position for position in range(len(vectors)) if pairings[facet, position] == 0) for facet in range(len(covectors))]
        generators = self.cell_stabilizer(cell)
        # Use precisely these generators to form the ray action.  Calling
        # point_permutations here repeats the stabilizer lookup and the
        # categorical basis-action conversion for every generator.
        positions = {tuple(int(entry) for entry in row): index for index, row in enumerate(integral_ray_rows.rows())}
        if len(positions) != len(vectors):
            raise ArithmeticError("a perfect-cell configuration contains repeated rays")
        permutations = []
        for generator in generators:
            action = matrix(SageZZ, [_coordinate_row(cell.lattice, generator(cell.lattice.module_generator(label)), SageZZ) for label in labels])
            try:
                permutations.append(tuple(positions[tuple(row)] for row in (integral_ray_rows * action).rows()))
            except KeyError as exc:
                raise ArithmeticError("a cell-stabilizer generator does not preserve its rays") from exc
        facet_set = set(facets)
        # Every stabilizer generator preserves the cone, hence must permute
        # its facets.  Verify this before forming Schreier transporters;
        # otherwise the unvisited-only traversal could silently lose images.
        for permutation in permutations:
            images = {frozenset(permutation[position] for position in facet) for facet in facet_set}
            if images != facet_set:
                raise ArithmeticError("a cell-stabilizer generator fails to permute the facets")
        unvisited = set(facets)
        orbits = []
        identity_permutation = tuple(range(len(vectors)))
        for start in facets:
            if start not in unvisited:
                continue
            unvisited.discard(start)
            orbit = [(start, identity_permutation)]
            for facet, element_permutation in orbit:
                for permutation in permutations:
                    image = frozenset(permutation[position] for position in facet)
                    if image in unvisited:
                        unvisited.discard(image)
                        composed = tuple(permutation[element_permutation[position]] for position in range(len(vectors)))
                        orbit.append((image, composed))
            orbits.append(tuple(orbit))
        return tuple(orbits)

    @cached_method
    def facet_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[tuple[FacetIncidence, LatticeIsometryMethods], ...], ...]:
        """Facet orbits with owned isometries, materialized for this public view."""
        labels = tuple(cell.lattice.module_generating_set())
        ray_rows = matrix(SageQQ, cell.rows)
        basis_positions = tuple(int(index) for index in ray_rows.transpose().pivots())
        if len(basis_positions) != len(labels):
            raise ArithmeticError("a perfect-cell configuration does not span its ambient lattice")
        inverse_basis = ray_rows.matrix_from_rows(basis_positions).inverse()
        result = []
        for orbit in self.facet_orbit_permutations(cell):
            elements = []
            for facet, permutation in orbit:
                image_basis = matrix(SageQQ, [ray_rows[permutation[position]] for position in basis_positions])
                action = inverse_basis * image_basis
                if any(entry.denominator() != 1 for entry in action.list()):
                    raise ArithmeticError("a Schreier facet permutation does not induce an integral lattice action")
                isometry = cell.lattice.O()(tuple(cell.lattice([int(entry) for entry in row]) for row in action.rows()))
                elements.append((facet, isometry))
            result.append(tuple(elements))
        return tuple(result)

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

    @cached_method
    def _integral_configuration_invariants(self, cell: LorentzianPerfectCell):
        r"""Exact pairing data of a configuration, computed once per cell.

        The quotient traversal compares each representative with many
        neighboring cells. The pairings and each vertex's pairing spectrum
        depend only on the configuration; the pairings are read from the
        n x rank rows on request, never stored as an n x n matrix. Cache these
        exact witnesses without replacing the integral-isometry decision by a
        hash.
        """
        rows = matrix(SageZZ, cell.rows)
        pairings = _ConfigurationPairings(cell.rows, _gram_matrix(cell.lattice))
        norms = pairings.diagonal()
        spectra = []
        for start, block in pairings.blocks():
            for offset, row in enumerate(block):
                index = start + offset
                spectra.append((int(norms[index]), gcd(*cell.rows[index]), _row_multiset(row, int(row[index]))))
        return rows, pairings, tuple(spectra), Counter(cell.rows), rows.elementary_divisors()

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

        A spanning subfamily of rank r determines an isometry uniquely from
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
        source_rows, source_pairings, source_spectra, source_multiset, source_smith = self._integral_configuration_invariants(source)
        target_rows, target_pairings, target_spectra, target_multiset, target_smith = self._integral_configuration_invariants(target)
        # This exact certificate subsumes any graph-isomorphism search.
        if not all_images and source.lattice is target.lattice and source_multiset == target_multiset:
            return source.lattice.Aut().identity()
        if source_smith != target_smith:
            return () if all_images else None
        if Counter(source_spectra) != Counter(target_spectra):
            return () if all_images else None
        colour_numbers = {signature: index for index, signature in enumerate(sorted(set(source_spectra)))}
        source_colours = tuple(colour_numbers[signature] for signature in source_spectra)
        target_colours = tuple(colour_numbers[signature] for signature in target_spectra)

        # Refine the two coloured, integral-weighted configuration graphs
        # (the vertex-coloured Gram graph of LinPolytopeIntegral_Isomorphism_GramMat)
        # simultaneously by their complete neighbour-colour distributions.
        # Each refinement signature consists entirely of integral-isometry
        # invariants: the preceding colour and the multiset of (pairing,
        # neighbour-colour) pairs. The partition only splits, hence reaches a
        # fixed point after finitely many steps. Every integral transporter
        # preserves all refinement stages; a discrepancy in the colour
        # multiplicities is an exact obstruction.
        while True:
            modulus = len(set(source_colours) | set(target_colours))
            source_refined = _refined_signatures(source_pairings, source_colours, modulus)
            target_refined = _refined_signatures(target_pairings, target_colours, modulus)
            refined_numbers = {signature: index for index, signature in enumerate(sorted(set(source_refined) | set(target_refined)))}
            refined_source_colours = tuple(refined_numbers[signature] for signature in source_refined)
            refined_target_colours = tuple(refined_numbers[signature] for signature in target_refined)
            if Counter(refined_source_colours) != Counter(refined_target_colours):
                return () if all_images else None
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

        # Any transporter carries the ordered source basis to an ordered
        # target basis with the same lattice index. This determinant
        # obstruction applies before constructing the candidate rational
        # action, and is stronger than testing the Smith invariants of
        # proper subsets alone.
        source_basis_index = abs(source_basis.det())

        source_labels = tuple(source.lattice.module_generating_set())
        target_labels = tuple(target.lattice.module_generating_set())
        source_submodule_invariants = {}
        target_submodule_invariants = {}
        # Smith data depend on the selected row subset, not its order.
        # Cache proper partial spans; full-rank leaves instead receive
        # the complete integral transporter certificate.
        # The source profile depends only on the chosen source anchors,
        # while a target profile depends only on their ordered images.
        # Backtracking revisits both across many different partial maps.
        source_anchor_profiles = {}
        target_anchor_profiles = {}
        stabilizer_elements: list[LatticeIsometryMethods] = []
        if generators_only:
            positions = {tuple(int(entry) for entry in row): index for index, row in enumerate(source_rows.rows())}
            if len(positions) != n:
                raise ArithmeticError("a perfect-cell stabilizer requires distinct configuration vectors")
            symmetric = SymmetricGroup(n)
            subgroup = symmetric.subgroup([])
            group_permutations = []

        def backtrack(images, certified_target_basis=None):
            nonlocal subgroup
            if images:
                # Every extension to a configuration isometry must preserve
                # the complete coloured multiset of pairings to the already
                # assigned vectors, not merely the pairings among the basis
                # images. In particular, a wrong candidate for one basis
                # image is often excluded before trying another basis image.
                anchors = tuple(sorted(images.items()))
                source_anchors = tuple(k for k, _ in anchors)
                target_anchors = tuple(v for _, v in anchors)
                if source_anchors not in source_anchor_profiles:
                    source_anchor_profiles[source_anchors] = _anchor_profile(source_pairings, source_colours, source_anchors)
                if target_anchors not in target_anchor_profiles:
                    target_anchor_profiles[target_anchors] = _anchor_profile(target_pairings, target_colours, target_anchors)
                source_profiles = source_anchor_profiles[source_anchors]
                target_profiles = target_anchor_profiles[target_anchors]
                if source_profiles != target_profiles:
                    return None
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
                if Counter(tuple(row) for row in (source_rows * action).rows()) != target_multiset:
                    return None
                witness = source.lattice.Isom(target.lattice)(
                    {
                        label: target.lattice.linear_combination(
                            {target_label: target.lattice.base_ring()(int(entry)) for target_label, entry in zip(target_labels, row, strict=True) if entry}
                        )
                        for label, row in zip(source_labels, action.rows(), strict=True)
                    }
                )
                if all_images:
                    if generators_only:
                        try:
                            permutation = symmetric([positions[tuple(int(entry) for entry in row)] + 1 for row in (source_rows * action).rows()])
                        except KeyError as exc:
                            raise ArithmeticError("a certified perfect-cell stabilizer fails to permute its configuration") from exc
                        if permutation not in subgroup:
                            stabilizer_elements.append(witness)
                            group_permutations.append(permutation)
                            subgroup = symmetric.subgroup(group_permutations)
                    else:
                        stabilizer_elements.append(witness)
                    return None
                return witness

            used = set(images.values())
            choice = None
            candidates_for_choice = None
            for i in basis:
                if i in images:
                    continue
                admissible_images = set(buckets[source_colours[i]]) - used
                for k, v in images.items():
                    key = (k, v, i)
                    if key not in compatible_images:
                        # The target neighbour partition depends on the
                        # image v alone, not on which source anchor k maps
                        # to v. Reuse it across distinct partial assignments.
                        if v not in target_neighbour_signatures:
                            neighbours = {}
                            for l, pairing in enumerate(target_pairings.column(v)):
                                neighbours.setdefault((target_colours[l], int(pairing)), set()).add(l)
                            target_neighbour_signatures[v] = {signature: frozenset(indices) for signature, indices in neighbours.items()}
                        signature = (source_colours[i], source_pairings[i, k])
                        compatible_images[key] = target_neighbour_signatures[v].get(signature, frozenset())
                    admissible_images.intersection_update(compatible_images[key])
                    if not admissible_images:
                        break
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
                        witness = backtrack(images, target_basis_integral)
                        if witness is not None:
                            return witness
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
                if target_selected not in target_submodule_invariants:
                    target_submodule_invariants[target_selected] = tuple(target_partial.elementary_divisors())
                if source_submodule_invariants[selected] == target_submodule_invariants[target_selected]:
                    witness = backtrack(images)
                    if witness is not None:
                        return witness
                del images[choice]
            return None

        first = backtrack({})
        return tuple(stabilizer_elements) if all_images else first


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
        self._traversal_representative_buckets: dict[PerfectFormKey, list[LorentzianPerfectCell]] = {}
        self._traversal_adjacencies: list[LorentzianCellAdjacency] = []
        self._traversal_processed: set[int] = set()
        self._traversal_facet_orbits = {}

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
        source_orbits = backend.facet_orbit_permutations(source)
        source_id = id(source)
        if sum(len(orbit) for orbit in source_orbits) != len({facet for orbit in source_orbits for facet, _ in orbit}):
            raise ArithmeticError("the cell facet-orbit decomposition contains repeated incidences")
        self._traversal_facet_orbits[source_id] = {orbit[0][0]: orbit for orbit in source_orbits}
        for orbit in source_orbits:
            facet = orbit[0][0]
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
            # LatticeDelaunays.g:145-176 records one adjacency for each
            # stabilizer orbit of facets, not an adjacency for every facet.
            # The other facets are carried to the selected one by the cell
            # stabilizer, whose generators are already retained separately.
            # Recording every orbit member would replicate the same quotient
            # edge and its side pairing, with an expensive polyhedral witness
            # check for each copy.
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
        target_key = _perfect_form_hash_key(target_cell)
        backend = self.local_backend()
        for source_cell in self._traversal_representative_buckets.get(target_key, ()):
            witness = backend.cell_transporter(source_cell, target_cell)
            if witness is not None:
                return witness
        while len(self._traversal_processed) < len(representatives):
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


type PerfectFormKey = tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...], tuple[int, ...], tuple[int, ...]]


def _perfect_form_hash_key(cell: LorentzianPerfectCell) -> PerfectFormKey:
    r"""ComputeInvariantPerfectForm (lorentzian_perfect.h:1212-1256), with two integral invariants.

    The source streams the configuration's Gram pairings and keeps the
    multiplicity of each diagonal value and of each off-diagonal value
    b(v_i, v_j), i < j.  The pairings are streamed here in blocks of rows, so
    no n x n matrix is formed.  The configuration's Smith invariants and the
    contents of its vectors are further necessary conditions for an integral
    isometry, read off the n x rank row matrix.  The exact isometry test
    remains mandatory.
    """
    exact_gram = _gram_matrix(cell.lattice)
    maximum_row = max(abs(int(entry)) for row in cell.rows for entry in row)
    maximum_gram = max(abs(int(entry)) for row in exact_gram.rows() for entry in row)
    bound = exact_gram.nrows() ** 2 * maximum_row ** 2 * maximum_gram
    if bound >= 2**62:
        raise ArithmeticError("perfect-cell pairings exceed exact 64-bit integer arithmetic")
    rows = np.array(cell.rows, dtype=np.int64)
    gram = np.array([[int(entry) for entry in row] for row in exact_gram.rows()], dtype=np.int64)
    covectors = rows @ gram
    diagonal = Counter(int(value) for value in np.einsum("ij,ij->i", covectors, rows))
    all_pairs: Counter[int] = Counter()
    for start in range(0, rows.shape[0], _PAIRING_BLOCK_ROWS):
        block = rows[start : start + _PAIRING_BLOCK_ROWS] @ covectors.T
        values, counts = np.unique(block, return_counts=True)
        all_pairs.update(dict(zip((int(value) for value in values), (int(count) for count in counts), strict=True)))
    # The full square counts each unordered pair twice and the diagonal once.
    off_diagonal = {value: (all_pairs[value] - diagonal.get(value, 0)) // 2 for value in all_pairs}
    smith_invariants = tuple(int(entry) for entry in matrix(SageZZ, cell.rows).elementary_divisors())
    vector_contents = tuple(sorted(gcd(*row) for row in cell.rows))
    return (
        tuple(sorted(diagonal.items())),
        tuple(sorted((value, count) for value, count in off_diagonal.items() if count)),
        smith_invariants,
        vector_contents,
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
    cell_cones = {id(cell): backend.cell_cone(cell) for cell in cells}
    # A supporting facet is the zero section of its defining covector.
    # Reuse the source cell's H-presentation for every member of its
    # stabilizer orbit instead of recomputing a V-to-H representation.
    facet_covectors_by_incidence = {}
    halfspace_rows_by_cell = {}
    for cell in cells:
        labels = tuple(cell.lattice.module_generating_set())
        rows = matrix(SageQQ, cell.rows)
        by_incidence = {}
        covectors = tuple(cell_cones[id(cell)].halfspace_covectors())
        halfspace_rows = []
        # The cone's stored H-presentation already is the minimized facet
        # family. `facet_covectors()` constructs every facet cone merely to
        # recover these same normals, which is prohibitively expensive on
        # high-facet cells.
        for covector in covectors:
            coordinate = covector.to_vector()
            normal = vector(SageQQ, (_rational(coordinate(label)) for label in labels))
            halfspace_rows.append(normal)
            incidence = frozenset(position for position, row in enumerate(rows.rows()) if row.dot_product(normal) == 0)
            if not incidence:
                raise ArithmeticError("a minimized H-presentation covector has no incident extreme rays")
            if incidence in by_incidence:
                raise ArithmeticError("distinct perfect-cell facets have the same extreme-ray incidence")
            by_incidence[incidence] = covector
        facet_covectors_by_incidence[id(cell)] = by_incidence
        halfspace_rows_by_cell[id(cell)] = matrix(SageQQ, halfspace_rows)
    ambient_labels = tuple(lattice.module_generating_set())
    ambient_dual = lattice.dual_module()
    dual_labels = tuple(ambient_dual.module_generating_set())
    ambient_gram = _gram_matrix(lattice)
    neighbor_cache = {}
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
        for incident, symmetry_permutation in orbit:
            # In row coordinates g(v)=vM. A covector column a cutting out
            # a>=0 on the target becomes M^{-1}a on its image. Since g is
            # integral unimodular, the transformed H-data remain integral.
            # This exact transformation avoids a new V-to-H conversion for
            # every member of a facet-stabilizer orbit.
            image_basis = matrix(
                SageQQ,
                [adjacency.source.rows[symmetry_permutation[position]] for position in basis_positions],
            )
            symmetry_action = inverse_source_basis * image_basis
            if any(entry.denominator() != 1 for entry in symmetry_action.list()):
                raise ArithmeticError("a Schreier stabilizer permutation does not induce an integral lattice action")
            symmetry_action = symmetry_action.change_ring(SageZZ)
            # In row coordinates, symmetry * transporter acts by the product
            # M_transporter M_symmetry.  The ray permutation determines the
            # exact integral action; check a nontrivial Schreier element in
            # each facet orbit against the owned categorical morphism.
            action = transporter_action * symmetry_action
            if symmetry_action * ambient_gram * symmetry_action.transpose() != ambient_gram or abs(symmetry_action.det()) != 1:
                raise ArithmeticError("a Schreier stabilizer permutation does not induce an integral lattice isometry")
            if action * ambient_gram * action.transpose() != ambient_gram or abs(action.det()) != 1:
                raise ArithmeticError("a Schreier side-pairing action is not an integral lattice isometry")
            # The Schreier permutation determines the exact action matrix.
            # Defer the owned categorical isometry until a caller requests the
            # side pairing, rather than constructing one for every orbit facet.
            target_to_neighbour_factory = (
                lambda symmetry_action=symmetry_action, transporter=adjacency.transporter, lattice=lattice:
                lattice.O()(tuple(lattice([int(entry) for entry in row]) for row in symmetry_action.rows())) * transporter
            )
            target_cone = cell_cones[id(adjacency.target)]
            target_key = id(adjacency.target)
            target_halfspace_rows = halfspace_rows_by_cell[target_key]
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
            try:
                facet_normal = facet_covectors_by_incidence[id(adjacency.source)][incident]
            except KeyError as exc:
                raise ArithmeticError("a Schreier facet has no defining covector in its source cone") from exc
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
            paired_facets[id(adjacency.source)].add(incident)
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
    for cell in cells:
        by_incidence = facet_covectors_by_incidence[id(cell)]
        rows = matrix(SageQQ, cell.rows)
        for incidence in by_incidence:
            if incidence in paired_facets[id(cell)]:
                continue
            facet_rows = rows.matrix_from_rows(sorted(incidence))
            pairings = facet_rows * ambient_gram * facet_rows.transpose()
            if not any(pairings.list()):
                continue
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
        # lie on its null extreme rays.  The configuration need not
        # consist of primitive generators; enumerate the primitive rays
        # of the actual cone rather than filtering those generators.
        zero = cell.lattice.base_ring().zero()
        marks = tuple(ray for ray in self.complex().local_backend().cell_cone(cell).primitive_rays() if ray.q() == zero)
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
        cone = backend.cell_cone(cell)
        rays = tuple(cone.primitive_rays())
        self.rows = tuple(_coordinate_row(lattice, ray, SageZZ) for ray in rays)
        self.covectors = tuple(vector(SageQQ, [_rational(covector.to_vector()(label)) for label in self.labels]) for covector in cone.facet_covectors())
        interior = sum(self.rows[1:], self.rows[0])
        if any(covector.dot_product(interior) <= 0 for covector in self.covectors):
            raise ArithmeticError("the facet covectors of a cell must be positive on its interior")
        self.incidences = frozenset(frozenset(index for index, row in enumerate(self.rows) if covector.dot_product(row) == 0) for covector in self.covectors)
        self.full_face = frozenset(range(len(self.rows)))
        self.group = libgap.Group([libgap.PermList([image + 1 for image in permutation]) for permutation in backend.point_permutations(cell, rays)], libgap.eval("()"))
        self._ranks: dict[frozenset[int], int] = {}

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
        chain: list[frozenset[int]] = []
        coefficients: list[Rational] = []
        while point != 0:
            tight = [covector for covector in self.covectors if covector.dot_product(point) == 0]
            face = frozenset(index for index, ray in enumerate(self.rows) if all(covector.dot_product(ray) == 0 for covector in tight))
            barycenter = vector(SageQQ, self.barycenter(face))
            coefficient = min(covector.dot_product(point) / covector.dot_product(barycenter) for covector in self.covectors if covector.dot_product(barycenter) > 0)
            chain.append(face)
            coefficients.append(coefficient)
            point -= coefficient * barycenter
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


def _negative_ldl_factorization(gram: Matrix_rational_dense) -> tuple[list[list[Rational]], list[Rational]]:
    """Exact ``G = L D L^T`` factorization of a positive-definite Gram matrix."""
    if gram.nrows() == 0 or gram.nrows() != gram.ncols() or gram != gram.transpose():
        raise ValueError("the fixed-distance LDL factorization requires a nonempty symmetric Gram matrix")
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
        return (((), bound),) if bound == 0 else ()
    lower, diagonal = _negative_ldl_factorization(gram)
    _validate_fixed_distance_ldl_factorization(gram, (lower, diagonal))
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


def _fixed_distance_dual_lll_solutions(
    gram: Matrix_rational_dense,
    target: tuple[Rational, ...],
    squared_distance: Rational,
    dual_transform: Matrix_integer_dense,
    enumeration_order_transform: Matrix_integer_dense,
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
        reduced_gram, transformed_target, squared_distance
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
    ellipsoid, the exact LDL recursion prunes a coordinate as soon as its
    accumulated square exceeds the prescribed distance.
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
    zero = cast(Rational, SageQQ.zero())
    lower, diagonal = ldl_factorization if ldl_factorization is not None else _negative_ldl_factorization(gram)
    _validate_fixed_distance_ldl_factorization(gram, (lower, diagonal))

    coordinates = [0] * rank
    solutions: list[tuple[int, ...]] = []

    def enumerate_coordinate(index: int, accumulated: Rational) -> None:
        if index < 0:
            if accumulated == bound:
                solutions.append(tuple(coordinates))
            return
        remaining = bound - accumulated
        if remaining < 0:
            return
        center = cast(
            Rational,
            reduced_target[index]
            - sum(
                (lower[later][index] * (coordinates[later] - reduced_target[later]) for later in range(index + 1, rank)),
                zero,
            ),
        )
        squared_radius = cast(Rational, remaining / diagonal[index])
        radius = isqrt(int(squared_radius.numerator() // squared_radius.denominator()))
        floor_center = int(center.numerator() // center.denominator())
        ceil_center = -int((-center).numerator() // (-center).denominator())
        # For r=floor(sqrt(squared_radius)), every integer in the true
        # interval lies in [ceil(center)-r-1, floor(center)+r+1]. Exact
        # comparisons below trim the endpoints. This also covers rational
        # centers and square or nonsquare squared radii uniformly.
        first_coordinate = ceil_center - radius - 1
        last_coordinate = floor_center + radius + 1
        # The conservative endpoints may include one extra integer when the
        # radius is nonintegral. Remove those endpoints exactly. The squared
        # distance is monotone away from the center, so this leaves precisely
        # the integer interval satisfying the remaining-distance inequality.
        while first_coordinate <= last_coordinate and diagonal[index] * (SageQQ(first_coordinate) - center) ** 2 > remaining:
            first_coordinate += 1
        while first_coordinate <= last_coordinate and diagonal[index] * (SageQQ(last_coordinate) - center) ** 2 > remaining:
            last_coordinate -= 1
        for coordinate in range(first_coordinate, last_coordinate + 1):
            term = cast(Rational, diagonal[index] * (SageQQ(coordinate) - center) ** 2)
            if term <= remaining:
                coordinates[index] = coordinate
                enumerate_coordinate(index - 1, cast(Rational, accumulated + term))

    enumerate_coordinate(rank - 1, zero)
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
        list[tuple[ComplexDoubleElement, list[FreeModuleElement[ComplexDoubleElement]], int]],
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
    integral_entries = tuple(SageZZ(entry) for entry in candidate)
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
            raise ArithmeticError("positive-vector perturbation is not an exact integral lattice automorphism") from error
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
        found = vector(SageZZ, found_raw)
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
        # Use the same exact PARI-backed CVP engine as close_vectors, but its
        # coordinate-level primitive avoids allocating an owned lattice vector
        # and indexed-family object for every intermediate shell candidate.
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
            )
        else:
            distances = self.cvp_engine.close_vector_coordinates(target, bound)
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
        # Lorentzian.g:918 evaluates the shortest shell, which may be empty
        # at an intermediate rational wall parameter.  An empty shell is
        # handled by EXIT 1 as a subset of CritSet, not as an exception.
        total = _first_positive_shell(lattice, tuple(test_direction_row), max_scal, allowed_squares)
        total_rows = set(total)
        # Lorentzian.g:922-924 rejects the impossible cardinality combination
        # CritSet subset ListTotal and |CritSet| > |ListTotal|.
        if critical_rows.issubset(total_rows) and len(critical_rows) > len(total_rows):
            raise ArithmeticError("the critical configuration cannot exceed a shell containing it")
        if total_rows == critical_rows:
            lower = middle
            continue
        # Lorentzian.g:924-930, EXIT 1: the shell lies inside CritSet.
        if total_rows.issubset(critical_rows):
            return total, normal, test_direction_row, max_scal
        # Lorentzian.g:922-930 enters the second scheme for every shell not
        # contained in CritSet. Do not add a transversality prefilter here:
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
        # Lorentzian.g:959-966, EXIT 2: ListTotal is contained in CritSet
        # together with the single vector selected above.
        if set(total).issubset(critical_rows | {tuple(int(entry) for entry in expanded[1:])}):
            return total, normal, test_direction_row, max_scal
