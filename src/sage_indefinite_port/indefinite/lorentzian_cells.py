"""Exact Lorentzian perfect-cell local backend."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from fractions import Fraction
from functools import cached_property
from itertools import product
from math import gcd, isqrt
from random import choice, randrange, shuffle
from typing import Literal, cast, overload

from dzack_research.preamble.all import FormModules, FramedFreeModules, Lattices
from dzack_research.preamble.categories.definite_lattices import LatticeReduction
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.pure.modules import _engine_matrix
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones
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
from dzack_research.preamble.categories.sets.set_categories import Sets
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.matrix.constructor import matrix
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
from sage.structure.element import Element, RingElement
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.indefinite_reduction import _indefinite_reduction
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
    generators = _binary_split_isometries(lattice, lattice)
    if not generators:
        raise ArithmeticError("the split binary isometry enumeration omitted the identity")
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


def lorentzian_orthogonal_group_generators(lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
    r"""Generators of O(L), with the pinned split rank-two branch."""
    if int(lattice.rank()) == 1:
        label = next(iter(lattice.module_generating_set()))
        negative = lattice.linear_combination({label: -lattice.base_ring().one()})
        return (lattice.Isom(lattice)({label: negative}),)
    if int(lattice.rank()) == 2 and lattice.is_isotropic():
        return _binary_split_isometries(lattice, lattice)
    return LorentzianPerfectComplex(lattice, NONNEGATIVE_SQUARES).full_orthogonal_group_generators()


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
        # Lorentzian.g:303-312 (LORENTZ_CheckCorrectnessVectorFamily):
        # the augmented vectors (1,v) span an affine hyperplane of rank n.
        # Both GetOnePerfect and DoFlipping certify this before returning.
        augmented = matrix(SageQQ, [[1, *row] for row in rows])
        if augmented.rank() != gram.nrows():
            raise ArithmeticError("a perfect-cell vector family must have augmented rank equal to the lattice rank")


@dataclass(frozen=True)
class LorentzianCellAdjacency:
    source: LorentzianPerfectCell
    facet: FacetIncidence
    target: LorentzianPerfectCell
    transporter: LatticeIsometryMethods | None


class LorentzianPerfectLocalBackend:
    @cached_method
    def cell_cone(self, cell: LorentzianPerfectCell):
        return RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)

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
    def cell_family(self, cell: LorentzianPerfectCell) -> Element:
        r"""The perfect-cell family with its intrinsic pairing-spectrum colouring.

        The colour of a vector is its norm together with the multiset of
        pairings against the other vectors. An isometry of configurations
        preserves these colours, and conversely every pairing-preserving
        bijection preserves them. Thus the coloured and uncoloured families
        have exactly the same integral stabilizers and transporters, while
        the labelled-graph engine can distinguish the colour classes before
        searching for a pairing-graph isomorphism.
        """
        vectors = cell.vector_configuration
        rows = matrix(SageZZ, [_coordinate_row(cell.lattice, point, SageZZ) for point in vectors])
        pairings = rows * _gram_matrix(cell.lattice) * rows.transpose()
        signatures = tuple(
            (
                int(pairings[i, i]),
                gcd(*(int(entry) for entry in rows.row(i))),
                tuple(sorted(int(pairings[i, j]) for j in range(len(vectors)) if i != j)),
            )
            for i in range(len(vectors))
        )
        classes = {signature: index for index, signature in enumerate(sorted(set(signatures)))}
        colours = Sets.Δ[len(classes) - 1]
        points = Sets().product((cell.lattice, colours))
        return Sets().Mor(Sets.Δ[len(vectors) - 1], points)(lambda index: points((vectors[int(index)], colours[classes[signatures[int(index)]]])))

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
    def facet_orbits(self, cell: LorentzianPerfectCell) -> tuple[tuple[tuple[FacetIncidence, LatticeIsometryMethods], ...], ...]:
        r"""Orbits of Aut((L, S)) on the facets spanned by S, each facet paired with an element carrying the orbit's first facet onto it.

        Lorentzian.g:1318-1322 (PermutationStabilizer := GRP_int) and
        OnSetsGroupFormalism: the stabilizer acts on facets through its
        permutation action on the cell's rays. The orbits and transversal are
        a Schreier traversal over the generators' ray permutations.
        """
        vectors = cell.vector_configuration
        cone = self.cell_cone(cell)
        covectors = tuple(cone.facet_covectors())
        if not covectors:
            return ()
        labels = tuple(cell.lattice.module_generating_set())
        ray_rows = matrix(SageQQ, [_coordinate_row(cell.lattice, point, SageQQ) for point in vectors])
        covector_rows = matrix(SageQQ, [[_rational(covector.to_vector()(label)) for label in labels] for covector in covectors])
        pairings = covector_rows * ray_rows.transpose()
        facets = [frozenset(position for position in range(len(vectors)) if pairings[facet, position] == 0) for facet in range(len(covectors))]
        generators = self.cell_stabilizer(cell)
        # Use precisely these generators to form the ray action.  Calling
        # point_permutations here repeats the stabilizer lookup and the
        # categorical basis-action conversion for every generator.
        positions = {tuple(int(entry) for entry in row): index for index, row in enumerate(ray_rows.rows())}
        if len(positions) != len(vectors):
            raise ArithmeticError("a perfect-cell configuration contains repeated rays")
        permutations = []
        for generator in generators:
            action = matrix(SageZZ, [_coordinate_row(cell.lattice, generator(cell.lattice.module_generator(label)), SageZZ) for label in labels])
            try:
                permutations.append(tuple(positions[tuple(int(entry) for entry in row)] for row in (ray_rows * action).rows()))
            except KeyError as exc:
                raise ArithmeticError("a cell-stabilizer generator does not preserve its rays") from exc
        identity = cell.lattice.Aut().identity()
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
        for start in facets:
            if start not in unvisited:
                continue
            unvisited.discard(start)
            orbit = [(start, identity)]
            for facet, element in orbit:
                for generator, permutation in zip(generators, permutations, strict=True):
                    image = frozenset(permutation[position] for position in facet)
                    if image in unvisited:
                        unvisited.discard(image)
                        orbit.append((image, generator * element))
            orbits.append(tuple(orbit))
        return tuple(orbits)

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
        normal_direction = _source_remove_fraction(vector(SageQQ, [0, *direction]))
        expanded = matrix(SageQQ, [[1, *_coordinate_row(normalized, item, SageQQ)] for item in vectors])
        base_nullspace = expanded.right_kernel_matrix()
        if base_nullspace.nrows() != 1:
            raise ArithmeticError("a perfect cell must have a unique affine supporting functional")
        base_normal = base_nullspace.row(0)
        if vector(SageQQ, base_normal[1:]).dot_product(outside) <= 0:
            base_normal = -base_normal
        base_normal = _source_remove_fraction(base_normal)
        critical = tuple(vectors[position] for position in sorted(facet))
        flipped, _normal, _direction, _max_scal = _kernel_flipping(normalized, critical, base_normal, normal_direction, cell.allowed_squares)
        returned = tuple(_transfer(cell.lattice, item) for item in flipped)
        return LorentzianPerfectCell(cell.lattice, returned, cell.allowed_squares)

    @cached_method
    def _integral_configuration_invariants(self, cell: LorentzianPerfectCell):
        r"""Exact Gram and determinantal-divisor data, computed once per cell.

        The quotient traversal compares each representative with many
        neighboring cells. Pairings and exterior-square contents depend
        only on the individual configuration. Cache these exact witnesses
        without replacing the integral-isometry decision by a hash.
        """
        rank = int(cell.lattice.module_rank())
        rows = matrix(SageZZ, [_coordinate_row(cell.lattice, point, SageZZ) for point in cell.vector_configuration])
        pairings = rows * _gram_matrix(cell.lattice) * rows.transpose()
        coordinates = tuple(tuple(int(entry) for entry in row) for row in rows.rows())
        n = rows.nrows()
        divisors = [[0] * n for _ in range(n)]
        for i in range(n):
            u = coordinates[i]
            for j in range(i + 1, n):
                v = coordinates[j]
                value = gcd(*(u[a] * v[b] - u[b] * v[a] for a in range(rank) for b in range(a + 1, rank)))
                divisors[i][j] = divisors[j][i] = value
        pair_divisors = tuple(tuple(row) for row in divisors)
        spectra = tuple(
            (
                int(pairings[i, i]),
                gcd(*coordinates[i]),
                tuple(sorted((int(pairings[i, j]), pair_divisors[i][j]) for j in range(n) if j != i)),
            )
            for i in range(n)
        )
        return (
            rows,
            pairings,
            pair_divisors,
            spectra,
            Counter(coordinates),
            rows.elementary_divisors(),
        )

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
        its r images. Enumerate all norm-, pairing-, and determinantal-divisor-
        compatible target images of a fixed source basis; at each leaf verify
        integrality, unimodularity, the entire Gram form, and the exact
        multiset of configuration vectors. No partial graph isomorphism is
        treated as a certificate, and no potential integral transporter is
        excluded by these necessary conditions.
        """
        if source.allowed_squares != target.allowed_squares:
            raise ValueError("an isometry of perfect cells must preserve their prescribed normalized-norm locus")
        if generators_only and (not all_images or source is not target):
            raise ValueError("stabilizer-generator search requires one identical source and target perfect cell")
        vectors_source = source.vector_configuration
        vectors_target = target.vector_configuration
        n = len(vectors_source)
        if n != len(vectors_target):
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
        (
            source_rows,
            source_pairings,
            source_pair_divisors,
            source_spectra,
            source_multiset,
            source_smith,
        ) = self._integral_configuration_invariants(source)
        (
            target_rows,
            target_pairings,
            target_pair_divisors,
            target_spectra,
            target_multiset,
            target_smith,
        ) = self._integral_configuration_invariants(target)
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
        # simultaneously by their complete neighbour-colour distributions.
        # Each refinement signature consists entirely of integral-isometry
        # invariants: the preceding colour and the multiset of (pairing,
        # wedge-content, neighbour-colour) triples. The partition only
        # splits, hence reaches a fixed point after finitely many steps.
        # Every integral transporter preserves all refinement stages; a
        # discrepancy in the colour multiplicities is an exact obstruction.
        while True:
            source_refined = tuple(
                (source_colours[i], tuple(sorted((int(source_pairings[i, j]), source_pair_divisors[i][j], source_colours[j]) for j in range(n) if i != j))) for i in range(n)
            )
            target_refined = tuple(
                (target_colours[i], tuple(sorted((int(target_pairings[i, j]), target_pair_divisors[i][j], target_colours[j]) for j in range(n) if i != j))) for i in range(n)
            )
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
        # k -> l by the exact Gram pairing and exterior-square content.
        # Only images of the selected basis occur as anchors. Index their
        # possible target neighbours by the exact pairing signature, and
        # materialize a compatibility fibre only when backtracking needs it.
        # Eagerly storing all (i,j,k) fibres costs cubic space in the
        # configuration size even when the search rejects its first branch.
        compatible_images = {}
        target_neighbour_signatures = {}
        # Each target anchor has the same exact pairing/wedge signature
        # partition regardless of the source vertex assigned to it.
        # Build its immutable fibres once; backtracking only intersects them.

        if source_rows.rank() != rank or target_rows.rank() != rank:
            raise ArithmeticError("a perfect-cell vector configuration does not span its lattice")
        # Source basis selection does not depend on the eventual transporter.
        # A rare spectrum minimizes the initial image search; subsequent
        # choices are constrained by the existing source basis.
        basis = []
        while len(basis) < rank:
            eligible = [i for i in range(n) if i not in basis and source_rows.matrix_from_rows([*basis, i]).rank() == len(basis) + 1]
            if not eligible:
                raise ArithmeticError("failed to select a basis of a spanning perfect-cell configuration")
            basis.append(min(eligible, key=lambda i: (len(buckets[source_colours[i]]), i)))
        basis = tuple(basis)
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
                    source_anchor_profiles[source_anchors] = Counter(
                        (source_colours[i], tuple((int(source_pairings[i, k]), source_pair_divisors[i][k]) for k in source_anchors)) for i in range(n)
                    )
                if target_anchors not in target_anchor_profiles:
                    target_anchor_profiles[target_anchors] = Counter(
                        (target_colours[j], tuple((int(target_pairings[j, v]), target_pair_divisors[j][v]) for v in target_anchors)) for j in range(n)
                    )
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
                            for l in range(n):
                                signature = (target_colours[l], target_pairings[l, v], target_pair_divisors[l][v])
                                neighbours.setdefault(signature, set()).add(l)
                            target_neighbour_signatures[v] = {signature: frozenset(indices) for signature, indices in neighbours.items()}
                        signature = (source_colours[i], source_pairings[i, k], source_pair_divisors[i][k])
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
        self._traversal_representative_buckets: dict[
            tuple[
                tuple[int, ...],
                tuple[int, ...],
                tuple[tuple[int, tuple[int, ...]], ...],
                tuple[int, ...],
                tuple[int, ...],
            ],
            list[LorentzianPerfectCell],
        ] = {}
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
            key=lambda index: (len(representatives[index].vector_configuration), index),
        )
        source = representatives[source_position]
        new_representatives: list[LorentzianPerfectCell] = []
        source_orbits = backend.facet_orbits(source)
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


def _perfect_form_hash_key(
    cell: LorentzianPerfectCell,
) -> tuple[
    tuple[int, ...],
    tuple[int, ...],
    tuple[tuple[int, tuple[int, ...]], ...],
    tuple[int, ...],
    tuple[int, ...],
]:
    r"""The upstream pairing multisets, refined by integral configuration invariants.

    ComputeInvariantPerfectForm hashes the multiplicities of diagonal and
    off-diagonal Gram pairings of the configuration. The upstream integral
    isomorphism test then canonicalizes its weighted Gram graph
    (PolytopeEquiStabInt.h:194-220). Each vertex's diagonal value and
    multiset of incident weights are preserved by any graph isomorphism;
    their sorted family therefore refines the initial hash without excluding
    any possible transporter. Further necessary integral conditions are the
    configuration's Smith invariants and the contents of its individual
    vectors: row permutations and integral unimodular coordinate changes
    preserve both. The exact isometry test remains mandatory.
    """
    rows = matrix(SageZZ, [_coordinate_row(cell.lattice, point, SageZZ) for point in cell.vector_configuration])
    pairings = rows * _gram_matrix(cell.lattice) * rows.transpose()
    diagonal = tuple(sorted(int(pairings[index, index]) for index in range(rows.nrows())))
    off_diagonal = tuple(sorted(int(pairings[left, right]) for left in range(rows.nrows()) for right in range(left + 1, rows.nrows())))
    vertex_spectra = tuple(
        sorted((int(pairings[index, index]), tuple(sorted(int(pairings[index, other]) for other in range(rows.nrows()) if other != index))) for index in range(rows.nrows()))
    )
    smith_invariants = tuple(int(entry) for entry in rows.elementary_divisors())
    vector_contents = tuple(sorted(gcd(*(int(entry) for entry in row)) for row in rows.rows()))
    return diagonal, off_diagonal, vertex_spectra, smith_invariants, vector_contents


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
    cells = complex_.quotient_cells()
    cell_cones = {id(cell): backend.cell_cone(cell) for cell in cells}
    # A supporting facet is the zero section of its defining covector.
    # Reuse the source cell's H-presentation for every member of its
    # stabilizer orbit instead of recomputing a V-to-H representation.
    facet_covectors_by_incidence = {}
    for cell in cells:
        labels = tuple(cell.lattice.module_generating_set())
        rows = matrix(SageQQ, [_coordinate_row(cell.lattice, point, SageQQ) for point in cell.vector_configuration])
        by_incidence = {}
        for covector in cell_cones[id(cell)].facet_covectors():
            coordinate = covector.to_vector()
            normal = vector(SageQQ, (_rational(coordinate(label)) for label in labels))
            incidence = frozenset(position for position, row in enumerate(rows.rows()) if row.dot_product(normal) == 0)
            if incidence in by_incidence:
                raise ArithmeticError("distinct perfect-cell facets have the same extreme-ray incidence")
            by_incidence[incidence] = covector
        facet_covectors_by_incidence[id(cell)] = by_incidence
    ambient_labels = tuple(lattice.module_generating_set())
    dual = lattice.dual_module()
    dual_labels = tuple(dual.module_generating_set())
    # Every target cone has one fixed H-presentation.  Schreier transporters
    # change its coordinates, not the source covectors themselves: cross the
    # H-rows from the categorical dual into exact integer columns just once.
    # In large facet orbits this avoids repeating a dual-coordinate
    # conversion for every side pairing.
    target_h_columns = {}
    for cell in cells:
        cone = cell_cones[id(cell)]
        inequalities = tuple(cone.halfspace_covectors())
        equations = tuple(cone.equation_covectors())
        columns = tuple(tuple(int(coordinate(label)) for label in dual_labels) for covector in (*inequalities, *equations) for coordinate in (covector.to_vector(),))
        # Batch the exact integral dual action on all H-rows in one matrix
        # product. Individual matrix-vector products for every facet of every
        # stabilizer orbit multiply the same inverse action unnecessarily.
        target_h_columns[id(cell)] = (
            len(inequalities),
            matrix(SageZZ, len(dual_labels), len(columns), lambda row, col, _columns=columns: _columns[col][row]),
        )
    neighbor_cache = {}
    identity_action = identity_matrix(SageZZ, len(ambient_labels))
    adjacencies = []
    for adjacency in complex_.adjacencies():
        # The GAP quotient (LatticeDelaunays.g:147-171) retains only one
        # edge per facet orbit.  The preamble's complete domain object also
        # requires a witness for each actual facet.  Recover those witnesses
        # using the Schreier transversal already computed for this orbit.
        orbit = complex_._traversal_facet_orbits[id(adjacency.source)][adjacency.facet]
        for incident, symmetry in orbit:
            target_to_neighbour = symmetry * adjacency.transporter
            # In row coordinates g(v)=vM. A covector column a cutting out
            # a>=0 on the target becomes M^{-1}a on its image. Since g is
            # integral unimodular, the transformed H-data remain integral.
            # This exact transformation avoids a new V-to-H conversion for
            # every member of a facet-stabilizer orbit.
            action = matrix(SageZZ, [_coordinate_row(lattice, target_to_neighbour(lattice.module_generator(label)), SageZZ) for label in ambient_labels])
            target_cone = cell_cones[id(adjacency.target)]
            action_key = (id(adjacency.target), tuple(int(entry) for entry in action.list()))
            neighbor = neighbor_cache.get(action_key)
            if neighbor is None:
                if action == identity_action:
                    neighbor = target_cone
                else:
                    inverse_action = action.change_ring(SageQQ).inverse()
                    if any(entry.denominator() != 1 for entry in inverse_action.list()):
                        raise ArithmeticError("a perfect-domain side pairing is not an integral automorphism")
                    inverse_action = inverse_action.change_ring(SageZZ)
                    inequality_count, columns = target_h_columns[id(adjacency.target)]
                    images = inverse_action * columns
                    transported = tuple(
                        dual.linear_combination({label: dual.base_ring()(int(entry)) for label, entry in zip(dual_labels, image, strict=True) if entry})
                        for image in images.columns()
                    )
                    neighbor = cones.from_halfspaces(
                        transported[:inequality_count],
                        equation_covectors=transported[inequality_count:],
                    )
                neighbor_cache[action_key] = neighbor
            # The Schreier orbit records the exact source facet incidence;
            # its cached supporting covector cuts out that face.
            try:
                facet_normal = facet_covectors_by_incidence[id(adjacency.source)][incident]
            except KeyError as exc:
                raise ArithmeticError("a Schreier facet has no defining covector in its source cone") from exc
            common_face = cell_cones[id(adjacency.source)].face_on_covector(facet_normal)
            adjacencies.append(
                PerfectDomainOrbitAdjacency(
                    cell_cones[id(adjacency.source)],
                    cell_cones[id(adjacency.target)],
                    neighbor,
                    common_face,
                    target_to_neighbour,
                )
            )
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
        self._local_mark_cache: dict[int, tuple[Lattices.ElementMethods, ...]] = {}

    def complex(self) -> LorentzianPerfectComplex:
        return self._complex

    def local_marks(self, cell: LorentzianPerfectCell) -> tuple[Lattices.ElementMethods, ...]:
        key = id(cell)
        if key in self._local_mark_cache:
            return self._local_mark_cache[key]
        if self._norm == 0:
            # Null vectors of a pointed cone in one closed Lorentzian sheet
            # lie on its null extreme rays.  The configuration need not
            # consist of primitive generators; enumerate the primitive rays
            # of the actual cone rather than filtering those generators.
            zero = cell.lattice.base_ring().zero()
            marks = tuple(ray for ray in self.complex().local_backend().cell_cone(cell).primitive_rays() if ray.q() == zero)
        else:
            cone = RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
            found: dict[tuple[int, ...], Lattices.ElementMethods] = {}
            for subcone in cone.triangulation().values():
                for mark in self._subcone_vectors_of_square(cone, subcone):
                    found[tuple(int(mark.to_vector()(label)) for label in cell.lattice.module_generating_set())] = mark
            marks = tuple(mark for mark in found.values() if mark.content() == 1)
        self._local_mark_cache[key] = marks
        return marks

    def _subcone_vectors_of_square(self, cone: Parent, subcone: Parent) -> tuple[Lattices.ElementMethods, ...]:
        r"""The lattice vectors of square n in one simplicial piece of a cell's cone.

        Let R be the matrix of the piece's primitive rays r_i and D = |det R|. A lattice
        vector of the piece is v = sum c_i r_i with C = D c in Z_{>=0}^k, and C mod D runs
        over the D residues of L / (sum Z r_i) (Smith form of R). The rays pair with one
        sign, so after that sign every pairing b_ij is nonnegative and
        D^2 q(v) = sum C_i C_j b_ij is monotone in each C_i. A timelike ray has
        C_i^2 b_ii <= D^2 |n|; a null ray needs a positively paired partner with C_j >= 1,
        so 2 C_i b_ij <= D^2 |n|. Each residue class is enumerated by depth-first search
        within these bounds, pruning a branch once its square exceeds |n|.
        """
        lattice = cone.ambient_lattice()
        # Search the timelike rays first: a null ray's coefficient raises
        # the square only once a positively paired coefficient is nonzero.
        rays = tuple(sorted(subcone.primitive_rays(), key=lambda ray: ray.q() == lattice.base_ring().zero()))
        if not rays:
            return ()
        labels = tuple(lattice.module_generating_set())
        rank = len(rays)
        if rank != len(labels):
            raise ArithmeticError("a simplicial piece of a full-dimensional cell cone must be full-dimensional")
        pairings = [[int(lattice.b(left, right)) for right in rays] for left in rays]
        # A simplicial face contained in a null ray has identically zero
        # quadratic form. It contributes no positive-norm marks, regardless
        # of the ambient lattice's signature. In particular, choosing the
        # sign from the zero pairing matrix would be mathematically invalid.
        if all(entry == 0 for row in pairings for entry in row):
            return ()
        sign = 1 if sum(sum(row) for row in pairings) > 0 else -1
        signed = [[sign * entry for entry in row] for row in pairings]
        if any(entry < 0 for row in signed for entry in row):
            raise ArithmeticError(f"the rays {rays} of a cell's cone do not pair with one sign, as rays of the closed positive cone do")
        target = sign * int(self._norm)
        if target <= 0:
            raise ArithmeticError("positive-square marks require a nonzero norm with the sign of the cone")
        ray_matrix = matrix(SageZZ, [[int(ray.to_vector()(label)) for label in labels] for ray in rays])
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
                    raise ArithmeticError(f"the isotropic ray {rays[position]} has no positive pairing in its simplicial piece")
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
        found: list[Lattices.ElementMethods] = []
        ring = lattice.base_ring()

        def search(scaled: list[int], paired: list[int], square: int, position: int) -> None:
            if position == rank:
                if square == scaled_target:
                    row = vector(SageZZ, scaled) * ray_matrix
                    if any(entry % determinant for entry in row):
                        raise ArithmeticError("a residue-class combination of the rays is not a lattice vector")
                    found.append(lattice.linear_combination({label: ring(int(entry) // determinant) for label, entry in zip(labels, row, strict=True) if entry}))
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

    def transport_marks(self, adjacency: LorentzianCellAdjacency) -> tuple[tuple[Lattices.ElementMethods, Lattices.ElementMethods], ...]:
        transporter = adjacency.transporter
        if transporter is None:
            raise ArithmeticError("a marked quotient adjacency needs a transporter")
        source = adjacency.source.vector_configuration
        facet = RationalPolyhedralCones(adjacency.source.lattice).from_rays(tuple(source[position] for position in sorted(adjacency.facet)))
        source_marks = tuple(mark for mark in self.local_marks(adjacency.source) if facet.contains(mark))
        target_marks = self.local_marks(adjacency.target)
        labels = adjacency.source.lattice.module_generating_set()
        source_by_coordinates = {tuple(int(mark.to_vector()(label)) for label in labels): mark for mark in source_marks}
        if len(source_by_coordinates) != len(source_marks):
            raise ArithmeticError("a facet contains duplicate marks")
        transported = []
        for target_mark in target_marks:
            source_mark = transporter(target_mark)
            key = tuple(int(source_mark.to_vector()(label)) for label in labels)
            if key in source_by_coordinates:
                transported.append((source_by_coordinates[key], target_mark))
        if len(transported) != len(source_marks):
            raise ArithmeticError("a facet mark has no unique mark in the target representative")
        return tuple(transported)

    def global_orbits(self) -> tuple[tuple[Lattices.ElementMethods, ...], ...]:
        cells = self.complex().quotient_cells()
        local = tuple(self.local_stabilizer_orbits(cell) for cell in cells)
        nodes = [(cell_position, orbit_position) for cell_position, cell_orbits in enumerate(local) for orbit_position in range(len(cell_orbits))]
        graph: dict[tuple[int, int], set[tuple[int, int]]] = {node: set() for node in nodes}
        local_mark_positions = tuple(
            {
                tuple(int(mark.to_vector()(label)) for label in cell.lattice.module_generating_set()): orbit_position
                for orbit_position, orbit in enumerate(cell_orbits)
                for mark in orbit
            }
            for cell, cell_orbits in zip(cells, local, strict=True)
        )

        def local_orbit_position(cell_position: int, mark: Lattices.ElementMethods) -> int:
            labels = cells[cell_position].lattice.module_generating_set()
            key = tuple(int(mark.to_vector()(label)) for label in labels)
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
                    labels = cells[cell_position].lattice.module_generating_set()
                    key = tuple(int(mark.to_vector()(label)) for label in labels)
                    if key not in seen_marks:
                        seen_marks.add(key)
                        marks.append(mark)
            global_orbits.append(tuple(marks))
        return tuple(global_orbits)


def _coordinate_row[S: RingElement](lattice: FramedFreeModules.ParentMethods, element: FramedFreeModules.ElementMethods, ring: Parent[S]) -> FreeModuleElement[S]:
    coordinates = element.to_vector()
    # A preamble scalar belongs to its own coefficient ring, not Sage's ZZ/QQ
    # parents.  Cross the private engine boundary by exact rational components;
    # Sage's coercion discovery cannot convert those categorical scalars.
    return vector(ring, tuple(ring(_rational(coordinates(label))) for label in lattice.module_generating_set()))


def _negative_ldl_factorization(gram: Matrix_rational_dense) -> tuple[list[list[Rational]], list[Rational]]:
    """Exact ``G = L D L^T`` factorization of a positive-definite Gram matrix."""
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
            raise ArithmeticError("the fixed-distance sphere requires a negative-definite lattice")
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
    target_row = cast("FreeModuleElement[Rational]", _coordinate_row(lattice, target, SageQQ))
    inverse = basis_inverse if basis_inverse is not None else cast(Matrix_rational_dense, basis_change.inverse())
    target_column = cast(Matrix_rational_dense, inverse * target_row.column())
    reduced_target = cast("FreeModuleElement[Rational]", vector(cast("Parent[Rational]", SageQQ), tuple(target_column[row, 0] for row in range(int(basis_change.nrows())))))
    if basis_change * target_column != target_row.column():
        raise ArithmeticError("the reduced fixed-distance center was not transported by the basis isomorphism")
    gram = cast(Matrix_rational_dense, -_gram_matrix(reduced).change_ring(SageQQ))
    rank = int(gram.nrows())
    zero = cast(Rational, SageQQ.zero())
    lower, diagonal = ldl_factorization if ldl_factorization is not None else _negative_ldl_factorization(gram)
    if len(lower) != rank or len(diagonal) != rank:
        raise ArithmeticError("the fixed-distance LDL factorization has incorrect rank")
    if any(len(row) != rank for row in lower):
        raise ArithmeticError("the fixed-distance LDL factorization has an invalid lower-triangular matrix")

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


def _owned_rational(lattice: Lattices.ParentMethods, value: Rational) -> RingElement:
    field = lattice.base_ring().fraction_field()
    return field(int(value.numerator())) / field(int(value.denominator()))


def _integral_positive_vector_allmeth(gram: Matrix_integer_dense) -> tuple[Integer, ...]:
    r"""Port of the V2 diagonal/eigen ApproxIterator branch in Positivity.h.

    IndefiniteReduction supplies the source's unimodular reduced Gram matrix.
    Exact rational diagonal directions are interleaved with numerical
    positive-eigenvalue directions as in GetIntegralVector_allmeth_V2. The
    exact diagonal family guarantees termination; every candidate is checked
    against the exact integral Gram matrix before acceptance.
    """
    if gram.nrows() == 0:
        raise ValueError("the Gram matrix has no positive direction")
    basis, reduced = _indefinite_reduction(gram)
    # SignatureSymmetric.h:142-200, NullspaceReduction followed by
    # DiagonalizeNonDegenerateSymmetricMatrix. Store basis vectors as rows.
    rank = int(reduced.nrows())
    qgram = cast(Matrix_rational_dense, reduced.change_ring(SageQQ))
    radical = cast(Matrix_rational_dense, qgram.left_kernel().basis_matrix())
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
    pivot_columns = set(radical.echelon_form().pivots())
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
            orthogonal_basis = cast(
                Matrix_rational_dense,
                cast(Matrix_rational_dense, nondegenerate * previous.transpose()).left_kernel().basis_matrix(),
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
    if any(diagonal_matrix[i, j] for i in range(rank) for j in range(rank) if i != j):
        raise ArithmeticError("pinned symmetric diagonalization left off-diagonal entries")
    if transform * reduced.change_ring(SageQQ) * transform.transpose() != diagonal_matrix:
        raise ArithmeticError("positive-direction diagonalization lost its congruence witness")
    if not transform.det():
        raise ArithmeticError("positive-direction diagonalization is singular")

    def make_approx_iterator(
        direction: tuple[Rational, ...] | tuple[float, ...],
    ) -> tuple[
        tuple[tuple[Rational | float, tuple[tuple[int, int], ...]], ...],
        list[Rational | float],
        list[int],
    ]:
        groups: dict[Rational | float, list[tuple[int, int]]] = {}
        for coordinate, entry in enumerate(direction):
            if entry:
                if isinstance(entry, float):
                    step: Rational | float = 1.0 / abs(entry)
                    sign = 1 if entry > 0 else -1
                else:
                    rational_entry = cast(Rational, entry)
                    step = SageQQ.one() / abs(rational_entry)
                    sign = 1 if rational_entry > 0 else -1
                groups.setdefault(step, []).append((coordinate, sign))
        ordered = tuple((step, tuple(group)) for step, group in sorted(groups.items()))
        nexts = [step / 2.0 if isinstance(step, float) else cast(Rational, step) / SageQQ(2) for step, _group in ordered]
        return ordered, nexts, [0] * rank

    iterators: list[
        tuple[
            tuple[tuple[Rational | float, tuple[tuple[int, int], ...]], ...],
            list[Rational | float],
            list[int],
        ]
    ] = []
    for index in range(rank):
        if diagonal_matrix[index, index] <= 0:
            continue
        row = cast("FreeModuleElement[Rational]", transform.row(index))
        norm = SageQQ.zero()
        for entry in row:
            norm += abs(entry)
        normalized = tuple(SageQQ(entry) / norm for entry in row)
        iterators.append(make_approx_iterator(normalized))
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
        eigenvector_data = cast(
            list[tuple[ComplexDoubleElement, list[FreeModuleElement[ComplexDoubleElement]], int]],
            double_gram.eigenvectors_right(),
        )
        for eigenvalue, eigenvectors, _multiplicity in eigenvector_data:
            if float(eigenvalue.real()) > 0:
                for eigenvector in eigenvectors:
                    direction = tuple(float(entry.real()) for entry in eigenvector)
                    eigen_norm = sum(abs(eigen_entry) for eigen_entry in direction)
                    if eigen_norm:
                        eigen_iterators.append(make_approx_iterator(tuple(eigen_entry / eigen_norm for eigen_entry in direction)))
    except ArithmeticError, OverflowError, RuntimeError, TypeError, ValueError:
        # This branch only accelerates the exact search. Overflow or failure
        # in the floating eigensolver must not disable its exact directions.
        eigen_iterators.clear()
    if not iterators:
        raise ValueError("the Gram matrix has no positive direction")
    while True:
        for ordered, nexts, approximation in iterators:
            least = min(nexts)
            for group_index, (step, group) in enumerate(ordered):
                if nexts[group_index] != least:
                    continue
                nexts[group_index] = cast(Rational, nexts[group_index]) + cast(Rational, step)
                for coordinate, sign in group:
                    approximation[coordinate] += sign
            candidate = vector(SageZZ, approximation)
            if (candidate * reduced * candidate.column())[0] > 0:
                original = basis.transpose() * candidate.column()
                result = tuple(SageZZ(original[row, 0]) for row in range(int(original.nrows())))
                original_vector = vector(SageZZ, result)
                original_norm = (original_vector * gram * original_vector.column())[0]
                reduced_norm = (candidate * reduced * candidate.column())[0]
                if original_norm != reduced_norm or original_norm <= 0:
                    raise ArithmeticError("positive-vector reduction lost its integral norm witness")
                return result
        for ordered, nexts, approximation in eigen_iterators:
            least = min(nexts)
            for group_index, (step, group) in enumerate(ordered):
                if nexts[group_index] != least:
                    continue
                nexts[group_index] = cast(float, nexts[group_index]) + cast(float, step)
                for coordinate, sign in group:
                    approximation[coordinate] += sign
            candidate = vector(SageZZ, approximation)
            if (candidate * reduced * candidate.column())[0] > 0:
                original = basis.transpose() * candidate.column()
                result = tuple(SageZZ(original[row, 0]) for row in range(int(original.nrows())))
                original_vector = vector(SageZZ, result)
                original_norm = (original_vector * gram * original_vector.column())[0]
                reduced_norm = (candidate * reduced * candidate.column())[0]
                if original_norm != reduced_norm or original_norm <= 0:
                    raise ArithmeticError("positive-vector reduction lost its integral norm witness")
                return result


def _positive_direction(lattice: Lattices.ParentMethods) -> tuple[Rational, ...]:
    r"""Pinned polyhedral_common@a55fcb7b src_isotropy/Positivity.h:565-673,
    ``INDEFINITE_GetShortPositiveVector``.

    The positive-vector finder is the independent preamble primitive. The
    Graver-basis descent, randomized unimodular reframing, 500-failure stop,
    and short-norm stop are the reference's *specialized* search, not a
    replacement one-shot reduction.
    """
    gram = _gram_matrix(lattice)
    rank = gram.nrows()
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

    perturbation = identity_matrix(SageZZ, rank)
    best_vector: FreeModuleElement[Integer] | None = None
    best_norm: Integer | None = None
    without_improvement = 0
    gram_l1 = sum(abs(entry) for entry in gram.list())
    ring = lattice.base_ring()
    module = lattice.unformed_module()
    bilinear_forms = module.bilinear_hom_space(ring)
    while True:
        perturbed_gram = perturbation * gram * perturbation.transpose()
        if perturbation == identity_matrix(SageZZ, rank):
            positive = lattice.positive_vector()
            found = _coordinate_row(lattice, positive, SageZZ)
        else:
            # The unimodular change of frame is applied to the *form* before
            # the general positive-vector primitive, exactly as upstream.
            form = bilinear_forms([[ring(int(perturbed_gram[i, j])) for j in range(rank)] for i in range(rank)])
            framed = FormModules(ring)(form)
            found = _coordinate_row(framed, framed.positive_vector(), SageZZ)
        current = perturbation.transpose() * found.column()
        improved = current.column(0)
        while True:
            changes = 0
            for direction in directions:
                current_norm = int((improved * gram * improved.column())[0])
                mixed = int((improved * gram * direction.column())[0])
                direction_norm = int((direction * gram * direction.column())[0])
                alpha = _quadratic_descent_prefix(current_norm, mixed, direction_norm)
                if alpha:
                    improved += alpha * direction
                    changes += alpha
            if changes == 0:
                break
        norm = (improved * gram * improved.column())[0]
        if best_norm is None or norm < best_norm:
            best_vector = improved
            best_norm = SageZZ(norm)
            without_improvement = 0
        else:
            without_improvement += 1
        if without_improvement == 500 or (best_norm is not None and best_norm < 10):
            if best_vector is None:
                raise ArithmeticError("the positive-vector search has no best vector")
            return tuple(SageQQ(entry) for entry in best_vector)
        # Lorentzian.g:1090-1093 resets the random frame only at 10000 times
        # the original entrywise L1 norm, not at 1000.
        if sum(abs(entry) for entry in perturbed_gram.list()) > 10000 * gram_l1:
            perturbation = identity_matrix(SageZZ, rank)
        modification = identity_matrix(SageZZ, rank)
        if choice((0, 1)) == 0:
            permutation = list(range(rank))
            shuffle(permutation)
            modification = matrix(SageZZ, rank, rank)
            for index, image in enumerate(permutation):
                modification[index, image] = choice((-1, 1))
        else:
            modification[randrange(rank), randrange(rank)] = choice((-1, 1))
        perturbation *= modification


def _primitive_integral_direction(lattice: Lattices.ParentMethods, coordinates: tuple[Rational, ...]) -> tuple[Lattices.ElementMethods, Rational]:
    r"""The primitive vector t of L on the ray of a rational direction u, and the scale s with u = t / s."""
    if len(coordinates) != int(lattice.module_rank()):
        raise ValueError("the rational direction must have the lattice rank")
    if not any(coordinates):
        raise ValueError("the rational direction must be nonzero")
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
        self.timelike = timelike
        self.timelike_row = _coordinate_row(lattice, timelike, SageZZ)
        self.timelike_covector = self.gram * self.timelike_row.column()
        self.scale = scale
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
        timelike_row = self.timelike_row
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

    @cached_method
    def vectors_at_multiplier(self, multiplier: int) -> tuple[Lattices.ElementMethods, ...]:
        target = self.base_target.parent().scalar_multiple(_owned_rational(self.lattice, SageQQ(multiplier)), self.base_target)
        bound = _owned_rational(self.lattice, -(SageQQ(multiplier) ** 2) * self.base_bound)
        # Preserve the pinned Lorentzian.g:369-373 closest-vector API here.
        # The specialized exact-sphere recursion is tested independently, but
        # replacing this reference path changed performance on the frozen
        # perfect-domain orbit-count fixtures.
        distances = self.kernel.close_vectors(target, bound)
        return self._ambient_vectors(multiplier, tuple(distances.index_set()))

    def _ambient_vectors(self, multiplier: int, close_vectors: tuple[Lattices.ElementMethods, ...]) -> tuple[Lattices.ElementMethods, ...]:
        shift = self.lattice.scalar_multiple(self.lattice.base_ring()(multiplier), self.partner)
        result: list[Lattices.ElementMethods] = []
        for close in close_vectors:
            candidate = shift + self.inclusion(close)
            row = _coordinate_row(self.lattice, candidate, SageZZ)
            candidate_square = int(row * self.gram * row)
            # Lorentzian.g:377-385 verifies both the scalar-product shell
            # and norm after reconstructing each closest-vector solution.
            pairing = (row * self.timelike_covector)[0]
            if pairing != multiplier * self.d or pairing <= 0 or pairing > self.scaled_max:
                raise ArithmeticError("a close-vector shell returned a vector outside its prescribed pairing fibre")
            if candidate_square < 0:
                raise ArithmeticError("a close-vector shell returned a negative-norm vector")
            # Lorentzian.g:382-388: the isotropic option keeps fNorm = 0 among
            # the closest-vector solutions, which lie within the distance.
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
    r"""The first shell along the direction, following LORENTZ_SearchInitialVectorOpt's successive integral bounds."""
    if allowed_squares == ZERO_SQUARES and not lattice.is_isotropic():
        raise ValueError(f"{lattice} is anisotropic: it has no isotropic perfect cells")
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    direction_row = vector(SageQQ, direction)
    max_scal = (direction_row * gram * direction_row.column())[0]
    enumerator = _PositiveVectorEnumerator(lattice, direction, max_scal, allowed_squares)
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
            vectors = enumerator.vectors_at_multiplier(multiplier)
            if vectors:
                return vectors
        next_multiplier = max(next_multiplier, maximum + 1)
        # Lorentzian.g:411–423 increases MaxScal by one, rather than
        # exponentially expanding the search interval.
        max_scal += 1
        enumerator.scaled_max += enumerator.scale


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
    basis = matrix(SageQQ, [base_vector.column(0), direction_vector.column(0)])
    restricted = basis * gram * basis.transpose()
    a, b, c = restricted[0, 0], restricted[0, 1], restricted[1, 1]
    discriminant = b * b - a * c
    isotropic_bound: Rational | None = None
    if discriminant > 0 and SageQQ(discriminant).is_square():
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
    gram = _gram_matrix(lattice).change_ring(SageQQ)
    perturbation = identity_matrix(SageZZ, spanning_rows.nrows())
    while True:
        rows = perturbation * spanning_rows
        restricted = rows * gram * rows.transpose()
        framed = Lattices(rational.base_ring())([[_owned_rational(lattice, restricted[i, j]) for j in range(restricted.ncols())] for i in range(restricted.nrows())])
        negative = framed.vector_of_sign(-1)
        if negative is None:
            raise ArithmeticError("the outside-ray span has no negative direction")
        coordinates = negative.to_vector()
        negative_row = vector(SageQQ, [_rational(coordinates(label)) for label in framed.module_generating_set()])
        direction = negative_row * rows
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
        # Lorentzian.g:920-923 / lorentzian_perfect.h:757: a shell
        # containing every critical point cannot have fewer distinct
        # points than the critical family.
        if critical_rows.issubset(total_rows) and len(critical) > len(total):
            raise ArithmeticError("a Lorentzian shell lost multiplicity of its critical configuration")
        if total_rows == critical_rows:
            lower = middle
            continue
        # Lorentzian.g:924-930, EXIT 1: the shell strictly contains CritSet.
        if critical_rows.issubset(total_rows):
            return total, normal, test_direction_row, max_scal
        break
    while True:
        if not total:
            raise ArithmeticError("perfect-cell flipping reached an empty positive-vector shell")
        # Lorentzian.g:944-966, second scheme.
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
        # Lorentzian.g:959, EXIT 2: the shell contains CritSet and eVect.
        if (critical_rows | _vector_rows((new_vector,))).issubset(_vector_rows(total)):
            return total, normal, test_direction_row, max_scal
