"""Allcock's edgewalk on a hyperbolic lattice, over preamble operations.

A translation of the walk in the pinned polyhedral_common edgewalk source. The walk along the
edges of the fundamental polyhedron of W(L), the cusp handling and the termination certificate
are this module's; Vinberg's chamber, cones, binary-form representations, affine lattice lines,
orders of induced actions and orbit closures are the preamble's. Edge planes are computed by
private rational linear algebra on coordinate rows, and every result is rebuilt as a lattice
element before it reaches a preamble operation.

Conventions: L has signature (n, 1); roots have positive norm; a vertex is a primitive k with
q(k) <= 0, and the walls through it are the roots r with b(r, k) = 0, all others pairing negatively.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key

from dzack_research.preamble.categories.sets.indexed_families import finite_indexed_family_from_values
from dzack_research.preamble.all import ZZ, HyperbolicLattices, Lattices
from dzack_research.preamble.categories.group.g_sets import FiniteGSets
from dzack_research.preamble.categories.hyperbolic_lattices import AllcockEdgewalkReport
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones
from dzack_research.preamble.categories.vector_configurations import (
    VectorConfigurations,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.edgewalk_extensions import (
    NormedDynkinExtension,
    finite_type_extensions,
    parabolic_extensions,
)
from sage_indefinite_port.indefinite.edgewalk_rank2 import anisotropic_cycle

type _Vector = Lattices.ElementMethods
type _Roots = tuple[_Vector, ...]
type _Key = tuple[int, ...]
type _RationalRow = FreeModuleElement[Rational]


def edgewalk_fundamental_domain(
    lattice: HyperbolicLattices.ParentMethods,
) -> AllcockEdgewalkReport:
    r"""Allcock's edgewalk on a hyperbolic lattice.

    Returns the preamble's ``AllcockEdgewalkReport`` of ``lattice``: the simple roots of a
    fundamental polyhedron of the reflection group W(L), as elements of ``lattice``, and whether
    ``lattice`` is reflective.
    """
    return _Edgewalk(lattice).walk()


def _key(vector_: _Vector) -> _Key:
    coordinates = vector_.to_vector()
    return tuple(
        int(coordinates(label)) for label in vector_.parent().module_generating_set()
    )


def _oriented_determinant(left: _Key, right: _Key) -> int:
    return left[0] * right[1] - left[1] * right[0]


def _orientation_compare(left: _Key, right: _Key) -> int:
    determinant = _oriented_determinant(left, right)
    if determinant > 0:
        return -1
    if determinant < 0:
        return 1
    return 0


def _signed_sqrt_compare(
    sign_left: int, square_left: Rational, sign_right: int, square_right: Rational
) -> int:
    if sign_left != sign_right:
        return -1 if sign_left < sign_right else 1
    if sign_left == 0:
        return 0
    if square_left == square_right:
        return 0
    if sign_left > 0:
        return -1 if square_left < square_right else 1
    return -1 if square_left > square_right else 1


@dataclass(frozen=True)
class _Vertex:
    generator: _Vector
    roots: _Roots


@dataclass(frozen=True)
class _Candidate:
    alpha: _Vector
    residual_norm: Rational
    root_norm: Rational
    vertex: _Vertex


@dataclass(frozen=True)
class _EdgePlane:
    r"""The rank-two lattice of an edge: the projection Pi of the root sublattice C onto the edge plane (span of the edge's roots)^perp.

    ``plane`` is Pi with its form scaled by ``scale`` to be integral, in the oriented basis ``basis``
    (rational rows of L tensor Q); ``r0`` is the edge's orientation vector in Pi; ``intersection``
    is the embedded sublattice C cap Pi whose quotient controls the period.
    """

    constraint: Parent
    basis: Matrix_rational_dense
    plane: Lattices.ParentMethods
    scale: Integer
    r0: _Vector
    intersection: Parent


class _Edgewalk:
    def __init__(self, lattice: HyperbolicLattices.ParentMethods) -> None:
        positive, negative = lattice.signature_pair()
        if int(negative) != 1 or int(positive) + 1 != int(lattice.rank()):
            raise ValueError(
                f"Allcock's edgewalk needs a nondegenerate lattice of signature (n, 1), but {lattice} has signature {lattice.signature_pair()}"
            )
        self.lattice = lattice
        self.rank = int(lattice.rank())
        self.labels = tuple(lattice.module_generating_set())
        generators = tuple(lattice.module_generators())
        self.gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageQQ)
        self.rational = lattice.vector_space()
        self.field = lattice.base_ring().fraction_field()
        self.norms = tuple(int(norm) for norm in lattice.possible_root_lengths())
        self.identity_key = tuple(_key(generator) for generator in generators)

    def element(self, row: _RationalRow | _Key) -> _Vector:
        ring = self.lattice.base_ring()
        return self.lattice.linear_combination(
            {
                label: ring(int(entry))
                for label, entry in zip(self.labels, row, strict=True)
                if entry
            }
        )

    def rational_element(self, row: _RationalRow) -> Lattices.ElementMethods:
        return self.rational.linear_combination(
            {
                label: self.field(int(SageQQ(entry).numerator()))
                / self.field(int(SageQQ(entry).denominator()))
                for label, entry in zip(self.labels, row, strict=True)
                if entry
            }
        )

    def row(self, vector_: _Vector) -> _RationalRow:
        return vector(SageQQ, _key(vector_))

    def rational_row(self, vector_: Lattices.ElementMethods) -> _RationalRow:
        coordinates = vector_.to_vector()
        return vector(
            SageQQ,
            [
                SageQQ(int(coordinates(label).numerator()))
                / SageQQ(int(coordinates(label).denominator()))
                for label in self.labels
            ],
        )

    def pair(self, left: _RationalRow, right: _RationalRow) -> Rational:
        r"""b on rational rows, for points of L tensor Q that are not lattice vectors."""
        return (left * self.gram * right.column())[0]

    def pairing(self, left: _Vector, right: _Vector) -> Rational:
        return self.pair(self.row(left), self.row(right))

    def norm(self, point: _Vector) -> Rational:
        row = self.row(point)
        return self.pair(row, row)

    def primitive(self, row: _RationalRow) -> _Vector:
        r"""The primitive vector of L on the ray of a nonzero rational row."""
        denominator = SageZZ.one()
        for entry in row:
            denominator = denominator.lcm(SageQQ(entry).denominator())
        scaled = self.element(tuple(SageZZ(denominator * entry) for entry in row))
        content = int(scaled.content())
        return self.element(tuple(entry // content for entry in _key(scaled)))

    def span_rank(self, vectors: _Roots) -> int:
        return (
            matrix(SageQQ, [self.row(item) for item in vectors]).rank()
            if vectors
            else 0
        )

    def perpendicular_rows(
        self, vectors: tuple[_RationalRow, ...]
    ) -> Matrix_rational_dense:
        return (
            matrix(
                SageQQ,
                len(vectors),
                self.rank,
                [entry for row in vectors for entry in row],
            )
            * self.gram
        ).right_kernel_matrix()

    def initial_isotropic_vertex(self, k: _Vector | None = None) -> _Vertex | None:
        r"""The ``isotropic_vinberg`` starting corner of pinned edgewalk.h:1880–2055.

        The source finds a primitive null vector, enumerates each admissible
        root-norm shell in the projection of C_N ∩ k-perp to k-perp/k,
        selects one simple cone, lifts its walls and completes the affine
        cusp roots. The preamble supplies the isotropic reduction and exact
        definite-lattice fixed-square enumeration.
        """
        if k is None:
            k = self.lattice.isotropic_vector()
        reduction = k.isotropic_reduction()
        perpendicular = reduction.orthogonal_complement()
        perpendicular_inclusion = perpendicular.inclusion()
        quotient_projection = reduction.projection()
        rank = int(reduction.rank())
        projected_roots: list[tuple[_Vector, _Vector]] = []

        for norm in self.norms:
            constraint = self.constraint(SageQQ(norm))
            intersection = constraint.intersection(perpendicular)
            intersection_inclusion = intersection.inclusion()

            projected_generators = {
                label: quotient_projection(
                    perpendicular_inclusion.lift(
                        intersection_inclusion(intersection.module_generator(label))
                    )
                )
                for label in intersection.module_generating_set()
            }
            projection = intersection.module_category().Mor(intersection, reduction)(
                projected_generators
            )
            image = reduction.subobject_on(tuple(projected_generators.values()))
            for root in image.vectors_of_square(reduction.base_ring()(norm)):
                projected = image.inclusion()(root)
                lifted = intersection_inclusion(projection.lift(projected))
                if self.norm(lifted) != norm:
                    raise ArithmeticError(
                        "the isotropic-reduction lift changed a prescribed root norm"
                    )
                projected_roots.append((projected, lifted))

        if not projected_roots:
            # get_simple_cone: an empty projected root system certifies
            # non-reflectivity, rather than an arithmetic failure.
            return None
        labels = tuple(reduction.module_generating_set())

        def quotient_row(point: _Vector) -> _RationalRow:
            coordinates = point.to_vector()
            return vector(SageQQ, [SageQQ(int(coordinates(label))) for label in labels])

        projected_rows = tuple(
            quotient_row(projected) for projected, _root in projected_roots
        )
        if matrix(SageQQ, projected_rows).rank() != rank:
            # edgewalk.h: the projected roots must span the rank-(n-2)
            # quotient for a reflective cusp.
            return None

        # GetFacetOneDomain_ListIdx (lorentzian_linalg.h:267–339):
        # regular splitting vector, positive half and its extreme rays.
        # A regular splitting is outside the finitely many hyperplanes
        # annihilating the nonzero projected roots. On the moment curve
        # (1,t,...,t^(rank-1)), each forbidden hyperplane has at most
        # rank-1 integer parameters t. Exhausting t=0,1,... therefore
        # constructs a regular splitting without probabilistic retries.
        if any(not any(row) for row in projected_rows):
            raise ArithmeticError("the projected root list contains a zero vector")
        parameter = 0
        while True:
            splitting = vector(SageQQ, [SageQQ(parameter) ** degree for degree in range(rank)])
            if all(splitting.dot_product(row) != 0 for row in projected_rows):
                break
            parameter += 1
        positive = [
            (projected, root)
            for (projected, root), row in zip(
                projected_roots, projected_rows, strict=True
            )
            if splitting.dot_product(row) > 0
        ]
        cone = RationalPolyhedralCones(reduction).from_rays(
            [projected for projected, _root in positive]
        )
        simple_roots: list[_Vector] = []
        simple_projected_rows: list[_RationalRow] = []
        for ray in cone.primitive_rays():
            ray_row = quotient_row(ray)
            for projected, root in positive:
                root_row = quotient_row(projected)
                if (
                    matrix(SageQQ, [ray_row, root_row]).rank() == 1
                    and reduction.b(ray, projected) > 0
                ):
                    simple_roots.append(root)
                    simple_projected_rows.append(root_row)
                    break
            else:
                raise ArithmeticError(
                    "a simple cone ray has no original admissible root lift"
                )
        # Source get_simple_cone accepts non-simplicial faces. Test rank in
        # the isotropic quotient, not in the ambient lattice of their lifts.
        if matrix(SageQQ, simple_projected_rows).rank() != rank:
            raise ArithmeticError(
                "the isotropic starting root cone does not span the projected lattice"
            )

        # edgewalk.h:2020–2049: choose a negative kP orthogonal to the
        # simple walls, orient by b(kP,k)<0 and complete the cusp by the
        # source's affine-line equation fibres.
        plane_rows = self.perpendicular_rows(
            tuple(self.row(root) for root in simple_roots)
        )
        plane = self.rational.vector_configuration(finite_indexed_family_from_values(None, 
            [self.rational_element(row) for row in plane_rows.rows()]
        ))
        negative = self.rational_row(plane.inclusion()(plane.vector_of_sign(-1)))
        if self.pair(negative, self.row(k)) > 0:
            negative = -negative
        previous = self.primitive(negative)
        return _Vertex(k, self.cusp_roots(tuple(simple_roots), k, previous))

    def initial_vertex(self) -> _Vertex | None:
        r"""Find a corner, first checking the finite root system at a timelike point.

        A full-rank perpendicular root system makes the point a corner of
        some chamber. Otherwise Vinberg's increasing prefixes locate a corner;
        no finite search bound is used as a nonexistence certificate.
        """
        if self.lattice.is_isotropic():
            return self.initial_isotropic_vertex()
        timelike = self.lattice.negative_vector()
        roots = self.ordinary_vertex_roots(timelike)
        if self.span_rank(roots) == self.rank - 1:
            return _Vertex(timelike, self.ordinary_vertex_roots(timelike, require_full=True))
        # vinberg.h FindRoots_Kernel / FindOneInitialRay: examine the accepted
        # root system after its initial finite cone and after each new root,
        # stopping on the first feasible nonpositive ray. The preamble's
        # VinAl backend supplies the actual root-decomposition enumeration.
        from vinal import GetIntegerPoints, VinAl

        engine = VinAl(self.gram, vector(SageZZ, _key(timelike)), False, None)
        # vinberg.h GetVinbergAux: enumerate translation classes of the
        # saturated perpendicular lattice together with v0, and reduce
        # each representative by the quotient of -<w,v0> and -v0^2.
        # VinAl sorts the unnormalized representatives, which changes the
        # source's decomposition sequence; reconstruct W before searching.
        coset_basis = matrix(SageZZ, tuple(engine.V1.gens()) + (engine.v0,))
        norm_shift = -engine.v0.inner_product(engine.v0)
        engine.W = []
        for coordinates in GetIntegerPoints(coset_basis):
            representative = engine.V(coordinates)
            quotient = (-representative.inner_product(engine.v0)) // norm_shift
            normalized = representative - quotient * engine.v0
            # GetVinbergAux uses QuoInt/ResInt to place each class in
            # the fundamental timelike strip, 0 <= -<w,v0> < -v0^2.
            height = -normalized.inner_product(engine.v0)
            if not 0 <= height < norm_shift:
                raise ArithmeticError("Vinberg coset representative leaves timelike strip")
            engine.W.append(normalized)
        # vinberg.h GetVinbergFromG supplies the caller's root-length set.
        # VinAl otherwise derives its own list from the form determinant,
        # which need not equal the admissible lengths of this lattice.
        engine.root_lengths = list(self.norms)
        # vinberg.h FundCone -> GetFacetOneDomain: select one member of
        # each perpendicular root pair by a regular linear functional,
        # and retain the original normed roots on the extreme cone rays.
        # Lexicographic sign is a deterministic regular functional for this
        # finite integral root family, unlike VinAl.FundCone's primitive
        # ray representatives.
        initial_roots = [
            root
            for norm in engine.root_lengths
            for root in engine.Roots_decomposed_into(engine.V([0] * engine.n), norm)
            if engine.IsRoot(root)
        ]
        # GetFacetOneDomain first reduces the root coordinate columns to
        # their actual span, and GetFacetOneDomain_ListIdx requires at least
        # twice that reduced dimension in the nonempty root family.
        if initial_roots:
            reduced_dimension = matrix(SageQQ, initial_roots).rank()
            if len(initial_roots) < 2 * reduced_dimension:
                raise ArithmeticError(
                    "initial Vinberg root family has fewer than twice its rank in roots"
                )
        positive = [
            root for root in initial_roots
            if next((entry for entry in root if entry), 0) > 0
        ]
        # GetFacetOneDomain_ListIdx allocates exactly half the roots on
        # the assumption that the family is centrally symmetric. A count
        # alone cannot establish the required opposite-root correspondence.
        if initial_roots:
            root_rows = [tuple(int(entry) for entry in root) for root in initial_roots]
            multiplicities = {}
            for row in root_rows:
                multiplicities[row] = multiplicities.get(row, 0) + 1
            if any(
                multiplicities.get(row, 0) != multiplicities.get(tuple(-entry for entry in row), 0)
                for row in multiplicities
            ):
                raise ArithmeticError("Vinberg perpendicular roots do not form opposite pairs")
        if positive:
            positive_cone = RationalPolyhedralCones(self.lattice).from_rays(
                self.element(tuple(int(entry) for entry in root)) for root in positive
            )
            rays = tuple(self.row(ray) for ray in positive_cone.primitive_rays())
            # GetFacetOneDomain_ListIdx selects nonredundant *indices*;
            # parallel root vectors cannot both survive its reduction.
            chosen = []
            used_rays = set()
            for root in positive:
                row = vector(SageQQ, root)
                for index, ray in enumerate(rays):
                    if (
                        index not in used_rays
                        and matrix(SageQQ, [row, ray]).rank() == 1
                        and row.dot_product(ray) > 0
                    ):
                        chosen.append(root)
                        used_rays.add(index)
                        break
            if len(used_rays) != len(rays):
                raise ArithmeticError("initial Vinberg cone ray has no normed root")
            engine.roots = chosen
        else:
            engine.roots = []

        def accepted_vertex() -> _Vertex | None:
            roots = tuple(self.element(tuple(int(entry) for entry in root)) for root in engine.roots)
            if not roots:
                return None
            equations = matrix(SageQQ, [self.row(root) * self.gram for root in roots])
            rank = equations.rank()
            if rank < self.rank - 1:
                return None
            if rank == self.rank - 1:
                generator = self.primitive(equations.right_kernel_matrix().row(0))
                if self.pairing(generator, timelike) > 0:
                    generator = -generator
                if self.norm(generator) > 0 or self.pairing(generator, timelike) >= 0:
                    return None
            else:
                # GetOneInteriorVertex: exact dual description of -R*G,
                # tested on each extreme ray in the selected time component.
                # GetOneInteriorVertex uses the covectors -G*r. The
                # lattice's algebraic correlation produces these directly,
                # without assuming that primal and dual generator labels agree.
                inequalities = tuple(-root.to_covector() for root in roots)
                cone = RationalPolyhedralCones(self.lattice).from_halfspaces(inequalities)
                generator = next((ray for ray in cone.primitive_rays()
                                  if self.norm(ray) <= 0 and self.pairing(ray, timelike) < 0), None)
                if generator is None:
                    return None
            if self.norm(generator) == 0:
                return self.initial_isotropic_vertex(generator)
            return _Vertex(generator, self.ordinary_vertex_roots(generator, require_full=True))

        initial = accepted_vertex()
        if initial is not None:
            return initial
        # vinberg.h FindRoots_Kernel tests its exit only after every root
        # belonging to one decomposition request has been inserted. VinAl's
        # NextRoot yields those roots individually; use the backend's called
        # decomposition helpers so an intermediate root cannot end the search.
        # vinberg.h IterateRootDecompositions::get_k (pinned 1592b246)
        # prioritizes each length's next coset by val/sqrt(val), with
        # val = -<v0, a>. VinAl instead uses val/sqrt(k), so its iterator
        # does not reproduce the C++ candidate sequence.
        candidates = {norm: 1 for norm in engine.root_lengths}
        while True:
            choices = {
                norm: engine.W[index % len(engine.W)]
                + (index // len(engine.W)) * engine.v0
                for norm, index in candidates.items()
            }
            norm = min(
                engine.root_lengths,
                key=lambda length: float(
                    -engine.v0.inner_product(choices[length])
                ) ** 0.5,
            )
            component = choices[norm]
            candidates[norm] += 1
            # vinberg.h FindRoot_filter/fct_CVP retains every norm-k
            # candidate satisfying all accepted wall inequalities. VinAl's
            # IsNewRoot additionally rejects linearly dependent candidates;
            # that is not the upstream filter and can discard actual walls.
            new_roots = [
                root for root in engine.Roots_decomposed_into(component, norm)
                if engine.IsRoot(root)
                and all(root.inner_product(wall) <= 0 for wall in engine.roots)
            ]
            if not new_roots:
                continue
            engine.roots.extend(new_roots)
            # vinberg.h ReduceListRoot removes generators redundant in the
            # accepted root cone. Keep the original normed root vectors,
            # using primitive rays only to identify their directions.
            ambient_roots = tuple(
                self.element(tuple(int(entry) for entry in root))
                for root in engine.roots
            )
            extreme = RationalPolyhedralCones(self.lattice).from_rays(
                ambient_roots
            ).primitive_rays()
            # ReduceListRoot/get_non_redundant_indices preserves the order
            # of the original family. Keep its first normed representative
            # of each extremal direction rather than reordering by rays.
            ray_rows = tuple(self.row(ray) for ray in extreme)
            reduced = []
            retained_directions = set()
            for root in engine.roots:
                row = vector(SageQQ, root)
                for index, ray_row in enumerate(ray_rows):
                    if (
                        index not in retained_directions
                        and matrix(SageQQ, [row, ray_row]).rank() == 1
                        and row.dot_product(ray_row) > 0
                    ):
                        reduced.append(root)
                        retained_directions.add(index)
                        break
            if len(retained_directions) != len(ray_rows):
                raise ArithmeticError("Vinberg cone ray has no accepted normed root")
            engine.roots = reduced
            initial = accepted_vertex()
            if initial is not None:
                return initial

    def ordinary_vertex_roots(self, timelike: _Vector, *, require_full: bool = False) -> _Roots:
        """A simple system of all ambient reflective roots perpendicular to k.

        Lexicographic positivity is a regular ordering on this finite root
        system. Its positive-root cone has precisely the simple-root rays.
        """
        inclusion = timelike.orthogonal_complement().inclusion()
        complement = inclusion.domain()
        positive: list[_Vector] = []
        for norm in self.norms:
            for point in complement.vectors_of_square(complement.base_ring()(norm)):
                root = inclusion(point)
                if next(entry for entry in _key(root) if entry) < 0:
                    continue
                if all(
                    (2 * int(self.pairing(root, basis))) % norm == 0
                    for basis in self.lattice.module_generators()
                ):
                    positive.append(root)
        if not positive:
            if require_full:
                # The upstream ordinary-cone constructor rejects a root
                # family whose facet count is not dimSpace (including zero).
                raise ArithmeticError("ordinary initial cone has incorrect simple-root count")
            return ()
        # edgewalk.h get_simple_cone_from_lattice selects actual normed
        # roots from FindFixedDistVectors, not primitive multiples of them.
        # Use geometric rays only to identify which original roots span
        # the selected cone, retaining each root's prescribed norm.
        rays = RationalPolyhedralCones(self.lattice).from_rays(positive).primitive_rays()
        if require_full and len(rays) != self.rank - 1:
            # edgewalk.h get_simple_cone_from_lattice requires exactly
            # dimSpace selected walls for an ordinary starting vertex.
            raise ArithmeticError("ordinary initial cone has incorrect simple-root count")
        return tuple(
            next(
                root for root in positive
                if matrix(SageQQ, [self.row(root), self.row(ray)]).rank() == 1
                and self.row(root).dot_product(self.row(ray)) > 0
            )
            for ray in rays
        )

    def constraint(self, norm: Rational) -> Parent:
        r"""C = {v in L : 2 b(v, L) / norm in Z}, the vectors whose reflection of that norm is integral."""
        integral_norm = int(norm)
        return self.constraint_for_divisor(
            integral_norm // 2 if integral_norm % 2 == 0 else integral_norm
        )

    @cached_method
    def constraint_for_divisor(self, divisor: int) -> Parent:
        return self.lattice.divisible_sublattice(divisor)

    def edge_plane(
        self, roots: _Roots, discarded: _Vector, k: _Vector, norm: Rational
    ) -> _EdgePlane:
        integral_norm = int(norm)
        return self.edge_plane_for_divisor(
            roots,
            discarded,
            k,
            integral_norm // 2 if integral_norm % 2 == 0 else integral_norm,
        )

    @cached_method
    def edge_plane_for_divisor(
        self, roots: _Roots, discarded: _Vector, k: _Vector, divisor: int
    ) -> _EdgePlane:
        root_rows = tuple(self.row(root) for root in roots)
        plane_basis = self.perpendicular_rows(root_rows)
        if plane_basis.nrows() != 2:
            raise ArithmeticError(
                "an edge must have a two-dimensional perpendicular plane"
            )
        r0_space = self.perpendicular_rows((*root_rows, self.row(k)))
        if r0_space.nrows() != 1:
            raise ArithmeticError("the edge orientation vector must be one-dimensional")
        r0 = r0_space.row(0)
        if self.pair(r0, self.row(discarded)) > 0:
            r0 = -r0
        if self.norm(k) < 0:
            oriented_basis = matrix(SageQQ, [r0, self.row(k)])
        else:
            plane = self.rational.vector_configuration(finite_indexed_family_from_values(None, 
                [self.rational_element(row) for row in plane_basis.rows()]
            ))
            positive = self.rational_row(plane.inclusion()(plane.vector_of_sign(1)))
            if self.pair(positive, self.row(k)) > 0:
                positive = -positive
            oriented_basis = matrix(SageQQ, [positive, self.row(k)])
            r0 = -self.row(k)

        constraint = self.constraint_for_divisor(divisor)
        constraint_rows = matrix(
            SageQQ,
            [
                self.row(constraint.inclusion()(generator))
                for generator in constraint.inclusion().domain().module_generators()
            ],
        )
        root_matrix = matrix(
            SageQQ,
            len(root_rows),
            self.rank,
            [entry for row in root_rows for entry in row],
        )
        root_gram = root_matrix * self.gram * root_matrix.transpose()
        projected = [
            row
            - (row * self.gram * root_matrix.transpose())
            * root_gram.inverse()
            * root_matrix
            for row in constraint_rows.rows()
        ]
        basis = matrix(
            SageQQ,
            list((SageQQ**self.rank).span(projected, SageZZ).basis_matrix().rows()),
        )
        if basis.nrows() != 2:
            raise ArithmeticError("the projected constraint lattice must have rank two")
        if (
            oriented_basis.transpose().solve_right(basis.transpose()).transpose().det()
            < 0
        ):
            basis.rescale_row(0, -1)
        work_gram = basis * self.gram * basis.transpose()
        scale = SageZZ.one()
        for entry in work_gram.list():
            scale = scale.lcm(entry.denominator())
        plane_lattice = Lattices(ZZ)(
            [
                [ZZ(int(scale * work_gram[row, column])) for column in range(2)]
                for row in range(2)
            ]
        )
        r0_coordinates = basis.transpose().solve_right(r0.column()).column(0)
        r0_vector = self.plane_primitive(plane_lattice, r0_coordinates)

        perpendicular = self.lattice.vector_configuration(finite_indexed_family_from_values(None, 
            list(roots)
        )).orthogonal_complement()
        # For divisor one C=L, so C cap P is already the perpendicular
        # sublattice. There is no second integral intersection to compute.
        intersection = (
            perpendicular if divisor == 1 else constraint.intersection(perpendicular)
        )
        intersection_rows = matrix(
            SageQQ,
            [
                self.row(intersection.inclusion()(generator))
                for generator in intersection.inclusion().domain().module_generators()
            ],
        )
        intersection_coordinates = (
            basis.transpose().solve_right(intersection_rows.transpose()).transpose()
        )
        if any(entry.denominator() != 1 for entry in intersection_coordinates.list()):
            raise ArithmeticError(
                "the perpendicular intersection must be integral in the projected lattice"
            )
        sublattice = plane_lattice.vector_configuration(finite_indexed_family_from_values(None, 
            [
                self.plane_element(plane_lattice, row)
                for row in intersection_coordinates.rows()
            ]
        ))
        return _EdgePlane(
            constraint,
            basis,
            plane_lattice,
            scale,
            r0_vector,
            sublattice,
        )

    @staticmethod
    def plane_element(
        plane: Lattices.ParentMethods, coordinates: _RationalRow
    ) -> _Vector:
        ring = plane.base_ring()
        return plane.linear_combination(
            {
                label: ring(int(entry))
                for label, entry in zip(
                    plane.module_generating_set(), coordinates, strict=True
                )
                if entry
            }
        )

    @staticmethod
    def plane_primitive(
        plane: Lattices.ParentMethods, coordinates: _RationalRow
    ) -> _Vector:
        r"""The primitive vector of the plane lattice on the ray of nonzero rational coordinates."""
        denominator = SageZZ.one()
        for entry in coordinates:
            denominator = denominator.lcm(SageQQ(entry).denominator())
        scaled = _Edgewalk.plane_element(
            plane, vector(SageQQ, [denominator * entry for entry in coordinates])
        )
        content = int(scaled.content())
        return _Edgewalk.plane_element(
            plane, vector(SageQQ, [entry // content for entry in _key(scaled)])
        )

    @cached_method
    def binary_candidates(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        norm: Rational,
        residual: Rational,
    ) -> _Roots:
        if residual <= 0:
            return ()
        edge = self.edge_plane(roots, discarded, k, norm)
        plane, scale, r0 = edge.plane, edge.scale, edge.r0
        work_vectors: list[_Vector] = []
        if plane.is_isotropic():
            square = scale * residual
            if square.denominator() != 1:
                return ()
            vectors = [
                item
                for item in plane.vectors_of_square(ZZ(int(square)))
                if _oriented_determinant(_key(r0), _key(item)) > 0
            ]
            if self.norm(k) < 0:
                vectors = [item for item in vectors if plane.b(r0, item) > 0]
            vectors.sort(key=lambda item: cmp_to_key(_orientation_compare)(_key(item)))
            work_vectors.extend(vectors)
        else:
            cycle = anisotropic_cycle(plane, scale * SageQQ(norm), r0)
            if cycle is None:
                return ()
            automorph, representatives = cycle
            multiples: list[_Vector] = []
            for item in representatives:
                square = SageQQ(int(item.q()))
                multiplier = 1
                while multiplier * multiplier * square <= scale * norm:
                    candidate = plane.scalar_multiple(
                        plane.base_ring()(multiplier), item
                    )
                    if SageQQ(int(candidate.q())) == scale * residual:
                        multiples.append(candidate)
                    multiplier += 1
            # lorentzian_linalg.h:70-143,
            # GetMatrixExponentSublattice_TrivClass: first find the least
            # positive power preserving C cap Pi, then the order of its
            # action on Pi/(C cap Pi). The source tests integral membership
            # in the sublattice and then the induced finite quotient action,
            # not the generally larger order in GL(2, Z/exp(Pi/S)).
            embedding = edge.intersection.inclusion()
            period = automorph
            preserving_order = 1
            while not all(
                embedding.is_in_image(period(embedding(point)))
                for point in edge.intersection.module_generators()
            ):
                period = automorph * period
                preserving_order += 1
            quotient_order = 1
            quotient_power = period
            while not all(
                embedding.is_in_image(quotient_power(point) - point)
                for point in plane.module_generators()
            ):
                quotient_power = period * quotient_power
                quotient_order += 1
            period_order = preserving_order * quotient_order
            work_vectors.extend(multiples)
            for _ in range(1, period_order):
                multiples = [automorph(item) for item in multiples]
                work_vectors.extend(multiples)

        return tuple(work_vectors)

    def extension_roots(
        self,
        roots: _Roots,
        discarded: _Vector,
        k: _Vector,
        extension: NormedDynkinExtension,
    ) -> _Roots:
        edge = self.edge_plane(roots, discarded, k, extension.norm)
        work_vectors = self.binary_candidates(
            roots, discarded, k, extension.norm, extension.residual_norm
        )
        component = vector(SageQQ, extension.u_component)
        answer: dict[_Key, _Vector] = {}
        for work_vector in work_vectors:
            ambient = component + edge.basis.transpose() * vector(
                SageQQ, _key(work_vector)
            )
            if any(entry.denominator() != 1 for entry in ambient):
                continue
            root = self.element(tuple(ambient))
            # edgewalk.h get_next_anisotropic/get_next_isotropic: every
            # integral candidate must belong to the norm-specific lattice.
            if not edge.constraint.inclusion().is_in_image(root):
                continue
            if self.pair(ambient, ambient) != extension.norm:
                raise ArithmeticError("rank-two reconstruction changed the root norm")
            answer.setdefault(_key(root), root)
        return tuple(answer.values())

    def candidate_compare(self, left: _Candidate, right: _Candidate, k: _Vector) -> int:
        left_scal = -SageQQ(int(self.pairing(left.alpha, k)))
        right_scal = -SageQQ(int(self.pairing(right.alpha, k)))
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

    def vertex_from_root(
        self, k: _Vector, roots: _Roots, discarded: _Vector, alpha: _Vector
    ) -> _Vertex | None:
        extended = (*roots, alpha)
        kernel = self.perpendicular_rows(tuple(self.row(root) for root in extended))
        if kernel.nrows() != 1:
            return None
        generator = self.primitive(kernel.row(0))
        old_pairing = self.pairing(generator, k)
        if old_pairing > 0:
            generator = -generator
            old_pairing = -old_pairing
        if old_pairing == 0 or self.pairing(discarded, generator) >= 0:
            return None
        # edgewalk.h:1032-1058 stores alpha itself in MatRoot. A scalar
        # multiple defines the same hyperplane, but not the same normed root
        # configuration; replacing alpha by its primitive part here changes
        # the reference algorithm's subsequent extensions and stabilizers.
        return _Vertex(generator, (*roots, alpha))

    def cusp_root(
        self, constraint: Parent, component: tuple[Rational, ...], k: _Vector
    ) -> _Vector | None:
        r"""The point u + c k of C on the line through the component u in direction k with the least c > 0."""
        points = constraint.affine_line_points(
            self.rational_element(vector(SageQQ, component)), k
        )
        if points is None:
            return None
        base, step = (constraint.inclusion()(point) for point in points)
        position = next(index for index, entry in enumerate(_key(k)) if entry)

        def coefficient(row: _RationalRow) -> Rational:
            return row[position] / _key(k)[position]

        offset = coefficient(self.row(base) - vector(SageQQ, component))
        increment = coefficient(self.row(step))
        threshold = -offset / increment
        shift = threshold.floor() + 1 if increment > 0 else threshold.ceil() - 1
        root = base + self.lattice.scalar_multiple(
            self.lattice.base_ring()(int(shift)), step
        )
        if offset + shift * increment <= 0:
            raise ArithmeticError("the cusp equation must return c > 0")
        return root

    def cusp_roots(self, roots: _Roots, k: _Vector, previous: _Vector) -> _Roots:
        candidates: list[tuple[int, Rational, Rational, _Vector]] = []
        for extension in parabolic_extensions(self.lattice, roots, self.norms):
            root = self.cusp_root(
                self.constraint(extension.norm), extension.u_component, k
            )
            if root is None:
                continue
            scalar = -SageQQ(int(self.pairing(previous, root)))
            sign = (scalar > 0) - (scalar < 0)
            candidates.append(
                (sign, scalar * scalar / extension.norm, extension.norm, root)
            )

        def compare(
            left: tuple[int, Rational, Rational, _Vector],
            right: tuple[int, Rational, Rational, _Vector],
        ) -> int:
            first = _signed_sqrt_compare(left[0], left[1], right[0], right[1])
            if first:
                return first
            if left[2] == right[2]:
                return 0
            return -1 if left[2] < right[2] else 1

        candidates.sort(key=cmp_to_key(compare))
        accepted = list(roots)
        for _sign, _quant, _norm, root in candidates:
            if all(self.pairing(old, root) <= 0 for old in accepted):
                # edgewalk.h:393-408, DetermineRootsCuspidalCase, appends the
                # actual chosen vector, without primitive normalization.
                accepted.append(root)
        return tuple(accepted)

    def isotropic_next_vertex(
        self, roots: _Roots, discarded: _Vector, k: _Vector
    ) -> _Vertex:
        perpendicular = self.lattice.vector_configuration(finite_indexed_family_from_values(None, 
            list(roots)
        )).orthogonal_complement()
        plane = perpendicular.inclusion().domain()
        choices: dict[_Key, _Vector] = {}
        for isotropic in plane.primitive_isotropic_vectors():
            candidate = perpendicular.inclusion()(isotropic)
            if (
                self.pairing(candidate, k) < 0
                and self.pairing(discarded, candidate) < 0
            ):
                choices.setdefault(_key(candidate), candidate)
        if len(choices) != 1:
            raise ArithmeticError(
                "the edge must determine one forward isotropic endpoint"
            )
        (generator,) = choices.values()
        return _Vertex(generator, self.cusp_roots(roots, generator, k))

    @cached_method
    def edge_step(self, vertex: _Vertex, roots: _Roots, discarded: _Vector) -> _Vertex:
        candidates: list[_Candidate] = []
        for extension in finite_type_extensions(self.lattice, roots, self.norms):
            # edgewalk.h:966-1061 calls get_next(extension) once: it takes
            # the first integral candidate in the binary search ordering.
            # Only that candidate is tested as a possible next vertex.
            chosen = self.extension_roots(roots, discarded, vertex.generator, extension)
            if chosen:
                alpha = chosen[0]
                next_vertex = self.vertex_from_root(
                    vertex.generator, roots, discarded, alpha
                )
                if next_vertex is not None:
                    candidates.append(
                        _Candidate(
                            alpha, extension.residual_norm, extension.norm, next_vertex
                        )
                    )
        if candidates:
            candidates.sort(
                key=cmp_to_key(
                    lambda left, right: self.candidate_compare(
                        left, right, vertex.generator
                    )
                )
            )
            return candidates[0].vertex
        return self.isotropic_next_vertex(roots, discarded, vertex.generator)

    @cached_method
    def edge_directions(self, vertex: _Vertex) -> tuple[tuple[_Roots, _Vector], ...]:
        r"""The walls along each outgoing edge and a wall orienting that edge."""
        if len(vertex.roots) == 1:
            return ()
        cone = RationalPolyhedralCones(self.lattice).from_rays(vertex.roots)
        directions: list[tuple[_Roots, _Vector]] = []
        for covector in cone.facet_covectors():
            included = tuple(
                root
                for root in vertex.roots
                if cone.evaluate_covector(covector, root) == 0
            )
            omitted = tuple(
                root
                for root in vertex.roots
                if cone.evaluate_covector(covector, root) != 0
            )
            if not omitted:
                raise ArithmeticError("a root-cone facet must omit a wall")
            # An ideal corner with several affine components can omit more
            # than one wall. Any omitted wall fixes the same forward side.
            directions.append((included, omitted[0]))
        return tuple(directions)

    @cached_method
    def vertex_configuration(self, vertex: _Vertex) -> tuple[_Roots, tuple[int, ...]]:
        r"""The coloured configuration of the vertex: its roots, the vertex, and at an ideal vertex its neighbours."""
        vectors = [*vertex.roots, vertex.generator]
        colours = [1] * len(vertex.roots) + [2]
        if self.norm(vertex.generator) == 0:
            adjacent: dict[_Key, _Vector] = {}
            for roots, discarded in self.edge_directions(vertex):
                neighbour = self.edge_step(vertex, roots, discarded).generator
                adjacent.setdefault(_key(neighbour), neighbour)
            vectors.extend(adjacent[key] for key in sorted(adjacent))
            colours.extend([3] * len(adjacent))
        return tuple(vectors), tuple(colours)

    def isometry_key(self, isometry: LatticeIsometryMethods) -> tuple[_Key, ...]:
        return tuple(
            _key(isometry(generator)) for generator in self.lattice.module_generators()
        )

    @staticmethod
    def common_fixed_space(
        invariant: Parent, isometry: LatticeIsometryMethods
    ) -> Parent | None:
        r"""The sublattice fixed by the group generated so far, or None once that group is infinite.

        A finite cone-preserving group fixes a timelike vector by averaging.
        Conversely its integral stabilizer is discrete in a compact real
        orthogonal group, hence finite. The common fixed lattice decides both
        directions, without separate order computations for the generators.
        """
        fixed = invariant.intersection(isometry.invariant_lattice())
        _positive, negative = fixed.inclusion().domain().signature_pair()
        return None if int(negative) == 0 else fixed

    def root_chamber_certificate(self, roots: _Roots) -> bool:
        r"""Upstream LORENTZ_RunEdgewalkAlgorithm dual-description ray test.

        The inequalities are -b(root, x) >= 0.  A pointed rational
        chamber is tested on its exact extreme rays, as in the pinned
        upstream dual-description termination condition.
        """
        if not roots:
            return False
        inequalities = tuple(-root.to_covector() for root in roots)
        chamber = RationalPolyhedralCones(self.lattice).from_halfspaces(inequalities)
        if not chamber.is_pointed():
            return False
        return all(self.norm(ray) <= 0 for ray in chamber.primitive_rays())

    def report(
        self,
        vertices: list[_Vertex],
        reflective: bool,
        simple_roots: _Roots,
        isometries: list[LatticeIsometryMethods],
    ) -> AllcockEdgewalkReport:
        cones = RationalPolyhedralCones(self.lattice)
        return AllcockEdgewalkReport(
            self.lattice,
            reflective=reflective,
            simple_roots=simple_roots,
            vertices=tuple(
                cones.from_rays([vertex.generator], wall_roots=vertex.roots)
                for vertex in vertices
            ),
            isometry_generators=tuple(isometries),
        )

    def edge_orbits(
        self, vertex: _Vertex, local_generators: tuple[LatticeIsometryMethods, ...]
    ) -> tuple[tuple[_Roots, _Vector], ...]:
        """Outgoing edges modulo the supplied full vertex stabilizer."""
        directions = self.edge_directions(vertex)
        direction_positions = {
            frozenset(_key(root) for root in roots): index
            for index, (roots, _discarded) in enumerate(directions)
        }

        def act_on_direction(isometry: LatticeIsometryMethods, index: int) -> int:
            return direction_positions[
                frozenset(_key(isometry(root)) for root in directions[index][0])
            ]

        local_group = self.lattice.O().subgroup(local_generators)
        direction_orbits = FiniteGSets(local_group)(
            tuple(range(len(directions))), act_on_direction
        ).orbits()
        return tuple(
            directions[int(orbit.representative())] for orbit in direction_orbits
        )

    def binary_walk(self) -> AllcockEdgewalkReport:
        """A binary chamber has two ordinary ends, or one wall in the split case."""
        root_locus = self.lattice.reflective_root_locus()
        # In the rationally split binary case, a nonzero norm fibre is
        # finite: after factoring the form over Q, its two integral linear
        # factors have a prescribed nonzero product. Enumerate that actual
        # fibre. In the anisotropic case it can be infinite, and the
        # preamble's finite orbit-representative computation is required.
        if self.lattice.is_isotropic():
            candidates = (
                point
                for norm in self.norms
                for point in self.lattice.vectors_of_square(ZZ(norm))
            )
        else:
            candidates = (
                point
                for norm in self.norms
                for point in self.lattice.binary_fixed_norm_representatives(ZZ(norm))
            )
        root = next((point for point in candidates if point in root_locus), None)
        if root is None:
            # An empty root system has no reflecting walls. In particular,
            # rational isotropy does not certify reflectivity: the trivial
            # reflection group has infinite-volume chamber in H^1.
            return self.report([], False, (), [])
        generator = self.primitive(self.perpendicular_rows((self.row(root),)).row(0))
        initial = _Vertex(generator, (root,))
        if self.lattice.is_isotropic():
            return self.report([initial], True, initial.roots, [])
        neighbor = self.edge_step(initial, (), root)
        category = VectorConfigurations(self.lattice.base_ring())
        symmetry = category.spanning_isometry(
            self.lattice,
            (root, generator),
            (1, 2),
            self.lattice,
            (*neighbor.roots, neighbor.generator),
            (1, 2),
        )
        return self.report(
            [initial, neighbor],
            True,
            (*initial.roots, *neighbor.roots),
            [] if symmetry is None else [symmetry],
        )

    def walk(self) -> AllcockEdgewalkReport:
        if self.rank == 2:
            return self.binary_walk()
        initial = self.initial_vertex()
        if initial is None:
            return self.report([], False, (), [])

        vertices = [initial]
        configurations = [self.vertex_configuration(initial)]
        seen: set[tuple[_Key, ...]] = {self.identity_key}
        isometries: list[LatticeIsometryMethods] = []
        invariant = self.lattice.vector_configuration(finite_indexed_family_from_values(None, 
            list(self.lattice.module_generators())
        ))

        def admit(isometry: LatticeIsometryMethods) -> bool:
            r"""Record a new symmetry of the polyhedron; False once the group they generate is infinite."""
            nonlocal invariant
            key = self.isometry_key(isometry)
            if key in seen:
                return True
            seen.add(key)
            isometries.append(isometry)
            fixed = self.common_fixed_space(invariant, isometry)
            if fixed is None:
                return False
            invariant = fixed
            return True

        def certified_roots() -> _Roots:
            # LORENTZ_RunEdgewalkAlgorithm::f_try_terminate: close the
            # accumulated walls under every discovered chamber symmetry.
            group = self.lattice.O().subgroup(tuple(isometries))
            seeds = tuple({
                _key(root): root for vertex in vertices for root in vertex.roots
            }.values())
            closure = FiniteGSets(group).from_finite_orbits(
                seeds, lambda isometry, root: isometry(root)
            )
            return tuple(sorted(closure, key=_key))

        certified_root_count = 0

        def try_reflectivity_certificate() -> AllcockEdgewalkReport | None:
            # LORENTZ_RunEdgewalkAlgorithm::f_try_terminate tests the dual
            # description only when the closed root orbit has grown.
            nonlocal certified_root_count
            roots = certified_roots()
            if len(roots) <= certified_root_count:
                return None
            certified_root_count = len(roots)
            if self.root_chamber_certificate(roots):
                return self.report(vertices, True, roots, isometries)
            return None

        # The source checks after insertion of the initial vertex and after
        # every newly discovered vertex, not only once the quotient closes.
        certificate = try_reflectivity_certificate()
        if certificate is not None:
            return certificate

        position = 0
        configuration_category = VectorConfigurations(self.lattice.base_ring())
        while position < len(vertices):
            vertex = vertices[position]
            configuration = configurations[position]
            local_generators = tuple(
                configuration_category.spanning_automorphism_generators(
                    self.lattice, *configuration
                )
            )
            for automorphism in local_generators:
                previous_generators = len(isometries)
                if not admit(automorphism):
                    return self.report(vertices, False, (), isometries)
                if len(isometries) > previous_generators:
                    certificate = try_reflectivity_certificate()
                    if certificate is not None:
                        return certificate
            for roots, discarded in self.edge_orbits(vertex, local_generators):
                neighbour = self.edge_step(vertex, roots, discarded)
                neighbour_configuration = self.vertex_configuration(neighbour)
                isometries_to_neighbour = (
                    configuration_category.spanning_isometry(
                        self.lattice, *known, self.lattice, *neighbour_configuration
                    )
                    for known in configurations
                )
                transporter = next(
                    (
                        isometry
                        for isometry in isometries_to_neighbour
                        if isometry is not None
                    ),
                    None,
                )
                if transporter is None:
                    vertices.append(neighbour)
                    configurations.append(neighbour_configuration)
                    certificate = try_reflectivity_certificate()
                    if certificate is not None:
                        return certificate
                else:
                    previous_generators = len(isometries)
                    if not admit(transporter):
                        return self.report(vertices, False, (), isometries)
                    if len(isometries) > previous_generators:
                        certificate = try_reflectivity_certificate()
                        if certificate is not None:
                            return certificate
            position += 1

        simple_roots = certified_roots()
        return self.report(vertices, self.root_chamber_certificate(simple_roots), simple_roots, isometries)


__all__ = ["edgewalk_fundamental_domain"]
