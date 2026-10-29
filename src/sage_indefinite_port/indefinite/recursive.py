r"""Dutour Sikirić's Eichler recursion over preamble operations (map §10, references/port-design.md).

For L of positive index at least 2 after the sign normalization of AttackProfile:

- generation step: if the A_L-translates of finitely many vectors r_i = t_i(v) cover the
  O(L)-orbit of a vector v, then O(L) = <A_L, Aut((L, v)), t_i>;
- Aut((L, v)): for q(v) != 0, the extensions across v^perp of generators of O(v^perp), then the
  subgroup of the group they generate preserving L; for q(v) = 0, the pointwise stabilizer of
  Zv in its parabolic subgroup;
- an element of Isom((L, v), (M, w)): for q(v) != 0, an isometry of the complements by
  recursion, its extension, and its correction by an element of the group above carrying the
  image of L onto M; for q(v) = 0, an integral lift of an isometry of the reductions, searched
  over the cosets of the Levi image;
- isometry(L, M): a splitting vector v of L, the covering vectors w of M of the same square
  and content, and an element of Isom((L, v), (M, w));
- primitive vector classes: covering vectors deduplicated by deciding Isom exactly.

Definite lattices are the preamble's; signature (1, n) is the Lorentzian leaf (#15).
"""

from dataclasses import dataclass

from dzack_research.preamble.all import ZZ, Lattices
from dzack_research.preamble.categories.group.groups import OwnedGroups
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.orbit_spans import (
    OrbitSpanConstructions,
)
from dzack_research.preamble.categories.modules.pure.modules import (
    Modules,
    ModulesOverIntegralDomains,
)
from dzack_research.preamble.categories.rational_integral_stabilizers import (
    IntegralStructureAction,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.rings.integer_ring import ZZ as SageZZ
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from sage_indefinite_port.indefinite.groupoid_cache import IsometryGroupoidCache
from sage_indefinite_port.indefinite.indefinite_reduction import (
    indefinite_reduced_presentation,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    lorentzian_isometry,
    lorentzian_orthogonal_group_generators,
)
from sage_indefinite_port.invariants import (
    AttackProfile,
    LatticePrefilter,
    VectorPrefilter,
)

_UNTWIST = Lattices(ZZ).twist_functor(-1)


def _untwisted(isometry: LatticeIsometryMethods, profile: AttackProfile) -> LatticeIsometryMethods:
    r"""An isometry of the signed view, as an isometry of the lattice itself."""
    return isometry if profile.sign == 1 else _UNTWIST(isometry)


def _extend_perpendicular_isometry(
    source_reduction,
    target_reduction,
    source_basis: tuple[Lattices.ElementMethods, ...],
    images: tuple[Lattices.ElementMethods, ...],
) -> tuple[LatticeIsometryMethods, int]:
    r"""Extend a supplied isometry of I^perp using the preamble Witt operation.

    source: Lorentzian.g:4-50 LORENTZ_ExtendOrthogonalIsotropicIsomorphism_Dim1
    (rank-one I, called from IndefiniteForms.g:1029 and :1638) and
    Lorentzian.g:115-249 LORENTZ_ExtendOrthogonalIsotropicIsomorphism
    get_one_transformation (rank-two I, IndefiniteForms.g:1097-1106, :1685-1686).
    The source's Subspace1/Subspace2 are a basis of N=I^perp and its images;
    here they are ``source_basis`` and ``images``. Its extension to L_Q -> M_Q, the Euclidean-nullspace
    complement (Fundamental.g:509-516), the pairing equations
    (SpecialEquationSolving, Lorentzian.g:53-99) and the prime-denominator
    elimination (Fundamental.g:751-794) belong to the preamble's
    ``witt_extension``. For rank-one I the extension is unique (an isometry
    fixing the hyperplane v^perp pointwise of a nondegenerate space with
    isotropic normal is trivial), so the choice of complement does not
    change the result. ``GetDenominatorMatrix`` (IndefiniteForms.g:1105) is
    the denominator ideal of the extension from L to M.
    """
    source = source_reduction.isotropic_embedding().codomain()
    target = target_reduction.isotropic_embedding().codomain()
    fraction_map = source.base_ring().fraction_field_map()
    field = fraction_map.codomain()
    rank = int(source_reduction.orthogonal_complement().inclusion().domain().module_rank())
    if len(source_basis) != rank or len(images) != rank:
        raise ValueError("the perpendicular partial isometry needs images of a full chosen basis")
    rational_source = source.base_change(fraction_map)
    rational_target = target.base_change(fraction_map)
    source_rationalize = ModulesOverIntegralDomains(source.base_ring())(source).generic_fibre_map()
    target_rationalize = ModulesOverIntegralDomains(target.base_ring())(target).generic_fibre_map()
    # The source's Subspace1 and Subspace2 are the rows of the two families.
    # Here they are the maps D -> L_Q and D -> M_Q out of the rational
    # lattice D on the basis, with the Gram pulled back along the first; the
    # partial isometry is the identity of D. No coordinate solve is needed.
    frame = field.free_module(rank)
    labels = tuple(frame.module_generating_set())
    gram = rational_source.gram_tensor().pullback(
        Modules(field).Mor(frame, rational_source)(dict(zip(labels, (source_rationalize(point).underlying_element() for point in source_basis), strict=True)))
    )
    span = Lattices(field)([[gram[row, column] for column in range(rank)] for row in range(rank)])
    span_labels = tuple(span.module_generating_set())
    source_inclusion = Modules(field).Mor(span, rational_source)(
        dict(zip(span_labels, (source_rationalize(point).underlying_element() for point in source_basis), strict=True))
    )
    target_inclusion = Modules(field).Mor(span, rational_target)(dict(zip(span_labels, (target_rationalize(point).underlying_element() for point in images), strict=True)))
    partial = span.Isom(span).identity()
    extension = partial.witt_extension(source_inclusion, target_inclusion)
    (denominator,) = tuple(extension.denominator_ideal(source, target).ideal_generators())
    return extension, int(denominator)


def _column_matrix(isometry: LatticeIsometryMethods) -> Matrix_integer_dense:
    r"""The integer matrix whose column ``j`` holds the coordinates of the image of generator ``j``."""
    domain, codomain = isometry.domain(), isometry.codomain()
    codomain_labels = tuple(codomain.module_generating_set())
    domain_labels = tuple(domain.module_generating_set())
    domain_basis = tuple(domain.module_generator(label) for label in domain_labels)
    images = tuple(isometry(generator).to_vector() for generator in domain_basis)
    return matrix(SageZZ, [[int(image(label)) for image in images] for label in codomain_labels])


def _from_column(lattice: Lattices.ParentMethods, column) -> Lattices.ElementMethods:
    r"""The element of ``lattice`` with the integer coordinates ``column``."""
    return lattice.linear_combination({label: int(entry) for label, entry in zip(lattice.module_generating_set(), column, strict=True) if entry})


def _conjugated(transport: LatticeIsometryMethods, generators: tuple[LatticeIsometryMethods, ...]) -> tuple[LatticeIsometryMethods, ...]:
    r"""The generators ``h`` of O(R) carried to L by ``transport: R -> L``, each built once from ``P H P^-1``.

    source: IndefiniteForms.g:1487 ``Inverse(RecRed.B) * eGen1 * RecRed.B``.
    Each conjugate is constructed from the images of the generators of L, so
    later evaluations do not compose the three maps again.
    """
    lattice = transport.codomain()
    change = _column_matrix(transport)
    inverse = change.inverse().change_ring(SageZZ)
    group = lattice.O()
    return tuple(group(tuple(_from_column(lattice, column) for column in (change * _column_matrix(generator) * inverse).columns())) for generator in generators)


def _transported(target_transport: LatticeIsometryMethods, witness: LatticeIsometryMethods, source_transport: LatticeIsometryMethods) -> LatticeIsometryMethods:
    r"""``target_transport * witness * source_transport^-1``, built once from its column matrix.

    source: IndefiniteForms.g:1389 ``Inverse(RecRed2.B) * test * RecRed1.B`` and
    :1419 (the degenerate extension ``Inverse(FullBasis2) * TheEquivTest * FullBasis1``).
    """
    source, target = source_transport.codomain(), target_transport.codomain()
    columns = _column_matrix(target_transport) * _column_matrix(witness) * _column_matrix(source_transport).inverse().change_ring(SageZZ)
    return source.Isom(target)(tuple(_from_column(target, column) for column in columns.columns()))


@dataclass(frozen=True)
class VectorClass:
    r"""An O(L)-orbit of primitive vectors: a representative, and the isometries of L found carrying it to the other covering vectors of its orbit."""

    representative: Lattices.ElementMethods
    witnesses: tuple[LatticeIsometryMethods, ...]


class EichlerRecursion:
    r"""The recursion; every general step is a preamble operation."""

    def __init__(self) -> None:
        # IndefiniteForms.g:1173-1310: the groupoid database is shared by
        # every recursive orthogonal-group and lattice-equivalence request.
        self._groupoid = IsometryGroupoidCache()
        self._models: dict[Lattices.ParentMethods, OrbitCoverModel] = {}

    def _degenerate_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""IndefiniteForms.g:853–901, 1496–1511: extend O(C) over the radical.

        The presentation of ``lattice`` is ``C ⊕ R`` with Gram ``C ⊕ 0``.
        Every integral isometry has block rows ``[A B; 0 D]``, where
        ``A`` preserves C, ``D`` belongs to GL(R), and B is arbitrary.
        The source extends the generators of O(C) and GL(R) and adjoining
        the first-radical-coordinate shears e_i -> e_i + z_1: the GL(R)
        action supplies the other coordinates of B.
        """
        gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
        size, nonsingular_rank = gram.nrows(), gram.rank()
        assert nonsingular_rank < size
        if any(gram[i, j] for i in range(nonsingular_rank, size) for j in range(size)):
            raise ArithmeticError("the radical must occupy the final Gram coordinates")
        matrices = []
        if nonsingular_rank:
            complement_gram = gram.matrix_from_rows_and_columns(range(nonsingular_rank), range(nonsingular_rank))
            complement = lattice.lattice_category()([[int(entry) for entry in row] for row in complement_gram.rows()])
            complement_labels = tuple(complement.module_generating_set())
            complement_basis = tuple(complement.module_generator(label) for label in complement_labels)
            for generator in self.orthogonal_group_generators(complement):
                # `gram` is expressed in module_generating_set() order. Build
                # the extended block in that same order rather than assuming
                # module_generators() enumerates the identical basis.
                top = tuple(
                    tuple(int(generator(point).to_vector()(label)) for label in complement_labels)
                    for point in complement_basis
                )
                block = identity_matrix(SageZZ, size)
                for row in range(nonsingular_rank):
                    for column in range(nonsingular_rank):
                        block[row, column] = top[row][column]
                matrices.append(block)
        radical_rank = size - nonsingular_rank
        general_linear = OwnedGroups().GL(radical_rank, lattice.base_ring())
        for generator in general_linear.group_generators():
            action = general_linear._to_engine(generator).matrix()
            block = identity_matrix(SageZZ, size)
            for i in range(radical_rank):
                for j in range(radical_rank):
                    block[nonsingular_rank + i, nonsingular_rank + j] = SageZZ(action[i, j])
            matrices.append(block)
        for row in range(nonsingular_rank):
            block = identity_matrix(SageZZ, size)
            block[row, nonsingular_rank] = 1
            matrices.append(block)
        return tuple(lattice.O()(tuple(_from_column(lattice, row) for row in block.rows())) for block in matrices)

    def _degenerate_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        r"""IndefiniteForms.g:1387–1423: split radicals, identify them, and extend the nondegenerate witness."""
        source_reduction = indefinite_reduced_presentation(source)
        target_reduction = indefinite_reduced_presentation(target)
        left, right = source_reduction.reduced, target_reduction.reduced
        left_gram = _engine_component_matrix(left.gram_tensor()).change_ring(SageZZ)
        right_gram = _engine_component_matrix(right.gram_tensor()).change_ring(SageZZ)
        size, rank = left_gram.nrows(), left_gram.rank()
        if rank != right_gram.rank():
            return None
        known, cached = self._groupoid.cached_isometry(left, right)
        if known:
            return None if cached is None else _transported(target_reduction.isometry, cached, source_reduction.isometry)
        extension = identity_matrix(SageZZ, size)
        if rank:
            left_complement_gram = left_gram.matrix_from_rows_and_columns(range(rank), range(rank))
            right_complement_gram = right_gram.matrix_from_rows_and_columns(range(rank), range(rank))
            left_complement = left.lattice_category()([[int(entry) for entry in row] for row in left_complement_gram.rows()])
            right_complement = right.lattice_category()([[int(entry) for entry in row] for row in right_complement_gram.rows()])
            complement_witness = self.isometry(left_complement, right_complement)
            if complement_witness is None:
                self._groupoid.record_isometry(left, right, None)
                return None
            left_labels = tuple(left_complement.module_generating_set())
            right_labels = tuple(right_complement.module_generating_set())
            left_basis = tuple(left_complement.module_generator(label) for label in left_labels)
            for row, point in enumerate(left_basis):
                image_coordinates = complement_witness(point).to_vector()
                for column, label in enumerate(right_labels):
                    extension[row, column] = int(image_coordinates(label))
        reduced_witness = left.Isom(right)(tuple(_from_column(right, row) for row in extension.rows()))
        self._groupoid.record_isometry(left, right, reduced_witness)
        return _transported(target_reduction.isometry, reduced_witness, source_reduction.isometry)

    def orthogonal_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        if lattice.is_definite():
            # IndefiniteForms.g:1275-1310 INDEF_FORM_AutomorphismGroup_Memoize
            # covers the definite kernel (:1435) too: a group computed
            # for an isometric lattice is transported along the recorded
            # isometry instead of being recomputed.
            cached = self._groupoid.cached_group(lattice)
            if cached is not None:
                return cached
            definite_generators = tuple(lattice.O().group_generators())
            self._groupoid.record_group(lattice, definite_generators)
            return definite_generators
        # IndefiniteForms.g:1474-1494: reduce before the memoized group
        # computation, then conjugate its generators back to the input form.
        # The reduction's isometry is directed from the reduced lattice to L.
        reduction = indefinite_reduced_presentation(lattice)
        cached = self._groupoid.cached_group(reduction.reduced)
        if cached is not None:
            return _conjugated(reduction.isometry, cached)
        if not reduction.reduced.is_nondegenerate():
            # IndefiniteForms.g:1496–1511 applies ExtendIsometryGroup to
            # the nondegenerate quotient before transporting the generators.
            reduced_generators = self._degenerate_group_generators(reduction.reduced)
        else:
            profile = AttackProfile.from_lattice(reduction.reduced)
            match profile.positive_index:
                case 1:
                    generators = lorentzian_orthogonal_group_generators(profile.signed_view)
                case _:
                    generators = self._generation_step(profile.signed_view)
            reduced_generators = tuple(_untwisted(generator, profile) for generator in generators)
        self._groupoid.record_group(reduction.reduced, reduced_generators)
        return _conjugated(reduction.isometry, reduced_generators)

    def _model(self, lattice: Lattices.ParentMethods) -> OrbitCoverModel:
        r"""The approximate model of ``lattice``, built once per lattice (IndefiniteForms.g:1338, :1346, :1450, :1530)."""
        if lattice not in self._models:
            self._models[lattice] = OrbitCoverModel(lattice)
        return self._models[lattice]

    def _generation_step(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        model = self._model(lattice)
        vector = model.splitting_vector()
        transporters = tuple(
            isometry for candidate in model.covering_vectors(int(vector.q()), int(vector.content())) if (isometry := self.pointed_isometry(vector, candidate)) is not None
        )
        return model.approximate_generators() + self.pointed_automorphism_generators(vector) + transporters

    def isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        # IndefiniteForms.g:1387–1423 splits the radical *before* the
        # nondegenerate genus/discriminant prefilters and memoized recursion.
        if source.is_nondegenerate() != target.is_nondegenerate():
            return None
        if not source.is_nondegenerate():
            if source.module_rank() != target.module_rank():
                return None
            return self._degenerate_isometry(source, target)
        if LatticePrefilter.from_lattice(source) != LatticePrefilter.from_lattice(target):
            return None
        if source.is_even() and source.genus() != target.genus():
            return None
        if not source.is_definite():
            source_reduction = indefinite_reduced_presentation(source)
            target_reduction = indefinite_reduced_presentation(target)
            reduced_witness = self._indefinite_isometry(source_reduction.reduced, target_reduction.reduced)
            if reduced_witness is None:
                return None
            return _transported(target_reduction.isometry, reduced_witness, source_reduction.isometry)
        # IndefiniteForms.g:1239-1273 INDEF_FORM_TestEquivalence_Memoize
        # records every exact verdict, the definite kernel's (:1326)
        # included, so later O(K) requests are transported along it.
        known, cached = self._groupoid.cached_isometry(source, target)
        if known:
            return cached
        homset = source.Isom(target)
        witness = None if homset.is_empty() else homset.an_element()
        self._groupoid.record_isometry(source, target, witness)
        return witness

    def _indefinite_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        known, cached = self._groupoid.cached_isometry(source, target)
        if known:
            return cached
        source_profile = AttackProfile.from_lattice(source)
        target_profile = AttackProfile.from_lattice(target)
        signed_source, signed_target = source_profile.signed_view, target_profile.signed_view
        match source_profile.positive_index:
            case 1:
                witness = lorentzian_isometry(signed_source, signed_target)
            case _:
                source_model = self._model(signed_source)
                target_model = self._model(signed_target)
                vector = source_model.splitting_vector()
                square, content = int(vector.q()), int(vector.content())
                # IndefiniteForms.g:1320-1359 and
                # ApproximateModels.h:GetCoveringOrbitRepresentatives: enumerate
                # the full target covering family at the source's first norm,
                # testing its pointed isometries in order. A separate
                # GetOneOrbitRepresentative trial is not part of this kernel:
                # it can demand an additional exact pointed equivalence and
                # does not replace any covering case.
                covering = target_model.covering_vectors(square, content)
                witness = next(
                    (isometry for candidate in covering if (isometry := self.pointed_isometry(vector, candidate)) is not None),
                    None,
                )
        result = None if witness is None else _untwisted(witness, source_profile)
        self._groupoid.record_isometry(source, target, result)
        return result

    def _lattice_action(self, vector: Lattices.ElementMethods) -> IntegralStructureAction:
        r"""The rational extensions acting at L, with a finite stable containing module.

        source: IndefiniteForms.g:988-1052 INDEF_FORM_GetVectorStructure for q(v) != 0:
        each generator of O(v^perp) is extended by the identity on v
        (``PmatInv * ExpandMatrix(eEndoRed) * Pmat``, :1019-1020).
        N=Zv+v-perp has finite index in L and is stable under every extension.
        If d annihilates L/N, then L is contained in the stable module N/d.
        This is the boundedness hypothesis for the finite residue computation.
        """
        lattice = vector.parent()
        perpendicular = vector.orthogonal_complement().inclusion()
        complement = perpendicular.domain()
        extensions = tuple(generator.extension_across(vector, vector) for generator in self.orthogonal_group_generators(complement))
        rational = lattice.base_change(lattice.base_ring().fraction_field_map())
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
        space = inclusion.codomain()
        split = lattice.subobject_on((vector, *(perpendicular(point) for point in complement.module_generators())))
        denominator = split.inclusion().cokernel().exponent()
        scale = rational.base_ring().one() / rational.base_ring()(denominator)
        containing = Modules(lattice.base_ring()).Mono(split, space)(
            lambda label: space.wrap(rational.scalar_multiple(scale, inclusion(split.inclusion()(split.module_generator(label))).underlying_element()))
        )
        return IntegralStructureAction(rational.O().subgroup(extensions), inclusion).with_stable_containing_bound(containing)

    def isotropic_stabilizer_action(
        self, isotropic: tuple[Lattices.ElementMethods, ...], radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...]
    ) -> IntegralStructureAction:
        r"""The stabilizer problem of the primitive totally isotropic J = span(isotropic), at L.

        source: IndefiniteForms.g:1555-1568 INDEF_FORM_StabilizerVector with :988-1052
        INDEF_FORM_GetVectorStructure (isotropic v), and :1060-1131
        INDEF_FORM_GetRec_IsotropicKplane MapOrthogonalSublatticeGroup with :1710-1723
        INDEF_FORM_Stabilizer_IsotropicKstuff_Kernel (k-planes and flags); the stabilizer
        of L is MatrixIntegral_Stabilizer.
        Aut(N) for N = J-perp is Hom(K, J) x| (O(K) x Gamma) with K = J-perp/J and Gamma the
        allowed action on J (trivial for a pointed vector, +-1 for a line, GL(J) for a plane, the
        flag stabilizer for a flag), given as images of the vectors of J. Its generators are
        extended by the preamble Witt operation, the kernel generators at the common
        denominator of those extensions are adjoined, and the stabilizer of L in the rational group
        they generate is the stabilizer of J (with its marking). Every element preserves N, so the
        orbit span of L is computed inside L/D for D = 2 d^2 e^2 a, where d is that denominator,
        e = |disc K| and a = |coker(L -> J^*)| (issue #18); leaving that module raises.
        """
        lattice = isotropic[0].parent()
        reduction = lattice.primitive_isotropic_subobject(*isotropic).isotropic_reduction()
        quotient = reduction.quotient_lattice()
        perpendicular = reduction.orthogonal_complement().inclusion()
        frame = reduction.coordinate_frame()
        labels = tuple(quotient.module_generating_set())
        lifts = {label: perpendicular(frame(label)) for label in labels}
        basis = (*isotropic, *(lifts[label] for label in labels))

        def lift(element):
            coordinates = element.to_vector()
            return sum((lattice.scalar_multiple(coordinates(label), lifts[label]) for label in labels if coordinates(label)), lattice.zero())

        # Images of the basis (J, then the frame lifts), in that order.
        frame = tuple(lifts[label] for label in labels)
        image_families = []
        for isometry in self.orthogonal_group_generators(quotient):
            image_families.append((*isotropic, *(lift(isometry(quotient.module_generator(label))) for label in labels)))
        for radical in radical_images:
            image_families.append((*radical, *frame))
        for sheared in range(len(labels)):
            for point in isotropic:
                image_families.append((*isotropic, *(lifted + point if position == sheared else lifted for position, lifted in enumerate(frame))))
        extended = [_extend_perpendicular_isometry(reduction, reduction, basis, images) for images in image_families]
        denominator = SageZZ.one()
        for _isometry, entry_denominator in extended:
            denominator = denominator.lcm(entry_denominator)
        # Lorentzian.g:115-249 supplies 1 + T/d for every generator T of the
        # saturated integral homogeneous module. The preamble supplies 1 + T
        # through the complete pointwise-perpendicular kernel. Since T kills
        # N=J^perp and has image in J contained in N, T_i T_j = 0: the scaled
        # generators give exactly the source's (1/d)-parameter module.
        field_map = lattice.base_ring().fraction_field_map()
        rational = lattice.base_change(field_map)
        inverse_denominator = rational.base_ring().one() / rational.base_ring()(int(denominator))
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
        kernel_generators = []
        for integral_generator in reduction.pointwise_perpendicular_kernel().group_generators():
            # T = generator - 1 is read on the integral basis and rationalized;
            # 1 + T/d is the rational map with those scaled images.
            def scaled_image(label, generator=integral_generator):
                point = lattice.module_generator(label)
                direction = inclusion(generator(point) - point).underlying_element()
                return rational.module_generator(label) + rational.scalar_multiple(inverse_denominator, direction)

            kernel_generators.append(rational.Isom(rational)(scaled_image))
        generators = (*(isometry for isometry, _entry_denominator in extended), *kernel_generators)
        space = inclusion.codomain()
        # Unlike the rational coordinate complement, the pairing L -> J^*
        # determines the invariant index a in D=2d^2e^2a intrinsically.
        isotropic_inclusion = reduction.isotropic_embedding()
        # The pairing L -> I^* is i^* composed with the metric map,
        # where i: I -> L. Dualization is represented on the opposite
        # category, so its input is the reversed arrow op(L) -> op(I).
        dualization = Modules(lattice.base_ring()).dualization()
        opposite = dualization.domain()
        opposite_inclusion = opposite.Mor(opposite.object(lattice), opposite.object(isotropic_inclusion.domain()))(isotropic_inclusion)
        pairing = dualization(opposite_inclusion) * lattice.metric_map()
        # Here pairing is a map between finite free ZZ-modules L -> J*. Its
        # cokernel is finite, and its order is the product of the nonzero Smith
        # entries of its coordinate matrix. Compute that exact index from these
        # presentations: the generic provider's index() currently delegates to
        # cardinality() on a represented torsion module whose construction does
        # not yet supply a decided finite cardinality.
        pairing_domain_labels = tuple(pairing.domain().module_generating_set())
        pairing_codomain_labels = tuple(pairing.codomain().module_generating_set())
        pairing_rows = tuple(
            tuple(pairing(pairing.domain().module_generator(label)).to_vector()(target_label) for target_label in pairing_codomain_labels)
            for label in pairing_domain_labels
        )
        pairing_matrix = matrix(SageZZ, pairing_rows)
        pairing_rank = int(pairing_matrix.rank())
        if pairing_rank != len(pairing_codomain_labels):
            raise ArithmeticError("the isotropic pairing L -> J* must have finite cokernel")
        smith = pairing_matrix.smith_form()[0]
        pairing_index = SageZZ.one()
        for position in range(pairing_rank):
            pairing_index *= abs(smith[position, position])
        bound = 2 * int(denominator) ** 2 * abs(int(quotient.discriminant())) ** 2 * int(pairing_index)
        scale = rational.base_ring().one() / rational.base_ring()(bound)
        containing = Modules(lattice.base_ring()).Mono(lattice, space)(
            lambda label: space.wrap(rational.scalar_multiple(scale, inclusion(lattice.module_generator(label)).underlying_element()))
        )
        action = IntegralStructureAction(rational.O().subgroup(tuple(generators)), inclusion)
        span = OrbitSpanConstructions(lattice.base_ring())(action.ambient_action_functor(), inclusion, containing=containing)
        return action.with_stable_containing_bound(span.inclusion())

    def line_stabilizer_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, Zv)) for primitive isotropic v: Gamma = {+-1} on the line."""
        lattice = vector.parent()
        action = self.isotropic_stabilizer_action((vector,), ((-vector,),))
        return tuple(lattice.integral_isometry(generator) for generator in action.stabilizer().group_generators())

    def levi_image(self, vector: Lattices.ElementMethods) -> Parent:
        r"""H, the image of Aut((L, Zv)) in O(K), K = v-perp/v, generated by the images of its computed generators."""
        reduction = vector.parent().primitive_isotropic_subobject(vector).isotropic_reduction()
        levi = reduction.levi_action()
        parabolic = reduction.parabolic_subgroup()
        return reduction.quotient_lattice().O().subgroup(tuple(levi(parabolic(generator)) for generator in self.line_stabilizer_generators(vector)))

    def pointed_automorphism_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, v)): the stabilizer of L in the rational extensions of Aut(v-perp) (IndefiniteForms.g:1533-1545)."""
        lattice = vector.parent()
        action = self.isotropic_stabilizer_action((vector,), ()) if vector.q() == lattice.base_ring().zero() else self._lattice_action(vector)
        return tuple(lattice.integral_isometry(generator) for generator in action.stabilizer().group_generators())

    def pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""An element of Isom((L, v), (M, w)), or None."""
        if VectorPrefilter.from_vector(source) != VectorPrefilter.from_vector(target):
            return None
        source_lattice, target_lattice = source.parent(), target.parent()
        if source_lattice is target_lattice and source == target:
            return source_lattice.O().identity()
        if source_lattice is not target_lattice and LatticePrefilter.from_lattice(source_lattice) != LatticePrefilter.from_lattice(target_lattice):
            return None
        if source.q() == source_lattice.base_ring().zero():
            if source_lattice.is_even() and source.div() == source_lattice.base_ring().one():
                return self._split_isotropic_isometry(source, target)
            return self._isotropic_pointed_isometry(source, target)
        source_complement = source.orthogonal_complement().inclusion().domain()
        target_complement = target.orthogonal_complement().inclusion().domain()
        reduced = self.isometry(source_complement, target_complement)
        if reduced is None:
            return None
        rational_transport = reduced.extension_across(source, target)
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)
        fraction_map = source_lattice.base_ring().fraction_field_map()
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()
        correction = self._lattice_action(target).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction.inverse() * rational_transport, source=source_lattice)

    def _reduced_isotropic_extension(self, source_reduction, target_reduction, source, target, reduced) -> LatticeIsometryMethods:
        r"""The rational extension L_Q -> M_Q of the integral isometry N -> N' that ``reduced`` determines.

        source: IndefiniteForms.g:1623-1641 (vectors) and :1673-1686 (k-planes and
        flags). The source's ``test`` is an integral isometry of the perpendicular
        lattices N=J^perp -> N'=J'^perp; here it sends the basis ``source`` of J to
        ``target`` and the frame lift of each k in K=N/J to the frame lift of
        ``reduced(k)``. Both families are Z-bases of N and N', so this is an
        integral isometry onto N' and carries the marking as given (the sign
        normalization of :1639-1641 is then the identity). Its extension is
        ``LORENTZ_ExtendOrthogonalIsotropicIsomorphism``.
        """
        target_lattice = target[0].parent()
        source_quotient = source_reduction.quotient_lattice()
        target_quotient = target_reduction.quotient_lattice()
        source_labels = tuple(source_quotient.module_generating_set())
        target_labels = tuple(target_quotient.module_generating_set())
        source_frame = source_reduction.coordinate_frame()
        target_frame = target_reduction.coordinate_frame()
        source_perpendicular = source_reduction.orthogonal_complement().inclusion()
        target_perpendicular = target_reduction.orthogonal_complement().inclusion()
        source_lifts = tuple(source_perpendicular(source_frame(label)) for label in source_labels)
        target_lifts = {label: target_perpendicular(target_frame(label)) for label in target_labels}

        def target_lift(element):
            coordinates = element.to_vector()
            return sum((target_lattice.scalar_multiple(coordinates(label), target_lifts[label]) for label in target_labels if coordinates(label)), target_lattice.zero())

        images = (*target, *(target_lift(reduced(source_quotient.module_generator(label))) for label in source_labels))
        rational_transport, _denominator = _extend_perpendicular_isometry(source_reduction, target_reduction, (*source, *source_lifts), images)
        return rational_transport

    def _isotropic_correction(
        self,
        rational_transport: LatticeIsometryMethods,
        source_lattice,
        target: tuple[Lattices.ElementMethods, ...],
        radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...],
    ) -> LatticeIsometryMethods | None:
        r"""source: IndefiniteForms.g:1569-1604 Kernel_Equivalence_Qmat.

        An integral ``EquivRat`` is returned as it is (:1578-1580). Otherwise
        the source searches its extended stabilizer group of J (:1581-1590,
        MatrixIntegral_Equivalence_Bis) for g with ``EquivRat * g`` integral.
        Here the same search is the integral-structure transporter of the
        conjugate group at J' (:meth:`isotropic_stabilizer_action`), with the
        allowed action ``radical_images`` on J'. The transport must carry N onto
        N' integrally (:meth:`_reduced_isotropic_extension`); every group element
        preserves N', so a correction exists exactly when an isometry with the
        marking exists.
        """
        target_lattice = target[0].parent()
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)
        fraction_map = source_lattice.base_ring().fraction_field_map()
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()
        correction = self.isotropic_stabilizer_action(tuple(target), radical_images).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction.inverse() * rational_transport, source=source_lattice)

    def isotropic_subspace_isometry(
        self, source: tuple[Lattices.ElementMethods, ...], target: tuple[Lattices.ElementMethods, ...], radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...]
    ) -> LatticeIsometryMethods | None:
        r"""An isometry L -> M carrying J = span(source) onto J' = span(target), compatible with the marking.

        source: IndefiniteForms.g:1660-1694 INDEF_FORM_Equivalence_IsotropicKstuff_Kernel.
        An isometry psi of the reductions K -> K' (``f_equiv``, :1673) gives an
        isometry of the perpendicular lattices, its rational extension
        (:meth:`_reduced_isotropic_extension`), and Kernel_Equivalence_Qmat
        (:meth:`_isotropic_correction`), whose group's allowed action
        ``radical_images`` on J' (images of ``target``) absorbs the choice of
        basis of J'.
        """
        source_lattice = source[0].parent()
        target_lattice = target[0].parent()
        source_reduction = source_lattice.primitive_isotropic_subobject(*source).isotropic_reduction()
        target_reduction = target_lattice.primitive_isotropic_subobject(*target).isotropic_reduction()
        reduced = self.isometry(source_reduction.quotient_lattice(), target_reduction.quotient_lattice())
        if reduced is None:
            return None
        rational_transport = self._reduced_isotropic_extension(source_reduction, target_reduction, source, target, reduced)
        return self._isotropic_correction(rational_transport, source_lattice, target, radical_images)

    def _isotropic_pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""Isom of primitive isotropic marked vectors.

        source: IndefiniteForms.g:1605-1651 INDEF_FORM_EquivalenceVector, isotropic
        branch. ``test`` (:1623) is a recursive isometry of the reductions; the
        rational transport extends the integral isometry of perpendicular
        lattices it determines and is integral already in the Ansatz case
        (:1578). Otherwise the preamble's marked-line lift torsor supplies the
        complete integral fibre over this reduction isometry (port-design
        "Isometry-extension fibres"). An empty fibre for this one reduction map
        does not make the pointed Isom empty: Kernel_Equivalence_Qmat then
        ranges over the full stabilizer of the marked vector, which acts on K
        through all of O(K).
        """
        source_lattice, target_lattice = source.parent(), target.parent()
        source_reduction = source_lattice.primitive_isotropic_subobject(source).isotropic_reduction()
        target_reduction = target_lattice.primitive_isotropic_subobject(target).isotropic_reduction()
        reduced = self.isometry(source_reduction.quotient_lattice(), target_reduction.quotient_lattice())
        if reduced is None:
            return None
        rational_transport = self._reduced_isotropic_extension(source_reduction, target_reduction, (source,), (target,), reduced)
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)

        torsor = source_reduction.rational_lifts(target_reduction, reduced)
        if torsor.integral_members().is_empty() is False:
            integral = torsor.integral_base_point()
            # The torsor fixes the reductions' chosen line generators; when one
            # is the negative of the supplied marking, the global sign of
            # IndefiniteForms.g:1639-1641 is applied against the actual vectors.
            fraction_map = source_lattice.base_ring().fraction_field_map()
            source_vector = ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()(source).underlying_element()
            target_vector = ModulesOverIntegralDomains(target_lattice.base_ring())(target_lattice).generic_fibre_map()(target).underlying_element()
            mapped = integral(source_vector)
            if mapped == -target_vector:
                target_rational = target_lattice.base_change(fraction_map)
                integral = target_rational.Isom(target_rational)(lambda label: -target_rational.module_generator(label)) * integral
            elif mapped != target_vector:
                raise ArithmeticError("the integral reduction lift does not carry the marked isotropic line onto the requested target")
            return target_lattice.integral_isometry(integral, source=source_lattice)
        return self._isotropic_correction(rational_transport, source_lattice, (target,), ())

    def _split_isotropic_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""For e of divisibility one in an even lattice, L = U + e^perp/e; an isometry is id_U plus one of the reductions."""
        source_splitting = source.parent().primitive_isotropic_subobject(source).isotropic_reduction().integral_hyperbolic_splitting()
        target_splitting = target.parent().primitive_isotropic_subobject(target).isotropic_reduction().integral_hyperbolic_splitting()
        hyperbolic_plane, source_complement = source_splitting.codomain().biproduct_factors()
        _hyperbolic_plane, target_complement = target_splitting.codomain().biproduct_factors()
        reduced = self.isometry(source_complement, target_complement)
        if reduced is None:
            return None
        split_map = hyperbolic_plane.identity_morphism().biproduct_map(reduced, source=source_splitting.codomain(), target=target_splitting.codomain())
        # The biproduct operation returns a module map.  Its summands are
        # isometries, so retain that structure before composing the splittings.
        split = source_splitting.codomain().Isom(target_splitting.codomain())(
            lambda label: split_map(source_splitting.codomain().module_generator(label))
        )
        return target_splitting.inverse() * split * source_splitting

    def primitive_vector_classes(self, lattice: Lattices.ParentMethods, square: int) -> tuple[VectorClass, ...]:
        r"""One class per O(L)-orbit of primitive vectors of the square, for L of positive index at least 2."""
        profile = AttackProfile.from_lattice(lattice)
        if profile.positive_index < 2:
            raise ValueError(
                f"{lattice} has positive index {profile.positive_index}: vector classes here need positive index at least 2; the Lorentzian case is marked cells (#15)"
            )
        # Real positive and negative index at least two does not imply
        # rational isotropy: diag(2, 2, -6, -6), for example, is anisotropic
        # by descent modulo 3.  Its zero-square fibre has no primitive
        # vectors, so decide that empty case before requesting the two
        # rational hyperbolic pairs of the Eichler envelope.
        if square == 0 and not lattice.is_isotropic():
            return ()
        buckets: dict[VectorPrefilter, list[tuple[Lattices.ElementMethods, list[LatticeIsometryMethods]]]] = {}
        for candidate in self._model(profile.signed_view).primitive_covering_vectors(profile.sign * int(square)):
            bucket = buckets.setdefault(VectorPrefilter.from_vector(candidate), [])
            match next(((witnesses, isometry) for kept, witnesses in bucket if (isometry := self.pointed_isometry(kept, candidate)) is not None), None):
                case None:
                    bucket.append((candidate, []))
                case (witnesses, isometry):
                    witnesses.append(isometry)
        return tuple(
            VectorClass(
                lattice.linear_combination({label: kept.to_vector()(label) for label in lattice.module_generating_set()}),
                tuple(_untwisted(isometry, profile) for isometry in witnesses),
            )
            for bucket in buckets.values()
            for kept, witnesses in bucket
        )


__all__ = ["EichlerRecursion", "VectorClass"]
