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
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from sage_indefinite_port.indefinite.groupoid_cache import IsometryGroupoidCache
from sage_indefinite_port.indefinite.indefinite_reduction import (
    indefinite_reduced_presentation,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    _coordinate_row,
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

    The specialized source, Indefinite.jl@374a5ebb Lorentzian.g:115–249,
    selects a rational extension of an isometry between coisotropic subspaces.
    The preamble owns that extension and its affine norm-correction equations.
    This adapter only expresses the supplied images of a basis of I^perp as
    the partial isometry on the preamble's actual perpendicular subobjects,
    then reads the extension's denominator ideal.
    """
    source = source_reduction.isotropic_embedding().codomain()
    target = target_reduction.isotropic_embedding().codomain()
    fraction_map = source.base_ring().fraction_field_map()
    rational_source = source.base_change(fraction_map)
    rational_target = target.base_change(fraction_map)
    source_inclusion = source_reduction.orthogonal_complement().inclusion().base_change(fraction_map)
    target_inclusion = target_reduction.orthogonal_complement().inclusion().base_change(fraction_map)
    source_perpendicular = source_inclusion.domain()
    target_perpendicular = target_inclusion.domain()
    if len(source_basis) != int(source_perpendicular.module_rank()) or len(images) != len(source_basis):
        raise ValueError("the perpendicular partial isometry needs images of a full chosen basis")

    basis = matrix(SageQQ, [_coordinate_row(source, vector, SageQQ) for vector in source_basis])
    if basis.rank() != len(source_basis):
        raise ValueError("the selected perpendicular generators are linearly dependent")
    target_labels = tuple(target.module_generating_set())
    target_field = rational_target.base_ring()

    def rationalize(point):
        coordinates = point.to_vector()
        return rational_target.linear_combination({label: target_field(int(coordinates(label))) for label in target_labels if coordinates(label)})

    partial_images = tuple(target_inclusion.lift(rationalize(point)) for point in images)

    def partial_image(label):
        source_point = source_inclusion(source_perpendicular.module_generator(label))
        coefficients = basis.solve_left(_coordinate_row(rational_source, source_point, SageQQ))
        return sum(
            (
                target_perpendicular.scalar_multiple(
                    target_field(int(coefficient.numerator())) / target_field(int(coefficient.denominator())),
                    image,
                )
                for coefficient, image in zip(coefficients, partial_images, strict=True)
                if coefficient
            ),
            target_perpendicular.zero(),
        )

    partial = source_perpendicular.Isom(target_perpendicular)(partial_image)
    extension = partial.witt_extension(source_inclusion, target_inclusion)
    (denominator,) = tuple(extension.denominator_ideal(source, target).ideal_generators())
    return extension, int(denominator)


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
        complement_gram = gram.matrix_from_rows_and_columns(range(nonsingular_rank), range(nonsingular_rank))
        complement = lattice.lattice_category()([[int(entry) for entry in row] for row in complement_gram.rows()])
        matrices = []
        if nonsingular_rank:
            for generator in self.orthogonal_group_generators(complement):
                top = tuple(tuple(int(value) for value in generator(point).to_vector()) for point in complement.module_generators())
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
        return tuple(lattice.O()(tuple(lattice(list(row)) for row in block.rows())) for block in matrices)

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
            return None if cached is None else target_reduction.isometry * cached * source_reduction.isometry.inverse()
        left_complement_gram = left_gram.matrix_from_rows_and_columns(range(rank), range(rank))
        right_complement_gram = right_gram.matrix_from_rows_and_columns(range(rank), range(rank))
        left_complement = left.lattice_category()([[int(entry) for entry in row] for row in left_complement_gram.rows()])
        right_complement = right.lattice_category()([[int(entry) for entry in row] for row in right_complement_gram.rows()])
        complement_witness = self.isometry(left_complement, right_complement) if rank else None
        if rank and complement_witness is None:
            self._groupoid.record_isometry(left, right, None)
            return None
        extension = identity_matrix(SageZZ, size)
        if complement_witness is not None:
            for row, point in enumerate(left_complement.module_generators()):
                coordinates = tuple(int(value) for value in complement_witness(point).to_vector())
                for column, coefficient in enumerate(coordinates):
                    extension[row, column] = coefficient
        reduced_witness = left.Isom(right)(tuple(right(list(row)) for row in extension.rows()))
        self._groupoid.record_isometry(left, right, reduced_witness)
        return target_reduction.isometry * reduced_witness * source_reduction.isometry.inverse()

    def orthogonal_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        if lattice.is_definite():
            return tuple(lattice.O().group_generators())
        # IndefiniteForms.g:1474-1494: reduce before the memoized group
        # computation, then conjugate its generators back to the input form.
        # The reduction's isometry is directed from the reduced lattice to L.
        reduction = indefinite_reduced_presentation(lattice)
        cached = self._groupoid.cached_group(reduction.reduced)
        if cached is not None:
            return tuple(reduction.isometry * generator * reduction.isometry.inverse() for generator in cached)
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
        return tuple(reduction.isometry * generator * reduction.isometry.inverse() for generator in reduced_generators)

    def _generation_step(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        model = OrbitCoverModel(lattice)
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
            return target_reduction.isometry * reduced_witness * source_reduction.isometry.inverse()
        if source.is_definite():
            homset = source.Isom(target)
            return None if homset.is_empty() else homset.an_element()
        return self._indefinite_isometry(source, target)

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
                source_model = OrbitCoverModel(signed_source)
                target_model = OrbitCoverModel(signed_target)
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

        IndefiniteForms.g:1533-1545 (vectors) and :1060-1130, :1711-1725 (k-planes and flags):
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
        kernel_generators = []
        for integral_generator in reduction.pointwise_perpendicular_kernel().group_generators():
            rational_generator = integral_generator.base_change(field_map)

            def scaled_image(label, generator=rational_generator):
                point = rational.module_generator(label)
                return point + rational.scalar_multiple(inverse_denominator, generator(point) - point)

            kernel_generators.append(rational.Isom(rational)(scaled_image))
        generators = (*(isometry for isometry, _entry_denominator in extended), *kernel_generators)
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
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
        bound = 2 * int(denominator) ** 2 * abs(int(quotient.discriminant())) ** 2 * int(pairing.index())
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

    def isotropic_subspace_isometry(
        self, source: tuple[Lattices.ElementMethods, ...], target: tuple[Lattices.ElementMethods, ...], radical_images: tuple[tuple[Lattices.ElementMethods, ...], ...]
    ) -> LatticeIsometryMethods | None:
        r"""An isometry L -> M carrying J = span(source) onto J' = span(target), compatible with the marking.

        INDEF_FORM_Equivalence_IsotropicKstuff_Kernel and INDEF_FORM_EquivalenceVector
        (IndefiniteForms.g:1607-1690). An isometry psi of the reductions K -> K' gives an isometry
        of the perpendicular lattices sending the basis of J to the basis of J' and the frame lift
        of k to the frame lift of psi(k); one extension of it
        (LORENTZ_ExtendOrthogonalIsotropicIsomorphism) is a rational isometry L -> M. When it is
        not integral, Kernel_Equivalence_Qmat looks for an element g of the extended stabilizer
        group of J' with g composed with it integral: the transporter of
        :meth:`isotropic_stabilizer_action` at J', with the same allowed action ``radical_images``
        on J' (images of ``target``), which absorbs the choice of basis of J'.
        """
        source_lattice, target_lattice = source[0].parent(), target[0].parent()
        source_reduction = source_lattice.primitive_isotropic_subobject(*source).isotropic_reduction()
        target_reduction = target_lattice.primitive_isotropic_subobject(*target).isotropic_reduction()
        source_quotient = source_reduction.quotient_lattice()
        target_quotient = target_reduction.quotient_lattice()
        reduced = self.isometry(source_quotient, target_quotient)
        if reduced is None:
            return None
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
        if rational_transport.is_integral_on(source_lattice, target_lattice):
            return target_lattice.integral_isometry(rational_transport, source=source_lattice)
        fraction_map = source_lattice.base_ring().fraction_field_map()
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()
        correction = self.isotropic_stabilizer_action(tuple(target), radical_images).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction.inverse() * rational_transport, source=source_lattice)

    def _isotropic_pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""Isom of primitive isotropic marked vectors, with the full target correction.

        IndefiniteForms.g:1605-1650 first selects a reduction isometry. The
        preamble's marked-line lift torsor supplies its entire integral fibre;
        if that fibre is empty, the source's integral-equivalence correction
        still ranges over the target's full reduced orthogonal-group action.
        """
        source_lattice, target_lattice = source.parent(), target.parent()
        source_reduction = source_lattice.primitive_isotropic_subobject(source).isotropic_reduction()
        target_reduction = target_lattice.primitive_isotropic_subobject(target).isotropic_reduction()
        reduced = self.isometry(source_reduction.quotient_lattice(), target_reduction.quotient_lattice())
        if reduced is None:
            return None
        torsor = source_reduction.rational_lifts(target_reduction, reduced)
        integral_locus = torsor.integral_members()
        empty = integral_locus.is_empty()

        fraction_map = source_lattice.base_ring().fraction_field_map()
        source_inclusion = ModulesOverIntegralDomains(source_lattice.base_ring())(source_lattice).generic_fibre_map()
        target_inclusion = ModulesOverIntegralDomains(target_lattice.base_ring())(target_lattice).generic_fibre_map()
        source_vector = source_inclusion(source).underlying_element()
        target_vector = target_inclusion(target).underlying_element()
        rational_transport = torsor.base_point()
        # The reduction's chosen primitive generator may differ by a global
        # sign from the supplied pointed marking (IndefiniteForms.g:1637-1642).
        mapped = rational_transport(source_vector)
        negate = mapped == -target_vector
        if not negate and mapped != target_vector:
            raise ArithmeticError("the rational reduction lift does not carry the marked isotropic line onto the requested target")
        if negate:
            target_rational = target_lattice.base_change(fraction_map)
            negative_identity = target_rational.Isom(target_rational)(lambda label: -target_rational.module_generator(label))
            rational_transport = negative_identity * rational_transport

        if empty is False:
            integral = torsor.integral_base_point()
            if negate:
                integral = negative_identity * integral
            return target_lattice.integral_isometry(integral, source=source_lattice)

        # Even when the lift-locus emptiness is undecided, the full target
        # stabilizer's integral-structure transporter decides whether a
        # correction of this rational transport carries L onto M.  Emptiness
        # for one fixed reduction map does not decide pointed isometry.
        restriction = Modules(fraction_map.codomain()).restriction_of_scalars(fraction_map)
        image = restriction(rational_transport) * source_inclusion
        correction = self.isotropic_stabilizer_action((target,), ()).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction.inverse() * rational_transport, source=source_lattice)

    def _split_isotropic_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""For e of divisibility one in an even lattice, L = U + e^perp/e; an isometry is id_U plus one of the reductions."""
        source_splitting = source.parent().primitive_isotropic_subobject(source).isotropic_reduction().integral_hyperbolic_splitting()
        target_splitting = target.parent().primitive_isotropic_subobject(target).isotropic_reduction().integral_hyperbolic_splitting()
        hyperbolic_plane, source_complement = source_splitting.codomain().biproduct_factors()
        _hyperbolic_plane, target_complement = target_splitting.codomain().biproduct_factors()
        reduced = self.isometry(source_complement, target_complement)
        if reduced is None:
            return None
        split = hyperbolic_plane.identity_morphism().biproduct_map(reduced, source=source_splitting.codomain(), target=target_splitting.codomain())
        return target_splitting.inverse() * split * source_splitting

    def primitive_vector_classes(self, lattice: Lattices.ParentMethods, square: int) -> tuple[VectorClass, ...]:
        r"""One class per O(L)-orbit of primitive vectors of the square, for L of positive index at least 2."""
        profile = AttackProfile.from_lattice(lattice)
        if profile.positive_index < 2:
            raise ValueError(
                f"{lattice} has positive index {profile.positive_index}: vector classes here need positive index at least 2; the Lorentzian case is marked cells (#15)"
            )
        buckets: dict[VectorPrefilter, list[tuple[Lattices.ElementMethods, list[LatticeIsometryMethods]]]] = {}
        for candidate in OrbitCoverModel(profile.signed_view).primitive_covering_vectors(profile.sign * int(square)):
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
