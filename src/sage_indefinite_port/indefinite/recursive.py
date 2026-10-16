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
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.orbit_spans import OrbitSpanConstructions
from dzack_research.preamble.categories.modules.pure.modules import Modules, ModulesOverIntegralDomains
from dzack_research.preamble.categories.rational_integral_stabilizers import IntegralStructureAction
from sage.matrix.constructor import matrix
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from sage_indefinite_port.indefinite.lorentzian_cells import _coordinate_row, _gram_matrix, lorentzian_isometry, lorentzian_orthogonal_group_generators
from sage_indefinite_port.invariants import AttackProfile, LatticePrefilter, VectorPrefilter

_UNTWIST = Lattices(ZZ).twist_functor(-1)


def _untwisted(isometry: LatticeIsometryMethods, profile: AttackProfile) -> LatticeIsometryMethods:
    r"""An isometry of the signed view, as an isometry of the lattice itself."""
    return isometry if profile.sign == 1 else _UNTWIST(isometry)


def _rational_isometry(source: Lattices.ParentMethods, target: Lattices.ParentMethods, transformation: Matrix_rational_dense) -> LatticeIsometryMethods:
    r"""Reconstruct the rational isometry L -> M whose row action, from coordinates of L to coordinates of M, is ``transformation``."""
    fraction_map = source.base_ring().fraction_field_map()
    rational_source = source.base_change(fraction_map)
    rational_target = target.base_change(fraction_map)
    source_labels = tuple(source.module_generating_set())
    target_labels = tuple(target.module_generating_set())
    field = rational_target.base_ring()

    def image(label):
        row = transformation.row(source_labels.index(label))
        return rational_target.linear_combination(
            {target: field(int(entry.numerator())) / field(int(entry.denominator())) for target, entry in zip(target_labels, row) if entry}
        )

    return rational_source.Isom(rational_target)(image)


type _ExtensionSystem = tuple[
    Lattices.ParentMethods, Matrix_rational_dense, Matrix_rational_dense, Matrix_rational_dense, Matrix_rational_dense, Matrix_rational_dense, Matrix_rational_dense
]
"""The target M, the image rows, the image rows of J, the candidate images of W, the pairing A and the defect B."""


class _IsotropicExtension:
    r"""Extensions to L of an isometry of N = J-perp, for J a primitive totally isotropic sublattice.

    Lorentzian.g:115-249, LORENTZ_ExtendOrthogonalIsotropicIsomorphism (and its codimension-one
    case :4-50). With a rational complement W of N, the image of each w in W is fixed modulo J_Q by
    its pairings with the image of N; writing the images as P + X J, the remaining conditions
    b(w'_i, w'_l) = b(w_i, w_l) read X A + A^T X^T = B with A = (b(J, P)), linear in X because J
    is totally isotropic. ``one`` selects one rational extension; the complete
    integral homogeneous kernel is supplied by the preamble's perpendicular kernel.
    """

    def __init__(self, lattice: Lattices.ParentMethods, isotropic: tuple[Lattices.ElementMethods, ...], perpendicular: tuple[Lattices.ElementMethods, ...]) -> None:
        self.lattice = lattice
        self.gram = _gram_matrix(lattice).change_ring(SageQQ)
        self.isotropic_rows = matrix(SageQQ, [_coordinate_row(lattice, point, SageQQ) for point in isotropic])
        self.perpendicular = perpendicular
        perpendicular_rows = matrix(SageQQ, [_coordinate_row(lattice, point, SageQQ) for point in perpendicular])
        if perpendicular_rows.rank() != perpendicular_rows.nrows():
            raise ValueError("the perpendicular sublattice must be given by a basis")
        if (perpendicular_rows * self.gram * self.isotropic_rows.transpose()) != 0 or (self.isotropic_rows * self.gram * self.isotropic_rows.transpose()) != 0:
            raise ValueError("J must be totally isotropic and perpendicular to N")
        self.perpendicular_rows = perpendicular_rows
        self.complement_rows = _rational_complement(perpendicular_rows)
        self.source_rows = perpendicular_rows.stack(self.complement_rows)

    def _solve(self, images: tuple[Lattices.ElementMethods, ...]) -> _ExtensionSystem:
        r"""The candidate images of W, the pairing A[m, l] = b(j'_m, P_l) and the defect B, in the target M of ``images``."""
        target = images[0].parent()
        target_gram = _gram_matrix(target).change_ring(SageQQ)
        image_rows = matrix(SageQQ, [_coordinate_row(target, image, SageQQ) for image in images])
        system = image_rows * target_gram
        candidates = matrix(SageQQ, [system.solve_right(self.perpendicular_rows * self.gram * row) for row in self.complement_rows])
        isotropic_images = image_rows.matrix_from_rows(range(self.isotropic_rows.nrows()))
        pairing = isotropic_images * target_gram * candidates.transpose()
        defect = self.complement_rows * self.gram * self.complement_rows.transpose() - candidates * target_gram * candidates.transpose()
        return target, target_gram, image_rows, isotropic_images, candidates, pairing, defect

    def one(self, images: tuple[Lattices.ElementMethods, ...]) -> tuple[LatticeIsometryMethods, int]:
        r"""One extension of the isometry of N sending the basis of N, in order, to ``images``, and the denominator of its matrix."""
        target, target_gram, image_rows, isotropic_images, candidates, pairing, defect = self._solve(images)
        # With A[m, l] = b(j'_m, P_l) the conditions read X A + (X A)^T = B; X = (B/2) A^-1 solves them.
        solution = (defect / 2) * pairing.inverse()
        target_rows = image_rows.stack(candidates + solution * isotropic_images)
        if target_rows * target_gram * target_rows.transpose() != self.source_rows * self.gram * self.source_rows.transpose():
            raise ArithmeticError("the isotropic extension does not preserve the form")
        transformation = self.source_rows.inverse() * target_rows
        return _rational_isometry(self.lattice, target, transformation), int(transformation.denominator())


def _rational_complement(rows: Matrix_rational_dense) -> Matrix_rational_dense:
    r"""Fundamental.g:509-515, SubspaceCompletionRational: Euclidean nullspace."""
    return rows.right_kernel_matrix()


@dataclass(frozen=True)
class VectorClass:
    r"""An O(L)-orbit of primitive vectors: a representative, and the isometries of L found carrying it to the other covering vectors of its orbit."""

    representative: Lattices.ElementMethods
    witnesses: tuple[LatticeIsometryMethods, ...]


class EichlerRecursion:
    r"""The recursion; every general step is a preamble operation."""

    def orthogonal_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        if lattice.is_definite():
            return tuple(lattice.O().group_generators())
        profile = AttackProfile.from_lattice(lattice)
        match profile.positive_index:
            case 1:
                generators = lorentzian_orthogonal_group_generators(profile.signed_view)
            case _:
                generators = self._generation_step(profile.signed_view)
        return tuple(_untwisted(generator, profile) for generator in generators)

    def _generation_step(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
        model = OrbitCoverModel(lattice)
        vector = model.splitting_vector()
        transporters = tuple(
            isometry for candidate in model.covering_vectors(int(vector.q()), int(vector.content())) if (isometry := self.pointed_isometry(vector, candidate)) is not None
        )
        return model.approximate_generators() + self.pointed_automorphism_generators(vector) + transporters

    def isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        if LatticePrefilter.from_lattice(source) != LatticePrefilter.from_lattice(target):
            return None
        if source.is_even() and source.genus() != target.genus():
            return None
        if source.is_definite():
            homset = source.Isom(target)
            return None if homset.is_empty() else homset.an_element()
        return self._indefinite_isometry(source, target)

    def _indefinite_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
        source_profile = AttackProfile.from_lattice(source)
        target_profile = AttackProfile.from_lattice(target)
        signed_source, signed_target = source_profile.signed_view, target_profile.signed_view
        match source_profile.positive_index:
            case 1:
                witness = lorentzian_isometry(signed_source, signed_target)
            case _:
                vector = OrbitCoverModel(signed_source).splitting_vector()
                covering = OrbitCoverModel(signed_target).covering_vectors(int(vector.q()), int(vector.content()))
                witness = next((isometry for candidate in covering if (isometry := self.pointed_isometry(vector, candidate)) is not None), None)
        return None if witness is None else _untwisted(witness, source_profile)

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
        extended to L (:class:`_IsotropicExtension`), the kernel generators at the common
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
        extension = _IsotropicExtension(lattice, tuple(isotropic), basis)

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
        extended = [extension.one(images) for images in image_families]
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

            def scaled_image(label):
                point = rational.module_generator(label)
                return point + rational.scalar_multiple(inverse_denominator, rational_generator(point) - point)

            kernel_generators.append(rational.Isom(rational)(scaled_image))
        generators = (*(isometry for isometry, _entry_denominator in extended), *kernel_generators)
        inclusion = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
        space = inclusion.codomain()
        # Unlike the rational coordinate complement, the pairing L -> J^*
        # determines the invariant index a in D=2d^2e^2a intrinsically.
        isotropic_lattice = reduction.isotropic_sublattice()
        isotropic_inclusion = reduction.isotropic_embedding()
        dual = isotropic_lattice.dual_module()
        isotropic_labels = tuple(isotropic_lattice.module_generating_set())
        pairing = lattice.module_category().Mor(lattice, dual)(
            lambda label: dual.linear_combination(
                {
                    isotropic_label: lattice.b(lattice.module_generator(label), isotropic_inclusion(isotropic_lattice.module_generator(isotropic_label)))
                    for isotropic_label in isotropic_labels
                }
            )
        )
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

        extension = _IsotropicExtension(source_lattice, tuple(source), (*source, *source_lifts))
        images = (*target, *(target_lift(reduced(source_quotient.module_generator(label))) for label in source_labels))
        rational_transport, _denominator = extension.one(images)
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
        if empty is not True and empty is not False:
            raise ValueError("the integral marked-line lift locus has undecided emptiness")

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

        # Emptiness for this one reduction map does not decide Isom((L,v),(M,w)).
        # Search its other target reduction actions using the complete rational
        # stabilizer and the integral-structure transporter, as in the source.
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
