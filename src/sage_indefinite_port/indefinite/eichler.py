r"""Dutour Sikirić's approximate model of the Eichler recursion (map §9), over preamble operations.

For a lattice L containing U + U, the preamble's 2U model of L supplies the Eichler
subgroup A of O(L) (SL2(Z) x SL2(Z) on U + U and the Eichler transvections) and, for each
discriminant class covering a square, one primitive vector with those invariants. A has
finitely many orbits on every set of primitive vectors with prescribed (q, div, class), so
those vectors cover every O(L)-orbit of primitive vectors of the square.

A lattice that does not split U + U is reached through Dutour Sikirić's envelope
(GetEichlerHyperplaneBasis): two hyperbolic pairs (v_1, w_1), (v_2, w_2) of L with orthogonal
complement K give a rational isometry L tensor Q -> E tensor Q onto E = U + U + K, and an
integral similarity sigma = N phi: L -> E. The approximate group of L is the subgroup of
A_E preserving phi(L), pulled back to L; its cosets in A_E move covering vectors of E onto
vectors of sigma(L).
"""

from __future__ import annotations

import weakref
from math import gcd
from random import choice, randrange, sample

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.categories.eichler_criterion import TwoUEichlerModel
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.framed.formed.discriminant_modules import (
    DiscriminantQuadraticModules,
)
from dzack_research.preamble.categories.rational_integral_stabilizers import (
    IntegralStructureAction,
)
from dzack_research.preamble.categories.rings.ring_foundation import (
    _owned_engine_element,
)
from dzack_research.preamble.categories.sets.finite_ordered_sets import (
    finite_ordered_set,
)
from dzack_research.preamble.categories.sets.indexed_families import (
    IndexedFamily,
    finite_indexed_family,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.graphs.graph import Graph
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector as vector_space_element
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.element import RingElement
from sage.structure.sage_object import SageObject

from sage_indefinite_port.indefinite.exact_isotropic import find_isotropic_exact
from sage_indefinite_port.indefinite.indefinite_reduction import (
    _compute_indefinite_reduction,
    _gram_l1,
    _indefinite_lll,
    _random_signed_permutation,
)


def _isometry_from_row_matrix(lattice: Lattices.ParentMethods, rows) -> LatticeIsometryMethods:
    """Build an exact owned isometry from source-defined row coordinates.

    Share the recursive module's checked row-matrix boundary so the pinned
    approximate-model matrices receive the same shape, Gram-identity, and
    declared-label-order guarantees as every other row-matrix consumer.
    The import is local because recursive owns OrbitCoverModel.
    """
    from sage_indefinite_port.indefinite.recursive import (
        _isometry_from_row_matrix as checked_row_isometry,
    )

    return checked_row_isometry(lattice, lattice, rows)



def square_divisors(integer: int) -> tuple[int, ...]:
    r"""The positive c with c^2 dividing the nonzero integer."""
    value = abs(int(integer))
    if value == 0:
        raise ValueError("zero has infinitely many square divisors: covers of norm 0 are by primitive isotropic vectors only")
    return tuple(int(divisor) for divisor in SageZZ(value).divisors() if value % (int(divisor) ** 2) == 0)


def _owned_integer_square(ring, value) -> RingElement:
    """Construct an owned integral norm without accepting lossy truncation.

    Public 2U consumers receive raw Sage ``ZZ`` values from some search
    callers. Convert those through Python's exact integer representation,
    while rejecting rational or inexact scalars rather than silently applying
    ``int`` to (for example) ``3/2`` or ``2.5``.
    """
    if isinstance(value, int):
        integer = int(value)
    elif isinstance(value, RingElement) and value.parent() in (ring, SageZZ):
        integer = int(value)
        if value != integer:
            raise ValueError(f"the prescribed square must be integral, got {value}")
    else:
        raise TypeError(f"the prescribed square must be an integer scalar over {ring}, got {value}")
    return ring(integer)


def _source_random_matrix_perturbation(rank: int):
    """GetRandomMatrixPerturbation (IndefiniteFormsFundamental.g:11-30).

    Random signed permutation, or one elementary entry replacement, exactly
    as in the pinned GAP source (including the diagonal replacement case).
    """
    step = identity_matrix(SageZZ, rank)
    if randrange(2) == 0:
        step = matrix(SageZZ, rank, rank)
        for row, column in enumerate(sample(range(rank), rank)):
            step[row, column] = choice((-1, 1))
    else:
        step[randrange(rank), randrange(rank)] = choice((-1, 1))
    return step


def _source_find_isotropic(gram):
    """INDEF_FindIsotropic/FindIsotropic, pinned ApproximateModels.h:130-166
    and Isotropic.h:76-198,252-377.

    The source first searches its two indefinite-LLL passes for an
    isotropic vector, then its isotropic-mode indefinite reduction, and
    finally the exact solver on the reduced form: FindIsotropicExact
    (Isotropic.h:279-306), translated in exact_isotropic for every rank.
    As in the source, the intermediate result is an integer coordinate
    vector. Only the returned hyperbolic pair becomes lattice elements.
    """
    rank = gram.nrows()
    if gram.rank() < rank:
        # INDEF_FindIsotropic: the first row of the integral nullspace.
        return vector_space_element(SageZZ, gram.right_kernel_matrix().row(0))

    def raise_vector(coordinates):
        if coordinates is None:
            return None
        rational = vector_space_element(SageQQ, coordinates)
        integral = rational * rational.denominator()
        entries = tuple(int(entry) for entry in integral)
        content = gcd(*entries)
        if content == 0:
            raise ArithmeticError("indefinite isotropic search returned the zero vector")
        result = vector_space_element(SageZZ, [entry // content for entry in entries])
        if result.dot_product(gram * result) != 0:
            raise ArithmeticError("indefinite isotropic search returned an invalid witness")
        return result

    if rank <= 2:
        # FindIsotropicExact: skip the LLL and reduction branches.
        return raise_vector(find_isotropic_exact(gram))

    def easy_vector(form):
        # get_isotropic_easy_method, Isotropic.h:76-112.
        for i in range(rank):
            if form[i, i] == 0:
                entries = [0] * rank
                entries[i] = 1
                return vector_space_element(SageZZ, entries)
        for i in range(rank):
            for j in range(i + 1, rank):
                if form[i, j] == 0 and form[i, i] + form[j, j] == 0:
                    entries = [0] * rank
                    entries[i] = entries[j] = 1
                    return vector_space_element(SageZZ, entries)
        return None

    transform = identity_matrix(SageZZ, rank)
    working = gram
    # FindIsotropic_LLL_nfixed makes exactly two attempts.
    for _ in range(2):
        step_transform = identity_matrix(SageZZ, rank)
        step_working = working
        previous_norm = _gram_l1(step_working)
        immediate = easy_vector(step_working)
        if immediate is not None:
            return raise_vector(transform.transpose() * immediate)
        while True:
            step, reduced, isotropic = _indefinite_lll(step_working)
            if isotropic is not None:
                return raise_vector(transform.transpose() * step_transform.transpose() * isotropic)
            step_transform = step * step_transform
            step_working = reduced
            immediate = easy_vector(step_working)
            if immediate is not None:
                return raise_vector(transform.transpose() * step_transform.transpose() * immediate)
            next_norm = _gram_l1(step_working)
            if next_norm >= previous_norm:
                break
            previous_norm = next_norm
            perturbation = _random_signed_permutation(rank)
            step_transform = perturbation * step_transform
            step_working = perturbation * step_working * perturbation.transpose()
        perturbation = _random_signed_permutation(rank)
        transform = perturbation * transform
        working = perturbation * working * perturbation.transpose()

    reduction, reduced_gram, isotropic = _compute_indefinite_reduction(gram, look_for_isotropic=True)
    if isotropic is not None:
        return raise_vector(isotropic)
    exact_coordinates = find_isotropic_exact(reduced_gram)
    if exact_coordinates is None:
        return None
    return raise_vector(reduction.transpose() * exact_coordinates)


def _source_hyperbolic_pair(lattice: Lattices.ParentMethods):
    """GetHyperbolicPlane, polyhedral_common ApproximateModels.h:996-1180.

    The pinned C++ revision of IndefiniteForms.g:465-534 maintains a set
    of sampled isotropic vectors, retrying an insertion at most 100 times
    for distinctness. For each sample it computes the least nonzero
    pairing with the set, retaining the first minimal-pairing pair. It
    exits on 10 subsequent samples at that level or 20 above it, with
    both counters reset by a strict improvement. These are the *C++*
    stopping branches, not the older GAP revision's 100 occurrences.

    The pair has positive integral pairing; unlike a U-splitting, that
    pairing need not be 1 or the divisibility of the first vector.
    """
    gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
    rank = gram.nrows()
    if rank < 2:
        raise ValueError("an isotropic hyperbolic pair needs rank at least two")
    labels = tuple(lattice.module_generating_set())
    witnesses = {}
    perturbation = identity_matrix(SageZZ, rank)
    best_pairing = None
    best_pair = None
    at_level = 0
    no_change = 0

    def insert_vector():
        nonlocal perturbation
        repeats = 0
        while True:
            perturbation *= _source_random_matrix_perturbation(rank)
            changed_gram = perturbation * gram * perturbation.transpose()
            selected = _source_find_isotropic(changed_gram)
            if selected is None:
                raise ValueError("the source isotropic search found no vector")
            original_row = selected * perturbation
            key = tuple(int(entry) for entry in original_row)
            if original_row.dot_product(gram * original_row) != 0:
                raise ArithmeticError("the perturbed isotropic vector is not isotropic")
            if key not in witnesses or repeats == 100:
                witnesses[key] = original_row
                return original_row
            repeats += 1

    while True:
        candidate = insert_vector()
        minimum = None
        partner = None
        # std::set orders vector coordinates lexicographically, and ties
        # retain the first partner: only strict improvements replace it.
        for key in sorted(witnesses):
            previous = witnesses[key]
            pairing = int(previous.dot_product(gram * candidate))
            if pairing and (minimum is None or abs(pairing) < minimum):
                minimum = abs(pairing)
                partner = previous if pairing > 0 else -previous
        if minimum is None:
            continue
        if best_pairing is None or minimum < best_pairing:
            best_pairing = minimum
            best_pair = (candidate, partner)
            at_level = 0
            no_change = 0
        elif minimum == best_pairing:
            at_level += 1
        else:
            no_change += 1
        if at_level == 10 or no_change == 20:
            result = tuple(lattice.linear_combination({label: lattice.base_ring()(int(value)) for label, value in zip(labels, row, strict=True) if value}) for row in best_pair)
            if any(point.q() != lattice.base_ring().zero() for point in result) or int(lattice.b(result[0], result[1])) != best_pairing:
                raise ArithmeticError("the sampled hyperbolic pair failed its Gram certificate")
            return result


class TwoUFullGeneratingDatum(SageObject):
    r"""Finite generators supplied by the ``2U`` approximate-model theorem.

    Let ``A`` be the finite-generator approximate subgroup, let ``v`` be the
    selected primitive base vector, and let ``r_i`` run through the approximate
    orbit representatives lying in the full ``O(L)``-orbit of ``v``.  Choose
    ``t_i`` with ``t_i(v)=r_i`` and generators of ``P_v=Stab(v)``.  Since the
    approximate-model oracle says every vector in that full orbit is
    ``A``-equivalent to some ``r_i``, the elementary orbit argument gives

    ``O(L) = < A, P_v, t_i >``.

    The object retains precisely those finite families.  It does not claim a
    word algorithm in the approximate subgroup beyond the sourced oracle.
    """

    def __init__(
        self,
        model,
        square,
        base_vector,
        approximate_generators,
        stabilizer_generators,
        orbit_representatives,
        orbit_transporters,
        generating_family,
    ) -> None:
        self._model = model
        self._square = square
        self._base_vector = base_vector
        self._approximate_generators = approximate_generators
        self._stabilizer_generators = stabilizer_generators
        self._orbit_representatives = orbit_representatives
        self._orbit_transporters = orbit_transporters
        self._generating_family = generating_family

    def model(self):
        return self._model

    def lattice(self):
        return self.model().lattice()

    def square(self):
        return self._square

    def base_vector(self):
        return self._base_vector

    def approximate_generators(self):
        return self._approximate_generators

    def stabilizer_generators(self):
        return self._stabilizer_generators

    def orbit_representatives(self):
        return self._orbit_representatives

    def orbit_transporters(self):
        return self._orbit_transporters

    def generating_family(self):
        return self._generating_family

    def __repr__(self) -> str:
        return f"Finite 2U orthogonal generating datum for {self.lattice()}"


class EichlerEnvelope:
    r"""The similarity sigma = N phi: L -> E = U + U + K of Dutour Sikirić's hyperplane-basis step.

    For an odd integral lattice, first use its canonical even index-two
    sublattice L_even = ker(q mod 2).  It has the same rational quadratic
    space; hence its Eichler hyperplane basis supplies the same rational
    isometry, whose restriction to L is cleared by the similarity.  The
    finite-index integral-structure descent is unchanged.

    source: ApproximateModels.h:1215-1289 GetEichlerHyperplaneBasis (two
    GetHyperbolicPlane calls, :995-1180, the second on the integral
    perpendicular of the first; EmbedSimpInv = diag(1/a1, 1, 1/a2, 1, ...),
    so b1 -> a1 e1, b2 -> f1, b3 -> a2 e2, b4 -> f2; Embed is the cleared
    inverse and scal its squared multiplier).  Deviation: the pinned source
    requires an even form (is_eichler_canonical, :246; EnumerateVectorOver-
    Discriminant halves norms).  For odd L the construction runs on the even
    sublattice above, which spans the same rational space.
    """

    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        construction_lattice = lattice
        even_inclusion = None
        if not lattice.is_even():
            labels = tuple(lattice.module_generating_set())
            generators = tuple(lattice.module_generator(label) for label in labels)
            parities = tuple(int(generator.q()) % 2 for generator in generators)
            odd_position = next(index for index, parity in enumerate(parities) if parity)
            odd_generator = generators[odd_position]
            even_subobject = lattice.subobject_on(
                tuple(
                    2 * generator if index == odd_position else generator - parity * odd_generator
                    for index, (generator, parity) in enumerate(zip(generators, parities, strict=True))
                )
            )
            even_inclusion = even_subobject.inclusion()
            construction_lattice = even_inclusion.domain()
            if not construction_lattice.is_even():
                raise ArithmeticError("the kernel of the norm-parity homomorphism must be even")

        first, first_partner = _source_hyperbolic_pair(construction_lattice)
        first_plane = construction_lattice.subobject_on([first, first_partner])
        complement = first_plane.orthogonal_complement()
        complement_lattice = complement.inclusion().domain()
        second, second_partner = _source_hyperbolic_pair(complement_lattice)
        embed = complement.inclusion()
        second, second_partner = embed(second), embed(second_partner)
        planes = construction_lattice.subobject_on([first, first_partner, second, second_partner])
        kernel = planes.orthogonal_complement()
        kernel_lattice = kernel.inclusion().domain()
        hyperbolic_plane = Lattices(lattice.base_ring())("U")
        self.lattice = lattice
        self.envelope = hyperbolic_plane + hyperbolic_plane + kernel_lattice
        self.model = self.envelope.two_u_eichler_model_from_represented_biproduct()
        e_1, f_1, e_2, f_2 = self.model.hyperbolic_basis()
        complement_embedding = self.model.complement_inclusion()
        first_pairing = construction_lattice.b(first, first_partner)
        second_pairing = construction_lattice.b(second, second_partner)
        source = (first, first_partner, second, second_partner) + tuple(kernel.inclusion()(vector) for vector in kernel_lattice.module_generators())
        target = (first_pairing * e_1, f_1, second_pairing * e_2, f_2) + tuple(complement_embedding(vector) for vector in kernel_lattice.module_generators())
        fraction_map = lattice.base_ring().fraction_field_map()
        # ApproximateModels.h:1256-1272: FullBasis (the rows of ``source``) is a
        # rational basis of L tensor Q, and preEmbed = Inverse(EmbedSimpInv *
        # FullBasis) writes each basis vector of L in the envelope basis. That
        # inverse is the coordinate solve below; no sublattice spanned by
        # ``source`` is constructed.
        labels = tuple(construction_lattice.module_generating_set())
        envelope_labels = tuple(self.envelope.module_generating_set())
        source_rows = matrix(SageQQ, [[int(vector.to_vector()(label)) for label in labels] for vector in source])
        target_rows = matrix(SageQQ, [[int(vector.to_vector()(label)) for label in envelope_labels] for vector in target])
        images = source_rows.inverse() * target_rows
        rational_construction = construction_lattice.base_change(fraction_map)
        rational_envelope = self.envelope.base_change(fraction_map)
        field = rational_envelope.base_ring()
        image_rows = {label: images.row(position) for position, label in enumerate(labels)}
        construction_to_envelope = rational_construction.Isom(rational_envelope)(
            lambda label: rational_envelope.linear_combination(
                {envelope_labels[column]: field(int(value.numerator())) / field(int(value.denominator())) for column, value in enumerate(image_rows[label]) if value}
            )
        )
        if even_inclusion is None:
            self.rational_isometry = construction_to_envelope
        else:
            rational_original = lattice.base_change(fraction_map)
            rational_even_inclusion = even_inclusion.base_change(fraction_map)
            even_to_original = rational_construction.Isom(rational_original)(lambda label: rational_even_inclusion(rational_construction.module_generator(label)))
            self.rational_isometry = construction_to_envelope * even_to_original.inverse()
        self.similarity = self.rational_isometry.integral_similarity(lattice, self.envelope)
        self.scale = SageZZ(self.similarity.scale())
        self.fraction_map = fraction_map

    @cached_method
    def stabilizer_of_lattice(self) -> IntegralStructureAction:
        r"""The action of the rational group generated by A_E on the integral structures of E tensor Q, at phi(L)."""
        rational_envelope = self.rational_isometry.codomain()
        # The rational extensions of the approximate generators, from their
        # integer matrices on E (source_generator_matrices) rather than by
        # base change of each integral isometry.
        generators = tuple(_isometry_from_row_matrix(rational_envelope, generator) for generator in OrbitCoverModel.source_generator_matrices(self.model))
        rational_lattice = self.rational_isometry.domain()
        image = rational_envelope.integral_structure(tuple(self.rational_isometry(rational_lattice.module_generator(label)) for label in self.lattice.module_generating_set()))
        # sigma = N phi is integral, so phi(L) lies in (1/N)E, which the integral
        # generators of A_E preserve: the stable containing module.
        inverse_scale = _owned_engine_element(
            rational_envelope.base_ring(), SageQQ(1) / SageQQ(self.scale)
        )
        labels = rational_envelope.module_generating_set()
        scaled_generators = tuple(rational_envelope.scalar_multiple(inverse_scale, rational_envelope.module_generator(label)) for label in labels)
        bound = rational_envelope.integral_structure(scaled_generators)
        return IntegralStructureAction(rational_envelope.O().subgroup(generators), image).with_stable_containing_bound(bound)

    @cached_method
    def right_coset_transversal(self) -> tuple[LatticeIsometryMethods, ...]:
        r"""Rational isometries t with A_E = union of G_L t: ListCoset of ApproximateModels.h:1318-1332.

        The pointed orbit G/G_L of phi(L) has the representative map
        p |-> t_p with t_p(phi(L)) = p, so G = union of t_p G_L and
        G = union of G_L t_p^-1.  The inverses are used only as rational
        isometries acting on vectors of E tensor Q (get_vector_representatives,
        :1374-1390); they are not re-admitted as words of the generated
        subgroup, whose membership the preamble does not decide.
        """
        coset_space = self.stabilizer_of_lattice().right_cosets()
        representative = coset_space.representative_map()
        return tuple(~representative(point) for point in coset_space)


_EASY_GENERATORS: dict[int, tuple[LatticeIsometryMethods, ...]] = {}
_GENERATOR_MATRICES: dict[int, tuple] = {}
"""The approximate generators of each live 2U model as integer matrices, keyed by the model's identity."""
"""GRPeasy per live 2U model, keyed by the model's identity (see source_easy_generators)."""


class OrbitCoverModel:
    r"""The approximate subgroup A_L of O(L) and the covering lists of the Eichler recursion."""

    @staticmethod
    def source_easy_generators(model: TwoUEichlerModel) -> tuple[LatticeIsometryMethods, ...]:
        r"""GRPeasy of ``model``, computed once per model.

        It serves both the approximate group and the class-orbit refinement,
        as the source computes it once (ApproximateModels.h:1353).  The model
        is not hashable, so the datum is held by the model's identity for the
        model's lifetime.
        """
        key = id(model)
        if key not in _EASY_GENERATORS:
            _EASY_GENERATORS[key] = OrbitCoverModel._computed_easy_generators(model)
            weakref.finalize(model, _EASY_GENERATORS.pop, key)
        return _EASY_GENERATORS[key]

    @staticmethod
    def _computed_easy_generators(model: TwoUEichlerModel) -> tuple[LatticeIsometryMethods, ...]:
        r"""Port GetEasyIsometries, IndefiniteForms.g:593–680.

        Extract the actual connected Gram blocks, extend positive-definite
        automorphisms by the identity and permute all exactly equal blocks.
        The resulting group is precisely GRPeasy used both to generate the
        approximate group and to refine its discriminant-class covering.
        """
        lattice = model.lattice()
        ring = lattice.base_ring()
        complement = model.orthogonal_complement()
        inclusion = model.complement_inclusion()
        basis = model.hyperbolic_basis() + tuple(inclusion(point) for point in complement.module_generators())
        gram = matrix(SageZZ, [[int(lattice.b(left, right)) for right in basis] for left in basis])

        graph = Graph(len(basis))
        graph.add_edges((first, second) for first in range(len(basis)) for second in range(first + 1, len(basis)) if gram[first, second])
        components = tuple(tuple(sorted(component)) for component in graph.connected_components())
        blocks = tuple(gram.matrix_from_rows_and_columns(component, component) for component in components)

        # Pinned polyhedral_common ApproximateModels.h:824-964: insert the
        # identity, full automorphisms of *both signs* of definite block,
        # the swap and negative identity of each aU block, and every
        # pairwise isometric definite/aU-block interchange. In particular,
        # equal Gram matrices are not the only interchangeable blocks.
        # The source inserts the identity first (ApproximateModels.h:842); it
        # generates nothing, and every generator costs one evaluation of each
        # downstream finite action, so it is omitted.
        easy_generators = []
        definite_blocks = []
        for component, block in zip(components, blocks, strict=True):
            definite = block.is_positive_definite() or (-block).is_positive_definite()
            block_lattice = Lattices(ring)([[ring(int(block[row, column])) for column in range(block.ncols())] for row in range(block.nrows())]) if definite else None
            definite_blocks.append(block_lattice)
            if definite:
                labels = tuple(block_lattice.module_generating_set())
                for isometry in block_lattice.O().group_generators():
                    images = list(basis)
                    for position, label in enumerate(labels):
                        coordinates = isometry(block_lattice.module_generator(label)).to_vector()
                        images[component[position]] = sum(
                            (
                                lattice.scalar_multiple(coordinates(inner_label), basis[component[inner_position]])
                                for inner_position, inner_label in enumerate(labels)
                                if coordinates(inner_label)
                            ),
                            lattice.zero(),
                        )
                    easy_generators.append(model._isometry_on_model_basis(tuple(images)))
            elif len(component) == 2 and block[0, 0] == 0 and block[1, 1] == 0:
                # INDEF_FORM_AutU, including the possibly nonunimodular aU.
                swap_images = list(basis)
                swap_images[component[0]], swap_images[component[1]] = (basis[component[1]], basis[component[0]])
                easy_generators.append(model._isometry_on_model_basis(tuple(swap_images)))
                negate_images = list(basis)
                for index in component:
                    negate_images[index] = -basis[index]
                easy_generators.append(model._isometry_on_model_basis(tuple(negate_images)))

        for first in range(len(blocks)):
            for second in range(first + 1, len(blocks)):
                left, right = blocks[first], blocks[second]
                if left.nrows() != right.nrows():
                    continue
                images = list(basis)
                if definite_blocks[first] is not None and definite_blocks[second] is not None:
                    # INDEF_FORM_TestEquivalence_PosNeg uses an exact definite
                    # Isom witness, not literal equality of two presentations.
                    left_positive = left.is_positive_definite()
                    if left_positive != right.is_positive_definite():
                        continue
                    left_lattice, right_lattice = definite_blocks[first], definite_blocks[second]
                    homset = left_lattice.Isom(right_lattice)
                    if homset.is_empty():
                        continue
                    witness = homset.an_element()
                    inverse = witness.inverse()
                    left_labels = tuple(left_lattice.module_generating_set())
                    right_labels = tuple(right_lattice.module_generating_set())
                    for position, label in enumerate(left_labels):
                        coordinates = witness(left_lattice.module_generator(label)).to_vector()
                        images[components[first][position]] = sum(
                            (
                                lattice.scalar_multiple(coordinates(other), basis[components[second][other_position]])
                                for other_position, other in enumerate(right_labels)
                                if coordinates(other)
                            ),
                            lattice.zero(),
                        )
                    for position, label in enumerate(right_labels):
                        coordinates = inverse(right_lattice.module_generator(label)).to_vector()
                        images[components[second][position]] = sum(
                            (
                                lattice.scalar_multiple(coordinates(other), basis[components[first][other_position]])
                                for other_position, other in enumerate(left_labels)
                                if coordinates(other)
                            ),
                            lattice.zero(),
                        )
                elif len(components[first]) == 2 and left[0, 0] == left[1, 1] == right[0, 0] == right[1, 1] == 0 and left[0, 1] == right[0, 1]:
                    # INDEF_FORM_IsU and its equal-scale identity transporter.
                    for position in range(2):
                        images[components[first][position]] = basis[components[second][position]]
                        images[components[second][position]] = basis[components[first][position]]
                else:
                    continue
                easy_generators.append(model._isometry_on_model_basis(tuple(images)))

        return tuple(easy_generators)

    @staticmethod
    def source_generating_family(model: TwoUEichlerModel) -> IndexedFamily[int, LatticeIsometryMethods]:
        r"""IndefiniteForms.g:593–680,175–320: GetEasyIsometries and GetApproximateGroup.

        Use the graph of nonzero Gram entries in the represented 2U + K basis.
        GetEasyIsometries inserts both signs of definite-block automorphisms,
        the aU-block involutions and the interchanges of isometric blocks.
        The source then adds both SL_2(ZZ) actions and
        transvections for the four coordinate isotropic vectors.

        source: ApproximateModels.h:406-529 GetApproximateGroup, :823-967
        GetEasyIsometries.  Deviation: both pinned revisions evaluate
        GetApproximateGroup before SetListClassesOrbitwise(GRPeasy) installs
        GRPstart (ApproximateModels.h:1309 before :1355; IndefiniteForms.g
        builds ApproximateGroup at record construction, :424), so their
        approximate group omits GRPeasy while their covering classes are
        reduced modulo GRPeasy.  A class c' = eps(c) dropped by that
        reduction is then reached only through eps, outside the group.  Here
        GRPeasy is a generating subfamily, which makes the GRPeasy-orbit
        reduction of :307-405 a valid covering reduction; every element is an
        integral isometry of the model, so the generated group is still a
        subgroup of O(L) containing the source's.
        """
        lattice = model.lattice()
        matrices = OrbitCoverModel.source_generator_matrices(model)
        return finite_indexed_family(
            finite_ordered_set(tuple(range(len(matrices)))),
            lambda index: _isometry_from_row_matrix(lattice, matrices[index]),
            name=f"Source-defined Eichler approximate generators in O({lattice})",
        )

    @staticmethod
    def source_generator_matrices(model: TwoUEichlerModel) -> tuple:
        r"""The approximate generators as integer matrices on the basis of the model lattice E, once per model.

        Row u of a matrix is the coordinate vector of the image of the u-th
        basis vector of E (the source's row convention, P Q P^T = Q).
        GRPeasy and the two SL_2(ZZ) actions are read off once from their
        isometries.  The transvections are the source's matrix formula
        INDEF_FORM_Eichler_Transvection (ApproximateModels.h:78-127),
        E_{f,x}(y) = y + (y,x) f - (x,x)/2 (y,f) f - (y,f) x, which with Gram
        matrix G is I + (G x) f - (x G x / 2)(G f) f - (G f) x in outer
        products.
        """
        key = id(model)
        if key not in _GENERATOR_MATRICES:
            _GENERATOR_MATRICES[key] = OrbitCoverModel._computed_generator_matrices(model)
            weakref.finalize(model, _GENERATOR_MATRICES.pop, key)
        return _GENERATOR_MATRICES[key]

    @staticmethod
    def _computed_generator_matrices(model: TwoUEichlerModel) -> tuple:
        lattice = model.lattice()
        labels = tuple(lattice.module_generating_set())
        gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
        rank = len(labels)

        def coordinates(vector: Lattices.ElementMethods):
            values = vector.to_vector()
            return vector_space_element(SageZZ, [int(values(label)) for label in labels])

        def isometry_matrix(isometry: LatticeIsometryMethods):
            return matrix(SageZZ, [coordinates(isometry(lattice.module_generator(label))) for label in labels])

        def transvection_matrix(point: Lattices.ElementMethods, vector: Lattices.ElementMethods):
            f, x = coordinates(point), coordinates(vector)
            half_norm = x.dot_product(gram * x) // 2
            gram_f, gram_x = gram * f, gram * x
            return identity_matrix(SageZZ, rank) + gram_x.column() * f.row() - half_norm * gram_f.column() * f.row() - gram_f.column() * x.row()

        complement = model.orthogonal_complement()
        inclusion = model.complement_inclusion()
        isotropic = model.hyperbolic_basis()
        e, f, e_prime, f_prime = isotropic
        embedded = tuple(inclusion(point) for point in complement.module_generators())
        # A basis of the isotropic point's perpendicular, less the point
        # itself: E_{f,f} is the identity (q(f) = 0 and (y,f)f - (y,f)f = 0),
        # and E_{f,x} is additive in x (ApproximateModels.h:62-76), so the
        # source's E_{f, basis of f-perp} (:505-520) generate the same group.
        perpendicular = ((e_prime, f_prime), (e_prime, f_prime), (e, f), (e, f))
        linear = tuple(model.special_linear_group().group_generators())
        matrices = [isometry_matrix(isometry) for isometry in OrbitCoverModel.source_easy_generators(model)]
        matrices += [isometry_matrix(model.left_action(element)) for element in linear]
        matrices += [isometry_matrix(model.right_action(element)) for element in linear]
        matrices += [transvection_matrix(point, vector) for point, subbasis in zip(isotropic, perpendicular, strict=True) for vector in (*subbasis, *embedded)]
        for generator in matrices:
            generator.set_immutable()
        return tuple(matrices)

    @staticmethod
    def split_generating_family(model: TwoUEichlerModel) -> IndexedFamily[int, LatticeIsometryMethods]:
        r"""The pinned split 2U approximate group, including the easy block actions."""
        return OrbitCoverModel.source_generating_family(model)

    @staticmethod
    def split_covering_vectors(model: TwoUEichlerModel, square: RingElement) -> IndexedFamily[DiscriminantQuadraticModules.ElementMethods, Lattices.ElementMethods]:
        r"""The source's orbit-refined covering vectors in the split 2U model.

        If x in A_K has order d and y in K dual lifts x, put k=d*y and
        m=(square-k^2)/(2d^2). Then d*e+d*m*f+k has the prescribed square
        and divided class. Its pairing with f is d, while the class order
        forces primitivity and divisibility exactly d. IndefiniteForms.g
        SetListClassesOrbitwise first keeps one class per orbit of GRPeasy:
        the other classes are already carried to these by the easy subgroup.

        source: ApproximateModels.h:548-717 EnumerateVectorOverDiscriminant
        and :718-758 GetCoveringOrbitRepresentatives (per content, below).
        The source's admissibility test A[e3]/2 = X/2 mod d followed by its
        w-correction accepts exactly the classes with
        q_A(x) = X/d^2 mod 2Z: e3 = d y and e3 + d w = d(y + w) have the same
        norm modulo 2d^2, since b(y, w) is integral and q(w) even, so the
        correction is solvable only with w = 0.  Those classes are
        ``covering_discriminant_classes(square)``, and then m = (X - q(k)) /
        (2 d^2) is integral.  Deviation: the source places d, d u in the
        second hyperbolic plane (e', f'); here they are in the first plane
        (e, f).  GRPeasy interchanges the two equal U blocks, so the two
        vectors are equivalent under the approximate group and the covering
        guarantee is the same.
        """
        lattice = model.lattice()
        ring = lattice.base_ring()
        # This public consumer can be reached with a raw Sage ZZ scalar from
        # the approximate-covering search.  Cross the scalar boundary through
        # Python's exact integer representation before constructing the owned
        # ring element; owned-ring coercion deliberately rejects raw Sage ZZ.
        square = _owned_integer_square(ring, square)
        complement = model.orthogonal_complement()
        discriminant = complement.discriminant_group()
        inclusion = model.complement_inclusion()
        e, f, _e_prime, _f_prime = model.hyperbolic_basis()

        # eClass2 = Ginv * eClass1 (:610-625): the correlation K -> K^# has
        # matrix G on the selected bases (lattices.py correlation_morphism), so
        # a dual lift with dual-basis coordinates c is y = G^-1 c in K tensor Q.
        complement_labels = tuple(complement.module_generating_set())
        gram_inverse = _engine_component_matrix(complement.gram_tensor()).change_ring(SageQQ).inverse()

        def scaled_lift(discriminant_class: DiscriminantQuadraticModules.ElementMethods) -> tuple[RingElement, Lattices.ElementMethods]:
            r"""The order d of the class and k = d*y in K for a dual lift y (eClass3 of :610-625)."""
            discriminant_class = discriminant(discriminant_class)
            order = ring(int(discriminant_class.additive_order()))
            dual_lift = discriminant.dual_lattice_lift(discriminant_class)
            dual_coordinates = dual_lift.to_vector()
            dual_labels = tuple(dual_lift.parent().module_generating_set())
            scaled = int(order) * gram_inverse * vector_space_element(SageQQ, [int(dual_coordinates(label)) for label in dual_labels])
            if scaled.denominator() != 1:
                raise ArithmeticError("the order of a discriminant class must clear the denominators of its dual lift")
            return order, complement.linear_combination({label: ring(int(value)) for label, value in zip(complement_labels, scaled, strict=True) if value})

        classes = tuple(complement.covering_discriminant_classes(square))
        lifts = {discriminant_class: scaled_lift(discriminant_class) for discriminant_class in classes}

        def representative(discriminant_class: DiscriminantQuadraticModules.ElementMethods) -> Lattices.ElementMethods:
            order, k = lifts[discriminant_class]
            numerator = square - complement.q(k)
            denominator = ring(2) * order * order
            coefficient = numerator // denominator
            if denominator * coefficient != numerator:
                raise ArithmeticError("a covering discriminant class must give an integral hyperbolic coefficient")
            return lattice.scalar_multiple(order, e) + lattice.scalar_multiple(order * coefficient, f) + inclusion(k)

        # SetListClassesOrbitwise(GRPeasy), ApproximateModels.h:313-405.  A
        # class y of K^#/K is extended by zero to E = 2U + K, and the class of
        # g(y) in E^#/E is read off as its coordinates modulo 1 (VectorModulo1
        # of QmatInv * v, :337-356).  A_K -> A_E is an isomorphism since 2U is
        # unimodular, so every image is one of the listed classes, even when a
        # GRPeasy element moves a hyperbolic block of K onto a distinguished
        # U.  With y = k/d this is the vector of E-coordinates of g(k) modulo d.
        labels = tuple(lattice.module_generating_set())

        def residue(vector: Lattices.ElementMethods, order: RingElement) -> tuple[SageQQ, ...]:
            coordinates = vector.to_vector()
            modulus = int(order)
            return tuple(SageQQ(int(coordinates(label)) % modulus) / modulus for label in labels)

        positions = {residue(inclusion(k), order): position for position, (order, k) in enumerate(lifts[point] for point in classes)}
        permutations = [
            [positions[residue(isometry(inclusion(k)), order)] + 1 for order, k in (lifts[point] for point in classes)]
            for isometry in OrbitCoverModel.source_easy_generators(model)
        ]
        # DecomposeOrbitPoint_FullRepr (:394-400): one class per orbit of the
        # permutation group of GRPeasy on the classes.
        orbits = PermutationGroup(permutations, domain=range(1, len(classes) + 1)).orbits() if classes else ()
        selected = finite_ordered_set(tuple(classes[orbit[0] - 1] for orbit in orbits))
        return finite_indexed_family(selected, representative, name=f"Eichler covering vectors of square {square} in {lattice}")

    @staticmethod
    def split_full_generating_data(model: TwoUEichlerModel, square: RingElement) -> TwoUFullGeneratingDatum:
        r"""Assemble O(L) from the approximate subgroup, one pointed stabilizer and exact witnesses.

        Every image of the base vector is approximately equivalent to one
        retained covering vector. Dividing an arbitrary isometry by that
        approximate element and the retained witness leaves an element of
        the base stabilizer. This proves generation by the three families.
        """
        from sage_indefinite_port.indefinite.recursive import EichlerRecursion

        recursion = EichlerRecursion()
        representatives = model.covering_vector_representatives(square)
        # The family is explicitly indexed by a finite ordered set. Testing
        # that index set avoids routing a decidable finite emptiness question
        # through cardinality() on its represented torsion-module values;
        # the pinned provider's generic cardinality path currently rejects
        # those values before reaching the family-level finite decision.
        representative_labels = tuple(representatives.index_set())
        if not representative_labels:
            raise ValueError("full-group assembly requires a nonempty primitive covering family")
        base_label = representative_labels[0]
        base_vector = representatives.value(base_label)
        witnesses: dict[DiscriminantQuadraticModules.ElementMethods, LatticeIsometryMethods] = {}
        for label in representative_labels:
            witness = recursion.pointed_isometry(base_vector, representatives.value(label))
            if witness is not None:
                witnesses[label] = witness
        labels = finite_ordered_set(tuple(witnesses))
        orbit_representatives = finite_indexed_family(labels, representatives.value)
        orbit_witnesses = finite_indexed_family(labels, witnesses.__getitem__)
        approximate = model.approximate_generating_family()
        stabilizer = finite_ordered_set(recursion.pointed_automorphism_generators(base_vector))
        generators = tuple(approximate) + tuple(stabilizer) + tuple(orbit_witnesses)
        family = finite_indexed_family(finite_ordered_set(tuple(range(len(generators)))), generators.__getitem__)
        # Match split_covering_vectors: public callers may supply a raw Sage
        # integer, which the owned coefficient ring intentionally rejects.
        # The prescribed norm is integral, so Python's exact integer is the
        # boundary representation before constructing the owned scalar.
        owned_square = _owned_integer_square(model.lattice().base_ring(), square)
        return TwoUFullGeneratingDatum(model, owned_square, base_vector, approximate, stabilizer, orbit_representatives, orbit_witnesses, family)

    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        self.lattice = lattice
        if lattice.is_even() and lattice.splits_two_hyperbolic_planes():
            self.envelope = None
            self.model = lattice.two_u_eichler_model_from_represented_biproduct()
        else:
            self.envelope = EichlerEnvelope(lattice)
            self.model = self.envelope.model

    def approximate_generators(self) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of A_L: the Eichler subgroup of the 2U model, preserving L when L is reached through the envelope.

        source: ApproximateModels.h:1309-1349 (LinearSpace_Stabilizer_RightCoset
        of Embed under the model group, conjugated by Embed) and :1401-1408.
        """
        if self.envelope is None:
            return tuple(self.model.approximate_generating_family())
        # :1337-1349 forms fGen = Embed * eGen * EmbedInv for each stabilizer
        # generator and requires it to be integral.  In rows, with R the
        # matrix of phi (row i the E-coordinates of phi(e_i)) and G that of a
        # generator, phi^-1 g phi has matrix R G R^-1.  Equal matrices give
        # equal generators and the identity generates nothing; a Schreier
        # family contains both, so each distinct matrix is kept once.
        envelope = self.envelope
        rational = envelope.rational_isometry
        rational_lattice, rational_envelope = rational.domain(), rational.codomain()
        lattice_labels = tuple(self.lattice.module_generating_set())
        envelope_labels = tuple(rational_envelope.module_generating_set())

        def coordinates(vector: Lattices.ElementMethods) -> list[SageQQ]:
            values = vector.to_vector()
            return [SageQQ(int(values(label).numerator())) / int(values(label).denominator()) for label in envelope_labels]

        embedding = matrix(SageQQ, [coordinates(rational(rational_lattice.module_generator(label))) for label in lattice_labels])
        embedding_inverse = embedding.inverse()
        identity = identity_matrix(SageQQ, len(lattice_labels))
        conjugates: dict[tuple[int, ...], LatticeIsometryMethods] = {}
        for generator in envelope.stabilizer_of_lattice().stabilizer().group_generators():
            conjugate = embedding * matrix(SageQQ, [coordinates(generator(rational_envelope.module_generator(label))) for label in envelope_labels]) * embedding_inverse
            if conjugate.denominator() != 1:
                raise ArithmeticError("a stabilizer generator of phi(L) must conjugate to an integral isometry of L")
            if conjugate == identity:
                continue
            key = tuple(int(entry) for entry in conjugate.list())
            if key not in conjugates:
                conjugates[key] = _isometry_from_row_matrix(self.lattice, conjugate)
        return tuple(conjugates.values())

    def _pullback_envelope_vector(self, point: Lattices.ElementMethods, coset: LatticeIsometryMethods) -> Lattices.ElementMethods | None:
        """Lift a rational coset image from the envelope through the integral similarity, if possible."""
        envelope = self.envelope
        if envelope is None:
            raise ValueError("an envelope coset requires a nonsplit model")
        envelope_lattice = self.model.lattice()
        rational_envelope = envelope.rational_isometry.codomain()
        labels = tuple(envelope_lattice.module_generating_set())
        coordinates = point.to_vector()
        rational_point = rational_envelope.linear_combination({label: rational_envelope.base_ring()(int(coordinates(label))) for label in labels if coordinates(label)})
        moved = coset(rational_point).to_vector()
        if any(moved(label).denominator() != 1 for label in labels):
            return None
        integral_point = envelope_lattice.linear_combination({label: envelope_lattice.base_ring()(int(moved(label))) for label in labels if moved(label)})
        if not envelope.similarity.is_in_image(integral_point):
            return None
        return envelope.similarity.lift(integral_point)

    def one_representative(self, square: int) -> Lattices.ElementMethods | None:
        r"""Port ``GetOneOrbitRepresentative`` (IndefiniteForms.g:416-423, 801-824).

        First try the distinguished hyperbolic representative and the envelope
        cosets. If it misses the original lattice, the source falls back to
        the complete covering list, including nonprimitive multiples when the
        square is nonzero. At the least positive represented norm, every
        representing vector is primitive.
        """
        if self.envelope is None:
            if square % 2:
                return None
            e, f, _e_prime, _f_prime = self.model.hyperbolic_basis()
            return e + self.lattice.scalar_multiple(self.lattice.base_ring()(square // 2), f)

        # source: ApproximateModels.h:1436-1472 GetOneOrbitRepresentative (envelope case),
        # calling :759-791 on the scaled norm Xscal = X * scal.
        envelope = self.envelope
        envelope_square = envelope.scale**2 * square
        if envelope_square % 2:
            # :766-776 and :1449-1451: the even envelope has no vector of odd
            # norm Xscal, so the source returns no representative and does not
            # reach the covering fallback.
            return None
        e, f, _e_prime, _f_prime = self.model.hyperbolic_basis()
        candidate = e + self.model.lattice().scalar_multiple(self.model.lattice().base_ring()(int(envelope_square // 2)), f)
        # :1456-1461 tests this representative against every coset of the
        # integral-structure stabilizer and keeps the first lift.  The
        # identity coset is necessarily in that transversal.  Test the exact
        # identity-coset equation before constructing the full finite
        # action; failure still follows the source's exhaustive search.
        if envelope.similarity.is_in_image(candidate):
            witness = envelope.similarity.lift(candidate)
            if int(witness.q()) != square:
                raise ArithmeticError("the identity-coset lift has the wrong norm")
            return witness
        for coset in envelope.right_coset_transversal():
            lifted = self._pullback_envelope_vector(candidate, coset)
            if lifted is not None:
                return lifted
        # :1466-1471: the first vector of the covering list of norm X.
        if square == 0:
            covers = self.primitive_covering_vectors(0)
            return covers[0] if covers else None
        for content in square_divisors(square):
            covers = self.primitive_covering_vectors(square // (content * content))
            if covers:
                return self.lattice.scalar_multiple(content, covers[0])
        return None

    def primitive_covering_vectors(self, square: int) -> tuple[Lattices.ElementMethods, ...]:
        r"""Primitive vectors of the square meeting every O(L)-orbit of primitive vectors of that square.

        source: ApproximateModels.h:1409-1435 GetCoveringOrbitRepresentatives
        with :1363-1392 get_vector_representatives.  The source lists vectors
        of every content in L, through every square divisor of Xscal
        (:735-753), but for X = 0 only E-primitive vectors (:729-731).  The
        frozen consumer asks for L-primitive vectors: sigma(v) = c p with p
        E-primitive and c the exact order of p in E / sigma(L), so c runs
        over the divisors of that exponent with c^2 | Xscal.  For X != 0 this
        is the source's list restricted to L-primitive vectors; for X = 0 it
        also covers L-primitive v whose image has E-content c > 1, which the
        source's X = 0 branch omits.
        """
        if self.envelope is None:
            return tuple(self.model.covering_vector_representatives(square))
        envelope = self.envelope
        # sigma(v) = c g p with g in A_E and p a covering representative; writing g = h t with
        # h in the stabilizer of phi(L) and t a right-coset representative, v is A_L-equivalent
        # to sigma^-1(c t p). So t ranges over a transversal of the right cosets G_L t.
        cosets = envelope.right_coset_transversal()
        envelope_square = envelope.scale**2 * square
        found: list[Lattices.ElementMethods] = []
        # sigma(v) for v primitive in L is c times a primitive vector of E with c dividing the
        # exponent of E / sigma(L). For a nonzero square c^2 also divides the square; for norm 0
        # only the exponent bounds c, so every divisor of it is a possible content.
        # Indeed, the class of the primitive E-vector in E/sigma(L) has exact
        # order c: a smaller multiple in sigma(L) would divide primitive v in L.
        exponent = int(envelope.similarity.cokernel().exponent())
        # Enumerate divisors of the finite quotient exponent, not divisors of
        # the potentially enormous norm.  The necessary square-divisibility
        # restriction then selects exactly the same possible contents.
        contents = tuple(int(divisor) for divisor in SageZZ(exponent).divisors() if envelope_square == 0 or envelope_square % (int(divisor) ** 2) == 0)
        # One covering list per reduced square Xscal / c^2: for X = 0 every
        # content reduces to the same square, whose list is computed once.
        representatives = {reduced: tuple(self.model.covering_vector_representatives(reduced)) for reduced in {envelope_square // content**2 for content in contents}}
        for content in contents:
            for primitive in representatives[envelope_square // content**2]:
                point = self.model.lattice().scalar_multiple(content, primitive)
                for representative in cosets:
                    vector = self._pullback_envelope_vector(point, representative)
                    if vector is not None and vector.is_primitive():
                        found.append(vector)
        return tuple(found)

    def covering_vectors(self, square: int, content: int) -> tuple[Lattices.ElementMethods, ...]:
        r"""Vectors of the square and content meeting every O(L)-orbit of such vectors: c times the primitive covers of norm square / c^2."""
        if square % (content * content):
            raise ValueError(f"no vector of content {content} has square {square}")
        return tuple(self.lattice.scalar_multiple(content, vector) for vector in self.primitive_covering_vectors(square // (content * content)))

    def splitting_vector(self) -> Lattices.ElementMethods:
        r"""The first represented positive norm, as in ``GetFirstNorm`` (IndefiniteForms.g:1131-1142).

        source: ApproximateModels.h:1482-1497 GetFirstNorm tries X = 1, 2, ...
        An even lattice represents no odd X, so its odd trials return nothing
        (one_representative returns None for odd Xscal) and are skipped.
        """
        square = 2 if self.lattice.is_even() else 1
        increment = 2 if self.lattice.is_even() else 1
        while True:
            representative = self.one_representative(square)
            if representative is not None:
                return representative
            square += increment


__all__ = [
    "EichlerEnvelope",
    "OrbitCoverModel",
    "TwoUFullGeneratingDatum",
    "square_divisors",
]
