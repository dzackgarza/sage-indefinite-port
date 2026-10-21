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

from math import gcd
from random import choice, randrange, sample

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.categories.eichler_criterion import TwoUEichlerModel
from dzack_research.preamble.categories.group.g_sets import FiniteGSets
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.framed.formed.discriminant_modules import (
    DiscriminantQuadraticModules,
)
from dzack_research.preamble.categories.rational_integral_stabilizers import (
    IntegralStructureAction,
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
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector as vector_space_element
from sage.quadratic_forms.qfsolve import qfsolve
from sage.rings.integer import Integer as SageInteger
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.element import RingElement
from sage.structure.sage_object import SageObject

from sage_indefinite_port.indefinite.indefinite_reduction import (
    _compute_indefinite_reduction,
    _gram_l1,
    _indefinite_lll,
    _random_signed_permutation,
)


def square_divisors(integer: int) -> tuple[int, ...]:
    r"""The positive c with c^2 dividing the nonzero integer."""
    value = abs(int(integer))
    if value == 0:
        raise ValueError("zero has infinitely many square divisors: covers of norm 0 are by primitive isotropic vectors only")
    return tuple(int(divisor) for divisor in SageZZ(value).divisors() if value % (int(divisor) ** 2) == 0)


def _source_random_matrix_perturbation(rank: int):
    """GetRandomMatrixPerturbation (IndefiniteFormsFundamental.g:11-30).

    Random signed permutation, or one elementary shear/diagonal sign.
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
    finally the exact solver on the reduced form. Its exact solver is
    Sage's qfsolve, the same PARI-backed operation used by the preamble's
    lattice.isotropic_vector(); no quadratic-solving algorithm is copied.
    As in the source, the intermediate result is an integer coordinate
    vector. Only the returned hyperbolic pair becomes lattice elements.
    """
    rank = gram.nrows()
    if gram.rank() < rank:
        # INDEF_FindIsotropic: the first row of the integral nullspace.
        return vector_space_element(SageZZ, gram.right_kernel_matrix().row(0))

    def raise_vector(coordinates):
        if isinstance(coordinates, SageInteger):
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
        return raise_vector(qfsolve(gram.change_ring(SageQQ)))

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
    exact_coordinates = qfsolve(reduced_gram.change_ring(SageQQ))
    if isinstance(exact_coordinates, SageInteger):
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
    """

    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        construction_lattice = lattice
        even_inclusion = None
        if not lattice.is_even():
            generators = tuple(lattice.module_generators())
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
        spanned = construction_lattice.subobject_on(list(source))
        # The span basis is echelonized; each basis vector is an integral combination of the
        # independent family source, and its image is the same combination of target.
        labels = tuple(construction_lattice.module_generating_set())
        source_rows = matrix(SageQQ, [[int(vector.to_vector()(label)) for label in labels] for vector in source])
        span_inclusion = spanned.inclusion()
        span_images = []
        for generator in span_inclusion.domain().module_generators():
            image = span_inclusion(generator)
            coefficients = source_rows.solve_left(vector_space_element([int(image.to_vector()(label)) for label in labels]))
            if any(entry.denominator() != 1 for entry in coefficients):
                raise ArithmeticError("a span basis vector must be an integral combination of the spanning family")
            span_images.append(
                sum(
                    (construction_lattice.base_ring()(int(entry)) * image_vector for entry, image_vector in zip(coefficients, target, strict=True)),
                    self.envelope.zero(),
                )
            )
        span_to_envelope = span_inclusion.domain().hom(span_images, codomain=self.envelope).base_change(fraction_map)
        rational_inclusion = spanned.inclusion().base_change(fraction_map)
        rational_construction = construction_lattice.base_change(fraction_map)
        construction_to_envelope = rational_construction.Isom(self.envelope.base_change(fraction_map))(
            lambda label: span_to_envelope(rational_inclusion.lift(rational_construction.module_generator(label)))
        )
        if even_inclusion is None:
            self.rational_isometry = construction_to_envelope
        else:
            rational_original = lattice.base_change(fraction_map)
            rational_even_inclusion = even_inclusion.base_change(fraction_map)
            even_to_original = rational_construction.Isom(rational_original)(lambda label: rational_even_inclusion(rational_construction.module_generator(label)))
            self.rational_isometry = construction_to_envelope * even_to_original.inverse()
        self.similarity = self.rational_isometry.integral_similarity(lattice, self.envelope)
        self.scale = self.similarity.scale()
        self.fraction_map = fraction_map

    @cached_method
    def stabilizer_of_lattice(self) -> IntegralStructureAction:
        r"""The action of the rational group generated by A_E on the integral structures of E tensor Q, at phi(L)."""
        rational_envelope = self.rational_isometry.codomain()
        generators = tuple(generator.base_change(self.fraction_map) for generator in self.model.approximate_generating_family())
        rational_lattice = self.rational_isometry.domain()
        image = rational_envelope.integral_structure(tuple(self.rational_isometry(rational_lattice.module_generator(label)) for label in self.lattice.module_generating_set()))
        # sigma = N phi is integral, so phi(L) lies in (1/N)E, which the integral
        # generators of A_E preserve: the stable containing module.
        inverse_scale = rational_envelope.base_ring().one() / rational_envelope.base_ring()(self.scale)
        labels = rational_envelope.module_generating_set()
        scaled_generators = tuple(rational_envelope.scalar_multiple(inverse_scale, rational_envelope.module_generator(label)) for label in labels)
        bound = rational_envelope.integral_structure(scaled_generators)
        return IntegralStructureAction(rational_envelope.O().subgroup(generators), image).with_stable_containing_bound(bound)


class OrbitCoverModel:
    r"""The approximate subgroup A_L of O(L) and the covering lists of the Eichler recursion."""

    @staticmethod
    def source_easy_generators(model: TwoUEichlerModel) -> tuple[LatticeIsometryMethods, ...]:
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
        easy_generators = [lattice.O().identity()]
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
        """
        lattice = model.lattice()
        complement = model.orthogonal_complement()
        inclusion = model.complement_inclusion()
        easy_generators = OrbitCoverModel.source_easy_generators(model)

        isotropic = model.hyperbolic_basis()
        e, f, e_prime, f_prime = isotropic
        embedded = tuple(inclusion(point) for point in complement.module_generators())
        perpendicular = ((e, e_prime, f_prime), (f, e_prime, f_prime), (e, f, e_prime), (e, f, f_prime))
        transvections = tuple((point, vector) for point, subbasis in zip(isotropic, perpendicular, strict=True) for vector in (*subbasis, *embedded))
        linear = tuple(model.special_linear_group().group_generators())
        easy_count, linear_count = len(easy_generators), len(linear)

        def generator(index: int) -> LatticeIsometryMethods:
            if index < easy_count:
                return easy_generators[index]
            index -= easy_count
            if index < linear_count:
                return model.left_action(linear[index])
            if index < 2 * linear_count:
                return model.right_action(linear[index - linear_count])
            point, vector = transvections[index - 2 * linear_count]
            return model.eichler_transvection(point, vector)

        return finite_indexed_family(
            finite_ordered_set(tuple(range(easy_count + 2 * linear_count + len(transvections)))),
            generator,
            name=f"Source-defined Eichler approximate generators in O({lattice})",
        )

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
        """
        lattice = model.lattice()
        ring = lattice.base_ring()
        square = ring(square)
        complement = model.orthogonal_complement()
        discriminant = complement.discriminant_group()
        correlation = complement.correlation_morphism()
        inclusion = model.complement_inclusion()
        e, f, _e_prime, _f_prime = model.hyperbolic_basis()

        def representative(discriminant_class: DiscriminantQuadraticModules.ElementMethods) -> Lattices.ElementMethods:
            discriminant_class = discriminant(discriminant_class)
            order = ring(int(discriminant_class.additive_order()))
            dual_lift = discriminant.dual_lattice_lift(discriminant_class)
            k = correlation.lift(dual_lift.parent().scalar_multiple(order, dual_lift))
            numerator = square - complement.q(k)
            denominator = ring(2) * order * order
            coefficient = numerator // denominator
            if denominator * coefficient != numerator:
                raise ArithmeticError("a covering discriminant class must give an integral hyperbolic coefficient")
            return lattice.scalar_multiple(order, e) + lattice.scalar_multiple(order * coefficient, f) + inclusion(k)

        classes = complement.covering_discriminant_classes(square)
        # IndefiniteForms.g:181-224, SetListClassesOrbitwise(GRPeasy).
        # The canonical inclusion A_K -> A_(2U+K) is an isomorphism, since
        # 2U is unimodular. GetEasyIsometries can permute a hyperbolic block
        # of K with either distinguished U, so its isometries need not
        # preserve K itself. Transfer their action through the *full*
        # discriminant group, rather than restricting their lattice maps to K.
        easy = OrbitCoverModel.source_easy_generators(model)
        discriminant_inclusion = inclusion.discriminant_inclusion()
        ambient_discriminant = lattice.discriminant_group()
        image = ambient_discriminant.orthogonal_group().subgroup(tuple(isometry.discriminant_morphism() for isometry in easy))
        embedded_classes = finite_ordered_set(tuple(discriminant_inclusion(point) for point in classes))
        action = FiniteGSets(image)(embedded_classes, lambda automorphism, point: automorphism(point))
        selected = finite_ordered_set(tuple(discriminant_inclusion.lift(orbit.representative()) for orbit in action.orbits()))
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
        if int(representatives.cardinality()) == 0:
            raise ValueError("full-group assembly requires a nonempty primitive covering family")
        base_label = next(iter(representatives.index_set()))
        base_vector = representatives.value(base_label)
        witnesses: dict[DiscriminantQuadraticModules.ElementMethods, LatticeIsometryMethods] = {}
        for label in representatives.index_set():
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
        return TwoUFullGeneratingDatum(model, model.lattice().base_ring()(square), base_vector, approximate, stabilizer, orbit_representatives, orbit_witnesses, family)

    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        self.lattice = lattice
        if lattice.is_even() and lattice.splits_two_hyperbolic_planes():
            self.envelope = None
            self.model = lattice.two_u_eichler_model_from_represented_biproduct()
        else:
            self.envelope = EichlerEnvelope(lattice)
            self.model = self.envelope.model

    def approximate_generators(self) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of A_L: the Eichler subgroup of the 2U model, preserving L when L is reached through the envelope."""
        if self.envelope is None:
            return tuple(self.model.approximate_generating_family())
        envelope = self.envelope
        stabilizer = envelope.stabilizer_of_lattice().stabilizer()
        inverse = envelope.rational_isometry.inverse()
        return tuple(self.lattice.integral_isometry(inverse * generator * envelope.rational_isometry) for generator in stabilizer.group_generators())

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

        envelope = self.envelope
        envelope_square = envelope.scale**2 * square
        if envelope_square % 2 == 0:
            e, f, _e_prime, _f_prime = self.model.hyperbolic_basis()
            candidate = e + self.model.lattice().scalar_multiple(self.model.lattice().base_ring()(envelope_square // 2), f)
            # ApproximateModels.h:1372-1400 tests this representative against
            # every coset of the integral-structure stabilizer.  Its identity
            # coset is necessarily in that transversal.  Test the exact
            # identity-coset equation before constructing the full finite
            # action; failure still follows the source's exhaustive search.
            if envelope.similarity.is_in_image(candidate):
                witness = envelope.similarity.lift(candidate)
                if int(witness.q()) != square:
                    raise ArithmeticError("the identity-coset lift has the wrong norm")
                return witness
            coset_space = envelope.stabilizer_of_lattice().right_cosets()
            section = coset_space.right_coset_section()
            for coset in coset_space:
                lifted = self._pullback_envelope_vector(candidate, section(coset))
                if lifted is not None:
                    return lifted
        if square == 0:
            covers = self.primitive_covering_vectors(0)
            return covers[0] if covers else None
        for content in square_divisors(square):
            covers = self.primitive_covering_vectors(square // (content * content))
            if covers:
                return self.lattice.scalar_multiple(content, covers[0])
        return None

    def primitive_covering_vectors(self, square: int) -> tuple[Lattices.ElementMethods, ...]:
        r"""Primitive vectors of the square meeting every O(L)-orbit of primitive vectors of that square."""
        if self.envelope is None:
            return tuple(self.model.covering_vector_representatives(square))
        envelope = self.envelope
        # sigma(v) = c g p with g in A_E and p a covering representative; writing g = h t with
        # h in the stabilizer of phi(L) and t a right-coset representative, v is A_L-equivalent
        # to sigma^-1(c t p). So t ranges over a transversal of the right cosets G_L t.
        coset_space = envelope.stabilizer_of_lattice().right_cosets()
        section = coset_space.right_coset_section()
        cosets = tuple(section(point) for point in coset_space)
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
        for content in contents:
            for primitive in self.model.covering_vector_representatives(envelope_square // content**2):
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
        r"""The first represented positive norm, as in ``GetFirstNorm`` (IndefiniteForms.g:1131-1142)."""
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
