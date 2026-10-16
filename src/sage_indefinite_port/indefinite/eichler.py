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

from math import isqrt

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.categories.eichler_criterion import TwoUEichlerModel
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.framed.formed.discriminant_modules import DiscriminantQuadraticModules
from dzack_research.preamble.categories.rational_integral_stabilizers import IntegralStructureAction
from dzack_research.preamble.categories.sets.finite_ordered_sets import finite_ordered_set
from dzack_research.preamble.categories.sets.indexed_families import IndexedFamily, finite_indexed_family
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector as vector_space_element
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.element import RingElement
from sage.structure.sage_object import SageObject


def square_divisors(integer: int) -> tuple[int, ...]:
    r"""The positive c with c^2 dividing the nonzero integer."""
    value = abs(int(integer))
    if value == 0:
        raise ValueError("zero has infinitely many square divisors: covers of norm 0 are by primitive isotropic vectors only")
    return tuple(divisor for divisor in range(1, isqrt(value) + 1) if value % (divisor * divisor) == 0)


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
    r"""The similarity sigma = N phi: L -> E = U + U + K of Dutour Sikirić's hyperplane-basis step."""

    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        first = lattice.isotropic_vector()
        first_partner = first.hyperbolic_partner()
        first_plane = lattice.subobject_on([first, first_partner])
        complement = first_plane.orthogonal_complement()
        complement_lattice = complement.inclusion().domain()
        second = complement_lattice.isotropic_vector()
        second_partner = second.hyperbolic_partner()
        embed = complement.inclusion()
        second, second_partner = embed(second), embed(second_partner)
        planes = lattice.subobject_on([first, first_partner, second, second_partner])
        kernel = planes.orthogonal_complement()
        kernel_lattice = kernel.inclusion().domain()
        hyperbolic_plane = Lattices(lattice.base_ring())("U")
        self.lattice = lattice
        self.envelope = hyperbolic_plane + hyperbolic_plane + kernel_lattice
        self.model = self.envelope.two_u_eichler_model_from_represented_biproduct()
        e_1, f_1, e_2, f_2 = self.model.hyperbolic_basis()
        complement_embedding = self.model.complement_inclusion()
        first_pairing = lattice.b(first, first_partner)
        second_pairing = lattice.b(second, second_partner)
        source = (first, first_partner, second, second_partner) + tuple(kernel.inclusion()(vector) for vector in kernel_lattice.module_generators())
        target = (first_pairing * e_1, f_1, second_pairing * e_2, f_2) + tuple(complement_embedding(vector) for vector in kernel_lattice.module_generators())
        fraction_map = lattice.base_ring().fraction_field_map()
        spanned = lattice.subobject_on(list(source))
        # The span basis is echelonized; each basis vector is an integral combination of the
        # independent family source, and its image is the same combination of target.
        labels = tuple(lattice.module_generating_set())
        source_rows = matrix(SageQQ, [[int(vector.to_vector()(label)) for label in labels] for vector in source])
        span_inclusion = spanned.inclusion()
        span_images = []
        for generator in span_inclusion.domain().module_generators():
            image = span_inclusion(generator)
            coefficients = source_rows.solve_left(vector_space_element([int(image.to_vector()(label)) for label in labels]))
            if any(entry.denominator() != 1 for entry in coefficients):
                raise ArithmeticError("a span basis vector must be an integral combination of the spanning family")
            span_images.append(sum((lattice.base_ring()(int(entry)) * image_vector for entry, image_vector in zip(coefficients, target, strict=True)), self.envelope.zero()))
        span_to_envelope = span_inclusion.domain().hom(span_images, codomain=self.envelope).base_change(fraction_map)
        rational_inclusion = spanned.inclusion().base_change(fraction_map)
        rational_lattice = lattice.base_change(fraction_map)
        self.rational_isometry = rational_lattice.Isom(self.envelope.base_change(fraction_map))(
            lambda label: span_to_envelope(rational_inclusion.lift(rational_lattice.module_generator(label)))
        )
        self.similarity = self.rational_isometry.integral_similarity(lattice, self.envelope)
        self.scale = self.similarity.scale()
        self.fraction_map = fraction_map

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
    def source_generating_family(model: TwoUEichlerModel):
        r"""Return the represented source-defined finite family inside ``O(2U+K)``.

        The family retains five mathematically different sources of
        isometries: the two determinant-model ``SL_2(ZZ)`` actions, the
        chosen generators of ``O(K)`` extended by the identity on ``2U``,
        the selected ``K``-direction Eichler transvections, and explicit
        lifts of generators of ``O(A_L)``.  Labels retain those sources even
        when two labels happen to determine the same lattice isometry.

        This is generating *data*, not a claim that the generated subgroup is
        the full orthogonal group.  That equality belongs to the recursive
        generation theorem and its parabolic-stabilizer hypotheses.
        """
        from dzack_research.preamble.categories.sets.finite_ordered_sets import (
            finite_ordered_set,
        )
        from dzack_research.preamble.categories.sets.indexed_families import (
            finite_indexed_family,
        )

        special_linear_generators = model.special_linear_group().group_generators()
        complement_generators = model.orthogonal_complement().O().group_generators()
        eichler = model.complement_eichler_transvections()
        discriminant_lifts = model.lattice().O().discriminant_orthogonal_group_generator_lifts()
        labels = finite_ordered_set(
            tuple(("left-SL2", generator) for generator in special_linear_generators)
            + tuple(("right-SL2", generator) for generator in special_linear_generators)
            + tuple(("O(K)", generator) for generator in complement_generators)
            + tuple(("Eichler", label) for label in eichler.index_set())
            + tuple(("discriminant-lift", label) for label in discriminant_lifts.index_set())
        )

        def generator(label):
            kind, datum = label
            match kind:
                case "left-SL2":
                    return model.left_action(datum)
                case "right-SL2":
                    return model.right_action(datum)
                case "O(K)":
                    return model.complement_action(datum)
                case "Eichler":
                    return eichler[datum]
                case "discriminant-lift":
                    return discriminant_lifts[datum]
                case _:
                    raise ValueError(f"{label!r} does not name a generator of O(2U + K) for {model}: its kind {kind!r} is not a known kind of generator")

        return finite_indexed_family(
            labels,
            generator,
            name=f"Source-defined Eichler generating family in O({model.lattice()})",
        )

    @staticmethod
    def split_generating_family(model: TwoUEichlerModel) -> IndexedFamily[int, LatticeIsometryMethods]:
        r"""The two SL2 actions and additive Eichler generators in a represented 2U model.

        For each of e, f, e', f', a basis of its perpendicular generates all
        its transvections by E(u,x) E(u,y) = E(u,x+y). Together with the left
        and right SL2 actions these give the finite approximate subgroup.
        """
        isotropic = model.hyperbolic_basis()
        e, f, e_prime, f_prime = isotropic
        complement = model.complement_inclusion()
        embedded = tuple(complement(point) for point in model.orthogonal_complement().module_generators())
        perpendicular = ((e, e_prime, f_prime), (f, e_prime, f_prime), (e, f, e_prime), (e, f, f_prime))
        transvections = tuple((point, vector) for point, basis in zip(isotropic, perpendicular, strict=True) for vector in (*basis, *embedded))
        linear = tuple(model.special_linear_group().group_generators())
        count = len(linear)

        def generator(index: int) -> LatticeIsometryMethods:
            if index < count:
                return model.left_action(linear[index])
            if index < 2 * count:
                return model.right_action(linear[index - count])
            point, vector = transvections[index - 2 * count]
            return model.eichler_transvection(point, vector)

        return finite_indexed_family(
            finite_ordered_set(tuple(range(2 * count + len(transvections)))), generator, name=f"Eichler approximate generators in O({model.lattice()})"
        )

    @staticmethod
    def split_covering_vectors(model: TwoUEichlerModel, square: RingElement) -> IndexedFamily[DiscriminantQuadraticModules.ElementMethods, Lattices.ElementMethods]:
        r"""One primitive vector per prescribed-invariant class in the 2U model.

        If x in A_K has order d and y in K dual lifts x, put k=d*y and
        m=(square-k^2)/(2d^2). Then d*e+d*m*f+k has the prescribed square
        and divided class. Its pairing with f is d, while the class order
        forces primitivity and divisibility exactly d.
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

        return finite_indexed_family(complement.covering_discriminant_classes(square), representative, name=f"Eichler covering vectors of square {square} in {lattice}")

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
        if lattice.splits_two_hyperbolic_planes():
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
        exponent = int(envelope.similarity.cokernel().exponent())
        contents = square_divisors(envelope_square) if envelope_square else tuple(divisor for divisor in range(1, exponent + 1) if exponent % divisor == 0)
        envelope_lattice = self.model.lattice()
        rational_envelope = envelope.rational_isometry.codomain()
        rational_ring = rational_envelope.base_ring()
        integers = envelope_lattice.base_ring()
        labels = tuple(envelope_lattice.module_generating_set())
        for content in contents:
            for primitive in self.model.covering_vector_representatives(envelope_square // content**2):
                point = self.model.lattice().scalar_multiple(content, primitive)
                rational_point = rational_envelope.linear_combination({label: rational_ring(int(point.to_vector()(label))) for label in labels if point.to_vector()(label)})
                for representative in cosets:
                    rational_moved = representative(rational_point).to_vector()
                    # Only an integral t p lies in E, and only a vector of E can be sigma(v).
                    if any(rational_moved(label).denominator() != 1 for label in labels):
                        continue
                    moved = envelope_lattice.linear_combination({label: integers(int(rational_moved(label))) for label in labels if rational_moved(label)})
                    if envelope.similarity.is_in_image(moved):
                        vector = envelope.similarity.lift(moved)
                        if vector.is_primitive():
                            found.append(vector)
        return tuple(found)

    def covering_vectors(self, square: int, content: int) -> tuple[Lattices.ElementMethods, ...]:
        r"""Vectors of the square and content meeting every O(L)-orbit of such vectors: c times the primitive covers of norm square / c^2."""
        if square % (content * content):
            raise ValueError(f"no vector of content {content} has square {square}")
        return tuple(self.lattice.scalar_multiple(content, vector) for vector in self.primitive_covering_vectors(square // (content * content)))

    def splitting_vector(self) -> Lattices.ElementMethods:
        r"""A primitive vector of positive norm: e + f of the first hyperbolic plane, of norm 2, through the envelope when needed."""
        e_1, f_1, _e_2, _f_2 = self.model.hyperbolic_basis()
        if self.envelope is None:
            return e_1 + f_1
        return self.lattice.positive_vector()


__all__ = [
    "EichlerEnvelope",
    "OrbitCoverModel",
    "TwoUFullGeneratingDatum",
    "square_divisors",
]
