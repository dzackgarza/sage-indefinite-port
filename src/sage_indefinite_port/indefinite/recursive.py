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
from dzack_research.preamble.categories.rational_integral_stabilizers import IntegralStructureAction

from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from sage_indefinite_port.indefinite.lorentzian_cells import lorentzian_isometry, lorentzian_orthogonal_group_generators
from sage_indefinite_port.invariants import AttackProfile, LatticePrefilter, VectorPrefilter

_UNTWIST = Lattices(ZZ).twist_functor(-1)


def _untwisted(isometry: LatticeIsometryMethods, profile: AttackProfile) -> LatticeIsometryMethods:
    r"""An isometry of the signed view, as an isometry of the lattice itself."""
    return isometry if profile.sign == 1 else _UNTWIST(isometry)


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
        reduction = lattice.indefinite_reduction()
        to_lattice = reduction.isometry
        return tuple(to_lattice * generator * to_lattice.inverse() for generator in self._reduced_orthogonal_group_generators(reduction.reduced))

    def _reduced_orthogonal_group_generators(self, lattice: Lattices.ParentMethods) -> tuple[LatticeIsometryMethods, ...]:
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
        source_reduction = source.indefinite_reduction()
        target_reduction = target.indefinite_reduction()
        reduced = self._reduced_isometry(source_reduction.reduced, target_reduction.reduced)
        return None if reduced is None else target_reduction.isometry * reduced * source_reduction.isometry.inverse()

    def _reduced_isometry(self, source: Lattices.ParentMethods, target: Lattices.ParentMethods) -> LatticeIsometryMethods | None:
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
        r"""For q(v) != 0: the group generated by the extensions across v^perp of generators of O(v^perp), acting on the integral structures of L tensor Q, at L."""
        lattice = vector.parent()
        complement = vector.orthogonal_complement().inclusion().domain()
        extensions = tuple(generator.extension_across(vector, vector) for generator in self.orthogonal_group_generators(complement))
        rational = lattice.base_change(lattice.base_ring().fraction_field_map())
        return IntegralStructureAction(rational.O().subgroup(extensions), rational.integral_structure(tuple(lattice.module_generators())))

    def pointed_automorphism_generators(self, vector: Lattices.ElementMethods) -> tuple[LatticeIsometryMethods, ...]:
        r"""Generators of Aut((L, v)): for q(v) = 0 the pointwise stabilizer of Zv in the parabolic of Zv (design: "the parabolic below")."""
        lattice = vector.parent()
        if vector.q() == lattice.base_ring().zero():
            return tuple(lattice.primitive_isotropic_subobject(vector).isotropic_reduction().pointwise_parabolic_subgroup().group_generators())
        return tuple(lattice.integral_isometry(generator) for generator in self._lattice_action(vector).stabilizer().group_generators())

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
        if rational_transport.is_integral_on(source_lattice):
            return target_lattice.integral_isometry(rational_transport)
        image = rational_transport.codomain().integral_structure(tuple(rational_transport(generator) for generator in source_lattice.module_generators()))
        correction = self._lattice_action(target).transporter(image)
        return None if correction is None else target_lattice.integral_isometry(correction * rational_transport)

    def _isotropic_pointed_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""An element of Isom((L, v), (M, w)) for isotropic v, w.

        Such an isometry induces psi: K -> K' on the reductions. With psi_0 one isometry K -> K'
        and H' the image of Aut((M, Zw)) in O(K'), the psi that lift integrally form a union of
        cosets H' c psi_0, so it suffices to look for an integral member among the lifts of
        +-c psi_0 for c in H' \ O(K') (the sign accounts for H' moving w to -w).
        """
        source_reduction = source.parent().primitive_isotropic_subobject(source).isotropic_reduction()
        target_reduction = target.parent().primitive_isotropic_subobject(target).isotropic_reduction()
        target_quotient = target_reduction.quotient_lattice()
        reduced = self.isometry(source_reduction.quotient_lattice(), target_quotient)
        if reduced is None:
            return None
        negation = target_quotient.Aut()({label: -target_quotient.module_generator(label) for label in target_quotient.module_generating_set()})
        group = target_quotient.O().subgroup(self.orthogonal_group_generators(target_quotient))
        for coset in group.right_cosets(target_reduction.levi_action().image()):
            for candidate in (coset * reduced, negation * coset * reduced):
                members = source_reduction.rational_lifts(target_reduction, candidate).integral_members()
                if not members.is_empty():
                    return target.parent().integral_isometry(members.an_element())
        return None

    def _split_isotropic_isometry(self, source: Lattices.ElementMethods, target: Lattices.ElementMethods) -> LatticeIsometryMethods | None:
        r"""For e of divisibility one in an even lattice, L = U + e^perp/e; an isometry is id_U plus one of the reductions."""
        source_splitting = source.parent().primitive_isotropic_subobject(source).isotropic_reduction().integral_hyperbolic_splitting()
        target_splitting = target.parent().primitive_isotropic_subobject(target).isotropic_reduction().integral_hyperbolic_splitting()
        hyperbolic_plane, source_complement = source_splitting.codomain().biproduct_factors()
        _hyperbolic_plane, target_complement = target_splitting.codomain().biproduct_factors()
        reduced = self.isometry(source_complement.inclusion().domain(), target_complement.inclusion().domain())
        if reduced is None:
            return None
        split = hyperbolic_plane.inclusion().domain().identity_morphism().biproduct_map(reduced, source=source_splitting.codomain(), target=target_splitting.codomain())
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
            VectorClass(lattice.linear_combination(kept.monomial_coefficients()), tuple(_untwisted(isometry, profile) for isometry in witnesses))
            for bucket in buckets.values()
            for kept, witnesses in bucket
        )


__all__ = ["EichlerRecursion", "VectorClass"]
