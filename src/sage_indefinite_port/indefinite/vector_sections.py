r"""Orthogonal sections of lattice vectors and their rational lifts."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass

from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeEmbeddingMethods,
    LatticeIsometryMethods,
)
from dzack_research.preamble.categories.lattices import IsotropicReductions, Lattices
from dzack_research.preamble.categories.rings.ring_foundation import OwnedRings

from sage_indefinite_port.indefinite.isotropic_lifts import (
    CodimensionOneIsotropicExtension,
    IsometryExtensionTorsor,
)


@dataclass(frozen=True)
class NonIsotropicVectorSection:
    r"""Orthogonal section of a nonisotropic vector."""

    vector: Lattices.ElementMethods
    perpendicular: Lattices.ParentMethods
    inclusion: LatticeEmbeddingMethods
    reduction: Lattices.ParentMethods

    def reduced_object(self) -> Lattices.ParentMethods:
        r"""Return the reduced lattice \(v^\perp\)."""
        return self.reduction

    def rational_lift(
        self,
        reduced_isometry: LatticeIsometryMethods,
        *,
        target: NonIsotropicVectorSection | None = None,
    ) -> IsometryExtensionTorsor:
        r"""Extend an isometry of perpendicular lattices uniquely over the fraction field."""
        match target:
            case None:
                target_section = self
            case _:
                target_section = target

        match (
            reduced_isometry.domain() is self.perpendicular,
            reduced_isometry.codomain() is target_section.perpendicular,
        ):
            case (True, True):
                pass
            case _:
                raise ValueError("the reduced isometry has the wrong endpoints")

        match self.vector.q() == target_section.vector.q():
            case True:
                pass
            case False:
                raise ValueError("source and target vectors must have the same norm")

        source_ambient = self.inclusion.codomain()
        target_ambient = target_section.inclusion.codomain()
        source_ring = source_ambient.base_ring()
        match source_ring is target_ambient.base_ring():
            case True:
                pass
            case False:
                raise ValueError("source and target lattices must have the same base ring")

        match source_ring:
            case OwnedRings.NoZeroDivisors.Commutative.ParentMethods() as domain:
                fraction_map = domain.fraction_field_map()
            case _:
                raise TypeError("rational lifting requires a lattice over an integral domain")
        source_rational: Lattices.ParentMethods = source_ambient.base_change(fraction_map)
        target_rational: Lattices.ParentMethods = target_ambient.base_change(fraction_map)

        def extend_source(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            result: Lattices.ElementMethods = source_rational.linear_combination(
                {
                    label: fraction_map(coordinates(label))
                    for label in source_ambient.module_generating_set()
                    if coordinates(label)
                }
            )
            return result

        def extend_target(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            result: Lattices.ElementMethods = target_rational.linear_combination(
                {
                    label: fraction_map(coordinates(label))
                    for label in target_ambient.module_generating_set()
                    if coordinates(label)
                }
            )
            return result

        source_perpendicular_inclusion = self.inclusion.base_change(fraction_map)
        target_perpendicular_inclusion = target_section.inclusion.base_change(fraction_map)
        reduced_rational = reduced_isometry.base_change(fraction_map)
        source_vector = extend_source(self.vector)
        target_vector = extend_target(target_section.vector)
        source_norm = source_rational.q(source_vector)

        def image(label: Hashable) -> Lattices.ElementMethods:
            source_generator = source_rational.module_generator(label)
            vector_coefficient = source_rational.b(source_generator, source_vector) / source_norm
            perpendicular_part = source_generator - source_rational.scalar_multiple(
                vector_coefficient,
                source_vector,
            )
            reduced_part = source_perpendicular_inclusion.lift(perpendicular_part)
            target_perpendicular_part = target_perpendicular_inclusion(reduced_rational(reduced_part))
            result: Lattices.ElementMethods = (
                target_perpendicular_part
                + target_rational.scalar_multiple(
                    vector_coefficient,
                    target_vector,
                )
            )
            return result

        return IsometryExtensionTorsor(source_rational.Isom(target_rational)(image), ())


@dataclass(frozen=True)
class IsotropicVectorSection:
    r"""Orthogonal section of a primitive isotropic vector."""

    vector: Lattices.ElementMethods
    perpendicular: Lattices.ParentMethods
    inclusion: LatticeEmbeddingMethods
    reduction: IsotropicReductions.ParentMethods

    def reduced_object(self) -> IsotropicReductions.ParentMethods:
        r"""Return the isotropic reduction \(v^\perp/Rv\)."""
        return self.reduction

    def rational_lift(
        self,
        reduced_isometry: LatticeIsometryMethods,
        *,
        target: IsotropicVectorSection | None = None,
    ) -> IsometryExtensionTorsor:
        r"""Lift a reduction isometry through \(v^\perp\) and over the fraction field."""
        match target:
            case None:
                target_section = self
            case _:
                target_section = target

        match (
            reduced_isometry.domain() is self.reduction,
            reduced_isometry.codomain() is target_section.reduction,
        ):
            case (True, True):
                pass
            case _:
                raise ValueError("the reduced isometry has the wrong endpoints")

        source_ambient = self.inclusion.codomain()
        target_ambient = target_section.inclusion.codomain()
        match source_ambient.base_ring() is target_ambient.base_ring():
            case True:
                pass
            case False:
                raise ValueError("source and target lattices must have the same base ring")

        source_reduction = self.reduction
        target_reduction = target_section.reduction
        source_projection = source_reduction.projection()
        source_line_in_perpendicular = source_reduction.isotropic_inclusion()
        target_line_in_perpendicular = target_reduction.isotropic_inclusion()
        source_line = source_reduction.isotropic_sublattice()
        target_line = target_reduction.isotropic_sublattice()
        source_line_vector = source_reduction.isotropic_embedding().lift(self.vector)
        target_line_vector = target_reduction.isotropic_embedding().lift(target_section.vector)
        (source_line_label,) = tuple(source_line.module_generating_set())
        source_vector_coordinate = source_line_vector.to_vector()(source_line_label)

        def reduction_lift(
            reduction: IsotropicReductions.ParentMethods,
            element: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            perpendicular = reduction.orthogonal_complement()
            coordinates = element.to_vector()
            lifts = reduction.coordinate_frame()
            result: Lattices.ElementMethods = sum(
                (
                    perpendicular.scalar_multiple(coordinates(label), lifts(label))
                    for label in reduction.module_generating_set()
                ),
                perpendicular.zero(),
            )
            return result

        def partial_image(label: Hashable) -> Lattices.ElementMethods:
            source_element = self.perpendicular.module_generator(label)
            quotient_element = source_projection(source_element)
            chosen_lift = reduction_lift(source_reduction, quotient_element)
            source_line_part = source_element - chosen_lift
            line_preimage = source_line_in_perpendicular.lift(source_line_part)
            line_coordinate = line_preimage.to_vector()(source_line_label)
            line_coefficient = line_coordinate / source_vector_coordinate
            target_line_element = target_line.scalar_multiple(
                line_coefficient,
                target_line_vector,
            )
            target_line_part = target_line_in_perpendicular(target_line_element)
            target_quotient_element: Lattices.ElementMethods = reduced_isometry(quotient_element)
            result: Lattices.ElementMethods = target_line_part + reduction_lift(
                target_reduction,
                target_quotient_element,
            )
            return result

        partial = self.perpendicular.Isom(target_section.perpendicular)(partial_image)

        source_ring = source_ambient.base_ring()
        match source_ring:
            case OwnedRings.NoZeroDivisors.Commutative.ParentMethods() as domain:
                fraction_map = domain.fraction_field_map()
            case _:
                raise TypeError("rational lifting requires a lattice over an integral domain")
        source_rational: Lattices.ParentMethods = source_ambient.base_change(fraction_map)
        target_rational: Lattices.ParentMethods = target_ambient.base_change(fraction_map)

        def extend_source(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            result: Lattices.ElementMethods = source_rational.linear_combination(
                {
                    label: fraction_map(coordinates(label))
                    for label in source_ambient.module_generating_set()
                    if coordinates(label)
                }
            )
            return result

        def extend_target(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            result: Lattices.ElementMethods = target_rational.linear_combination(
                {
                    label: fraction_map(coordinates(label))
                    for label in target_ambient.module_generating_set()
                    if coordinates(label)
                }
            )
            return result

        source_perpendicular_rational = source_rational.subobject_on(tuple(extend_source(self.inclusion(generator)) for generator in self.perpendicular.module_generators()))
        target_perpendicular_rational = target_rational.subobject_on(
            tuple(extend_target(target_section.inclusion(generator)) for generator in target_section.perpendicular.module_generators())
        )
        source_extended: Lattices.ParentMethods = self.perpendicular.base_change(fraction_map)
        target_extended: Lattices.ParentMethods = target_section.perpendicular.base_change(fraction_map)
        source_to_subspace = source_extended.Isom(source_perpendicular_rational)(
            lambda label: source_perpendicular_rational.inclusion().lift(extend_source(self.inclusion(self.perpendicular.module_generator(label))))
        )
        target_to_subspace = target_extended.Isom(target_perpendicular_rational)(
            lambda label: target_perpendicular_rational.inclusion().lift(extend_target(target_section.inclusion(target_section.perpendicular.module_generator(label))))
        )
        partial_extended = partial.base_change(fraction_map)
        subspace_to_source = source_to_subspace.inverse()
        rational_partial = source_perpendicular_rational.Isom(target_perpendicular_rational)(
            lambda label: target_to_subspace(partial_extended(subspace_to_source(source_perpendicular_rational.module_generator(label))))
        )
        extension = CodimensionOneIsotropicExtension(
            source_perpendicular_rational,
            target_perpendicular_rational,
            rational_partial,
        ).rational_extension()
        return IsometryExtensionTorsor(extension.extension, ())


type VectorOrthogonalSection = NonIsotropicVectorSection | IsotropicVectorSection


def orthogonal_section(
    vector: Lattices.ElementMethods,
) -> VectorOrthogonalSection:
    r"""Return the orthogonal section of a nonzero vector."""
    match vector.is_isotropic():
        case False:
            perpendicular = vector.orthogonal_complement()
            return NonIsotropicVectorSection(
                vector=vector,
                perpendicular=perpendicular,
                inclusion=perpendicular.inclusion(),
                reduction=perpendicular,
            )
        case True:
            reduction = vector.isotropic_reduction()
            perpendicular = reduction.orthogonal_complement()
            return IsotropicVectorSection(
                vector=vector,
                perpendicular=perpendicular,
                inclusion=perpendicular.inclusion(),
                reduction=reduction,
            )
        case _:
            raise TypeError("is_isotropic() must return a Boolean value")
