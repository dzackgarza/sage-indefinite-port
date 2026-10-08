r"""Orthogonal sections of lattice vectors and their rational lifts."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass
from functools import cached_property

from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeEmbeddingMethods,
    LatticeIsometryMethods,
)
from dzack_research.preamble.categories.lattices import IsotropicReductions, Lattices
from dzack_research.preamble.categories.rings.ring_foundation import (
    OwnedRings,
    RingMorphism,
)

from sage_indefinite_port.indefinite.isotropic_lifts import (
    CodimensionOneIsotropicExtension,
    IsometryExtensionTorsor,
)

type IsotropicLiftContext = tuple[
    RingMorphism,
    Lattices.ParentMethods,
    Lattices.ParentMethods,
    Lattices.ParentMethods,
    LatticeIsometryMethods,
    LatticeIsometryMethods,
]
"""Fraction map, rational ambient, rational perpendicular subspace, base-changed perpendicular, and the isometry to that subspace with its inverse."""


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
        r"""Extend the perpendicular isometry over the fraction field."""
        target_section = self if target is None else target
        if reduced_isometry.domain() is not self.perpendicular or reduced_isometry.codomain() is not target_section.perpendicular:
            raise ValueError("the reduced isometry has the wrong endpoints")
        extension = reduced_isometry.extend_from_perpendicular(self.vector, target_section.vector)
        return IsometryExtensionTorsor(extension, ())

    def integral_lift(
        self,
        reduced_isometry: LatticeIsometryMethods,
        *,
        target: NonIsotropicVectorSection | None = None,
    ) -> LatticeIsometryMethods | None:
        r"""Return the extension when it preserves both integral lattices."""
        target_section = self if target is None else target
        extension = self.rational_lift(reduced_isometry, target=target_section).rational_isometry()
        return extension.integral_restriction(self.inclusion.codomain(), target_section.inclusion.codomain())


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

    @cached_property
    def _rational_lift_context(self) -> IsotropicLiftContext:
        source_ambient = self.inclusion.codomain()
        source_ring = source_ambient.base_ring()
        match source_ring:
            case OwnedRings.NoZeroDivisors.Commutative.ParentMethods() as domain:
                fraction_map = domain.fraction_field_map()
            case _:
                raise TypeError("rational lifting requires a lattice over an integral domain")
        source_rational: Lattices.ParentMethods = source_ambient.base_change(fraction_map)

        def extend_source(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            return source_rational.linear_combination({label: fraction_map(coordinates(label)) for label in source_ambient.module_generating_set() if coordinates(label)})

        source_perpendicular_rational = source_rational.subobject_on(tuple(extend_source(self.inclusion(generator)) for generator in self.perpendicular.module_generators()))
        source_extended: Lattices.ParentMethods = self.perpendicular.base_change(fraction_map)
        source_to_subspace = source_extended.Isom(source_perpendicular_rational)(
            lambda label: source_perpendicular_rational.inclusion().lift(extend_source(self.inclusion(self.perpendicular.module_generator(label))))
        )
        return (
            fraction_map,
            source_rational,
            source_perpendicular_rational,
            source_extended,
            source_to_subspace,
            source_to_subspace.inverse(),
        )

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
                (perpendicular.scalar_multiple(coordinates(label), lifts(label)) for label in reduction.module_generating_set()),
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

        partial_images = tuple(partial_image(label) for label in self.perpendicular.module_generating_set())
        partial = self.perpendicular.Isom(target_section.perpendicular)(partial_images)

        (
            fraction_map,
            _source_rational,
            source_perpendicular_rational,
            _source_extended,
            _source_to_subspace,
            subspace_to_source,
        ) = self._rational_lift_context
        (
            _target_fraction_map,
            _target_rational,
            target_perpendicular_rational,
            _target_extended,
            target_to_subspace,
            _subspace_to_target,
        ) = target_section._rational_lift_context
        partial_extended = partial.base_change(fraction_map)
        rational_partial_images = tuple(
            target_to_subspace(partial_extended(subspace_to_source(source_perpendicular_rational.module_generator(label))))
            for label in source_perpendicular_rational.module_generating_set()
        )
        rational_partial = source_perpendicular_rational.Isom(target_perpendicular_rational)(rational_partial_images)
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
