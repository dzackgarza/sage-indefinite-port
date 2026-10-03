r"""Orthogonal sections of lattice vectors and their rational lifts."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass

from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeEmbeddingMethods,
    LatticeIsometryMethods,
)
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.categories.modules.framed.framed_free_modules import (
    FramedFreeModules,
)

from sage_indefinite_port.indefinite.isotropic_lifts import IsometryExtensionTorsor


@dataclass(frozen=True)
class NonIsotropicVectorSection:
    r"""Orthogonal section of a nonisotropic vector."""

    vector: FramedFreeModules.ElementMethods
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
        match source_ambient.base_ring() is target_ambient.base_ring():
            case True:
                pass
            case False:
                raise ValueError("source and target lattices must have the same base ring")

        fraction_map = source_ambient.base_ring().fraction_field_map()
        source_rational = source_ambient.base_change(fraction_map)
        target_rational = target_ambient.base_change(fraction_map)

        def extend_source(
            vector: FramedFreeModules.ElementMethods,
        ) -> FramedFreeModules.ElementMethods:
            coordinates = vector.to_vector()
            return source_rational.linear_combination({label: fraction_map(coordinates(label)) for label in source_ambient.module_generating_set() if coordinates(label)})

        def extend_target(
            vector: FramedFreeModules.ElementMethods,
        ) -> FramedFreeModules.ElementMethods:
            coordinates = vector.to_vector()
            return target_rational.linear_combination({label: fraction_map(coordinates(label)) for label in target_ambient.module_generating_set() if coordinates(label)})

        source_perpendicular_inclusion = self.inclusion.base_change(fraction_map)
        target_perpendicular_inclusion = target_section.inclusion.base_change(fraction_map)
        reduced_rational = reduced_isometry.base_change(fraction_map)
        source_vector = extend_source(self.vector)
        target_vector = extend_target(target_section.vector)
        source_norm = source_rational.q(source_vector)

        def image(label: Hashable) -> FramedFreeModules.ElementMethods:
            source_generator = source_rational.module_generator(label)
            vector_coefficient = source_rational.b(source_generator, source_vector) / source_norm
            perpendicular_part = source_generator - source_rational.scalar_multiple(
                vector_coefficient,
                source_vector,
            )
            reduced_part = source_perpendicular_inclusion.lift(perpendicular_part)
            target_perpendicular_part = target_perpendicular_inclusion(reduced_rational(reduced_part))
            return target_perpendicular_part + target_rational.scalar_multiple(
                vector_coefficient,
                target_vector,
            )

        return IsometryExtensionTorsor(source_rational.Isom(target_rational)(image), ())


def orthogonal_section(
    vector: FramedFreeModules.ElementMethods,
) -> NonIsotropicVectorSection:
    r"""Return the orthogonal section of a nonisotropic vector."""
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
            raise ValueError("the isotropic branch uses IsotropicReduction")
