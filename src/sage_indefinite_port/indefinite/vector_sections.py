r"""Orthogonal sections of lattice vectors and their rational lifts."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass
from functools import cached_property

from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeEmbeddingMethods,
    LatticeIsometryMethods,
    _module_matrix,
)
from dzack_research.preamble.categories.lattices import IsotropicReductions, Lattices
from dzack_research.preamble.categories.modules.pure.modules import _engine_matrix
from dzack_research.preamble.categories.rings.ring_foundation import (
    OwnedRings,
    _engine_element,
)
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ as SageQQ

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

    @cached_property
    def _rational_lift_context(self):
        source_ambient = self.inclusion.codomain()
        source_ring = source_ambient.base_ring()
        match source_ring:
            case OwnedRings.NoZeroDivisors.Commutative.ParentMethods() as domain:
                fraction_map = domain.fraction_field_map()
            case _:
                raise TypeError("rational lifting requires a lattice over an integral domain")
        source_rational: Lattices.ParentMethods = source_ambient.base_change(fraction_map)
        source_perpendicular_inclusion = self.inclusion.base_change(fraction_map)
        coordinates = self.vector.to_vector()
        source_vector = source_rational.linear_combination({label: fraction_map(coordinates(label)) for label in source_ambient.module_generating_set() if coordinates(label)})
        source_norm = source_rational.q(source_vector)
        source_decomposition = {}
        for label in source_rational.module_generating_set():
            source_generator = source_rational.module_generator(label)
            vector_coefficient = source_rational.b(source_generator, source_vector) / source_norm
            perpendicular_part = source_generator - source_rational.scalar_multiple(
                vector_coefficient,
                source_vector,
            )
            source_decomposition[label] = (
                vector_coefficient,
                source_perpendicular_inclusion.lift(perpendicular_part),
            )
        return (
            fraction_map,
            source_rational,
            source_perpendicular_inclusion,
            source_vector,
            source_norm,
            source_decomposition,
        )

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

        if target_section is self:
            (
                fraction_map,
                source_rational,
                source_perpendicular_inclusion,
                source_vector,
                source_norm,
                source_decomposition,
            ) = self._rational_lift_context
            target_rational = source_rational
            target_perpendicular_inclusion = source_perpendicular_inclusion
        else:
            match source_ring:
                case OwnedRings.NoZeroDivisors.Commutative.ParentMethods() as domain:
                    fraction_map = domain.fraction_field_map()
                case _:
                    raise TypeError("rational lifting requires a lattice over an integral domain")
            source_rational = source_ambient.base_change(fraction_map)
            target_rational = target_ambient.base_change(fraction_map)
            source_perpendicular_inclusion = self.inclusion.base_change(fraction_map)
            source_vector = None
            source_norm = None
            source_decomposition = None

        def extend_source(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            result: Lattices.ElementMethods = source_rational.linear_combination(
                {label: fraction_map(coordinates(label)) for label in source_ambient.module_generating_set() if coordinates(label)}
            )
            return result

        def extend_target(
            vector: Lattices.ElementMethods,
        ) -> Lattices.ElementMethods:
            coordinates = vector.to_vector()
            result: Lattices.ElementMethods = target_rational.linear_combination(
                {label: fraction_map(coordinates(label)) for label in target_ambient.module_generating_set() if coordinates(label)}
            )
            return result

        if target_section is not self:
            target_perpendicular_inclusion = target_section.inclusion.base_change(fraction_map)
        reduced_rational = reduced_isometry.base_change(fraction_map)
        if source_vector is None:
            source_vector = extend_source(self.vector)
            source_norm = source_rational.q(source_vector)
        target_vector = source_vector if target_section is self else extend_target(target_section.vector)

        def image(label: Hashable) -> Lattices.ElementMethods:
            if source_decomposition is None:
                source_generator = source_rational.module_generator(label)
                vector_coefficient = source_rational.b(source_generator, source_vector) / source_norm
                perpendicular_part = source_generator - source_rational.scalar_multiple(
                    vector_coefficient,
                    source_vector,
                )
                reduced_part = source_perpendicular_inclusion.lift(perpendicular_part)
            else:
                vector_coefficient, reduced_part = source_decomposition[label]
            target_perpendicular_part = target_perpendicular_inclusion(reduced_rational(reduced_part))
            result: Lattices.ElementMethods = target_perpendicular_part + target_rational.scalar_multiple(
                vector_coefficient,
                target_vector,
            )
            return result

        return IsometryExtensionTorsor(source_rational.Isom(target_rational)(image), ())

    def integral_lift(
        self,
        reduced_isometry: LatticeIsometryMethods,
        *,
        target: NonIsotropicVectorSection | None = None,
    ) -> LatticeIsometryMethods | None:
        r"""Lift a reduced isometry directly when its unique rational extension is integral."""
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

        source_labels = tuple(source_ambient.module_generating_set())
        target_labels = tuple(target_ambient.module_generating_set())
        source_coordinates = self.vector.to_vector()
        target_coordinates = target_section.vector.to_vector()
        source_vector_row = matrix(
            SageQQ,
            1,
            len(source_labels),
            [
                SageQQ(
                    _engine_element(
                        source_ambient.base_ring(),
                        source_coordinates(label),
                    )
                )
                for label in source_labels
            ],
        )
        target_vector_row = matrix(
            SageQQ,
            1,
            len(target_labels),
            [
                SageQQ(
                    _engine_element(
                        target_ambient.base_ring(),
                        target_coordinates(label),
                    )
                )
                for label in target_labels
            ],
        )
        source_perpendicular_rows = _engine_matrix(
            _module_matrix(self.inclusion)
        ).transpose().change_ring(SageQQ)
        target_perpendicular_rows = _engine_matrix(
            _module_matrix(target_section.inclusion)
        ).transpose().change_ring(SageQQ)
        reduced_action = reduced_isometry.parent()._row_action_matrix(
            reduced_isometry
        ).change_ring(SageQQ)
        source_basis = source_vector_row.stack(source_perpendicular_rows)
        image_rows = target_vector_row.stack(
            reduced_action * target_perpendicular_rows
        )
        ambient_action = source_basis.inverse() * image_rows
        if any(entry.denominator() != 1 for entry in ambient_action.list()):
            return None
        lift = source_ambient.Isom(target_ambient)._isometry_from_column_matrix(
            ambient_action.transpose()
        )
        match lift(self.vector) == target_section.vector:
            case True:
                pass
            case False:
                raise ArithmeticError("the integral lift does not carry the selected source vector to the target vector")
        return lift


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
    def _rational_lift_context(self):
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
            source_to_subspace,
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
