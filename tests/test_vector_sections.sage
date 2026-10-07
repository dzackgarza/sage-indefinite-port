r"""Acceptance specimens for vector-orthogonal sections and rational lifts."""

from dzack_research.preamble.all import Lattices, ZZ as OwnedZZ

from sage_indefinite_port.indefinite.vector_sections import (
    IsotropicVectorSection,
    NonIsotropicVectorSection,
    orthogonal_section,
)


def test_nonisotropic_section_in_U_lifts_minus_one_to_the_reflection() -> None:
    plane = Lattices(OwnedZZ)("U")
    e, f = plane.module_generators()
    vector = e + f
    section = orthogonal_section(vector)

    assert isinstance(section, NonIsotropicVectorSection)
    assert section.inclusion.domain() is section.perpendicular
    assert section.reduced_object() is section.perpendicular

    (perpendicular_generator,) = section.perpendicular.module_generators()
    reduced_minus_identity = section.perpendicular.O()((-perpendicular_generator,))
    torsor = section.rational_lift(reduced_minus_identity, target=section)
    integral_locus = torsor.integral_parameters(plane, plane)

    assert integral_locus is not None
    lifted = torsor.one_integral_extension()
    reflection = plane.reflection(e - f)
    assert all(
        lifted(generator) == reflection(generator)
        for generator in plane.module_generators()
    )
    assert lifted(vector) == vector


def test_nonisotropic_section_in_U2_has_a_nonempty_integral_locus() -> None:
    two = OwnedZZ.one() + OwnedZZ.one()
    plane = Lattices(OwnedZZ)("U").twist(two)
    e, f = plane.module_generators()
    section = orthogonal_section(e + f)
    assert isinstance(section, NonIsotropicVectorSection)
    (perpendicular_generator,) = section.perpendicular.module_generators()
    reduced_minus_identity = section.perpendicular.O()((-perpendicular_generator,))

    torsor = section.rational_lift(reduced_minus_identity, target=section)

    assert torsor.integral_parameters(plane, plane) is not None


def test_nonisotropic_integral_lift_between_distinct_U_presentations() -> None:
    source = Lattices(OwnedZZ)("U")
    target = Lattices(OwnedZZ)("U")
    source_e, source_f = source.module_generators()
    target_e, target_f = target.module_generators()
    source_vector = source_e + source_f
    target_vector = target_e + target_f
    source_section = orthogonal_section(source_vector)
    target_section = orthogonal_section(target_vector)
    assert isinstance(source_section, NonIsotropicVectorSection)
    assert isinstance(target_section, NonIsotropicVectorSection)
    reduced_isometry = source_section.perpendicular.Isom(
        target_section.perpendicular
    )(
        tuple(target_section.perpendicular.module_generators())
    )

    lifted = source_section.integral_lift(
        reduced_isometry,
        target=target_section,
    )

    assert lifted is not None
    assert lifted.domain() is source
    assert lifted.codomain() is target
    assert lifted(source_vector) == target_vector
    assert lifted(source_e - source_f) == target_e - target_f


def test_isotropic_section_in_U_delegates_to_isotropic_reduction() -> None:
    plane = Lattices(OwnedZZ)("U")
    e, _f = plane.module_generators()
    section = orthogonal_section(e)

    assert isinstance(section, IsotropicVectorSection)
    assert section.reduced_object() is section.reduction
    assert section.perpendicular is section.reduction.orthogonal_complement()
    assert section.inclusion is section.perpendicular.inclusion()


def test_isotropic_section_lifts_the_reduction_identity() -> None:
    lattice = Lattices(OwnedZZ)("U") + Lattices(OwnedZZ)("A1")
    e, _f, _a = lattice.module_generators()
    section = orthogonal_section(e)
    assert isinstance(section, IsotropicVectorSection)
    (reduced_generator,) = section.reduction.module_generators()
    reduced_identity = section.reduction.O()((reduced_generator,))
    torsor = section.rational_lift(reduced_identity, target=section)

    assert torsor.integral_parameters(lattice, lattice) is not None
    lifted = torsor.one_integral_extension()
    assert lifted(e) == e
