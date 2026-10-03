r"""Acceptance specimens for vector-orthogonal sections and rational lifts."""

from dzack_research.preamble.all import Lattices, ZZ as OwnedZZ

from sage_indefinite_port.indefinite.vector_sections import (
    NonIsotropicVectorSection,
    orthogonal_section,
)


def test_nonisotropic_section_in_U_lifts_minus_one_to_the_reflection() -> None:
    plane = Lattices(OwnedZZ)("U")
    e, f = plane.module_generators()
    vector = e + f
    section = orthogonal_section(vector)

    assert isinstance(section, NonIsotropicVectorSection)
    assert section.inclusion is section.perpendicular.inclusion()
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
    (perpendicular_generator,) = section.perpendicular.module_generators()
    reduced_minus_identity = section.perpendicular.O()((-perpendicular_generator,))

    torsor = section.rational_lift(reduced_minus_identity, target=section)

    assert torsor.integral_parameters(plane, plane) is not None
