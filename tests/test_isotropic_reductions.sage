r"""Acceptance specimens for the owned isotropic reduction K_I = I^perp / I.

The reduction itself belongs to dzack_research.preamble.  This port consumes
that public construction rather than introducing a second lattice or quotient
owner.
"""

import pytest

from dzack_research.preamble.all import *


def test_primitive_isotropic_line_of_U_plus_E8_reduces_to_E8() -> None:
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)("E8")
    root_lattice = Lattices(ZZ)("E8")
    isotropic = lattice.module_generators()[0]

    assert isotropic.is_isotropic()
    assert isotropic.is_primitive()
    assert isotropic.isotropic_reduction().is_isometric(root_lattice)


def test_reduction_equality_excludes_the_coordinate_complement() -> None:
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)("A2")
    isotropic = lattice.module_generators()[0]
    reduction = isotropic.isotropic_reduction()
    labels = tuple(reduction.module_generating_set())
    frame = reduction.coordinate_frame()
    line_generator = reduction.isotropic_sublattice().module_generators()[0]
    isotropic_in_perpendicular = reduction.isotropic_inclusion()(line_generator)
    shifted_frame = {
        label: frame[label] + isotropic_in_perpendicular
        if label == labels[0]
        else frame[label]
        for label in labels
    }

    shifted = reduction.with_coordinate_frame(shifted_frame)

    assert shifted.coordinate_frame()[labels[0]] != frame[labels[0]]
    assert shifted == reduction
    assert hash(shifted) == hash(reduction)


def test_nonprimitive_isotropic_vector_retains_torsion_only_in_formed_quotient() -> None:
    isotropic = Lattices(ZZ)("U").module_generators()[0]
    nonprimitive = 2 * isotropic

    with pytest.raises(NotPrimitiveError):
        nonprimitive.isotropic_reduction()

    quotient = nonprimitive.isotropic_quotient()

    assert quotient in FormModules(ZZ)
    assert not quotient.is_free()
    assert tuple(quotient.invariant_factors()) == (ZZ(2),)
