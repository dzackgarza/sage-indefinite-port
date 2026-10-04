r"""Acceptance specimens for isotropic flag types."""

from dzack_research.preamble.all import Lattices, ZZ as OwnedZZ

from sage_indefinite_port.indefinite.isotropic_reductions import FlagType


def test_plane_and_complete_flag_types_have_expected_blocks() -> None:
    assert FlagType.plane(2).dimensions == (2,)
    assert FlagType.plane(2).block_sizes == (2,)
    assert FlagType.complete(3).dimensions == (1, 2, 3)
    assert FlagType.complete(3).block_sizes == (1, 1, 1)


def test_two_step_isotropic_flag_recovers_its_flag_type() -> None:
    lattice = Lattices(OwnedZZ)("U") + Lattices(OwnedZZ)("U")
    first = lattice.module_generators()[0]
    second = lattice.module_generators()[2]
    flag = lattice.isotropic_flag(first, second)

    assert flag.flag_type() == FlagType((1, 2))
    assert flag.flag_type().block_sizes == (1, 1)
