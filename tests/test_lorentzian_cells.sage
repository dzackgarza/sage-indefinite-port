r"""Acceptance specimens for Bliss-backed configuration equivalence."""

import pytest
from dzack_research.preamble.all import ZZ as OwnedZZ
from dzack_research.preamble.all import Lattices
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ as SageZZ

from sage_indefinite_port.backends.canonization import (
    CellConfiguration,
    cell_stabilizer,
    cell_transporter,
    presentation_bucket_key,
)
from tests.fixtures.oracle_fixtures import (
    LorentzianEquivalenceCase,
    load_lorentzian_equivalence_cases,
)


def _identity_rows(rank: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(int(row == column) for column in range(rank))
        for row in range(rank)
    )


def _transported_basis_rows(
    case: LorentzianEquivalenceCase,
) -> tuple[tuple[int, ...], ...]:
    witness_rows = case["transporter_witness"]
    if witness_rows is None:
        raise ValueError(f'{case["id"]} has no transporter witness')
    witness = matrix(SageZZ, witness_rows)
    inverse = witness.inverse()
    return tuple(
        tuple(int(inverse[row, column]) for column in range(inverse.ncols()))
        for row in range(inverse.nrows())
    )


_EQUIVALENCE_CASES = load_lorentzian_equivalence_cases()


@pytest.mark.parametrize(
    "case",
    _EQUIVALENCE_CASES,
    ids=[case["id"] for case in _EQUIVALENCE_CASES],
)
def test_lorentzian_equivalence_corpus_recovers_verified_transporters(
    case: LorentzianEquivalenceCase,
) -> None:
    source_lattice = Lattices(OwnedZZ)(case["mat1"])
    target_lattice = Lattices(OwnedZZ)(case["mat2"])
    source = CellConfiguration.from_coordinate_rows(
        source_lattice,
        _identity_rows(case["dimension"]),
    )
    target = CellConfiguration.from_coordinate_rows(
        target_lattice,
        _transported_basis_rows(case),
    )

    transporter = cell_transporter(source, target)

    assert transporter is not None
    assert transporter.domain() is source_lattice
    assert transporter.codomain() is target_lattice
    assert {
        transporter(vector)
        for vector in source.vectors
    } == set(target.vectors)


def test_U_and_U2_have_distinct_presentation_buckets_and_no_transporter() -> None:
    plane = Lattices(OwnedZZ)("U")
    two = OwnedZZ.one() + OwnedZZ.one()
    doubled_plane = plane.twist(two)
    source = CellConfiguration.from_coordinate_rows(plane, _identity_rows(2))
    target = CellConfiguration.from_coordinate_rows(doubled_plane, _identity_rows(2))

    assert presentation_bucket_key(plane) != presentation_bucket_key(doubled_plane)
    assert cell_transporter(source, target) is None


def test_U_basis_configuration_stabilizer_is_realized_by_integral_isometries() -> None:
    plane = Lattices(OwnedZZ)("U")
    configuration = CellConfiguration.from_coordinate_rows(plane, _identity_rows(2))

    stabilizer = cell_stabilizer(configuration)

    assert stabilizer.generators()
    for generator in stabilizer.generators():
        assert {
            generator(vector)
            for vector in configuration.vectors
        } == set(configuration.vectors)
