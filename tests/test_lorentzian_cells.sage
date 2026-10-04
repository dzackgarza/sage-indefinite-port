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
from sage_indefinite_port.backends.polyhedral import (
    configuration_cone,
    configuration_facets,
    facet_orbits,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectComplex,
    LorentzianPerfectLocalBackend,
    MarkedCellOrbitAlgorithm,
    NoLocalMarkTheoremError,
    perfect_domain_traversal,
)
from tests.fixtures.oracle_fixtures import (
    LorentzianEquivalenceCase,
    load_lorentzian_equivalence_cases,
    load_lorentzian_perfect_domains,
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


def test_rank_three_perfect_cell_facets_and_flip_stay_in_the_same_complex() -> None:
    """Recover the rank-three source cell and its unique facet orbit."""
    case = next(
        case
        for case in load_lorentzian_perfect_domains()
        if case["id"] == "lor_perf_dim3_case04"
    )
    lattice = Lattices(OwnedZZ)(case["gram"])
    source = CellConfiguration.from_coordinate_rows(
        lattice,
        (
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, -1),
            (2, 1, -2),
            (1, 2, -2),
        ),
    )
    neighbor = CellConfiguration.from_coordinate_rows(
        lattice,
        (
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 1),
            (2, 1, 2),
            (1, 2, 2),
        ),
    )

    source_cone = configuration_cone(source)
    neighbor_cone = configuration_cone(neighbor)
    facets = configuration_facets(source)
    orbits = facet_orbits(source, cell_stabilizer(source))

    assert source_cone.primitive_rays().cardinality() == 4
    assert len(facets) == 4
    assert sum(len(orbit) for orbit in orbits) == len(facets)
    assert len(orbits) == 1
    assert frozenset({0, 1}) in orbits[0]
    assert source_cone.is_adjacent_to(neighbor_cone)
    transporter = cell_transporter(source, neighbor)
    assert transporter is not None
    assert {
        transporter(vector)
        for vector in source.vectors
    } == set(neighbor.vectors)
    assert case["total_count"] == 1


def test_native_rank_three_local_backend_builds_and_flips_both_cell_modes() -> None:
    case = next(
        case
        for case in load_lorentzian_perfect_domains()
        if case["id"] == "lor_perf_dim3_case01"
    )
    lattice = Lattices(OwnedZZ)(case["gram"])
    backend = LorentzianPerfectLocalBackend()

    isotropic = backend.initial_cell(lattice, "isotropic")
    total = backend.initial_cell(lattice, "total")
    orbits = backend.facet_orbits(total)
    neighbor = backend.flip_across(total, orbits[0][0])
    transporter = backend.cell_transporter(total, neighbor)

    assert all(vector.q() == 0 for vector in isotropic.vector_configuration)
    assert sum(len(orbit) for orbit in orbits) == 4
    assert transporter is not None
    assert {
        transporter(vector)
        for vector in total.vector_configuration
    } == set(neighbor.vector_configuration)


def test_rank_three_quotient_traversal_and_groups_match_frozen_counts() -> None:
    case = next(
        case
        for case in load_lorentzian_perfect_domains()
        if case["id"] == "lor_perf_dim3_case01"
    )
    lattice = Lattices(OwnedZZ)(case["gram"])

    total = LorentzianPerfectComplex(lattice, "total")
    isotropic = LorentzianPerfectComplex(lattice, "isotropic")

    total_cells = total.quotient_cells()
    isotropic_cells = isotropic.quotient_cells()
    assert len(total_cells) == case["total_count"]
    assert len(isotropic_cells) == case["isotropic_count"]

    for adjacency in total.adjacencies():
        assert adjacency.transporter is not None
        transported = {
            adjacency.transporter(vector)
            for vector in adjacency.target.vector_configuration
        }
        flipped = total.local_backend().flip_across(adjacency.source, adjacency.facet)
        assert transported == set(flipped.vector_configuration)

    component = total.component_preserving_group()
    full = total.full_orthogonal_group()
    assert full.generators()[:-1] == component.generators()
    negation = full.generators()[-1]
    assert all(
        negation(generator) == -generator
        for generator in lattice.module_generators()
    )

    records = perfect_domain_traversal(case["gram"], "total")
    assert len(records) == case["total_count"]
    assert all(record["x"]["EXT"] for record in records)
    assert all(record["ListAdj"] for record in records)


def test_rank_three_marked_isotropic_orbits_follow_cell_adjacencies() -> None:
    case = next(
        case
        for case in load_lorentzian_perfect_domains()
        if case["id"] == "lor_perf_dim3_case01"
    )
    lattice = Lattices(OwnedZZ)(case["gram"])
    complex_ = LorentzianPerfectComplex(lattice, "total")
    algorithm = MarkedCellOrbitAlgorithm(complex_)

    orbits = algorithm.global_orbits()

    assert len(orbits) == 1
    assert orbits[0]
    assert all(vector.q() == 0 for vector in orbits[0])
    with pytest.raises(NoLocalMarkTheoremError):
        MarkedCellOrbitAlgorithm(complex_, norm=2)
