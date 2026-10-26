"""Lorentzian perfect-domain traversal against polyhedral_common's recorded cell counts.

28B_LorentzianPerfStabEqui records, for 40 Lorentzian lattices of rank 3 to 6, the number
of perfect-domain orbits in the "total" and "isotropic" modes. The total mode is reached
through the preamble's reduction complex, whose generators must generate O(L).
"""

from dataclasses import replace

import pytest

from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectLocalBackend,
    NONNEGATIVE_SQUARES,
    perfect_complex,
    perfect_domain_traversal,
)
from tests.acceptance.consumer import PERFECT_DOMAINS, lattice, require
from tests.fixtures.oracle_fixtures import load_lorentzian_perfect_domains


@pytest.mark.parametrize("case", load_lorentzian_perfect_domains(), ids=lambda case: case["id"])
def test_perfect_domain_orbit_counts_match_polyhedral_common(case) -> None:
    require(PERFECT_DOMAINS)

    traversal = lattice(case["gram"]).lorentzian_reduction_complex()

    assert traversal.is_complete()
    assert len(traversal.cells()) == case["total_count"]
    assert traversal.generates_orthogonal_group()


@pytest.mark.parametrize("case", load_lorentzian_perfect_domains(), ids=lambda case: case["id"])
def test_isotropic_mode_counts_match_polyhedral_common(case) -> None:
    """The "isotropic" mode is a reference-implementation option with no preamble consumer,
    so the case calls the port's entry point directly, on the preamble lattice."""
    require(PERFECT_DOMAINS)

    assert len(perfect_domain_traversal(lattice(case["gram"]), "isotropic").cells()) == case["isotropic_count"]


def test_lazy_neighbor_cones_retain_exact_transformed_rays() -> None:
    """Materializing a deferred side-pairing cone checks its exact ray certificate."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    traversal = lattice(case["gram"]).lorentzian_reduction_complex()
    adjacencies = tuple(traversal.adjacencies())
    assert adjacencies
    for adjacency in adjacencies:
        neighbor = adjacency.neighbor()
        common_face = adjacency.common_face()
        assert neighbor.ambient_lattice() is traversal.lattice()
        assert common_face.ambient_lattice() is traversal.lattice()
        assert int(common_face.dimension()) == int(traversal.lattice().rank()) - 1
        assert common_face.is_face_of(adjacency.source())
        assert common_face.is_face_of(neighbor)
        assert adjacency.target_representative().transport(adjacency.target_to_neighbor()).is_equal_to(neighbor)


def test_facet_orbit_permutation_view_materializes_exact_isometries() -> None:
    """The public owned-isometry view agrees with the internal Schreier permutations."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(lattice(case["gram"]), NONNEGATIVE_SQUARES)
    permutation_orbits = backend.facet_orbit_permutations(cell)
    isometry_orbits = backend.facet_orbits(cell)

    assert tuple(tuple(facet for facet, _ in orbit) for orbit in permutation_orbits) == tuple(
        tuple(facet for facet, _ in orbit) for orbit in isometry_orbits
    )
    for permutation_orbit, isometry_orbit in zip(permutation_orbits, isometry_orbits, strict=True):
        for (facet, permutation), (isometry_facet, isometry) in zip(permutation_orbit, isometry_orbit, strict=True):
            assert facet == isometry_facet
            for position, image_position in enumerate(permutation):
                assert isometry(cell.vector_configuration[position]) == cell.vector_configuration[image_position]


def test_large_cell_cone_uses_exact_normaliz_h_presentation() -> None:
    """The large-configuration backend preserves the owned exact cone interface."""
    require(PERFECT_DOMAINS)

    case = next(case for case in load_lorentzian_perfect_domains() if case["id"] == "Result_Enumeration_1")
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(lattice(case["gram"]), NONNEGATIVE_SQUARES)
    # Repeating the same rays keeps the cone unchanged while forcing the
    # large-configuration Normaliz path without a costly high-rank fixture.
    large_cell = replace(cell, rows=cell.rows * 4_000)
    assert len(large_cell.rows) >= 10_000
    large_cone = backend.cell_cone(large_cell)
    small_cone = backend.cell_cone(cell)
    labels = tuple(cell.lattice.module_generating_set())
    assert large_cone._engine_polyhedron().backend() == "normaliz"

    def covector_rows(cone):
        return {
            tuple(int(covector.to_vector()(label)) for label in labels)
            for covector in cone.halfspace_covectors()
        }

    assert covector_rows(large_cone) == covector_rows(small_cone)
    assert large_cone.dimension() == small_cone.dimension()
    assert set(large_cone.primitive_rays()) == set(small_cone.primitive_rays())
    assert len(large_cone.faces(int(large_cone.dimension()) - 1)) == len(small_cone.faces(int(small_cone.dimension()) - 1))
