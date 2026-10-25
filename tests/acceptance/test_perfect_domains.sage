"""Lorentzian perfect-domain traversal against polyhedral_common's recorded cell counts.

28B_LorentzianPerfStabEqui records, for 40 Lorentzian lattices of rank 3 to 6, the number
of perfect-domain orbits in the "total" and "isotropic" modes. The total mode is reached
through the preamble's reduction complex, whose generators must generate O(L).
"""

import pytest

from sage_indefinite_port.indefinite.lorentzian_cells import perfect_domain_traversal
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
