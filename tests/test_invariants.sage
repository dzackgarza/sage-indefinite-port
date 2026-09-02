"""Prefilters checked against upstream's equivalent Lorentzian pairs."""

import pytest
from dzack_research.preamble.categories.lattices import Lattices
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.invariants import lattice_prefilter
from tests.fixtures.oracle_fixtures import LorentzianEquivalenceCase, load_lorentzian_equivalence_cases

INTEGRAL_LATTICES = Lattices(ZZ)


@pytest.mark.parametrize("case", load_lorentzian_equivalence_cases(), ids=lambda case: case["id"])
def test_lattice_prefilter_agrees_on_upstream_equivalent_pairs(case: LorentzianEquivalenceCase) -> None:
    """Upstream ships each pair with an exact transporter; no isometry invariant may separate them."""
    first = lattice_prefilter(INTEGRAL_LATTICES(case["mat1"]))
    second = lattice_prefilter(INTEGRAL_LATTICES(case["mat2"]))
    assert first.rank == second.rank
    assert first.signature == second.signature
    assert first.parity == second.parity
    assert first.discriminant == second.discriminant
    assert first.discriminant_elementary_divisors == second.discriminant_elementary_divisors
    assert first == second
