"""Prefilters checked against cited equivalences and genera."""

import pytest
from dzack_research.preamble.categories.lattices import Lattices
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.invariants import lattice_prefilter
from tests.fixtures.oracle_fixtures import (
    LorentzianEquivalenceCase,
    load_conway_sloane_cases,
    load_lorentzian_equivalence_cases,
)

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
    assert first == second


def test_lattice_prefilter_agrees_on_the_conway_sloane_same_genus_pair() -> None:
    """Conway--Sloane, Chapter 15 §11: (51a) and (51b) lie in one genus, I_{2,1}(2 x 64),
    so every genus invariant agrees although the forms are not integrally equivalent."""
    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first = lattice_prefilter(INTEGRAL_LATTICES(pair["form_a"]["gram"]))
    second = lattice_prefilter(INTEGRAL_LATTICES(pair["form_b"]["gram"]))
    assert first == second
    assert list(first.signature[:2]) == pair["signature"]


def test_lattice_prefilter_accepts_owned_isotropic_reduction() -> None:
    plane = Lattices(ZZ)("U")
    reduction = plane.module_generators()[0].isotropic_reduction()
    prefilter = lattice_prefilter(reduction)
    assert int(prefilter.rank) == int(reduction.rank())
    assert prefilter.discriminant == 1
