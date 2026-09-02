"""Attack profiles and prefilters checked against the oracle corpus."""

import pytest
from dzack_research.preamble.categories.lattices import Lattices
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.invariants import (
    ZeroVectorError,
    attack_profile,
    content,
    lattice_prefilter,
    vector_prefilter,
)
from tests.fixtures.lattice_corpus import CORPUS, LatticeEntry
from tests.fixtures.oracle_fixtures import LorentzianEquivalenceCase, load_lorentzian_equivalence_cases

INTEGRAL_LATTICES = Lattices(ZZ)


@pytest.mark.parametrize("entry", CORPUS, ids=lambda entry: entry["name"])
def test_attack_profile_positive_index_is_the_smaller_signature_entry(entry: LatticeEntry) -> None:
    positive, negative, _degenerate = entry["signature"]
    profile = attack_profile(INTEGRAL_LATTICES(entry["gram"]))
    assert profile.positive_index == min(positive, negative)
    assert profile.negative_index == max(positive, negative)
    assert profile.sign == (1 if positive <= negative else -1)
    assert profile.signed_view.signature_pair() == (min(positive, negative), max(positive, negative))
    if profile.sign == 1:
        assert profile.signed_view is profile.lattice


@pytest.mark.parametrize("case", load_lorentzian_equivalence_cases(), ids=lambda case: case["id"])
def test_lattice_prefilter_agrees_on_upstream_equivalent_pairs(case: LorentzianEquivalenceCase) -> None:
    first = lattice_prefilter(INTEGRAL_LATTICES(case["mat1"]))
    second = lattice_prefilter(INTEGRAL_LATTICES(case["mat2"]))
    assert first == second


def test_lattice_prefilter_separates_discriminant_parity_and_signature() -> None:
    plane = lattice_prefilter(INTEGRAL_LATTICES("U"))
    scaled = lattice_prefilter(INTEGRAL_LATTICES("U").twist(2))
    odd = lattice_prefilter(INTEGRAL_LATTICES([[1, 0], [0, -1]]))
    negative_definite = lattice_prefilter(INTEGRAL_LATTICES([[-2]]))
    positive_definite = lattice_prefilter(INTEGRAL_LATTICES([[2]]))

    assert plane.discriminant_elementary_divisors == ()
    assert scaled.discriminant_elementary_divisors == (2, 2)
    assert plane.parity == "even" and odd.parity == "odd"
    assert plane.signature == (1, 1, 0)
    assert negative_definite.signature == (0, 1, 0) and positive_definite.signature == (1, 0, 0)
    assert len({plane, scaled, odd, negative_definite, positive_definite}) == 5


def test_vector_prefilter_separates_content_from_divisor() -> None:
    scaled = INTEGRAL_LATTICES("U").twist(2)
    generator = scaled.module_generator(0)
    isotropic = vector_prefilter(scaled, generator)
    doubled = vector_prefilter(scaled, 2 * generator)
    assert (isotropic.norm, isotropic.content, isotropic.divisor) == (0, 1, 2)
    assert (doubled.norm, doubled.content, doubled.divisor) == (0, 2, 4)

    plane = INTEGRAL_LATTICES("U")
    unimodular = vector_prefilter(plane, plane.module_generator(0) + plane.module_generator(1))
    assert (unimodular.norm, unimodular.content, unimodular.divisor) == (2, 1, 1)


def test_content_of_the_zero_vector_is_rejected() -> None:
    plane = INTEGRAL_LATTICES("U")
    with pytest.raises(ZeroVectorError):
        content(0 * plane.module_generator(0))
