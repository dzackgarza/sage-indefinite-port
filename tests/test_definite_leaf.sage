r"""Acceptance specimens for the delegated definite leaf."""

import pytest

from dzack_research.preamble.all import Lattices, ZZ as OwnedZZ

from sage_indefinite_port.backends.definite import (
    definite_isometry,
    definite_orthogonal_group,
)
from tests.fixtures.oracle_fixtures import load_conway_sloane_cases


_ROOT_CASES = load_conway_sloane_cases()["root_lattice_automorphism_orders"]["cases"]


@pytest.mark.parametrize(
    "case",
    _ROOT_CASES,
    ids=[case["id"] for case in _ROOT_CASES],
)
def test_conway_sloane_definite_orthogonal_group_orders(case) -> None:
    lattice = Lattices(OwnedZZ)(case["gram"])
    group = definite_orthogonal_group(lattice)
    assert int(group.cardinality()) == case["order"]


def test_E8_isometry_witness_preserves_forms_and_composes_with_inverse() -> None:
    standard = next(case for case in _ROOT_CASES if case["id"] == "cs99_E8")
    upstream = next(case for case in _ROOT_CASES if case["id"] == "indefinite_jl_E8")
    source = Lattices(OwnedZZ)(standard["gram"])
    target = Lattices(OwnedZZ)(upstream["gram"])

    witness = definite_isometry(source, target)
    inverse = ~witness

    assert witness.domain() is source
    assert witness.codomain() is target
    for left in source.module_generators():
        for right in source.module_generators():
            assert source.b(left, right) == target.b(witness(left), witness(right))
        assert inverse(witness(left)) == left
