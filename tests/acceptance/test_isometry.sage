"""Isometry decisions and witnesses between indefinite lattices.

Verdicts come from certified corpus pairs (Indefinite.jl and Hecke tests, polyhedral_common's
145 equivalent Lorentzian pairs, the Conway--Sloane same-genus pair). A returned witness is
accepted when it carries one Gram matrix to the other. polyhedral_common's 19_IndefiniteComp
suite is adapted as it is: a lattice and its image under a unimodular change of basis are
isometric.
"""

import pytest

from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import (
    ISOMETRY,
    is_witness,
    lattice,
    require,
    rows_of,
)
from tests.fixtures.oracle_fixtures import (
    load_ci_indefinite_comp,
    load_conway_sloane_cases,
    load_indefinite_isometry_pairs,
    load_lorentzian_equivalence_cases,
)


def _assert_isometric(source_gram, target_gram) -> None:
    source, target = lattice(source_gram), lattice(target_gram)

    assert source.is_isometric(target)
    witness = source.Isom(target).an_element()
    assert is_witness(rows_of(witness), source_gram, target_gram)


@pytest.mark.parametrize("case", load_indefinite_isometry_pairs(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="O(L), isometry, vector stabilizers: owned by #18", raises=UnfinishedCapability, strict=True)
def test_isometry_verdicts_match_certified_pairs(case) -> None:
    require(ISOMETRY)

    if case["isometric"]:
        _assert_isometric(case["gram1"], case["gram2"])
    else:
        assert not lattice(case["gram1"]).is_isometric(lattice(case["gram2"]))
        assert lattice(case["gram1"]).Isom(lattice(case["gram2"])).is_empty()


@pytest.mark.parametrize("case", load_lorentzian_equivalence_cases(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="O(L), isometry, vector stabilizers: owned by #18", raises=UnfinishedCapability, strict=True)
def test_lorentzian_pairs_are_isometric(case) -> None:
    require(ISOMETRY)

    assert is_witness(case["transporter_witness"], case["mat2"], case["mat1"])

    _assert_isometric(case["mat1"], case["mat2"])


@pytest.mark.xfail(reason="O(L), isometry, vector stabilizers: owned by #18", raises=UnfinishedCapability, strict=True)
def test_conway_sloane_same_genus_pair_is_not_isometric() -> None:
    require(ISOMETRY)

    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first, second = lattice(pair["form_a"]["gram"]), lattice(pair["form_b"]["gram"])

    assert first.genus() == second.genus()
    assert not first.is_isometric(second)
    assert first.Isom(second).is_empty()


def _unimodular_change_of_basis(rank: int) -> list[list[int]]:
    """The upper unitriangular all-ones matrix: unimodular and mixes every basis vector."""
    return [[1 if j >= i else 0 for j in range(rank)] for i in range(rank)]


@pytest.mark.parametrize("case", load_ci_indefinite_comp(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="O(L), isometry, vector stabilizers: owned by #18", raises=UnfinishedCapability, strict=True)
def test_change_of_basis_gives_an_isometric_lattice(case) -> None:
    require(ISOMETRY)

    gram = matrix(ZZ, case["gram"])
    change = matrix(ZZ, _unimodular_change_of_basis(gram.nrows()))
    conjugate = change * gram * change.transpose()

    _assert_isometric(case["gram"], [list(row) for row in conjugate.rows()])
