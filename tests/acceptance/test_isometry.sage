"""Isometry decisions and witnesses between indefinite lattices.

Verdicts come from certified corpus pairs (Indefinite.jl and Hecke tests, polyhedral_common's
145 equivalent Lorentzian pairs, the Conway--Sloane same-genus pair). A returned witness is
accepted when it carries one Gram matrix to the other. polyhedral_common's 19_IndefiniteComp
suite is adapted as it is: a lattice and its image under a unimodular change of basis are
isometric.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_ISOMETRY,
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
def test_isometry_verdicts_match_certified_pairs(case) -> None:
    require(MARKED_ISOMETRY)

    if case["isometric"]:
        _assert_isometric(case["gram1"], case["gram2"])
    else:
        assert not lattice(case["gram1"]).is_isometric(lattice(case["gram2"]))
        assert lattice(case["gram1"]).Isom(lattice(case["gram2"])).is_empty()


@pytest.mark.parametrize("case", load_lorentzian_equivalence_cases(), ids=lambda case: case["id"])
def test_lorentzian_pairs_are_isometric(case) -> None:
    require(MARKED_ISOMETRY)

    assert is_witness(case["transporter_witness"], case["mat2"], case["mat1"])

    _assert_isometric(case["mat1"], case["mat2"])


def test_equal_gram_groupoid_identification_preserves_coordinate_labels() -> None:
    """Memoized equal-Gram transport is coordinate identity, not iterator-position identity."""
    from sage_indefinite_port.indefinite.groupoid_cache import coordinate_identification

    source = lattice([[2, 1], [1, 2]])
    target = lattice([[2, 1], [1, 2]])
    witness = coordinate_identification(source, target)
    source_labels = tuple(source.module_generating_set())
    target_labels = tuple(target.module_generating_set())

    for index, source_label in enumerate(source_labels):
        image = witness(source.module_generator(source_label)).to_vector()
        assert tuple(int(image(label)) for label in target_labels) == tuple(int(index == j) for j in range(len(target_labels)))


def test_degenerate_isometry_witness_respects_complement_coordinates() -> None:
    """A non-diagonal nonsingular summand plus radical survives a basis change."""
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    require(MARKED_ISOMETRY)

    gram = matrix(SageZZ, [[2, 1, 0], [1, 2, 0], [0, 0, 0]])
    change = matrix(SageZZ, [[1, 1, 0], [0, 1, 1], [0, 0, 1]])
    conjugate = change * gram * change.transpose()
    _assert_isometric([list(row) for row in gram.rows()], [list(row) for row in conjugate.rows()])


def test_conway_sloane_same_genus_pair_is_not_isometric() -> None:
    require(MARKED_ISOMETRY)

    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first, second = lattice(pair["form_a"]["gram"]), lattice(pair["form_b"]["gram"])

    assert first.genus() == second.genus()
    assert not first.is_isometric(second)
    assert first.Isom(second).is_empty()


def _unimodular_change_of_basis(rank: int) -> list[list[int]]:
    """The upper unitriangular all-ones matrix: unimodular and mixes every basis vector."""
    return [[1 if j >= i else 0 for j in range(rank)] for i in range(rank)]


@pytest.mark.parametrize("case", load_ci_indefinite_comp(), ids=lambda case: case["id"])
def test_change_of_basis_gives_an_isometric_lattice(case) -> None:
    require(MARKED_ISOMETRY)

    gram = matrix(ZZ, case["gram"])
    change = matrix(ZZ, _unimodular_change_of_basis(gram.nrows()))
    conjugate = change * gram * change.transpose()

    _assert_isometric(case["gram"], [list(row) for row in conjugate.rows()])
