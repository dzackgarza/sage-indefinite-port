"""Exact integral Gram pairings of hyperbolic-plane configurations."""

from __future__ import annotations

import pytest

from sage.all import ZZ, matrix

from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectCell,
    NONNEGATIVE_SQUARES,
    ZERO_SQUARES,
    _ConfigurationPairings,
    _configuration_pairing_block_size,
    _perfect_form_hash_key,
)
from tests.acceptance.consumer import lattice


def test_hyperbolic_plane_pairings_preserve_integral_change_of_basis() -> None:
    """In U, b((a,b),(c,d)) = ad+bc, in every integral basis."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    large = 2**80 + 1
    rows = ((1, 0), (0, 1), (large, 1), (1, -large))
    expected = matrix(ZZ, [[a * d + b * c for c, d in rows] for a, b in rows])
    shear = matrix(ZZ, [[1, large], [0, 1]])
    for coordinates, form in (
        (matrix(ZZ, rows), gram),
        (matrix(ZZ, rows) * shear, shear.inverse() * gram * shear.inverse().transpose()),
    ):
        pairings = _ConfigurationPairings(tuple(tuple(row) for row in coordinates), form)
        assert len(pairings) == len(rows)
        assert tuple(pairings.diagonal()) == tuple(expected.diagonal())
        for i in range(len(rows)):
            assert tuple(pairings.column(i)) == tuple(expected.column(i))
            for j in range(len(rows)):
                assert pairings[i, j] == expected[i, j]
        for start, block in pairings.blocks():
            assert matrix(ZZ, block.tolist()) == expected[start : start + len(block), :]

    # The backtracking cache must stay bounded even if every column is read;
    # otherwise streaming columns eventually retains the full Gram matrix.
    rows = tuple((i, 1) for i in range(80))
    streaming = _ConfigurationPairings(rows, gram)
    for i in range(len(rows)):
        assert tuple(streaming.column(i)) == tuple(i + j for j in range(len(rows)))
    assert not streaming._column_cache

    cached_pairings = _ConfigurationPairings(rows, gram, cache_columns=True)
    for start, block in cached_pairings.blocks():
        assert 0 < len(block) <= 16
        assert matrix(ZZ, block.tolist()) == matrix(
            ZZ, [[i + j for j in range(len(rows))] for i in range(start, start + len(block))]
        )
    assert not cached_pairings._column_cache

    pairings = _ConfigurationPairings(rows, gram, cache_columns=True)
    for i in range(len(rows)):
        assert tuple(pairings.column(i)) == tuple(i + j for j in range(len(rows)))
    assert len(pairings._column_cache) == 64
    # An evicted column is recomputed exactly and re-enters the bounded cache.
    assert tuple(pairings.column(0)) == tuple(range(len(rows)))
    assert len(pairings._column_cache) == 64


def test_configuration_pairing_blocks_have_a_bounded_exact_working_set() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        _configuration_pairing_block_size(-1)
    assert _configuration_pairing_block_size(0) == 16
    assert _configuration_pairing_block_size(49_861) == 16
    assert _configuration_pairing_block_size(100_000) == 10
    assert _configuration_pairing_block_size(2_000_000) == 1


def test_isotropic_cell_invariants_preserve_integral_change_of_basis() -> None:
    """The primitive null basis of U has norms zero and mutual pairing one."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    shear = matrix(ZZ, [[1, 2**80 + 1], [0, 1]])
    changed_gram = shear.inverse() * gram * shear.inverse().transpose()
    rows = matrix(ZZ, [[1, 0], [0, 1]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), tuple(tuple(row) for row in rows), ZERO_SQUARES)
    changed_cell = LorentzianPerfectCell(
        lattice(changed_gram.rows()), tuple(tuple(row) for row in rows * shear), ZERO_SQUARES
    )
    expected = (((0, 2),), ((1, 1),), (1, 1), ((1, 2),), ((1, 2),))
    assert _perfect_form_hash_key(cell) == expected
    assert _perfect_form_hash_key(changed_cell) == expected


def test_large_coarse_hash_preserves_gcd_colours_under_integral_isometry(monkeypatch) -> None:
    """Linear-cost content colours are invariant under an integral basis change."""
    from sage_indefinite_port.indefinite import lorentzian_cells

    gram = matrix(ZZ, [[0, 1], [1, 0]])
    shear = matrix(ZZ, [[1, 2**80 + 1], [0, 1]])
    changed_gram = shear.inverse() * gram * shear.inverse().transpose()
    rows = matrix(ZZ, [[2, 0], [0, 3]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), tuple(tuple(row) for row in rows), ZERO_SQUARES)
    changed_cell = LorentzianPerfectCell(
        lattice(changed_gram.rows()), tuple(tuple(row) for row in rows * shear), ZERO_SQUARES
    )
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)

    key = _perfect_form_hash_key(cell)

    assert key == _perfect_form_hash_key(changed_cell)
    assert key[1] == ()
    assert sum(count for _value, count in key[3]) == len(cell.rows)
    assert sum(count for _value, count in key[4]) == len(cell.rows)


def test_nonnegative_square_face_has_the_exact_pairing_multisets() -> None:
    """The face x=1 of x^2-y^2 has norms (0,1,0) and pairings (1,2,1)."""
    cell = LorentzianPerfectCell(
        lattice([[1, 0], [0, -1]]), ((1, -1), (1, 0), (1, 1)), NONNEGATIVE_SQUARES
    )
    assert _perfect_form_hash_key(cell) == (
        ((0, 2), (1, 1)), ((1, 2), (2, 1)), (1, 1), ((1, 3),), ((1, 3),)
    )
