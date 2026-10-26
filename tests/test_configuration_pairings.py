"""Exact integral Gram pairings of hyperbolic-plane configurations."""

from __future__ import annotations

from sage.all import ZZ, matrix

from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectCell,
    ZERO_SQUARES,
    _ConfigurationPairings,
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


def test_isotropic_cell_invariants_preserve_integral_change_of_basis() -> None:
    """The primitive null basis of U has norms zero and mutual pairing one."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    shear = matrix(ZZ, [[1, 2**80 + 1], [0, 1]])
    changed_gram = shear.inverse() * gram * shear.inverse().transpose()
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    changed_cell = LorentzianPerfectCell(
        lattice(changed_gram.rows()), tuple(tuple(row) for row in shear), ZERO_SQUARES
    )
    expected = (((0, 2),), ((1, 1),), (1, 1), (1, 1))
    assert _perfect_form_hash_key(cell) == expected
    assert _perfect_form_hash_key(changed_cell) == expected
