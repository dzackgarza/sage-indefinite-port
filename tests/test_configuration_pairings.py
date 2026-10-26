"""Exact integral Gram pairings of hyperbolic-plane configurations."""

from __future__ import annotations

from sage.all import ZZ, matrix

from sage_indefinite_port.indefinite.lorentzian_cells import _ConfigurationPairings


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
