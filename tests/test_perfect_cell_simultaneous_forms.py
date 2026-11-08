"""Both forms of the pinned configuration isomorphism determine transport.

For diag(5,-1,-2), the first positive pairing with (1,0,0) is 5.
Its critical vectors are (1,x,y) with x^2+2*y^2 <= 5. The unique
norm-five vector fixes the time axis. The integral spatial group of
diag(-1,-2) consists of the four sign changes, so these are all cell
automorphisms. The inverse second-moment form adds independent data.
"""

from itertools import product

from sage.all import ZZ, diagonal_matrix, elementary_matrix, matrix

from sage_indefinite_port.indefinite import lorentzian_cells
from tests.acceptance.consumer import lattice, rows_of


def test_asymmetric_moment_cell_has_its_complete_integral_transporters(monkeypatch):
    gram = diagonal_matrix(ZZ, (5, -1, -2))
    ambient = lattice(gram)
    rows = tuple(
        (1, x, y) for x in range(-2, 3) for y in range(-1, 2)
        if x * x + 2 * y * y <= 5
    )
    enumerator = lorentzian_cells._PositiveVectorEnumerator(
        ambient, (1, 0, 0), 5, lorentzian_cells.NONNEGATIVE_SQUARES,
    )
    first = enumerator.first_shell()
    assert first is not None and first[0] == 1
    assert set(first[1]) == set(rows)
    moment = matrix(ZZ, rows).transpose() * matrix(ZZ, rows)
    assert moment == diagonal_matrix(ZZ, (11, 14, 6))
    auxiliary = moment.adjugate()
    assert auxiliary[1, 1] * gram[2, 2] != auxiliary[2, 2] * gram[1, 1]

    change = elementary_matrix(ZZ, 3, row1=1, row2=2, scale=1)
    source_gram = change * gram * change.transpose()
    source_rows = tuple(tuple(row) for row in -matrix(ZZ, rows) * change.inverse())
    source = lorentzian_cells.LorentzianPerfectCell(
        lattice(source_gram), source_rows, lorentzian_cells.NONNEGATIVE_SQUARES,
    )
    target = lorentzian_cells.LorentzianPerfectCell(
        ambient, rows, lorentzian_cells.NONNEGATIVE_SQUARES,
    )
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_CANONICAL_LABEL_LIMIT", 0)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    backend = lorentzian_cells.LorentzianPerfectLocalBackend()
    backend.cell_stabilizer(source)
    witness = backend.cell_transporter(source, target)
    complete = backend.cell_transporter(source, target, all_images=True)
    expected = {
        tuple(map(tuple, -change * diagonal_matrix(ZZ, (1, x, y))))
        for x, y in product((-1, 1), repeat=2)
    }
    assert witness is not None
    assert tuple(map(tuple, rows_of(witness))) in expected
    assert {tuple(map(tuple, rows_of(arrow))) for arrow in complete} == expected
