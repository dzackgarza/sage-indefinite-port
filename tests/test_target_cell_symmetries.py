"""Integral transport between the two octahedral null configurations.

In 3 A1 + A1(-1), the six vectors t +/- e_i are null. Their sum
determines t, so their integral stabilizer is the signed permutation
group on the three negative axes, of order 2^3 * 3! = 48.
Negation transports the opposite configuration onto this one.
"""

import pytest
from sage.all import ZZ, elementary_matrix, identity_matrix, matrix
from sage_extensions.all import Lattices, ZZ as PreambleZZ

from sage_indefinite_port.indefinite import lorentzian_cells
from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
    ZERO_SQUARES,
)
from tests.acceptance.consumer import lattice, rows_of


@pytest.mark.parametrize("shear", (False, True), ids=("canonical_basis", "sheared_basis"))
def test_octahedral_transport_has_one_witness_and_all_signed_permutations(monkeypatch, shear):
    negative = Lattices(PreambleZZ)("A1")
    ambient = negative + negative + negative + negative.twist(PreambleZZ(-1))
    first, second, third, time = tuple(ambient.module_generators())
    target_points = tuple(
        point for root in (first, second, third) for point in (time + root, time - root)
    )
    labels = tuple(ambient.module_generating_set())
    target_rows = tuple(
        tuple(int(point.to_vector()(label)) for label in labels) for point in target_points
    )
    assert all(ambient.b(point, point) == PreambleZZ(0) for point in target_points)
    target_gram = lorentzian_cells._gram_matrix(ambient)
    basis_change = elementary_matrix(ZZ, 4, row1=0, row2=1, scale=1) if shear else identity_matrix(ZZ, 4)
    source_gram = basis_change * target_gram * basis_change.transpose()
    source_rows = tuple(tuple(row) for row in -matrix(ZZ, target_rows) * basis_change.inverse())
    source = LorentzianPerfectCell(lattice(source_gram.rows()), source_rows, ZERO_SQUARES)
    target = LorentzianPerfectCell(ambient, target_rows, ZERO_SQUARES)
    assert set(map(tuple, matrix(ZZ, source_rows) * -basis_change)) == set(target_rows)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_CANONICAL_LABEL_LIMIT", 0)
    backend = LorentzianPerfectLocalBackend()

    witness = backend.cell_transporter(source, target)
    complete = backend.cell_transporter(source, target, all_images=True)

    assert witness is not None
    actions = {
        tuple(tuple(row) for row in rows_of(isometry)) for isometry in complete
    }
    assert len(actions) == 48
    assert tuple(tuple(row) for row in rows_of(witness)) in actions
    for action_rows in actions:
        action = matrix(ZZ, action_rows)
        assert abs(action.det()) == 1
        assert action * target_gram * action.transpose() == source_gram
        assert set(map(tuple, matrix(ZZ, source_rows) * action)) == set(target_rows)
