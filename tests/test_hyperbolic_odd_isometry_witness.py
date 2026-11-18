"""U + <-1>^2 is isometric to [2 1; 1 0] + <-1>^2 (#18).

Independent reason: in U with basis e, f (e^2 = f^2 = 0, e.f = 1) the basis e + f, e
has Gram matrix [2 1; 1 0]. The returned witness must carry one Gram matrix to the
other in the row convention of tests/acceptance/consumer.py.
"""

from tests.acceptance.consumer import is_witness, lattice, rows_of

U_PLUS_TWO_MINUS_ONES = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]
SHEARED_PLUS_TWO_MINUS_ONES = [[2, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]


def test_hand_change_of_basis_is_a_witness():
    assert is_witness([[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], SHEARED_PLUS_TWO_MINUS_ONES, U_PLUS_TWO_MINUS_ONES)


def test_returned_witness_carries_one_gram_matrix_to_the_other():
    source, target = lattice(U_PLUS_TWO_MINUS_ONES), lattice(SHEARED_PLUS_TWO_MINUS_ONES)

    assert source.is_isometric(target)
    witness = source.Isom(target).an_element()
    assert is_witness(rows_of(witness), U_PLUS_TWO_MINUS_ONES, SHEARED_PLUS_TWO_MINUS_ONES)
