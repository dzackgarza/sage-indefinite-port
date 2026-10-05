from dzack_research.preamble.all import Lattices, ZZ

from sage_indefinite_port.indefinite.recursive import IndefiniteOrthogonalAlgorithm
from tests.fixtures.oracle_fixtures import load_conway_sloane_cases


def test_conway_sloane_51a_51b_public_isometry_rejects() -> None:
    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first = Lattices(ZZ)(pair["form_a"]["gram"])
    second = Lattices(ZZ)(pair["form_b"]["gram"])

    witness = IndefiniteOrthogonalAlgorithm().isometry(first, second)

    assert pair["integrally_equivalent"] is False
    assert witness is None
