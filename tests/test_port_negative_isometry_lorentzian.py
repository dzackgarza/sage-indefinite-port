"""The port rejects a same-genus non-isometric pair of signature (2, 1).

EichlerRecursion.isometry calls _indefinite_isometry, whose genus and spinor-genus
decision applies to integral lattices of rank at least three, including odd forms.
The Conway--Sloane forms (51a) and (51b) are rejected there before the Lorentzian
witness search. This tests the port entry point's negative decision.

The departure from the source's exhaustive search is recorded at
https://github.com/dzackgarza/sage-indefinite-port/issues/68#issuecomment-6821405364.
By Eichler--Kneser, each indefinite spinor genus of rank at least three contains
one integral isometry class. The independent negative expectation is
Conway--Sloane, SPLAG 3rd ed., Ch. 15, section 11: (51a) and (51b) lie in one
genus and are not integrally equivalent.
"""

from tests.acceptance.consumer import lattice
from tests.fixtures.oracle_fixtures import load_conway_sloane_cases
from sage_indefinite_port.indefinite.recursive import EichlerRecursion


def test_eichler_recursion_rejects_the_conway_sloane_pair_at_index_one():
    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first, second = lattice(pair["form_a"]["gram"]), lattice(pair["form_b"]["gram"])
    assert not first.is_even()
    assert tuple(int(index) for index in first.signature_pair()) == (2, 1)

    assert EichlerRecursion().isometry(first, second) is None
