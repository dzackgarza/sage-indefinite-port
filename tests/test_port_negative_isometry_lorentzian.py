"""A same-genus non-isometric pair decided by the port's Lorentzian branch (#68).

The public Isom of the provider can answer from its spinor-genus prefilter
(sage-extensions lattice_morphisms.py, the indefinite rank >= 3 case). This test calls
EichlerRecursion.isometry directly. The Conway--Sloane forms (51a) and (51b) are odd, so
its even-genus check does not apply; the pair has signature (2, 1), so the decision is
made by the positive-index-1 route of _indefinite_isometry (lorentzian_isometry on the
signed view). Expected answer: Conway--Sloane, SPLAG 3rd ed., Ch. 15, section 11,
(51a) and (51b) lie in one genus and are not integrally equivalent.
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
