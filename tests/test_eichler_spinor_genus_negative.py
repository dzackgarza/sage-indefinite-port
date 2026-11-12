"""Same-genus, distinct-spinor-genus pair at signature (2, 2) decided without the covering search (#68).

L = diag(1, -2, -64, 128) has rational Witt index 2, so EichlerRecursion takes the
Eichler route. M is the second class of Genus(L). By Eichler-Kneser (Cassels,
Rational Quadratic Forms, Ch. 11, Thm. 1.3) an indefinite genus of rank at least 3
has one class per spinor genus, so L and M lie in different spinor genera and are
not isometric.

Expectation: the call returns None within a unit-test budget. Before the
spinor-genus decision in EichlerRecursion._indefinite_isometry, the call did not
finish: EichlerEnvelope.right_coset_transversal expanded more than 10^6 cosets
(10 h 15 min observed without a result).
"""

from tests.acceptance.consumer import lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion

FIRST = [[1, 0, 0, 0], [0, -2, 0, 0], [0, 0, -64, 0], [0, 0, 0, 128]]
SECOND = [[697, 0, 0, 896], [0, -2, 0, 0], [0, 0, -64, 0], [896, 0, 0, 1152]]


def test_distinct_spinor_genera_return_none():
    first, second = lattice(FIRST), lattice(SECOND)
    assert tuple(int(index) for index in first.signature_pair()) == (2, 2)
    assert first.genus() == second.genus()

    assert EichlerRecursion().isometry(first, second) is None
