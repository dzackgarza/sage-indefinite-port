"""Direct Allcock Promised branch regression (arXiv:2110.03784, Algorithm Promised)."""

import pytest

from sage_indefinite_port.indefinite.edgewalk_rank2 import _OrientedBinaryCycle


def test_promised_go_left_replaces_both_search_vectors():
    # With Gram [[2,3],[3,-1]], the first middle m=(1,1) has
    # q(m)=7>1, b(m,l)=2>=0 and q(l)=-1. Thus the source takes
    # [Go left], which recurses on (m,r), not (m,l).
    cycle = object.__new__(_OrientedBinaryCycle)
    cycle.form = (2, 3, -1)
    pairing = cycle.pairing
    seen = []

    def check(left, right):
        seen.append((left, right))
        if len(seen) == 4:
            raise StopIteration
        return pairing(left, right)

    cycle.pairing = check
    with pytest.raises(StopIteration):
        cycle.promised_step(1, (1, 0), (0, 1))
    # On the second iteration, [Find middle] evaluates q((m+r)).
    assert seen[3] == ((2, 1), (2, 1))
