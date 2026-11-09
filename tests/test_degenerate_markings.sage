r"""Stabilizers of a line in L = <2> + 0, x generating <2> and r the radical.

O(L) is (x, r) -> (e x + b r, d r). The setwise stabilizer of Zx has b = 0, of order 4; the
pointwise one has also e = 1, of order 2. For S = Z(x + r), g(x + r) = e x + (b + d) r must be
+-(x + r), so e = 1 gives b = 1 - d and e = -1 gives b = -1 - d: order 4 again, and it includes
(e, b, d) = (1, 2, -1).
"""

from sage_extensions.all import *
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.integer_ring import ZZ as SageZZ
from sage_indefinite_port.indefinite.degenerate_markings import degenerate_marking_generators


def _closure(generators):
    one = sage_identity(SageZZ, int(2))
    known = {one.str(): one}
    frontier = [one]
    while frontier:
        current = frontier.pop()
        for g in generators:
            image = g * current
            if image.str() not in known:
                known[image.str()] = image
                frontier.append(image)
    return known


def _setup():
    L = Lattices(ZZ)([[2, 0], [0, 0]])
    x, r = tuple(L.module_generators())
    return L, x, r


def test_setwise_and_pointwise_stabilizers_of_x():
    L, x, _r = _setup()
    assert len(_closure(degenerate_marking_generators(L, (x,), False))) == int(4)
    assert len(_closure(degenerate_marking_generators(L, (x,), True))) == int(2)


def test_setwise_stabilizer_of_a_tilted_line():
    L, x, r = _setup()
    group = _closure(degenerate_marking_generators(L, (x + r,), False))
    assert len(group) == int(4)
    assert all(g * (sage_identity(SageZZ, int(2)).column(int(0)) + sage_identity(SageZZ, int(2)).column(int(1))) in (sage_identity(SageZZ, int(2)).column(int(0)) + sage_identity(SageZZ, int(2)).column(int(1)), -(sage_identity(SageZZ, int(2)).column(int(0)) + sage_identity(SageZZ, int(2)).column(int(1)))) for g in group.values())


def test_public_stabilizers_on_a_degenerate_lattice():
    L, x, r = _setup()
    embedding = L.subobject_on((x,)).inclusion()
    setwise = tuple(L.O().setwise_stabilizer(embedding).group_generators())
    pointwise = tuple(L.O().pointwise_stabilizer(embedding).group_generators())
    assert setwise and pointwise
    assert all(g(x) in (x, -x) and g(r) in (r, -r) for g in setwise)
    assert all(g(x) == x for g in pointwise)
