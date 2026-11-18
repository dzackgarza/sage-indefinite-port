r"""Stabilizers of a line in L = <2> + 0, x generating <2> and r the radical.

O(L) is (x, r) -> (e x + b r, d r). The setwise stabilizer of Zx has b = 0, of order 4; the
pointwise one has also e = 1, of order 2. For S = Z(x + r), g(x + r) = e x + (b + d) r must be
+-(x + r), so e = 1 gives b = 1 - d and e = -1 gives b = -1 - d: order 4 again, and it includes
(e, b, d) = (1, 2, -1).
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.modules.free_module_element import vector as sage_vector
from sage.rings.finite_rings.finite_field_constructor import GF as SageGF
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


def _order_mod_three(generators, size):
    field = SageGF(int(3))
    one = sage_identity(field, int(size))
    one.set_immutable()
    known = set([one])
    frontier = [one]
    reduced = [g.change_ring(field) for g in generators]
    while frontier:
        current = frontier.pop()
        for g in reduced:
            image = current * g
            image.set_immutable()
            if image not in known:
                known.add(image)
                frontier.append(image)
    return len(known)


def test_marking_inside_the_radical():
    r"""L = <2> + 0 + 0, S = Z r1. Modulo 3, (e, b, D): e = +-1, b in F_3^2, D upper triangular.
    Setwise: D diag +-1 with a shear gives 12, so the order is 2 * 9 * 12 = 216.
    Pointwise: the corner is 1 and the other diagonal is +-1, so 2 * 9 * 6 = 108."""
    L = Lattices(ZZ)([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    x, r1, r2 = tuple(L.module_generators())
    assert _order_mod_three(degenerate_marking_generators(L, (r1,), False), 3) == int(216)
    assert _order_mod_three(degenerate_marking_generators(L, (r1,), True), 3) == int(108)


def test_marking_meeting_the_radical_and_the_quotient():
    r"""S = <r1, x + r2> in the same L: every generator preserves S, and pointwise ones fix it."""
    L = Lattices(ZZ)([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    x, r1, r2 = tuple(L.module_generators())
    span = sage_matrix(SageZZ, [[0, 1, 0], [1, 0, 1]]).row_module()
    for g in degenerate_marking_generators(L, (r1, x + r2), False):
        assert all(sage_vector(SageZZ, list(g * v)) in span for v in span.basis())
    fixed = [sage_vector(SageZZ, [0, 1, 0]), sage_vector(SageZZ, [1, 0, 1])]
    for g in degenerate_marking_generators(L, (r1, x + r2), True):
        assert all(g * v == v for v in fixed)
