r"""L = U + U with Gamma = Stab(e1). v_m = e2 + m f1 has b(e1, v_m) = m, so distinct m are
distinct orbits, an infinite family in one shell. The Eichler map E(e1, a), a = e2 + f2,
fixes e1 and sends v_1 to w = -2 e1 + f1 + 2 e2 + f2.
"""

from sage_extensions.all import *
from sage_indefinite_port.indefinite.pairing_slices import pointwise_transporter, slice_invariant


def _setup():
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e1, f1, e2, f2 = tuple(L.module_generators())
    return L, e1, f1, e2, f2


def test_infinite_family_is_separated_by_the_pairing():
    L, e1, f1, e2, _f2 = _setup()
    family = [e2 + L.scalar_multiple(m, f1) for m in range(int(6))]
    assert len(set(int(slice_invariant(L, (e1,), v)[0][0]) for v in family)) == int(6)
    assert all(L.q(v) == 0 for v in family)
    assert pointwise_transporter(L, (e1,), family[1], family[2]) is None


def test_transporter_in_the_stabilizer():
    L, e1, f1, e2, f2 = _setup()
    v = e2 + f1
    w = L.scalar_multiple(-2, e1) + f1 + L.scalar_multiple(2, e2) + f2
    assert L.q(w) == 0
    g = pointwise_transporter(L, (e1,), v, w)
    assert g is not None and g(e1) == e1 and g(v) == w
