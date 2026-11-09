r"""L = U + U with Gamma = Stab(e1). v_m = e2 + m f1 has b(e1, v_m) = m, so distinct m are
distinct orbits, an infinite family in one shell. The Eichler map E(e1, a), a = e2 + f2,
fixes e1 and sends v_1 to w = -2 e1 + f1 + 2 e2 + f2.
"""

from sage_extensions.all import *
from sage_indefinite_port.indefinite.pairing_slices import PointwiseStabilizerOrbits, line_transporter, pointwise_transporter, slice_invariant


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


def test_orbit_stream_on_the_plane_has_the_four_primitive_null_vectors():
    r"""Stab(e) in O(U) is trivial, so each of +-e, +-f is its own orbit."""
    L = Lattices(ZZ)("U")
    e, f = tuple(L.module_generators())
    orbits = iter(PointwiseStabilizerOrbits(L, (e,), 0))
    found = [next(orbits) for _ in range(int(4))]
    assert set(str(v) for v in found) == set(str(v) for v in (e, -e, f, -f))


def test_orbit_stream_on_two_planes_is_free_of_duplicates():
    r"""Pairings with e1 are invariants, and the first orbits realize several values."""
    L, e1, _f1, _e2, _f2 = _setup()
    orbits = PointwiseStabilizerOrbits(L, (e1,), 0)
    stream = iter(orbits)
    found = [next(stream) for _ in range(int(6))]
    assert all(v in orbits for v in found)
    assert all(not orbits.equivalent(found[i], found[j]) for i in range(int(6)) for j in range(i + int(1), int(6)))
    assert len(set(int(L.b(e1, v)) for v in found)) >= int(2)


def test_setwise_line_stabilizer_identifies_signs():
    r"""The setwise stabilizer of Ze in O(U) is {+-1}, so e ~ -e and e is not f."""
    L = Lattices(ZZ)("U")
    e, f = tuple(L.module_generators())
    assert line_transporter(L, e, f, -f) is not None
    assert pointwise_transporter(L, (e,), f, -f) is None
    assert line_transporter(L, e, e, f) is None


def test_tuple_slice_of_two_vectors():
    r"""On U + U, the swap of planes moves (e2, f2) to (e1, f1), but it does not fix e1."""
    L, e1, f1, e2, f2 = _setup()
    assert pointwise_transporter(L, (), (e2, f2), (e1, f1)) is not None
    assert pointwise_transporter(L, (e1,), (e2, f2), (e1, f1)) is None
