r"""L = U + U with Gamma = Stab(e1). v_m = e2 + m f1 has b(e1, v_m) = m, so distinct m are
distinct orbits, an infinite family in one shell. The Eichler map E(e1, a), a = e2 + f2,
fixes e1 and sends v_1 to w = -2 e1 + f1 + 2 e2 + f2.
"""

from sage_extensions.all import *
from sage_indefinite_port.indefinite.pairing_slices import PointwiseStabilizerOrbits, SubgroupVectorOrbits, line_transporter, setwise_transporter, pointwise_transporter, slice_invariant


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


def test_setwise_on_a_marking_with_a_radical():
    r"""U + A1 + A1, S = A1 + A1, marking (x1, x2, e + x1), which has radical e. The swap sends it to (x2, x1, e + x2)."""
    L = Lattices(ZZ)("U") + Lattices(ZZ)("A1") + Lattices(ZZ)("A1")
    e, _f, x1, x2 = tuple(L.module_generators())
    g = setwise_transporter(L, (x1, x2), e + x1, e + x2)
    assert g is not None and g(e + x1) == e + x2


def test_setwise_rank_two_stabilizer_exchanges_summands():
    r"""L = U + A1 + A1, S = A1 + A1. The swap preserves S and sends h + x1 to h + x2, h = e + f.
    Pointwise, the pairing with x1 differs. The vector h pairs to 0 with S, so it is not h + x1."""
    L = Lattices(ZZ)("U") + Lattices(ZZ)("A1") + Lattices(ZZ)("A1")
    e, f, x1, x2 = tuple(L.module_generators())
    h = e + f
    left, right = h + x1, h + x2
    assert pointwise_transporter(L, (x1, x2), left, right) is None
    g = setwise_transporter(L, (x1, x2), left, right)
    assert g is not None and g(left) == right
    assert set(str(g(x)) for x in (x1, x2)) <= set(str(v) for v in (x1, x2, -x1, -x2))
    assert setwise_transporter(L, (x1, x2), left, h) is None


def test_mixed_marking_is_decided_in_the_complement():
    r"""L = U + A1, x the generator of A1. (x, e + x) spans a marking with radical e. O(U) swaps e, f,
    so (x, e + x) moves to (x, f + x). The tuple (x, x + x) has a different Gram matrix.
    """
    L = Lattices(ZZ)("U") + Lattices(ZZ)("A1")
    e, f, x = tuple(L.module_generators())
    g = pointwise_transporter(L, (x,), e + x, f + x)
    assert g is not None and g(x) == x and g(e + x) == f + x
    assert pointwise_transporter(L, (x,), e + x, x + x) is None


def test_congruence_subgroup_orbit_stream():
    r"""O(U) is faithful modulo three, so Gamma(3) is trivial and each primitive null vector is its own orbit."""
    L = Lattices(ZZ)("U")
    e, f = tuple(L.module_generators())
    gamma = L.O().congruence_preimage(3)
    stream = iter(SubgroupVectorOrbits(L, gamma, 0))
    found = [next(stream) for _ in range(int(4))]
    assert set(str(v) for v in found) == set(str(v) for v in (e, -e, f, -f))


def test_ordered_and_nonprimitive_null_tuples_on_two_planes():
    r"""Issue 61, items 2 and 3. U + U swaps the planes, so (e1, e2) goes to (e2, e1). The tuple (2 e1, e2)
    spans a sublattice of index two in its saturation and (e1, e2) has index one, so no isometry relates them."""
    from sage_indefinite_port.indefinite.pairing_slices import marked_transporter

    L, e1, f1, e2, f2 = _setup()
    g = marked_transporter(L, (e1, e2), (e2, e1))
    assert g is not None and g(e1) == e2 and g(e2) == e1
    double = L.scalar_multiple(2, e1)
    assert marked_transporter(L, (double, e2), (e1, e2)) is None
    assert marked_transporter(L, (double, e2), (double, e2)) is not None


def test_equal_gram_with_only_a_rational_transporter():
    r"""Issue 61, item 9. (e1) and (2 e1) both have Gram (0), and no isometry maps a
    primitive vector to a nonprimitive one."""
    from sage_indefinite_port.indefinite.pairing_slices import marked_transporter

    L, e1, _f1, _e2, _f2 = _setup()
    assert marked_transporter(L, (e1,), (L.scalar_multiple(2, e1),)) is None


def test_setwise_and_based_stabilizers_of_a_null_plane():
    r"""Issue 61, item 1. For I = <e1, e2> in U + U, based elements fix e1 and e2; setwise ones preserve I."""
    L, e1, f1, e2, f2 = _setup()
    embedding = L.primitive_isotropic_subobject(e1, e2).inclusion()
    based = tuple(L.O().pointwise_stabilizer(embedding).group_generators())
    setwise = tuple(L.O().setwise_stabilizer(embedding).group_generators())
    assert based and setwise
    assert all(g(e1) == e1 and g(e2) == e2 for g in based)
    for g in setwise:
        for v in (e1, e2):
            assert int(L.b(g(v), e1)) == int(0) and int(L.b(g(v), e2)) == int(0)
    assert any(g(e1) != e1 or g(e2) != e2 for g in setwise)


def test_congruence_transporter_is_in_the_subgroup():
    r"""Issue 61, item 8. A transporter under Gamma(3) reduces to the identity modulo three."""
    L, e1, f1, _e2, _f2 = _setup()
    gamma = L.O().congruence_preimage(3)
    g = gamma.vector_equivalence_witness(e1, e1 + L.scalar_multiple(3, _e2))
    assert g is not None and g in gamma and g(e1) == e1 + L.scalar_multiple(3, _e2)


def test_skipping_flag_tuples_in_three_planes():
    r"""Issue 61, item 6. In U + U + U, the flag e1 in <e1, e2, e3> is the ordered tuple (e1, e2, e3) up to the
    upper triangular source group, and the transporter carries the tuple and hence each face."""
    from sage_indefinite_port.indefinite.pairing_slices import marked_transporter

    L = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e1, f1, e2, f2, e3, f3 = tuple(L.module_generators())
    g = marked_transporter(L, (e1, e2, e3), (e2, e1, e3))
    assert g is not None and (g(e1), g(e2), g(e3)) == (e2, e1, e3)
    assert marked_transporter(L, (e1, e2, e3), (e1, e2, e2 + e3)) is not None


def test_public_orbit_set_interface():
    r"""O(U) has one orbit of primitive null vectors, and Gamma(3) = 1 has four."""
    L = Lattices(ZZ)("U")
    e, f = tuple(L.module_generators())
    full = L.O().vector_orbit_set(0)
    assert full.equivalent(e, -f)
    level = L.O().vector_orbit_set(0, L.O().congruence_preimage(3))
    stream = iter(level)
    assert len([next(stream) for _ in range(int(4))]) == int(4)
    assert not level.equivalent(e, f)


def test_setwise_indefinite_sublattice_respects_the_glue():
    r"""L = U, S = Z(e + 3f), K = S^perp = Z(e - 3f), glued along Z/6.

    O(U) = {1, -1, s, -s}, with s the swap. The swap sends e + 3f to 3e + f, so the setwise
    stabilizer of S is {1, -1} and the orbit of e is {e, -e}. The vector 3f has the same
    projections up to the sign of the K part, so only the glue rejects it: (1, -1) on S + K
    acts by different signs on the two copies of Z/6.
    """
    L = Lattices(ZZ)("U")
    e, f = tuple(L.module_generators())
    s = e + f + f + f
    g = setwise_transporter(L, (s,), e, -e)
    assert g is not None and g(e) == -e and g(s) == -s
    assert setwise_transporter(L, (s,), e, f + f + f) is None


def test_setwise_indefinite_unimodular_summand():
    r"""L = U + U, S the first U. Its setwise stabilizer is O(U) x O(U).

    The swap of e1 and f1 moves e1 + e2 to f1 + e2. The square of the projection to S is
    invariant, and it is 2 for e1 + f1 + e2 and 0 for e1 + e2 + f2.
    """
    L, e1, f1, e2, f2 = _setup()
    g = setwise_transporter(L, (e1, f1), e1 + e2, f1 + e2)
    assert g is not None and g(e1 + e2) == f1 + e2
    assert set(str(g(x)) for x in (e1, f1)) <= set(str(v) for v in (e1, f1, -e1, -f1))
    assert setwise_transporter(L, (e1, f1), e1 + f1 + e2, e1 + e2 + f2) is None


def _pell_plus_two():
    r"""L = <1> + <-2> + <2> with basis (e1, e2, e3) and f = (3 + 2 sqrt 2) + 1, of infinite order.

    C(f) = (+-f^Z) x {+-1}: the eigenspaces of f split L with no glue, and the centralizer of
    multiplication by a unit of Z[sqrt 2] in O(<1> + <-2>) is the norm-one group +-(3 + 2 sqrt 2)^Z.
    """
    from sage_indefinite_port.indefinite.pairing_slices import CentralizerSubgroup

    L = Lattices(ZZ)([[1, 0, 0], [0, -2, 0], [0, 0, 2]])
    e1, e2, e3 = tuple(L.module_generators())
    labels = tuple(L.module_generating_set())
    images = (L.scalar_multiple(3, e1) + L.scalar_multiple(2, e2), L.scalar_multiple(4, e1) + L.scalar_multiple(3, e2), e3)
    f = L.O()(dict(zip(labels, images, strict=True)))
    return L, e1, e2, e3, f, CentralizerSubgroup(L, (f,))


def test_centralizer_vector_witness_and_obstruction():
    r"""f(e1) = 3 e1 + 2 e2 is in the orbit of e1 (witness f). v = e1 + e2 + e3 has norm 1 but is not:
    g(e1) = v forces g(f e1) = f v, and the saturated span of e1, f e1 is <e1, e2>, whose image would be
    spanned by v and 2 e1 + e2 - e3 of norm 4, not -2."""
    L, e1, e2, e3, f, group = _pell_plus_two()
    w = f(e1)
    g = group.vector_equivalence_witness(e1, w)
    assert g is not None and g(e1) == w and g in group
    v = e1 + e2 + e3
    assert L.q(v) == L.q(e1)
    assert group.vector_equivalence_witness(e1, v) is None
    assert group.vector_equivalence_witness(e1, -e1) is not None


def test_centralizer_orbit_stream_separates_the_glued_shell():
    r"""Norm-one vectors with coordinates in {-1, 0, 1}: (+-1, 0, 0), one orbit, and (+-1, +-1, +-1), one orbit,
    since x + y sqrt 2 of norm -1 is +- (1 + sqrt 2) (3 + 2 sqrt 2)^n and e3 -> -e3 is in C(f)."""
    from sage_indefinite_port.indefinite.pairing_slices import SubgroupVectorOrbits

    L, e1, e2, e3, f, group = _pell_plus_two()
    stream = iter(SubgroupVectorOrbits(L, group, 1))
    first, second = next(stream), next(stream)
    assert set(int(abs(c)) for c in first) == set([int(1)])
    assert list(second)[int(1):] == [int(0), int(0)]
    assert group.vector_equivalence_witness(first, e1 + e2 + e3) is not None


def test_centralizer_stabilizers_of_the_fixed_line():
    r"""Stab(e3) = (+-f^Z) x 1, infinite. The setwise stabilizer of Z e3 adds 1 + (-1)."""
    L, e1, e2, e3, f, group = _pell_plus_two()
    pointwise = group.vector_stabilizer_generators(e3)
    assert all(g(e3) == e3 and g in group for g in pointwise)
    assert any(g(e1) not in (e1, -e1) for g in pointwise)
    setwise = group.setwise_stabilizer_generators((e3,))
    assert all(g(e3) in (e3, -e3) and g in group for g in setwise)
    assert any(g(e3) == -e3 for g in setwise)


def test_centralizer_setwise_stabilizer_of_a_nonprimitive_invariant_sublattice():
    r"""Lambda = <3 e1, e2> has Lambda_Q = <e1, e2>, which f preserves. f = [[3, 4], [2, 3]] is [[0, 1], [2, 0]]
    modulo 3, so f(e2) = 4 e1 + 3 e2 leaves Lambda, while f^2 = -1 modulo 3 preserves it. Hence
    Stab(Lambda) = (+-f^{2Z}) x {+-1 on e3}, of index 2 in C(f), and it contains elements of infinite order."""
    L, e1, e2, e3, f, group = _pell_plus_two()
    basis = (L.scalar_multiple(3, e1), e2)
    generators = group.setwise_stabilizer_generators(basis)
    from sage.matrix.constructor import matrix as sage_matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    span = lambda vectors: sage_matrix(SageZZ, [[int(c) for c in v] for v in vectors]).row_module()
    assert span([f(v) for v in basis]) != span(basis)
    assert span([f(f(v)) for v in basis]) == span(basis)
    assert all(g in group and span([g(v) for v in basis]) == span(basis) for g in generators)
    assert any(g(e3) == -e3 for g in generators)
    assert any(not all(p.is_cyclotomic() for p, _e in group._matrix(g).charpoly().factor()) for g in generators)


def _swap_of_copies(gram):
    r"""L = N + N for the Gram matrix ``gram`` of N, with F the swap of the copies, in the centralizer subgroup."""
    from sage_indefinite_port.indefinite.pairing_slices import CentralizerSubgroup

    n = len(gram)
    rows = [[int(gram[i][j]) if i < n and j < n else int(gram[i - n][j - n]) if i >= n and j >= n else int(0)
             for j in range(int(2) * n)] for i in range(int(2) * n)]
    L = Lattices(ZZ)(rows)
    frame = tuple(L.module_generators())
    labels = tuple(L.module_generating_set())
    images = tuple(frame[(i + n) % (int(2) * n)] for i in range(int(2) * n))
    F = L.O()(dict(zip(labels, images, strict=True)))
    return L, frame, F, CentralizerSubgroup(L, (F,))


def test_centralizer_setwise_stabilizer_of_a_summand_moved_by_the_swap():
    r"""L = U + U, F the swap of the copies, Lambda = U + 0, indefinite and not F-invariant (F Lambda = 0 + U).

    g in C(F) preserving Lambda also preserves F Lambda = 0 + U, so g = P + Q, and gF = Fg forces Q = P.
    Hence Stab(Lambda) = {P + P : P in O(U)} = {1, -1, s + s, -s - s}, s the swap of e and f, of order 4."""
    from sage.matrix.constructor import matrix as sage_matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    L, frame, F, group = _swap_of_copies([[0, 1], [1, 0]])
    generators = [group._matrix(g) for g in group.setwise_stabilizer_generators(frame[:int(2)])]
    one = sage_matrix.identity(SageZZ, int(4))
    seen = {one.str(): one}
    frontier = [one]
    while frontier:
        current = frontier.pop()
        for g in generators:
            image = g * current
            if image.str() not in seen:
                seen[image.str()] = image
                frontier.append(image)
    s = sage_matrix(SageZZ, [[0, 1], [1, 0]])
    two = sage_matrix.identity(SageZZ, int(2))
    expected = [sage_matrix.block_diagonal([P, P]) for P in (two, -two, s, -s)]
    assert set(seen) == set(m.str() for m in expected)


def test_centralizer_setwise_stabilizer_of_a_moved_pell_summand_is_infinite():
    r"""L = P + P, P = <1> + <-2>, F the swap and Lambda = P + 0. As above Stab(Lambda) = {h + h : h in O(P)}, and
    O(P) contains -1 and f = 3 + 2 sqrt 2 of infinite order."""
    L, frame, F, group = _swap_of_copies([[1, 0], [0, -2]])
    generators = [group._matrix(g) for g in group.setwise_stabilizer_generators(frame[:int(2)])]
    for g in generators:
        assert g.submatrix(int(0), int(2), int(2), int(2)).is_zero() and g.submatrix(int(2), int(0), int(2), int(2)).is_zero()
        assert g.submatrix(int(0), int(0), int(2), int(2)) == g.submatrix(int(2), int(2), int(2), int(2))
    assert any(not all(p.is_cyclotomic() for p, _e in g.charpoly().factor()) for g in generators)
    assert any(g[int(0), int(0)] < 0 for g in generators)
