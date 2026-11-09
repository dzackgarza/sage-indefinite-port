r"""Identities stated in section 7 of Allcock-Carlson-Toledo.

(7.4) adjacent R_i braid and nonadjacent commute on the affine E6 diagram.
(7.11) R_2^{-1} R_1 = T_{omega,0,0;theta/2}; [T, T'] = T_{0,0,0;theta}.
(7.8.1) composition law of translations.
"""

from sage_extensions.all import *
from sage.modules.free_module_element import vector as sage_vector
from sage_indefinite_port.indefinite.act_eisenstein import (
    AUT_PLUS_GENERATORS, BAR, BIFLECTION, FIELD, FORM, FULL_GENERATORS, OMEGA, ONE, R, RHO, ROOTS, THETA, ZERO,
    h, null_transporter, reduce_short_root, root_transporter, is_aut_plus, is_gamma_theta, trace_generators, triflection, is_isometry, is_primitive, lattice, reduce_null_vector, spinor_character, translation,
)

EDGES = ((1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (6, 7))


def _t(*values):
    return sage_vector(FIELD, [v if hasattr(v, 'vector') else FIELD(int(v)) for v in values])


def test_roots_and_generators():
    assert all(h(r, r) == ONE for r in ROOTS.values())
    assert all(is_isometry(g) for g in FULL_GENERATORS)
    assert len(AUT_PLUS_GENERATORS) == int(7)
    assert R[1] * ROOTS[1] == -OMEGA * ROOTS[1]


def test_affine_e6_braid_relations():
    edges = {(int(a), int(b)) for a, b in EDGES}
    labels = sorted(int(i) for i in R)
    for i in labels:
        for j in labels:
            if i < j:
                if (i, j) in edges:
                    assert R[i] * R[j] * R[i] == R[j] * R[i] * R[j]
                else:
                    assert R[i] * R[j] == R[j] * R[i]


def test_source_translations():
    half = THETA / FIELD(int(2))
    assert ~R[2] * R[1] == translation(_t(OMEGA, 0, 0), half)
    first = translation(_t(OMEGA, 0, 0), half)
    second = translation(_t(-BAR, 0, 0), half)
    assert first * second * ~first * ~second == translation(_t(0, 0, 0), THETA)


def test_translation_composition_law():
    half = THETA / FIELD(int(2))
    a, b = _t(OMEGA, 0, 0), _t(0, ONE, 0)
    composed = translation(a, half) * translation(b, half)
    pairing = b.dot_product(a.conjugate())
    assert composed == translation(a + b, half + half + (pairing - pairing.conjugate()) / FIELD(int(2)))


def test_null_reduction_from_source_shape():
    v = _t(THETA, 0, 0, 1, OMEGA - ONE)
    assert h(v, v) == ZERO
    g = reduce_null_vector(v)
    assert is_isometry(g)
    assert g * v == RHO


def test_reduction_of_a_group_image_of_rho():
    g = R[3] * R[6] * translation(_t(ONE, OMEGA, 0), ZERO) * R[2] * R[3]
    v = g * RHO
    assert is_primitive(v)
    assert v[3] != ZERO
    back = reduce_null_vector(v)
    assert is_isometry(back)
    assert back * v == RHO


def test_orthogonal_height_one_root_case():
    r"""v = (1,1,1;theta,1-omega): Sum |v_i|^2 = 3 = |H|^2 and r = (0,0,0;1,-omega) is orthogonal to it."""
    v = _t(ONE, ONE, ONE, THETA, ONE - OMEGA)
    assert h(v, v) == ZERO
    assert h(v, ROOTS[3]) == ZERO
    g = reduce_null_vector(v)
    assert is_isometry(g) and g * v == RHO


def test_unit_multiples_of_rho_reach_rho():
    for unit in (-ONE, OMEGA, -OMEGA, BAR, -BAR):
        g = reduce_null_vector(unit * RHO)
        assert is_isometry(g) and g * (unit * RHO) == RHO


def test_eisenstein_primitivity_is_not_integer_primitivity():
    assert not is_primitive(THETA * RHO)


def test_trace_lattice_rank_ten_signature():
    trace = lattice().trace_lattice()
    assert int(trace.module_rank()) == int(10)
    assert tuple(int(x) for x in trace.signature_pair()) == (int(8), int(2))


def test_spinor_character_separates_aut_plus():
    r"""Hexflections reduce to reflections in vectors of norm 1, a square; B to one of norm 2 = -1, a nonsquare."""
    assert all(is_aut_plus(g) for g in AUT_PLUS_GENERATORS)
    assert spinor_character(BIFLECTION) == int(-1)
    assert not is_aut_plus(BIFLECTION)


def test_spinor_character_is_multiplicative():
    elements = (R[1], R[3], R[6], BIFLECTION, translation(_t(ONE, OMEGA, 0), ZERO))
    assert all(
        spinor_character(a * b) == spinor_character(a) * spinor_character(b)
        for a in elements for b in elements
    )


def test_generators_are_trace_isometries_commuting_with_omega():
    hermitian = lattice()
    action = hermitian.order_action(OMEGA)
    trace = hermitian.trace_lattice()
    for isometry in trace_generators():
        assert all(isometry(action(p)) == action(isometry(p)) for p in trace.module_generators())


def test_gamma_theta_contains_triflections_and_not_hexflections():
    r"""omega = 1 modulo theta, so an omega-reflection is trivial on V; -omega acts as -1 on the root."""
    assert all(is_gamma_theta(triflection(r)) for r in ROOTS.values())
    assert all(triflection(r) == R[i] ** int(4) for i, r in ROOTS.items())
    assert not any(is_gamma_theta(R[i]) for i in R)


def test_null_transporter_carries_source_vectors():
    v = _t(ONE, ONE, ONE, THETA, ONE - OMEGA)
    w = _t(THETA, 0, 0, 1, OMEGA - ONE)
    g = null_transporter(v, w)
    assert is_isometry(g) and g * v == w and is_aut_plus(g)


def test_short_roots_reduce_to_r7():
    for root in (ROOTS[1], ROOTS[2], ROOTS[3], ROOTS[4], ROOTS[6], -ROOTS[5], OMEGA * ROOTS[1]):
        g = reduce_short_root(root)
        assert is_isometry(g) and g * root == ROOTS[7]


def test_root_transporter_on_a_group_image():
    start = ROOTS[3]
    target = R[1] * R[2] * R[3] * R[6] * start
    g = root_transporter(start, target)
    assert is_isometry(g) and g * start == target


def test_nonroot_is_rejected():
    from pytest import raises

    with raises(ValueError):
        reduce_short_root(RHO)


def test_public_centralizer_of_omega_is_generated_by_the_source_group():
    r"""U(Lambda, h) is the centralizer of omega in O(L), and (7.23) gives its generators R_i and B."""
    trace = lattice().trace_lattice()
    omega = lattice().unit_isometry(OMEGA)
    generators = tuple(trace.O().centralizer(omega).group_generators())
    assert len(generators) == len(FULL_GENERATORS)
    assert all(g * omega == omega * g for g in generators)
