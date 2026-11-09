r"""Identities stated in section 7 of Allcock-Carlson-Toledo.

(7.4) adjacent R_i braid and nonadjacent commute on the affine E6 diagram.
(7.11) R_2^{-1} R_1 = T_{omega,0,0;theta/2}; [T, T'] = T_{0,0,0;theta}.
(7.8.1) composition law of translations.
"""

from sage_extensions.all import *
from sage.modules.free_module_element import vector as sage_vector
from sage_indefinite_port.indefinite.act_eisenstein import (
    AUT_PLUS_GENERATORS, BAR, BIFLECTION, FIELD, FORM, FULL_GENERATORS, OMEGA, ONE, R, RHO, ROOTS, THETA, ZERO,
    h, is_isometry, is_primitive, lattice, reduce_null_vector, translation,
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
    g, unit = reduce_null_vector(v)
    assert is_isometry(g)
    assert g * v == unit * RHO


def test_reduction_of_a_group_image_of_rho():
    g = R[3] * R[6] * translation(_t(ONE, OMEGA, 0), ZERO) * R[2] * R[3]
    v = g * RHO
    assert is_primitive(v)
    assert v[3] != ZERO
    back, unit = reduce_null_vector(v)
    assert back * v == unit * RHO


def test_eisenstein_primitivity_is_not_integer_primitivity():
    assert not is_primitive(THETA * RHO)


def test_trace_lattice_rank_ten_signature():
    trace = lattice().trace_lattice()
    assert int(trace.module_rank()) == int(10)
    assert tuple(int(x) for x in trace.signature_pair()) == (int(8), int(2))
