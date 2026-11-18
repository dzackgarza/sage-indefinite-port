r"""The order ZZ<1,i,j,k> of (-1,3)_QQ with q = 2 Nrd has Gram diag(2,2,-6,-6).

Nrd(2+j) = 4-3 = 1. By hand: (2+j)1 = 2+j, (2+j)i = 2i-k, (2+j)j = 3+2j and
(2+j)k = -3i+2k, while on the right i(2+j) = 2i+k and k(2+j) = 3i+2k.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.quaternion_lattices import QuaternionNormLattice


def _model():
    return QuaternionNormLattice(int(-1), int(3), int(2))


def _unit(model):
    return model.element(int(2), int(0), int(1), int(0))


def test_the_norm_lattice_has_the_expected_gram_matrix():

    assert [[int(_model().lattice().gram_matrix()[i, j]) for j in range(int(4))] for i in range(int(4))] == [[2 if i == j and i < 2 else -6 if i == j else 0 for j in range(4)] for i in range(4)]


def test_left_multiplication_by_two_plus_j_has_the_hand_matrix():
    model = _model()
    expected = sage_matrix(SageQQ, [[int(x) for x in row] for row in [[2, 0, 3, 0], [0, 2, 0, -3], [1, 0, 2, 0], [0, -1, 0, 2]]])
    assert model.left_multiplication_matrix(_unit(model)) == expected


def test_right_multiplication_differs_from_left_on_i_and_k():
    model = _model()
    expected = sage_matrix(SageQQ, [[int(x) for x in row] for row in [[2, 0, 3, 0], [0, 2, 0, 3], [1, 0, 2, 0], [0, 1, 0, 2]]])
    assert model.right_multiplication_matrix(_unit(model)) == expected


def test_the_isometry_preserves_the_form_and_inverts_through_the_conjugate():
    model = _model()
    unit = _unit(model)
    forward = model.isometry(unit, "left")
    backward = model.isometry(unit.conjugate(), "left")
    points = tuple(model.lattice().module_generators())
    for x in points:
        assert backward(forward(x)) == x
        for y in points:
            assert model.lattice().b(forward(x), forward(y)) == model.lattice().b(x, y)
    one = points[0]
    assert forward(one) == int(2) * points[0] + points[2]


def test_conjugation_is_an_isometry_outside_every_pair_action():
    r"""Conjugation fixes 1 and negates i, j, k, so it has determinant -1; a pair action has determinant 1."""
    model = _model()
    unit = _unit(model)
    assert model.conjugation_matrix().det() == int(-1)
    assert model.pair_matrix(unit, unit).det() == int(1)
    image = model.conjugation_isometry()
    points = tuple(model.lattice().module_generators())
    assert image(points[0]) == points[0] and image(points[1]) == -points[1]


def test_a_pair_of_norm_minus_one_units_has_nontrivial_spinor_class():
    r"""Nrd(1+i+j) = 1+1-3 = -1, so the pair is an isometry with spinor class -1, unlike norm-one pairs."""
    model = _model()
    u = model.element(int(1), int(1), int(1), int(0))
    forward = model.pair_isometry(u, u)
    points = tuple(model.lattice().module_generators())
    for x in points:
        for y in points:
            assert model.lattice().b(forward(x), forward(y)) == model.lattice().b(x, y)
    assert model.pair_spinor_class(u) == int(-1)
    assert model.pair_spinor_class(_unit(model)) == int(1)


def test_the_normalizer_quotient_has_exactly_twelve_classes():
    r"""N(Lambda) / QQ^x Lambda^x embeds in the local groups, as QQ has class number one.

    At 3 the algebra is ramified and Lambda is maximal, giving order 2. At 2 it is
    ramified with Lambda = Z_2 + P, and A^x / Q^x Lambda^x has order 2 * |F_4^x / F_2^x| = 6.
    So there are at most 12 classes, and twelve distinct classes in a box prove equality.
    """
    assert len(_model().normalizer_classes_in_box(int(3))) == int(12)


def test_conjugation_by_one_plus_i_exchanges_j_and_k():
    r"""(1+i) j (1+i)^{-1} = (j+k)(1-i)/2 = k, since ji = -k and ki = j; likewise k goes to -j."""
    model = _model()
    element = model.element(int(1), int(1), int(0), int(0))
    points = tuple(model.lattice().module_generators())
    isometry = model.conjugation_by_isometry(element)
    assert isometry(points[0]) == points[0] and isometry(points[1]) == points[1]
    assert isometry(points[2]) == points[3] and isometry(points[3]) == -points[2]

