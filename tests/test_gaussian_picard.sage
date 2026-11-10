r"""Falbel-Francsics-Lax-Parker matrices and their reduction.

N_{(1+i,1)}^2 = N_{(2+2i,4)} by the Heisenberg law a_2 + b_2 + conj(a_1) b_1 = 1 + 1 + 2.
A word in the generators re-evaluates to the same matrix after reduction.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage_indefinite_port.indefinite.gaussian_picard import (
    CUSP_STABILIZER_GENERATORS, FIELD, FORM, GENERATORS, I, ONE, SCALAR, ZERO,
    _translation, evaluate, is_element, trace_isometry, word,
)


def test_four_source_matrices_lie_in_the_group():
    assert all(is_element(g) for g in GENERATORS.values())
    assert not is_element(SCALAR)
    assert SCALAR.det() == -I


def test_heisenberg_square_of_the_source_translation():
    two = FIELD(int(2))
    assert GENERATORS["N1"] ** int(2) == _translation(two + two * I, FIELD(int(4)))


def test_cusp_element_word():
    three = FIELD(int(3))
    element = _translation(ONE + I, ONE + three * I) * GENERATORS["M"]
    assert is_element(element)
    assert element[0, 2] == ZERO
    assert evaluate(word(element)) == element
    assert all(g[0, 2] == ZERO for g in CUSP_STABILIZER_GENERATORS)


def test_reduction_recovers_long_products():
    names = ("N1", "J", "N0", "N1", "J", "M", "N1", "J", "N0")
    product = evaluate([(name, int(1)) for name in names])
    assert product[0, 2] != ZERO
    reduced = word(product)
    assert evaluate(reduced) == product


def test_non_element_is_rejected_by_membership():
    half = _translation(ONE, ZERO)
    assert not is_element(half)


def test_trace_coordinates_are_order_linear_isometries():
    from sage_indefinite_port.indefinite.gaussian_picard import lattice
    hermitian = lattice()
    trace = hermitian.trace_lattice()
    assert int(trace.module_rank()) == int(6)
    assert tuple(int(x) for x in trace.signature_pair()) == (int(4), int(2))
    for g in GENERATORS.values():
        assert hermitian.is_order_linear(trace_isometry(g))


def test_generators_in_a_changed_basis_remain_order_linear_isometries():
    r"""P in GL_3(Z[i]) changes the O-basis: the form becomes P^* H P and g becomes P^{-1} g P.

    Every transformed generator is an integral isometry of the new rank-six trace lattice that
    commutes with multiplication by i.
    """
    from sage_extensions.categories.hermitian_lattices import HermitianLattice

    P = sage_matrix(FIELD, int(3), int(3), [ONE, ONE + I, ZERO, ZERO, ONE, ZERO, ZERO, I, ONE])
    assert P.det() in (ONE, -ONE, I, -I)
    changed = P.conjugate_transpose() * FORM * P
    hermitian = HermitianLattice((ONE, I), changed.transpose())
    trace = hermitian.trace_lattice()
    for g in GENERATORS.values():
        moved = ~P * g * P
        assert moved.conjugate_transpose() * changed * moved == changed
        isometry = hermitian.isometry_from_matrix(moved)
        assert hermitian.is_order_linear(isometry)
        assert all(trace.b(isometry(x), isometry(y)) == trace.b(x, y) for x in trace.module_generators() for y in trace.module_generators())
