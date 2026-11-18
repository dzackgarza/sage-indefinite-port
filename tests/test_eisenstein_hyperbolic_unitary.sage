r"""U(H + I_m) over Z[omega] for m = 1, 2 (#51).

The reflection in the norm-one vector r = e + b_1 + theta f with root of unity -omega,
x -> x - (1 + omega) h(x, r) r, is integral and unitary, and h(e, r) = conj(theta) is nonzero,
so it moves e and its word needs the swap.
"""

from sage_extensions.all import *
from sage.modules.free_module_element import vector as sage_vector
from sage_indefinite_port.indefinite.eisenstein_hyperbolic_unitary import (
    FIELD, OMEGA, ONE, THETA, ZERO, evaluate, form, generators, is_element, word,
)


def _reflection(m, root, xi):
    J = form(m)
    norm = (root.conjugate() * J * root)
    columns = []
    for e in J.parent().identity_matrix().columns():
        pairing = root.conjugate() * J * e
        columns.append(e - (ONE - xi) * pairing / norm * root)
    from sage.matrix.constructor import matrix
    return matrix(FIELD, columns).transpose()


def test_generators_are_integral_unitary():
    for m in (int(1), int(2)):
        assert all(is_element(m, g) for g in generators(m).values())


def test_words_of_products_of_generators():
    for m in (int(1), int(2)):
        table = generators(m)
        g = table["W"] * table["T0"] * table["W"] * table["Tw0"] ** int(2) * table["D"] * table["W"]
        assert evaluate(m, word(m, g)) == g


def test_word_of_a_reflection_that_moves_e():
    for m in (int(1), int(2)):
        root = sage_vector(FIELD, [ONE, ONE] + [ZERO] * (m - int(1)) + [THETA])
        g = _reflection(m, root, -OMEGA)
        assert is_element(m, g)
        assert g.column(int(0)) != sage_vector(FIELD, [ONE] + [ZERO] * (m + int(1)))
        assert evaluate(m, word(m, g)) == g


def test_trace_generators_commute_with_omega():
    from sage_indefinite_port.indefinite.eisenstein_hyperbolic_unitary import lattice, trace_generators

    for m in (int(1), int(2)):
        hermitian = lattice(m)
        action = hermitian.order_action(OMEGA)
        trace = hermitian.trace_lattice()
        assert trace.rank() == int(2) * m + int(4)
        for isometry in trace_generators(m):
            assert all(isometry(action(p)) == action(isometry(p)) for p in trace.module_generators())
