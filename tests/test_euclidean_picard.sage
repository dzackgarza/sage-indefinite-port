r"""Zhao's five matrices for d = 2, 7, 11 in the antidiagonal form.

Cusp relations are those of the source presentations of each d. The trace form on O_d^3 has signature (4,2).
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage_indefinite_port.indefinite.euclidean_picard import (
    field, form, full_unitary_generators, generators, is_element, lattice,
)


def _identity(d):
    return generators(d)["R1"] ** int(2)


def test_source_matrices_are_integral_unitary():
    for d in (int(2), int(7), int(11)):
        assert all(is_element(d, g) for g in generators(d).values())
        assert len(full_unitary_generators(d)) == int(6)


def test_source_cusp_presentations():
    two = int(2)
    g = generators(two)
    one = _identity(two)
    assert all(g[n] ** two == one and g["T"] * g[n] == g[n] * g["T"] for n in ("R1", "R2", "R3"))
    assert (g["T"] ** two * g["R1"] * g["R3"] * g["R2"]) ** two == one

    seven = int(7)
    g = generators(seven)
    one = _identity(seven)
    assert g["R1"] ** two == one and g["R3"] ** two == one
    assert g["T"] * g["R1"] == g["R1"] * g["T"] and g["T"] * g["R3"] == g["R3"] * g["T"]
    assert g["R2"] ** two == g["T"]
    assert (g["R1"] * g["R3"] * g["R2"]) ** two == one

    eleven = int(11)
    g = generators(eleven)
    one = _identity(eleven)
    assert g["R1"] ** two == one and g["T"] * g["R1"] == g["R1"] * g["T"]
    assert g["R2"] ** two == g["T"] and g["R3"] ** two == g["T"]
    assert (g["R1"] * g["R3"] * g["R2"]) ** two == g["T"]


def test_maximal_order_is_used_at_seven():
    g = generators(int(7))["R2"]
    a = field(int(7)).gen()
    assert a not in tuple(g.list())
    assert g[0, 2].denominator() == int(2)
    assert is_element(int(7), g)


def test_trace_form_signature():
    for d in (int(2), int(7), int(11)):
        assert tuple(int(x) for x in lattice(d).trace_lattice().signature_pair()) == (int(4), int(2))
