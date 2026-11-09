r"""U(1,1) of the hyperbolic plane. Issue 51's fixture is [[1, i], [0, 1]], with distinct powers."""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage_indefinite_port.indefinite.binary_hyperbolic_unitary import evaluate, field, form, generators, is_element, word


def test_issue_fixture_is_unipotent_of_infinite_order():
    d = int(1)
    E = field(d)
    g = sage_matrix(E, int(2), int(2), [E.one(), E.gen(), E.zero(), E.one()])
    assert is_element(d, g) and g ** int(5) != g ** int(2)
    assert evaluate(d, word(d, g)) == g


def test_generators_are_unitary():
    for d in (int(1), int(2), int(3)):
        assert all(is_element(d, g) for g in generators(d).values())


def test_words_recover_long_products():
    from random import Random

    rng = Random(int(11))
    for d in (int(1), int(2), int(3)):
        table = generators(d)
        names = tuple(table)
        for _trial in range(int(8)):
            product = table["W"] ** int(2)
            for _letter in range(int(16)):
                product = product * table[rng.choice(names)]
            assert evaluate(d, word(d, product)) == product


def test_non_element_has_no_word():
    from pytest import raises

    d = int(1)
    with raises(ValueError):
        word(d, field(d)(int(2)) * form(d))
