r"""Witt-index-one fixed-norm classes that need long runs (#40)."""

from math import gcd
from sage_extensions.all import *
from sage.modules.free_module_element import vector as sage_vector
from sage.rings.integer_ring import ZZ as SageZZ
from sage_indefinite_port.indefinite.witt_one_hilbert import HilbertWittOneLattice
from sage_indefinite_port.indefinite.witt_one_fixed_norm_classes import (
    HilbertFixedNormClasses,
    witt_one_fixed_norm_vector_classes,
)


def _is_transporter(h, x, y, gram):
    return h is not None and h * x == y and h.transpose() * gram * h == gram


def _sqrt_ten_classes():
    return HilbertFixedNormClasses(HilbertWittOneLattice(int(10)))


def _v(*entries):
    return sage_vector(SageZZ, [SageZZ(int(c)) for c in entries])


def test_sqrt_ten_norm_two_transporter_and_divisor_obstruction():
    classes = _sqrt_ten_classes()
    gram = HilbertWittOneLattice(int(10)).gram()
    x, y = _v(1, -1, 0, 0), _v(19, -19, 0, -6)
    assert int(y * gram * y) == int(2)
    assert _is_transporter(classes.transporter(x, y), x, y, gram)
    assert classes.transporter(x, _v(0, 0, 1, 0)) is None
