"""
Shared utilities for sage_indefinite_port.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, matrix
from dzack_research.preamble.categories.lattices import Lattices
from dzack_research.preamble.tensors import tensor


def signature_pair(gram):
    """
    Compute (p, q) signature via preamble L.signature_pair().
    """
    C = Lattices(ZZ)
    n = len(gram)
    rows = tuple(tuple(int(x) for x in row) for row in gram)
    g = tensor(ZZ, (), (n, n), rows)
    L = C(g)
    return L.signature_pair()


def det_abs(gram):
    """Compute |det(gram)| via preamble L.discriminant()."""
    C = Lattices(ZZ)
    n = len(gram)
    rows = tuple(tuple(int(x) for x in row) for row in gram)
    g = tensor(ZZ, (), (n, n), rows)
    L = C(g)
    return abs(int(L.discriminant()))
