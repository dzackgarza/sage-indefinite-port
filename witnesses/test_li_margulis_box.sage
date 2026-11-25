r"""Li-Margulis isotropic reductions whose box searches need long runs (#56)."""

from sage_extensions.all import *
import pytest
from sage_indefinite_port.indefinite.li_margulis import isotropic_vector_orbits_in_box


def test_isotropic_orbits_of_u_plus_u2():
    r"""Issue 56, E. L = U + U(2), A_L = u(2) with isotropic elements (1,0), (0,1) swapped by O(L).

    L contains U, so by Eichler primitive isotropic orbits are the pairs (div v, v/div v in A_L):
    div 1 one orbit, div 2 one orbit. Representatives e_1 (with f_1) and e_2 (with f_2) lie in
    the box of size 1, which is why that box suffices here.
    """
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U").twist(2)
    orbits = isotropic_vector_orbits_in_box(L, int(1))
    assert sorted(int(v.divisor()) for v in orbits) == [int(1), int(2)]


def test_isotropic_orbits_of_u_plus_u():
    r"""Unimodular L = U + U: by Eichler one orbit of primitive isotropic vectors, represented by e_1."""
    L = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    assert len(isotropic_vector_orbits_in_box(L, int(1))) == int(1)
