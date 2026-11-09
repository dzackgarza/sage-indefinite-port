r"""diag(1, -3) over Z[i]: anisotropic since (-1, 3)_3 = -1, with SU(h) the norm-one units of Z<1,i,j,k>.

Every generator is unitary, the SU generators have determinant 1, and i * I is a scalar of U(h)
that dies in PU(h) and is not the identity.
"""

from sage_extensions.all import *
from sage.matrix.special import identity_matrix as sage_identity
from sage_indefinite_port.indefinite.anisotropic_binary_hermitian import (
    FIELD, is_anisotropic, is_element, scalar_kernel, special_unitary_generators, unitary_generators,
)


def test_the_form_is_anisotropic():
    assert is_anisotropic()


def test_generators_are_integral_and_unitary():
    assert all(is_element(g) for g in unitary_generators())
    assert all(g.det() == FIELD.one() for g in special_unitary_generators())
    assert any(g.det() != FIELD.one() for g in unitary_generators())


def test_the_special_unitary_group_is_infinite():
    r"""Some generator has infinite order, since a trace of absolute value above 2 is hyperbolic."""
    assert any(abs(g.trace().vector()[int(0)]) > int(2) for g in special_unitary_generators())


def test_scalar_kernel_is_the_fourth_roots_of_unity():
    z = scalar_kernel()
    assert is_element(z) and z != sage_identity(FIELD, int(2)) and z ** int(4) == sage_identity(FIELD, int(2))
