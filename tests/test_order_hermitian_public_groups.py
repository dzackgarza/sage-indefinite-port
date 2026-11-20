"""Complete unitary group on the original nonmaximal/nonfree integral module."""

import pytest

from sage.matrix.constructor import block_matrix, diagonal_matrix, identity_matrix, matrix
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF

from sage_extensions.algebras.order_hermitian_modules import OrderHermitianModule


def _module(kind):
    one = identity_matrix(QQ, 2)
    i = matrix(QQ, [[0, -2], [QQ(1) / 2, 0]])
    if kind == "rank_one_nonmaximal":
        return OrderHermitianModule(diagonal_matrix(ZZ, [2, 8]), [2*i], [one, i])
    if kind == "rank_two_proper":
        standard = matrix(QQ, [[0, -1], [1, 0]])
        zero = matrix(QQ, 2, 2)
        maximal = block_matrix(QQ, [[standard, zero], [zero, i]], subdivide=False)
        return OrderHermitianModule(
            diagonal_matrix(ZZ, [2, 2, 2, 8]), [2*maximal],
            [identity_matrix(QQ, 4), maximal],
        )
    if kind == "nonfree_ideal":
        root = matrix(ZZ, [[-1, -3], [2, 1]])  # square -5, nonprincipal ideal.
        gram = matrix(ZZ, [[8, 4], [4, 12]])
        assert root.transpose() * gram * root == 5*gram
        return OrderHermitianModule(gram, [root], [one, root])
    raise ValueError(kind)


@pytest.mark.parametrize("kind,expected_order,conductor", (
    ("rank_one_nonmaximal", 2, 2),
    ("rank_two_proper", 8, 2),
    ("nonfree_ideal", 2, 1),
))
def test_arithmetic_unitary_generation_respects_original_order_and_ideal(kind, expected_order, conductor):
    r"""All integral unitaries are recovered, with exact finite group order."""
    module = _module(kind)
    assert int(module.conductor_exponent) == conductor
    generators = tuple(module.unitary_group_generators())
    assert generators
    assert all(module.is_unitary(g) for g in generators)
    assert int(MatrixGroup([g.change_ring(GF(3)) for g in generators]).order()) == expected_order
    if kind == "rank_two_proper":
        assert expected_order == 8 < 32  # Saturating to Z[i]^2 changes the group.


def test_nonmaximal_indefinite_hyperbolic_order_has_full_infinite_integral_generators():
    r"""U(H;Z[2i]) contains an explicit infinite-order integral unipotent.

    The algorithm computes the complete integer commutant in the original
    basis, not the larger U(H;Z[i]) or only a finite congruence image.
    The finite mod-2 image of the computed family has order eight, while
    its integral group remains infinite.
    """
    i = matrix(QQ, [[0, -2], [QQ(1)/2, 0]])
    zero = matrix(QQ, 2, 2)
    i_action = block_matrix(QQ, [[i, zero], [zero, i]], subdivide=False)
    gram = matrix(ZZ, [[0, 0, 2, 0], [0, 0, 0, 8], [2, 0, 0, 0], [0, 8, 0, 0]])
    module = OrderHermitianModule(gram, [2*i_action], [identity_matrix(QQ, 4), i_action])
    assert int(module.conductor_exponent) == 2
    generators = tuple(module.unitary_group_generators())
    assert generators and all(module.is_unitary(g) for g in generators)
    one = identity_matrix(ZZ, 4)
    assert any(g != one and (g-one)**2 == 0 for g in generators)
    assert int(MatrixGroup([g.change_ring(GF(2)) for g in generators]).order()) == 8

    # An independently exhibited nonzero unipotent of the exact order R:
    # the upper-right block is multiplication by 2i on (1,2i).
    shear = matrix(ZZ, one)
    shear.set_block(0, 2, (2*i).change_ring(ZZ))
    assert module.is_unitary(shear)
    assert (shear-one)**2 == 0 and shear != one
