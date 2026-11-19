"""Exact full split/direct-product Hermitian Levi image."""

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.rings.number_field.number_field import CyclotomicField

from sage_extensions.categories.hermitian_lattices import HermitianLattice, HermitianWittFrame


def test_split_direct_product_eisenstein_line_levi_image():
    r"""Every one of the 36 pairs lifts in the split H+<1> lattice."""
    field = CyclotomicField(3)
    omega = field.gen()
    theta = omega - omega.conjugate()
    one, zero = field.one(), field.zero()

    split_gram = matrix(field, 3, 3, [zero, zero, one, zero, one, zero, one, zero, zero])
    split_lattice = HermitianLattice((one, omega), split_gram)
    split = HermitianWittFrame(split_lattice, identity_matrix(field, 3),
                               matrix(field, 1, 1, [one]), 1)
    candidates = split.levi_candidates()
    assert len(candidates) == 36
    assert len(split.levi_image()) == 36
    for A, U in candidates:
        Z, T = split.levi_lift(A, U)
        assert Z.is_zero() and T.is_zero()
        assert split.preserves_lattice(split.levi(A, U))
