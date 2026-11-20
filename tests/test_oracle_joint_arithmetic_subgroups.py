"""Elementary finite orthogonal groups and exact arithmetic preimages.

No expected group order below is imported from an earlier test or from the
port. For U=[[0,1],[1,0]], an integral isometry A=[[a,b],[c,d]] satisfies
ac=bd=0 and ad+bc=1. If a is nonzero, c=0, whence ad=1;
then b=0 and d=a^{-1}=a=+/-1, giving +/-I. If a=0, then bc=1,
so b,c are units; bd=0 forces d=0, giving +/-S. These are all cases.

For I_2, each column of an integral isometry has squared length one, hence
is one of ±e_1,±e_2; orthogonality forces the other column to be the other
coordinate vector. Thus O(I_2) is the eight signed permutation matrices.
Modulo 3 this representation is faithful because entries in {-1,0,1}
are distinct modulo 3. The swap subgroup has two elements, so its inverse
image has index 8/2=4. It is nonnormal: conjugating the swap by a sign
change gives the negative swap, which is not in the subgroup.
"""

from sage.all import GF, MatrixGroup, ZZ, matrix

from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import isometry_from_rows, lattice


def test_hyperbolic_plane_joint_character_and_level_image():
    L = lattice(((0, 1), (1, 0)))
    recursion = EichlerRecursion()
    e = isometry_from_rows(L, ((1, 0), (0, 1)))
    minus_e = isometry_from_rows(L, ((-1, 0), (0, -1)))
    swap = isometry_from_rows(L, ((0, 1), (1, 0)))
    minus_swap = isometry_from_rows(L, ((0, -1), (-1, 0)))
    four = (e, minus_e, swap, minus_swap)

    # Direct determinant and Cartan--Dieudonne spinor signs. The reflection
    # in e-f has norm -2 and normalized sn_{-b/2}=+1; e+f has norm +2
    # and normalized sn_{-b/2}=-1. Positive-cone sign is read on e+f.
    spinor = L.O_plus().preimage_morphism()
    expected = ((1, 1, 1), (1, -1, -1), (-1, 1, 1), (-1, -1, -1))
    actual = []
    for g in four:
        actual.append((int(g.determinant()), 1 if spinor(g).is_one() else -1,
                       1 if g.preserves_positive_cone() else -1))
    assert tuple(actual) == expected

    # Actual finite image of all requested simultaneous factors.
    rep = recursion.joint_discriminant_character_representation(
        L, determinant=True, real_spinor=True, positive_cone=True,
        level_modulus=2,
    )
    assert int(rep.codomain().cardinality()) == 4
    images = tuple(rep(g) for g in four)
    assert len(set(images)) == 4
    assert images[0] == rep.codomain().one()
    joint_image = recursion.joint_discriminant_character_image(
        L, determinant=True, real_spinor=True, positive_cone=True,
        level_modulus=2,
    )
    kernel_index, kernel_generators = joint_image.preimage_generators(
        (joint_image.identity,)
    )
    assert int(kernel_index) == 4
    kernel = recursion.joint_discriminant_character_preimage_subgroup(
        L, (joint_image.identity,), determinant=True, real_spinor=True,
        positive_cone=True, level_modulus=2,
    )
    assert all((g in kernel) == (g == e) for g in four)


def test_positive_definite_binary_level_three_nonnormal_preimage():
    L = lattice(((1, 0), (0, 1)))
    recursion = EichlerRecursion()
    # Elementary column-norm argument above proves the ambient group is D8.
    # Its mod-3 action is faithful, so the actual finite image has order 8.
    reduction = recursion.congruence_image(L, 3)
    reduced_generators = tuple(reduction.reduce(g) for g in reduction.generators)
    image = MatrixGroup([matrix(GF(3), 2, 2, g.list()) for g in reduced_generators])
    assert int(image.order()) == 8

    swap = matrix(GF(3), ((0, 1), (1, 0)))
    allowed_group = MatrixGroup([swap])
    assert int(allowed_group.order()) == 2
    swap_lift = matrix(ZZ, ((0, 1), (1, 0)))
    allowed = (reduction.image_one, reduction.reduce(swap_lift))
    index, schreier = reduction.preimage_generators(allowed)
    assert int(index) == 4
    assert schreier

    # The negative swap is outside <swap>; hence the inverse image is
    # genuinely nonnormal. Its full preimage generators are not merely
    # the lifts of the two subgroup elements.
    negative_swap = -swap
    assert negative_swap not in allowed_group
    sign_change = matrix(GF(3), ((-1, 0), (0, 1)))
    assert sign_change * swap * sign_change.inverse() == negative_swap
    assert int(MatrixGroup([matrix(ZZ, g) for g in schreier]).order()) == 2
