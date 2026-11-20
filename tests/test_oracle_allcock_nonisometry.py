"""Allcock's published indefinite same-genus, nonisometric pair.

Primary source: the pinned Allcock rk3.tex and its machine-readable RK3_all
table, revision 63a9067a. The paper identifies genus
1_7 8_1 256_1 as containing L_{341.1} with Weyl group W_341 and
L_{342.1} with Weyl group W_342 (rk3.tex:196--204, 2311--2314).
It explicitly says that these lattices have the same genus but different
Weyl groups (rk3.tex:567--573). The two cross-reference entries at
rk3.tex:129658--129662 and 129851--129855 independently state that each
lattice shares its genus with the other.

The paper's notation W_n is one-based. The corresponding zero-based final
fields in RK3_all are 340 and 341, respectively (records at lines 8261 and
8263); field four is the number of simple roots, i.e. walls of the
fundamental Weyl chamber. The fixture's two Gram matrices reproduce the
first fields of those exact records, and its source metadata points to the
same RK3_all lines. Do not substitute record (342,1): it has a different
genus.
"""

from tests.acceptance.consumer import MARKED_ISOMETRY, lattice, require


def test_allcock_rank_three_same_genus_different_weyl_chambers():
    r"""The two Allcock lattices have 8 and 6 intrinsic chamber walls.

    An integral isometry sends roots to roots, hence conjugates the full
    root-reflection subgroup in O(2,1). It therefore identifies the
    associated mirror arrangements and sends chambers to chambers. The
    number of walls of a chamber is invariant under this identification.
    Allcock's classification identifies the Weyl groups as distinct and
    records 8 simple roots for W_341 and 6 for W_342, so these lattices
    cannot be isometric. This is independent of the port's isometry code.
    """
    require(MARKED_ISOMETRY)
    lattice_w341_l341_1 = lattice(
        [
            [-239360, -11008, -11520],
            [-11008, -504, -528],
            [-11520, -528, -553],
        ]
    )
    lattice_w342_l342_1 = lattice(
        [
            [13951232, 6751744, 425216],
            [6751744, 3267528, 205784],
            [425216, 205784, 12959],
        ]
    )

    assert lattice_w341_l341_1.genus() == lattice_w342_l342_1.genus()
    assert not lattice_w341_l341_1.is_isometric(lattice_w342_l342_1)
    assert lattice_w341_l341_1.Isom(lattice_w342_l342_1).is_empty()
