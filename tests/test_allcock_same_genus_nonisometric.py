"""Two reflective rank-3 lattices of one genus are not isometric (#18).

Allcock's table (references/vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/rk3.tex)
records that the genus 1^{1}_{7} 8^{1}_{1} 256^{1}_{1} holds two reflective lattices,
L_{341.1} and L_{342.1}, with different Weyl groups W_{341} and W_{342}
(rk3.tex:193-204, :2204, :2278; "shares genus with", :129658-129662 and :129851-129855;
the genus is "not a strong enough invariant", :567-573). The Weyl group of the reflections
of L is an isometry invariant, so the two lattices are not isometric. In RK3_all
(lines 8259 and 8263) and in the fixture their Weyl-group ids are 340 and 341, counted from 0.
"""

from tests.acceptance.consumer import lattice
from tests.fixtures.oracle_fixtures import load_allcock_rank3_reflective

L_341_1 = [[-239360, -11008, -11520], [-11008, -504, -528], [-11520, -528, -553]]
L_342_1 = [[13951232, 6751744, 425216], [6751744, 3267528, 205784], [425216, 205784, 12959]]


def _gram(weyl_group_id):
    (case,) = (case for case in load_allcock_rank3_reflective() if (case["weyl_group_id"], case["lattice_id"]) == (weyl_group_id, 1))
    return case


def test_fixture_holds_the_two_printed_lattices_with_equal_elementary_divisors():
    first, second = _gram(340), _gram(341)

    assert first["gram"] == L_341_1
    assert second["gram"] == L_342_1
    assert first["elementary_divisors"] == second["elementary_divisors"] == [256, 8, 1]
    assert first["num_simple_roots"] != second["num_simple_roots"]


def test_same_genus_lattices_with_different_weyl_groups_are_not_isometric():
    first, second = lattice(L_341_1), lattice(L_342_1)

    assert first.genus() == second.genus()
    assert not first.is_isometric(second)
    assert first.Isom(second).is_empty()
