"""Independently sourced orthogonal-group properties for indefinite lattices.

Recorded outputs come from Hecke's QuadBin tests and OSCAR's spinor-norm tests.
Finite reduction *orders* alone do not establish equality of finite image
subgroups, let alone completeness of generators of an infinite integral
orthogonal group. The former finite-order matching tests have been removed
as mathematical acceptance assertions.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    discriminant_image_order,
    lattice,
    require,
)
from tests.fixtures.oracle_fixtures import (
    load_binary_form_automorphisms,
    load_oscar_lattice_oracles,
)


@pytest.mark.parametrize(
    "case",
    load_binary_form_automorphisms()["improper_automorphisms"],
    ids=lambda case: case["id"],
)
def test_ambiguous_binary_form_has_an_improper_automorphism(case) -> None:
    require(MARKED_AUTOMORPHISMS)

    computed = lattice(case["gram"]).O().group_generators()

    assert any(int(generator.determinant()) == -1 for generator in computed)


def test_orthogonal_group_of_the_hyperbolic_plane_has_order_four() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["isometry_groups"] if case["id"] == "O_U_order")

    assert lattice(case["gram"]).O().order() == case["orthogonal_group_order"]


_IMAGE_CASES = [case for case in load_oscar_lattice_oracles()["discriminant_images"] if "O_qL_order" in case]


@pytest.mark.parametrize("case", _IMAGE_CASES, ids=lambda case: case["id"])
def test_discriminant_images_match_oscar(case) -> None:
    """|O(q_L)|, the image of O(L), and the image of O^+(L) = ker sn_R (OSCAR's image_in_Oq_signed)."""
    require(MARKED_AUTOMORPHISMS)

    lattice_ = lattice(case["gram"])

    assert int(lattice_.discriminant_group().orthogonal_group().cardinality()) == case["O_qL_order"]
    assert int(lattice_.O().discriminant_image().cardinality()) == case["image_in_Oq_order"]
    assert discriminant_image_order(lattice_, lattice_.O_plus().group_generators()) == case["image_in_Oq_signed_order"]


def test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["discriminant_images"] if case["id"] == "A7_plus_diag_1_1_-1")
    lattice_ = lattice(case["gram"])

    assert lattice_.O().discriminant_image().cardinality() == lattice_.discriminant_group().orthogonal_group().cardinality()


def test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["discriminant_images"] if case["id"] == "U_plus_-3")
    lattice_ = lattice(case["gram"])
    signed = discriminant_image_order(lattice_, lattice_.O_plus().group_generators())

    assert signed == case["image_in_Oq_signed_order"]
    assert signed == int(lattice_.discriminant_group().orthogonal_group().cardinality())
