"""Generators of O(L) for indefinite L, against recorded generators and finite images.

A computed generating set is accepted when its group has the finite images of the
recorded group: the image in O(q_L) and the reductions modulo the primes of
``comparison_primes``. Recorded generating sets come from Mertens (arXiv:1303.3478),
Hecke's QuadBin tests and polyhedral_common's 28B suite; recorded image orders come from
OSCAR's spinor-norm tests.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    assert_same_finite_images,
    discriminant_image_order,
    is_isometry,
    lattice,
    require,
    rows_of,
    transpose,
)
from tests.fixtures.oracle_fixtures import (
    load_binary_form_automorphisms,
    load_lorentzian_stabilizer_cases,
    load_mertens_generators,
    load_oscar_lattice_oracles,
)


def test_mertens_omega_and_minus_one_generate_o_l() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = load_mertens_generators()
    gram = case["gram"]
    minus_one = [[-1 if i == j else 0 for j in range(3)] for i in range(3)]
    recorded = [generator["matrix"] for generator in case["omega_generators"]] + [minus_one]

    computed = lattice(gram).O().group_generators()

    assert_same_finite_images(gram, computed, recorded)


@pytest.mark.parametrize(
    "case",
    load_binary_form_automorphisms()["explicit_generators"],
    ids=lambda case: case["id"],
)
def test_binary_form_generators_match_hecke(case) -> None:
    require(MARKED_AUTOMORPHISMS)

    recorded = [transpose(generator) for generator in case["automorphism_group_generators"]]

    computed = lattice(case["gram"]).O().group_generators()

    assert_same_finite_images(case["gram"], computed, recorded)


@pytest.mark.parametrize(
    "case",
    load_binary_form_automorphisms()["improper_automorphisms"],
    ids=lambda case: case["id"],
)
def test_ambiguous_binary_form_has_an_improper_automorphism(case) -> None:
    require(MARKED_AUTOMORPHISMS)

    computed = lattice(case["gram"]).O().group_generators()

    assert any(int(generator.determinant()) == -1 for generator in computed)


@pytest.mark.parametrize("case", load_lorentzian_stabilizer_cases(), ids=lambda case: case["id"])
def test_lorentzian_generators_match_polyhedral_common(case) -> None:
    require(MARKED_AUTOMORPHISMS)

    computed = lattice(case["gram"]).O().group_generators()

    assert_same_finite_images(case["gram"], computed, case["generators"])


def test_degenerate_orthogonal_group_contains_radical_shear_and_gl_action() -> None:
    """For <2> plus a rank-one radical, O(L) includes both signs and every shear."""
    require(MARKED_AUTOMORPHISMS)

    gram = [[2, 0], [0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(tuple(tuple(row) for row in rows_of(generator)) for generator in computed)
    expected = {
        ((-1, 0), (0, 1)),
        ((1, 0), (0, -1)),
        ((1, 1), (0, 1)),
    }

    assert all(is_isometry(rows_of(generator), gram) for generator in computed)
    assert expected <= set(rows)


def test_degenerate_orthogonal_group_generates_rank_two_radical_action() -> None:
    """The radical GL_2(Z) action and the nonsingular-to-radical shear survive assembly."""
    require(MARKED_AUTOMORPHISMS)

    gram = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(tuple(tuple(row) for row in rows_of(generator)) for generator in computed)

    assert all(is_isometry(rows_of(generator), gram) for generator in computed)
    assert any(generator_rows[0][1] != 0 for generator_rows in rows)
    assert any(generator_rows[1][2] != 0 or generator_rows[2][1] != 0 for generator_rows in rows)
    assert any(
        tuple(tuple(generator_rows[i][j] for j in (1, 2)) for i in (1, 2)) != ((1, 0), (0, 1))
        for generator_rows in rows
    )


def test_degenerate_orthogonal_group_embeds_non_diagonal_complement() -> None:
    """A non-diagonal nonsingular summand embeds in the radical extension's coordinate frame."""
    require(MARKED_AUTOMORPHISMS)

    gram = [[2, 1, 0], [1, 2, 0], [0, 0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(rows_of(generator) for generator in computed)

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert any(generator_rows[0][2] != 0 for generator_rows in rows)
    assert any(generator_rows[1][2] != 0 for generator_rows in rows)
    assert any(generator_rows[2][2] == -1 for generator_rows in rows)


def test_zero_form_orthogonal_group_is_the_full_general_linear_group() -> None:
    """When the whole lattice is radical, assembly is exactly GL_2(Z)."""
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    require(MARKED_AUTOMORPHISMS)

    gram = [[0, 0], [0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(rows_of(generator) for generator in computed)

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert all(abs(int(matrix(SageZZ, generator_rows).determinant())) == 1 for generator_rows in rows)
    assert any(generator_rows[0][1] != 0 or generator_rows[1][0] != 0 for generator_rows in rows)
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
