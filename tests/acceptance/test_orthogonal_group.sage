"""Generators of O(L) for indefinite L, against recorded generators and finite images.

A computed generating set is accepted when its group has the finite images of the
recorded group: the image in O(q_L) and the reductions modulo the primes of
``comparison_primes``. Recorded generating sets come from Mertens (arXiv:1303.3478),
Hecke's QuadBin tests and polyhedral_common's 28B suite; recorded image orders come from
OSCAR's spinor-norm tests.
"""

import pytest

from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import (
    O_L,
    assert_same_finite_images,
    discriminant_image_order,
    lattice,
    require,
    transpose,
)
from tests.fixtures.oracle_fixtures import (
    load_binary_form_automorphisms,
    load_lorentzian_stabilizer_cases,
    load_mertens_generators,
    load_oscar_lattice_oracles,
)


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_mertens_omega_and_minus_one_generate_o_l() -> None:
    require(O_L)

    case = load_mertens_generators()
    gram = case["gram"]
    minus_one = [[-1 if i == j else 0 for j in range(3)] for i in range(3)]
    recorded = [generator["matrix"] for generator in case["omega_generators"]] + [minus_one]

    computed = lattice(gram).O().framing().group_generators()

    assert_same_finite_images(gram, computed, recorded)


@pytest.mark.parametrize(
    "case",
    load_binary_form_automorphisms()["explicit_generators"],
    ids=lambda case: case["id"],
)
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_binary_form_generators_match_hecke(case) -> None:
    require(O_L)

    recorded = [transpose(generator) for generator in case["automorphism_group_generators"]]

    computed = lattice(case["gram"]).O().framing().group_generators()

    assert_same_finite_images(case["gram"], computed, recorded)


@pytest.mark.parametrize(
    "case",
    load_binary_form_automorphisms()["improper_automorphisms"],
    ids=lambda case: case["id"],
)
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_ambiguous_binary_form_has_an_improper_automorphism(case) -> None:
    require(O_L)

    computed = lattice(case["gram"]).O().framing().group_generators()

    assert any(generator.determinant() == -1 for generator in computed)


@pytest.mark.parametrize("case", load_lorentzian_stabilizer_cases(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_lorentzian_generators_match_polyhedral_common(case) -> None:
    require(O_L)

    computed = lattice(case["gram"]).O().framing().group_generators()

    assert_same_finite_images(case["gram"], computed, case["generators"])


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_orthogonal_group_of_the_hyperbolic_plane_has_order_four() -> None:
    require(O_L)

    case = next(case for case in load_oscar_lattice_oracles()["isometry_groups"] if case["id"] == "O_U_order")

    assert lattice(case["gram"]).O().order() == case["orthogonal_group_order"]


_IMAGE_CASES = [case for case in load_oscar_lattice_oracles()["discriminant_images"] if "O_qL_order" in case]


@pytest.mark.parametrize("case", _IMAGE_CASES, ids=lambda case: case["id"])
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_discriminant_images_match_oscar(case) -> None:
    """|O(q_L)|, the image of O(L), and the image of O^+(L) = ker sn_R (OSCAR's image_in_Oq_signed)."""
    require(O_L)

    lattice_ = lattice(case["gram"])

    assert int(lattice_.discriminant_group().orthogonal_group().cardinality()) == case["O_qL_order"]
    assert int(lattice_.O().discriminant_image().cardinality()) == case["image_in_Oq_order"]
    assert discriminant_image_order(lattice_, lattice_.O_plus().group_generators()) == case["image_in_Oq_signed_order"]


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q() -> None:
    require(O_L)

    case = next(case for case in load_oscar_lattice_oracles()["discriminant_images"] if case["id"] == "A7_plus_diag_1_1_-1")
    lattice_ = lattice(case["gram"])

    assert lattice_.O().discriminant_image().cardinality() == lattice_.discriminant_group().orthogonal_group().cardinality()


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q() -> None:
    require(O_L)

    case = next(case for case in load_oscar_lattice_oracles()["discriminant_images"] if case["id"] == "U_plus_-3")
    lattice_ = lattice(case["gram"])
    signed = discriminant_image_order(lattice_, lattice_.O_plus().group_generators())

    assert signed == case["image_in_Oq_signed_order"]
    assert signed == int(lattice_.discriminant_group().orthogonal_group().cardinality())
