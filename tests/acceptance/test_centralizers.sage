"""Centralizers of finite-order isometries of indefinite lattices.

Recorded data: OSCAR's QuadFormAndIsom tests. For an isometry f of L, the image in O(q_L)
of the centralizer Z_{O(L)}(f) has the recorded order (or is all of O(q_L)). The image is
computed from the generators of Z_{O(L)}(f) that this port returns, never from OSCAR's
own centralizer computation.
"""

import pytest

from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import (
    CENTRALIZER,
    isometry_from_rows,
    lattice,
    require,
)
from tests.fixtures.oracle_fixtures import (
    load_isometry_centralizers,
)

_DATA = load_isometry_centralizers()


@pytest.mark.parametrize("case", _DATA["centralizer_cases"], ids=lambda case: case["id"])
@pytest.mark.xfail(reason="centralizers of indefinite lattices: owned by #23", raises=UnfinishedCapability, strict=True)
def test_centralizer_images_match_oscar(case) -> None:
    require(CENTRALIZER)

    lattice_ = lattice(case["gram"])
    isometry = isometry_from_rows(lattice_, case["isometry"])
    assert isometry**case["isometry_order"] == lattice_.Aut().identity()

    target = lattice_.discriminant_group().orthogonal_group()
    centralizer = lattice_.O().centralizer(isometry)
    image = target.subgroup_on(tuple(generator.discriminant_morphism() for generator in centralizer.group_generators()))

    if "centralizer_image_order" in case:
        assert int(image.cardinality()) == case["centralizer_image_order"]
    else:
        assert image.cardinality() == lattice_.discriminant_group().orthogonal_group().cardinality()


@pytest.mark.xfail(reason="research-owned: Genus.equivariant_classes is missing in the preamble; the case moves to research's suite under #33", raises=AttributeError, strict=True)
def test_involution_classes_in_the_signature_1_9_genus() -> None:
    (case,) = _DATA["involution_classes"]
    x = polygen(ZZ, "x")

    classes = lattice(case["genus_representative_gram"]).genus().equivariant_classes((x - 1) ** 4 * (x + 1) ** 6)

    assert len(classes) == case["classes_in_genus"]


@pytest.mark.xfail(reason="research-owned: Genus.equivariant_classes is missing in the preamble; the case moves to research's suite under #33", raises=AttributeError, strict=True)
def test_u_plus_e8_has_one_class_of_order_30() -> None:
    case = next(case for case in _DATA["lattice_class_counts"] if case["id"] == "U_plus_E8_order_30")
    x = polygen(ZZ, "x")

    classes = lattice(case["gram"]).genus().equivariant_classes((x - 1) ** 2 * cyclotomic_polynomial(30, x))

    assert len(classes) == case["class_count"]
    (representative,) = classes
    assert representative.isometry().invariant_lattice().discriminant() == -1


@pytest.mark.xfail(reason="research-owned: Genus.equivariant_classes is missing in the preamble; the case moves to research's suite under #33", raises=AttributeError, strict=True)
def test_4u_has_three_hermitian_classes_of_order_5() -> None:
    """representatives_of_hermitian_type(L, 5): f of order 5 without fixed vectors, so chi_f = Phi_5^2."""
    case = next(case for case in _DATA["lattice_class_counts"] if case["id"] == "4U_hermitian_order_5")
    x = polygen(ZZ, "x")

    classes = lattice(case["gram"]).genus().equivariant_classes(cyclotomic_polynomial(5, x) ** 2)

    assert len(classes) == case["class_count"]
