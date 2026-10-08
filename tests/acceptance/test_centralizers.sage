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
