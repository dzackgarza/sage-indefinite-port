"""Centralizers of finite-order isometries of indefinite lattices.

Recorded data: OSCAR's QuadFormAndIsom tests. For an isometry f of L, the image in O(q_L)
of the centralizer Z_{O(L)}(f) has the recorded order (or is all of O(q_L)). The image is
f.centralizer_discriminant_image(): hermitian Miranda--Morrison theory (Brandhorst--Hofmann
§6, dzackgarza/research#410) where f has a quadratic cyclotomic minimal polynomial on an
indefinite lattice, and otherwise the image of generators of Z_{O(L)}(f) glued from its
cyclotomic components (dzackgarza/research#397). Indefinite unitary groups have no
generator algorithm in the literature, so the image is not read off generators there.
OSCAR's own centralizer computation is the recorded oracle only, never the route.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    isometry_from_rows,
    lattice,
    require,
)
from tests.fixtures.oracle_fixtures import (
    load_isometry_centralizers,
)

_DATA = load_isometry_centralizers()


@pytest.mark.parametrize("case", _DATA["centralizer_cases"], ids=lambda case: case["id"])
def test_centralizer_images_match_oscar(case) -> None:
    require(MARKED_AUTOMORPHISMS)

    lattice_ = lattice(case["gram"])
    isometry = isometry_from_rows(lattice_, case["isometry"])
    assert isometry**case["isometry_order"] == lattice_.Aut().identity()

    image = isometry.centralizer_discriminant_image()

    if "centralizer_image_order" in case:
        assert int(image.cardinality()) == case["centralizer_image_order"]
    else:
        assert image.cardinality() == lattice_.discriminant_group().orthogonal_group().cardinality()
