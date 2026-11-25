"""Automorphisms and isometries of pointed indefinite lattices.

The stabilizer of a vector v is Aut((L, v)) for the pointed lattice (L, <q(v)> -> L), and
an isometry carrying v to w is an element of Isom((L, v), (L, w)). OSCAR records an
automorphism group of order 2. Brandhorst--Gonzalez-Alonso record, for every O(E10)-orbit
of primitive vectors of norm 0..10, the index of the image of its automorphism group in
O(E10 tensor F_2), the orthogonal group of even type O^+_10(2). Isometries of pointed
lattices are checked against isometries the corpus provides independently: the
reflections in Dutour Sikirić--Hulek's simple roots of U + E8(-1).
"""

from collections import Counter

import pytest
from sage.all import MatrixGroup, ZZ, matrix

from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    MARKED_CLASSES,
    automorphisms,
    element,
    is_isometry,
    lattice,
    marked_vector,
    pointed,
    pointed_category,
    reduction_image_order,
    require,
    rows_of,
    underlying_isometry,
    vector_classes,
)
from tests.fixtures.oracle_fixtures import (
    load_e10_fundamental_domain,
    load_e10_vector_orbits,
    load_oscar_lattice_oracles,
)

_E10 = load_e10_fundamental_domain()
_ORDER_O_PLUS_10_2 = 2**21 * 3**5 * 5**2 * 7 * 17 * 31


def _automorphism_rows(marked):
    category = pointed_category(marked_vector(marked).q())
    return [rows_of(underlying_isometry(generator)) for generator in automorphisms(category, marked).group_generators()]


def test_vector_stabilizer_has_oscars_order() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["isometry_groups"] if case["id"] == "stabilizer_of_100")
    lattice_ = lattice(case["gram"])

    rows = _automorphism_rows(pointed(lattice_, element(lattice_, case["vector"])))

    assert all(is_isometry(generator, case["gram"]) for generator in rows)
    assert MatrixGroup([matrix(ZZ, generator) for generator in rows]).order() == case["vector_stabilizer_order"]


@pytest.mark.parametrize("square", sorted(Counter(orbit["h_squared"] for orbit in load_e10_vector_orbits()["orbits"])))
def test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso(square) -> None:
    require(MARKED_CLASSES, MARKED_AUTOMORPHISMS)

    recorded = Counter(
        orbit["stabilizer_image_index_in_O_E10_F2"] for orbit in load_e10_vector_orbits()["orbits"] if orbit["h_squared"] == square
    )

    image_orders = [
        reduction_image_order(_E10["gram"], _automorphism_rows(marked), 2)
        for marked in vector_classes(lattice(_E10["gram"]), square)
    ]
    assert all(_ORDER_O_PLUS_10_2 % order == 0 for order in image_orders)
    computed = Counter(_ORDER_O_PLUS_10_2 // order for order in image_orders)

    assert computed == recorded
