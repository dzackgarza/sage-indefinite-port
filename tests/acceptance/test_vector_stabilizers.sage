"""Stabilizers and equivalence witnesses of vectors in indefinite lattices.

OSCAR records a vector stabilizer of order 2. Brandhorst--Gonzalez-Alonso record, for every
O(E10)-orbit of primitive vectors of norm 0..10, the index of the image of its stabilizer
in O(E10 tensor F_2), the orthogonal group of even type O^+_10(2). Equivalence witnesses
are checked against isometries the corpus provides independently: the reflections in Dutour
Sikirić--Hulek's simple roots of U + E8(-1).
"""

from collections import Counter

import pytest

from tests.acceptance.consumer import coordinates, element, is_isometry, lattice, primitive, reduction_image_order, rows_of
from tests.fixtures.oracle_fixtures import (
    load_e10_fundamental_domain,
    load_e10_vector_orbits,
    load_oscar_lattice_oracles,
)

_E10 = load_e10_fundamental_domain()
_ORDER_O_PLUS_10_2 = 2**21 * 3**5 * 5**2 * 7 * 17 * 31


def test_vector_stabilizer_has_oscars_order() -> None:
    case = next(case for case in load_oscar_lattice_oracles()["isometry_groups"] if case["id"] == "stabilizer_of_100")
    lattice_ = lattice(case["gram"])

    generators = lattice_.O().vector_stabilizer_generators(element(lattice_, case["vector"]))

    assert all(is_isometry(rows_of(generator), case["gram"]) for generator in generators)
    assert MatrixGroup([matrix(ZZ, rows_of(generator)) for generator in generators]).order() == case["vector_stabilizer_order"]


def test_o_plus_10_2_has_the_order_dutour_sikiric_hulek_state() -> None:
    assert GO(10, GF(2), e=1).order() == _ORDER_O_PLUS_10_2


@pytest.mark.parametrize("square", sorted(Counter(orbit["h_squared"] for orbit in load_e10_vector_orbits()["orbits"])))
def test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso(square) -> None:
    group = lattice(_E10["gram"]).O()
    recorded = Counter(
        orbit["stabilizer_image_index_in_O_E10_F2"] for orbit in load_e10_vector_orbits()["orbits"] if orbit["h_squared"] == square
    )

    computed = Counter()
    for representative in primitive(group.vector_orbit_representatives(square)):
        generators = [rows_of(generator) for generator in group.vector_stabilizer_generators(representative)]
        computed[_ORDER_O_PLUS_10_2 // reduction_image_order(generators, 2)] += 1

    assert computed == recorded


def _reflection_rows(gram, root) -> list[list[int]]:
    """Row matrix of the reflection in a norm -2 root: v -> v + b(v, r) r."""
    form = matrix(ZZ, gram)
    r = vector(ZZ, root)
    return [list(basis_vector + (basis_vector * form * r) * r) for basis_vector in identity_matrix(ZZ, len(root)).rows()]


@pytest.mark.parametrize("root", _E10["simple_roots"], ids=lambda root: f"r_{root['label']}")
def test_equivalence_witness_carries_a_representative_to_its_reflection(root) -> None:
    gram = _E10["gram"]
    group = lattice(gram).O()
    reflection = matrix(ZZ, _reflection_rows(gram, root["vector"]))
    assert is_isometry(reflection.rows(), gram)
    for representative in _E10["orbit_representatives"]:
        source = vector(ZZ, representative["vector"])
        target = source * reflection

        witness = group.vector_equivalence_witness(element(group.domain(), source), element(group.domain(), target))

        assert is_isometry(rows_of(witness), gram)
        assert vector(ZZ, coordinates(witness(element(group.domain(), source)))) == target
