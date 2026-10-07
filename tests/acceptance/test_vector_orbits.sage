"""O(L)-orbits of primitive vectors of a fixed square in indefinite lattices.

Recorded orbit counts and representatives come from Dutour Sikirić--Hulek (U + E8(-1),
norms 2..72, and representatives of norm <= 30), Brandhorst--Gonzalez-Alonso (E10, norms
0..10), Allcock (I_{2,10}, norm -1) and Gritsenko--Hulek--Sankaran's Eichler criterion
(the K3 lattice). Imprimitive representatives are discarded: primitivity is an orbit
invariant, and every recorded count is a count of primitive orbits.
"""

from collections import Counter

import pytest

from tests.acceptance.consumer import element, lattice, primitive
from tests.fixtures.oracle_fixtures import (
    load_allcock_i_2_10_orbits,
    load_e10_fundamental_domain,
    load_e10_vector_orbits,
    load_enriques_polarization_orbits,
    load_k3_modular_strata,
)

_E10 = load_e10_fundamental_domain()


def _e10():
    return lattice(_E10["gram"])


@pytest.mark.parametrize("case", load_enriques_polarization_orbits(), ids=lambda case: f"2d_{case['two_d']}")
def test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek(case) -> None:
    representatives = primitive(_e10().O().vector_orbit_representatives(case["two_d"]))

    assert len(representatives) == case["primitive_vector_orbits"]


@pytest.mark.parametrize("degree", sorted({rep["degree"] for rep in _E10["orbit_representatives"]}))
def test_published_representatives_are_one_per_orbit(degree) -> None:
    group = _e10().O()
    published = [element(group.domain(), rep["vector"]) for rep in _E10["orbit_representatives"] if rep["degree"] == degree]
    computed = primitive(group.vector_orbit_representatives(degree))

    assert len(computed) == len(published)
    for representative in computed:
        assert sum(group.vectors_are_equivalent(representative, listed) for listed in published) == 1


@pytest.mark.parametrize("square", sorted(Counter(orbit["h_squared"] for orbit in load_e10_vector_orbits()["orbits"])))
def test_e10_orbit_counts_match_brandhorst_gonzalez_alonso(square) -> None:
    recorded = sum(orbit["h_squared"] == square for orbit in load_e10_vector_orbits()["orbits"])

    assert len(primitive(_e10().O().vector_orbit_representatives(square))) == recorded


def test_i_2_10_is_transitive_on_norm_minus_one_vectors() -> None:
    case = load_allcock_i_2_10_orbits()

    representatives = primitive(lattice(case["gram"]).O().vector_orbit_representatives(-1))

    assert len(representatives) == case["norm_minus_one_vector_orbits"]


def _k3_gram() -> list[list[int]]:
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    e8 = -matrix(ZZ, CartanMatrix(["E", 8]))
    gram = block_diagonal_matrix([hyperbolic, hyperbolic, hyperbolic, e8, e8])
    return [list(row) for row in gram.rows()]


@pytest.mark.parametrize("square", load_k3_modular_strata()["k3_unimodular_lattice"]["tested_represented_norms"])
def test_k3_lattice_has_one_primitive_orbit_per_norm(square) -> None:
    representatives = primitive(lattice(_k3_gram()).O().vector_orbit_representatives(square))

    assert len(representatives) == 1
