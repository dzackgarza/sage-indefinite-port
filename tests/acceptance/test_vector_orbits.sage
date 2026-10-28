"""O(L)-orbits of primitive vectors of a fixed square in indefinite lattices.

An orbit of primitive vectors of square n is an isomorphism class of pointed lattices
(L, <n> -> L) with L fixed and the marking primitive. Recorded orbit counts and
representatives come from Dutour Sikirić--Hulek (U + E8(-1), norms 2..72, and
representatives of norm <= 30), Brandhorst--Gonzalez-Alonso (E10, norms 0..10), Allcock
(I_{2,10}, norm -1) and Gritsenko--Hulek--Sankaran's Eichler criterion (the K3 lattice).
"""

from collections import Counter

import pytest

from tests.acceptance.consumer import (
    MARKED_CLASSES,
    MARKED_ISOMETRY,
    are_isomorphic,
    element,
    lattice,
    pointed,
    pointed_category,
    require,
    vector_classes,
)
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
    require(MARKED_CLASSES)

    assert len(vector_classes(_e10(), case["two_d"])) == case["primitive_vector_orbits"]


@pytest.mark.parametrize("degree", sorted({rep["degree"] for rep in _E10["orbit_representatives"]}))
def test_published_representatives_are_one_per_orbit(degree) -> None:
    require(MARKED_CLASSES, MARKED_ISOMETRY)

    lattice_ = _e10()
    category = pointed_category(degree)
    published = [pointed(lattice_, element(lattice_, rep["vector"])) for rep in _E10["orbit_representatives"] if rep["degree"] == degree]
    computed = vector_classes(lattice_, degree)

    assert len(computed) == len(published)
    for representative in computed:
        assert sum(are_isomorphic(category, representative, listed) for listed in published) == 1


@pytest.mark.parametrize("square", sorted(Counter(orbit["h_squared"] for orbit in load_e10_vector_orbits()["orbits"])))
def test_e10_orbit_counts_match_brandhorst_gonzalez_alonso(square) -> None:
    require(MARKED_CLASSES)

    recorded = sum(orbit["h_squared"] == square for orbit in load_e10_vector_orbits()["orbits"])

    assert len(vector_classes(_e10(), square)) == recorded


def test_i_2_10_is_transitive_on_norm_minus_one_vectors() -> None:
    require(MARKED_CLASSES)

    case = load_allcock_i_2_10_orbits()

    assert len(vector_classes(lattice(case["gram"]), -1)) == case["norm_minus_one_vector_orbits"]


def test_anisotropic_binary_vector_classes_use_primitive_embeddings() -> None:
    """An anisotropic (1,1) plane has no perfect-cell cusp to enumerate."""
    require(MARKED_CLASSES, MARKED_ISOMETRY)

    binary = lattice([[2, 0], [0, -3]])
    assert not binary.is_isotropic()
    assert vector_classes(binary, 0) == ()
    # 2*x^2 - 3*y^2 = 1 has no integer solution modulo 3.
    assert vector_classes(binary, 1) == ()
    for square in (2, -3):
        classes = vector_classes(binary, square)
        assert classes
        category = pointed_category(square)
        canonical = pointed(binary, element(binary, (1, 0) if square == 2 else (0, 1)))
        assert sum(are_isomorphic(category, marked, canonical) for marked in classes) == 1
        witnesses = []
        for marked in classes:
            generators = tuple(marked.arrow()(g) for g in marked.arrow().domain().module_generators())
            assert len(generators) == 1
            witness, = generators
            assert witness.q() == square
            assert witness.content() == 1
            assert all(not are_isomorphic(category, marked, previous) for previous in witnesses)
            witnesses.append(marked)


def _k3_gram() -> list[list[int]]:
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    e8 = -matrix(ZZ, CartanMatrix(["E", 8]))
    gram = block_diagonal_matrix([hyperbolic, hyperbolic, hyperbolic, e8, e8])
    return [list(row) for row in gram.rows()]


@pytest.mark.parametrize("square", load_k3_modular_strata()["k3_unimodular_lattice"]["tested_represented_norms"])
def test_k3_lattice_has_one_primitive_orbit_per_norm(square) -> None:
    require(MARKED_CLASSES)

    assert len(vector_classes(lattice(_k3_gram()), square)) == 1
