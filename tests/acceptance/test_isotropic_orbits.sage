"""Orbits of isotropic lines, planes and flags in indefinite lattices.

Recorded data: the existence of isotropic vectors for polyhedral_common's 8,821 forms;
Dutour Sikirić--Hulek's orbits in N = U + U(2) + E8(-2) under O^+(N) and the stable group
tilde O^+(N), with representatives and the index of each stabilizer's image in O(D(N));
Allcock's two orbits of isotropic vectors and of isotropic planes in I_{2,10}, told apart
by the parity of v^perp/v and V^perp/V. O^+ is the kernel of the real spinor norm
(Gritsenko--Hulek--Sankaran), the group of these sources.
"""

import pytest

from tests.acceptance.consumer import (
    FLAG_ORBITS,
    ISOTROPIC_ORBITS,
    ISOTROPIC_STABILIZER,
    ISOTROPIC_WITNESS,
    O_L,
    SPLIT_ORBIT,
    SUBGROUP_GENERATORS,
    element,
    lattice,
    require,
)
from tests.fixtures.oracle_fixtures import (
    load_allcock_i_2_10_orbits,
    load_isotropic_cases,
    load_unpolarized_enriques,
)


@pytest.mark.parametrize("case", load_isotropic_cases(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="vector and primitive isotropic orbits: owned by #19", strict=True)
def test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors(case) -> None:
    require(ISOTROPIC_ORBITS)

    lines = lattice(case["gram"]).O().isotropic_orbit_representatives(1)

    assert (len(lines) > 0) == case["has_isotropic"]


_N = load_unpolarized_enriques()


def _n():
    return lattice(_N["lattice"]["gram"])


@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", strict=True)
def test_o_plus_n_has_two_line_and_two_plane_orbits() -> None:
    require(ISOTROPIC_ORBITS, SPLIT_ORBIT)

    group = _n().O_plus()
    counts = _N["component_preserving_group"]["counts"]

    assert len(group.isotropic_orbit_representatives(1)) == counts["line_orbits"]
    assert len(group.isotropic_orbit_representatives(2)) == counts["plane_orbits"]


@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", strict=True)
def test_stable_o_plus_n_orbits_match_dutour_sikiric_hulek() -> None:
    require(ISOTROPIC_ORBITS, FLAG_ORBITS, SPLIT_ORBIT)

    n = _n()
    group = n.stable_orthogonal_group().intersection(n.O_plus())
    counts = _N["stable_component_preserving_group"]["counts"]

    assert len(group.isotropic_orbit_representatives(1)) == counts["line_orbits"]
    assert len(group.isotropic_orbit_representatives(2)) == counts["plane_orbits"]
    assert len(group.isotropic_orbit_representatives(2, flag=True)) == counts["flag_orbits"]


@pytest.mark.xfail(reason="O(L), isometry, vector stabilizers: owned by #18", strict=True)
def test_o_n_maps_onto_o_of_the_discriminant_form() -> None:
    require(O_L)

    n = _n()
    image = _N["discriminant_image"]

    assert int(n.discriminant_group().orthogonal_group().cardinality()) == image["O_qN_order"]
    assert int(n.O().discriminant_image().cardinality()) == image["O_qN_order"]


def _isotropic_sublattice(n, basis):
    return n.vector_configuration([element(n, vector) for vector in basis])


@pytest.mark.parametrize("rank", [1, 2])
@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", strict=True)
def test_published_representatives_lie_in_distinct_o_plus_orbits(rank) -> None:
    require(ISOTROPIC_ORBITS, ISOTROPIC_WITNESS, SPLIT_ORBIT)

    n = _n()
    group = n.O_plus()
    key = "line_representatives" if rank == 1 else "plane_representatives"
    published = [_isotropic_sublattice(n, rep["basis"]) for rep in _N["component_preserving_group"][key]]

    computed = group.isotropic_orbit_representatives(rank)

    for representative in computed:
        assert sum(group.isotropic_equivalence_witness(representative, listed) is not None for listed in published) == 1


@pytest.mark.parametrize(
    "representative",
    _N["component_preserving_group"]["line_representatives"] + _N["component_preserving_group"]["plane_representatives"],
    ids=lambda rep: rep["id"],
)
@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", strict=True)
def test_stabilizer_images_have_dutour_sikiric_hulek_indices(representative) -> None:
    require(ISOTROPIC_STABILIZER, SUBGROUP_GENERATORS)

    n = _n()
    group = n.O_plus()
    sublattice = _isotropic_sublattice(n, representative["basis"])
    target = n.discriminant_group().orthogonal_group()

    stabilizer = group.isotropic_stabilizer_generators(sublattice)
    image = target.subgroup_on(tuple(generator.discriminant_morphism() for generator in stabilizer))

    assert int(target.cardinality()) // int(image.cardinality()) == representative["discriminant_stabilizer_index"]


_ALLCOCK = load_allcock_i_2_10_orbits()
_PARITY = {"I_{1,9}": False, "II_{1,9}": True, "E8(-1)": True, "I_{0,8}": False}


@pytest.mark.xfail(reason="vector and primitive isotropic orbits: owned by #19", strict=True)
def test_i_2_10_isotropic_vectors_form_an_odd_and_an_even_orbit() -> None:
    require(ISOTROPIC_ORBITS)

    group = lattice(_ALLCOCK["gram"]).O()

    lines = group.isotropic_orbit_representatives(1)

    recorded = sorted(_PARITY[orbit["isomorphism_class"]] for orbit in _ALLCOCK["primitive_isotropic_vector_orbits"])
    assert sorted(line.isotropic_reduction().is_even() for line in lines) == recorded


@pytest.mark.xfail(reason="isotropic sublattice and flag orbits: owned by #21", strict=True)
def test_i_2_10_isotropic_planes_form_an_odd_and_an_even_orbit() -> None:
    require(ISOTROPIC_ORBITS)

    group = lattice(_ALLCOCK["gram"]).O()

    planes = group.isotropic_orbit_representatives(2)

    recorded = sorted(_PARITY[orbit["isomorphism_class"]] for orbit in _ALLCOCK["isotropic_plane_orbits"])
    assert sorted(plane.isotropic_reduction().is_even() for plane in planes) == recorded
