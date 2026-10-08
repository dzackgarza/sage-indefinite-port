"""Isotropic lines, planes and flags in indefinite lattices, as marked lattices.

An orbit of primitive totally isotropic sublattices of rank k is an isomorphism class of
marked lattices (L, I -> L) in Ar(Lattices) with L fixed, and a flag orbit is a class in
Fun([2], Lattices). Orbits under a subgroup of O(L) split those classes. Recorded data: the
existence of isotropic vectors for polyhedral_common's 8,821 forms; Dutour Sikirić--Hulek's
orbits in N = U + U(2) + E8(-2) under O^+(N) and the stable group tilde O^+(N), with
representatives and the index of each stabilizer's image in O(D(N)); Allcock's two orbits of
isotropic vectors and of isotropic planes in I_{2,10}, told apart by the parity of v^perp/v
and V^perp/V. O^+ is the kernel of the real spinor norm (Gritsenko--Hulek--Sankaran), the
group of these sources.
"""

import pytest

from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    MARKED_CLASSES,
    MARKED_ISOMETRY,
    element,
    flag_classes,
    isomorphism_in,
    isotropic_classes,
    lattice,
    marked_sublattice,
    require,
    subgroup_classes,
    sublattice_category,
    sublattice_of,
)
from tests.fixtures.oracle_fixtures import (
    load_allcock_i_2_10_orbits,
    load_isotropic_cases,
    load_unpolarized_enriques,
)


@pytest.mark.parametrize("case", load_isotropic_cases(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors(case) -> None:
    require(MARKED_CLASSES)

    assert (len(isotropic_classes(lattice(case["gram"]), 1)) > 0) == case["has_isotropic"]


_N = load_unpolarized_enriques()


def _n():
    return lattice(_N["lattice"]["gram"])


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_o_plus_n_has_two_line_and_two_plane_orbits() -> None:
    require(MARKED_CLASSES)

    n = _n()
    counts = _N["component_preserving_group"]["counts"]

    assert len(subgroup_classes(n.O_plus(), isotropic_classes(n, 1))) == counts["line_orbits"]
    assert len(subgroup_classes(n.O_plus(), isotropic_classes(n, 2))) == counts["plane_orbits"]


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_stable_o_plus_n_orbits_match_dutour_sikiric_hulek() -> None:
    require(MARKED_CLASSES)

    n = _n()
    group = n.stable_orthogonal_group().intersection(n.O_plus())
    counts = _N["stable_component_preserving_group"]["counts"]

    assert len(subgroup_classes(group, isotropic_classes(n, 1))) == counts["line_orbits"]
    assert len(subgroup_classes(group, isotropic_classes(n, 2))) == counts["plane_orbits"]
    assert len(subgroup_classes(group, flag_classes(n))) == counts["flag_orbits"]


@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_o_n_maps_onto_o_of_the_discriminant_form() -> None:
    require(MARKED_AUTOMORPHISMS)

    n = _n()
    image = _N["discriminant_image"]

    assert int(n.discriminant_group().orthogonal_group().cardinality()) == image["O_qN_order"]
    assert int(n.O().discriminant_image().cardinality()) == image["O_qN_order"]


def _published(n, key):
    return [marked_sublattice(n, [element(n, vector) for vector in rep["basis"]]) for rep in _N["component_preserving_group"][key]]


@pytest.mark.parametrize("rank", [1, 2])
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_published_representatives_lie_in_distinct_o_plus_orbits(rank) -> None:
    require(MARKED_CLASSES, MARKED_ISOMETRY)

    n = _n()
    group = n.O_plus()
    published = _published(n, "line_representatives" if rank == 1 else "plane_representatives")

    computed = subgroup_classes(group, isotropic_classes(n, rank))

    for representative in computed:
        assert sum(isomorphism_in(group, sublattice_category(), representative, listed) is not None for listed in published) == 1


@pytest.mark.parametrize(
    "representative",
    _N["component_preserving_group"]["line_representatives"] + _N["component_preserving_group"]["plane_representatives"],
    ids=lambda rep: rep["id"],
)
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_stabilizer_images_have_dutour_sikiric_hulek_indices(representative) -> None:
    require(MARKED_AUTOMORPHISMS)

    n = _n()
    sublattice = sublattice_of(marked_sublattice(n, [element(n, vector) for vector in representative["basis"]]))
    target = n.discriminant_group().orthogonal_group()

    stabilizer = n.O_plus().stabilizer(sublattice, action="setwise")
    image = target.subgroup_on(tuple(generator.discriminant_morphism() for generator in stabilizer.group_generators()))

    assert int(target.cardinality()) // int(image.cardinality()) == representative["discriminant_stabilizer_index"]


_ALLCOCK = load_allcock_i_2_10_orbits()
_PARITY = {"I_{1,9}": False, "II_{1,9}": True, "E8(-1)": True, "I_{0,8}": False}


@pytest.mark.parametrize("rank, key", [(1, "primitive_isotropic_vector_orbits"), (2, "isotropic_plane_orbits")])
@pytest.mark.xfail(reason="recursive indefinite algorithm: owned by #18", raises=UnfinishedCapability, strict=True)
def test_i_2_10_isotropic_classes_are_one_odd_and_one_even(rank, key) -> None:
    require(MARKED_CLASSES)

    classes = isotropic_classes(lattice(_ALLCOCK["gram"]), rank)

    recorded = sorted(_PARITY[orbit["isomorphism_class"]] for orbit in _ALLCOCK[key])
    assert sorted(sublattice_of(marked).isotropic_reduction().is_even() for marked in classes) == recorded
