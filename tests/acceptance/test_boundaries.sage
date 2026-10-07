"""Baily--Borel boundaries: orbits of isotropic lines (points), planes (curves) and flags (edges).

For an arithmetic group Gamma acting on a lattice of signature (2, n), the Type III points,
Type II curves and their incidences are the Gamma-orbits of isotropic lines, isotropic
planes and line-in-plane flags. Recorded boundaries: Dawes's Tits buildings
(arXiv:2205.10601, 2108.06236), Laza--O'Grady's D-tower (arXiv:1801.04845 with the group of
arXiv:1607.01324), Scattone's K3 boundaries of degree 2 and 4 (via arXiv:1205.3144 and
2502.04301), and Dutour Sikirić--Hulek's 87 Enriques modular groups (arXiv:2302.01679).
"""

import pytest

from sage_indefinite_port.readiness import UnfinishedCapability
from tests.acceptance.consumer import (
    EQUIVARIANT_LATTICE,
    FLAG_ORBITS,
    ISOTROPIC_ORBITS,
    SPLIT_ORBIT,
    SUBGROUP_GENERATORS,
    basis,
    element,
    lattice,
    require,
)
from tests.fixtures.oracle_fixtures import (
    load_dawes_buildings,
    load_dtower_boundaries,
    load_e10_fundamental_domain,
    load_enriques_87_polarizations,
    load_enriques_classical_indices,
    load_k3_modular_strata,
)


def _counts(group) -> tuple[int, int, int]:
    return (
        len(group.isotropic_orbit_representatives(1)),
        len(group.isotropic_orbit_representatives(2)),
        len(group.isotropic_orbit_representatives(2, flag=True)),
    )


def _stable_o_plus(lattice_):
    return lattice_.stable_orthogonal_group().intersection(lattice_.O_plus())


def _o_plus_fixing(lattice_, discriminant_class):
    """{g in O^+(L) : g fixes the class in A_L}."""
    target = lattice_.discriminant_group().orthogonal_group()
    preimage = lattice_.Aut().discriminant_preimage(target.stabilizer_of_element(discriminant_class))
    return lattice_.O_plus().intersection(preimage)


_DAWES = load_dawes_buildings()


def _dawes_groups(case):
    lattice_ = lattice(case["gram"])
    match case["id"]:
        case "dawes_2U_A2_stable" | "dawes_2U2_A2_stable":
            return [_stable_o_plus(lattice_)]
        case "dawes_2U_minus2_minus6":
            return [lattice_.O_plus(), _stable_o_plus(lattice_)]
        case "dawes_L2_gamma2":
            v = basis(lattice_)[4]
            assert v.q() == -2
            return [_o_plus_fixing(lattice_, v.divided_discriminant_class())]
    raise AssertionError(f"no group for {case['id']}")


@pytest.mark.parametrize("case", _DAWES["buildings"], ids=lambda case: case["id"])
@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", raises=UnfinishedCapability, strict=True)
def test_dawes_buildings_have_the_published_counts(case) -> None:
    require(ISOTROPIC_ORBITS, FLAG_ORBITS, SPLIT_ORBIT)

    building = case["building"]

    for group in _dawes_groups(case):
        assert _counts(group) == (building["points"], building["curves"], building["edges"])


@pytest.mark.parametrize("graph", _DAWES["incidence_graphs"], ids=lambda graph: graph["building_id"])
@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", raises=UnfinishedCapability, strict=True)
def test_dawes_building_incidence_graph_is_the_published_figure(graph) -> None:
    require(ISOTROPIC_ORBITS, FLAG_ORBITS, SPLIT_ORBIT)

    case = next(case for case in _DAWES["buildings"] if case["id"] == graph["building_id"])
    (group,) = _dawes_groups(case)
    lines, planes = list(group.cusps(1)), list(group.cusps(2))

    computed = _incidence_graph(
        range(len(lines)),
        range(len(planes)),
        [(lines.index(incidence.line_cusp()), planes.index(incidence.plane_cusp())) for incidence in group.tits_building_incidence()],
    )
    published = _incidence_graph(graph["points"], graph["curves"], graph["point_curve_incidences"])

    assert computed == published


@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", raises=UnfinishedCapability, strict=True)
def test_dawes_last_index_in_the_chain_is_the_discriminant_image_of_o_plus() -> None:
    require(SUBGROUP_GENERATORS)

    (chain,) = _DAWES["index_chains"]
    case = next(case for case in _DAWES["buildings"] if case["id"] == "dawes_2U_A2_stable")
    lattice_ = lattice(case["gram"])
    target = lattice_.discriminant_group().orthogonal_group()

    image = target.subgroup_on(tuple(g.discriminant_morphism() for g in lattice_.O_plus().group_generators()))

    assert int(image.cardinality()) == chain["indices"][-1]


def _incidence_graph(points, curves, incidences):
    """The point-curve incidence graph up to isomorphism preserving points and curves."""
    point_vertices = [("p", point) for point in points]
    curve_vertices = [("c", curve) for curve in curves]
    graph = Graph([point_vertices + curve_vertices, [(("p", point), ("c", curve)) for point, curve in incidences]], format="vertices_and_edges")
    return graph.canonical_label(partition=[point_vertices, curve_vertices]).copy(immutable=True)


_DTOWER = load_dtower_boundaries()


def _lambda_n(n: int) -> list[list[int]]:
    """U^2 + D_{N-2}, D negative definite, D_1 = <-4> (Laza--O'Grady line 625).

    Laza--O'Grady define D_m for m >= 3 and m = 1; for N = 4 the case uses D_2 = A_1 + A_1.
    """
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    match n - 2:
        case 1:
            d = matrix(ZZ, [[-4]])
        case 2:
            d = matrix(ZZ, [[-2, 0], [0, -2]])
        case rank:
            d = -matrix(ZZ, CartanMatrix(["D", rank]))
    return [list(row) for row in block_diagonal_matrix([hyperbolic, hyperbolic, d]).rows()]


def _gamma(lattice_):
    """Gamma(N) = {phi in O^+(Lambda_N) : phi(xi) = xi} for a decoration xi, q(xi) = 1 mod 2Z."""
    form = lattice_.discriminant_group()
    one = form.quadratic_value_module()(1)
    decorations = [x for x in form.elements() if form.q(x) == one]
    assert decorations, f"{lattice_} has no decoration"
    return _o_plus_fixing(lattice_, decorations[0])


@pytest.mark.parametrize("case", _DTOWER["type_counts"], ids=lambda case: f"N_{case['N']}")
@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", raises=UnfinishedCapability, strict=True)
def test_dtower_boundary_component_counts(case) -> None:
    require(ISOTROPIC_ORBITS, SPLIT_ORBIT)

    group = _gamma(lattice(_lambda_n(case["N"])))

    assert len(group.isotropic_orbit_representatives(1)) == case["type_iii_components"]
    assert len(group.isotropic_orbit_representatives(2)) == case["type_ii_components"]


@pytest.mark.parametrize("picture", _DTOWER["boundary_pictures"], ids=lambda picture: picture["id"])
@pytest.mark.xfail(reason="finite-index subgroups and orbit splitting: owned by #22", raises=UnfinishedCapability, strict=True)
def test_dtower_boundary_pictures(picture) -> None:
    require(ISOTROPIC_ORBITS, FLAG_ORBITS, SPLIT_ORBIT)

    group = _gamma(lattice(picture["gram"]))

    assert _counts(group) == (picture["type_iii_points"], picture["type_ii_curves"], len(picture["incidences"]))


_K3 = load_k3_modular_strata()


def _lambda_2k(k: int) -> list[list[int]]:
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    e8 = -matrix(ZZ, CartanMatrix(["E", 8]))
    return [list(row) for row in block_diagonal_matrix([hyperbolic, hyperbolic, e8, e8, matrix(ZZ, [[-2 * k]])]).rows()]


@pytest.mark.xfail(reason="isotropic sublattice and flag orbits: owned by #21", raises=UnfinishedCapability, strict=True)
def test_degree_two_k3_type_ii_components_and_root_types() -> None:
    require(ISOTROPIC_ORBITS)

    case = _K3["degree_two_polarized_k3"]
    group = lattice(_lambda_2k(1)).O()

    planes = group.isotropic_orbit_representatives(2)

    assert len(planes) == case["baily_borel_boundary"]["type_ii_curves"]
    assert sorted(plane.isotropic_reduction().norm_two_root_types() for plane in planes) == sorted(case["type_ii_root_types"])


@pytest.mark.xfail(reason="isotropic sublattice and flag orbits: owned by #21", raises=UnfinishedCapability, strict=True)
def test_degree_four_k3_has_nine_type_ii_components() -> None:
    require(ISOTROPIC_ORBITS)

    case = _K3["degree_four_polarized_k3"]
    group = lattice(_lambda_2k(2)).O()

    assert len(group.isotropic_orbit_representatives(2)) == case["baily_borel_boundary"]["type_ii_curves"]


_ENRIQUES = load_enriques_87_polarizations()
_E10 = load_e10_fundamental_domain()
_ORDER_O_PLUS_10_2 = 2**21 * 3**5 * 5**2 * 7 * 17 * 31


def _gamma_h(face_polarization):
    """Gamma_h^+ = O^+(N) cap pi_N^{-1}(pi_M(O(M, h))) for h in M(1/2) = U + E8(-1).

    The K3 lattice with the Enriques involution carries M = U(2) + E8(-2) as invariant and
    N as anti-invariant lattice. h is transported from Dutour Sikirić--Hulek's coordinates
    on U + E8(-1) into M along an isometry (U + E8(-1))(2) -> M.
    """
    from dzack_research.preamble.catalogue import Involutions

    involution = Involutions.I_En
    invariant = involution.primitive_extension().invariant
    transport = lattice(_E10["gram"]).twist(2).Isom(invariant).an_element()
    polarization = invariant.inclusion()(transport(element(transport.domain(), face_polarization)))
    anti_invariant_group = involution.polarized(polarization).coinvariant_extension_subgroup()
    anti_invariant = anti_invariant_group.supergroup().domain()
    return anti_invariant_group.intersection(anti_invariant.O_plus())


@pytest.mark.parametrize("case", _ENRIQUES, ids=lambda case: f"case_{case['case']}")
@pytest.mark.xfail(reason="EquivariantLattice and centralizers: owned by #23", raises=UnfinishedCapability, strict=True)
def test_enriques_modular_groups_match_dutour_sikiric_hulek(case) -> None:
    require(EQUIVARIANT_LATTICE, ISOTROPIC_ORBITS, FLAG_ORBITS, SPLIT_ORBIT)

    group = _gamma_h(case["face_polarization"])
    anti_invariant = group.supergroup().domain()
    target = anti_invariant.discriminant_group().orthogonal_group()

    image = target.subgroup_on(tuple(g.discriminant_morphism() for g in group.group_generators()))

    assert int(image.cardinality()) == case["group_order"]
    assert _counts(group) == (case["line_orbits"], case["plane_orbits"], case["flag_orbits"])


@pytest.mark.parametrize(
    "case",
    load_enriques_classical_indices()["consistent_with_table"] + load_enriques_classical_indices()["inconsistent_with_table"],
    ids=lambda case: f"degree_{case['degree']}",
)
@pytest.mark.xfail(reason="EquivariantLattice and centralizers: owned by #23", raises=UnfinishedCapability, strict=True)
def test_classical_polarization_indices(case) -> None:
    require(EQUIVARIANT_LATTICE, SUBGROUP_GENERATORS)

    row = next(row for row in _ENRIQUES if row["case"] == case["case"])
    group = _gamma_h(row["face_polarization"])
    anti_invariant = group.supergroup().domain()
    target = anti_invariant.discriminant_group().orthogonal_group()

    image = target.subgroup_on(tuple(g.discriminant_morphism() for g in group.group_generators()))

    assert int(image.cardinality()) == case["gamma_h_over_stable_order"]
    assert _ORDER_O_PLUS_10_2 // int(image.cardinality()) == case["index_in_O_plus_N"]
