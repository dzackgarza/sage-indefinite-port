import pytest

from dzack_research.preamble.all import Lattices, ZZ

from sage_indefinite_port.indefinite.recursive import IndefiniteOrthogonalAlgorithm
from tests.fixtures.oracle_fixtures import (
    load_ci_indefinite_comp,
    load_conway_sloane_cases,
    load_unpolarized_enriques,
)

def _fixture_component(name):
    two = ZZ.one() + ZZ.one()
    match name:
        case "U":
            return Lattices(ZZ)("U")
        case "2U":
            return Lattices(ZZ)("U").twist(two)
        case "2E8":
            return Lattices(ZZ)("E8").twist(two)
        case _:
            return Lattices(ZZ)(name)


def _ci_lattice(case, *, changed_first_u=False):
    components = tuple(case["components"])
    assert components and components[0] == "U"
    match changed_first_u:
        case True:
            lattice = Lattices(ZZ)([[2, 1], [1, 0]])
        case False:
            lattice = _fixture_component("U")
    for component in components[1:]:
        lattice = lattice + _fixture_component(component)
    return lattice


def _discriminant_image(lattice, generators):
    target = lattice.discriminant_group().orthogonal_group()
    return target.subgroup_on(
        tuple(generator.discriminant_morphism() for generator in generators)
    )


def test_conway_sloane_51a_51b_public_isometry_rejects() -> None:
    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first = Lattices(ZZ)(pair["form_a"]["gram"])
    second = Lattices(ZZ)(pair["form_b"]["gram"])

    witness = IndefiniteOrthogonalAlgorithm().isometry(first, second)

    assert pair["integrally_equivalent"] is False
    assert witness is None


@pytest.mark.parametrize(
    "case",
    load_ci_indefinite_comp(),
    ids=lambda case: case["id"],
)
@pytest.mark.parametrize(
    "square",
    (ZZ.zero(), ZZ.one() + ZZ.one()),
    ids=("q0", "q2"),
)
def test_ci_indefinite_comp_public_basis_covariance_and_orbit_counts(case, square) -> None:
    lattice = _ci_lattice(case)
    changed = _ci_lattice(case, changed_first_u=True)
    algorithm = IndefiniteOrthogonalAlgorithm()

    witness = algorithm.isometry(lattice, changed)
    assert witness is not None
    assert witness.domain() is lattice
    assert witness.codomain() is changed

    source_orbits = algorithm.vector_orbit_representatives(lattice, square)
    target_orbits = algorithm.vector_orbit_representatives(changed, square)
    assert len(source_orbits) == len(target_orbits)


def test_enriques_public_discriminant_stabilizer_indices() -> None:
    fixture = load_unpolarized_enriques()
    lattice = Lattices(ZZ)(fixture["lattice"]["gram"])
    algorithm = IndefiniteOrthogonalAlgorithm()
    full_group = algorithm.orthogonal_group(lattice)
    full_image = _discriminant_image(lattice, full_group.generators())
    minus_identity = lattice.Aut()(
        tuple(-generator for generator in lattice.module_generators())
    )
    labels = tuple(lattice.module_generating_set())

    for representative in fixture["component_preserving_group"]["line_representatives"]:
        coordinates = representative["basis"][0]
        vector = lattice.linear_combination(
            {
                label: lattice.base_ring()(coefficient)
                for label, coefficient in zip(labels, coordinates, strict=True)
                if coefficient
            }
        )
        assert vector.q() == lattice.base_ring().zero()
        stabilizer = algorithm.vector_stabilizer(vector)
        stabilizer_image = _discriminant_image(
            lattice,
            (*stabilizer.generators(), minus_identity),
        )
        index = int(full_image.cardinality()) // int(stabilizer_image.cardinality())
        assert index == representative["discriminant_stabilizer_index"]
