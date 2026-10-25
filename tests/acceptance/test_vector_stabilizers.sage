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

from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    MARKED_CLASSES,
    MARKED_ISOMETRY,
    automorphisms,
    coordinates,
    element,
    is_isometry,
    isomorphisms,
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
from sage_indefinite_port.indefinite.recursive import EichlerRecursion

_E10 = load_e10_fundamental_domain()
_ORDER_O_PLUS_10_2 = 2**21 * 3**5 * 5**2 * 7 * 17 * 31


def test_isotropic_pointed_automorphism_generators_preserve_marking_and_form() -> None:
    """Exercise the rank-one isotropic Witt-extension and kernel assembly on 2U."""
    lattice_ = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    marked = element(lattice_, [1, 0, 0, 0])
    assert int(marked.q()) == 0
    assert int(marked.content()) == 1

    recursion = EichlerRecursion()
    line_generators = recursion.line_stabilizer_generators(marked)
    pointed_generators = recursion.pointed_automorphism_generators(marked)
    basis = tuple(lattice_.module_generators())
    assert line_generators and pointed_generators
    for generator in line_generators:
        assert generator(marked) in (marked, -marked)
        assert all(
            lattice_.b(generator(left), generator(right)) == lattice_.b(left, right)
            for left in basis
            for right in basis
        )
    for generator in pointed_generators:
        assert generator(marked) == marked
        assert all(
            lattice_.b(generator(left), generator(right)) == lattice_.b(left, right)
            for left in basis
            for right in basis
        )


def test_divisibility_two_isotropic_pointed_isometry_has_exact_witness() -> None:
    """Exercise the general isotropic pointed-isometry route beyond the div-one split case."""
    lattice_ = lattice([[0, 2, 0, 0], [2, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    source = element(lattice_, [1, 0, 0, 0])
    target = element(lattice_, [0, 0, 1, 0])
    assert int(source.q()) == int(target.q()) == 0
    assert int(source.content()) == int(target.content()) == 1
    assert int(source.div()) == int(target.div()) == 2

    witness = EichlerRecursion().pointed_isometry(source, target)
    assert witness is not None
    assert witness(source) == target
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_divisibility_four_isotropic_stabilizer_constructs_exact_containing_bound() -> None:
    """The pairing L -> (Zv)^* has finite cokernel of order four in this example."""
    lattice_ = lattice([[0, 4, 0], [4, 0, 0], [0, 0, -8]])
    marked = element(lattice_, [1, 0, 0])
    assert int(marked.q()) == 0
    assert int(marked.content()) == 1
    assert int(marked.div()) == 4

    action = EichlerRecursion().isotropic_stabilizer_action((marked,), ())
    assert action.stable_containing_bound() is not None


def test_divisibility_four_isotropic_plane_stabilizer_constructs_exact_containing_bound() -> None:
    """The rank-two pairing L -> J* has cokernel order sixteen in U(4) + U(4)."""
    lattice_ = lattice(
        [
            [0, 4, 0, 0],
            [4, 0, 0, 0],
            [0, 0, 0, 4],
            [0, 0, 4, 0],
        ]
    )
    plane = (element(lattice_, [1, 0, 0, 0]), element(lattice_, [0, 0, 1, 0]))
    assert all(int(vector.q()) == 0 and int(vector.content()) == 1 and int(vector.div()) == 4 for vector in plane)
    assert lattice_.b(plane[0], plane[1]) == 0

    action = EichlerRecursion().isotropic_stabilizer_action(plane, ())
    assert action.stable_containing_bound() is not None


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

    computed = Counter(_ORDER_O_PLUS_10_2 // reduction_image_order(_automorphism_rows(marked), 2) for marked in vector_classes(lattice(_E10["gram"]), square))

    assert computed == recorded


def _reflection_rows(gram, root) -> list[list[int]]:
    """Row matrix of the reflection in a norm -2 root: v -> v + b(v, r) r."""
    form = matrix(ZZ, gram)
    r = vector(ZZ, root)
    return [list(basis_vector + (basis_vector * form * r) * r) for basis_vector in identity_matrix(ZZ, len(root)).rows()]


@pytest.mark.parametrize("root", _E10["simple_roots"], ids=lambda root: f"r_{root['label']}")
def test_equivalence_witness_carries_a_representative_to_its_reflection(root) -> None:
    require(MARKED_ISOMETRY)

    gram = _E10["gram"]
    lattice_ = lattice(gram)
    reflection = matrix(ZZ, _reflection_rows(gram, root["vector"]))
    assert is_isometry(reflection.rows(), gram)
    for representative in _E10["orbit_representatives"]:
        source = vector(ZZ, representative["vector"])
        target = source * reflection
        source_vector = element(lattice_, source)

        isometry = underlying_isometry(
            isomorphisms(pointed_category(source_vector.q()), pointed(lattice_, source_vector), pointed(lattice_, element(lattice_, target))).an_element()
        )

        assert is_isometry(rows_of(isometry), gram)
        assert vector(ZZ, coordinates(isometry(source_vector))) == target
