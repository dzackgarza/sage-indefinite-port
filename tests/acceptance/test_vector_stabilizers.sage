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
from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction

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
    lattice_ = lattice([[0, 2], [2, 1]])
    source = element(lattice_, [1, 0])
    target = element(lattice_, [1, -4])
    assert int(source.q()) == int(target.q()) == 0
    assert int(source.content()) == int(target.content()) == 1
    assert int(source.div()) == int(target.div()) == 2

    witness = EichlerRecursion().pointed_isometry(source, target)
    assert witness is not None
    assert witness(source) == target
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_odd_lattice_isotropic_pointed_isometry_uses_general_divisibility_two_route() -> None:
    """The generic isotropic route also handles odd lattices and divisibility two."""
    lattice_ = lattice([[0, 2], [2, 1]])
    source = element(lattice_, [1, 0])
    target = element(lattice_, [1, -4])
    assert not lattice_.is_even()
    assert int(source.q()) == int(target.q()) == 0
    assert int(source.content()) == int(target.content()) == 1
    assert int(source.div()) == int(target.div()) == 2

    witness = EichlerRecursion().pointed_isometry(source, target)
    assert witness is not None
    assert witness(source) == target
    assert witness(element(lattice_, [0, 1])) == -element(lattice_, [0, 1])
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_general_isotropic_pointed_isometry_carries_negative_marking() -> None:
    """The general (divisibility-two) route must normalize against the actual marked vector."""
    lattice_ = lattice([[0, 2, 0, 0], [2, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    source = element(lattice_, [1, 0, 0, 0])
    target = -source

    assert int(source.q()) == int(target.q()) == 0
    assert int(source.content()) == int(target.content()) == 1
    assert int(source.div()) == int(target.div()) == 2
    witness = EichlerRecursion().pointed_isometry(source, target)
    assert witness is not None
    assert witness(source) == target
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_empty_integral_fibre_falls_back_to_full_marked_vector_stabilizer(monkeypatch) -> None:
    """An empty lift fibre for one reduction map must not imply pointed Isom is empty."""
    require(MARKED_ISOMETRY)
    lattice_ = lattice([[0, 2, 0, 0], [2, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    source = element(lattice_, [1, 0, 0, 0])
    target = element(lattice_, [0, 0, 1, 0])
    recursion = EichlerRecursion()
    source_reduction = lattice_.primitive_isotropic_subobject(source).isotropic_reduction()
    reduction_type = type(source_reduction)
    correction_calls = []
    original_correction = recursion._isotropic_correction

    class EmptyIntegralFibre:
        def is_empty(self):
            return True

    class EmptyTorsor:
        def integral_members(self):
            return EmptyIntegralFibre()

    def no_integral_lifts(self, target_reduction, reduced_isometry):
        return EmptyTorsor()

    def record_correction(*args, **kwargs):
        correction_calls.append(True)
        return original_correction(*args, **kwargs)

    monkeypatch.setattr(reduction_type, "rational_lifts", no_integral_lifts)
    monkeypatch.setattr(recursion, "_isotropic_correction", record_correction)
    witness = recursion.pointed_isometry(source, target)

    assert correction_calls
    assert witness is not None
    assert witness(source) == target
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_isotropic_pointed_isometry_transports_between_changed_basis_parents() -> None:
    """The recursive route returns an owned witness for distinct, isometric lattice parents."""
    gram = matrix(ZZ, [[0, 2, 0, 0], [2, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]])
    change = matrix(ZZ, [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    source_lattice = lattice(gram)
    target_lattice = lattice(change.transpose() * gram * change)
    source = element(source_lattice, [1, 0, 0, 0])
    target = element(target_lattice, [1, 0, 0, 0])

    assert int(source.q()) == int(target.q()) == 0
    assert int(source.div()) == int(target.div()) == 2
    witness = EichlerRecursion().pointed_isometry(source, target)
    assert witness is not None
    assert witness.domain() is source_lattice
    assert witness.codomain() is target_lattice
    assert witness(source) == target
    source_basis = tuple(source_lattice.module_generators())
    target_basis = tuple(target_lattice.module_generators())
    assert all(
        target_lattice.b(witness(left), witness(right)) == source_lattice.b(left, right)
        for left in source_basis
        for right in source_basis
    )


def test_split_divisibility_one_pointed_isometry_carries_actual_negative_marking() -> None:
    lattice_ = lattice([[0, 1], [1, 0]])
    source = element(lattice_, [1, 0])
    target = -source
    assert int(source.q()) == int(target.q()) == 0
    assert int(source.content()) == int(target.content()) == 1
    assert int(source.div()) == int(target.div()) == 1

    witness = EichlerRecursion().pointed_isometry(source, target)
    assert witness is not None
    assert witness(source) == target
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_isotropic_pointed_isometry_rejects_different_divisibility() -> None:
    """Divisibility is preserved by every integral isometry of a lattice."""
    lattice_ = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, -8]])
    source = element(lattice_, [1, 0, 0, 0])
    target = element(lattice_, [0, 2, 2, 1])

    assert int(source.q()) == int(target.q()) == 0
    assert int(source.content()) == int(target.content()) == 1
    assert int(source.div()) == 1
    assert int(target.div()) == 2
    assert EichlerRecursion().pointed_isometry(source, target) is None


def test_isotropic_plane_isometry_carries_changed_basis_plane() -> None:
    """Plane equivalence is setwise: a unimodular basis change does not change the plane."""
    lattice_ = lattice([[0, 4, 0, 0], [4, 0, 0, 0], [0, 0, 0, 4], [0, 0, 4, 0]])
    first = element(lattice_, [1, 0, 0, 0])
    second = element(lattice_, [0, 0, 1, 0])
    changed_first = first + second
    source_plane = lattice_.primitive_isotropic_subobject(first, second)
    target_plane = lattice_.primitive_isotropic_subobject(changed_first, second)

    assert lattice_.b(first, second) == 0
    assert lattice_.b(changed_first, second) == 0
    witness = IsotropicInduction(EichlerRecursion()).plane_isometry(source_plane, target_plane)
    assert witness is not None
    target_basis_matrix = matrix(ZZ, [coordinates(changed_first), coordinates(second)]).transpose()
    first_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(first))))
    second_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(second))))
    induced_basis_change = matrix(ZZ, [
        [first_coordinates[0], second_coordinates[0]],
        [first_coordinates[1], second_coordinates[1]],
    ])
    assert abs(induced_basis_change.det()) == 1
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_isotropic_plane_isometry_transports_between_changed_basis_parents() -> None:
    """The plane-isometry witness retains distinct source and target lattice parents."""
    gram = matrix(ZZ, [[0, 4, 0, 0], [4, 0, 0, 0], [0, 0, 0, 4], [0, 0, 4, 0]])
    change = matrix(ZZ, [[1, 0, 0, 0], [0, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 1]])
    source_lattice = lattice(gram)
    target_lattice = lattice(change.transpose() * gram * change)
    source_first = element(source_lattice, [1, 0, 0, 0])
    source_second = element(source_lattice, [0, 0, 1, 0])
    target_first = element(target_lattice, [1, 0, 0, 0])
    target_second = element(target_lattice, [0, 0, 1, 0])
    source_plane = source_lattice.primitive_isotropic_subobject(source_first, source_second)
    target_plane = target_lattice.primitive_isotropic_subobject(target_first, target_second)

    witness = IsotropicInduction(EichlerRecursion()).plane_isometry(source_plane, target_plane)
    assert witness is not None
    assert witness.domain() is source_lattice
    assert witness.codomain() is target_lattice
    target_basis_matrix = matrix(ZZ, [coordinates(target_first), coordinates(target_second)]).transpose()
    first_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(source_first))))
    second_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(source_second))))
    induced_basis_change = matrix(ZZ, [
        [first_coordinates[0], second_coordinates[0]],
        [first_coordinates[1], second_coordinates[1]],
    ])
    assert abs(induced_basis_change.det()) == 1
    source_basis = tuple(source_lattice.module_generators())
    assert all(
        target_lattice.b(witness(left), witness(right)) == source_lattice.b(left, right)
        for left in source_basis
        for right in source_basis
    )


def test_isotropic_flag_isometry_preserves_line_inside_plane() -> None:
    """A flag transporter preserves the marked line, not merely the ambient plane."""
    lattice_ = lattice([[0, 4, 0, 0], [4, 0, 0, 0], [0, 0, 0, 4], [0, 0, 4, 0]])
    first = element(lattice_, [1, 0, 0, 0])
    second = element(lattice_, [0, 0, 1, 0])
    changed_line = first + second
    plane = lattice_.primitive_isotropic_subobject(first, second)
    induction = IsotropicInduction(EichlerRecursion())

    witness = induction.flag_isometry((first, plane), (changed_line, plane))
    assert witness is not None
    assert witness(first) in (changed_line, -changed_line)
    target_basis = induction.flag_basis((changed_line, plane))
    target_basis_matrix = matrix(ZZ, [coordinates(target_basis[0]), coordinates(target_basis[1])]).transpose()
    first_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(first))))
    second_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(second))))
    induced_basis_change = matrix(ZZ, [
        [first_coordinates[0], second_coordinates[0]],
        [first_coordinates[1], second_coordinates[1]],
    ])
    assert abs(induced_basis_change.det()) == 1
    basis = tuple(lattice_.module_generators())
    assert all(lattice_.b(witness(left), witness(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_isotropic_flag_isometry_transports_between_changed_basis_parents() -> None:
    """The flag transporter preserves both the marked line and plane across presentations."""
    gram = matrix(ZZ, [[0, 4, 0, 0], [4, 0, 0, 0], [0, 0, 0, 4], [0, 0, 4, 0]])
    change = matrix(ZZ, [[1, 0, 0, 0], [0, 1, 0, 0], [1, 0, 1, 0], [0, 0, 0, 1]])
    source_lattice = lattice(gram)
    target_lattice = lattice(change.transpose() * gram * change)
    source_first = element(source_lattice, [1, 0, 0, 0])
    source_second = element(source_lattice, [0, 0, 1, 0])
    target_first = element(target_lattice, [1, 0, 0, 0])
    target_second = element(target_lattice, [0, 0, 1, 0])
    source_plane = source_lattice.primitive_isotropic_subobject(source_first, source_second)
    target_plane = target_lattice.primitive_isotropic_subobject(target_first, target_second)
    induction = IsotropicInduction(EichlerRecursion())

    witness = induction.flag_isometry((source_first, source_plane), (target_first, target_plane))
    assert witness is not None
    assert witness.domain() is source_lattice
    assert witness.codomain() is target_lattice
    assert witness(source_first) in (target_first, -target_first)
    target_basis = induction.flag_basis((target_first, target_plane))
    target_basis_matrix = matrix(ZZ, [coordinates(target_basis[0]), coordinates(target_basis[1])]).transpose()
    first_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(source_first))))
    second_coordinates = target_basis_matrix.solve_right(vector(ZZ, coordinates(witness(source_second))))
    induced_basis_change = matrix(ZZ, [
        [first_coordinates[0], second_coordinates[0]],
        [first_coordinates[1], second_coordinates[1]],
    ])
    assert abs(induced_basis_change.det()) == 1
    source_basis = tuple(source_lattice.module_generators())
    assert all(
        target_lattice.b(witness(left), witness(right)) == source_lattice.b(left, right)
        for left in source_basis
        for right in source_basis
    )


def test_isotropic_plane_and_flag_automorphism_generators_preserve_their_marks() -> None:
    """The plane stabilizer preserves J setwise; the flag stabilizer also preserves Zu."""
    lattice_ = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    first = element(lattice_, [1, 0, 0, 0])
    second = element(lattice_, [0, 0, 1, 0])
    plane = lattice_.primitive_isotropic_subobject(first, second)
    flag_line = first + second
    induction = IsotropicInduction(EichlerRecursion())
    plane_inclusion = plane.inclusion()
    plane_basis = tuple(plane_inclusion(generator) for generator in plane_inclusion.domain().module_generators())
    ambient_basis = tuple(lattice_.module_generators())

    line_generators = induction.line_automorphism_generators(lattice_, first)
    assert line_generators
    for generator in line_generators:
        assert generator(first) in (first, -first)
        assert all(
            lattice_.b(generator(left), generator(right)) == lattice_.b(left, right)
            for left in ambient_basis
            for right in ambient_basis
        )

    for generators, preserve_line, expected_mod_three_order in (
        (induction.plane_automorphism_generators(plane), False, 48),
        (induction.flag_automorphism_generators((flag_line, plane)), True, 12),
    ):
        assert generators
        induced_actions = []
        for generator in generators:
            if preserve_line:
                assert generator(flag_line) in (flag_line, -flag_line)
            images = tuple(generator(point) for point in plane_basis)
            image_coordinates = tuple(plane_inclusion.lift(point).to_vector() for point in images)
            coordinate_matrix = matrix(ZZ, [
                [coordinate_vector(label) for label in plane_inclusion.domain().module_generating_set()]
                for coordinate_vector in image_coordinates
            ])
            assert abs(coordinate_matrix.det()) == 1
            induced_actions.append(matrix(GF(3), coordinate_matrix))
            assert all(
                lattice_.b(generator(left), generator(right)) == lattice_.b(left, right)
                for left in ambient_basis
                for right in ambient_basis
            )
        # The plane action is all of GL(2,Z); the flag action is its line
        # stabilizer. Their exact finite images have orders 48 and 12 over F_3.
        assert MatrixGroup(induced_actions).order() == expected_mod_three_order


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
    first = element(lattice_, [1, 0, 0, 0])
    second = element(lattice_, [0, 0, 1, 0])
    for plane in ((first, second), (first + second, second)):
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
