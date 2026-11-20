"""Independently sourced orthogonal-group properties for indefinite lattices.

Recorded outputs come from Hecke's QuadBin tests and OSCAR's spinor-norm tests.
Finite reduction *orders* alone do not establish equality of finite image
subgroups, let alone completeness of generators of an infinite integral
orthogonal group. The former finite-order matching tests have been removed
as mathematical acceptance assertions.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_AUTOMORPHISMS,
    discriminant_image_order,
    is_isometry,
    lattice,
    require,
    rows_of,
)
from tests.fixtures.oracle_fixtures import (
    load_binary_form_automorphisms,
    load_oscar_lattice_oracles,
)


@pytest.mark.parametrize(
    "case",
    load_binary_form_automorphisms()["improper_automorphisms"],
    ids=lambda case: case["id"],
)
def test_ambiguous_binary_form_has_an_improper_automorphism(case) -> None:
    require(MARKED_AUTOMORPHISMS)

    computed = lattice(case["gram"]).O().group_generators()

    assert any(int(generator.determinant()) == -1 for generator in computed)


def test_degenerate_orthogonal_group_contains_radical_shear_and_gl_action() -> None:
    """For <2> plus a rank-one radical, O(L) includes both signs and every shear."""
    require(MARKED_AUTOMORPHISMS)

    gram = [[2, 0], [0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(tuple(tuple(row) for row in rows_of(generator)) for generator in computed)
    expected = {
        ((-1, 0), (0, 1)),
        ((1, 0), (0, -1)),
        ((1, 1), (0, 1)),
    }

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert expected <= set(rows)


def test_degenerate_orthogonal_group_respects_permuted_radical_coordinates() -> None:
    """The radical extension generators are expressed in the declared coordinate frame."""
    from sage_extensions.all import ZZ as PreambleZZ
    from sage_extensions.all import Lattices

    require(MARKED_AUTOMORPHISMS)
    labels = ("radical", "nonsingular")
    gram = [[0, 0], [0, 2]]
    lattice_ = Lattices(PreambleZZ)(gram, module_generators=labels)
    assert tuple(lattice_.module_generating_set()) == labels

    from sage_indefinite_port.indefinite.recursive import _column_matrix

    computed = lattice_.O().group_generators()
    rows = tuple(
        tuple(tuple(int(entry) for entry in row) for row in _column_matrix(generator).transpose().rows())
        for generator in computed
    )
    expected = {
        ((-1, 0), (0, 1)),
        ((1, 0), (0, -1)),
        ((1, 0), (1, 1)),
    }

    assert all(is_isometry(rows_of(generator), gram) for generator in computed)
    assert expected <= set(rows)


def test_degenerate_orthogonal_group_generates_rank_two_radical_action() -> None:
    """The radical GL_2(Z) action and the nonsingular-to-radical shear survive assembly."""
    require(MARKED_AUTOMORPHISMS)

    gram = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(tuple(tuple(row) for row in rows_of(generator)) for generator in computed)

    assert all(is_isometry(rows_of(generator), gram) for generator in computed)
    assert any(generator_rows[0][1] != 0 for generator_rows in rows)
    assert any(generator_rows[1][2] != 0 or generator_rows[2][1] != 0 for generator_rows in rows)
    assert any(
        tuple(tuple(generator_rows[i][j] for j in (1, 2)) for i in (1, 2)) != ((1, 0), (0, 1))
        for generator_rows in rows
    )


def test_degenerate_orthogonal_group_generates_permuted_rank_two_radical_action() -> None:
    """Rank-two radical GL and shears survive when the nonsingular coordinate is last."""
    from sage_extensions.all import ZZ as PreambleZZ
    from sage_extensions.all import Lattices

    require(MARKED_AUTOMORPHISMS)
    labels = ("radical_b", "radical_a", "nonsingular")
    gram = [[0, 0, 0], [0, 0, 0], [0, 0, 2]]
    lattice_ = Lattices(PreambleZZ)(gram, module_generators=labels)
    assert tuple(lattice_.module_generating_set()) == labels

    from sage_indefinite_port.indefinite.recursive import _column_matrix

    computed = lattice_.O().group_generators()
    rows = tuple(
        tuple(tuple(int(entry) for entry in row) for row in _column_matrix(generator).transpose().rows())
        for generator in computed
    )

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert ((1, 0, 0), (0, 1, 0), (1, 0, 1)) in set(rows)
    assert any(generator_rows[2][0] != 0 or generator_rows[2][1] != 0 for generator_rows in rows)
    assert any(
        tuple(tuple(generator_rows[i][j] for j in (0, 1)) for i in (0, 1)) != ((1, 0), (0, 1))
        for generator_rows in rows
    )
    assert any(generator_rows[2][2] == -1 for generator_rows in rows)


def test_degenerate_orthogonal_group_embeds_non_diagonal_complement() -> None:
    """A non-diagonal nonsingular summand embeds in the radical extension's coordinate frame."""
    require(MARKED_AUTOMORPHISMS)

    gram = [[2, 1, 0], [1, 2, 0], [0, 0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(rows_of(generator) for generator in computed)

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert any(generator_rows[0][2] != 0 for generator_rows in rows)
    assert any(generator_rows[1][2] != 0 for generator_rows in rows)
    assert any(generator_rows[2][2] == -1 for generator_rows in rows)


def test_zero_form_orthogonal_group_is_the_full_general_linear_group() -> None:
    """When the whole lattice is radical, assembly is exactly GL_2(Z)."""
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    require(MARKED_AUTOMORPHISMS)

    gram = [[0, 0], [0, 0]]
    computed = lattice(gram).O().group_generators()
    rows = tuple(tuple(tuple(row) for row in rows_of(generator)) for generator in computed)

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert all(abs(int(matrix(SageZZ, generator_rows).determinant())) == 1 for generator_rows in rows)
    assert any(generator_rows[0][1] != 0 or generator_rows[1][0] != 0 for generator_rows in rows)
    assert any(
        generator_rows in {((1, 1), (0, 1)), ((1, 0), (1, 1))}
        for generator_rows in rows
    )
    assert any(
        generator_rows in {
            ((0, 1), (1, 0)),
            ((0, -1), (1, 0)),
            ((0, 1), (-1, 0)),
            ((0, -1), (-1, 0)),
        }
        for generator_rows in rows
    )
    assert any(int(generator.determinant()) == -1 for generator in computed)


def test_zero_form_orthogonal_group_uses_permuted_coordinate_labels() -> None:
    """The full GL_2(Z) action is returned in the declared radical coordinate order."""
    from sage_extensions.all import ZZ as PreambleZZ
    from sage_extensions.all import Lattices
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    require(MARKED_AUTOMORPHISMS)
    labels = ("z_second", "a_first")
    gram = [[0, 0], [0, 0]]
    lattice_ = Lattices(PreambleZZ)(gram, module_generators=labels)
    assert tuple(lattice_.module_generating_set()) == labels

    from sage_indefinite_port.indefinite.recursive import _column_matrix

    computed = lattice_.O().group_generators()
    rows = tuple(
        tuple(tuple(int(entry) for entry in row) for row in _column_matrix(generator).transpose().rows())
        for generator in computed
    )

    assert all(is_isometry(generator_rows, gram) for generator_rows in rows)
    assert all(abs(int(matrix(SageZZ, generator_rows).determinant())) == 1 for generator_rows in rows)
    assert any(generator_rows[0][1] != 0 or generator_rows[1][0] != 0 for generator_rows in rows)
    assert any(int(generator.determinant()) == -1 for generator in computed)


def test_orthogonal_group_of_the_hyperbolic_plane_has_order_four() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["isometry_groups"] if case["id"] == "O_U_order")

    assert lattice(case["gram"]).O().order() == case["orthogonal_group_order"]


_IMAGE_CASES = [case for case in load_oscar_lattice_oracles()["discriminant_images"] if "O_qL_order" in case]


@pytest.mark.parametrize("case", _IMAGE_CASES, ids=lambda case: case["id"])
def test_discriminant_images_match_oscar(case) -> None:
    """|O(q_L)|, the image of O(L), and the image of O^+(L) = ker sn_R (OSCAR's image_in_Oq_signed)."""
    require(MARKED_AUTOMORPHISMS)

    lattice_ = lattice(case["gram"])

    assert int(lattice_.discriminant_group().orthogonal_group().cardinality()) == case["O_qL_order"]
    assert int(lattice_.O().discriminant_image().cardinality()) == case["image_in_Oq_order"]
    assert discriminant_image_order(lattice_, lattice_.O_plus().group_generators()) == case["image_in_Oq_signed_order"]


def test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["discriminant_images"] if case["id"] == "A7_plus_diag_1_1_-1")
    lattice_ = lattice(case["gram"])

    assert lattice_.O().discriminant_image().cardinality() == lattice_.discriminant_group().orthogonal_group().cardinality()


def test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q() -> None:
    require(MARKED_AUTOMORPHISMS)

    case = next(case for case in load_oscar_lattice_oracles()["discriminant_images"] if case["id"] == "U_plus_-3")
    lattice_ = lattice(case["gram"])
    signed = discriminant_image_order(lattice_, lattice_.O_plus().group_generators())

    assert signed == case["image_in_Oq_signed_order"]
    assert signed == int(lattice_.discriminant_group().orthogonal_group().cardinality())
