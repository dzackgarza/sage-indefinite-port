"""Isometry decisions and witnesses between indefinite lattices.

Verdicts come from certified corpus pairs (Indefinite.jl and Hecke tests, polyhedral_common's
145 equivalent Lorentzian pairs, the Conway--Sloane same-genus pair). A returned witness is
accepted when it carries one Gram matrix to the other. polyhedral_common's 19_IndefiniteComp
suite is adapted as it is: a lattice and its image under a unimodular change of basis are
isometric.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_ISOMETRY,
    is_witness,
    lattice,
    require,
    rows_of,
)
from tests.fixtures.oracle_fixtures import (
    load_ci_indefinite_comp,
    load_conway_sloane_cases,
    load_indefinite_isometry_pairs,
    load_lorentzian_equivalence_cases,
)


def _assert_isometric(source_gram, target_gram) -> None:
    source, target = lattice(source_gram), lattice(target_gram)

    assert source.is_isometric(target)
    witness = source.Isom(target).an_element()
    assert is_witness(rows_of(witness), source_gram, target_gram)


@pytest.mark.parametrize("case", load_indefinite_isometry_pairs(), ids=lambda case: case["id"])
def test_isometry_verdicts_match_certified_pairs(case) -> None:
    require(MARKED_ISOMETRY)

    if case["isometric"]:
        _assert_isometric(case["gram1"], case["gram2"])
    else:
        assert not lattice(case["gram1"]).is_isometric(lattice(case["gram2"]))
        assert lattice(case["gram1"]).Isom(lattice(case["gram2"])).is_empty()


@pytest.mark.parametrize("case", load_lorentzian_equivalence_cases(), ids=lambda case: case["id"])
def test_lorentzian_pairs_are_isometric(case) -> None:
    require(MARKED_ISOMETRY)

    assert is_witness(case["transporter_witness"], case["mat2"], case["mat1"])

    _assert_isometric(case["mat1"], case["mat2"])


def test_row_matrix_morphism_roundtrip_uses_coordinate_labels() -> None:
    """The approximate-model row convention is translated at the owned O() boundary."""
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    from sage_indefinite_port.indefinite.eichler import _isometry_from_row_matrix

    lattice_ = lattice([[2, 1], [1, 2]])
    labels = tuple(lattice_.module_generating_set())
    rows = matrix(SageZZ, [[1, 0], [1, 1]])
    morphism = _isometry_from_row_matrix(lattice_, rows)

    for generator in lattice_.module_generators():
        coordinates = generator.to_vector()
        source_coordinates = tuple(int(coordinates(label)) for label in labels)
        expected = tuple(sum(source_coordinates[row] * rows[row, column] for row in range(len(labels))) for column in range(len(labels)))
        actual = morphism(generator).to_vector()
        assert tuple(int(actual(label)) for label in labels) == expected


def test_column_matrix_morphism_roundtrip_uses_coordinate_labels() -> None:
    """Matrix conventions survive different source/target generator iteration orders."""
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    from sage_indefinite_port.indefinite.eichler import (
        _isometry_from_row_matrix as eichler_row_isometry,
    )
    from sage_indefinite_port.indefinite.recursive import (
        _column_matrix,
        _isometry_from_column_matrix,
        _isometry_from_row_matrix,
    )

    source_labels = ("right", "left")
    target_labels = ("left", "right")
    gram = [[2, 1], [1, 2]]
    source = Lattices(SageZZ)(gram, module_generators=source_labels)
    target = Lattices(SageZZ)(gram, module_generators=target_labels)
    assert tuple(source.module_generating_set()) == source_labels
    assert tuple(target.module_generating_set()) == target_labels

    def generator_coordinate_order(lattice_, labels):
        return tuple(
            next(index for index, label in enumerate(labels) if generator.to_vector()(label))
            for generator in lattice_.module_generators()
        )

    assert generator_coordinate_order(source, source_labels) != generator_coordinate_order(target, target_labels)

    columns = matrix(SageZZ, [[-1, 0], [1, 1]])
    column_morphism = _isometry_from_column_matrix(source, target, columns)
    assert is_witness(_column_matrix(column_morphism).transpose().rows(), gram, gram)
    assert _column_matrix(column_morphism) == columns
    for generator in source.module_generators():
        coordinates = generator.to_vector()
        expected = tuple(
            sum(
                columns[row, column] * int(coordinates(source_labels[column]))
                for column in range(len(source_labels))
            )
            for row in range(len(source_labels))
        )
        actual = column_morphism(generator).to_vector()
        assert tuple(int(actual(label)) for label in target_labels) == expected

    row_matrix = matrix(SageZZ, [[-1, 1], [0, 1]])
    row_morphism = _isometry_from_row_matrix(source, target, row_matrix)
    assert is_witness(_column_matrix(row_morphism).transpose().rows(), gram, gram)
    for generator in source.module_generators():
        coordinates = generator.to_vector()
        source_coordinates = tuple(int(coordinates(label)) for label in source_labels)
        expected = tuple(
            sum(source_coordinates[row] * row_matrix[row, column] for row in range(len(source_labels)))
            for column in range(len(source_labels))
        )
        actual = row_morphism(generator).to_vector()
        assert tuple(int(actual(label)) for label in target_labels) == expected

    eichler_morphism = eichler_row_isometry(source, row_matrix)
    same_lattice_morphism = _isometry_from_row_matrix(source, source, row_matrix)
    assert is_witness(_column_matrix(eichler_morphism).transpose().rows(), gram, gram)
    for generator in source.module_generators():
        assert eichler_morphism(generator) == same_lattice_morphism(generator)


def test_equal_gram_groupoid_identification_preserves_coordinate_labels() -> None:
    """Equal Gram presentations with different labels are identified by coordinates."""
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ

    from sage_indefinite_port.indefinite.groupoid_cache import coordinate_identification

    source_labels = ("zeta", "alpha")
    target_labels = ("beta", "omega")
    source = Lattices(SageZZ)([[2, 1], [1, 2]], module_generators=source_labels)
    target = Lattices(SageZZ)([[2, 1], [1, 2]], module_generators=target_labels)
    assert source is not target
    assert tuple(source.module_generating_set()) == source_labels
    assert tuple(target.module_generating_set()) == target_labels

    def generator_coordinate_order(lattice, labels):
        return tuple(
            next(index for index, label in enumerate(labels) if generator.to_vector()(label))
            for generator in lattice.module_generators()
        )

    source_order = generator_coordinate_order(source, source_labels)
    target_order = generator_coordinate_order(target, target_labels)
    assert source_order != target_order
    witness = coordinate_identification(source, target)

    for index, source_label in enumerate(source_labels):
        image = witness(source.module_generator(source_label)).to_vector()
        assert tuple(int(image(label)) for label in target_labels) == tuple(int(index == j) for j in range(len(target_labels)))


def test_degenerate_isometry_witness_respects_complement_coordinates() -> None:
    """A non-diagonal nonsingular summand plus radical survives a basis change."""
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ
    from sage.all import matrix

    require(MARKED_ISOMETRY)

    gram = matrix(SageZZ, [[2, 1, 0], [1, 2, 0], [0, 0, 0]])
    change = matrix(SageZZ, [[1, 1, 0], [0, 1, 1], [0, 0, 1]])
    conjugate = change * gram * change.transpose()
    source_labels = ("zeta", "alpha", "mu")
    target_labels = ("omega", "beta", "null")
    source = Lattices(SageZZ)(gram, module_generators=source_labels)
    target = Lattices(SageZZ)(conjugate, module_generators=target_labels)
    assert tuple(source.module_generating_set()) == source_labels
    assert tuple(target.module_generating_set()) == target_labels

    assert source.is_isometric(target)
    witness = source.Isom(target).an_element()
    assert is_witness(rows_of(witness), [list(row) for row in gram.rows()], [list(row) for row in conjugate.rows()])


def test_eichler_generator_families_preserve_nonlexical_complement() -> None:
    """GRPeasy and full generator matrices preserve a model with non-lexical complement labels."""
    from dzack_research.preamble.all import ZZ as PreambleZZ
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ

    require(MARKED_ISOMETRY)
    category = Lattices(PreambleZZ)
    complement = category([[2, 0], [0, 4]], module_generators=("zeta", "alpha"))
    lattice_ = category.orthogonal_sum((category("U"), category("U"), complement))
    model = lattice_.two_u_eichler_model_from_represented_biproduct()
    selected_complement = model.orthogonal_complement()
    labels = tuple(selected_complement.module_generating_set())
    assert labels == ("zeta", "alpha")

    from sage_indefinite_port.indefinite.eichler import OrbitCoverModel

    lattice_basis = tuple(model.lattice().module_generators())
    model_basis = model.hyperbolic_basis() + model._embedded_complement_basis()
    easy_generators = OrbitCoverModel.source_easy_generators(model)
    assert easy_generators
    assert all(
        model.lattice().b(generator(left), generator(right)) == model.lattice().b(left, right)
        for generator in easy_generators
        for left in lattice_basis
        for right in lattice_basis
    )
    for selected_label in labels:
        selected_vector = model.complement_inclusion()(selected_complement.module_generator(selected_label))
        expected_images = tuple(
            -point if point == selected_vector else point
            for point in model_basis
        )
        assert any(
            all(generator(point) == image for point, image in zip(model_basis, expected_images, strict=True))
            for generator in easy_generators
        )
    from sage.all import matrix as sage_matrix

    model_gram = sage_matrix(
        SageZZ,
        [
            [model.lattice().b(left, right) for right in model_basis]
            for left in model_basis
        ],
    )
    full_matrices = OrbitCoverModel.source_generator_matrices(model)
    assert full_matrices
    assert all(
        generator_matrix * model_gram * generator_matrix.transpose() == model_gram
        for generator_matrix in full_matrices
    )


def test_odd_eichler_envelope_uses_declared_coordinate_order() -> None:
    """The even-sublattice envelope construction follows the declared basis labels."""
    from dzack_research.preamble.all import ZZ as PreambleZZ
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ

    require(MARKED_ISOMETRY)
    labels = ("z_odd", "a_first", "m_second", "b_radical", "c_extra")
    gram = [
        [1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1],
        [0, 0, 0, 1, 0],
    ]
    permutation = (3, 0, 4, 1, 2)
    permuted_gram = [[gram[i][j] for j in permutation] for i in permutation]
    permuted_labels = tuple(labels[index] for index in permutation)
    multiple_odd_labels = ("z_odd", "a_odd", "m_first", "b_first", "c_second", "d_second")
    multiple_odd_gram = [
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0],
    ]
    multiple_odd_permutation = (5, 1, 3, 0, 4, 2)
    multiple_odd_permuted_gram = [
        [multiple_odd_gram[i][j] for j in multiple_odd_permutation]
        for i in multiple_odd_permutation
    ]
    multiple_odd_permuted_labels = tuple(multiple_odd_labels[index] for index in multiple_odd_permutation)
    presentations = (
        (gram, labels),
        (permuted_gram, permuted_labels),
        (multiple_odd_gram, multiple_odd_labels),
        (multiple_odd_permuted_gram, multiple_odd_permuted_labels),
    )

    from sage_indefinite_port.indefinite.eichler import EichlerEnvelope

    for gram_matrix, coordinate_labels in presentations:
        lattice_ = Lattices(PreambleZZ)(gram_matrix, module_generators=coordinate_labels)
        assert not lattice_.is_even()
        assert tuple(lattice_.module_generating_set()) == coordinate_labels
        envelope = EichlerEnvelope(lattice_)
        assert envelope.scale == 2
        assert envelope.similarity.multiplier() == 4
        assert envelope.envelope.is_even()
        images = tuple(envelope.similarity(lattice_.module_generator(label)) for label in coordinate_labels)
        for i, left in enumerate(images):
            for j, right in enumerate(images):
                assert envelope.envelope.b(left, right) == envelope.scale**2 * lattice_.b(
                    lattice_.module_generator(coordinate_labels[i]), lattice_.module_generator(coordinate_labels[j])
                )
        envelope_labels = tuple(envelope.envelope.module_generating_set())
        assert all(
            all(image.to_vector()(label) in SageZZ for label in envelope_labels)
            for image in images
        )


def test_conway_sloane_same_genus_pair_is_not_isometric() -> None:
    require(MARKED_ISOMETRY)

    pair = load_conway_sloane_cases()["spinor_genus_pair_determinant_minus_128"]
    first, second = lattice(pair["form_a"]["gram"]), lattice(pair["form_b"]["gram"])

    assert first.genus() == second.genus()
    assert not first.is_isometric(second)
    assert first.Isom(second).is_empty()


def _unimodular_change_of_basis(rank: int) -> list[list[int]]:
    """The upper unitriangular all-ones matrix: unimodular and mixes every basis vector."""
    return [[1 if j >= i else 0 for j in range(rank)] for i in range(rank)]


@pytest.mark.parametrize("case", load_ci_indefinite_comp(), ids=lambda case: case["id"])
def test_change_of_basis_gives_an_isometric_lattice(case) -> None:
    from sage.all import ZZ, matrix

    require(MARKED_ISOMETRY)

    gram = matrix(ZZ, case["gram"])
    change = matrix(ZZ, _unimodular_change_of_basis(gram.nrows()))
    conjugate = change * gram * change.transpose()

    _assert_isometric(case["gram"], [list(row) for row in conjugate.rows()])
