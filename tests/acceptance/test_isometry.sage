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
    # This is an actual A2 isometry: R G R^T = G.  The former shear
    # [[1, 0], [1, 1]] did not preserve this Gram form.
    rows = matrix(SageZZ, [[-1, 1], [0, 1]])
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

    source_labels = ("a", "z")
    target_labels = ("z", "a")
    gram = [[2, 1], [1, 2]]
    source = Lattices(SageZZ)(gram, module_generators=source_labels)
    target = Lattices(SageZZ)(gram, module_generators=target_labels)
    assert tuple(source.module_generating_set()) == source_labels
    assert tuple(target.module_generating_set()) == target_labels

    def generator_coordinate_order(lattice_, labels):
        order = []
        for generator in lattice_.module_generators():
            coordinates = generator.to_vector()
            nonzero = tuple(index for index, label in enumerate(labels) if coordinates(label))
            assert len(nonzero) == 1
            index = nonzero[0]
            assert int(coordinates(labels[index])) == 1
            assert all(not coordinates(label) for label in labels if label != labels[index])
            order.append(index)
        assert sorted(order) == list(range(len(labels)))
        return tuple(order)

    assert generator_coordinate_order(source, source_labels) == (0, 1)
    assert generator_coordinate_order(target, target_labels) == (0, 1)
    assert source_labels != target_labels

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

    try:
        eichler_row_isometry(source, matrix(SageZZ, [[1, 0], [0, 2]]))
    except ValueError as error:
        assert "does not define an isometry" in str(error)
    else:
        raise AssertionError("Eichler row-matrix constructor accepted a non-isometry")

    try:
        _isometry_from_row_matrix(source, target, matrix(SageZZ, [[1, 0], [0, 2]]))
    except ValueError as error:
        assert "does not define an isometry" in str(error)
    else:
        raise AssertionError("row-matrix constructor accepted a non-isometry")

    try:
        _isometry_from_row_matrix(source, target, matrix(SageZZ, [[1, 0, 0]]))
    except ValueError as error:
        assert "shape" in str(error)
    else:
        raise AssertionError("row-matrix constructor accepted a wrong-shaped matrix")

    from sage.all import QQ as SageQQ

    rational_target = Lattices(SageQQ)(gram, module_generators=target_labels)
    try:
        _isometry_from_row_matrix(source, rational_target, row_matrix)
    except ValueError as error:
        assert "common base ring" in str(error)
    else:
        raise AssertionError("row-matrix constructor accepted different base rings")

    try:
        _isometry_from_column_matrix(source, target, matrix(SageZZ, [[1], [0]]))
    except ValueError as error:
        assert "shape" in str(error)
    else:
        raise AssertionError("column-matrix constructor accepted a wrong-shaped matrix")

    try:
        _isometry_from_column_matrix(source, target, matrix(SageZZ, [[1, 0], [0, 2]]))
    except ValueError as error:
        assert "does not define an isometry" in str(error)
    else:
        raise AssertionError("column-matrix constructor accepted a non-isometry")

    try:
        _isometry_from_column_matrix(source, rational_target, columns)
    except ValueError as error:
        assert "common base ring" in str(error)
    else:
        raise AssertionError("column-matrix constructor accepted different base rings")


def test_equal_gram_groupoid_identification_preserves_coordinate_labels() -> None:
    """Equal Gram presentations with different labels are identified by coordinates."""
    from dzack_research.preamble.all import Lattices
    from sage.all import ZZ as SageZZ

    from sage_indefinite_port.indefinite.groupoid_cache import (
        IsometryGroupoidCache,
        coordinate_identification,
    )

    source_labels = ("zeta", "alpha")
    target_labels = ("beta", "omega")
    source = Lattices(SageZZ)([[2, 1], [1, 2]], module_generators=source_labels)
    target = Lattices(SageZZ)([[2, 1], [1, 2]], module_generators=target_labels)
    assert source is not target
    assert tuple(source.module_generating_set()) == source_labels
    assert tuple(target.module_generating_set()) == target_labels

    def generator_coordinate_order(lattice, labels):
        order = []
        for generator in lattice.module_generators():
            coordinates = generator.to_vector()
            nonzero = tuple(index for index, label in enumerate(labels) if coordinates(label))
            assert len(nonzero) == 1
            index = nonzero[0]
            assert int(coordinates(labels[index])) == 1
            assert all(not coordinates(label) for label in labels if label != labels[index])
            order.append(index)
        assert sorted(order) == list(range(len(labels)))
        return tuple(order)

    source_order = generator_coordinate_order(source, source_labels)
    target_order = generator_coordinate_order(target, target_labels)
    assert source_order == (0, 1)
    assert target_order == (0, 1)
    assert source_labels != target_labels
    witness = coordinate_identification(source, target)

    for index, source_label in enumerate(source_labels):
        image = witness(source.module_generator(source_label)).to_vector()
        assert tuple(int(image(label)) for label in target_labels) == tuple(int(index == j) for j in range(len(target_labels)))

    different_gram = Lattices(SageZZ)([[2, 0], [0, 4]], module_generators=target_labels)
    try:
        coordinate_identification(source, different_gram)
    except ValueError as error:
        assert "identical Gram matrices" in str(error)
    else:
        raise AssertionError("coordinate identification accepted unequal Gram presentations")

    different_rank = Lattices(SageZZ)([[2]], module_generators=("only",))
    try:
        coordinate_identification(source, different_rank)
    except ValueError as error:
        assert "equal coordinate ranks" in str(error)
    else:
        raise AssertionError("coordinate identification accepted unequal coordinate ranks")

    from sage.all import QQ as SageQQ

    different_ring = Lattices(SageQQ)([[2, 1], [1, 2]], module_generators=target_labels)
    try:
        coordinate_identification(source, different_ring)
    except ValueError as error:
        assert "common base ring" in str(error)
    else:
        raise AssertionError("coordinate identification accepted different base rings")

    cache = IsometryGroupoidCache()
    try:
        cache._key(different_ring)
    except TypeError as error:
        assert "integral ZZ-lattices" in str(error)
    else:
        raise AssertionError("the integer-keyed groupoid cache accepted a rational lattice")
    assert cache._representatives == {}, "rejected rational lattices must not mutate the integer-keyed cache"


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
    source = Lattices(SageZZ)([[int(entry) for entry in row] for row in gram.rows()], module_generators=source_labels)
    target = Lattices(SageZZ)([[int(entry) for entry in row] for row in conjugate.rows()], module_generators=target_labels)
    assert tuple(source.module_generating_set()) == source_labels
    assert tuple(target.module_generating_set()) == target_labels

    def generator_coordinate_order(lattice_, labels):
        order = []
        for generator in lattice_.module_generators():
            coordinates = generator.to_vector()
            nonzero = tuple(index for index, label in enumerate(labels) if coordinates(label))
            assert len(nonzero) == 1
            index = nonzero[0]
            assert int(coordinates(labels[index])) == 1
            assert all(not coordinates(label) for label in labels if label != labels[index])
            order.append(index)
        assert sorted(order) == list(range(len(labels)))
        return tuple(order)

    assert generator_coordinate_order(source, source_labels) == (0, 1, 2)
    assert generator_coordinate_order(target, target_labels) == (0, 1, 2)
    assert source_labels != target_labels
    assert source.is_isometric(target)
    from sage_indefinite_port.indefinite.recursive import _column_matrix

    witness = source.Isom(target).an_element()
    witness_rows = _column_matrix(witness).transpose().rows()
    assert is_witness(witness_rows, [list(row) for row in gram.rows()], [list(row) for row in conjugate.rows()])


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


def test_two_u_full_generating_data_accepts_raw_sage_integer_norm() -> None:
    """The public full-generation consumer returns its norm in the owned ring."""
    from dzack_research.preamble.all import Lattices, ZZ as PreambleZZ
    from sage.all import QQ as SageQQ
    from sage.all import ZZ as SageZZ

    category = Lattices(PreambleZZ)
    lattice_ = category.orthogonal_sum((category("U"), category("U"), category([[2]])))
    model = lattice_.two_u_eichler_model_from_represented_biproduct()

    try:
        model.covering_vector_representatives(SageQQ(3) / 2)
    except ValueError as error:
        assert "must be integral" in str(error)
    else:
        raise AssertionError("the covering-vector consumer truncated a nonintegral rational square")

    datum = model.full_generating_data(SageZZ(2))

    two = PreambleZZ.one() + PreambleZZ.one()
    assert datum.square() == two
    assert datum.base_vector().q() == two
    assert tuple(datum.generating_family())


def test_odd_eichler_envelope_uses_declared_coordinate_order() -> None:
    """The even-sublattice envelope construction follows the declared basis labels."""
    from dzack_research.preamble.all import ZZ as PreambleZZ
    from dzack_research.preamble.all import Lattices

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
        assert int(envelope.similarity.multiplier()) == 4
        assert envelope.envelope.is_even()
        if envelope._even_inclusion is not None:
            # The norm-parity kernel has index two, so doubling every
            # original basis vector must land in the even sublattice.
            even_inclusion = envelope._even_inclusion
            assert even_inclusion.domain().is_even()
            assert int(even_inclusion.index()) == 2
            fraction_map = lattice_.base_ring().fraction_field_map()
            rational_inclusion = even_inclusion.base_change(fraction_map)
            rational_even = rational_inclusion.domain()
            rational_original = rational_inclusion.codomain()
            rational_inverse = rational_inclusion.inverse()
            even_labels = tuple(rational_even.module_generating_set())
            for coordinate_label in coordinate_labels:
                preimage = rational_inverse(rational_original.module_generator(coordinate_label))
                two = lattice_.base_ring().one() + lattice_.base_ring().one()
                doubled = rational_even.scalar_multiple(two, preimage)
                assert all(doubled.to_vector()(label).denominator() == 1 for label in even_labels)
        images = tuple(envelope.similarity(lattice_.module_generator(label)) for label in coordinate_labels)
        for i, left in enumerate(images):
            for j, right in enumerate(images):
                assert envelope.envelope.b(left, right) == envelope.scale**2 * lattice_.b(
                    lattice_.module_generator(coordinate_labels[i]), lattice_.module_generator(coordinate_labels[j])
                )
        envelope_labels = tuple(envelope.envelope.module_generating_set())
        coordinate_denominators = tuple(
            tuple(int(image.to_vector()(label).denominator()) for label in envelope_labels)
            for image in images
        )
        coordinate_values = tuple(
            tuple(str(image.to_vector()(label)) for label in envelope_labels)
            for image in images
        )
        assert all(all(denominator == 1 for denominator in row) for row in coordinate_denominators), (
            coordinate_labels,
            coordinate_denominators,
            coordinate_values,
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
