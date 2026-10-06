"""Integral-structure actions for rational isometry groups."""

from collections.abc import Callable

import pytest

from dzack_research.preamble.all import Lattices, Modules, QQ, RestrictedScalarsModules, ZZ
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.module_morphisms.module_morphisms import (
    ModuleEmbeddingMethods,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.graphs.graph import Graph
from sage.groups.perm_gps.permgroup import PermutationGroup_generic
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.libs.gap.element import GapElement, GapElement_Permutation
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import matrix
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ

from sage_indefinite_port.groups.integral_structures import (
    ArithmeticSubgroup,
    CentralizerSubgroup,
    DoubleCosetDecomposition,
    FinitePermutationRepresentation,
    FinitePreimageSubgroup,
    GeneratedSubgroup,
    IntegralStructureAction,
    IntersectionSubgroup,
    KernelSubgroup,
    RationalMatrixGroup,
    StabilizerSubgroup,
)
from tests.fixtures.oracle_fixtures import (
    DoubleCosetCase,
    PermutationImages,
    PermutationGroupData,
    RootSystemCase,
    load_classification_simplices,
    load_double_coset_cases,
    load_root_systems,
)


def _standard_hyperbolic_lattice() -> tuple[
    Lattices.ParentMethods,
    RestrictedScalarsModules.ParentMethods,
    ModuleEmbeddingMethods,
    Lattices.ElementMethods,
    Lattices.ElementMethods,
]:
    plane = Lattices(QQ)("U")
    restriction = Modules(QQ).restriction_of_scalars(
        ZZ.Mor(QQ)(lambda element: QQ(element))
    )
    space = restriction(plane)
    e, f = plane.module_generators()
    standard_domain = ZZ.free_module(2)
    standard = standard_domain.Mono(space)(
        {0: space.wrap(e), 1: space.wrap(f)},
    )
    return plane, space, standard, e, f


def _fixture_gap_group(
    data: PermutationGroupData,
) -> tuple[GapElement, tuple[GapElement_Permutation, ...]]:
    generators = tuple(
        libgap.PermList([image + 1 for image in permutation])
        for permutation in data["generators"]
    )
    return libgap.Group(list(generators)), generators


def _permutation_images(
    permutation: GapElement,
    degree: int,
) -> PermutationImages:
    return tuple(
        int(image) - 1
        for image in libgap.ListPerm(permutation, degree).sage()
    )


def _basis_permutation(
    action: IntegralStructureAction,
    automorphism: LatticeIsometryMethods,
) -> GapElement_Permutation:
    transformation = action._ambient_action_matrix(automorphism)
    images = []
    for column in range(transformation.ncols()):
        matching = tuple(
            row + 1
            for row in range(transformation.nrows())
            if transformation[row, column] == 1
        )
        assert len(matching) == 1
        assert all(
            transformation[row, column] in (0, 1)
            for row in range(transformation.nrows())
        )
        images.append(matching[0])
    return libgap.PermList(images)


def _double_coset_fixture_action(
    case: DoubleCosetCase,
) -> tuple[
    IntegralStructureAction,
    GapElement,
    GapElement,
    GapElement,
    Callable[[GapElement], LatticeIsometryMethods],
]:
    big, _big_generators = _fixture_gap_group(case["group_h"])
    small, _small_generators = _fixture_gap_group(case["group_g"])
    assert int(big.Size()) == case["big_group_order"]
    assert int(small.Size()) == case["small_group_order"]
    assert bool(libgap.IsSubgroup(big, small))

    right_cosets = libgap.RightCosets(big, small)
    index = case["big_group_order"] // case["small_group_order"]
    assert int(libgap.Length(right_cosets)) == index
    coset_action = libgap.ActionHomomorphism(big, right_cosets, libgap.OnRight)
    image_group = libgap.Image(coset_action)
    image_generators = tuple(libgap.SmallGeneratingSet(image_group))

    rational_lattice = Lattices(QQ)(index)
    labels = tuple(rational_lattice.module_generating_set())
    basis = tuple(rational_lattice.module_generator(label) for label in labels)
    automorphisms = rational_lattice.Aut()
    live_isometries: dict[PermutationImages, LatticeIsometryMethods] = {}

    def live_isometry(permutation: GapElement) -> LatticeIsometryMethods:
        coset_images = _permutation_images(permutation, index)
        match coset_images in live_isometries:
            case True:
                return live_isometries[coset_images]
            case False:
                pass
        transformation = matrix(
            SageQQ,
            [
                [
                    1 if row == coset_images[column] else 0
                    for column in range(index)
                ]
                for row in range(index)
            ],
        )
        result = automorphisms._isometry_from_column_matrix(transformation)
        live_isometries[coset_images] = result
        return result

    group = RationalMatrixGroup(
        rational_lattice,
        tuple(live_isometry(generator) for generator in image_generators),
    )
    restriction = Modules(QQ).restriction_of_scalars(
        ZZ.Mor(QQ)(lambda element: QQ(element))
    )
    space = restriction(rational_lattice)
    integral_module = ZZ.free_module(index)
    integral_labels = tuple(integral_module.module_generating_set())
    selected = integral_module.Mono(space)(
        {
            integral_labels[position]: space.wrap(
                rational_lattice.scalar_multiple(
                    QQ(ZZ.one() + ZZ.one()) if position == 0 else QQ.one(),
                    basis[position],
                )
            )
            for position in range(index)
        },
    )
    return (
        IntegralStructureAction(group, selected),
        big,
        small,
        coset_action,
        live_isometry,
    )


def _configuration_group(
    rows: list[list[int]],
) -> tuple[Matrix_rational_dense, Matrix_rational_dense, PermutationGroup_generic]:
    configuration = matrix(SageQQ, rows)
    gram = configuration.transpose() * configuration
    projector = configuration * gram.inverse() * configuration.transpose()
    points = range(configuration.nrows())
    graph = Graph()
    graph.add_vertices(points)
    for left in points:
        for right in range(left + 1, configuration.nrows()):
            graph.add_edge(left, right, projector[left, right])
    diagonal_parts: dict[Rational, list[int]] = {}
    for point in points:
        diagonal_parts.setdefault(projector[point, point], []).append(point)
    group = graph.automorphism_group(
        partition=list(diagonal_parts.values()),
        edge_labels=True,
        algorithm="sage",
    )
    return configuration, gram, group


def _configuration_matrix(
    configuration: Matrix_rational_dense,
    gram_inverse: Matrix_rational_dense,
    permutation: PermutationGroupElement,
) -> Matrix_rational_dense:
    target = matrix(
        SageQQ,
        [configuration.row(permutation(point)) for point in range(configuration.nrows())],
    )
    return target.transpose() * configuration * gram_inverse


def _configuration_action(
    rows: list[list[int]],
) -> tuple[
    IntegralStructureAction,
    Matrix_rational_dense,
    Matrix_rational_dense,
    PermutationGroup_generic,
    Callable[[PermutationGroupElement], LatticeIsometryMethods],
]:
    configuration, gram, permutation_group = _configuration_group(rows)
    rational_gram = gram.inverse()
    rank = configuration.ncols()
    rational_lattice = Lattices(QQ)(
        [
            [
                QQ(ZZ(int(rational_gram[row, column].numerator())))
                / QQ(ZZ(int(rational_gram[row, column].denominator())))
                for column in range(rank)
            ]
            for row in range(rank)
        ]
    )
    labels = tuple(rational_lattice.module_generating_set())

    def live_isometry(
        permutation: PermutationGroupElement,
    ) -> LatticeIsometryMethods:
        transformation = _configuration_matrix(configuration, rational_gram, permutation)
        assert transformation.transpose() * rational_gram * transformation == rational_gram
        return rational_lattice.Aut()._isometry_from_column_matrix(transformation)

    group = RationalMatrixGroup(
        rational_lattice,
        tuple(live_isometry(generator) for generator in permutation_group.gens()),
    )
    restriction = Modules(QQ).restriction_of_scalars(
        ZZ.Mor(QQ)(lambda element: QQ(element))
    )
    space = restriction(rational_lattice)
    integral_module = ZZ.free_module(rank)
    integral_labels = tuple(integral_module.module_generating_set())
    basis = tuple(rational_lattice.module_generator(label) for label in labels)
    selected = integral_module.Mono(space)(
        {
            integral_labels[position]: space.wrap(basis[position])
            for position in range(rank)
        },
    )
    return (
        IntegralStructureAction(group, selected),
        configuration,
        rational_gram,
        permutation_group,
        live_isometry,
    )


def _configuration_permutation(
    action: IntegralStructureAction,
    automorphism: LatticeIsometryMethods,
    configuration: Matrix_rational_dense,
) -> GapElement_Permutation:
    transformation = action._ambient_action_matrix(automorphism)
    moved = configuration * transformation.transpose()
    images = []
    for row in moved.rows():
        matching = tuple(
            position
            for position, candidate in enumerate(configuration.rows())
            if row == candidate
        )
        assert len(matching) == 1
        images.append(matching[0] + 1)
    return libgap.PermList(images)


def _integral_configuration_subgroup(
    configuration: Matrix_rational_dense,
    gram_inverse: Matrix_rational_dense,
    permutation_group: PermutationGroup_generic,
) -> GapElement:
    integral_permutations = []
    for permutation in permutation_group:
        transformation = _configuration_matrix(configuration, gram_inverse, permutation)
        if all(entry.denominator() == 1 for entry in transformation.list()):
            integral_permutations.append(
                libgap.PermList(
                    [
                        int(permutation(point)) + 1
                        for point in range(configuration.nrows())
                    ]
                )
            )
    return libgap.Group(integral_permutations)


def _assert_configuration_stabilizer(
    rows: list[list[int]],
) -> tuple[IntegralStructureAction, Matrix_rational_dense, GapElement, int]:
    action, configuration, gram_inverse, permutation_group, _live_isometry = _configuration_action(rows)
    expected = _integral_configuration_subgroup(
        configuration,
        gram_inverse,
        permutation_group,
    )

    expected_order = int(expected.Size())
    rational_order = int(permutation_group.order())
    assert rational_order % expected_order == 0

    stabilizer = action.lattice_stabilizer()
    returned = libgap.Group(
        [
            _configuration_permutation(action, generator, configuration)
            for generator in stabilizer.generators()
        ]
    )
    assert int(returned.Size()) == expected_order
    assert bool(libgap.IsSubgroup(expected, returned))
    return action, configuration, expected, rational_order


def _assert_configuration_double_cosets(
    action: IntegralStructureAction,
    configuration: Matrix_rational_dense,
    expected: GapElement,
    rational_order: int,
) -> None:
    decomposition = action.double_cosets(action.lattice_stabilizer())
    source_representatives: list[GapElement_Permutation] = []
    source_double_coset_sizes: list[int] = []
    for representative in decomposition.representatives():
        source = _configuration_permutation(action, representative, configuration)
        source_double_coset = libgap.DoubleCoset(expected, source, expected)
        assert all(
            previous not in source_double_coset
            for previous in source_representatives
        )
        source_representatives.append(source)
        source_double_coset_sizes.append(int(source_double_coset.Size()))
    assert sum(source_double_coset_sizes) == rational_order
    assert sum(
        piece.double_coset_size() for piece in decomposition.intersections()
    ) == decomposition.finite_ambient_order()


def test_an_integral_group_keeps_the_selected_lattice() -> None:
    plane, _space, standard, e, f = _standard_hyperbolic_lattice()
    swap = plane.Aut()({0: f, 1: e})
    action = IntegralStructureAction(RationalMatrixGroup(plane, (swap,)), standard)

    assert action.invariant_overlattice() is standard
    assert int(action.quotient_exponent()) == 1


def test_a_rational_involution_closes_to_the_smallest_invariant_overlattice() -> None:
    plane, space, standard, e, f = _standard_hyperbolic_lattice()
    involution = plane.Aut()(
        {
            0: plane.scalar_multiple(QQ(ZZ.one() + ZZ.one()), f),
            1: plane.scalar_multiple(QQ.one() / QQ(ZZ.one() + ZZ.one()), e),
        }
    )
    assert involution * involution == plane.Aut().identity()

    action = IntegralStructureAction(
        RationalMatrixGroup(plane, (involution,)),
        standard,
    )
    invariant = action.invariant_overlattice()
    half_e = space.wrap(
        plane.scalar_multiple(QQ.one() / QQ(ZZ.one() + ZZ.one()), e)
    )

    assert invariant.codomain() is space
    assert invariant.is_in_image(space.wrap(e))
    assert invariant.is_in_image(space.wrap(f))
    assert invariant.is_in_image(half_e)
    assert int(action.quotient_exponent()) == 2

    for basis_vector in invariant.domain().module_generators():
        embedded = invariant(basis_vector).underlying_element()
        moved = space.wrap(involution(embedded))
        assert invariant.is_in_image(moved)


def test_the_finite_module_orbit_recovers_the_lattice_stabilizer() -> None:
    plane, _space, standard, e, f = _standard_hyperbolic_lattice()
    involution = plane.Aut()(
        {
            0: plane.scalar_multiple(QQ(ZZ.one() + ZZ.one()), f),
            1: plane.scalar_multiple(QQ.one() / QQ(ZZ.one() + ZZ.one()), e),
        }
    )
    group = RationalMatrixGroup(plane, (involution,))
    action = IntegralStructureAction(group, standard)
    finite = action.finite_representation()

    assert int(action.finite_module().cardinality()) == 4
    assert len(finite.orbit()) == 2
    assert finite.image_order() == 2
    assert len(finite.orbit_witnesses()) == 2
    assert finite.orbit_witnesses()[0] == plane.Aut().identity()
    assert finite.orbit_witnesses()[1](e) in (involution(e), (~involution)(e))
    assert finite.orbit_images(e) == tuple(witness(e) for witness in finite.orbit_witnesses())
    assert all(
        finite._submodule_key(submodule) == key
        for submodule, key in zip(
            finite.orbit(),
            finite._orbit_keys,
            strict=True,
        )
    )

    stabilizer = action.lattice_stabilizer()
    assert stabilizer.supergroup() is group
    assert all(action.preserves_selected_lattice(generator) for generator in stabilizer.generators())
    assert not action.preserves_selected_lattice(involution)


def test_matrix_backed_rational_group_preserves_the_finite_action() -> None:
    plane, _space, standard, e, f = _standard_hyperbolic_lattice()
    involution = plane.Aut()(
        {
            0: plane.scalar_multiple(QQ(ZZ.one() + ZZ.one()), f),
            1: plane.scalar_multiple(QQ.one() / QQ(ZZ.one() + ZZ.one()), e),
        }
    )
    matrix_ = plane.Aut()._row_action_matrix(involution).transpose()
    group = RationalMatrixGroup._from_column_matrices(plane, (matrix_,))
    action = IntegralStructureAction(group, standard)

    assert group._generator_column_matrices() == (matrix_,)
    assert group.generators() == (involution,)
    finite = action.finite_representation()
    assert finite.image_order() == 2
    matrix_keys = {
        tuple(tuple(entry for entry in row) for row in matrix_.rows())
        for matrix_ in finite._lattice_stabilizer_column_matrices
    }
    live_keys = {
        tuple(tuple(entry for entry in row) for row in plane.Aut()._row_action_matrix(generator).transpose().rows())
        for generator in finite.lattice_stabilizer().generators()
    }
    assert matrix_keys == live_keys


def test_transporter_and_right_subgroup_cosets_retain_their_actual_sides() -> None:
    plane, space, standard, e, f = _standard_hyperbolic_lattice()
    involution = plane.Aut()(
        {
            0: plane.scalar_multiple(QQ(ZZ.one() + ZZ.one()), f),
            1: plane.scalar_multiple(QQ.one() / QQ(ZZ.one() + ZZ.one()), e),
        }
    )
    group = RationalMatrixGroup(plane, (involution,))
    action = IntegralStructureAction(group, standard)
    target_domain = ZZ.free_module(2)
    target = target_domain.Mono(space)(
        {
            0: space.wrap(involution(e)),
            1: space.wrap(involution(f)),
        },
    )

    witness = action.transporter(standard, target)
    assert witness is not None
    for label in standard.domain().module_generating_set():
        vector = standard(standard.domain().module_generator(label)).underlying_element()
        assert target.is_in_image(space.wrap(witness(vector)))

    reverse = action.transporter(target, standard)
    assert reverse is not None
    for label in target.domain().module_generating_set():
        vector = target(target.domain().module_generator(label)).underlying_element()
        assert standard.is_in_image(space.wrap(reverse(vector)))

    cosets = action.right_cosets()
    assert cosets.ambient_group() is group
    assert cosets.right_subgroup() is action.lattice_stabilizer()
    assert cosets.cardinality() == 2
    assert len(cosets.representatives()) == 2
    assert cosets.representatives()[0] == plane.Aut().identity()


def test_double_cosets_name_both_subgroups_and_exhaust_the_finite_image() -> None:
    plane, _space, standard, e, f = _standard_hyperbolic_lattice()
    involution = plane.Aut()(
        {
            0: plane.scalar_multiple(QQ(ZZ.one() + ZZ.one()), f),
            1: plane.scalar_multiple(QQ.one() / QQ(ZZ.one() + ZZ.one()), e),
        }
    )
    group = RationalMatrixGroup(plane, (involution,))
    action = IntegralStructureAction(group, standard)
    trivial = ArithmeticSubgroup(group, (plane.Aut().identity(),))

    decomposition = action.double_cosets(trivial)
    assert decomposition.left_subgroup() is trivial
    assert decomposition.ambient_group() is group
    assert decomposition.right_subgroup() is action.lattice_stabilizer()
    assert decomposition.cardinality() == 2
    assert sum(piece.double_coset_size() for piece in decomposition.intersections()) == decomposition.finite_ambient_order()
    assert all(piece.finite_image_order() == 1 for piece in decomposition.intersections())

    whole = ArithmeticSubgroup(group, (involution,))
    one_double_coset = action.double_cosets(whole)
    assert one_double_coset.cardinality() == 1
    assert one_double_coset.intersections()[0].double_coset_size() == one_double_coset.finite_ambient_order()


@pytest.mark.parametrize(
    "case",
    load_double_coset_cases(),
    ids=lambda case: case["id"],
)
def test_upstream_double_coset_cases_run_through_integral_structure_action(
    case: DoubleCosetCase,
) -> None:
    """Reproduce every ``GRP_DoubleCoset`` DBL source check through the public T2 API.

    The auxiliary coset block realizes ``G/H`` as integral structures: one
    coordinate is doubled, so its lattice stabilizer is exactly the stabilizer
    of the distinguished coset.  The original source group independently
    supplies the expected double-coset count for every face stabilizer.
    """
    action, big, small, coset_action, live_isometry = _double_coset_fixture_action(case)
    expected_index = case["big_group_order"] // case["small_group_order"]

    assert int(action.quotient_exponent()) == 2
    assert action.right_cosets().cardinality() == expected_index
    image_group = libgap.Image(coset_action)
    right_image = libgap.Stabilizer(image_group, 1)
    small_image = libgap.Image(coset_action, small)
    assert int(small_image.Size()) == int(right_image.Size())
    assert bool(libgap.IsSubgroup(right_image, small_image))
    for generator in libgap.SmallGeneratingSet(small_image):
        assert action.preserves_selected_lattice(live_isometry(generator))

    decompositions: dict[
        frozenset[PermutationImages],
        tuple[ArithmeticSubgroup, DoubleCosetDecomposition, int],
    ] = {}
    for vector in case["vectors"]:
        face = [position + 1 for position, entry in enumerate(vector) if entry]
        left_gap = libgap.Stabilizer(big, face, libgap.OnSets)
        left_image = libgap.Image(coset_action, left_gap)
        left_key = frozenset(
            _permutation_images(element, expected_index)
            for element in libgap.AsList(left_image)
        )
        match left_key in decompositions:
            case True:
                left, decomposition, source_cardinality = decompositions[left_key]
            case False:
                left = ArithmeticSubgroup(
                    action.rational_group(),
                    tuple(
                        live_isometry(generator)
                        for generator in libgap.SmallGeneratingSet(left_image)
                    ),
                )
                decomposition = action.double_cosets(left)
                # ``big/small`` is exactly the set on which ``coset_action``
                # acts, and its kernel is ``core_big(small) <= small``.
                # Therefore the orbit set ``left_gap \ big / small`` depends
                # only on ``left_image``.  Compute the independent source
                # oracle once for each exact finite left image rather than
                # repeating the same GAP double-coset computation for every
                # face with that image.
                source_cardinality = int(
                    libgap.Length(libgap.DoubleCosets(big, left_gap, small))
                )
                decompositions[left_key] = (
                    left,
                    decomposition,
                    source_cardinality,
                )

                representative_permutations = [
                    _basis_permutation(action, representative)
                    for representative in decomposition.representatives()
                ]

                double_coset_sizes = []
                for position, representative in enumerate(representative_permutations):
                    source_double_coset = libgap.DoubleCoset(
                        left_image,
                        representative,
                        right_image,
                    )
                    assert all(
                        previous not in source_double_coset
                        for previous in representative_permutations[:position]
                    )
                    double_coset_sizes.append(int(source_double_coset.Size()))
                assert sum(double_coset_sizes) == int(image_group.Size())

        assert decomposition.left_subgroup() is left
        assert decomposition.ambient_group() is action.rational_group()
        assert decomposition.right_subgroup().supergroup() is action.rational_group()
        returned_right_image = libgap.Group(
            [
                _basis_permutation(action, generator)
                for generator in decomposition.right_subgroup().generators()
            ]
        )
        assert int(returned_right_image.Size()) == int(right_image.Size())
        assert bool(libgap.IsSubgroup(right_image, returned_right_image))
        assert sum(
            piece.double_coset_size() for piece in decomposition.intersections()
        ) == decomposition.finite_ambient_order()

        assert decomposition.cardinality() == source_cardinality


_CLASSIFIED_SIMPLEX_FIXTURE = load_classification_simplices()
_CLASSIFIED_SIMPLEX_CASES = [
    (f"{family}_{position}", simplex)
    for family, data in (
        ("dim5", _CLASSIFIED_SIMPLEX_FIXTURE["dim5"]),
        ("dim6", _CLASSIFIED_SIMPLEX_FIXTURE["dim6"]),
        ("dim7", _CLASSIFIED_SIMPLEX_FIXTURE["dim7"]),
    )
    for position, simplex in enumerate(data["simplices"])
]


@pytest.mark.parametrize(
    "_case_id,rows",
    _CLASSIFIED_SIMPLEX_CASES,
    ids=[case_id for case_id, _rows in _CLASSIFIED_SIMPLEX_CASES],
)
def test_rat_int_automorphy_classified_simplices_run_through_public_action(
    _case_id: str,
    rows: list[list[int]],
) -> None:
    """Run every classified-simplex example used by ``01_RatIntAutomorphy``."""
    action, configuration, expected, rational_order = _assert_configuration_stabilizer(rows)
    _assert_configuration_double_cosets(
        action,
        configuration,
        expected,
        rational_order,
    )


@pytest.mark.parametrize(
    "case",
    load_root_systems(),
    ids=lambda case: case["id"],
)
def test_rat_int_automorphy_root_systems_run_through_public_action(
    case: RootSystemCase,
) -> None:
    """Run every simple-root-system integral-stabilizer example from ``01_RatIntAutomorphy``."""
    _assert_configuration_stabilizer(case["roots"])



def _root_reflection_group(name: str) -> tuple[Lattices.ParentMethods, RationalMatrixGroup]:
    lattice = Lattices(ZZ)(name)
    basis = tuple(lattice.module_generators())
    gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageQQ)
    automorphisms = lattice.Aut()
    generators = []
    for position in range(len(basis)):
        norm = gram[position, position]
        transformation = matrix(
            SageQQ,
            len(basis),
            len(basis),
            lambda row, column: 1 if row == column else 0,
        )
        for column in range(len(basis)):
            transformation[position, column] -= 2 * gram[column, position] / norm
        generators.append(automorphisms._isometry_from_column_matrix(transformation))
    generators = tuple(generators)
    return lattice, RationalMatrixGroup(lattice, generators)


def test_construction_aware_subgroup_carriers_retain_defining_data() -> None:
    lattice, group = _root_reflection_group("A2")
    roots = tuple(lattice.roots())
    representation = FinitePermutationRepresentation(group, roots, lambda isometry, root: isometry(root))
    assert representation.is_faithful()

    generated = GeneratedSubgroup(group, group.generators()[:1])
    preimage = FinitePreimageSubgroup(representation, group.generators()[:1])
    stabilizer = StabilizerSubgroup(representation, 0)
    centralizer = CentralizerSubgroup(representation, group.generators()[0])
    intersection = IntersectionSubgroup(representation, (stabilizer, centralizer))

    assert generated.selected_generators() == group.generators()[:1]
    assert preimage.finite_image_generators() == group.generators()[:1]
    assert stabilizer.point() == roots[0]
    assert centralizer.centralizing_element() == group.generators()[0]
    assert intersection.subgroups() == (stabilizer, centralizer)
    assert all(generator(root) == root for generator in stabilizer.generators() for root in (roots[0],))
    assert all(
        generator * group.generators()[0] == group.generators()[0] * generator
        for generator in centralizer.generators()
    )
    assert all(
        generator(roots[0]) == roots[0]
        and generator * group.generators()[0] == group.generators()[0] * generator
        for generator in intersection.generators()
    )

    determinant = FinitePermutationRepresentation(
        group,
        (0, 1),
        lambda isometry, sign: sign if isometry.determinant() == ZZ.one() else 1 - sign,
    )
    kernel = KernelSubgroup(determinant)
    assert kernel.index() == 2
    assert all(generator.determinant() == ZZ.one() for generator in kernel.generators())


def test_E8_root_action_word_lift_kernel_and_stabilizer_indices() -> None:
    lattice = Lattices(ZZ)("E8")
    group = RationalMatrixGroup(
        lattice,
        tuple(lattice.Aut().framing().group_generators()),
    )
    roots = tuple(lattice.roots())
    assert len(roots) == 240

    labels = tuple(lattice.module_generating_set())
    root_coordinates = {
        root: tuple(int(root.to_vector()(label)) for label in labels)
        for root in roots
    }
    root_by_coordinates = {
        coordinates: root
        for root, coordinates in root_coordinates.items()
    }
    automorphisms = lattice.Aut()
    action_matrices = {
        id(generator): automorphisms._row_action_matrix(generator)
        for generator in group.generators()
    }

    root_action = FinitePermutationRepresentation(
        group,
        roots,
        lambda isometry, root: root_by_coordinates[
            tuple(
                (
                    matrix(SageQQ, [root_coordinates[root]])
                    * action_matrices[id(isometry)]
                ).row(0)
            )
        ],
    )
    assert root_action.image_order() == 696729600
    word = ((0, 1), (1, -1), (2, 1))
    lifted = root_action.generator_word_lift(word)
    expected = group.generators()[2] * (~group.generators()[1]) * group.generators()[0]
    assert lifted == expected

    determinant = FinitePermutationRepresentation(
        group,
        (0, 1),
        lambda isometry, sign: sign if isometry.determinant() == ZZ.one() else 1 - sign,
    )
    determinant_kernel = KernelSubgroup(determinant)
    assert determinant_kernel.index() == 2
    assert all(generator.determinant() == ZZ.one() for generator in determinant_kernel.generators())

    root_stabilizer = StabilizerSubgroup(root_action, 0)
    assert root_stabilizer.index() == 240
    assert all(generator(roots[0]) == roots[0] for generator in root_stabilizer.generators())
