r"""Acceptance specimens for exact isotropic extension equations."""

from sage.all import MatrixSpace, QQ
import pytest

from dzack_research.preamble.all import Lattices, ZZ as OwnedZZ

from sage_indefinite_port.indefinite.isotropic_lifts import (
    CodimensionOneIsotropicExtension,
    IsometryExtensionTorsor,
    NoIntegralExtensionError,
    solve_isotropic_extension_equation,
)


def test_extension_equation_has_expected_affine_dimension() -> None:
    for rank in range(1, 5):
        matrices = MatrixSpace(QQ, rank, rank)
        A = matrices.random_element()
        while not A.is_invertible():
            A = matrices.random_element()
        expected = matrices.random_element()
        B = expected * A + A.transpose() * expected.transpose()

        solution = solve_isotropic_extension_equation(A, B)

        X = solution.particular
        assert X * A + A.transpose() * X.transpose() == B
        assert solution.homogeneous_space.dimension() == rank * (rank - 1) // 2


def test_codimension_one_extension_retains_subspace_inclusions() -> None:
    plane = Lattices(OwnedZZ)("U").base_change(OwnedZZ.fraction_field_map())
    isotropic, _partner = plane.module_generators()
    line = plane.subobject_on((isotropic,))
    partial = line.identity_morphism()

    result = CodimensionOneIsotropicExtension(line, line, partial).rational_extension()

    assert result.source_inclusion is line.inclusion()
    assert result.target_inclusion is line.inclusion()
    assert result.extension(isotropic) == isotropic
    assert result.extension in plane.Isom(plane)


def test_extension_torsor_reports_an_empty_integral_locus() -> None:
    integral_plane = Lattices(OwnedZZ)("U")
    rational_plane = integral_plane.base_change(OwnedZZ.fraction_field_map())
    field = rational_plane.base_ring()
    e, f = rational_plane.module_generators()
    rational_isometry = rational_plane.Isom(rational_plane)(
        (
            rational_plane.scalar_multiple(field(2), e),
            rational_plane.scalar_multiple(field(1) / field(2), f),
        )
    )
    torsor = IsometryExtensionTorsor(rational_isometry, ())

    assert torsor.integral_parameters(integral_plane, integral_plane) is None
    with pytest.raises(NoIntegralExtensionError):
        torsor.one_integral_extension()


def test_every_E8_orthogonal_generator_lifts_through_U_plus_E8() -> None:
    plane = Lattices(OwnedZZ)("U")
    root_lattice = Lattices(OwnedZZ)("E8")
    lattice = plane + root_lattice
    ambient_labels = tuple(lattice.module_generating_set())
    ambient_generators = tuple(lattice.module_generators())
    root_generators = tuple(root_lattice.module_generators())
    isotropic = ambient_generators[0]

    for generator in root_lattice.Aut().group_generators():
        images = list(ambient_generators[:2])
        for root_generator in root_generators:
            coordinates = generator(root_generator).to_vector()
            images.append(
                lattice.linear_combination(
                    {
                        ambient_labels[index + 2]: coordinates(root_label)
                        for index, root_label in enumerate(root_lattice.module_generating_set())
                        if coordinates(root_label)
                    }
                )
            )
        lifted = lattice.O()(tuple(images))
        assert lifted in lattice.O()
        assert lifted(isotropic) == isotropic
