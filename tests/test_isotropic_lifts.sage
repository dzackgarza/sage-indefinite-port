r"""Acceptance specimens for exact isotropic extension equations."""

from sage.all import MatrixSpace, QQ

from dzack_research.preamble.all import Lattices, ZZ as OwnedZZ

from sage_indefinite_port.indefinite.isotropic_lifts import (
    CodimensionOneIsotropicExtension,
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
