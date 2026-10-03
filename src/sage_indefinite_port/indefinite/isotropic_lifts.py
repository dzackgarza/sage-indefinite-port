r"""Exact extension problems along isotropic subspaces."""

from dataclasses import dataclass

from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeEmbedding,
    LatticeIsometry,
)
from dzack_research.preamble.categories.rings.ring_foundation import (
    _engine_element,
    _owned_engine_element,
)
from sage.matrix.constructor import matrix
from sage.matrix.matrix0 import Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module_element import vector
from sage.modules.module import Module
from sage.rings.rational_field import QQ


@dataclass(frozen=True)
class MatrixEquationSolution:
    r"""Affine solutions of \(XA+A^T X^T=B\) over \(\mathbf Q\)."""

    particular: Matrix
    homogeneous_space: Module
    homogeneous_lattice: Module | None


@dataclass(frozen=True)
class CodimensionOneIsotropicExtensionResult:
    r"""A rational extension and the two retained subspace inclusions."""

    extension: LatticeIsometry
    source_inclusion: LatticeEmbedding
    target_inclusion: LatticeEmbedding


class CodimensionOneIsotropicExtension:
    r"""Extend an isometry of codimension-one isotropic subspaces over \(\mathbf Q\)."""

    def __init__(self, source_subspace, target_subspace, partial_isometry) -> None:
        self._source_subspace = source_subspace
        self._target_subspace = target_subspace
        self._partial_isometry = partial_isometry

    def rational_extension(self) -> CodimensionOneIsotropicExtensionResult:
        r"""Return the exact norm-corrected rational extension."""
        source = self._source_subspace
        target = self._target_subspace
        partial = self._partial_isometry
        source_inclusion = source.inclusion()
        target_inclusion = target.inclusion()
        source_ambient = source_inclusion.codomain()
        target_ambient = target_inclusion.codomain()
        if int(source.module_rank()) + 1 != int(source_ambient.module_rank()):
            raise ValueError("the source subspace must have codimension one")
        if int(target.module_rank()) + 1 != int(target_ambient.module_rank()):
            raise ValueError("the target subspace must have codimension one")
        if partial.domain() is not source or partial.codomain() is not target:
            raise ValueError("the partial isometry must map the source subspace to the target subspace")

        source_complement = next(generator for generator in source_ambient.module_generators() if generator not in source)
        source_generators = tuple(source.module_generators())
        target_ambient_generators = tuple(target_ambient.module_generators())
        target_ring = target_ambient.base_ring()
        target_images = tuple(target_inclusion(partial(generator)) for generator in source_generators)
        pairing_matrix = matrix(
            QQ,
            (tuple(_engine_element(target_ring, target_ambient.b(generator, image)) for generator in target_ambient_generators) for image in target_images),
        )
        prescribed_pairings = vector(
            QQ,
            tuple(
                _engine_element(
                    source_ambient.base_ring(),
                    source_ambient.b(source_complement, source_inclusion(generator)),
                )
                for generator in source_generators
            ),
        )
        target_coordinates = pairing_matrix.solve_right(prescribed_pairings)
        target_complement = target_ambient.linear_combination(
            {
                label: _owned_engine_element(target_ring, target_coordinates[position])
                for position, label in enumerate(target_ambient.module_generating_set())
                if target_coordinates[position]
            }
        )

        radical = target_inclusion.orthogonal_complement()
        if int(radical.module_rank()) != 1:
            raise ValueError("the target subspace must have a one-dimensional isotropic radical")
        radical_vector = radical.inclusion()(radical.module_generators()[0])
        denominator = 2 * target_ambient.b(target_complement, radical_vector)
        if denominator == 0:
            raise ValueError("the target radical does not pair nontrivially with the complementary vector")
        correction = (source_ambient.q(source_complement) - target_ambient.q(target_complement)) / denominator
        target_complement += target_ambient.scalar_multiple(correction, radical_vector)

        source_frame = source_ambient.subobject_on(tuple(source_inclusion(generator) for generator in source_generators) + (source_complement,))
        source_frame_inclusion = source_frame.inclusion()
        target_frame = target_images + (target_complement,)

        def image(label):
            coordinates = source_frame_inclusion.lift(source_ambient.module_generator(label)).to_vector()
            target_labels = tuple(target_ambient.module_generating_set())
            return target_ambient.linear_combination(
                {
                    target_label: coefficient
                    for target_label in target_labels
                    if (
                        coefficient := sum(
                            coordinates(frame_label) * target_frame[index].to_vector()(target_label) for index, frame_label in enumerate(source_frame.module_generating_set())
                        )
                    )
                }
            )

        extension = source_ambient.Isom(target_ambient)(image)
        return CodimensionOneIsotropicExtensionResult(extension, source_inclusion, target_inclusion)


def solve_isotropic_extension_equation(A: Matrix, B: Matrix) -> MatrixEquationSolution:
    r"""Solve \(XA+A^T X^T=B\) exactly over \(\mathbf Q\).

    The map is constructed as a Sage module morphism from the full matrix
    space to its symmetric submodule. Sage supplies its kernel and matrix; the
    latter is used only to obtain one affine preimage of ``B``.
    """
    if A.nrows() != A.ncols():
        raise ValueError("A must be square")
    rank = A.nrows()
    if B.nrows() != rank or B.ncols() != rank:
        raise ValueError("B must have the same shape as A")
    if B != B.transpose():
        raise ValueError("B must be symmetric")

    matrices = MatrixSpace(QQ, rank, rank)
    A = matrices(A)
    B = matrices(B)
    matrix_basis = matrices.basis()
    symmetric_generators = [matrix_basis[(i, i)] for i in range(rank)]
    symmetric_generators.extend(matrix_basis[(i, j)] + matrix_basis[(j, i)] for i in range(rank) for j in range(i + 1, rank))
    symmetric_matrices = matrices.submodule(symmetric_generators)
    phi = matrices.module_morphism(
        on_basis=lambda index: symmetric_matrices.retract(matrix_basis[index] * A + A.transpose() * matrix_basis[index].transpose()),
        codomain=symmetric_matrices,
    )
    target = symmetric_matrices.retract(B)
    target_coefficients = target.monomial_coefficients()
    target_coordinates = vector(
        QQ,
        (target_coefficients.get(label, QQ.zero()) for label in symmetric_matrices.basis().keys()),
    )
    particular_coordinates = phi.matrix().solve_right(target_coordinates)
    particular = matrices(tuple(particular_coordinates))

    return MatrixEquationSolution(
        particular=particular,
        homogeneous_space=phi.kernel(),
        homogeneous_lattice=None,
    )
