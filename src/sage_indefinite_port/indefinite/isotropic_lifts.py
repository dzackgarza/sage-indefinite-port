r"""Exact extension problems along isotropic subspaces."""

from dataclasses import dataclass

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
