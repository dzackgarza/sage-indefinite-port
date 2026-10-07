r"""Exact extension problems along isotropic subspaces."""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass

from dzack_research.preamble.categories.lattice_morphisms import (
    LatticeEmbeddingMethods,
    LatticeIsometryMethods,
    LatticeIsometryMor,
)
from dzack_research.preamble.categories.lattices import IsotropicReductions, Lattices
from dzack_research.preamble.categories.modules.framed.framed_free_modules import (
    FramedFreeModules,
)
from dzack_research.preamble.categories.modules.pure.modules import ModuleSubobjects
from dzack_research.preamble.categories.rings.ring_foundation import (
    _engine_element,
    _engine_ring,
    _owned_engine_element,
)
from dzack_research.preamble.categories.sets.set_categories import OwnedSetMorphism, Sets
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import FreeModule, FreeModule_generic
from sage.modules.free_module_element import FreeModuleElement
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ


@dataclass(frozen=True)
class MatrixEquationSolution:
    r"""Affine solutions of \(XA+A^T X^T=B\) over \(\mathbf Q\)."""

    particular: Matrix_rational_dense
    homogeneous_space: FreeModule_generic[Rational]
    homogeneous_lattice: FreeModule_generic[Integer]


@dataclass(frozen=True)
class CodimensionOneIsotropicExtensionResult:
    r"""A rational extension and the two retained subspace inclusions."""

    extension: LatticeIsometryMethods
    source_inclusion: LatticeEmbeddingMethods
    target_inclusion: LatticeEmbeddingMethods


class CodimensionOneIsotropicExtension:
    r"""Extend an isometry of codimension-one isotropic subspaces over \(\mathbf Q\)."""

    def __init__(
        self,
        source_subspace: ModuleSubobjects.ParentMethods,
        target_subspace: ModuleSubobjects.ParentMethods,
        partial_isometry: LatticeIsometryMethods,
    ) -> None:
        self._source_subspace = source_subspace
        self._target_subspace = target_subspace
        self._partial_isometry = partial_isometry

    def rational_extension(self) -> CodimensionOneIsotropicExtensionResult:
        r"""Return the exact norm-corrected rational extension."""
        source = self._source_subspace
        target = self._target_subspace
        partial = self._partial_isometry
        raw_source_inclusion: object = source.inclusion()
        raw_target_inclusion: object = target.inclusion()
        if not isinstance(raw_source_inclusion, LatticeEmbeddingMethods):
            raise TypeError("the source subspace must carry a lattice embedding")
        if not isinstance(raw_target_inclusion, LatticeEmbeddingMethods):
            raise TypeError("the target subspace must carry a lattice embedding")
        source_inclusion: LatticeEmbeddingMethods = raw_source_inclusion
        target_inclusion: LatticeEmbeddingMethods = raw_target_inclusion
        source_ambient = source_inclusion.codomain()
        target_ambient = target_inclusion.codomain()
        if len(tuple(source.module_generators())) + 1 != int(source_ambient.module_rank()):
            raise ValueError("the source subspace must have codimension one")
        if len(tuple(target.module_generators())) + 1 != int(target_ambient.module_rank()):
            raise ValueError("the target subspace must have codimension one")
        if partial.domain() is not source_inclusion.domain() or partial.codomain() is not target_inclusion.domain():
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

        source_frame = tuple(source_inclusion(generator) for generator in source_generators) + (source_complement,)
        target_frame = target_images + (target_complement,)
        source_labels = tuple(source_ambient.module_generating_set())
        source_ring = source_ambient.base_ring()
        source_frame_matrix = matrix(
            QQ,
            [[_engine_element(source_ring, element.to_vector()(label)) for label in source_labels] for element in source_frame],
        )
        if source_frame_matrix.rank() != len(source_labels):
            raise ArithmeticError("the selected source frame does not span the ambient rational lattice")

        def image(label: Hashable) -> FramedFreeModules.ElementMethods:
            generator = source_ambient.module_generator(label)
            generator_coordinates = vector(
                QQ,
                (_engine_element(source_ring, generator.to_vector()(source_label)) for source_label in source_labels),
            )
            coordinates = source_frame_matrix.transpose().solve_right(generator_coordinates)
            target_labels = tuple(target_ambient.module_generating_set())
            return target_ambient.linear_combination(
                {
                    target_label: coefficient
                    for target_label in target_labels
                    if (
                        coefficient := sum(
                            _owned_engine_element(target_ring, coordinates[index]) * target_frame[index].to_vector()(target_label) for index in range(len(target_frame))
                        )
                    )
                }
            )

        extension = source_ambient.Isom(target_ambient)(image)
        return CodimensionOneIsotropicExtensionResult(extension, source_inclusion, target_inclusion)


class NoIntegralExtensionError(ValueError):
    r"""Raised when a rational extension torsor has no integral member."""


@dataclass(frozen=True)
class IntegralParameterCoset:
    r"""A coset \(\lambda_0+\Lambda\subseteq\mathbf Q^r\)."""

    particular: FreeModuleElement[Rational]
    lattice: FreeModule_generic[Integer]


class IsometryExtensionTorsor:
    r"""An affine family \(T_0+\sum_i\lambda_iT_i\) of rational lattice maps."""

    def __init__(
        self,
        particular: LatticeIsometryMethods,
        homogeneous_directions: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        self._particular = particular
        self._directions = tuple(homogeneous_directions)
        self._integral_locus: IntegralParameterCoset | None = None
        self._integral_source: Lattices.ParentMethods | None = None
        self._integral_target: Lattices.ParentMethods | None = None

    def rational_isometry(self) -> LatticeIsometryMethods:
        r"""Return the selected rational point of this extension torsor."""
        return self._particular

    def integral_parameters(
        self,
        source_lattice: Lattices.ParentMethods,
        target_lattice: Lattices.ParentMethods,
    ) -> IntegralParameterCoset | None:
        r"""Return the exact parameter coset giving integral maps, or ``None``."""
        domain = self._particular.domain()
        codomain = self._particular.codomain()
        domain_labels = tuple(domain.module_generating_set())
        codomain_labels = tuple(codomain.module_generating_set())
        if int(source_lattice.module_rank()) != len(domain_labels) or int(target_lattice.module_rank()) != len(codomain_labels):
            raise ValueError("integral lattices must have the ranks of the rational endpoints")

        rows = []
        for source_position in range(len(domain_labels)):
            rational_generator = domain.module_generator(domain_labels[source_position])
            base_image = self._particular(rational_generator).to_vector()
            direction_images = tuple(direction(rational_generator).to_vector() for direction in self._directions)
            for target_label in codomain_labels:
                rows.append(
                    (
                        _engine_element(codomain.base_ring(), base_image(target_label)),
                        tuple(_engine_element(codomain.base_ring(), image(target_label)) for image in direction_images),
                    )
                )
        constant = vector(QQ, (entry[0] for entry in rows))
        coefficients = matrix(QQ, (entry[1] for entry in rows)) if self._directions else matrix(QQ, len(rows), 0)
        annihilator = coefficients.left_kernel().basis_matrix()
        rhs = annihilator * constant
        denominators = [QQ(entry).denominator() for entry in annihilator.list()] + [QQ(entry).denominator() for entry in rhs.list()]
        common_denominator: Integer = ZZ(1)
        for denominator in denominators:
            common_denominator = common_denominator.lcm(ZZ(denominator))
        equations = matrix(
            ZZ,
            annihilator.nrows(),
            annihilator.ncols(),
            tuple(ZZ(common_denominator * QQ(entry)) for entry in annihilator.list()),
        )
        integral_rhs = vector(
            ZZ,
            tuple(ZZ(common_denominator * QQ(entry)) for entry in rhs.list()),
        )
        smith_data = equations.smith_form()
        if not isinstance(smith_data, tuple):
            raise ArithmeticError("Smith form did not return transformation matrices")
        diagonal, left_change, right_change = smith_data
        transformed_rhs = left_change * integral_rhs
        smith_coordinates = vector(ZZ, [0] * len(rows))
        diagonal_rank = min(diagonal.nrows(), diagonal.ncols())
        for position in range(diagonal.nrows()):
            diagonal_entry = diagonal[position, position] if position < diagonal_rank else 0
            value = transformed_rhs[position]
            if diagonal_entry:
                if value % diagonal_entry:
                    self._integral_locus = None
                    return None
                smith_coordinates[position] = value // diagonal_entry
            elif value:
                self._integral_locus = None
                return None
        integral_point = right_change * smith_coordinates
        parameter_target = vector(
            QQ,
            (
                (
                    QQ(integral_point[position]).numerator() * QQ(constant[position]).denominator()
                    - QQ(constant[position]).numerator() * QQ(integral_point[position]).denominator()
                )
                / (QQ(integral_point[position]).denominator() * QQ(constant[position]).denominator())
                for position in range(len(rows))
            ),
        )
        parameter_point = coefficients.solve_right(parameter_target)
        integral_image_lattice = coefficients.column_space().intersection(FreeModule(ZZ, len(rows)))
        parameter_directions = tuple(coefficients.solve_right(vector(QQ, lattice_vector.list())) for lattice_vector in integral_image_lattice.gens())
        parameter_ambient = FreeModule(QQ, len(self._directions))
        parameter_lattice = parameter_ambient.span(parameter_directions, ZZ)
        locus = IntegralParameterCoset(parameter_point, parameter_lattice)
        self._integral_locus = locus
        self._integral_source = source_lattice
        self._integral_target = target_lattice
        return locus

    def one_integral_extension(self) -> LatticeIsometryMethods:
        r"""Return one integral isometry in the last computed integral locus."""
        if self._integral_locus is None or self._integral_source is None or self._integral_target is None:
            raise NoIntegralExtensionError("the extension torsor has no computed integral member")
        parameters = self._integral_locus.particular
        source = self._integral_source
        target = self._integral_target
        rational_domain = self._particular.domain()
        rational_codomain = self._particular.codomain()
        rational_labels = tuple(rational_domain.module_generating_set())
        target_labels = tuple(target.module_generating_set())
        target_ring = target.base_ring()

        def image(label: Hashable) -> FramedFreeModules.ElementMethods:
            position = tuple(source.module_generating_set()).index(label)
            rational_generator = rational_domain.module_generator(rational_labels[position])
            value = self._particular(rational_generator)
            for coefficient, direction in zip(parameters.list(), self._directions, strict=True):
                value += rational_codomain.scalar_multiple(
                    _owned_engine_element(rational_codomain.base_ring(), coefficient),
                    direction(rational_generator),
                )
            coordinates = value.to_vector()
            return target.linear_combination(
                {
                    target_labels[index]: _owned_engine_element(
                        target_ring,
                        _engine_element(rational_codomain.base_ring(), coordinates(rational_label)),
                    )
                    for index, rational_label in enumerate(rational_codomain.module_generating_set())
                    if coordinates(rational_label)
                }
            )

        extension = source.Isom(target)(image)
        if not isinstance(extension, LatticeIsometryMethods):
            raise ArithmeticError("the integral extension constructor did not return a lattice isometry")
        return extension


class PointwisePerpendicularKernel:
    r"""Integral isometries acting identically on \(I^\perp\).

    For an isotropic reduction \(R=I^\perp/I\), the represented quotient
    \(L/I^\perp\) supplies integral complementary lifts.  If \(U\) is
    their pairing matrix with the chosen framing of \(I\), an integral
    matrix \(H\) changes those lifts by \(H I\).  The form is preserved
    exactly when \(H U^T + U H^T=0\).

    The integral homogeneous lattice returned by the T1 extension-equation
    solver is therefore the parameter lattice.  Addition of parameters maps
    to composition because every image fixes \(I^\perp\) pointwise.
    """

    def __init__(self, reduction: IsotropicReductions.ParentMethods) -> None:
        self._reduction = reduction
        embedding = reduction.isotropic_embedding()
        ambient = embedding.codomain()
        ring = ambient.base_ring()
        if _engine_ring(ring) is not ZZ:
            raise ValueError("the pointwise perpendicular kernel is implemented for integral lattices over ZZ")

        isotropic = reduction.isotropic_sublattice()
        perpendicular = reduction.orthogonal_complement()
        quotient_projection = perpendicular.inclusion().cokernel_projection()
        quotient = quotient_projection.codomain()
        quotient_trivialization = quotient.finite_free_trivialization()
        quotient_free = quotient_trivialization.forward().codomain()
        quotient_section = quotient_projection.section()
        quotient_labels = tuple(quotient_free.module_generating_set())
        isotropic_labels = tuple(isotropic.module_generating_set())
        if len(quotient_labels) != len(isotropic_labels):
            raise ArithmeticError("L/I^perp and I have different ranks, so the nondegenerate pairing between them was not represented correctly")

        complement = tuple(quotient_section(quotient_trivialization.inverse()(quotient_free.module_generator(label))) for label in quotient_labels)
        embedded_isotropic = tuple(embedding(isotropic.module_generator(label)) for label in isotropic_labels)
        pairing = matrix(
            ZZ,
            (tuple(_engine_element(ring, ambient.b(vector, isotropic_vector)) for isotropic_vector in embedded_isotropic) for vector in complement),
        )
        rank = len(isotropic_labels)
        equation = solve_isotropic_extension_equation(
            matrix(QQ, pairing.transpose()),
            matrix(QQ, rank, rank, 0),
        )
        integral_kernel = equation.homogeneous_lattice
        if integral_kernel is None:
            raise ArithmeticError("the homogeneous extension equation did not return its integral solution lattice")
        directions = tuple(matrix(ZZ, rank, rank, basis_vector.list()) for basis_vector in integral_kernel.gens())

        self._ambient = ambient
        self._embedding = embedding
        self._isotropic = isotropic
        self._quotient_projection = quotient_projection
        self._quotient_trivialization = quotient_trivialization
        self._quotient_labels = quotient_labels
        self._isotropic_labels = isotropic_labels
        self._directions = directions
        raw_parameter_lattice: object = ring.free_module(len(directions))
        if not isinstance(raw_parameter_lattice, FramedFreeModules.ParentMethods):
            raise ArithmeticError("the parameter module is not a framed free module")
        self._parameter_lattice: FramedFreeModules.ParentMethods = raw_parameter_lattice

    def parameter_lattice(self) -> FramedFreeModules.ParentMethods:
        r"""Return the free integral lattice of homogeneous solutions \(H\)."""
        return self._parameter_lattice

    def _direction_matrix(self, parameter: FramedFreeModules.ElementMethods) -> Matrix_integer_dense:
        parameter = self._parameter_lattice(parameter)
        coordinates = parameter.to_vector()
        labels = tuple(self._parameter_lattice.module_generating_set())
        rank = len(self._isotropic_labels)
        result = matrix(ZZ, rank, rank, 0)
        for position, label in enumerate(labels):
            coefficient = coordinates(label)
            if coefficient:
                result += _engine_element(self._parameter_lattice.base_ring(), coefficient) * self._directions[position]
        return result

    def _isometry(
        self,
        parameter: FramedFreeModules.ElementMethods,
        target: LatticeIsometryMor,
    ) -> LatticeIsometryMethods:
        parameter = self._parameter_lattice(parameter)
        direction = self._direction_matrix(parameter)
        if not direction:
            return target.identity()

        ambient = self._ambient
        ring = ambient.base_ring()
        quotient_to_free = self._quotient_trivialization.forward()

        def image(label: Hashable) -> FramedFreeModules.ElementMethods:
            source = ambient.module_generator(label)
            quotient_coordinates = quotient_to_free(self._quotient_projection(source)).to_vector()
            shift = self._isotropic.linear_combination(
                {
                    isotropic_label: coefficient
                    for column, isotropic_label in enumerate(self._isotropic_labels)
                    if (
                        coefficient := sum(
                            (
                                quotient_coordinates(quotient_label) * _owned_engine_element(ring, direction[row, column])
                                for row, quotient_label in enumerate(self._quotient_labels)
                            ),
                            ring.zero(),
                        )
                    )
                }
            )
            image_value = source + self._embedding(shift)
            if not isinstance(image_value, Lattices.ElementMethods):
                raise ArithmeticError("the pointwise-kernel image is not a lattice element")
            return image_value

        isometry = target(image)
        if not isinstance(isometry, LatticeIsometryMethods):
            raise ArithmeticError("the pointwise-kernel constructor did not return a lattice isometry")
        return isometry

    def embedding_into(self, target: LatticeIsometryMor) -> OwnedSetMorphism:
        r"""Return the injective parameter map into the ambient orthogonal group."""
        ambient_orthogonal_group = self._ambient.O()
        if target is not ambient_orthogonal_group:
            raise ValueError(f"the pointwise perpendicular kernel embeds in {ambient_orthogonal_group}, not in {target}")
        embedding = Sets().Mor(self._parameter_lattice, target)(lambda parameter: self._isometry(parameter, target))
        if not isinstance(embedding, OwnedSetMorphism):
            raise ArithmeticError("the parameter embedding is not a represented set morphism")
        return embedding

    def gens(self) -> tuple[LatticeIsometryMethods, ...]:
        r"""Return the isometries attached to a basis of the parameter lattice."""
        target = self._ambient.O()
        embedding = self.embedding_into(target)
        return tuple(embedding(generator) for generator in self._parameter_lattice.module_generators())


def pointwise_perpendicular_kernel(
    reduction: IsotropicReductions.ParentMethods,
) -> PointwisePerpendicularKernel:
    r"""Return the exact kernel fixing the reduction's \(I^\perp\) pointwise."""
    return PointwisePerpendicularKernel(reduction)


def solve_isotropic_extension_equation(
    A: Matrix_rational_dense,
    B: Matrix_rational_dense,
) -> MatrixEquationSolution:
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
    A_matrix: Matrix_rational_dense = matrices.matrix(A)
    B_matrix: Matrix_rational_dense = matrices.matrix(B)
    A_transpose: Matrix_rational_dense = matrices.matrix(tuple(A_matrix[column, row] for row in range(rank) for column in range(rank)))
    source_positions = tuple((row, column) for row in range(rank) for column in range(rank))
    symmetric_positions = tuple((i, i) for i in range(rank)) + tuple((i, j) for i in range(rank) for j in range(i + 1, rank))

    def matrix_unit(index: tuple[int, int]) -> Matrix_rational_dense:
        row, column = index
        entries = [QQ.zero()] * (rank * rank)
        entries[row * rank + column] = QQ.one()
        return matrices.matrix(entries)

    image_columns = tuple(
        vector(
            QQ,
            (QQ(image[row, column]) for row, column in symmetric_positions),
        )
        for index in source_positions
        for image in (matrix_unit(index) * A_matrix + A_transpose * matrix_unit((index[1], index[0])),)
    )
    relation_matrix = matrix(
        QQ,
        len(symmetric_positions),
        len(source_positions),
        tuple(QQ(image_columns[column][row]) for row in range(len(symmetric_positions)) for column in range(len(source_positions))),
    )
    if not isinstance(relation_matrix, Matrix_rational_dense):
        raise ArithmeticError("the extension-equation relation matrix is not rational dense")
    target_coordinates = vector(
        QQ,
        (QQ(B_matrix[row, column]) for row, column in symmetric_positions),
    )
    particular_coordinates = relation_matrix.solve_right(target_coordinates)
    particular = matrix(QQ, rank, rank, particular_coordinates.list())
    if not isinstance(particular, Matrix_rational_dense):
        raise ArithmeticError("the extension-equation particular solution is not rational dense")
    homogeneous_space = relation_matrix.right_kernel()
    common_denominator: Integer = ZZ(1)
    for coefficient in relation_matrix.list():
        common_denominator = common_denominator.lcm(QQ(coefficient).denominator())
    integral_relation_matrix = matrix(
        ZZ,
        relation_matrix.nrows(),
        relation_matrix.ncols(),
        tuple(ZZ(common_denominator * QQ(coefficient)) for coefficient in relation_matrix.list()),
    )
    homogeneous_lattice = integral_relation_matrix.right_kernel()

    return MatrixEquationSolution(
        particular=particular,
        homogeneous_space=homogeneous_space,
        homogeneous_lattice=homogeneous_lattice,
    )
