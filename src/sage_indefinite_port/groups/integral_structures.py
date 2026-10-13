"""Integral structures for finitely generated rational isometry groups.

The public objects in this module retain preamble lattices, embeddings, and
isometries.  Sage matrices are used only privately to normalize a finite
``ZZ``-span inside a rational vector space.
"""

from __future__ import annotations

from collections.abc import Callable, Hashable, Iterable, Sequence
from functools import cached_property

from dzack_research.preamble.all import (
    ZZ,
    Lattices,
    Modules,
    ModuleSubobjects,
    RestrictedScalarsModules,
)
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.modules.module_morphisms.module_morphisms import (
    ModuleEmbeddingMethods,
    ModuleMorphismMethods,
)
from dzack_research.preamble.categories.rings.ring_foundation import (
    _engine_element,
    _owned_engine_element,
)
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_rational_dense import Matrix_rational_dense
from sage.modules.vector_integer_dense import Vector_integer_dense
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ
from sage.structure.element import RingElement
from sage.structure.sage_object import SageObject

type RowLatticeKey = tuple[tuple[int, ...], ...]
type IntegralRows = Sequence[Sequence[int | Integer] | Vector_integer_dense]
type RationalMatrixKey = tuple[tuple[Rational, ...], ...]


class RationalMatrixGroup(SageObject):
    """A finitely generated rational isometry group on one preamble lattice.

    The selected generators are actual elements of ``V.Aut()``.  No raw matrix
    is part of the public representation; matrix lowering belongs to private
    computation seams such as the finite integral-structure action below.
    """

    def __init__(
        self,
        rational_lattice: Lattices.ParentMethods,
        generators: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        selected = tuple(generators)
        automorphisms = rational_lattice.Aut()
        match any(generator.parent() is not automorphisms for generator in selected):
            case True:
                raise ValueError("rational-group generators must be automorphisms of one lattice")
            case False:
                pass
        self._rational_lattice = rational_lattice
        self._generators: tuple[LatticeIsometryMethods, ...] | None = selected
        self._selected_column_matrices: tuple[Matrix_rational_dense, ...] | None = None

    @classmethod
    def _from_column_matrices(
        cls,
        rational_lattice: Lattices.ParentMethods,
        column_matrices: tuple[Matrix_rational_dense, ...],
    ) -> RationalMatrixGroup:
        r"""Construct privately from exact generator matrices, materializing isometries lazily."""
        selected = tuple(matrix_.change_ring(SageQQ) for matrix_ in column_matrices)
        rank = int(rational_lattice.module_rank())
        gram = _engine_component_matrix(rational_lattice.gram_tensor()).change_ring(SageQQ)
        for matrix_ in selected:
            if matrix_.nrows() != rank or matrix_.ncols() != rank:
                raise ValueError("a rational-group generator matrix must be square of the lattice rank")
            if matrix_.transpose() * gram * matrix_ != gram:
                raise ValueError("a rational-group generator matrix must preserve the lattice form")
        result = cls.__new__(cls)
        result._rational_lattice = rational_lattice
        result._generators = None
        result._selected_column_matrices = selected
        return result

    def rational_lattice(self) -> Lattices.ParentMethods:
        """Return the rational lattice acted on by this group."""
        return self._rational_lattice

    def generators(self) -> tuple[LatticeIsometryMethods, ...]:
        """Return the selected live isometry generators."""
        match self._generators:
            case None:
                automorphisms = self.rational_lattice().Aut()
                self._generators = tuple(automorphisms._isometry_from_column_matrix(matrix_) for matrix_ in self._generator_column_matrices())
            case _:
                pass
        return tuple(self._generators)

    def _generator_column_matrices(self) -> tuple[Matrix_rational_dense, ...]:
        """Return private exact generator column matrices without forcing live isometries."""
        match self._selected_column_matrices:
            case None:
                automorphisms = self.rational_lattice().Aut()
                self._selected_column_matrices = tuple(automorphisms._row_action_matrix(generator).transpose().change_ring(SageQQ) for generator in self.generators())
            case _:
                pass
        return tuple(self._selected_column_matrices)

    def generators_and_inverses(self) -> tuple[LatticeIsometryMethods, ...]:
        """Return the selected generators and their inverses."""
        generators = self.generators()
        return generators + tuple(~generator for generator in generators)

    def _repr_(self) -> str:
        return f"Rational isometry group on {self.rational_lattice()} generated by {len(self._generator_column_matrices())} elements"


def _embedded_basis(
    inclusion: ModuleEmbeddingMethods,
) -> tuple[RestrictedScalarsModules.ElementMethods, ...]:
    domain = inclusion.domain()
    return tuple(inclusion(generator) for generator in domain.module_generators())


def _underlying_coordinates(
    space: RestrictedScalarsModules.ParentMethods,
    element: RestrictedScalarsModules.ElementMethods,
) -> tuple[Rational, ...]:
    ambient = space.module_over_extension()
    underlying = space(element).underlying_element()
    coefficients = ambient._framing_lift(underlying)
    base_ring = ambient.base_ring()
    return tuple(
        SageQQ(
            _engine_element(
                base_ring,
                coefficients(label),
            )
        )
        for label in ambient.module_generating_set()
    )


def _span_embedding(
    space: RestrictedScalarsModules.ParentMethods,
    elements: Iterable[RestrictedScalarsModules.ElementMethods],
) -> ModuleEmbeddingMethods:
    """Return the exact ``ZZ``-span of finitely many elements of ``Res(V)``."""
    elements = tuple(elements)
    match elements:
        case ():
            raise ValueError("an integral lattice span requires at least one generator")
        case _:
            pass

    rows = tuple(_underlying_coordinates(space, element) for element in elements)
    denominator = SageZZ.one()
    for row in rows:
        for entry in row:
            denominator = denominator.lcm(entry.denominator())
    integral_rows = tuple(tuple(SageZZ(denominator * entry) for entry in row) for row in rows)
    integer_span = matrix(SageZZ, integral_rows).row_module()
    basis = integer_span.basis_matrix()

    ambient = space.module_over_extension()
    ambient_labels = tuple(ambient.module_generating_set())
    domain = ZZ.free_module(basis.nrows())

    def image(position: int) -> RestrictedScalarsModules.ElementMethods:
        row = basis[int(position)]
        coefficients = {}
        for column, label in enumerate(ambient_labels):
            match row[column] == 0:
                case True:
                    pass
                case False:
                    coefficients[label] = ambient.base_ring()(int(row[column])) / ambient.base_ring()(int(denominator))
        return space.wrap(ambient.linear_combination(coefficients))

    return domain.Mono(space)(
        {label: image(label) for label in domain.module_generating_set()},
    )


class ArithmeticSubgroup(RationalMatrixGroup):
    """A generated arithmetic subgroup retaining its ambient rational group."""

    def __init__(
        self,
        supergroup: RationalMatrixGroup,
        generators: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        self._supergroup = supergroup
        super().__init__(supergroup.rational_lattice(), generators)

    def supergroup(self) -> RationalMatrixGroup:
        return self._supergroup


class FinitePermutationRepresentation[FinitePointT: Hashable](SageObject):
    r"""A finite permutation quotient with exact lifts to live isometries."""

    def __init__(
        self,
        group: RationalMatrixGroup,
        points: tuple[FinitePointT, ...],
        action: Callable[[LatticeIsometryMethods, FinitePointT], FinitePointT],
        *,
        generator_permutations: tuple[GapElement, ...] | None = None,
    ) -> None:
        if not points:
            raise ValueError("a finite permutation representation needs at least one point")
        self._group = group
        self._points = tuple(points)
        self._position_by_point = {point: position for position, point in enumerate(self._points)}
        if len(self._position_by_point) != len(self._points):
            raise ValueError("a finite permutation representation needs distinct represented points")
        self._action = action
        self._matrix_isometry_cache: dict[RationalMatrixKey, LatticeIsometryMethods] = {}
        match generator_permutations:
            case None:
                self._generator_permutations = tuple(self._permutation_of(generator, check_image=False) for generator in group.generators())
            case _:
                self._generator_permutations = tuple(generator_permutations)
                if len(self._generator_permutations) != len(group._generator_column_matrices()):
                    raise ValueError("finite generator permutations must match the rational-group generator framing")
        self._free_group = libgap.FreeGroup(len(self._generator_permutations))
        self._permutation_group = libgap.Group(list(self._generator_permutations))
        self._homomorphism = libgap.GroupHomomorphismByImages(
            self._free_group,
            self._permutation_group,
            self._free_group.GeneratorsOfGroup(),
            list(self._generator_permutations),
        )

    def group(self) -> RationalMatrixGroup:
        return self._group

    def points(self) -> tuple[FinitePointT, ...]:
        return self._points

    def is_faithful(self) -> bool:
        r"""Return whether the finite action is faithful on the generated matrix group."""
        automorphisms = self.group().rational_lattice().Aut()
        engine_group = automorphisms._engine_subgroup_from_generators(self.group().generators())
        try:
            generated_order = int(engine_group.order())
        except TypeError, ValueError, OverflowError:
            return False
        return generated_order == self.image_order()

    def image_order(self) -> int:
        return int(self._permutation_group.Size())

    @cached_property
    def _generator_action_matrices(self) -> tuple[Matrix_rational_dense, ...]:
        return self.group()._generator_column_matrices()

    def _generator_word_matrix(
        self,
        word: tuple[tuple[int, int], ...],
    ) -> Matrix_rational_dense:
        r"""Evaluate one generator word in the private exact column matrices."""
        automorphisms = self.group().rational_lattice().Aut()
        matrices = self._generator_action_matrices
        matrix_result = matrices[0].parent().one() if matrices else automorphisms._row_action_matrix(automorphisms.identity()).transpose()
        for generator_position, exponent in word:
            if generator_position < 0 or generator_position >= len(matrices):
                raise ValueError("a generator word refers to a generator outside the selected framing")
            matrix_result = (matrices[generator_position] ** int(exponent)) * matrix_result
        return matrix_result

    def generator_word_lift(
        self,
        word: tuple[tuple[int, int], ...],
    ) -> LatticeIsometryMethods:
        r"""Lift an exponent word in the selected generators to a live isometry."""
        automorphisms = self.group().rational_lattice().Aut()
        matrix_result = self._generator_word_matrix(word)
        matrix_key = tuple(tuple(entry for entry in row) for row in matrix_result.rows())
        cached = self._matrix_isometry_cache.get(matrix_key)
        if cached is not None:
            return cached
        result = automorphisms._isometry_from_column_matrix(matrix_result)
        self._matrix_isometry_cache[matrix_key] = result
        return result

    def _permutation_of(
        self,
        automorphism: LatticeIsometryMethods,
        *,
        check_image: bool = True,
    ) -> GapElement:
        images = []
        for point in self.points():
            moved = self._action(automorphism, point)
            position = self._position_by_point.get(moved)
            if position is None:
                raise ArithmeticError("a finite action does not preserve the represented point set")
            images.append(position + 1)
        permutation = libgap.PermList(images)
        if check_image and permutation not in self._permutation_group:
            raise ValueError("an isometry does not lie in the represented finite image")
        return permutation

    def _gap_word_lift(self, word: GapElement) -> LatticeIsometryMethods:
        representation = tuple(int(entry) for entry in libgap.ExtRepOfObj(word).sage())
        exponent_word = tuple((representation[position] - 1, representation[position + 1]) for position in range(0, len(representation), 2))
        return self.generator_word_lift(exponent_word)

    def _gap_word_matrix(self, word: GapElement) -> Matrix_rational_dense:
        representation = tuple(int(entry) for entry in libgap.ExtRepOfObj(word).sage())
        exponent_word = tuple((representation[position] - 1, representation[position + 1]) for position in range(0, len(representation), 2))
        return self._generator_word_matrix(exponent_word)

    def _preimage_generators(
        self,
        finite_subgroup: GapElement,
        *,
        faithful: bool = False,
    ) -> tuple[LatticeIsometryMethods, ...]:
        if faithful:
            return tuple(self._gap_word_lift(libgap.PreImagesRepresentative(self._homomorphism, generator)) for generator in finite_subgroup.GeneratorsOfGroup())
        preimage = libgap.PreImage(self._homomorphism, finite_subgroup)
        return tuple(self._gap_word_lift(word) for word in preimage.GeneratorsOfGroup())

    def _image_subgroup(
        self,
        generators: tuple[LatticeIsometryMethods, ...],
    ) -> GapElement:
        return libgap.Subgroup(
            self._permutation_group,
            [self._permutation_of(generator) for generator in generators],
        )

    def _subgroup_index(self, finite_subgroup: GapElement) -> int:
        subgroup_order = int(finite_subgroup.Size())
        image_order = self.image_order()
        if image_order % subgroup_order:
            raise ArithmeticError("a finite subgroup order does not divide its ambient image order")
        return image_order // subgroup_order


class GeneratedSubgroup(ArithmeticSubgroup):
    r"""An arithmetic subgroup retaining the selected generating family."""

    def __init__(
        self,
        supergroup: RationalMatrixGroup,
        generators: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        self._selected_generators = tuple(generators)
        super().__init__(supergroup, self._selected_generators)

    def selected_generators(self) -> tuple[LatticeIsometryMethods, ...]:
        return self._selected_generators


class FinitePreimageSubgroup[FinitePointT: Hashable](ArithmeticSubgroup):
    r"""The exact preimage of a represented subgroup of a finite quotient."""

    def __init__(
        self,
        representation: FinitePermutationRepresentation[FinitePointT],
        finite_image_generators: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        self._representation = representation
        self._finite_image_generators = tuple(finite_image_generators)
        self._finite_subgroup = representation._image_subgroup(self._finite_image_generators)
        super().__init__(
            representation.group(),
            representation._preimage_generators(self._finite_subgroup),
        )

    def representation(self) -> FinitePermutationRepresentation[FinitePointT]:
        return self._representation

    def finite_image_generators(self) -> tuple[LatticeIsometryMethods, ...]:
        return self._finite_image_generators

    def index(self) -> int:
        return self.representation()._subgroup_index(self._finite_subgroup)


class KernelSubgroup[FinitePointT: Hashable](ArithmeticSubgroup):
    r"""The exact kernel of a finite permutation representation."""

    def __init__(
        self,
        representation: FinitePermutationRepresentation[FinitePointT],
    ) -> None:
        self._representation = representation
        self._finite_subgroup = libgap.TrivialSubgroup(representation._permutation_group)
        super().__init__(
            representation.group(),
            representation._preimage_generators(self._finite_subgroup),
        )

    def representation(self) -> FinitePermutationRepresentation[FinitePointT]:
        return self._representation

    def index(self) -> int:
        return self.representation()._subgroup_index(self._finite_subgroup)


class StabilizerSubgroup[FinitePointT: Hashable](ArithmeticSubgroup):
    r"""The exact preimage of a point stabilizer in a finite action."""

    def __init__(
        self,
        representation: FinitePermutationRepresentation[FinitePointT],
        point_position: int,
    ) -> None:
        position = int(point_position)
        if position < 0 or position >= len(representation.points()):
            raise ValueError("a finite-action stabilizer needs a represented point position")
        self._representation = representation
        self._point_position = position
        self._finite_subgroup = libgap.Stabilizer(
            representation._permutation_group,
            position + 1,
        )
        faithful = representation.is_faithful()
        super().__init__(
            representation.group(),
            representation._preimage_generators(
                self._finite_subgroup,
                faithful=faithful,
            ),
        )

    def representation(self) -> FinitePermutationRepresentation[FinitePointT]:
        return self._representation

    def point(self) -> FinitePointT:
        return self.representation().points()[self._point_position]

    def index(self) -> int:
        return self.representation()._subgroup_index(self._finite_subgroup)


class CentralizerSubgroup[FinitePointT: Hashable](ArithmeticSubgroup):
    r"""The centralizer of an isometry recovered from a faithful finite action."""

    def __init__(
        self,
        representation: FinitePermutationRepresentation[FinitePointT],
        element: LatticeIsometryMethods,
    ) -> None:
        if not representation.is_faithful():
            raise ValueError("a finite-image centralizer requires a faithful representation")
        self._representation = representation
        self._centralizing_element = element
        self._finite_subgroup = libgap.Centralizer(
            representation._permutation_group,
            representation._permutation_of(element),
        )
        super().__init__(
            representation.group(),
            representation._preimage_generators(
                self._finite_subgroup,
                faithful=True,
            ),
        )

    def representation(self) -> FinitePermutationRepresentation[FinitePointT]:
        return self._representation

    def centralizing_element(self) -> LatticeIsometryMethods:
        return self._centralizing_element

    def index(self) -> int:
        return self.representation()._subgroup_index(self._finite_subgroup)


class IntersectionSubgroup[FinitePointT: Hashable](ArithmeticSubgroup):
    r"""The intersection of arithmetic subgroups in a faithful finite image."""

    def __init__(
        self,
        representation: FinitePermutationRepresentation[FinitePointT],
        subgroups: tuple[ArithmeticSubgroup, ...],
    ) -> None:
        if not representation.is_faithful():
            raise ValueError("finite-image intersection requires a faithful representation")
        selected = tuple(subgroups)
        if not all(subgroup.supergroup() is representation.group() for subgroup in selected):
            raise ValueError("intersected arithmetic subgroups must lie in the represented group")
        self._representation = representation
        self._subgroups = selected
        finite_subgroups = tuple(representation._image_subgroup(subgroup.generators()) for subgroup in selected)
        if finite_subgroups:
            finite = finite_subgroups[0]
            for subgroup in finite_subgroups[1:]:
                finite = libgap.Intersection(finite, subgroup)
            self._finite_subgroup = finite
        else:
            self._finite_subgroup = representation._permutation_group
        super().__init__(
            representation.group(),
            representation._preimage_generators(
                self._finite_subgroup,
                faithful=True,
            ),
        )

    def representation(self) -> FinitePermutationRepresentation[FinitePointT]:
        return self._representation

    def subgroups(self) -> tuple[ArithmeticSubgroup, ...]:
        return self._subgroups

    def index(self) -> int:
        return self.representation()._subgroup_index(self._finite_subgroup)


class RightCosetDecomposition(SageObject):
    """A finite decomposition ``G/H`` with the subgroup explicitly on the right."""

    def __init__(
        self,
        ambient_group: RationalMatrixGroup,
        right_subgroup: ArithmeticSubgroup,
        representatives: tuple[LatticeIsometryMethods, ...],
    ) -> None:
        self._ambient_group = ambient_group
        self._right_subgroup = right_subgroup
        self._representatives = tuple(representatives)

    def ambient_group(self) -> RationalMatrixGroup:
        return self._ambient_group

    def right_subgroup(self) -> ArithmeticSubgroup:
        return self._right_subgroup

    def representatives(self) -> tuple[LatticeIsometryMethods, ...]:
        return self._representatives

    def cardinality(self) -> int:
        return len(self._representatives)

    def _repr_(self) -> str:
        return f"Right-subgroup cosets {self.ambient_group()}/{self.right_subgroup()} with {self.cardinality()} representatives"


class DoubleCosetIntersection(SageObject):
    """Finite-image intersection data for one double coset ``V g H``."""

    def __init__(
        self,
        representative: LatticeIsometryMethods,
        finite_image_order: int,
        double_coset_size: int,
    ) -> None:
        self._representative = representative
        self._finite_image_order = int(finite_image_order)
        self._double_coset_size = int(double_coset_size)

    def representative(self) -> LatticeIsometryMethods:
        return self._representative

    def finite_image_order(self) -> int:
        return self._finite_image_order

    def double_coset_size(self) -> int:
        return self._double_coset_size


class DoubleCosetDecomposition(SageObject):
    r"""A finite decomposition ``V \ G / H`` retaining all three group sides."""

    def __init__(
        self,
        left_subgroup: ArithmeticSubgroup,
        ambient_group: RationalMatrixGroup,
        right_subgroup: ArithmeticSubgroup,
        representatives: tuple[LatticeIsometryMethods, ...],
        intersections: tuple[DoubleCosetIntersection, ...],
        finite_ambient_order: int,
        finite_left_order: int,
        finite_right_order: int,
    ) -> None:
        self._left_subgroup = left_subgroup
        self._ambient_group = ambient_group
        self._right_subgroup = right_subgroup
        self._representatives = tuple(representatives)
        self._intersections = tuple(intersections)
        self._finite_ambient_order = int(finite_ambient_order)
        self._finite_left_order = int(finite_left_order)
        self._finite_right_order = int(finite_right_order)

    def left_subgroup(self) -> ArithmeticSubgroup:
        return self._left_subgroup

    def ambient_group(self) -> RationalMatrixGroup:
        return self._ambient_group

    def right_subgroup(self) -> ArithmeticSubgroup:
        return self._right_subgroup

    def representatives(self) -> tuple[LatticeIsometryMethods, ...]:
        return self._representatives

    def intersections(self) -> tuple[DoubleCosetIntersection, ...]:
        return self._intersections

    def cardinality(self) -> int:
        return len(self._representatives)

    def finite_ambient_order(self) -> int:
        return self._finite_ambient_order

    def finite_left_order(self) -> int:
        return self._finite_left_order

    def finite_right_order(self) -> int:
        return self._finite_right_order

    def _repr_(self) -> str:
        return rf"Double cosets {self.left_subgroup()} \ {self.ambient_group()} / {self.right_subgroup()} with {self.cardinality()} representatives"


class FiniteIntegralRepresentation(SageObject):
    """The finite action controlling the integral structures in one commensurability class.

    The public orbit consists of actual submodules of ``M/dM``.  libGAP owns
    the private permutation group and subgroup preimage calculation; every
    resulting free-group word is evaluated back in the selected live rational
    isometries before it crosses the public boundary.
    """

    def __init__(self, action: IntegralStructureAction) -> None:
        self._action = action
        quotient_exponent = SageZZ(_engine_element(ZZ, action.quotient_exponent()))
        self._prime_modulus = quotient_exponent if quotient_exponent.is_prime() else None
        self._rational_permutation_cache: dict[int, tuple[LatticeIsometryMethods, GapElement]] = {}
        self._generator_ambient_matrices = action.rational_group()._generator_column_matrices()
        self._generator_restricted_matrices = tuple(action._restricted_matrix_from_ambient(matrix_) for matrix_ in self._generator_ambient_matrices)
        (
            self._orbit_keys,
            self._orbit_witness_matrices,
            self._orbit_position_by_key,
        ) = self._compute_submodule_orbit()
        self._orbit: tuple[ModuleSubobjects.ParentMethods, ...] | None = None
        self._orbit_witnesses: tuple[LatticeIsometryMethods, ...] | None = None
        self._finite_action = FinitePermutationRepresentation(
            action.rational_group(),
            tuple(self._orbit_keys),
            self._act_on_key,
            generator_permutations=tuple(self._permutation_of_restricted_matrix(matrix_) for matrix_ in self._generator_restricted_matrices),
        )
        self._generator_permutations = self._finite_action._generator_permutations
        self._permutation_group = self._finite_action._permutation_group
        self._homomorphism = self._finite_action._homomorphism

    def action(self) -> IntegralStructureAction:
        return self._action

    def orbit(self) -> tuple[ModuleSubobjects.ParentMethods, ...]:
        """Return the finite orbit of ``S=L/dM`` as represented submodules."""
        match self._orbit:
            case None:
                self._orbit = tuple(self._submodule_from_key(key) for key in self._orbit_keys)
            case _:
                pass
        return tuple(self._orbit)

    def orbit_witnesses(self) -> tuple[LatticeIsometryMethods, ...]:
        """Return live rational isometries carrying ``S`` to the orbit members."""
        match self._orbit_witnesses:
            case None:
                automorphisms = self.action().rational_group().rational_lattice().Aut()
                self._orbit_witnesses = tuple(automorphisms._isometry_from_column_matrix(witness_matrix) for witness_matrix in self._orbit_witness_matrices)
            case _:
                pass
        return tuple(self._orbit_witnesses)

    def orbit_images(self, element: Lattices.ElementMethods) -> tuple[Lattices.ElementMethods, ...]:
        """Apply every finite-orbit witness to one semantic lattice element.

        The finite action already retains the exact witness column matrices.
        Applying those matrices directly avoids constructing full forward and
        inverse categorical isometries when only their values on one element
        are required. Sage matrices remain private and the boundary returns
        owned lattice elements.
        """
        lattice = self.action().rational_group().rational_lattice()
        labels = tuple(lattice.module_generating_set())
        ring = lattice.base_ring()
        return tuple(
            lattice.linear_combination({label: _owned_engine_element(ring, image_coordinates[row, 0]) for row, label in enumerate(labels) if image_coordinates[row, 0]})
            for image_coordinates in self._orbit_image_coordinate_columns(element)
        )

    def _orbit_image_coordinate_columns(
        self,
        element: Lattices.ElementMethods,
    ) -> tuple[Matrix_rational_dense, ...]:
        r"""Return finite-orbit images as private exact coordinate columns."""
        lattice = self.action().rational_group().rational_lattice()
        source = element if element.parent() is lattice else lattice(element)
        labels = tuple(lattice.module_generating_set())
        coordinates = source.to_vector()
        coordinate_column = matrix(
            SageQQ,
            len(labels),
            1,
            [
                SageQQ(
                    _engine_element(
                        lattice.base_ring(),
                        coordinates(label),
                    )
                )
                for label in labels
            ],
        )
        return tuple(witness_matrix * coordinate_column for witness_matrix in self._orbit_witness_matrices)

    def image_order(self) -> int:
        """Return the order of the finite permutation image."""
        return self._finite_action.image_order()

    @staticmethod
    def _row_lattice_key(
        rows: IntegralRows,
    ) -> RowLatticeKey:
        basis = matrix(SageZZ, rows).row_module().basis_matrix()
        return tuple(tuple(int(entry) for entry in row) for row in basis.rows())

    def _finite_submodule_key(
        self,
        rows: IntegralRows,
    ) -> RowLatticeKey:
        r"""Return the exact canonical key of a submodule of M/dM.

        When d=p is prime, submodules of M/pM are vector subspaces of M/pM.
        Their reduced row spaces over GF(p) are canonical and avoid recomputing
        an integral Hermite basis on every orbit edge. Composite exponents
        retain the integral preimage-lattice key.
        """
        selected_rows = tuple(tuple(entry for entry in row) for row in rows)
        match self._prime_modulus:
            case None:
                return self._row_lattice_key(selected_rows)
            case prime:
                pass
        match selected_rows:
            case ():
                return ()
            case _:
                pass
        modulus = int(prime)
        reduced = [[int(entry) % modulus for entry in row] for row in selected_rows]
        column_count = len(reduced[0])
        pivot_row = 0
        for column in range(column_count):
            pivot = next(
                (position for position in range(pivot_row, len(reduced)) if reduced[position][column]),
                None,
            )
            match pivot:
                case None:
                    continue
                case _:
                    pass
            reduced[pivot_row], reduced[pivot] = reduced[pivot], reduced[pivot_row]
            inverse = pow(reduced[pivot_row][column], -1, modulus)
            reduced[pivot_row] = [(inverse * entry) % modulus for entry in reduced[pivot_row]]
            for position, row in enumerate(reduced):
                if position == pivot_row or row[column] == 0:
                    continue
                coefficient = row[column]
                reduced[position] = [
                    (entry - coefficient * pivot_entry) % modulus
                    for entry, pivot_entry in zip(
                        row,
                        reduced[pivot_row],
                        strict=True,
                    )
                ]
            pivot_row += 1
            if pivot_row == len(reduced):
                break
        return tuple(tuple(row) for row in reduced[:pivot_row])

    def _submodule_key(
        self,
        submodule: ModuleSubobjects.ParentMethods,
    ) -> RowLatticeKey:
        r"""Return the canonical integral preimage lattice of a finite submodule.

        For F = Z^n/R and S <= F generated by rows G, the inverse image of S
        in Z^n is R + <G>.  Its Hermite basis is therefore an exact equality
        key, without enumerating residue vectors.
        """
        inclusion = submodule.inclusion()
        ambient = inclusion.codomain()
        framing = ambient.framing_source()
        framing_labels = tuple(framing.module_generating_set())
        base_ring = framing.base_ring()

        def coordinates(element: Modules.ElementMethods) -> tuple[Integer, ...]:
            vector = framing(element).to_vector()
            return tuple(SageZZ(_engine_element(base_ring, vector(label))) for label in framing_labels)

        presentation = ambient.presentation()
        rows = [coordinates(presentation(generator)) for generator in presentation.domain().module_generators()]
        augmentation = ambient.framing_morphism()
        rows.extend(coordinates(augmentation.lift(inclusion(generator))) for generator in inclusion.domain().module_generators())
        return self._finite_submodule_key(rows)

    @staticmethod
    def _image_key(
        key: RowLatticeKey,
        action_matrix: Matrix_integer_dense,
    ) -> RowLatticeKey:
        moved = matrix(SageZZ, key) * action_matrix
        return FiniteIntegralRepresentation._row_lattice_key(moved.rows())

    def _finite_image_key(
        self,
        key: RowLatticeKey,
        action_matrix: Matrix_integer_dense,
    ) -> RowLatticeKey:
        match key:
            case ():
                return ()
            case _:
                pass
        moved = matrix(SageZZ, key) * action_matrix
        return self._finite_submodule_key(moved.rows())

    def _submodule_from_key(
        self,
        key: RowLatticeKey,
    ) -> ModuleSubobjects.ParentMethods:
        invariant = self.action().invariant_overlattice().domain()
        invariant_labels = tuple(invariant.module_generating_set())
        source = ZZ.free_module(len(key))
        inclusion = source.Mono(invariant)(
            {
                label: invariant.linear_combination({invariant_labels[column]: ZZ(entry) for column, entry in enumerate(key[int(label)]) if entry})
                for label in source.module_generating_set()
            }
        )
        return (self.action().finite_projection() * inclusion).image()

    def _orbit_position(
        self,
        candidate: ModuleSubobjects.ParentMethods,
    ) -> int | None:
        key = self._submodule_key(candidate)
        return self._orbit_position_by_key.get(key)

    def _compute_submodule_orbit(
        self,
    ) -> tuple[
        list[RowLatticeKey],
        list[Matrix_rational_dense],
        dict[RowLatticeKey, int],
    ]:
        orbit_keys = [self._finite_submodule_key(self.action()._selected_invariant_coordinates.rows())]
        orbit_position_by_key = {orbit_keys[0]: 0}
        witness_matrices = [self.action()._invariant_basis_matrix.parent().one()]
        steps = tuple(zip(self._generator_restricted_matrices, self._generator_ambient_matrices, strict=True))
        frontier = [0]
        while frontier:
            source_position = frontier.pop()
            source_key = orbit_keys[source_position]
            source_witness_matrix = witness_matrices[source_position]
            for action_matrix, ambient_action_matrix in steps:
                candidate_key = self._finite_image_key(source_key, action_matrix)
                if candidate_key not in orbit_position_by_key:
                    candidate_position = len(orbit_keys)
                    orbit_keys.append(candidate_key)
                    orbit_position_by_key[candidate_key] = candidate_position
                    witness_matrices.append(ambient_action_matrix * source_witness_matrix)
                    frontier.append(candidate_position)
        return orbit_keys, witness_matrices, orbit_position_by_key

    def _permutation_of_restricted_matrix(
        self,
        action_matrix: Matrix_integer_dense,
    ) -> GapElement:
        images = []
        for key in self._orbit_keys:
            image_key = self._finite_image_key(key, action_matrix)
            position = self._orbit_position_by_key.get(image_key)
            if position is None:
                raise ArithmeticError("a finite-module generator left the computed orbit")
            images.append(position + 1)
        return libgap.PermList(images)

    def _act_on_key(
        self,
        automorphism: LatticeIsometryMethods,
        key: RowLatticeKey,
    ) -> RowLatticeKey:
        action_matrix = self.action()._restricted_matrix(automorphism)
        image_key = self._finite_image_key(key, action_matrix)
        if image_key not in self._orbit_position_by_key:
            raise ArithmeticError("the finite-module generator left the computed orbit")
        return image_key

    def _evaluate_free_word(self, word: GapElement) -> LatticeIsometryMethods:
        return self._finite_action._gap_word_lift(word)

    @cached_property
    def _lattice_stabilizer(self) -> ArithmeticSubgroup:
        point_stabilizer = libgap.Stabilizer(self._permutation_group, 1)
        preimage = libgap.PreImage(self._homomorphism, point_stabilizer)
        lifted_by_identity = {}
        for word in preimage.GeneratorsOfGroup():
            isometry = self._evaluate_free_word(word)
            lifted_by_identity[id(isometry)] = isometry
        lifted = tuple(lifted_by_identity.values())
        match all(self.action().preserves_selected_lattice(generator) for generator in lifted):
            case False:
                raise ArithmeticError("a lifted finite stabilizer generator does not preserve the lattice")
            case True:
                pass
        return ArithmeticSubgroup(self.action().rational_group(), lifted)

    @cached_property
    def _lattice_stabilizer_column_matrices(self) -> tuple[Matrix_rational_dense, ...]:
        r"""Return exact private column matrices for the finite lattice stabilizer generators."""
        point_stabilizer = libgap.Stabilizer(self._permutation_group, 1)
        preimage = libgap.PreImage(self._homomorphism, point_stabilizer)
        by_key = {}
        for word in preimage.GeneratorsOfGroup():
            matrix_ = self._finite_action._gap_word_matrix(word)
            key = tuple(tuple(entry for entry in row) for row in matrix_.rows())
            by_key[key] = matrix_
        return tuple(by_key.values())

    def lattice_stabilizer(self) -> ArithmeticSubgroup:
        """Return the exact preimage of the stabilizer of ``S=L/dM``."""
        return self._lattice_stabilizer

    def transporter_to(
        self,
        target_submodule: ModuleSubobjects.ParentMethods,
    ) -> LatticeIsometryMethods | None:
        """Return a live rational isometry carrying ``S`` to ``target_submodule``."""
        position = self._orbit_position(target_submodule)
        match position:
            case None:
                return None
            case _:
                return self.orbit_witnesses()[position]

    def transporter_to_key(
        self,
        target_key: RowLatticeKey,
    ) -> LatticeIsometryMethods | None:
        """Return a live rational isometry carrying ``S`` to a keyed finite submodule."""
        position = self._orbit_position_by_key.get(target_key)
        match position:
            case None:
                return None
            case _:
                return self.orbit_witnesses()[position]

    @cached_property
    def _right_cosets(self) -> RightCosetDecomposition:
        return RightCosetDecomposition(
            self.action().rational_group(),
            self.lattice_stabilizer(),
            self.orbit_witnesses(),
        )

    def right_cosets(self) -> RightCosetDecomposition:
        """Return ``G/H`` where ``H`` is the selected lattice stabilizer."""
        return self._right_cosets

    def _permutation_of_rational_isometry(
        self,
        automorphism: LatticeIsometryMethods,
    ) -> GapElement:
        cached = self._rational_permutation_cache.get(id(automorphism))
        match cached:
            case (cached_automorphism, permutation) if cached_automorphism is automorphism:
                return permutation
            case _:
                pass
        permutation = self._finite_action._permutation_of(automorphism)
        self._rational_permutation_cache[id(automorphism)] = (automorphism, permutation)
        return permutation

    def double_cosets(self, left_subgroup: ArithmeticSubgroup) -> DoubleCosetDecomposition:
        r"""Return ``left_subgroup \ G / H`` in the finite integral representation."""
        match left_subgroup.supergroup() is self.action().rational_group():
            case False:
                raise ValueError("a double-coset left subgroup must lie in the selected rational group")
            case True:
                pass
        left_image = libgap.Subgroup(
            self._permutation_group,
            [self._permutation_of_rational_isometry(generator) for generator in left_subgroup.generators()],
        )
        right_image = libgap.Stabilizer(self._permutation_group, 1)
        representatives = []
        intersections = []
        for double_coset in libgap.DoubleCosets(
            self._permutation_group,
            left_image,
            right_image,
        ):
            finite_representative = double_coset.Representative()
            word = libgap.PreImagesRepresentative(
                self._homomorphism,
                finite_representative,
            )
            match str(word) == "fail":
                case True:
                    raise ArithmeticError("a finite double-coset representative has no rational lift")
                case False:
                    pass
            representative = self._evaluate_free_word(word)
            double_coset_size = int(double_coset.Size())
            numerator = int(left_image.Size()) * int(right_image.Size())
            match numerator % double_coset_size:
                case 0:
                    intersection_order = numerator // double_coset_size
                case _:
                    raise ArithmeticError("the finite double-coset size violates the orbit-stabilizer formula")
            representatives.append(representative)
            intersections.append(
                DoubleCosetIntersection(
                    representative,
                    intersection_order,
                    double_coset_size,
                )
            )
        return DoubleCosetDecomposition(
            left_subgroup,
            self.action().rational_group(),
            self.lattice_stabilizer(),
            tuple(representatives),
            tuple(intersections),
            int(self._permutation_group.Size()),
            int(left_image.Size()),
            int(right_image.Size()),
        )


class IntegralStructureAction(SageObject):
    """A rational group acting on lattices commensurable with one selected lattice.

    ``lattice_inclusion`` is the actual monomorphism ``L -> Res_Z^Q(V)``.  The
    invariant-overlattice procedure has the explicit precondition from the
    source algorithm: the selected rational group stabilizes some lattice
    commensurable with ``L``.  Under that hypothesis the ascending sequence of
    finite ``ZZ``-spans stabilizes.
    """

    def __init__(
        self,
        rational_group: RationalMatrixGroup,
        lattice_inclusion: ModuleEmbeddingMethods,
    ) -> None:
        space = lattice_inclusion.codomain()
        match space.module_over_extension() is rational_group.rational_lattice():
            case False:
                raise ValueError("the integral structure must lie in the rational group's space")
            case True:
                pass
        match lattice_inclusion.domain().base_ring() is space.base_ring():
            case False:
                raise ValueError("the lattice and restricted rational space require one base ring")
            case True:
                pass
        self._rational_group = rational_group
        self._lattice_inclusion = lattice_inclusion
        self._ambient_action_matrix_cache: dict[int, tuple[LatticeIsometryMethods, Matrix_rational_dense]] = {}
        self._restricted_matrix_cache: dict[int, tuple[LatticeIsometryMethods, Matrix_integer_dense]] = {}

    def rational_group(self) -> RationalMatrixGroup:
        return self._rational_group

    def lattice_inclusion(self) -> ModuleEmbeddingMethods:
        """Return the selected lattice embedding ``L -> Res(V)``."""
        return self._lattice_inclusion

    def _ambient_action_matrix(
        self,
        automorphism: LatticeIsometryMethods,
    ) -> Matrix_rational_dense:
        cached = self._ambient_action_matrix_cache.get(id(automorphism))
        match cached:
            case (cached_automorphism, action_matrix) if cached_automorphism is automorphism:
                return action_matrix
            case _:
                pass

        space = self.lattice_inclusion().codomain()
        ambient = space.module_over_extension()
        labels = tuple(ambient.module_generating_set())
        basis = tuple(ambient.module_generator(label) for label in labels)
        columns = tuple(
            _underlying_coordinates(
                space,
                space.wrap(automorphism(vector)),
            )
            for vector in basis
        )
        result = matrix(
            SageQQ,
            [[columns[column][row] for column in range(len(columns))] for row in range(len(labels))],
        )
        self._ambient_action_matrix_cache[id(automorphism)] = (automorphism, result)
        return result

    def _restricted_matrix_from_ambient(
        self,
        ambient_action_matrix: Matrix_rational_dense,
    ) -> Matrix_integer_dense:
        basis = self._invariant_basis_matrix
        moved = basis * ambient_action_matrix.transpose()
        coordinates = moved * basis.inverse()
        match all(entry.denominator() == 1 for entry in coordinates.list()):
            case False:
                raise ArithmeticError("a rational-group generator does not preserve the computed invariant over-lattice")
            case True:
                pass
        return matrix(
            SageZZ,
            [[SageZZ(entry) for entry in row] for row in coordinates.rows()],
        )

    @cached_property
    def _invariant_overlattice(self) -> ModuleEmbeddingMethods:
        selected = self.lattice_inclusion()
        space = selected.codomain()
        ambient = space.module_over_extension()
        ambient_labels = tuple(ambient.module_generating_set())
        basis = matrix(
            SageQQ,
            tuple(_underlying_coordinates(space, element) for element in _embedded_basis(selected)),
        )
        match basis.nrows() == basis.ncols() and basis.rank() == basis.ncols():
            case False:
                raise ValueError("an integral structure must be a full-rank lattice in the selected rational space")
            case True:
                pass

        generator_matrices = self.rational_group()._generator_column_matrices()
        action_matrices = generator_matrices + tuple(action_matrix.inverse() for action_matrix in generator_matrices)

        changed = False
        while True:
            inverse_basis: Matrix_rational_dense = basis.inverse()
            spanning_rows = list(basis.rows())
            stable = True
            for action_matrix in action_matrices:
                moved = basis * action_matrix.transpose()
                spanning_rows.extend(moved.rows())
                for row in moved.rows():
                    coordinates = row * inverse_basis
                    if any(entry.denominator() != 1 for entry in coordinates):
                        stable = False

            match stable:
                case True:
                    match changed:
                        case False:
                            return selected
                        case True:
                            elements = []
                            for row in basis.rows():
                                coefficients = {}
                                for column, label in enumerate(ambient_labels):
                                    entry = row[column]
                                    match entry == 0:
                                        case True:
                                            pass
                                        case False:
                                            coefficients[label] = ambient.base_ring()(int(entry.numerator())) / ambient.base_ring()(int(entry.denominator()))
                                elements.append(space.wrap(ambient.linear_combination(coefficients)))
                            return _span_embedding(space, elements)
                case False:
                    pass

            denominator = SageZZ.one()
            for row in spanning_rows:
                for entry in row:
                    denominator = denominator.lcm(entry.denominator())
            integral_rows = tuple(tuple(SageZZ(denominator * entry) for entry in row) for row in spanning_rows)
            integral_basis = matrix(SageZZ, integral_rows).row_module().basis_matrix()
            basis = matrix(
                SageQQ,
                [[SageQQ(entry) / SageQQ(denominator) for entry in row] for row in integral_basis.rows()],
            )
            changed = True

    def invariant_overlattice(self) -> ModuleEmbeddingMethods:
        """Return the smallest represented over-lattice stable under the generators.

        This is the source operation ``MatrixIntegral_GetInvariantSpace``:
        repeatedly take the ``ZZ``-span of the current lattice and all generator
        and inverse images.  No saturation is taken inside the divisible
        rational ambient space.
        """
        return self._invariant_overlattice

    @cached_property
    def _invariant_basis_matrix(self) -> Matrix_rational_dense:
        space = self.invariant_overlattice().codomain()
        return matrix(
            SageQQ,
            tuple(_underlying_coordinates(space, element) for element in _embedded_basis(self.invariant_overlattice())),
        )

    @cached_property
    def _selected_invariant_coordinates(self) -> Matrix_integer_dense:
        space = self.lattice_inclusion().codomain()
        selected_basis = matrix(
            SageQQ,
            tuple(_underlying_coordinates(space, element) for element in _embedded_basis(self.lattice_inclusion())),
        )
        coordinates = selected_basis * self._invariant_basis_matrix.inverse()
        match all(entry.denominator() == 1 for entry in coordinates.list()):
            case False:
                raise ArithmeticError("the selected lattice is not integral in its invariant over-lattice")
            case True:
                pass
        return matrix(
            SageZZ,
            [[SageZZ(entry) for entry in row] for row in coordinates.rows()],
        )

    @cached_property
    def _quotient_exponent(self) -> RingElement:
        diagonal, _left, _right = self._selected_invariant_coordinates.smith_form()
        exponent = SageZZ.one()
        for entry in diagonal.diagonal():
            match entry:
                case 0:
                    raise ArithmeticError("the selected lattice must have finite index in its invariant over-lattice")
                case _:
                    exponent = exponent.lcm(abs(SageZZ(entry)))
        return ZZ(int(exponent))

    def quotient_exponent(self) -> RingElement:
        """Return the exponent of ``M/L`` for ``M`` the invariant over-lattice."""
        return self._quotient_exponent

    @cached_property
    def _scaling_morphism(self) -> ModuleMorphismMethods:
        invariant = self.invariant_overlattice().domain()
        modulus = self.quotient_exponent()
        return invariant.Mor(invariant)({label: invariant.scalar_multiple(modulus, invariant.module_generator(label)) for label in invariant.module_generating_set()})

    def finite_module(self) -> Modules.ParentMethods:
        """Return ``F=M/dM`` as the actual cokernel of multiplication by ``d``."""
        return self._scaling_morphism.cokernel()

    def finite_projection(self) -> ModuleMorphismMethods:
        """Return the quotient map ``M -> M/dM``."""
        return self._scaling_morphism.cokernel_projection()

    def intermediate_image(
        self,
        lattice_inclusion: ModuleEmbeddingMethods,
    ) -> ModuleSubobjects.ParentMethods:
        """Return ``L'/dM`` inside ``M/dM`` for ``dM <= L' <= M``."""
        invariant = self.invariant_overlattice()
        match lattice_inclusion.codomain() is invariant.codomain():
            case False:
                raise ValueError("an intermediate lattice lies in the selected rational space")
            case True:
                pass
        into_invariant = lattice_inclusion.factor_through(invariant)
        scaled_invariant = invariant * self._scaling_morphism
        scaled_invariant.factor_through(lattice_inclusion)
        return (self.finite_projection() * into_invariant).image()

    def _intermediate_image_key(
        self,
        lattice_inclusion: ModuleEmbeddingMethods,
    ) -> RowLatticeKey:
        r"""Return the canonical preimage-row-lattice key of ``L'/dM <= M/dM``."""
        invariant = self.invariant_overlattice()
        match lattice_inclusion.codomain() is invariant.codomain():
            case False:
                raise ValueError("an intermediate lattice lies in the selected rational space")
            case True:
                pass
        space = invariant.codomain()
        basis = matrix(
            SageQQ,
            tuple(_underlying_coordinates(space, element) for element in _embedded_basis(lattice_inclusion)),
        )
        coordinates = basis * self._invariant_basis_matrix.inverse()
        match all(entry.denominator() == 1 for entry in coordinates.list()):
            case False:
                raise ValueError("an intermediate lattice is not contained in the invariant over-lattice")
            case True:
                pass
        integral_coordinates = matrix(
            SageZZ,
            [[SageZZ(entry) for entry in row] for row in coordinates.rows()],
        )
        exponent = SageZZ(_engine_element(ZZ, self.quotient_exponent()))
        scaled = exponent * matrix.identity(SageZZ, integral_coordinates.ncols())
        return self.finite_representation()._finite_submodule_key(tuple(integral_coordinates.rows()) + tuple(scaled.rows()))

    @cached_property
    def _selected_submodule(self) -> ModuleSubobjects.ParentMethods:
        return self.intermediate_image(self.lattice_inclusion())

    def selected_submodule(self) -> ModuleSubobjects.ParentMethods:
        """Return ``S=L/dM <= M/dM`` for the selected lattice."""
        return self._selected_submodule

    def _restricted_matrix(
        self,
        automorphism: LatticeIsometryMethods,
    ) -> Matrix_integer_dense:
        cached = self._restricted_matrix_cache.get(id(automorphism))
        match cached:
            case (cached_automorphism, restricted_matrix) if cached_automorphism is automorphism:
                return restricted_matrix
            case _:
                pass

        result = self._restricted_matrix_from_ambient(self._ambient_action_matrix(automorphism))
        self._restricted_matrix_cache[id(automorphism)] = (automorphism, result)
        return result

    def preserves_selected_lattice(self, automorphism: LatticeIsometryMethods) -> bool:
        """Return whether the selected rational isometry carries ``L`` onto itself."""
        selected_key = FiniteIntegralRepresentation._row_lattice_key(self._selected_invariant_coordinates.rows())
        image_key = FiniteIntegralRepresentation._image_key(
            selected_key,
            self._restricted_matrix(automorphism),
        )
        return image_key == selected_key

    @cached_property
    def _finite_representation(self) -> FiniteIntegralRepresentation:
        return FiniteIntegralRepresentation(self)

    def finite_representation(self) -> FiniteIntegralRepresentation:
        """Return the finite submodule-orbit representation controlling integrality."""
        return self._finite_representation

    def lattice_stabilizer(self) -> ArithmeticSubgroup:
        """Return ``{g in G : g(L)=L}`` from the finite representation."""
        return self.finite_representation().lattice_stabilizer()

    def _carries_lattice(
        self,
        automorphism: LatticeIsometryMethods,
        source_inclusion: ModuleEmbeddingMethods,
        target_inclusion: ModuleEmbeddingMethods,
    ) -> bool:
        match source_inclusion.codomain() is target_inclusion.codomain():
            case False:
                return False
            case True:
                pass
        space = source_inclusion.codomain()
        source_basis = matrix(
            SageQQ,
            tuple(_underlying_coordinates(space, element) for element in _embedded_basis(source_inclusion)),
        )
        target_basis = matrix(
            SageQQ,
            tuple(_underlying_coordinates(space, element) for element in _embedded_basis(target_inclusion)),
        )
        moved_basis = source_basis * self._ambient_action_matrix(automorphism).transpose()
        denominator = SageZZ.one()
        for entry in tuple(moved_basis.list()) + tuple(target_basis.list()):
            denominator = denominator.lcm(entry.denominator())
        moved_key = FiniteIntegralRepresentation._row_lattice_key(tuple(tuple(SageZZ(denominator * entry) for entry in row) for row in moved_basis.rows()))
        target_key = FiniteIntegralRepresentation._row_lattice_key(tuple(tuple(SageZZ(denominator * entry) for entry in row) for row in target_basis.rows()))
        return moved_key == target_key

    def transporter(
        self,
        source_inclusion: ModuleEmbeddingMethods,
        target_inclusion: ModuleEmbeddingMethods,
    ) -> LatticeIsometryMethods | None:
        """Return one ``g`` in the rational group with ``g(source)=target``.

        The finite submodule ``source/dM`` determines the source lattice among
        the intermediate lattices ``dM <= source <= M``.  An orbit witness in
        ``M/dM`` therefore lifts to a genuine lattice transporter, which is
        checked again on both actual embeddings before being returned.
        """
        selected = self.lattice_inclusion()
        match source_inclusion is selected, target_inclusion is selected:
            case True, _:
                try:
                    target_submodule = self.intermediate_image(target_inclusion)
                except ValueError:
                    return None
                witness = self.finite_representation().transporter_to(target_submodule)
                match witness:
                    case None:
                        return None
                    case _:
                        pass
                match self._carries_lattice(witness, source_inclusion, target_inclusion):
                    case False:
                        raise ArithmeticError("a finite-module transporter does not carry the actual lattice")
                    case True:
                        return witness
            case _, True:
                try:
                    source_key = self._intermediate_image_key(source_inclusion)
                except ValueError:
                    return None
                selected_to_source = self.finite_representation().transporter_to_key(source_key)
                match selected_to_source:
                    case None:
                        return None
                    case _:
                        pass
                witness = ~selected_to_source
                match self._carries_lattice(witness, source_inclusion, target_inclusion):
                    case False:
                        raise ArithmeticError("a finite-module transporter does not carry the actual lattice")
                    case True:
                        return witness
            case _:
                pass

        source_action = IntegralStructureAction(self.rational_group(), source_inclusion)
        try:
            target_submodule = source_action.intermediate_image(target_inclusion)
        except ValueError:
            return None
        witness = source_action.finite_representation().transporter_to(target_submodule)
        match witness:
            case None:
                return None
            case _:
                pass
        match self._carries_lattice(witness, source_inclusion, target_inclusion):
            case False:
                raise ArithmeticError("a finite-module transporter does not carry the actual lattice")
            case True:
                return witness

    def right_cosets(self) -> RightCosetDecomposition:
        """Return the finite ``G/H`` decomposition for the selected stabilizer ``H``."""
        return self.finite_representation().right_cosets()

    def double_cosets(
        self,
        left_subgroup: ArithmeticSubgroup,
    ) -> DoubleCosetDecomposition:
        r"""Return ``left_subgroup \ G / G_L`` with all sides retained."""
        return self.finite_representation().double_cosets(left_subgroup)

    def double_cosets_from_generators(
        self,
        left_generators: tuple[LatticeIsometryMethods, ...],
    ) -> DoubleCosetDecomposition:
        r"""Return ``V \ G / G_L`` for the subgroup generated by live ``V`` generators."""
        return self.double_cosets(ArithmeticSubgroup(self.rational_group(), tuple(left_generators)))

    def _repr_(self) -> str:
        return f"Integral-structure action of {self.rational_group()} on {self.lattice_inclusion().domain()}"


def integral_structure_action(
    rational_lattice: Lattices.ParentMethods,
    generators: Iterable[LatticeIsometryMethods],
    lattice_inclusion: ModuleEmbeddingMethods,
) -> IntegralStructureAction:
    """Build the T2 integral-structure action from live preamble objects."""
    return IntegralStructureAction(
        RationalMatrixGroup(rational_lattice, tuple(generators)),
        lattice_inclusion,
    )


def integral_structure_action_for_group(
    rational_group: RationalMatrixGroup,
    lattice_inclusion: ModuleEmbeddingMethods,
) -> IntegralStructureAction:
    """Build the T2 action while retaining the selected rational-group parent."""
    return IntegralStructureAction(rational_group, lattice_inclusion)


__all__ = [
    "ArithmeticSubgroup",
    "CentralizerSubgroup",
    "DoubleCosetDecomposition",
    "DoubleCosetIntersection",
    "FiniteIntegralRepresentation",
    "FinitePermutationRepresentation",
    "FinitePreimageSubgroup",
    "GeneratedSubgroup",
    "IntegralStructureAction",
    "IntersectionSubgroup",
    "KernelSubgroup",
    "integral_structure_action",
    "integral_structure_action_for_group",
    "RationalMatrixGroup",
    "RightCosetDecomposition",
    "StabilizerSubgroup",
]
