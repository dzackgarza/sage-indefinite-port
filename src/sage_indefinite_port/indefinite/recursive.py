"""Public recursive indefinite-lattice operations on live preamble objects."""

from random import Random

from dzack_research.preamble.all import QQ, ZZ, Lattices, Modules
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.modules.free_module_element import vector as sage_vector
from sage.rings.integer_ring import ZZ as SageZZ
from sage.sets.disjoint_set import DisjointSet as SageDisjointSet

from sage_indefinite_port.backends.canonization import (
    _identity_rows,
    _rational_lattice_with_integral_structure,
)
from sage_indefinite_port.groups.integral_structures import (
    GeneratedSubgroup,
    IntegralStructureAction,
    RationalMatrixGroup,
)
from sage_indefinite_port.indefinite.eichler import OrbitCoverModel
from sage_indefinite_port.indefinite.lorentzian_cells import (
    IndefiniteOrthogonalAlgorithm as _RecursiveBackend,
)
from sage_indefinite_port.indefinite.vector_sections import (
    IsotropicVectorSection,
    NonIsotropicVectorSection,
    orthogonal_section,
)
from sage_indefinite_port.readiness import unfinished
from sage_indefinite_port.invariants import (
    LatticePrefilter,
    VectorPrefilter,
    lattice_prefilter,
    vector_content,
)


class IndefiniteOrthogonalAlgorithm(_RecursiveBackend):
    r"""Exact recursion above the definite and Lorentzian leaf backends."""

    def __init__(self) -> None:
        super().__init__()
        self._orbit_cover_model_cache: dict[
            Lattices.ParentMethods,
            OrbitCoverModel,
        ] = {}
        self._reduced_framing_cache: dict[
            tuple[tuple[int, ...], ...],
            Lattices.ParentMethods,
        ] = {}
        self._simple_reduction_cache: dict[
            int,
            tuple[
                Lattices.ParentMethods,
                Lattices.ParentMethods,
                Matrix_integer_dense,
            ],
        ] = {}
        self._vector_orbit_cache: dict[
            tuple[Lattices.ParentMethods, int],
            tuple[Lattices.ElementMethods, ...],
        ] = {}
        self._section_cache: dict[
            int,
            tuple[
                Lattices.ElementMethods,
                IsotropicVectorSection | NonIsotropicVectorSection,
            ],
        ] = {}
        self._lattice_prefilter_cache: dict[
            int,
            tuple[Lattices.ParentMethods, LatticePrefilter],
        ] = {}
        self._isometry_cache: dict[
            tuple[int, int],
            tuple[
                Lattices.ParentMethods,
                Lattices.ParentMethods,
                LatticeIsometryMethods | None,
            ],
        ] = {}
        self._transporter_cache: dict[
            tuple[int, int],
            tuple[
                Lattices.ElementMethods,
                Lattices.ElementMethods,
                LatticeIsometryMethods | None,
            ],
        ] = {}
        self._split_isotropic_data_cache: dict[
            int,
            tuple[
                Lattices.ElementMethods,
                Lattices.ParentMethods,
                Matrix_integer_dense,
            ],
        ] = {}

    def _section(
        self,
        vector: Lattices.ElementMethods,
    ) -> IsotropicVectorSection | NonIsotropicVectorSection:
        cached = self._section_cache.get(id(vector))
        if cached is not None and cached[0] is vector:
            return cached[1]
        section = orthogonal_section(vector)
        self._section_cache[id(vector)] = (vector, section)
        return section

    def _lattice_prefilter(
        self,
        lattice: Lattices.ParentMethods,
    ) -> LatticePrefilter:
        cached = self._lattice_prefilter_cache.get(id(lattice))
        if cached is not None and cached[0] is lattice:
            return cached[1]
        result = lattice_prefilter(lattice)
        self._lattice_prefilter_cache[id(lattice)] = (lattice, result)
        return result

    def _orbit_cover_model(
        self,
        lattice: Lattices.ParentMethods,
    ) -> OrbitCoverModel:
        match self._orbit_cover_model_cache.get(lattice):
            case OrbitCoverModel() as model:
                return model
            case None:
                model = OrbitCoverModel.from_lattice(lattice)
                self._orbit_cover_model_cache[lattice] = model
                return model

    def _simple_indefinite_reduction_data(
        self,
        lattice: Lattices.ParentMethods,
    ) -> tuple[Lattices.ParentMethods, Matrix_integer_dense]:
        r"""Apply the source SimpleIndefiniteReduction elementary moves.

        Each accepted move is an integral unimodular row operation
        P = I + c E_ij which strictly decreases the L1 norm of the Gram
        matrix under G -> P G P^T.  The returned matrix is the reference
        algorithm's accumulated basis change from the reduced framing to the
        input framing.
        """
        cached = self._simple_reduction_cache.get(id(lattice))
        if cached is not None and cached[0] is lattice:
            return cached[1], cached[2]
        ring = lattice.base_ring()
        gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(SageZZ)
        rank = gram.nrows()
        transformation = identity_matrix(SageZZ, rank)
        work = gram
        indices = list(range(rank))
        search_rng = Random(0)

        while True:
            selected = None
            first = search_rng.randrange(rank)
            second = search_rng.randrange(rank)
            if first != second:
                indices[first], indices[second] = indices[second], indices[first]
            for row in indices:
                if selected is not None:
                    break
                current_affected_norm = abs(work[row, row]) + 2 * sum(abs(work[row, column]) for column in range(rank) if column != row)
                for source_row in indices:
                    if row == source_row:
                        continue
                    for direction in (-1, 1):
                        coefficient = direction
                        best_coefficient = None
                        best_delta = SageZZ.zero()
                        while True:
                            candidate_diagonal = work[row, row] + 2 * coefficient * work[row, source_row] + coefficient * coefficient * work[source_row, source_row]
                            candidate_affected_norm = abs(candidate_diagonal) + 2 * sum(
                                abs(work[row, column] + coefficient * work[source_row, column]) for column in range(rank) if column != row
                            )
                            delta = current_affected_norm - candidate_affected_norm
                            if delta <= best_delta:
                                break
                            best_delta = delta
                            best_coefficient = coefficient
                            coefficient += direction
                        if best_coefficient is not None:
                            selected = (row, source_row, best_coefficient)
                            break
                    if selected is not None:
                        break
            if selected is None:
                break
            row, source_row, coefficient = selected
            elementary = identity_matrix(SageZZ, rank)
            elementary[row, source_row] = coefficient
            transformation = elementary * transformation
            work = elementary * work * elementary.transpose()

        reduction_key = tuple(tuple(int(work[row, column]) for column in range(rank)) for row in range(rank))
        reduced = self._reduced_framing_cache.get(reduction_key)
        match reduced:
            case None:
                match transformation.is_one():
                    case True:
                        reduced = lattice
                    case False:
                        reduced = Lattices(ring)([[ring(int(work[row, column])) for column in range(rank)] for row in range(rank)])
                self._reduced_framing_cache[reduction_key] = reduced
            case _:
                pass
        self._simple_reduction_cache[id(lattice)] = (
            lattice,
            reduced,
            transformation,
        )
        return reduced, transformation

    def _simple_indefinite_reduction(
        self,
        lattice: Lattices.ParentMethods,
    ) -> tuple[Lattices.ParentMethods, LatticeIsometryMethods | None]:
        reduced, transformation = self._simple_indefinite_reduction_data(lattice)
        if reduced is lattice:
            return lattice, None

        ring = lattice.base_ring()
        source_labels = tuple(lattice.module_generating_set())
        reduced_labels = tuple(reduced.module_generating_set())

        def image(label):
            position = reduced_labels.index(label)
            return lattice.linear_combination(
                {source_label: ring(int(transformation[position, column])) for column, source_label in enumerate(source_labels) if transformation[position, column]}
            )

        reduced_to_source = reduced.Isom(lattice)(image)
        return reduced, reduced_to_source

    @staticmethod
    def _retarget_isometry(
        isometry: LatticeIsometryMethods,
        source: Lattices.ParentMethods,
        target: Lattices.ParentMethods,
    ) -> LatticeIsometryMethods:
        match source is target:
            case True:
                expected_parent = source.Aut()
            case False:
                expected_parent = source.Isom(target)
        match (
            isometry.domain() is source,
            isometry.codomain() is target,
            isometry.parent() is expected_parent,
        ):
            case (True, True, True):
                return isometry
            case _:
                pass
        source_labels = tuple(source.module_generating_set())
        target_labels = tuple(target.module_generating_set())
        rational_source_labels = tuple(isometry.domain().module_generating_set())
        rational_target_labels = tuple(isometry.codomain().module_generating_set())

        def image(label):
            position = source_labels.index(label)
            moved = isometry(isometry.domain().module_generator(rational_source_labels[position]))
            coordinates = moved.to_vector()
            return target.linear_combination(
                {
                    target_label: target.base_ring()(coordinates(rational_label))
                    for target_label, rational_label in zip(target_labels, rational_target_labels, strict=True)
                    if coordinates(rational_label)
                }
            )

        match source is target:
            case True:
                return expected_parent(image)
            case False:
                return expected_parent(image)

    @staticmethod
    def _selected_generators(group) -> tuple[LatticeIsometryMethods, ...]:
        match group:
            case RationalMatrixGroup():
                return tuple(group.generators())
            case _:
                return tuple(group.framing().group_generators())

    def orthogonal_group(self, lattice: Lattices.ParentMethods) -> GeneratedSubgroup:
        if lattice.is_definite():
            from sage_indefinite_port.backends.definite import definite_orthogonal_group

            group = definite_orthogonal_group(lattice)
            generators = self._selected_generators(group)
            if generators and generators[0].domain() is not lattice:
                generators = tuple(self._retarget_isometry(generator, lattice, lattice) for generator in generators)
            return GeneratedSubgroup(RationalMatrixGroup(lattice, generators), generators)
        profile = self.attack_profile(lattice)
        reduced, reduced_to_lattice = self._simple_indefinite_reduction(lattice)
        if reduced is not lattice:
            assert reduced_to_lattice is not None
            reduced_group = self.orthogonal_group(reduced)
            generators = tuple(reduced_to_lattice * generator * ~reduced_to_lattice for generator in reduced_group.generators())
            return GeneratedSubgroup(
                RationalMatrixGroup(lattice, generators),
                generators,
            )
        profile = self.attack_profile(lattice)
        signed_lattice = profile.signed_view
        match profile.positive_index:
            case 0 | 1:
                group = super().orthogonal_group(signed_lattice)
                generators = self._selected_generators(group)
            case _:
                model = self._orbit_cover_model(signed_lattice)
                vector = model.choose_splitting_vector()
                approximate = tuple(model.subgroup().generators())
                stabilizer = tuple(self.vector_stabilizer(vector).generators())
                transporters = tuple(
                    witness
                    for candidate in model.covering_representatives(
                        vector.q(),
                        primitive=vector.is_primitive(),
                    )
                    if (witness := self.vector_transporter(vector, candidate)) is not None
                )
                generators = approximate + stabilizer + transporters
        if generators and generators[0].domain() is not lattice:
            generators = tuple(self._retarget_isometry(generator, lattice, lattice) for generator in generators)
        return GeneratedSubgroup(RationalMatrixGroup(lattice, generators), generators)

    def isometry(
        self,
        source: Lattices.ParentMethods,
        target: Lattices.ParentMethods,
    ) -> LatticeIsometryMethods | None:
        cache_key = (id(source), id(target))
        cached = self._isometry_cache.get(cache_key)
        if cached is not None and cached[0] is source and cached[1] is target:
            return cached[2]
        if self._lattice_prefilter(source) != self._lattice_prefilter(target):
            self._isometry_cache[cache_key] = (source, target, None)
            return None
        reduced_source, source_reduction = self._simple_indefinite_reduction_data(source)
        reduced_target, target_reduction = self._simple_indefinite_reduction_data(target)
        if reduced_source is reduced_target:
            row_action = (source_reduction.inverse() * target_reduction).change_ring(SageZZ)
            source_labels = tuple(source.module_generating_set())
            target_labels = tuple(target.module_generating_set())

            def image(label):
                row = source_labels.index(label)
                return target.linear_combination(
                    {target_label: target.base_ring()(int(row_action[row, column])) for column, target_label in enumerate(target_labels) if row_action[row, column]}
                )

            witness = source.Isom(target)(image)
            self._isometry_cache[cache_key] = (source, target, witness)
            return witness
        reduced_source_to_source = None
        reduced_target_to_target = None
        if reduced_source is not source:
            _, reduced_source_to_source = self._simple_indefinite_reduction(source)
        if reduced_target is not target:
            _, reduced_target_to_target = self._simple_indefinite_reduction(target)
        if reduced_source is not source or reduced_target is not target:
            reduced_witness = self.isometry(reduced_source, reduced_target)
            if reduced_witness is None:
                self._isometry_cache[cache_key] = (source, target, None)
                return None
            row_action = reduced_witness.parent()._row_action_matrix(reduced_witness)
            if reduced_source_to_source is not None:
                source_reduction_action = reduced_source_to_source.parent()._row_action_matrix(reduced_source_to_source)
                row_action = source_reduction_action.inverse() * row_action
            if reduced_target_to_target is not None:
                target_reduction_action = reduced_target_to_target.parent()._row_action_matrix(reduced_target_to_target)
                row_action = row_action * target_reduction_action
            witness = source.Isom(target)._isometry_from_column_matrix(row_action.transpose())
            self._isometry_cache[cache_key] = (source, target, witness)
            return witness
        match source.gram_tensor() == target.gram_tensor():
            case True:
                witness = source.Isom(target)(tuple(target.module_generators()))
                self._isometry_cache[cache_key] = (source, target, witness)
                return witness
            case False:
                pass

        source_profile = self.attack_profile(source)
        target_profile = self.attack_profile(target)
        if source_profile.sign != target_profile.sign:
            self._isometry_cache[cache_key] = (source, target, None)
            return None
        signed_source = source_profile.signed_view
        signed_target = target_profile.signed_view

        match source_profile.positive_index:
            case 0 | 1:
                witness = super().isometry(signed_source, signed_target)
            case _:
                target_model = self._orbit_cover_model(signed_target)
                basis = tuple(signed_source.module_generators())
                candidates = list(basis)
                for left_position, left in enumerate(basis):
                    for right in basis[left_position + 1 :]:
                        candidates.append(left + right)
                        candidates.append(left - right)
                positive = tuple(vector for vector in candidates if vector.q() > signed_source.base_ring().zero())
                match positive:
                    case ():
                        raise ValueError("the higher-Witt source lattice has no positive vector in its framing span")
                    case _:
                        pass
                labels = tuple(signed_source.module_generating_set())
                source_vector = min(
                    positive,
                    key=lambda vector: (
                        abs(int(vector.q())),
                        tuple(int(vector.to_vector()(label)) for label in labels),
                    ),
                )
                witness = None
                for target_vector in target_model.covering_representatives(
                    source_vector.q(),
                    primitive=source_vector.is_primitive(),
                ):
                    witness = self.vector_transporter(source_vector, target_vector)
                    if witness is not None:
                        break

        if witness is None:
            self._isometry_cache[cache_key] = (source, target, None)
            return None
        if witness.domain() is source and witness.codomain() is target:
            self._isometry_cache[cache_key] = (source, target, witness)
            return witness
        witness = self._retarget_isometry(witness, source, target)
        self._isometry_cache[cache_key] = (source, target, witness)
        return witness

    def _integral_stabilizer_from_rational_lifts(
        self,
        vector: Lattices.ElementMethods,
        lifts: tuple[LatticeIsometryMethods, ...],
    ) -> GeneratedSubgroup:
        ambient = vector.parent()
        fraction_map = ambient.base_ring().fraction_field_map()
        rational_ambient = ambient.base_change(fraction_map)
        rational_lifts = tuple(self._retarget_isometry(lift, rational_ambient, rational_ambient) for lift in lifts)
        rational_group = RationalMatrixGroup(rational_ambient, rational_lifts)
        try:
            integral_generators = tuple(self._retarget_isometry(generator, ambient, ambient) for generator in rational_lifts)
        except TypeError, ValueError:
            integral_generators = ()
        if len(integral_generators) == len(rational_lifts):
            if any(generator(vector) != vector for generator in integral_generators):
                raise ArithmeticError("an integral vector-stabilizer generator moves the selected vector")
            return GeneratedSubgroup(RationalMatrixGroup(ambient, integral_generators), integral_generators)

        restriction = Modules(QQ).restriction_of_scalars(ZZ.Mor(QQ)(lambda element: QQ(element)))
        integral_structure_space = restriction(rational_ambient)
        inclusion = _rational_lattice_with_integral_structure(
            ambient,
            rational_ambient,
            integral_structure_space,
            _identity_rows(int(ambient.module_rank())),
        )
        rational_stabilizer = IntegralStructureAction(rational_group, inclusion).lattice_stabilizer()
        generators = tuple(self._retarget_isometry(generator, ambient, ambient) for generator in rational_stabilizer.generators())
        if any(generator(vector) != vector for generator in generators):
            raise ArithmeticError("an integral vector-stabilizer generator moves the selected vector")
        return GeneratedSubgroup(RationalMatrixGroup(ambient, generators), generators)

    def _rational_stabilizer_lifts(
        self,
        section,
        reduced_group: GeneratedSubgroup,
    ) -> tuple[LatticeIsometryMethods, ...]:
        lifts = tuple(section.rational_lift(generator, target=section).rational_isometry() for generator in reduced_group.generators())
        match section:
            case IsotropicVectorSection():
                fraction_map = section.vector.parent().base_ring().fraction_field_map()
                match section.vector.parent().is_even():
                    case True:
                        line = section.vector.parent().primitive_isotropic_subobject(section.vector)
                        unipotent = line.unipotent_group_generators()
                        return lifts + tuple(unipotent[label].base_change(fraction_map) for label in unipotent.index_set())
                    case False:
                        kernel = section.reduction.pointwise_perpendicular_kernel()
                        return lifts + tuple(generator.base_change(fraction_map) for generator in kernel.gens())
            case NonIsotropicVectorSection():
                return lifts

    def _split_unimodular_isotropic_transporter(
        self,
        source_vector: Lattices.ElementMethods,
        target_vector: Lattices.ElementMethods,
    ) -> LatticeIsometryMethods | None:
        r"""Use the integral split-U extension when both isotropic vectors have divisibility one."""
        source = source_vector.parent()
        target = target_vector.parent()
        match (
            source_vector.q() == source.base_ring().zero(),
            target_vector.q() == target.base_ring().zero(),
            source.is_even(),
            target.is_even(),
            source_vector.is_primitive(),
            target_vector.is_primitive(),
            source_vector.div() == source.base_ring().one(),
            target_vector.div() == target.base_ring().one(),
        ):
            case (True, True, True, True, True, True, True, True):
                pass
            case _:
                return None

        source_complement, source_basis = self._split_unimodular_isotropic_data(source_vector)
        target_complement, target_basis = self._split_unimodular_isotropic_data(target_vector)
        complement_isometry = self.isometry(source_complement, target_complement)
        match complement_isometry:
            case None:
                return None
            case _:
                pass

        match source is target:
            case True:
                target_parent = source.Aut()
            case False:
                target_parent = source.Isom(target)
        complement_row_action = complement_isometry.parent()._row_action_matrix(complement_isometry).change_ring(SageZZ)
        rank = source_basis.nrows()
        split_action = identity_matrix(SageZZ, rank)
        for row in range(complement_row_action.nrows()):
            for column in range(complement_row_action.ncols()):
                split_action[row + 2, column + 2] = complement_row_action[
                    row,
                    column,
                ]
        ambient_row_action = source_basis.inverse() * split_action * target_basis
        witness = target_parent._isometry_from_column_matrix(ambient_row_action.transpose())
        match witness(source_vector) == target_vector:
            case True:
                return witness
            case False:
                raise ArithmeticError("the split-U isotropic transporter does not carry the selected vector")

    def _split_unimodular_isotropic_data(
        self,
        vector: Lattices.ElementMethods,
    ) -> tuple[Lattices.ParentMethods, Matrix_integer_dense]:
        cached = self._split_isotropic_data_cache.get(id(vector))
        if cached is not None and cached[0] is vector:
            return cached[1], cached[2]

        ambient = vector.parent()
        ring = ambient.base_ring()
        labels = tuple(ambient.module_generating_set())
        gram = _engine_component_matrix(ambient.gram_tensor()).change_ring(SageZZ)
        coordinate_function = vector.to_vector()
        isotropic_coordinates = sage_vector(
            SageZZ,
            tuple(SageZZ(int(coordinate_function(label))) for label in labels),
        )
        pairings = isotropic_coordinates * gram
        gcd_value = SageZZ.zero()
        coefficients = []
        for pairing in pairings:
            new_gcd, old_coefficient, new_coefficient = gcd_value.xgcd(SageZZ(pairing))
            coefficients = [old_coefficient * coefficient for coefficient in coefficients]
            coefficients.append(new_coefficient)
            gcd_value = new_gcd
        if gcd_value == -SageZZ.one():
            coefficients = [-coefficient for coefficient in coefficients]
            gcd_value = -gcd_value
        if gcd_value != SageZZ.one():
            raise ArithmeticError("a divisibility-one isotropic vector has no Bezout partner")
        bezout_coordinates = sage_vector(SageZZ, coefficients)
        if (isotropic_coordinates * gram * bezout_coordinates.column())[0] != 1:
            raise ArithmeticError("the computed Bezout partner does not pair to one")
        bezout_square = (bezout_coordinates * gram * bezout_coordinates.column())[0]
        if bezout_square % 2:
            raise ArithmeticError("an even-lattice Bezout partner has odd square")
        partner_coordinates = bezout_coordinates - SageZZ(bezout_square // 2) * isotropic_coordinates
        match (
            (partner_coordinates * gram * partner_coordinates.column())[0] == 0,
            (isotropic_coordinates * gram * partner_coordinates.column())[0] == 1,
        ):
            case (True, True):
                pass
            case _:
                raise ArithmeticError("the divisibility-one Witt partner does not span a hyperbolic plane")
        perpendicular_pairings = matrix(
            SageZZ,
            [
                isotropic_coordinates * gram,
                partner_coordinates * gram,
            ],
        )
        complement_rows = perpendicular_pairings.right_kernel_matrix()
        complement_gram = complement_rows * gram * complement_rows.transpose()
        complement = Lattices(ring)([[ring(int(complement_gram[row, column])) for column in range(complement_gram.ncols())] for row in range(complement_gram.nrows())])
        split_basis = matrix(
            SageZZ,
            [
                tuple(isotropic_coordinates),
                tuple(partner_coordinates),
                *(tuple(row) for row in complement_rows.rows()),
            ],
        )
        if abs(split_basis.det()) != 1:
            raise ArithmeticError("a divisibility-one hyperbolic plane does not split the ambient lattice unimodularly")
        self._split_isotropic_data_cache[id(vector)] = (
            vector,
            complement,
            split_basis,
        )
        return complement, split_basis

    def vector_stabilizer(self, vector: Lattices.ElementMethods) -> GeneratedSubgroup:
        section = self._section(vector)
        reduced_group = self.orthogonal_group(section.reduced_object())
        rational_lifts = self._rational_stabilizer_lifts(section, reduced_group)
        return self._integral_stabilizer_from_rational_lifts(vector, rational_lifts)

    def vector_transporter(
        self,
        source_vector: Lattices.ElementMethods,
        target_vector: Lattices.ElementMethods,
    ) -> LatticeIsometryMethods | None:
        match (
            source_vector.q() == target_vector.q(),
            vector_content(source_vector) == vector_content(target_vector),
            source_vector.div() == target_vector.div(),
        ):
            case (True, True, True):
                pass
            case _:
                return None
        if source_vector.parent() is target_vector.parent() and source_vector == target_vector:
            return source_vector.parent().O().identity()

        split_witness = self._split_unimodular_isotropic_transporter(
            source_vector,
            target_vector,
        )
        match split_witness:
            case None:
                pass
            case _:
                return split_witness

        source_section = self._section(source_vector)
        target_section = self._section(target_vector)
        if isinstance(source_section, NonIsotropicVectorSection) != isinstance(target_section, NonIsotropicVectorSection):
            return None
        reduced_isometry = self.isometry(source_section.reduced_object(), target_section.reduced_object())
        if reduced_isometry is None:
            return None

        source = source_vector.parent()
        target = target_vector.parent()
        match (source_section, target_section):
            case (
                NonIsotropicVectorSection(),
                NonIsotropicVectorSection(),
            ):
                integral_witness = source_section.integral_lift(
                    reduced_isometry,
                    target=target_section,
                )
                match integral_witness:
                    case None:
                        pass
                    case _:
                        return integral_witness
            case _:
                pass

        base_torsor = source_section.rational_lift(
            reduced_isometry,
            target=target_section,
        )
        match base_torsor.integral_parameters(source, target):
            case None:
                pass
            case _:
                witness = base_torsor.one_integral_extension()
                match witness(source_vector) == target_vector:
                    case True:
                        return witness
                    case False:
                        raise ArithmeticError("an integral vector transporter does not carry the selected source vector to the target")

        base_transport = base_torsor.rational_isometry()
        target_reduced_group = self.orthogonal_group(target_section.reduced_object())
        target_rational = base_transport.codomain()
        rational_lifts = tuple(
            self._retarget_isometry(
                lift,
                target_rational,
                target_rational,
            )
            for lift in self._rational_stabilizer_lifts(
                target_section,
                target_reduced_group,
            )
        )
        rational_group = RationalMatrixGroup(target_rational, rational_lifts)
        restriction = Modules(QQ).restriction_of_scalars(ZZ.Mor(QQ)(lambda element: QQ(element)))
        integral_structure_space = restriction(target_rational)
        target_inclusion = _rational_lattice_with_integral_structure(
            target,
            target_rational,
            integral_structure_space,
            _identity_rows(int(target.module_rank())),
        )
        base_action = base_transport.parent()._row_action_matrix(base_transport)
        source_image = _rational_lattice_with_integral_structure(
            source,
            target_rational,
            integral_structure_space,
            tuple(tuple(base_action[row, column] for column in range(base_action.ncols())) for row in range(base_action.nrows())),
        )
        correction = IntegralStructureAction(
            rational_group,
            target_inclusion,
        ).transporter(
            source_image,
            target_inclusion,
        )
        match correction:
            case None:
                return None
            case _:
                pass
        witness = self._retarget_isometry(
            correction * base_transport,
            source,
            target,
        )
        match witness(source_vector) == target_vector:
            case True:
                return witness
            case False:
                raise ArithmeticError("an integral vector transporter correction does not carry the selected source vector to the target")

    def vector_orbit_representatives(
        self,
        lattice: Lattices.ParentMethods,
        square,
    ) -> tuple[Lattices.ElementMethods, ...]:
        r"""Return exact higher-Witt O(L)-orbit representatives of a square."""
        reduced, reduced_to_lattice = self._simple_indefinite_reduction(lattice)
        if reduced is not lattice:
            assert reduced_to_lattice is not None
            return tuple(
                reduced_to_lattice(representative)
                for representative in self.vector_orbit_representatives(
                    reduced,
                    square,
                )
            )
        cache_key = (lattice, int(square))
        cached = self._vector_orbit_cache.get(cache_key)
        match cached:
            case tuple() as representatives:
                return representatives
            case _ if cached is not None:
                return cached
            case None:
                pass
        profile = self.attack_profile(lattice)
        if profile.positive_index < 2:
            raise NotImplementedError("recursive vector-orbit representatives currently require higher Witt index")
        signed_lattice = profile.signed_view
        signed_square = signed_lattice.base_ring()(profile.sign * int(square))
        model = self._orbit_cover_model(signed_lattice)
        cover = model.covering_representatives(
            signed_square,
            primitive=signed_square == signed_lattice.base_ring().zero(),
        )
        labels = tuple(signed_lattice.module_generating_set())
        cover_coordinates = tuple(tuple(int(candidate.to_vector()(label)) for label in labels) for candidate in cover)
        coordinate_positions = {coordinates: position for position, coordinates in enumerate(cover_coordinates)}
        classes = SageDisjointSet(len(cover_coordinates))
        known_generator_rows = set()
        generator_matrices = []
        for matrix_ in model._subgroup_column_matrices():
            integral_matrix = matrix_.change_ring(SageZZ)
            rows = tuple(tuple(int(integral_matrix[row, column]) for column in range(integral_matrix.ncols())) for row in range(integral_matrix.nrows()))
            if rows not in known_generator_rows:
                known_generator_rows.add(rows)
                generator_matrices.append(integral_matrix)

        if cover_coordinates and generator_matrices:
            rank = len(labels)
            stacked_actions = matrix(
                SageZZ,
                tuple(row for generator_matrix in generator_matrices for row in generator_matrix.rows()),
            )
            candidate_columns = matrix(
                SageZZ,
                rank,
                len(cover_coordinates),
                lambda row, column: cover_coordinates[column][row],
            )
            image_columns = stacked_actions * candidate_columns
            for generator_position in range(len(generator_matrices)):
                row_offset = generator_position * rank
                for candidate_position in range(len(cover_coordinates)):
                    image = tuple(int(image_columns[row_offset + row, candidate_position]) for row in range(rank))
                    image_position = coordinate_positions.get(image)
                    if image_position is not None:
                        classes.union(candidate_position, image_position)

        simplified_coordinates = tuple(
            sorted(
                min(
                    (cover_coordinates[position] for position in component),
                )
                for component in classes.root_to_elements_dict().values()
            )
        )
        cover = tuple(
            signed_lattice.linear_combination(
                {
                    label: signed_lattice.base_ring()(coordinate)
                    for label, coordinate in zip(
                        labels,
                        coordinates,
                        strict=True,
                    )
                    if coordinate
                }
            )
            for coordinates in simplified_coordinates
        )
        buckets: dict[VectorPrefilter, list[Lattices.ElementMethods]] = {}
        for candidate in cover:
            bucket = buckets.setdefault(
                VectorPrefilter.from_vector(
                    candidate,
                    include_orthogonal_reduction=False,
                ),
                [],
            )
            for representative in bucket:
                if self.vector_transporter(representative, candidate) is not None:
                    break
            else:
                bucket.append(candidate)

        representatives = tuple(representative for bucket in buckets.values() for representative in bucket)
        if signed_lattice is lattice:
            self._vector_orbit_cache[cache_key] = representatives
            return representatives

        lattice_labels = tuple(lattice.module_generating_set())
        signed_labels = tuple(signed_lattice.module_generating_set())
        transported = tuple(
            lattice.linear_combination(
                {
                    lattice_label: lattice.base_ring()(representative.to_vector()(signed_label))
                    for lattice_label, signed_label in zip(
                        lattice_labels,
                        signed_labels,
                        strict=True,
                    )
                    if representative.to_vector()(signed_label)
                }
            )
            for representative in representatives
        )
        self._vector_orbit_cache[cache_key] = transported
        return transported


@unfinished(18)
def orthogonal_group_generators(homset) -> tuple[LatticeIsometryMethods, ...]:
    lattice = homset.domain()
    group = IndefiniteOrthogonalAlgorithm().orthogonal_group(lattice)
    return tuple(group.generators())


@unfinished(18)
def isometry(source: Lattices.ParentMethods, target: Lattices.ParentMethods):
    return IndefiniteOrthogonalAlgorithm().isometry(source, target)


@unfinished(18)
def vector_equivalence_witness(homset, left, right):
    if left.parent() is not homset.domain() or right.parent() is not homset.domain():
        raise ValueError("vector equivalence requires vectors in the homset lattice")
    return IndefiniteOrthogonalAlgorithm().vector_transporter(left, right)


@unfinished(18)
def vector_stabilizer_generators(homset, element) -> tuple[LatticeIsometryMethods, ...]:
    if element.parent() is not homset.domain():
        raise ValueError("vector stabilizer requires an element of the homset lattice")
    subgroup = IndefiniteOrthogonalAlgorithm().vector_stabilizer(element)
    return tuple(subgroup.generators())


__all__ = [
    "IndefiniteOrthogonalAlgorithm",
    "isometry",
    "orthogonal_group_generators",
    "vector_equivalence_witness",
    "vector_stabilizer_generators",
]
