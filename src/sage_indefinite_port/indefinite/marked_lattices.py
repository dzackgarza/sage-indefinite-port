r"""Isomorphisms and automorphisms of marked indefinite lattices.

A marked lattice is a diagram of lattice morphisms ending at a lattice L
(references/port-design.md, "Marked lattices"):

- a plain lattice L, an object of ``Lattices``;
- a pointed lattice ``<n> -> L``, an object of the coslice ``<n>/Lattices``;
- a primitive totally isotropic sublattice ``I -> L`` of any attainable rank, an object of ``Ar(Lattices)``;
- a partial or complete isotropic flag ``I_1 -> ... -> I_t -> L``, an object of ``Fun([t], Lattices)``.

The research preamble's core homsets of these categories call these three entry points, with
the category, when L is indefinite. An isomorphism of marked lattices is returned as its
isometry of the underlying lattices: the markings are monomorphisms, so the isometry determines
the other components.
"""

from itertools import pairwise
from math import gcd

from sage_extensions.all import ZZ, Cat, Lattices
from sage_extensions.categories.abstract_categories.arrow_categories import (
    CosliceCategory,
)
from sage_extensions.categories.abstract_categories.products import (
    FiniteOrdinalCategory,
    _FiniteSequenceFunctorCategory,
)
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage_extensions.categories.sets.finite_ordered_sets import (
    finite_ordered_set,
)
from sage.categories.category import Category
from sage.categories.morphism import Morphism
from sage.structure.parent import Parent

from sage_indefinite_port.indefinite.indefinite_reduction import (
    indefinite_reduced_presentation,
)
from sage_indefinite_port.indefinite.isotropic_induction import (
    ARROWS,
    Flag,
    IsotropicInduction,
    UnsupportedRankError,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    MarkedCellOrbitAlgorithm,
    _transfer,
    binary_split_primitive_vector_classes,
    perfect_complex,
)
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _gram_in_label_order
from sage_indefinite_port.invariants import AttackProfile

LATTICES = Lattices(ZZ)
FLAGS = Cat().Mor(FiniteOrdinalCategory(3), LATTICES)
"""Fun([2], Lattices): flags I_1 -> I_2 -> L."""

RECURSION = EichlerRecursion()
INDUCTION = IsotropicInduction(RECURSION)


def _marked_vector(marked: Parent) -> Lattices.ElementMethods:
    arrow = marked.arrow()
    return arrow(next(iter(arrow.domain().module_generators())))


def _isotropic_sublattice(arrow: Morphism) -> Parent:
    r"""The primitive totally isotropic sublattice that a marking of Ar(Lattices) embeds, as a subobject of L."""
    lattice = arrow.codomain()
    vectors = tuple(arrow(generator) for generator in arrow.domain().module_generators())
    if any(lattice.b(left, right) != 0 for left in vectors for right in vectors):
        raise ValueError(f"{arrow} does not mark a totally isotropic sublattice: Isom and Aut of other sublattice markings are the preamble's (research#399)")
    if not arrow.is_primitive():
        raise ValueError(f"{arrow} does not mark a primitive sublattice")
    return lattice.primitive_isotropic_subobject(*vectors)


def _is_totally_isotropic(arrow: Morphism) -> bool:
    lattice = arrow.codomain()
    vectors = tuple(arrow(generator) for generator in arrow.domain().module_generators())
    return all(lattice.b(left, right) == 0 for left in vectors for right in vectors)


def _is_nondegenerate_primitive(arrow: Morphism) -> bool:
    return arrow.is_primitive() and arrow.domain().is_nondegenerate() is True


def _saturated_arrow(arrow: Morphism) -> Morphism:
    r"""The inclusion of the saturation of the image of ``arrow``."""
    return arrow.saturation().inclusion()


def _saturation_exponent(arrow: Morphism, inclusion: Morphism) -> int:
    return int(arrow.factor_through(inclusion).cokernel().exponent())


def _nonprimitive_setwise_generators(arrow: Morphism, generators: tuple[LatticeIsometryMethods, ...]) -> tuple[LatticeIsometryMethods, ...]:
    r"""Schreier generators of the stabilizer of ``im(arrow)`` in the setwise stabilizer of its saturation.

    The stabilizer acts on ``S'/nS'``. Reidemeister--Schreier on the finite orbit of the span of ``S``.
    """
    lattice = arrow.codomain()
    inclusion = _saturated_arrow(arrow)
    modulus = _saturation_exponent(arrow, inclusion)
    base = [inclusion.lift(arrow(point)) for point in arrow.domain().module_generators()]

    def span_key(points):
        from sage.matrix.constructor import matrix
        from sage.rings.finite_rings.integer_mod_ring import Zmod

        rows = [[int(point.to_vector()(label)) % modulus for label in inclusion.domain().module_generating_set()] for point in points]
        return tuple(tuple(row) for row in matrix(Zmod(modulus), rows).echelon_form().rows() if any(row))

    def act(generator, points):
        return [inclusion.lift(generator(inclusion(point))) for point in points]

    transversal = {span_key(base): (lattice.O().one(), base)}
    frontier = [span_key(base)]
    schreier = []
    while frontier:
        key = frontier.pop()
        element, points = transversal[key]
        for generator in generators:
            moved = act(generator, points)
            moved_key = span_key(moved)
            if moved_key not in transversal:
                transversal[moved_key] = (generator * element, moved)
                frontier.append(moved_key)
            else:
                schreier.append(~transversal[moved_key][0] * generator * element)
    return tuple(element for element in schreier if element != lattice.O().one())


def _product(values: tuple[int, ...]) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def nonprimitive_sublattice_classes(saturated: Morphism, divisors: tuple[int, ...]) -> tuple[Morphism, ...]:
    r"""Inclusions of the classes of ``J <= S`` with elementary divisors ``divisors``, under ``Stab(S)``.

    #48: ``S`` is primitive and each ``J`` contains ``mS``, so the classes are
    the orbits of the actual stabilizer image on the finite set of such ``J``.
    """
    from itertools import product

    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    ambient = saturated.codomain()
    module = saturated.domain()
    labels = tuple(module.module_generating_set())
    rank = len(labels)
    if len(divisors) != rank:
        raise ValueError(f"{divisors} must have one elementary divisor for each of the {rank} generators of {module}")
    wanted = sorted(int(divisor) for divisor in divisors)
    index = _product(tuple(wanted))
    candidates: list[tuple[tuple[int, ...], ...]] = []
    for diagonal in product(*[[d for d in range(1, index + 1) if index % d == 0]] * rank):
        if _product(diagonal) != index:
            continue
        free = [(i, j) for i in range(rank) for j in range(i + 1, rank)]
        for values in product(*[range(diagonal[j]) for _i, j in free]):
            rows = [[0] * rank for _ in range(rank)]
            for i, value in enumerate(diagonal):
                rows[i][i] = value
            for (i, j), value in zip(free, values, strict=True):
                rows[i][j] = value
            if sorted(int(d) for d in matrix(SageZZ, rows).elementary_divisors()) == wanted:
                candidates.append(tuple(tuple(row) for row in rows))
    generators = tuple(automorphism_group_generators(ARROWS, ARROWS.object(saturated)))

    def act(generator, rows):
        moved = []
        for row in rows:
            point = saturated(module.linear_combination({label: module.base_ring()(entry) for label, entry in zip(labels, row, strict=True) if entry}))
            coordinates = saturated.lift(generator(point)).to_vector()
            moved.append([int(coordinates(label)) for label in labels])
        return tuple(tuple(int(entry) for entry in row) for row in matrix(SageZZ, moved).hermite_form().rows())

    remaining = set(candidates)
    classes = []
    while remaining:
        start = min(remaining)
        orbit = {start}
        frontier = [start]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                moved = act(generator, current)
                if moved not in orbit:
                    orbit.add(moved)
                    frontier.append(moved)
        remaining -= orbit
        classes.append(start)
    return tuple(
        ambient.subobject_on(tuple(
            saturated(module.linear_combination({label: module.base_ring()(entry) for label, entry in zip(labels, row, strict=True) if entry}))
            for row in rows
        )).inclusion()
        for rows in classes
    )


def _flag(marked: Parent) -> Flag:
    first, second = _flag_arrows(marked)
    plane = _isotropic_sublattice(second)
    return second(first(next(iter(first.domain().module_generators())))), plane


def _flag_arrows(marked: Parent) -> tuple[Morphism, ...]:
    r"""Read adjacent arrows from the retained finite-ordinal diagram."""
    diagram = marked.functor()
    shape = diagram.domain()
    return tuple(
        diagram(shape.Mor(source, target).unique())
        for source, target in pairwise(shape.objects())
    )


def _flag_action(marked: Parent) -> tuple[tuple[int, ...], tuple[Lattices.ElementMethods, ...], tuple[tuple[Lattices.ElementMethods, ...], ...]]:
    r"""Read a primitive isotropic chain and its radical flag action.

    CombinedAlgorithms.h:655-679 f_get_list_spaces and :682-717
    GetAutomorphismOfFlag, expressed as subobjects of the last isotropic term.
    """
    arrows = _flag_arrows(marked)
    space = _isotropic_sublattice(arrows[-1])
    terms = []
    composite = None
    for arrow in reversed(arrows[:-1]):
        composite = arrow if composite is None else composite * arrow
        terms.append(_isotropic_sublattice(composite))
    terms.reverse()
    dimensions = tuple(int(term.inclusion().domain().rank()) for term in terms) + (int(space.inclusion().domain().rank()),)
    if any(left >= right for left, right in zip(dimensions, dimensions[1:])):
        raise ValueError("an isotropic flag must have strictly increasing ranks")
    basis, images = INDUCTION.flag_action_data(space, tuple(terms))
    return dimensions, basis, images


def isometry(category: Category, source: Parent, target: Parent) -> LatticeIsometryMethods | None:
    r"""The isometry of an isomorphism of marked lattices ``source -> target`` in ``category``, or ``None`` when there is none."""
    match category:
        case CosliceCategory():
            rank = int(category.base_object().rank())
            if rank == 0:
                return RECURSION.isometry(source.arrow().codomain(), target.arrow().codomain())
            if rank == 1:
                source_vector, target_vector = _marked_vector(source), _marked_vector(target)
                if source_vector.content() != target_vector.content():
                    return None
                # An isometry carries v=c*v0 to w=c*w0 exactly when it carries v0 to w0 (#48).
                return RECURSION.pointed_isometry(source_vector.primitive_part(), target_vector.primitive_part())
            if not _is_totally_isotropic(source.arrow()) and _is_nondegenerate_primitive(source.arrow()):
                from sage_extensions.categories.lattice_embeddings import _nondegenerate_marking_transporter

                return _nondegenerate_marking_transporter(source.arrow(), target.arrow(), pointwise=True)
            _isotropic_sublattice(source.arrow())
            _isotropic_sublattice(target.arrow())
            labels = tuple(category.base_object().module_generating_set())
            basis = tuple(category.base_object().module_generator(label) for label in labels)
            # IndefiniteForms.g:1660-1694 with trivial action on the
            # source: every supplied generator is fixed by the marking.
            return RECURSION.isotropic_subspace_isometry(
                tuple(source.arrow()(point) for point in basis),
                tuple(target.arrow()(point) for point in basis),
                (),
            )
        case _ if category is LATTICES:
            return RECURSION.isometry(source, target)
        case _ if category is ARROWS and not source.arrow().is_primitive():
            inclusion = _saturated_arrow(source.arrow())
            target_inclusion = _saturated_arrow(target.arrow())
            transporter = isometry(ARROWS, ARROWS.object(inclusion), ARROWS.object(target_inclusion))
            if transporter is None:
                return None
            return _nonprimitive_transporter(source.arrow(), target.arrow(), transporter, inclusion)
        case _ if category is ARROWS:
            source_rank = int(source.arrow().domain().rank())
            target_rank = int(target.arrow().domain().rank())
            if source_rank == 0 or target_rank == 0:
                return RECURSION.isometry(source.arrow().codomain(), target.arrow().codomain()) if source_rank == target_rank else None
            if not _is_totally_isotropic(source.arrow()) and _is_nondegenerate_primitive(source.arrow()):
                from sage_extensions.categories.lattice_embeddings import _nondegenerate_marking_transporter

                return _nondegenerate_marking_transporter(source.arrow(), target.arrow(), pointwise=False)
            source_sublattice = _isotropic_sublattice(source.arrow())
            target_sublattice = _isotropic_sublattice(target.arrow())
            match int(source_sublattice.inclusion().domain().rank()), int(target_sublattice.inclusion().domain().rank()):
                case 1, 1:
                    return INDUCTION.line_isometry(INDUCTION.plane_vectors(source_sublattice)[0], INDUCTION.plane_vectors(target_sublattice)[0])
                case source_rank, target_rank if source_rank == target_rank:
                    # CombinedAlgorithms.h:1401-1404: unequal isotropic k-plane invariants
                    # (INDEF_FORM_Invariant_IsotropicKplane, :2424-2428) decide inequivalence.
                    if INDUCTION.flag_invariant(source.arrow().codomain(), INDUCTION.plane_vectors(source_sublattice), (source_rank,)) != INDUCTION.flag_invariant(
                        target.arrow().codomain(), INDUCTION.plane_vectors(target_sublattice), (target_rank,)
                    ):
                        return None
                    return INDUCTION.plane_isometry(source_sublattice, target_sublattice)
                case _:
                    return None
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            source_dimensions, source_basis, _source_images = _flag_action(source)
            target_dimensions, target_basis, target_images = _flag_action(target)
            if source_dimensions != target_dimensions:
                return None
            # CombinedAlgorithms.h:1401-1404: unequal isotropic flag invariants
            # (INDEF_FORM_Invariant_IsotropicKflag, :2429-2433) decide inequivalence.
            if INDUCTION.flag_invariant(source_basis[0].parent(), source_basis, source_dimensions) != INDUCTION.flag_invariant(
                target_basis[0].parent(), target_basis, target_dimensions
            ):
                return None
            return RECURSION.isotropic_subspace_isometry(source_basis, target_basis, target_images)
    raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")


def _nonprimitive_transporter(source: Morphism, target: Morphism, transporter: LatticeIsometryMethods, inclusion: Morphism) -> LatticeIsometryMethods | None:
    r"""Correct ``transporter`` by the stabilizer of the target saturation, or return ``None``.

    BFS in the finite orbit of the moved finite-index sublattice in ``S'/nS'``.
    """
    from sage.matrix.constructor import matrix
    from sage.rings.finite_rings.integer_mod_ring import Zmod

    lattice = source.codomain()
    target_inclusion = _saturated_arrow(target)
    modulus = _saturation_exponent(target, target_inclusion)
    generators = tuple(automorphism_group_generators(ARROWS, ARROWS.object(target_inclusion)))

    def span_key(points):
        rows = [[int(target_inclusion.lift(point).to_vector()(label)) % modulus for label in target_inclusion.domain().module_generating_set()] for point in points]
        return tuple(tuple(row) for row in matrix(Zmod(modulus), rows).echelon_form().rows() if any(row))

    target_key = span_key([target(point) for point in target.domain().module_generators()])
    start_points = [transporter(source(point)) for point in source.domain().module_generators()]
    seen = {span_key(start_points): (transporter, start_points)}
    frontier = [span_key(start_points)]
    while frontier and target_key not in seen:
        element, points = seen[frontier.pop()]
        for generator in generators:
            moved = [generator(point) for point in points]
            key = span_key(moved)
            if key not in seen:
                seen[key] = (generator * element, moved)
                frontier.append(key)
    return seen[target_key][0] if target_key in seen else None


def automorphism_group_generators(category: Category, marked: Parent) -> Parent:
    r"""Generators of the automorphism group of a marked lattice in ``category``, as a finite ordered set of isometries of L."""
    if isinstance(category, CosliceCategory) or category is ARROWS:
        arrow = marked.arrow()
        ambient = arrow.codomain()
        source = arrow.domain()
        if int(source.rank()) > 0 and not ambient.is_nondegenerate():
            from sage_indefinite_port.indefinite.degenerate_markings import degenerate_marking_generators

            frame = tuple(ambient.module_generators())
            labels = tuple(ambient.module_generating_set())
            vectors = tuple(arrow(point) for point in source.module_generators())
            matrices = degenerate_marking_generators(ambient, vectors, isinstance(category, CosliceCategory))
            return finite_ordered_set(
                tuple(
                    ambient.O()(
                        {
                            label: sum(
                                (ambient.scalar_multiple(int(c), g) for c, g in zip(m.column(j), frame, strict=True)),
                                ambient.zero(),
                            )
                            for j, label in enumerate(labels)
                        }
                    )
                    for m in matrices
                )
            )
    match category:
        case CosliceCategory():
            rank = int(category.base_object().rank())
            if rank == 0:
                generators = RECURSION.orthogonal_group_generators(marked.arrow().codomain())
            elif rank == 1:
                generators = RECURSION.pointed_automorphism_generators(_marked_vector(marked).primitive_part())
            elif not _is_totally_isotropic(marked.arrow()) and _is_nondegenerate_primitive(marked.arrow()):
                from sage_extensions.categories.lattice_embeddings import _nondegenerate_marking_stabilizer_generators

                generators = _nondegenerate_marking_stabilizer_generators(marked.arrow(), pointwise=True)
            else:
                _isotropic_sublattice(marked.arrow())
                # The same perpendicular extension as the setwise case,
                # with trivial radical action (IndefiniteForms.g:1710-1723).
                basis = tuple(marked.arrow()(point) for point in category.base_object().module_generators())
                generators = INDUCTION._stabilizer_generators(basis, ())
        case _ if category is LATTICES:
            generators = RECURSION.orthogonal_group_generators(marked)
        case _ if category is ARROWS and not marked.arrow().is_primitive():
            inclusion = _saturated_arrow(marked.arrow())
            generators = _nonprimitive_setwise_generators(marked.arrow(), tuple(automorphism_group_generators(ARROWS, ARROWS.object(inclusion))))
        case _ if category is ARROWS:
            if int(marked.arrow().domain().rank()) == 0:
                generators = RECURSION.orthogonal_group_generators(marked.arrow().codomain())
            elif not _is_totally_isotropic(marked.arrow()) and _is_nondegenerate_primitive(marked.arrow()):
                from sage_extensions.categories.lattice_embeddings import _nondegenerate_marking_stabilizer_generators

                generators = _nondegenerate_marking_stabilizer_generators(marked.arrow(), pointwise=False)
            else:
                sublattice = _isotropic_sublattice(marked.arrow())
                match int(sublattice.inclusion().domain().rank()):
                    case 1:
                        generators = INDUCTION.line_automorphism_generators(marked.arrow().codomain(), INDUCTION.plane_vectors(sublattice)[0])
                    case _:
                        generators = INDUCTION.plane_automorphism_generators(sublattice)
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            _dimensions, basis, images = _flag_action(marked)
            generators = INDUCTION._stabilizer_generators(basis, images)
        case _:
            raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")
    return finite_ordered_set(tuple(generators))


def _radical_residues(content: int, radical_rank: int) -> tuple[tuple[int, ...], ...]:
    r"""Representatives of ``{r in R/dR : gcd(content(r), d) = 1} / GL(R)``.

    Rank at least two: ``SL_t(Z) -> SL_t(Z/d)`` is onto and transitive on
    unimodular residues. Rank one: the units modulo ``+-1``.
    """
    if content == 1:
        return ((0,) * radical_rank,)
    if radical_rank >= 2:
        return ((1,) + (0,) * (radical_rank - 1),)
    return tuple((unit,) for unit in range(1, content // 2 + 1) if gcd(unit, content) == 1)


def _degenerate_pointed_classes(category: Category, lattice: Lattices.ParentMethods, square: int) -> tuple[Parent, ...]:
    r"""Primitive vector classes of nonzero square on ``L = C + R`` with ``R`` the radical (#41).

    ``O(L)`` has blocks ``[A B; 0 D]`` with ``A in O(C)``, ``D in GL(R)`` and
    ``B`` arbitrary. For ``v = (c, r)`` with ``c`` of content ``d``, the shears
    translate ``r`` by exactly ``dR``, and ``v`` is primitive when ``r`` is
    unimodular modulo ``d``. So classes are pairs of an ``O(C)``-class of ``c``
    and a ``GL(R)``-class of ``r`` modulo ``d``.
    """
    reduction = indefinite_reduced_presentation(lattice)
    reduced = reduction.reduced
    gram = _gram_in_label_order(reduced)
    size, rank = int(gram.nrows()), int(gram.rank())
    labels = tuple(reduced.module_generating_set())
    complement = lattice.lattice_category()([[int(gram[i, j]) for j in range(rank)] for i in range(rank)])
    complement_labels = tuple(complement.module_generating_set())
    ring = reduced.base_ring()
    vectors = []
    for point in complement.O().vector_orbit_representatives(square):
        coordinates = tuple(int(point.to_vector()(label)) for label in complement_labels)
        content = gcd(*coordinates)
        for residue in _radical_residues(content, size - rank):
            vectors.append(reduction.isometry(reduced.linear_combination({
                label: ring(value) for label, value in zip(labels, coordinates + residue, strict=True) if value
            })))
    base = category.base_object()
    return tuple(category.object(base.hom([vector_], codomain=lattice)) for vector_ in vectors)


class DegenerateNullClasses:
    r"""Primitive square-zero classes of ``L = C + R``, in blocks of the content ``d`` of the ``C``-part (#41).

    Block zero is the radical class ``(0, e_1)``. Block ``d >= 1`` has one class
    per primitive isotropic class ``c'`` of ``C`` and per residue of
    :func:`_radical_residues` modulo ``d``.
    """

    def __init__(self, lattice: Lattices.ParentMethods) -> None:
        self.reduction = indefinite_reduced_presentation(lattice)
        reduced = self.reduction.reduced
        gram = _gram_in_label_order(reduced)
        self.size, self.rank = int(gram.nrows()), int(gram.rank())
        self.labels = tuple(reduced.module_generating_set())
        self.ring = reduced.base_ring()
        self.primitive: tuple[tuple[int, ...], ...] = ()
        if self.rank:
            complement = lattice.lattice_category()([[int(gram[i, j]) for j in range(self.rank)] for i in range(self.rank)])
            names = tuple(complement.module_generating_set())
            self.primitive = tuple(
                tuple(int(point.to_vector()(name)) for name in names)
                for point in complement.O().vector_orbit_representatives(0, content=1)
            )

    def is_finite(self) -> bool:
        return not self.primitive

    def _point(self, coordinates: tuple[int, ...]) -> Lattices.ElementMethods:
        reduced = self.reduction.reduced
        return self.reduction.isometry(reduced.linear_combination({
            label: self.ring(value) for label, value in zip(self.labels, coordinates, strict=True) if value
        }))

    def _block(self, content: int) -> tuple[tuple[int, ...], ...]:
        return tuple(
            tuple(content * entry for entry in primitive) + residue
            for primitive in self.primitive
            for residue in _radical_residues(content, self.size - self.rank)
        )

    def coordinates(self):
        yield (0,) * self.rank + (1,) + (0,) * (self.size - self.rank - 1)
        if not self.primitive:
            return
        content = 1
        while True:
            yield from self._block(content)
            content += 1

    def __iter__(self):
        return map(self._point, self.coordinates())

    def index_of(self, vector: Lattices.ElementMethods) -> int | None:
        inverse = self.reduction.isometry.inverse()
        target = tuple(int(inverse(vector).to_vector()(label)) for label in self.labels)
        bound = gcd(*target[: self.rank])
        for index, candidate in enumerate(self.coordinates()):
            if candidate == target:
                return index
            if gcd(*candidate[: self.rank]) > bound > 0 or (bound == 0 and index > 0):
                return None
        return None


def _pointed_classes(category: Category, lattice: Lattices.ParentMethods, square: int) -> tuple[Parent, ...]:
    if not lattice.is_nondegenerate() and square != 0:
        return _degenerate_pointed_classes(category, lattice, square)
    profile = AttackProfile.from_lattice(lattice)
    if profile.positive_index == 1 and int(lattice.rank()) == 2 and tuple(map(int, lattice.signature_pair())) == (1, 1) and not lattice.is_isotropic():
        # An anisotropic binary plane has no isotropic cusp and therefore no
        # perfect-cell fundamental domain. Use the research owner's exact
        # nonsplit binary fixed-norm orbit set, whose finite transversal
        # includes imprimitive vectors; retain only primitive markings here.
        # An anisotropic plane has no nonzero vector of square zero.
        if square == 0:
            return ()
        orbit_set = lattice.binary_fixed_norm_representatives(square)
        section = orbit_set.section()
        vectors = tuple(
            vector_
            for orbit in orbit_set
            if (vector_ := section(orbit)).content() == 1
        )
        base = category.base_object()
        return tuple(category.object(base.hom([vector_], codomain=lattice)) for vector_ in vectors)
    if profile.positive_index == 1 and int(lattice.rank()) == 2 and tuple(map(int, lattice.signature_pair())) == (1, 1) and lattice.is_isotropic():
        # The pinned binary split dispatcher covers both timelike and
        # spacelike squares; perfect cells cover only the closed positive cone.
        vectors = binary_split_primitive_vector_classes(lattice, square)
    else:
        match profile.positive_index, square:
            case 0, _:
                raise ValueError(f"{lattice} is definite: its vector classes are the preamble's")
            case 1, 0:
                vectors = INDUCTION.isotropic_vector_classes(lattice)
            case 1, _ if profile.sign * square < 0:
                # These vectors lie outside the positive cone covered by
                # perfect cells.  Classify their rank-one primitive
                # embeddings by the generic Nikulin genus/gluing algorithm
                # at its research owner (pinned OSCAR embeddings.jl,
                # d135b70b, primitive_embeddings(::ZZGenus, ::ZZLat)).
                base = category.base_object()
                return tuple(
                    category.object(embedding)
                    for embedding in base.Emb(lattice).primitive_embedding_class_representatives(
                        classification="emb"
                    )
                )
            case 1, _:
                signed_square = profile.sign * square
                reduction = indefinite_reduced_presentation(profile.signed_view)
                vectors = tuple(
                    _transfer(lattice, reduction.isometry(orbit[0]))
                    for orbit in MarkedCellOrbitAlgorithm(perfect_complex(reduction.reduced, NONNEGATIVE_SQUARES), signed_square).global_orbits()
                )
            case _:
                vectors = tuple(vector_class.representative for vector_class in RECURSION.primitive_vector_classes(lattice, square))
    base = category.base_object()
    return tuple(category.object(base.hom([vector_], codomain=lattice)) for vector_ in vectors)


def primitive_marking_classes(category: Category, source: Parent | Morphism, lattice: Parent) -> Parent:
    r"""One marked lattice for each isomorphism class of primitive markings of ``lattice`` of shape ``source``.

    ``category`` is the category of the marking: the coslice under ``<n>`` for pointed
    lattices, ``Ar(Lattices)`` for sublattices, ``Fun([t], Lattices)`` for finite flags. ``source``
    is the shape: the lattice ``<n>``, a totally isotropic lattice of rank ``k`` (zero form),
    or the chain of such lattices for a flag.
    """
    match category:
        case CosliceCategory():
            if int(source.rank()) == 0:
                classes = (category.object(source.hom([], codomain=lattice)),)
            elif int(source.rank()) == 1:
                classes = _pointed_classes(category, lattice, int(source.q(next(iter(source.module_generators())))))
            else:
                source_basis = tuple(source.module_generators())
                if any(source.b(left, right) != 0 for left in source_basis for right in source_basis):
                    # #46: fixed-source classes of a nondegenerate source are genus and glue data.
                    return finite_ordered_set(tuple(
                        category.object(embedding)
                        for embedding in source.Emb(lattice).primitive_embedding_class_representatives(classification="emb")
                    ))
                classes = tuple(
                    category.object(source.hom(list(images), codomain=lattice))
                    for basis in INDUCTION.isotropic_flag_bases(lattice, (int(source.rank()),))
                    for images in INDUCTION.basis_marking_classes(lattice, basis)
                )
        case _ if category is ARROWS:
            rank = int(source.rank())
            if rank == 0:
                classes = (ARROWS.object(source.hom([], codomain=lattice)),)
            else:
                classes = tuple(ARROWS.object(lattice.primitive_isotropic_subobject(*basis).inclusion()) for basis in INDUCTION.isotropic_flag_bases(lattice, (rank,)))
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            if isinstance(source, Morphism):
                shape_arrows = (source,)
            else:
                shape_arrows = _flag_arrows(source)
            dimensions = (int(shape_arrows[0].domain().rank()), *(int(arrow.codomain().rank()) for arrow in shape_arrows))
            classes = tuple(category(_flag_marking_from_basis(lattice, basis, dimensions)) for basis in INDUCTION.isotropic_flag_bases(lattice, dimensions))
        case _:
            raise ValueError(f"{category} is not a category of marked lattices this port classifies: a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")
    return finite_ordered_set(classes)


def _flag_marking_from_basis(lattice: Lattices.ParentMethods, basis: tuple[Lattices.ElementMethods, ...], dimensions: tuple[int, ...]) -> tuple[Morphism, ...]:
    r"""The chain of primitive isotropic prefix sublattices of an adapted basis."""
    space = lattice.primitive_isotropic_subobject(*basis)
    composite = space.inclusion()
    arrows = [composite]
    for rank in reversed(dimensions[:-1]):
        term = composite.domain().primitive_isotropic_subobject(*(composite.lift(point) for point in basis[:rank]))
        arrows.insert(0, term.inclusion())
        composite = composite * term.inclusion()
    return tuple(arrows)


def _flag_marking(line: Lattices.ElementMethods, plane: Parent) -> tuple[Morphism, Morphism]:
    r"""The chain Zu -> J -> L of a flag, as the composable pair of embeddings."""
    inclusion = plane.inclusion()
    line_in_plane = inclusion.domain().primitive_isotropic_subobject(inclusion.lift(line))
    return line_in_plane.inclusion(), inclusion


def isotropic_stabilizer_right_cosets(category: Category, marked: Parent) -> Parent:
    r"""A transversal of ``G/H`` for a primitive isotropic sublattice or flag.

    CombinedAlgorithms.h:1565-1609 INDEF_FORM_RightCosets_IsotropicKstuff_Reduced,
    specialized at :2461-2471. ``G`` is the rational group extending the
    source's stabilizer of ``I^perp`` and ``H = G`` meets ``GL(L)``.
    """
    match category:
        case _ if category is ARROWS:
            sublattice = _isotropic_sublattice(marked.arrow())
            basis, images = INDUCTION.plane_vectors(sublattice), INDUCTION.sublattice_action_images(sublattice)
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            _dimensions, basis, images = _flag_action(marked)
        case _:
            raise ValueError(f"{category} is not Ar(Lattices) or Fun([t], Lattices): the source defines these cosets for k-planes and k-flags")
    cosets = RECURSION.isotropic_stabilizer_action(basis, images).right_cosets()
    representative = cosets.representative_map()
    return finite_ordered_set(tuple(~representative(coset) for coset in cosets))


__all__ = [
    "automorphism_group_generators",
    "isometry",
    "isotropic_stabilizer_right_cosets",
    "primitive_marking_classes",
]
