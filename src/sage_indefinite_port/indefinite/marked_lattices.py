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
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
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



def _homothetic_isotropic_saturation(arrow: Morphism) -> Parent | None:
    r"""Primitive saturation of an isotropic ``J=cS``, or ``None`` otherwise.

    The content of the inclusion ``J -> S`` is the gcd of all its
    coordinates.  If this is ``c``, then ``J`` is contained in ``cS``;
    equality holds precisely when ``[S:J]=c**rank``.  This is an exact
    finite-index test, not an assumption that every sublattice is homothetic.
    """
    from math import gcd

    lattice = arrow.codomain()
    points = tuple(arrow(point) for point in arrow.domain().module_generators())
    if not points or any(lattice.b(x, y) != 0 for x in points for y in points):
        return None
    saturation = arrow.saturation()
    inclusion = saturation.inclusion()
    lifts = tuple(inclusion.lift(point) for point in points)
    labels = tuple(inclusion.domain().module_generating_set())
    content = 0
    for point in lifts:
        for label in labels:
            content = gcd(content, int(point.to_vector()(label)))
    if content == 0:
        return None
    subobject = inclusion.domain().subobject_on(lifts)
    if int(subobject.inclusion().cokernel().cardinality()) != content ** len(points):
        return None
    return lattice.primitive_isotropic_subobject(*(inclusion(point) for point in inclusion.domain().module_generators()))


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
                return RECURSION.pointed_isometry(_marked_vector(source), _marked_vector(target))
            _isotropic_sublattice(source.arrow())
            _isotropic_sublattice(target.arrow())
            labels = tuple(category.base_object().module_generating_set())
            basis = tuple(category.base_object().module_generator(label) for label in labels)
            source_images = tuple(source.arrow()(point) for point in basis)
            target_images = tuple(target.arrow()(point) for point in basis)
            # An integral isometry preserves the content/divisibility of each
            # *marked* generator.  Reject an impossible based transporter
            # before constructing its rational reduction and correction group.
            if any(left.divisor() != right.divisor() for left, right in zip(source_images, target_images, strict=True)):
                return None
            # IndefiniteForms.g:1660-1694 with trivial action on the
            # source: every supplied generator is fixed by the marking.
            return RECURSION.isotropic_subspace_isometry(source_images, target_images, ())
        case _ if category is LATTICES:
            return RECURSION.isometry(source, target)
        case _ if category is ARROWS:
            source_rank = int(source.arrow().domain().rank())
            target_rank = int(target.arrow().domain().rank())
            if source_rank == 0 or target_rank == 0:
                return RECURSION.isometry(source.arrow().codomain(), target.arrow().codomain()) if source_rank == target_rank else None
            if source_rank == target_rank == 1:
                source_arrow, target_arrow = source.arrow(), target.arrow()
                source_point = source_arrow(next(iter(source_arrow.domain().module_generators())))
                target_point = target_arrow(next(iter(target_arrow.domain().module_generators())))
                if source_point.content() != target_point.content():
                    return None
                # The source automorphism group of a nonzero rank-one
                # lattice is {+/-1}, whether its square is zero or not.
                # Cancelling content reduces both permitted transports to
                # primitive pointed isometries, preserving the original norm.
                if not source_point.content():
                    return RECURSION.isometry(source_point.parent(), target_point.parent())
                return INDUCTION.line_isometry(source_point.primitive_part(), target_point.primitive_part())
            source_arrow, target_arrow = source.arrow(), target.arrow()
            source_sublattice = _isotropic_sublattice(source_arrow) if source_arrow.is_primitive() else _homothetic_isotropic_saturation(source_arrow)
            target_sublattice = _isotropic_sublattice(target_arrow) if target_arrow.is_primitive() else _homothetic_isotropic_saturation(target_arrow)
            if source_sublattice is None or target_sublattice is None:
                raise ValueError("nonhomothetic nonprimitive moving sublattices require finite-module stabilizer orbits")
            source_point = source_arrow(next(iter(source_arrow.domain().module_generators())))
            target_point = target_arrow(next(iter(target_arrow.domain().module_generators())))
            if (not source_arrow.is_primitive() or not target_arrow.is_primitive()) and source_point.content() != target_point.content():
                return None
            match int(source_sublattice.inclusion().domain().rank()), int(target_sublattice.inclusion().domain().rank()):
                case 1, 1:
                    return INDUCTION.line_isometry(INDUCTION.plane_vectors(source_sublattice)[0], INDUCTION.plane_vectors(target_sublattice)[0])
                case source_rank, target_rank if source_rank == target_rank:
                    return INDUCTION.plane_isometry(source_sublattice, target_sublattice)
                case _:
                    return None
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            source_dimensions, source_basis, _source_images = _flag_action(source)
            target_dimensions, target_basis, target_images = _flag_action(target)
            if source_dimensions != target_dimensions:
                return None
            return RECURSION.isotropic_subspace_isometry(source_basis, target_basis, target_images)
    raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")


def automorphism_group_generators(category: Category, marked: Parent) -> Parent:
    r"""Generators of the automorphism group of a marked lattice in ``category``, as a finite ordered set of isometries of L."""
    match category:
        case CosliceCategory():
            rank = int(category.base_object().rank())
            if rank == 0:
                generators = RECURSION.orthogonal_group_generators(marked.arrow().codomain())
            elif rank == 1:
                generators = RECURSION.pointed_automorphism_generators(_marked_vector(marked))
            else:
                _isotropic_sublattice(marked.arrow())
                # The same perpendicular extension as the setwise case,
                # with trivial radical action (IndefiniteForms.g:1710-1723).
                basis = tuple(marked.arrow()(point) for point in category.base_object().module_generators())
                generators = INDUCTION._stabilizer_generators(basis, ())
        case _ if category is LATTICES:
            generators = RECURSION.orthogonal_group_generators(marked)
        case _ if category is ARROWS:
            if int(marked.arrow().domain().rank()) == 0:
                generators = RECURSION.orthogonal_group_generators(marked.arrow().codomain())
            else:
                arrow = marked.arrow()
                if int(arrow.domain().rank()) == 1:
                    # For any nonzero v, Stab(Zv) is generated by Stab(v)
                    # and, when it exists, one isometry taking v to -v.
                    # Both groups depend only on the primitive saturation.
                    lattice = arrow.codomain()
                    point = arrow(next(iter(arrow.domain().module_generators())))
                    if not point.content():
                        generators = RECURSION.orthogonal_group_generators(lattice)
                    elif point.q() == lattice.base_ring().zero():
                        generators = INDUCTION.line_automorphism_generators(lattice, point.primitive_part())
                    else:
                        primitive = point.primitive_part()
                        fixed = RECURSION.pointed_automorphism_generators(primitive)
                        reversal = RECURSION.pointed_isometry(primitive, -primitive)
                        generators = (*fixed, *((reversal,) if reversal is not None else ()))
                else:
                    sublattice = _isotropic_sublattice(arrow) if arrow.is_primitive() else _homothetic_isotropic_saturation(arrow)
                    if sublattice is None:
                        raise ValueError("nonhomothetic nonprimitive isotropic stabilizers require the finite-module action on saturation")
                    match int(sublattice.inclusion().domain().rank()):
                        case 1:
                            generators = INDUCTION.line_automorphism_generators(arrow.codomain(), INDUCTION.plane_vectors(sublattice)[0])
                        case _:
                            generators = INDUCTION.plane_automorphism_generators(sublattice)
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            _dimensions, basis, images = _flag_action(marked)
            generators = INDUCTION._stabilizer_generators(basis, images)
        case _:
            raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")
    return finite_ordered_set(tuple(generators))


def _pointed_classes(category: Category, lattice: Lattices.ParentMethods, square: int) -> tuple[Parent, ...]:
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



def pointed_content_classes(category: Category, lattice: Parent, square: int, content: int) -> Parent:
    r"""O(L)-classes of rank-one markings of square ``square`` and exact positive content.

    Every vector in this stratum is uniquely ``content * primitive_part``.
    Isometries commute with scaling, so its orbit set is the primitive
    square-``square/content**2`` orbit set.  In particular, specifying
    content makes the isotropic stratum finite when the primitive orbit
    algorithm is finite, without pretending that all nonzero isotropic
    vectors have finitely many orbits across unbounded contents.
    """
    if not isinstance(category, CosliceCategory) or int(category.base_object().rank()) != 1:
        raise ValueError("content strata require a rank-one pointed coslice category")
    source = category.base_object()
    source_generator = next(iter(source.module_generators()))
    if int(source.q(source_generator)) != int(square):
        raise ValueError("the pointed source form must have the requested square")
    if content < 1:
        raise ValueError("nonzero pointed-vector content must be positive")
    if square % (content * content):
        return finite_ordered_set(())
    primitive_square = square // (content * content)
    primitive_source = lattice.lattice_category()([[primitive_square]])
    primitive_category = CosliceCategory(lattice.lattice_category(), primitive_source)
    primitive_classes = _pointed_classes(primitive_category, lattice, primitive_square)
    classes = (
        category.object(source.hom([lattice.scalar_multiple(lattice.base_ring()(content),
            _marked_vector(marking))], codomain=lattice))
        for marking in primitive_classes
    )
    return finite_ordered_set(tuple(classes))



def pointed_nonzero_norm_classes(category: Category, lattice: Parent, square: int) -> Parent:
    r"""All pointed O(L)-orbits of nonzero square, partitioned by exact content.

    The finite list of possible contents is precisely the positive integers
    ``c`` with ``c*c`` dividing ``abs(square)``.  Each stratum uses the
    existing primitive orbit classification, without losing nonprimitive
    representatives or identifying different contents.
    """
    from math import isqrt

    if square == 0:
        raise ValueError("isotropic contents are unbounded; request pointed_content_classes with a prescribed positive content")
    classes = (
        marking
        for content in range(1, isqrt(abs(int(square))) + 1)
        if square % (content * content) == 0
        for marking in pointed_content_classes(category, lattice, square, content)
    )
    return finite_ordered_set(tuple(classes))



def iter_pointed_isotropic_classes(category: Category, lattice: Parent, *, include_zero: bool = False):
    r"""Lazily enumerate every nonzero isotropic pointed orbit, by increasing content.

    For a nondegenerate indefinite lattice admitting primitive isotropic
    vectors, each positive content contributes the primitive zero-norm
    orbit set, scaled by that content.  No finite bound on the content is
    imposed and no finite cardinality is claimed.  With ``include_zero``,
    the zero marking is emitted first as its own orbit.  If there is no
    nonzero rational isotropic vector, only the optional zero orbit remains.
    """
    from itertools import count

    if not isinstance(category, CosliceCategory) or int(category.base_object().rank()) != 1:
        raise ValueError("isotropic content enumeration requires a rank-one pointed coslice")
    source = category.base_object()
    if source.q(next(iter(source.module_generators()))) != source.base_ring().zero():
        raise ValueError("the source of an isotropic pointed marking must have zero square")
    if include_zero:
        yield category.object(source.hom([lattice.zero()], codomain=lattice))
    if not lattice.is_isotropic():
        return
    for content in count(1):
        yield from pointed_content_classes(category, lattice, 0, content)



def nonprimitive_isotropic_line_classes(lattice: Parent, content: int) -> Parent:
    r"""Moving rank-one isotropic sublattice orbits of prescribed index in saturation.

    A rank-one sublattice of index ``content`` in its saturation is
    ``content * S`` for a unique primitive isotropic line ``S``.  The
    operation commutes with integral isometries, so the primitive moving
    line classes determine these orbits without a second classification.
    """
    if content < 1:
        raise ValueError("rank-one sublattice index in saturation must be positive")
    source = lattice.lattice_category()([[0]])
    primitive = primitive_marking_classes(ARROWS, source, lattice)
    classes = (
        ARROWS.object(source.hom([
            lattice.scalar_multiple(lattice.base_ring()(content),
                marking.arrow()(next(iter(marking.arrow().domain().module_generators()))))
        ], codomain=lattice))
        for marking in primitive
    )
    return finite_ordered_set(tuple(classes))



def nonisotropic_line_classes(lattice: Parent, square: int) -> Parent:
    r"""Moving rank-one sublattice orbits of nonzero generator square.

    First enumerate every pointed vector orbit, including imprimitive
    contents.  Passing from pointed vectors to rank-one sublattices
    identifies exactly the orbits related by sign, as the source
    automorphism group is ``{+1,-1}``.  Decide that identification with
    the existing exact moving-arrow transporter, not with a heuristic
    comparison of norms or divisibilities.
    """
    if not square:
        raise ValueError("use nonprimitive_isotropic_line_classes with prescribed index for square zero")
    source = lattice.lattice_category()([[square]])
    category = CosliceCategory(lattice.lattice_category(), source)
    candidates = pointed_nonzero_norm_classes(category, lattice, square)
    kept = []
    for candidate in candidates:
        moving = ARROWS.object(candidate.arrow())
        if all(isometry(ARROWS, moving, previous) is None for previous in kept):
            kept.append(moving)
    return finite_ordered_set(tuple(kept))



def homothetic_isotropic_sublattice_classes(lattice: Parent, rank: int, multiplier: int) -> Parent:
    r"""Moving isotropic sublattice classes ``J=multiplier*S`` for primitive rank-``rank`` ``S``.

    For a free rank-``rank`` saturation, ``S/J`` has elementary divisors
    ``(multiplier,)*rank``.  The map ``S -> multiplier*S`` is equivariant
    under every integral isometry, so it induces a bijection of orbit sets.
    Other elementary-divisor types require the actual finite-module
    stabilizer action and are not silently conflated with this family.
    """
    if rank < 1 or multiplier < 1:
        raise ValueError("rank and saturation multiplier must be positive")
    source = lattice.lattice_category()([[0] * rank for _ in range(rank)])
    primitive = primitive_marking_classes(ARROWS, source, lattice)
    classes = []
    for marking in primitive:
        arrow = marking.arrow()
        images = tuple(
            lattice.scalar_multiple(lattice.base_ring()(multiplier), arrow(point))
            for point in arrow.domain().module_generators()
        )
        classes.append(ARROWS.object(source.hom(list(images), codomain=lattice)))
    return finite_ordered_set(tuple(classes))


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
                    raise ValueError("fixed-source class induction requires a totally isotropic source")
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


__all__ = [
    "automorphism_group_generators",
    "isometry",
    "primitive_marking_classes",
    "pointed_content_classes",
    "pointed_nonzero_norm_classes",
    "iter_pointed_isotropic_classes",
    "nonprimitive_isotropic_line_classes",
    "nonisotropic_line_classes",
    "homothetic_isotropic_sublattice_classes",
]
