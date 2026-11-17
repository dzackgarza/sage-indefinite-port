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



def _finite_index_isotropic_transporter(source_arrow: Morphism, target_arrow: Morphism) -> LatticeIsometryMethods | None:
    r"""Transport two finite-index isotropic markings using the actual saturated stabilizer.

    First carry the primitive saturations onto each other.  Then enumerate
    the finite orbit of the transported source submodule under the integral
    Levi image, retaining an ambient isometry word for each visited member.
    """
    from sage.all import ZZ as SageZZ, matrix

    source_lattice, target_lattice = source_arrow.codomain(), target_arrow.codomain()
    source_saturation = source_arrow.saturation()
    target_saturation = target_arrow.saturation()
    source_basis = tuple(source_saturation.inclusion()(x) for x in source_saturation.module_generators())
    target_basis = tuple(target_saturation.inclusion()(x) for x in target_saturation.module_generators())
    source_space = source_lattice.primitive_isotropic_subobject(*source_basis)
    target_space = target_lattice.primitive_isotropic_subobject(*target_basis)
    transport = INDUCTION.plane_isometry(source_space, target_space)
    if transport is None:
        return None
    inclusion = target_space.inclusion()
    labels = tuple(inclusion.domain().module_generating_set())

    def hnf(points):
        rows = []
        for point in points:
            coordinates = inclusion.lift(point).to_vector()
            rows.append([int(coordinates(label)) for label in labels])
        return matrix(SageZZ, rows).hermite_form()

    source_points = tuple(transport(source_arrow(x)) for x in source_arrow.domain().module_generators())
    target_points = tuple(target_arrow(x) for x in target_arrow.domain().module_generators())
    start, goal = hnf(source_points), hnf(target_points)
    key = lambda h: tuple(tuple(int(v) for v in row) for row in h.rows())
    if start.det() != goal.det() or start.elementary_divisors() != goal.elementary_divisors():
        return None
    goal_key = key(goal)
    if key(start) == goal_key:
        return transport
    actions = []
    for generator in INDUCTION.plane_automorphism_generators(target_space):
        rows = []
        for point in target_basis:
            coordinates = inclusion.lift(generator(point)).to_vector()
            rows.append([int(coordinates(label)) for label in labels])
        actions.append((matrix(SageZZ, rows), generator))
    identity = target_lattice.O().identity()
    seen = {key(start)}
    frontier = [(start, identity)]
    while frontier:
        current, word = frontier.pop()
        for action, generator in actions:
            moved = (current * action).hermite_form()
            moved_key = key(moved)
            if moved_key == goal_key:
                return generator * word * transport
            if moved_key not in seen:
                seen.add(moved_key)
                frontier.append((moved, generator * word))
    return None



def _finite_index_isotropic_stabilizer_generators(arrow: Morphism) -> tuple[LatticeIsometryMethods, ...]:
    r"""Schreier generators of Stab(J) inside Stab(S), for S the isotropic saturation.

    The orbit of the finite-index submodule J in S is finite.  Retain an
    ambient word sending J to each orbit member; Schreier's lemma gives
    generators ``w(gJ)^-1 * g * w(J)`` for its stabilizer.  This uses the
    *actual* integral Levi image, never the full abstract GL(S).
    """
    from sage.all import ZZ as SageZZ, matrix

    lattice = arrow.codomain()
    saturation = arrow.saturation()
    inclusion = saturation.inclusion()
    basis = tuple(inclusion(x) for x in inclusion.domain().module_generators())
    space = lattice.primitive_isotropic_subobject(*basis)
    labels = tuple(inclusion.domain().module_generating_set())
    points = tuple(arrow(x) for x in arrow.domain().module_generators())

    def key(h):
        return tuple(tuple(int(v) for v in row) for row in h.rows())

    start = matrix(SageZZ, [
        [int(inclusion.lift(point).to_vector()(label)) for label in labels]
        for point in points
    ]).hermite_form()
    actions = []
    for generator in INDUCTION.plane_automorphism_generators(space):
        rows = []
        for point in basis:
            coordinates = inclusion.lift(generator(point)).to_vector()
            rows.append([int(coordinates(label)) for label in labels])
        actions.append((matrix(SageZZ, rows), generator))
    identity = lattice.O().identity()
    words = {key(start): identity}
    frontier = [start]
    schreier = []
    while frontier:
        current = frontier.pop()
        source_word = words[key(current)]
        for action, generator in actions:
            moved = (current * action).hermite_form()
            target_key = key(moved)
            if target_key not in words:
                words[target_key] = generator * source_word
                frontier.append(moved)
            schreier.append(words[target_key].inverse() * generator * source_word)
    return tuple(schreier)



def _split_nondegenerate_embedding_isometry(source_arrow: Morphism, target_arrow: Morphism, *, based: bool, source_action=None, complement_action=None) -> LatticeIsometryMethods | None:
    r"""Extend source and perpendicular isometries across an integral orthogonal direct sum.

    For L=S+S^perp as *integral* orthogonal summands, independent isometries
    of S and S^perp extend uniquely.  If the decomposition has nontrivial
    glue, do not infer a nonexistent extension from rational splitting.
    """
    from sage.all import ZZ as SageZZ, matrix

    source_lattice, target_lattice = source_arrow.codomain(), target_arrow.codomain()
    if source_action is None and complement_action is None and source_lattice is target_lattice and source_arrow.domain() is target_arrow.domain() and all(
        source_arrow(point) == target_arrow(point) for point in source_arrow.domain().module_generators()
    ):
        return source_lattice.O().identity()
    source_sub = source_lattice.subobject_on(tuple(source_arrow(x) for x in source_arrow.domain().module_generators()))
    target_sub = target_lattice.subobject_on(tuple(target_arrow(x) for x in target_arrow.domain().module_generators()))
    if not source_sub.inclusion().is_primitive() or not target_sub.inclusion().is_primitive():
        raise ValueError("the split nondegenerate extension requires primitive embeddings")
    source_complement = source_lattice.orthogonal_complement(source_sub)
    target_complement = target_lattice.orthogonal_complement(target_sub)
    source_points = tuple(source_arrow(x) for x in source_arrow.domain().module_generators())
    target_points = tuple(target_arrow(x) for x in target_arrow.domain().module_generators())
    source_cpoints = tuple(source_complement.inclusion()(x) for x in source_complement.module_generators())
    target_cpoints = tuple(target_complement.inclusion()(x) for x in target_complement.module_generators())
    for lattice, points in ((source_lattice, (*source_points, *source_cpoints)),
                            (target_lattice, (*target_points, *target_cpoints))):
        if int(lattice.subobject_on(points).inclusion().cokernel().cardinality()) != 1:
            raise ValueError("the primitive embedding has nontrivial complement glue")
    source_iso = source_action if source_action is not None else (None if based else RECURSION.isometry(source_arrow.domain(), target_arrow.domain()))
    if not based and source_iso is None:
        return None
    complement_iso = complement_action if complement_action is not None else RECURSION.isometry(source_complement.inclusion().domain(), target_complement.inclusion().domain())
    if complement_iso is None:
        return None
    images = (
        tuple(target_arrow(x if based else source_iso(x)) for x in source_arrow.domain().module_generators())
        + tuple(target_complement.inclusion()(complement_iso(x)) for x in source_complement.module_generators())
    )
    labels = tuple(source_lattice.module_generating_set())
    source_matrix = matrix(SageZZ, [
        [int(point.to_vector()(label)) for label in labels]
        for point in (*source_points, *source_cpoints)
    ])
    inverse = source_matrix.inverse().change_ring(SageZZ)
    target_basis = target_lattice.base_ring()
    def image(label):
        column = inverse.row(labels.index(label))
        return sum((target_lattice.scalar_multiple(target_basis(int(coefficient)), point)
                    for coefficient, point in zip(column, images, strict=True) if coefficient), target_lattice.zero())
    return source_lattice.Isom(target_lattice)(image)



def _split_nondegenerate_stabilizer_generators(arrow: Morphism, *, based: bool) -> tuple[LatticeIsometryMethods, ...]:
    r"""The split stabilizer O(S) x O(S^perp), or its based kernel O(S^perp)."""
    lattice = arrow.codomain()
    subobject = lattice.subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators()))
    complement = lattice.orthogonal_complement(subobject).inclusion().domain()
    complement_identity = complement.O().identity()
    source_identity = arrow.domain().O().identity()
    family = []
    if not based:
        for action in RECURSION.orthogonal_group_generators(arrow.domain()):
            family.append(_split_nondegenerate_embedding_isometry(
                arrow, arrow, based=False, source_action=action, complement_action=complement_identity
            ))
    for action in RECURSION.orthogonal_group_generators(complement):
        family.append(_split_nondegenerate_embedding_isometry(
            arrow, arrow, based=based, source_action=source_identity,
            complement_action=action
        ))
    return tuple(family)



def _finite_glued_nondegenerate_stabilizer_generators(arrow: Morphism, *, based: bool) -> tuple[LatticeIsometryMethods, ...]:
    r"""Glue-compatible stabilizer when both orthogonal factors are definite.

    The port supplies the complete finite factor groups; the lattice
    embedding owner performs the integral glue extension.
    """
    from sage_extensions.categories.lattice_embeddings import (
        _glue_compatible_integral_extensions,
    )

    lattice = arrow.codomain()
    subobject = lattice.subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators()))
    if not subobject.inclusion().is_primitive() or not subobject.inclusion().domain().is_nondegenerate():
        raise ValueError("finite glue stabilization requires a primitive nondegenerate source")
    complement = lattice.orthogonal_complement(subobject).inclusion().domain()
    if not arrow.domain().is_definite() or not complement.is_definite():
        raise ValueError("infinite orthogonal factors require the discriminant-image glue stabilizer")
    source_actions = (arrow.domain().O().identity(),) if based else tuple(arrow.domain().O())
    return _glue_compatible_integral_extensions(arrow, arrow, source_actions, tuple(complement.O()))


def _nondegenerate_stabilizer_generators(arrow: Morphism, *, based: bool) -> tuple[LatticeIsometryMethods, ...]:
    lattice = arrow.codomain()
    subobject = lattice.subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators()))
    complement = lattice.orthogonal_complement(subobject)
    basis = tuple(arrow(x) for x in arrow.domain().module_generators())
    complementary = tuple(complement.inclusion()(x) for x in complement.module_generators())
    combined = lattice.subobject_on((*basis, *complementary))
    if int(combined.inclusion().cokernel().cardinality()) == 1:
        return _split_nondegenerate_stabilizer_generators(arrow, based=based)
    return _finite_glued_nondegenerate_stabilizer_generators(arrow, based=based)



def _nondegenerate_embedding_isometry(source_arrow: Morphism, target_arrow: Morphism, *, based: bool) -> LatticeIsometryMethods | None:
    r"""Split extension or exact finite-factor glue-compatible transporter."""
    from sage_extensions.categories.lattice_embeddings import (
        _glue_compatible_integral_extensions,
    )

    source_lattice, target_lattice = source_arrow.codomain(), target_arrow.codomain()
    if source_lattice is target_lattice and source_arrow.domain() is target_arrow.domain() and all(
        source_arrow(x) == target_arrow(x) for x in source_arrow.domain().module_generators()
    ):
        return source_lattice.O().identity()
    source_sub = source_lattice.subobject_on(tuple(source_arrow(x) for x in source_arrow.domain().module_generators()))
    target_sub = target_lattice.subobject_on(tuple(target_arrow(x) for x in target_arrow.domain().module_generators()))
    source_comp = source_lattice.orthogonal_complement(source_sub)
    target_comp = target_lattice.orthogonal_complement(target_sub)
    source_points = tuple(source_arrow(x) for x in source_arrow.domain().module_generators())
    target_points = tuple(target_arrow(x) for x in target_arrow.domain().module_generators())
    source_cpoints = tuple(source_comp.inclusion()(x) for x in source_comp.module_generators())
    target_cpoints = tuple(target_comp.inclusion()(x) for x in target_comp.module_generators())
    source_index = int(source_lattice.subobject_on((*source_points, *source_cpoints)).inclusion().cokernel().cardinality())
    target_index = int(target_lattice.subobject_on((*target_points, *target_cpoints)).inclusion().cokernel().cardinality())
    if source_index == target_index == 1:
        return _split_nondegenerate_embedding_isometry(source_arrow, target_arrow, based=based)
    if not (source_arrow.domain().is_definite() and target_arrow.domain().is_definite()
            and source_comp.inclusion().domain().is_definite() and target_comp.inclusion().domain().is_definite()):
        raise ValueError("nontrivial glue with an infinite orthogonal factor requires the finite discriminant-image transporter")
    source_iso = None if based else RECURSION.isometry(source_arrow.domain(), target_arrow.domain())
    if not based and source_iso is None:
        return None
    complement_iso = RECURSION.isometry(source_comp.inclusion().domain(), target_comp.inclusion().domain())
    if complement_iso is None:
        return None
    source_actions = (
        (source_arrow.domain().O().identity(),) if based
        else tuple(g * source_iso for g in target_arrow.domain().O())
    )
    complement_actions = tuple(g * complement_iso for g in target_comp.inclusion().domain().O())
    extensions = _glue_compatible_integral_extensions(
        source_arrow, target_arrow, source_actions, complement_actions
    )
    return extensions[0] if extensions else None


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
            source_arrow, target_arrow = source.arrow(), target.arrow()
            if source_arrow.domain().is_nondegenerate():
                return _nondegenerate_embedding_isometry(source_arrow, target_arrow, based=True)
            source_saturation = _isotropic_sublattice(source_arrow) if source_arrow.is_primitive() else _homothetic_isotropic_saturation(source_arrow)
            target_saturation = _isotropic_sublattice(target_arrow) if target_arrow.is_primitive() else _homothetic_isotropic_saturation(target_arrow)
            if source_saturation is None or target_saturation is None:
                raise ValueError("nonhomothetic based isotropic markings require an integral source-module action")
            labels = tuple(category.base_object().module_generating_set())
            basis = tuple(category.base_object().module_generator(label) for label in labels)
            source_images = tuple(source.arrow()(point) for point in basis)
            target_images = tuple(target.arrow()(point) for point in basis)
            # An integral isometry preserves the content/divisibility of each
            # *marked* generator.  Reject an impossible based transporter
            # before constructing its rational reduction and correction group.
            if any(left.divisor() != right.divisor() for left, right in zip(source_images, target_images, strict=True)):
                return None
            source_contents = tuple(int(point.content()) for point in source_images)
            target_contents = tuple(int(point.content()) for point in target_images)
            if source_contents != target_contents:
                return None
            # For homothetic J=cS the marked generators become a primitive
            # basis of S after cancelling c.  No source GL action is allowed.
            source_primitive = tuple(point.primitive_part() for point in source_images)
            target_primitive = tuple(point.primitive_part() for point in target_images)
            return RECURSION.isotropic_subspace_isometry(source_primitive, target_primitive, ())
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
            if source_arrow.domain().is_nondegenerate() and target_arrow.domain().is_nondegenerate():
                return _nondegenerate_embedding_isometry(source_arrow, target_arrow, based=False)
            source_sublattice = _isotropic_sublattice(source_arrow) if source_arrow.is_primitive() else _homothetic_isotropic_saturation(source_arrow)
            target_sublattice = _isotropic_sublattice(target_arrow) if target_arrow.is_primitive() else _homothetic_isotropic_saturation(target_arrow)
            if source_sublattice is None or target_sublattice is None:
                return _finite_index_isotropic_transporter(source_arrow, target_arrow)
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
                arrow = marked.arrow()
                if arrow.domain().is_nondegenerate():
                    generators = _nondegenerate_stabilizer_generators(arrow, based=True)
                    return finite_ordered_set(tuple(generators))
                saturation = _isotropic_sublattice(arrow) if arrow.is_primitive() else _homothetic_isotropic_saturation(arrow)
                if saturation is None:
                    raise ValueError("nonhomothetic based isotropic stabilizers require an integral source-module action")
                basis = tuple(arrow(point).primitive_part() for point in category.base_object().module_generators())
                generators = INDUCTION._stabilizer_generators(basis, ())
        case _ if category is LATTICES:
            generators = RECURSION.orthogonal_group_generators(marked)
        case _ if category is ARROWS:
            if int(marked.arrow().domain().rank()) == 0:
                generators = RECURSION.orthogonal_group_generators(marked.arrow().codomain())
            else:
                arrow = marked.arrow()
                if int(arrow.domain().rank()) > 1 and arrow.domain().is_nondegenerate():
                    generators = _nondegenerate_stabilizer_generators(arrow, based=False)
                    return finite_ordered_set(tuple(generators))
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
                        generators = _finite_index_isotropic_stabilizer_generators(arrow)
                    else:
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
    # The binary split primitive algorithm is even-form based.  For an odd
    # split plane, double the pairing and transport full-group orbits back;
    # never transport the doubled discriminant-kernel subgroup.
    if not lattice.is_even() and int(lattice.rank()) == 2 and lattice.is_isotropic():
        return tuple(odd_primitive_pointed_classes_via_doubling(category, lattice, square))
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




def odd_primitive_pointed_classes_via_doubling(category: Category, lattice: Parent, square: int) -> Parent:
    r"""Primitive full-O(L) pointed classes on an odd lattice via its even double.

    The underlying free module and the *full* integral isometry group are
    unchanged by doubling the pairing.  Classify square ``2*square`` on
    ``(L,2b)`` and transport each vector's unchanged coordinates back to
    ``(L,b)``.  No discriminant-kernel, stable or congruence subgroup is
    transported from the doubled lattice: those must be imposed on L.
    """
    if lattice.is_even():
        raise ValueError("this reduction is specifically for odd integral lattices")
    if not isinstance(category, CosliceCategory) or int(category.base_object().rank()) != 1:
        raise ValueError("odd primitive pointed classes require a rank-one coslice")
    source = category.base_object()
    if int(source.q(next(iter(source.module_generators())))) != square:
        raise ValueError("the source marking square must match the requested square")
    doubled, _twist = lattice.doubled_form()
    doubled_source = doubled.lattice_category()([[2 * square]])
    doubled_category = CosliceCategory(doubled.lattice_category(), doubled_source)
    labels = tuple(lattice.module_generating_set())
    classes = []
    for marking in _pointed_classes(doubled_category, doubled, 2 * square):
        coordinates = _marked_vector(marking).to_vector()
        point = lattice.linear_combination({label: coordinates(label) for label in labels if coordinates(label)})
        classes.append(category.object(source.hom([point], codomain=lattice)))
    return finite_ordered_set(tuple(classes))


def pointed_content_classes(category: Category, lattice: Parent, square: int, content: int, *, acting_group=None) -> Parent:
    r"""O(L)-classes of rank-one markings of square ``square`` and exact positive content.

    Every vector in this stratum is uniquely ``content * primitive_part``.
    Isometries commute with scaling, so its orbit set is the primitive
    square-``square/content**2`` orbit set.  In particular, specifying
    content makes the isotropic stratum finite when the primitive orbit
    algorithm is finite, without pretending that all nonzero isotropic
    vectors have finitely many orbits across unbounded contents.  An
    ``acting_group`` splits each primitive orbit under that *same* subgroup
    before scaling; it is never silently replaced by the full O(L).
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
    if acting_group is not None and acting_group is not lattice.O():
        primitive_classes = tuple(
            split for marking in primitive_classes
            for split in lattice.O().split_orbit(marking, acting_group)
        )
    classes = (
        category.object(source.hom([lattice.scalar_multiple(lattice.base_ring()(content),
            _marked_vector(marking))], codomain=lattice))
        for marking in primitive_classes
    )
    return finite_ordered_set(tuple(classes))



def pointed_nonzero_norm_classes(category: Category, lattice: Parent, square: int, *, acting_group=None) -> Parent:
    r"""All pointed O(L)-orbits of nonzero square, partitioned by exact content.

    The finite list of possible contents is precisely the positive integers
    ``c`` with ``c*c`` dividing ``abs(square)``.  Each stratum uses the
    existing primitive orbit classification, without losing nonprimitive
    representatives or identifying different contents.  The optional
    ``acting_group`` is retained across every content stratum.
    """
    from math import isqrt

    if square == 0:
        raise ValueError("isotropic contents are unbounded; request pointed_content_classes with a prescribed positive content")
    classes = (
        marking
        for content in range(1, isqrt(abs(int(square))) + 1)
        if square % (content * content) == 0
        for marking in pointed_content_classes(category, lattice, square, content, acting_group=acting_group)
    )
    return finite_ordered_set(tuple(classes))



def iter_pointed_isotropic_classes(category: Category, lattice: Parent, *, include_zero: bool = False, acting_group=None):
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
        yield from pointed_content_classes(category, lattice, 0, content, acting_group=acting_group)



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



def nonprimitive_sublattice_classes(saturated: Parent | Morphism, divisors: tuple[int, ...]) -> tuple[Morphism, ...]:
    r"""Finite-index sublattice classes under the actual ambient stabilizer image.

    The port supplies the integral action on the primitive saturation.  The
    generic finite-index submodule orbit construction belongs to the module
    owner and returns represented subobjects, not matrix surrogates.
    """
    if isinstance(saturated, Morphism):
        lattice = saturated.codomain()
        saturated = lattice.subobject_on(tuple(saturated(point) for point in saturated.domain().module_generators()))
    inclusion = saturated.inclusion()
    lattice = inclusion.codomain()
    if not inclusion.is_primitive():
        raise ValueError("the supplied saturation must be primitive")
    basis = tuple(inclusion(point) for point in inclusion.domain().module_generators())
    isotropic = all(lattice.b(x, y) == 0 for x in basis for y in basis)
    if isotropic:
        ambient_generators = INDUCTION.plane_automorphism_generators(saturated)
        module = inclusion.domain()
        actions = tuple(module.Isom(module)(
            lambda label, generator=g: inclusion.lift(generator(inclusion(module.module_generator(label))))
        ) for g in ambient_generators)
    else:
        if not inclusion.domain().is_nondegenerate():
            raise ValueError("nonisotropic degenerate saturation requires its radical stabilizer action")
        complement = lattice.orthogonal_complement(saturated)
        combined = lattice.subobject_on((*basis, *(complement.inclusion()(x) for x in complement.module_generators())))
        if int(combined.inclusion().cokernel().cardinality()) != 1:
            raise ValueError("nonsplit nondegenerate saturation requires the complement-glue stabilizer")
        # Every O(S)-element extends by identity on its orthogonal complement.
        actions = RECURSION.orthogonal_group_generators(inclusion.domain())
    submodules = inclusion.domain().finite_index_submodule_orbits(divisors, actions)
    return tuple(lattice.subobject_on(tuple(inclusion(subobject(point))
        for point in subobject.domain().module_generators())).inclusion() for subobject in submodules)


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
                if source.is_nondegenerate():
                    # The owned primitive-embedding classification='emb'
                    # fixes the source marking pointwise.  Classification
                    # 'sub' would silently quotient by O(source) instead.
                    classes = tuple(category.object(embedding) for embedding in
                        source.Emb(lattice).primitive_embedding_class_representatives(classification="emb"))
                else:
                    if any(source.b(left, right) != 0 for left in source_basis for right in source_basis):
                        raise ValueError("mixed degenerate Gram markings require their radical/glue classification")
                    classes = tuple(
                        category.object(source.hom(list(images), codomain=lattice))
                        for basis in INDUCTION.isotropic_flag_bases(lattice, (int(source.rank()),))
                        for images in INDUCTION.basis_marking_classes(lattice, basis)
                    )
        case _ if category is ARROWS:
            rank = int(source.rank())
            if rank == 0:
                classes = (ARROWS.object(source.hom([], codomain=lattice)),)
            elif source.is_nondegenerate():
                # Moving sublattices allow the full source isometry group.
                # The owned 'sub' classification quotients by both O(source)
                # and O(target), unlike the based 'emb' classification.
                classes = tuple(ARROWS.object(embedding) for embedding in
                    source.Emb(lattice).primitive_embedding_class_representatives(classification="sub"))
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
    "odd_primitive_pointed_classes_via_doubling",
    "pointed_nonzero_norm_classes",
    "iter_pointed_isotropic_classes",
    "nonprimitive_isotropic_line_classes",
    "nonisotropic_line_classes",
    "homothetic_isotropic_sublattice_classes",
    "nonprimitive_sublattice_classes",
]
