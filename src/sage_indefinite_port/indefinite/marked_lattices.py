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


def _zero_form_orthogonal_generators(ambient, linear_generators):
    r"""Carry integral GL generators into the owned O(L) homset for q=0."""
    from sage.all import ZZ as SageZZ, matrix
    from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

    labels = tuple(ambient.module_generating_set())
    return tuple(_isometry_from_row_matrix(ambient, ambient, matrix(SageZZ, [
        [int(generator(point).to_vector()(label)) for label in labels]
        for point in ambient.module_generators()
    ])) for generator in linear_generators)


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



def _finite_index_sublattice_transporter(source_arrow: Morphism, target_arrow: Morphism) -> LatticeIsometryMethods | None:
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
    source_space = source_saturation
    target_space = target_saturation
    if source_space.inclusion().domain().is_nondegenerate() and target_space.inclusion().domain().is_nondegenerate():
        transport = _nondegenerate_embedding_isometry(source_space.inclusion(), target_space.inclusion(), based=False)
    else:
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
    target_generators = (
        _nondegenerate_stabilizer_generators(target_space.inclusion(), based=False)
        if target_space.inclusion().domain().is_nondegenerate()
        else INDUCTION.plane_automorphism_generators(target_space)
    )
    for generator in target_generators:
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



def _finite_index_sublattice_stabilizer_generators(arrow: Morphism) -> tuple[LatticeIsometryMethods, ...]:
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
    space = (saturation if inclusion.domain().is_nondegenerate()
             or any(lattice.b(x, y) != 0 for x in basis for y in basis)
             else lattice.primitive_isotropic_subobject(*basis))
    labels = tuple(inclusion.domain().module_generating_set())
    points = tuple(arrow(x) for x in arrow.domain().module_generators())

    def key(h):
        return tuple(tuple(int(v) for v in row) for row in h.rows())

    start = matrix(SageZZ, [
        [int(inclusion.lift(point).to_vector()(label)) for label in labels]
        for point in points
    ]).hermite_form()
    actions = []
    if all(lattice.b(x, y) == 0 for x in lattice.module_generators() for y in lattice.module_generators()):
        # O(L)=GL(L).  The saturation stabilizer is the integral
        # parabolic, even when J has nonconstant elementary divisors.
        # Schreier below intersects it with the *actual* finite-index
        # submodule stabilizer; a Witt-index induction is inapplicable.
        linear = lattice.general_linear_group()
        preserving = (linear.group_generators() if int(inclusion.domain().rank()) == int(lattice.rank())
                      else linear.flag_stabilizer((saturation,)).group_generators())
        ambient_generators = _zero_form_orthogonal_generators(lattice, preserving)
    else:
        ambient_generators = (
            _nondegenerate_stabilizer_generators(space.inclusion(), based=False)
            if space.inclusion().domain().is_nondegenerate()
            else (_degenerate_radical_flag_generators(lattice, (len(basis),), basis)
                  if not lattice.is_nondegenerate() and all(
                      RECURSION.vector_lies_in_radical(point) for point in basis
                  )
                  else (_split_radical_saturation_stabilizer_generators(space)
                        if not all(lattice.b(x, y) == 0 for x in basis for y in basis)
                        else INDUCTION.plane_automorphism_generators(space)))
        )
    for generator in ambient_generators:
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



def _nondegenerate_block_group(arrow: Morphism, *, based: bool):
    r"""Full glue stabilizer from exact finite-image preimage of block O(S)xO(K).

    The block group may be infinite.  The finite-index embedding's owned
    ``extended_subgroup`` computes the *whole* preimage of the glue
    stabilizer in the actual discriminant image, including Schreier kernel
    generators.  Filtering block generators individually is insufficient.
    """
    from sage.all import QQ as SageQQ, matrix

    lattice = arrow.codomain()
    subobject = lattice.subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators()))
    complement = lattice.orthogonal_complement(subobject)
    inclusion = lattice.subobject_on((
        *(arrow(x) for x in arrow.domain().module_generators()),
        *(complement.inclusion()(x) for x in complement.module_generators()),
    )).inclusion()
    module = inclusion.domain()
    source_basis = tuple(arrow(x) for x in arrow.domain().module_generators())
    comp_basis = tuple(complement.inclusion()(x) for x in complement.module_generators())
    frame = (*source_basis, *comp_basis)
    labels = tuple(lattice.module_generating_set())
    coordinates = matrix(SageQQ, [[int(point.to_vector()(label)) for label in labels]
        for point in frame]).inverse()
    mlabels = tuple(module.module_generating_set())

    def block(source_action, complement_action):
        images = (
            *(arrow(source_action(x)) for x in arrow.domain().module_generators()),
            *(complement.inclusion()(complement_action(x)) for x in complement.module_generators()),
        )
        def image(label):
            point = inclusion(module.module_generator(label))
            coefficients = point.to_vector()
            row = matrix(SageQQ, [[int(coefficients(key)) for key in labels]]) * coordinates
            result = sum((lattice.scalar_multiple(lattice.base_ring()(int(value)), element)
                for value, element in zip(row.row(0), images, strict=True) if value), lattice.zero())
            return inclusion.lift(result)
        return module.O()(image)

    source_identity = arrow.domain().O().identity()
    comp_domain = complement.inclusion().domain()
    comp_identity = comp_domain.O().identity()
    blocks = []
    if not based:
        blocks.extend(block(g, comp_identity) for g in RECURSION.orthogonal_group_generators(arrow.domain()))
    blocks.extend(block(source_identity, g) for g in RECURSION.orthogonal_group_generators(comp_domain))
    return inclusion, module.O().subgroup(tuple(blocks))



def _glued_nondegenerate_stabilizer_generators(arrow: Morphism, *, based: bool) -> tuple[LatticeIsometryMethods, ...]:
    inclusion, blocks = _nondegenerate_block_group(arrow, based=based)
    return tuple(inclusion.extended_subgroup(blocks).group_generators())


def _nondegenerate_stabilizer_generators(arrow: Morphism, *, based: bool) -> tuple[LatticeIsometryMethods, ...]:
    lattice = arrow.codomain()
    subobject = lattice.subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators()))
    complement = lattice.orthogonal_complement(subobject)
    basis = tuple(arrow(x) for x in arrow.domain().module_generators())
    complementary = tuple(complement.inclusion()(x) for x in complement.module_generators())
    combined = lattice.subobject_on((*basis, *complementary))
    if int(combined.inclusion().cokernel().cardinality()) == 1:
        return _split_nondegenerate_stabilizer_generators(arrow, based=based)
    return _glued_nondegenerate_stabilizer_generators(arrow, based=based)




def _glued_nondegenerate_embedding_isometry(source_arrow: Morphism, target_arrow: Morphism, *, based: bool, source_action=None) -> LatticeIsometryMethods | None:
    r"""Glue-compatible transporter via the actual finite target block action.

    A base isometry of the two orthogonal sums transports source glue to a
    subgroup of the target discriminant module.  The target block group
    acts on the finite set of glue subgroups; its *actual* transporter
    supplies the missing arithmetic factor word.  The resulting rational
    map extends integrally exactly when the glue subgroups agree.
    """
    from sage.all import QQ as SageQQ, matrix
    from sage_extensions.categories.group.g_sets import FiniteGSets
    from sage_extensions.categories.lattice_embeddings import _finite_overlattice_isometry_extension

    source_lattice, target_lattice = source_arrow.codomain(), target_arrow.codomain()
    source_sub = source_lattice.subobject_on(tuple(source_arrow(x) for x in source_arrow.domain().module_generators()))
    target_sub = target_lattice.subobject_on(tuple(target_arrow(x) for x in target_arrow.domain().module_generators()))
    source_comp = source_lattice.orthogonal_complement(source_sub)
    target_comp = target_lattice.orthogonal_complement(target_sub)
    source_iso = (source_action if source_action is not None else
                  (source_arrow.domain().O().identity() if based else RECURSION.isometry(source_arrow.domain(), target_arrow.domain())))
    comp_iso = RECURSION.isometry(source_comp.inclusion().domain(), target_comp.inclusion().domain())
    if source_iso is None or comp_iso is None:
        return None
    source_basis = (
        *(source_arrow(x) for x in source_arrow.domain().module_generators()),
        *(source_comp.inclusion()(x) for x in source_comp.module_generators()),
    )
    target_images = (
        *(target_arrow(source_iso(x)) for x in source_arrow.domain().module_generators()),
        *(target_comp.inclusion()(comp_iso(x)) for x in source_comp.module_generators()),
    )
    source_inclusion = source_lattice.subobject_on(source_basis).inclusion()
    target_inclusion, target_blocks = _nondegenerate_block_group(target_arrow, based=based)
    source_module, target_module = source_inclusion.domain(), target_inclusion.domain()
    labels = tuple(source_lattice.module_generating_set())
    frame = matrix(SageQQ, [[int(point.to_vector()(label)) for label in labels]
        for point in source_basis]).inverse()

    def base_image(label):
        coordinates = source_inclusion(source_module.module_generator(label)).to_vector()
        row = matrix(SageQQ, [[int(coordinates(key)) for key in labels]]) * frame
        result = sum((target_lattice.scalar_multiple(target_lattice.base_ring()(int(value)), point)
            for value, point in zip(row.row(0), target_images, strict=True) if value), target_lattice.zero())
        return target_inclusion.lift(result)

    base = source_module.Isom(target_module)(base_image)
    source_glue = source_inclusion.gluing_subgroup().inclusion()
    target_glue = target_inclusion.gluing_subgroup().inclusion()
    target_form = target_module.discriminant_group()
    points = target_form.elements()
    action = FiniteGSets(target_blocks)(
        points, lambda element, point: element.discriminant_isometry()(point)
    )
    source_discriminant_map = base.discriminant_isometry()
    moved_source = tuple(source_discriminant_map(source_glue(x)) for x in source_glue.domain().elements())
    target_members = tuple(target_glue(x) for x in target_glue.domain().elements())
    correction = action.setwise_transporter_witness(moved_source, target_members)
    if correction is None:
        return None
    return _finite_overlattice_isometry_extension(source_inclusion, target_inclusion, correction * base)



def _nonprimitive_nondegenerate_based_isometry(source_arrow: Morphism, target_arrow: Morphism) -> LatticeIsometryMethods | None:
    r"""Extend the unique rational saturation map determined by a fixed source basis.

    A based isometry must send every generator of J to its corresponding
    target generator.  Since J has full rank in its primitive saturation,
    this determines a unique rational map S_Q -> T_Q.  It extends only
    if its matrix on S,T is integral unimodular, and its ambient glue
    admits a compatible complement action.
    """
    from sage.all import QQ as SageQQ, ZZ as SageZZ, matrix

    source_saturation = source_arrow.saturation().inclusion()
    target_saturation = target_arrow.saturation().inclusion()
    source_module, target_module = source_saturation.domain(), target_saturation.domain()
    source_labels = tuple(source_module.module_generating_set())
    target_labels = tuple(target_module.module_generating_set())
    source_points = tuple(source_arrow(x) for x in source_arrow.domain().module_generators())
    target_points = tuple(target_arrow(x) for x in target_arrow.domain().module_generators())
    if len(source_points) != len(target_points):
        return None
    a = matrix(SageQQ, [[int(source_saturation.lift(x).to_vector()(label)) for label in source_labels]
        for x in source_points])
    b = matrix(SageQQ, [[int(target_saturation.lift(x).to_vector()(label)) for label in target_labels]
        for x in target_points])
    if not a.is_square() or not a.det():
        return None
    action = a.inverse() * b
    if abs(action.det()) != 1 or not all(entry in SageZZ for entry in action.list()):
        return None
    integral = action.change_ring(SageZZ)
    def image(label):
        return target_module.linear_combination({target_label: target_module.base_ring()(int(value))
            for target_label, value in zip(target_labels, integral.row(source_labels.index(label)), strict=True) if value})
    saturation_isometry = source_module.Isom(target_module)(image)
    return _glued_nondegenerate_embedding_isometry(source_saturation, target_saturation,
        based=True, source_action=saturation_isometry)


def _nondegenerate_embedding_isometry(source_arrow: Morphism, target_arrow: Morphism, *, based: bool) -> LatticeIsometryMethods | None:
    r"""Split extension or exact finite-factor glue-compatible transporter."""
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
    return _glued_nondegenerate_embedding_isometry(source_arrow, target_arrow, based=based)




def _degenerate_radical_based_isometry(source_arrow: Morphism, target_arrow: Morphism):
    r"""Exact based integral radical-frame transporter by column Hermite equivalence.

    For row frames A,B in the saturated radical R, an ambient isometry
    carrying each marked generator to its counterpart exists exactly when
    A D = B for some D in GL(R).  This is equivalent to equality of the
    column lattices of A and B in the source frame coordinates.  Hermite
    transformations produce D and the reduced split extends it by the
    identity on the nondegenerate quotient.
    """
    from sage.all import ZZ as SageZZ, identity_matrix, matrix
    from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

    source_ambient, target_ambient = source_arrow.codomain(), target_arrow.codomain()
    ambient_transport = (source_ambient.O().identity() if source_ambient is target_ambient
                         else RECURSION.isometry(source_ambient, target_ambient))
    if ambient_transport is None:
        return None
    reduction, n, k, _quotient = _degenerate_quotient_and_radical(target_ambient)
    labels = tuple(reduction.reduced.module_generating_set())
    inverse = reduction.isometry.inverse()

    def rows(arrow, *, transport=None):
        vectors = tuple(arrow(point) for point in arrow.domain().module_generators())
        if transport is not None:
            vectors = tuple(transport(point) for point in vectors)
        coordinates = tuple(tuple(int(inverse(point).to_vector()(label)) for label in labels)
                            for point in vectors)
        if any(any(row[:n]) for row in coordinates):
            return None
        return matrix(SageZZ, [row[n:] for row in coordinates])

    left = rows(source_arrow, transport=ambient_transport)
    right = rows(target_arrow)
    if left is None or right is None or left.nrows() != right.nrows():
        return None
    left_hnf, left_change = left.transpose().hermite_form(transformation=True)
    right_hnf, right_change = right.transpose().hermite_form(transformation=True)
    if left_hnf != right_hnf:
        return None
    radical_action = left_change.transpose() * right_change.transpose().inverse()
    if left * radical_action != right or abs(radical_action.det()) != 1:
        raise ArithmeticError("the Hermite transporter failed its integral frame certificate")
    block = identity_matrix(SageZZ, n + k)
    for i in range(k):
        for j in range(k):
            block[n + i, n + j] = radical_action[i, j]
    reduced_witness = _isometry_from_row_matrix(reduction.reduced, reduction.reduced, block)
    correction = reduction.isometry * reduced_witness * reduction.isometry.inverse()
    return correction * ambient_transport


def _degenerate_full_quotient_based_isometry(source_arrow: Morphism, target_arrow: Morphism):
    r"""Exact based mixed-Gram transport for a full unimodular quotient frame.

    Write L=C+R in the reduced integral chart.  Suppose the marked
    nonradical rows project to a Z-basis of C and the marked radical rows
    form a primitive free direct summand of R.  Their images determine the
    quotient action A and an integral radical completion D; the unique
    shear B is X^-1(Y_rad-X_rad D).  The block [[A,B],[0,D]] is integral
    unimodular and preserves the form precisely when the marked Gram
    tensors agree.  This gives an actual ambient isometry, not a finite
    search or a rational extension.
    """
    from sage.matrix.constructor import identity_matrix, matrix
    from sage.rings.integer_ring import ZZ as SageZZ
    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
    from sage_indefinite_port.indefinite.recursive import _gram_in_label_order, _isometry_from_row_matrix, _conjugated

    lattice = source_arrow.codomain()
    if target_arrow.codomain() is not lattice:
        transport = RECURSION.isometry(lattice, target_arrow.codomain())
        if transport is None:
            return None
        moved = source_arrow.domain().hom([transport(source_arrow(x)) for x in source_arrow.domain().module_generators()],
                                          codomain=target_arrow.codomain())
        correction = _degenerate_full_quotient_based_isometry(moved, target_arrow)
        return None if correction is None else correction * transport
    reduction = indefinite_reduced_presentation(lattice)
    reduced = reduction.reduced
    gram = _gram_in_label_order(reduced)
    n, k = int(gram.rank()), int(gram.nrows() - gram.rank())
    labels = tuple(reduced.module_generating_set())
    inverse = reduction.isometry.inverse()
    left = tuple(source_arrow(x) for x in source_arrow.domain().module_generators())
    right = tuple(target_arrow(x) for x in target_arrow.domain().module_generators())
    if len(left) != len(right):
        return None
    radical_indices = tuple(i for i, point in enumerate(left) if RECURSION.vector_lies_in_radical(point))
    nonradical_indices = tuple(i for i in range(len(left)) if i not in radical_indices)
    r = len(radical_indices)
    if len(nonradical_indices) != n or not 1 <= r <= k:
        raise ValueError("mixed based Gram transport requires a full quotient basis and primitive radical frame")
    if any(not RECURSION.vector_lies_in_radical(right[i]) for i in radical_indices):
        return None

    def coordinates(points):
        return matrix(SageZZ, [[int(inverse(point).to_vector()(label)) for label in labels] for point in points])

    left_non = coordinates(tuple(left[i] for i in nonradical_indices))
    right_non = coordinates(tuple(right[i] for i in nonradical_indices))
    left_rad = coordinates(tuple(left[i] for i in radical_indices))[:, n:]
    right_rad = coordinates(tuple(right[i] for i in radical_indices))[:, n:]
    x, y = left_non[:, :n], right_non[:, :n]
    if not x.det() or not y.det():
        raise ValueError("the nonradical marked vectors must span the rational nondegenerate quotient")
    if (abs(int(x.det())) != 1 or abs(int(y.det())) != 1) and r != k:
        raise ValueError("nonprimitive quotient with a proper radical frame requires constrained integral stabilizer action")
    a_rational = x.inverse() * y
    if abs(a_rational.det()) != 1 or any(value not in SageZZ for value in a_rational.list()):
        return None
    a = a_rational.change_ring(SageZZ)
    left_smith, u_left, v_left = left_rad.smith_form()
    right_smith, u_right, v_right = right_rad.smith_form()
    if left_smith != right_smith:
        return None
    # E M = U_left U_right^-1 E for E=diag(d_1,...,d_r).
    # The top-left block of M is therefore uniquely forced to be
    # diag(d)^-1 U_left U_right^-1 diag(d).  Integrality of that block
    # is the exact obstruction to a based integral radical-frame map;
    # equality of Smith invariants alone is insufficient.
    from sage.rings.rational_field import QQ as SageQQ
    frame_change = u_left * u_right.inverse().change_ring(SageZZ)
    diagonal = matrix(SageQQ, r, r, lambda i, j: int(left_smith[i, i]) if i == j else 0)
    required = diagonal.inverse() * frame_change.change_ring(SageQQ) * diagonal
    if any(value not in SageZZ for value in required.list()):
        return None
    middle = identity_matrix(SageZZ, k)
    for i in range(r):
        for j in range(r):
            middle[i, j] = SageZZ(required[i, j])
    d = v_left * middle * v_right.inverse().change_ring(SageZZ)
    if left_rad * d != right_rad or abs(int(d.det())) != 1:
        raise ArithmeticError("the integral radical completion does not carry the marked frame")
    b_rational = x.inverse() * (right_non[:, n:] - left_non[:, n:] * d)
    if any(value not in SageZZ for value in b_rational.list()):
        return None
    b = b_rational.change_ring(SageZZ)
    block = identity_matrix(SageZZ, n + k)
    for i in range(n):
        for j in range(n):
            block[i, j] = a[i, j]
        for j in range(k):
            block[i, n + j] = b[i, j]
    for i in range(k):
        for j in range(k):
            block[n + i, n + j] = d[i, j]
    if block * gram * block.transpose() != gram:
        return None
    witness = _conjugated(reduction.isometry, (_isometry_from_row_matrix(reduced, reduced, block),))[0]
    if any(witness(first) != second for first, second in zip(left, right, strict=True)):
        raise ArithmeticError("the full-quotient mixed Gram isometry fails its based marking")
    return witness


def _degenerate_rank_two_based_isometry(source_arrow: Morphism, target_arrow: Morphism):
    r"""Based mixed rank-two markings with a radical generator and primitive quotient image.

    Transport the radical generator first.  A primitive nonradical quotient
    vector has complement content one, so its remaining radical difference
    is corrected by a shear with D=1, fixing the already transported
    radical generator.  Both factors are actual integral ambient maps.
    """
    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
    from sage_indefinite_port.indefinite.recursive import _gram_in_label_order
    source_lattice, target_lattice = source_arrow.codomain(), target_arrow.codomain()
    if source_lattice is not target_lattice:
        ambient = RECURSION.isometry(source_lattice, target_lattice)
        if ambient is None:
            return None
        source_arrow = source_arrow.domain().hom(
            [ambient(source_arrow(x)) for x in source_arrow.domain().module_generators()],
            codomain=target_lattice)
        correction = _degenerate_rank_two_based_isometry(source_arrow, target_arrow)
        return None if correction is None else correction * ambient
    lattice = source_lattice
    left = tuple(source_arrow(x) for x in source_arrow.domain().module_generators())
    right = tuple(target_arrow(x) for x in target_arrow.domain().module_generators())
    radical = [i for i, (a, b) in enumerate(zip(left, right, strict=True))
               if RECURSION.vector_lies_in_radical(a) and RECURSION.vector_lies_in_radical(b)]
    if len(radical) != 1:
        raise ValueError("mixed degenerate rank-two source requires exactly one marked radical generator")
    r = radical[0]
    p = 1 - r
    reduction = indefinite_reduced_presentation(lattice)
    gram = _gram_in_label_order(reduction.reduced)
    n = int(gram.rank())
    labels = tuple(reduction.reduced.module_generating_set())[:n]
    coordinates = reduction.isometry.inverse()(left[p]).to_vector()
    from math import gcd
    content = 0
    for label in labels:
        content = gcd(content, int(coordinates(label)))
    radical_rank = int(gram.nrows()) - n
    first = RECURSION.pointed_isometry(left[r], right[r])
    if first is None:
        return None
    if content == 1:
        correction = RECURSION.pointed_isometry(first(left[p]), right[p])
    else:
        inverse = reduction.isometry.inverse()
        all_labels = tuple(reduction.reduced.module_generating_set())
        fixed = tuple(int(inverse(right[r]).to_vector()(label)) for label in all_labels[n:])
        start = tuple(int(inverse(first(left[p])).to_vector()(label)) for label in all_labels[n:])
        goal = tuple(int(inverse(right[p]).to_vector()(label)) for label in all_labels[n:])
        d = RECURSION.radical_stabilizer_congruence_image(lattice, content, fixed).transporter(start, goal)
        correction = (None if d is None else RECURSION._mixed_degenerate_pointed_isometry(
            first(left[p]), right[p], radical_action=d
        ))
    if correction is None:
        return None
    result = correction * first
    if result(left[r]) != right[r] or result(left[p]) != right[p]:
        raise ArithmeticError("the mixed rank-two transporter must fix both marked generators")
    return result


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
            ambient = source_arrow.codomain()
            # Pure zero-form frames are included in the split-radical
            # Hermite transporter below (with zero-dimensional quotient).
            # No unbounded enumeration of GL_n(ZZ) is necessary.
            if not ambient.is_nondegenerate() and all(
                RECURSION.vector_lies_in_radical(source_arrow(point))
                for point in source_arrow.domain().module_generators()
            ):
                if target_arrow.codomain().is_nondegenerate():
                    return None
                return _degenerate_radical_based_isometry(source_arrow, target_arrow)
            if source_arrow.domain().is_nondegenerate():
                return (_nondegenerate_embedding_isometry(source_arrow, target_arrow, based=True)
                        if source_arrow.is_primitive() and target_arrow.is_primitive()
                        else _nonprimitive_nondegenerate_based_isometry(source_arrow, target_arrow))
            if (not source_arrow.domain().is_nondegenerate()
                    and any(source_arrow.domain().b(x, y) != 0
                            for x in source_arrow.domain().module_generators()
                            for y in source_arrow.domain().module_generators())):
                if int(source_arrow.domain().rank()) == 2:
                    return _degenerate_rank_two_based_isometry(source_arrow, target_arrow)
                return _degenerate_full_quotient_based_isometry(source_arrow, target_arrow)
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
            if source_rank != target_rank:
                return None
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
            ambient = source_arrow.codomain()
            if (int(ambient.rank()) == int(target_arrow.codomain().rank())
                    and all(ambient.b(x, y) == 0 for x in ambient.module_generators() for y in ambient.module_generators())
                    and all(target_arrow.codomain().b(x, y) == 0 for x in target_arrow.codomain().module_generators()
                            for y in target_arrow.codomain().module_generators())):
                from sage.all import ZZ as SageZZ, matrix
                from sage_indefinite_port.indefinite.recursive import _fair_unimodular_matrices, _isometry_from_row_matrix

                source_labels = tuple(ambient.module_generating_set())
                target_labels = tuple(target_arrow.codomain().module_generating_set())
                left = tuple(source_arrow(point) for point in source_arrow.domain().module_generators())
                right = tuple(target_arrow(point) for point in target_arrow.domain().module_generators())
                left_rows = matrix(SageZZ, [[int(point.to_vector()(label)) for label in source_labels] for point in left])
                right_rows = matrix(SageZZ, [[int(point.to_vector()(label)) for label in target_labels] for point in right])
                if left_rows.elementary_divisors() != right_rows.elementary_divisors():
                    return None
                target_subobject = target_arrow.codomain().subobject_on(right)
                for rows in _fair_unimodular_matrices(int(ambient.rank())):
                    witness = _isometry_from_row_matrix(ambient, target_arrow.codomain(), rows)
                    if all(target_subobject.inclusion().is_in_image(witness(point)) for point in left):
                        return witness
            if source_arrow.domain().is_nondegenerate() and target_arrow.domain().is_nondegenerate():
                return (_nondegenerate_embedding_isometry(source_arrow, target_arrow, based=False)
                        if source_arrow.is_primitive() and target_arrow.is_primitive()
                        else _finite_index_sublattice_transporter(source_arrow, target_arrow))
            source_sublattice = _isotropic_sublattice(source_arrow) if source_arrow.is_primitive() else _homothetic_isotropic_saturation(source_arrow)
            target_sublattice = _isotropic_sublattice(target_arrow) if target_arrow.is_primitive() else _homothetic_isotropic_saturation(target_arrow)
            if source_sublattice is None or target_sublattice is None:
                return _finite_index_sublattice_transporter(source_arrow, target_arrow)
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
            source_arrows = _flag_arrows(source)
            target_arrows = _flag_arrows(target)
            if len(source_arrows) == len(target_arrows) == 2 and all(
                tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2)
                for arrows in (source_arrows, target_arrows)
            ):
                source_ambient = source_arrows[-1].codomain()
                target_ambient = target_arrows[-1].codomain()
                if all(
                    not ambient.is_nondegenerate()
                    and int(arrows[-1].domain().rank()) == int(ambient.radical().rank())
                    and all(RECURSION.vector_lies_in_radical(arrows[-1](point))
                            for point in arrows[-1].domain().module_generators())
                    for ambient, arrows in ((source_ambient, source_arrows),
                                            (target_ambient, target_arrows))
                ):
                    source_line = source_arrows[-1](source_arrows[0](
                        next(iter(source_arrows[0].domain().module_generators()))
                    ))
                    target_line = target_arrows[-1](target_arrows[0](
                        next(iter(target_arrows[0].domain().module_generators()))
                    ))
                    # The terminal planes are the whole radicals, hence
                    # preserved by every ambient isometry. The flag
                    # transporter is exactly a transporter of the two
                    # primitive radical lines, up to source sign.
                    witness = RECURSION.pointed_isometry(source_line, target_line)
                    return witness if witness is not None else RECURSION.pointed_isometry(source_line, -target_line)
            source_dimensions, source_basis, _source_images = _flag_action(source)
            target_dimensions, target_basis, target_images = _flag_action(target)
            if source_dimensions != target_dimensions:
                return None
            ambient = source_basis[0].parent()
            target_ambient = target_basis[0].parent()
            if (int(ambient.rank()) == int(target_ambient.rank())
                    and all(ambient.b(x, y) == 0 for x in ambient.module_generators() for y in ambient.module_generators())
                    and all(target_ambient.b(x, y) == 0 for x in target_ambient.module_generators()
                            for y in target_ambient.module_generators())):
                # On a zero-form lattice GL_n(ZZ) acts transitively on
                # primitive flags of any fixed rank profile. A fair
                # unimodular search therefore terminates positively; each
                # candidate is certified against every marked flag term.
                from sage_indefinite_port.indefinite.recursive import _fair_unimodular_matrices, _isometry_from_row_matrix

                target_terms = tuple(target_ambient.primitive_isotropic_subobject(*target_basis[:rank])
                                     for rank in source_dimensions)
                for rows in _fair_unimodular_matrices(int(ambient.rank())):
                    witness = _isometry_from_row_matrix(ambient, target_ambient, rows)
                    if all(term.inclusion().is_in_image(witness(point))
                           for rank, term in zip(source_dimensions, target_terms, strict=True)
                           for point in source_basis[:rank]):
                        return witness
            return RECURSION.isotropic_subspace_isometry(source_basis, target_basis, target_images)
    raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")


def automorphism_group_generators(category: Category, marked: Parent) -> Parent:
    r"""Generators of the automorphism group of a marked lattice in ``category``, as a finite ordered set of isometries of L."""
    match category:
        case CosliceCategory():
            rank = int(category.base_object().rank())
            if rank == 0:
                ambient = marked.arrow().codomain()
                generators = (_zero_form_orthogonal_generators(ambient, ambient.general_linear_group().group_generators())
                              if all(ambient.b(x, y) == 0 for x in ambient.module_generators()
                                     for y in ambient.module_generators())
                              else RECURSION.orthogonal_group_generators(ambient))
            elif rank == 1:
                generators = RECURSION.pointed_automorphism_generators(_marked_vector(marked))
            else:
                arrow = marked.arrow()
                if all(
                    arrow.codomain().b(x, y) == 0
                    for x in arrow.codomain().module_generators()
                    for y in arrow.codomain().module_generators()
                ):
                    from sage.all import ZZ as SageZZ, identity_matrix
                    from sage_indefinite_port.indefinite.recursive import _fair_unimodular_matrices, _isometry_from_row_matrix

                    ambient = arrow.codomain()
                    n, k = int(ambient.rank()), int(arrow.domain().rank())
                    labels = tuple(ambient.module_generating_set())
                    # Pointwise fixation of an integral frame is equivalent
                    # to pointwise fixation of its primitive saturation:
                    # if m*x is fixed, torsion-freeness forces x fixed.
                    frame_arrow = arrow if arrow.is_primitive() else arrow.saturation().inclusion()
                    fixed = tuple(frame_arrow(point) for point in frame_arrow.domain().module_generators())
                    frame = tuple(tuple(int(point.to_vector()(label)) for label in labels) for point in fixed)
                    completion = next(rows for rows in _fair_unimodular_matrices(n)
                                      if tuple(tuple(int(entry) for entry in row) for row in rows.rows()[:k]) == frame)
                    inverse = completion.inverse()
                    # In the completed frame the pointwise stabilizer is
                    # Z^((n-k)k) semidirect GL_(n-k)(Z): lower-left shears
                    # and lower-right unimodular generators suffice.
                    # Conjugation by the integral completion transports
                    # these exact generators to the original coordinates.
                    blocks = []
                    for i in range(k, n):
                        for j in range(k):
                            shear = identity_matrix(SageZZ, n)
                            shear[i, j] = 1
                            blocks.append(shear)
                    if n > k:
                        complement = ambient.lattice_category()([[0] * (n - k) for _ in range(n - k)])
                        complement_labels = tuple(complement.module_generating_set())
                        for generator in complement.general_linear_group().group_generators():
                            block = identity_matrix(SageZZ, n)
                            for i, point in enumerate(complement.module_generators()):
                                for j, label in enumerate(complement_labels):
                                    block[k + i, k + j] = int(generator(point).to_vector()(label))
                            blocks.append(block)
                    generators = tuple(_isometry_from_row_matrix(ambient, ambient, inverse * block * completion)
                                       for block in blocks) or (ambient.O().identity(),)
                    return finite_ordered_set(tuple(generators))
                if not arrow.codomain().is_nondegenerate() and all(
                    RECURSION.vector_lies_in_radical(arrow(point))
                    for point in arrow.domain().module_generators()
                ):
                    # Pointwise fixation of a full-rank finite-index frame
                    # is equivalent, over QQ, to pointwise fixation of its
                    # primitive saturation; torsion-freeness gives the
                    # integral equivalence.  No homothetic-index assumption.
                    frame = arrow if arrow.is_primitive() else arrow.saturation().inclusion()
                    basis = tuple(frame(point) for point in frame.domain().module_generators())
                    generators = _degenerate_radical_flag_generators(
                        arrow.codomain(), (int(frame.domain().rank()),), basis, pointwise=True
                    )
                    return finite_ordered_set(tuple(generators))
                elif arrow.domain().is_nondegenerate():
                    generators = (_nondegenerate_stabilizer_generators(arrow, based=True)
                                  if arrow.is_primitive() else _nondegenerate_stabilizer_generators(arrow.saturation().inclusion(), based=True))
                    return finite_ordered_set(tuple(generators))
                saturation = _isotropic_sublattice(arrow) if arrow.is_primitive() else _homothetic_isotropic_saturation(arrow)
                if saturation is None:
                    raise ValueError("nonhomothetic based isotropic stabilizers require an integral source-module action")
                basis = tuple(arrow(point).primitive_part() for point in category.base_object().module_generators())
                generators = INDUCTION._stabilizer_generators(basis, ())
        case _ if category is LATTICES:
            generators = (_zero_form_orthogonal_generators(marked, marked.general_linear_group().group_generators())
                          if all(marked.b(x, y) == 0 for x in marked.module_generators()
                                 for y in marked.module_generators())
                          else RECURSION.orthogonal_group_generators(marked))
        case _ if category is ARROWS:
            if int(marked.arrow().domain().rank()) == 0:
                ambient = marked.arrow().codomain()
                generators = (_zero_form_orthogonal_generators(ambient, ambient.general_linear_group().group_generators())
                              if all(ambient.b(x, y) == 0 for x in ambient.module_generators()
                                     for y in ambient.module_generators())
                              else RECURSION.orthogonal_group_generators(ambient))
            else:
                arrow = marked.arrow()
                if all(
                    arrow.codomain().b(x, y) == 0
                    for x in arrow.codomain().module_generators()
                    for y in arrow.codomain().module_generators()
                ):
                    ambient = arrow.codomain()
                    if not arrow.is_primitive():
                        return finite_ordered_set(tuple(_finite_index_sublattice_stabilizer_generators(arrow)))
                    subobject = ambient.primitive_isotropic_subobject(*(
                        arrow(point) for point in arrow.domain().module_generators()
                    ))
                    generators = _zero_form_orthogonal_generators(
                        ambient, ambient.general_linear_group().flag_stabilizer((subobject,)).group_generators()
                    )
                    return finite_ordered_set(tuple(generators))
                if not arrow.codomain().is_nondegenerate() and arrow.is_primitive() and all(
                    RECURSION.vector_lies_in_radical(arrow(point))
                    for point in arrow.domain().module_generators()
                ):
                    basis = tuple(arrow(point) for point in arrow.domain().module_generators())
                    generators = _degenerate_radical_flag_generators(
                        arrow.codomain(), (int(arrow.domain().rank()),), basis
                    )
                    return finite_ordered_set(tuple(generators))
                if not arrow.codomain().is_nondegenerate() and arrow.is_primitive():
                    from sage.all import identity_matrix
                    from sage.rings.integer_ring import ZZ as SageZZ
                    from sage_indefinite_port.indefinite.recursive import _conjugated, _isometry_from_row_matrix

                    ambient = arrow.codomain()
                    reduction, n, k, quotient = _degenerate_quotient_and_radical(ambient)
                    if quotient.is_isotropic() and int(quotient.witt_index()) == 1 and int(arrow.domain().rank()) == k + 1:
                        reduced = reduction.reduced
                        labels = tuple(reduced.module_generating_set())
                        quotient_labels = tuple(quotient.module_generating_set())
                        inverse = reduction.isometry.inverse()
                        images = tuple(arrow(point) for point in arrow.domain().module_generators())
                        projected = tuple(quotient.linear_combination({
                            quotient_labels[i]: quotient.base_ring()(int(inverse(point).to_vector()(labels[i])))
                            for i in range(n)
                        }) for point in images)
                        nonzero = next((point for point in projected if point != quotient.zero()), None)
                        if nonzero is not None and all(quotient.b(nonzero, point) == 0 for point in projected):
                            radical_basis = tuple(reduction.isometry(reduced.module_generator(label))
                                                  for label in labels[n:])
                            generators = [generator for generator in _degenerate_radical_flag_generators(
                                ambient, (k,), radical_basis
                            ) if all(arrow.is_in_image(generator(point)) for point in images)]
                            line = quotient.primitive_isotropic_subobject(nonzero.primitive_part())
                            quotient_marked = ARROWS.object(line.inclusion())
                            for action in automorphism_group_generators(ARROWS, quotient_marked):
                                block = identity_matrix(SageZZ, n + k)
                                for i, label in enumerate(quotient_labels):
                                    coordinates = action(quotient.module_generator(label)).to_vector()
                                    for j, target_label in enumerate(quotient_labels):
                                        block[i, j] = int(coordinates(target_label))
                                generators.extend(_conjugated(reduction.isometry, (
                                    _isometry_from_row_matrix(reduced, reduced, block),
                                )))
                            return finite_ordered_set(tuple(generators))
                if int(arrow.domain().rank()) > 1 and arrow.domain().is_nondegenerate():
                    generators = (_nondegenerate_stabilizer_generators(arrow, based=False)
                                  if arrow.is_primitive() else _finite_index_sublattice_stabilizer_generators(arrow))
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
                elif not arrow.domain().is_nondegenerate() and any(
                    arrow.domain().b(x, y) != 0 for x in arrow.domain().module_generators()
                    for y in arrow.domain().module_generators()
                ):
                    generators = (_split_radical_saturation_stabilizer_generators(
                        arrow.codomain().subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators())))
                        if arrow.is_primitive() else _finite_index_sublattice_stabilizer_generators(arrow))
                else:
                    sublattice = _isotropic_sublattice(arrow) if arrow.is_primitive() else _homothetic_isotropic_saturation(arrow)
                    if sublattice is None:
                        generators = _finite_index_sublattice_stabilizer_generators(arrow)
                    else:
                        match int(sublattice.inclusion().domain().rank()):
                            case 1:
                                generators = INDUCTION.line_automorphism_generators(arrow.codomain(), INDUCTION.plane_vectors(sublattice)[0])
                            case _:
                                generators = INDUCTION.plane_automorphism_generators(sublattice)
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            arrows = _flag_arrows(marked)
            if len(arrows) == 2 and tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2):
                first_arrow, last_arrow = arrows
                ambient = last_arrow.codomain()
                if not ambient.is_nondegenerate() and all(
                    RECURSION.vector_lies_in_radical(last_arrow(point))
                    for point in last_arrow.domain().module_generators()
                ):
                    from sage.all import ZZ as SageZZ, matrix

                    plane = last_arrow.domain()
                    labels = tuple(plane.module_generating_set())
                    line = first_arrow(next(iter(first_arrow.domain().module_generators())))
                    row = tuple(int(line.to_vector()(label)) for label in labels)
                    _, _, change = matrix(SageZZ, [row]).smith_form()
                    adapted = change.inverse()
                    basis = tuple(last_arrow(plane.linear_combination({
                        label: plane.base_ring()(int(adapted[i, j]))
                        for j, label in enumerate(labels)
                    })) for i in range(2))
                    generators = _degenerate_radical_flag_generators(ambient, (1, 2), basis)
                    return finite_ordered_set(tuple(generators))
            dimensions, basis, images = _flag_action(marked)
            ambient = basis[0].parent()
            if all(ambient.b(x, y) == 0 for x in ambient.module_generators() for y in ambient.module_generators()):
                terms = tuple(ambient.primitive_isotropic_subobject(*basis[:rank]) for rank in dimensions)
                generators = _zero_form_orthogonal_generators(
                    ambient, ambient.general_linear_group().flag_stabilizer(terms).group_generators()
                )
            elif not ambient.is_nondegenerate() and all(
                RECURSION.vector_lies_in_radical(point) for point in basis
            ):
                generators = _degenerate_radical_flag_generators(ambient, dimensions, basis)
            else:
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




def _degenerate_quotient_and_radical(lattice: Parent):
    r"""The reduced nondegenerate quotient and integral radical coordinates.

    The pinned unimodular reduction orders the nondegenerate block before
    the saturated radical.  All marked classifiers use this same framing.
    """
    from sage_indefinite_port.indefinite.recursive import _gram_in_label_order

    reduction = indefinite_reduced_presentation(lattice)
    gram = _gram_in_label_order(reduction.reduced)
    n = int(gram.rank())
    quotient = (lattice.lattice_category()([
        [int(gram[i, j]) for j in range(n)] for i in range(n)
    ]) if n else None)
    return reduction, n, int(gram.nrows()) - n, quotient


def _degenerate_radical_flag_generators(ambient: Parent, dimensions: tuple[int, ...], basis: tuple, *, pointwise: bool = False):
    r"""Exact stabilizer of a primitive flag inside the saturated radical.

    In the reduced split C + R, O(L) is (O(C) × GL(R)) ⋉ Hom(C,R).
    The flag imposes precisely its GL(R)-stabilizer, leaving O(C) and
    every shear unrestricted.  The pinned reduction transports the result.
    """
    from sage.all import identity_matrix, matrix
    from sage.rings.integer_ring import ZZ as SageZZ
    from sage_indefinite_port.indefinite.recursive import _conjugated, _isometry_from_row_matrix
    from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int

    reduction, n, k, _quotient = _degenerate_quotient_and_radical(ambient)
    reduced = reduction.reduced
    labels = tuple(reduced.module_generating_set())
    radical = ambient.lattice_category()([[0] * k for _ in range(k)])
    radical_labels = tuple(radical.module_generating_set())
    inverse = reduction.isometry.inverse()
    radical_basis = tuple(radical.linear_combination({
        label: radical.base_ring()(int(inverse(point).to_vector()(labels[n + i])))
        for i, label in enumerate(radical_labels)
    }) for point in basis)
    radical_identity = identity_matrix(SageZZ, k)
    first = reduced.module_generator(labels[n])
    blocks = [block for block in RECURSION._radical_blocks(reduced, first)
              if block.matrix_from_rows_and_columns(range(n, n + k), range(n, n + k))
              == radical_identity]
    basis_rows = matrix(SageZZ, [
        [int(point.to_vector()(label)) for label in radical_labels]
        for point in radical_basis
    ])
    adapted = basis_rows.stack(_source_subspace_completion_int(basis_rows, k))
    if abs(adapted.det()) != 1:
        raise ArithmeticError("the primitive radical flag does not admit an integral adapted completion")
    if pointwise:
        if len(dimensions) != 1 or dimensions[0] != len(radical_basis):
            raise ValueError("a based radical marking fixes its entire ordered frame")
        rank = dimensions[0]
        standard_matrices = []
        for i in range(rank, k):
            for j in range(rank):
                shear = identity_matrix(SageZZ, k)
                shear[i, j] = 1
                standard_matrices.append(shear)
        if k > rank:
            quotient = radical.lattice_category()([[0] * (k - rank) for _ in range(k - rank)])
            quotient_labels = tuple(quotient.module_generating_set())
            for action in quotient.general_linear_group().group_generators():
                block = identity_matrix(SageZZ, k)
                for i, label in enumerate(quotient_labels):
                    image = action(quotient.module_generator(label)).to_vector()
                    for j, target_label in enumerate(quotient_labels):
                        block[rank + i, rank + j] = int(image(target_label))
                standard_matrices.append(block)
    else:
        standard_basis = tuple(radical.module_generator(label) for label in radical_labels)
        standard_terms = tuple(radical.primitive_isotropic_subobject(*standard_basis[:rank])
                               for rank in dimensions)
        standard_matrices = [matrix(SageZZ, [
            [int(action(radical.module_generator(label)).to_vector()(target_label))
             for target_label in radical_labels]
            for label in radical_labels
        ]) for action in radical.general_linear_group().flag_stabilizer(standard_terms).group_generators()]
    radical_actions = tuple(adapted.inverse() * action * adapted for action in standard_matrices)
    for action in radical_actions:
        block = identity_matrix(SageZZ, n + k)
        for i in range(k):
            for j in range(k):
                block[n + i, n + j] = action[i, j]
        blocks.append(block)
    return _conjugated(reduction.isometry, tuple(
        _isometry_from_row_matrix(reduced, reduced, block) for block in blocks
    ))


def _degenerate_radical_residues(lattice: Parent, radical_rank: int, modulus: int):
    """Canonical orbit representatives of (Z/modulus)^k under the actual GL_k(Z) image."""
    from itertools import product

    if modulus == 1:
        return ((0,) * radical_rank,)
    image = RECURSION.radical_congruence_image(lattice, modulus)
    unseen = set(product(range(modulus), repeat=radical_rank))
    representatives = []
    while unseen:
        representative = min(unseen)
        representatives.append(representative)
        unseen.difference_update(tuple(int(value) % modulus for value in row)
            for row in image.orbit_of(representative))
    return tuple(representatives)


def _degenerate_isotropic_pointed_iterator(category: Category, lattice: Parent, *, acting_group=None):
    r"""Fair infinite enumeration by complement content and finite radical residue orbits."""
    from itertools import count
    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
    from sage_indefinite_port.indefinite.recursive import _gram_in_label_order

    reduction = indefinite_reduced_presentation(lattice)
    reduced = reduction.reduced
    gram = _gram_in_label_order(reduced)
    n, k = int(gram.rank()), int(gram.nrows() - gram.rank())
    labels = tuple(reduced.module_generating_set())
    primitive = ()
    if n:
        complement = reduced.lattice_category()([[int(gram[i, j]) for j in range(n)] for i in range(n)])
        if not complement.is_definite() and complement.is_isotropic():
            primitive = tuple(complement.O().vector_orbit_representatives(0, content=1))
        c_labels = tuple(complement.module_generating_set())
    source = category.base_object()

    def emit(point):
        marked = category.object(source.hom([reduction.isometry(point)], codomain=lattice))
        if acting_group is None:
            yield marked
        else:
            yield from lattice.O().split_orbit(marked, acting_group)

    for content in count(1):
        # The pure radical has exactly one GL(R)-orbit at each content.
        radical = reduced.scalar_multiple(reduced.base_ring()(content), reduced.module_generator(labels[n]))
        yield from emit(radical)
        if not primitive:
            continue
        for point in primitive:
            coordinates = point.to_vector()
            for residue in _degenerate_radical_residues(lattice, k, content):
                values = (
                    *(content * int(coordinates(label)) for label in c_labels),
                    *residue,
                )
                reduced_point = reduced.linear_combination({label: reduced.base_ring()(value)
                    for label, value in zip(labels, values, strict=True) if value})
                yield from emit(reduced_point)


def _degenerate_nonzero_pointed_classes(category: Category, lattice: Parent, square: int, *, acting_group=None) -> Parent:
    r"""All nonzero-norm pointed orbits through C+rad(L) and exact finite residues.

    A complement vector x of content c determines a finite radical
    residue r modulo cR.  Shears change r by cR, while GL(R) acts through
    its *actual* integral congruence image.  The full vector is primitive
    precisely when gcd(c,r)=1; nonprimitive vectors are retained as well.
    """
    from math import isqrt
    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
    from sage_indefinite_port.indefinite.recursive import _gram_in_label_order

    if not square:
        raise ValueError("zero norm requires content-indexed lazy radical and isotropic families")
    reduction = indefinite_reduced_presentation(lattice)
    reduced = reduction.reduced
    gram = _gram_in_label_order(reduced)
    n, k = int(gram.rank()), int(gram.nrows() - gram.rank())
    if not n:
        return finite_ordered_set(())
    complement = reduced.lattice_category()([[int(gram[i, j]) for j in range(n)] for i in range(n)])
    labels = tuple(reduced.module_generating_set())
    classes = []
    for content in range(1, isqrt(abs(square)) + 1):
        if square % (content * content):
            continue
        complement_square = square // (content * content)
        primitive = (
            tuple(point for point in complement.O().vector_orbit_representatives(complement_square)
                  if int(point.content()) == 1)
            if complement.is_definite()
            else complement.O().vector_orbit_representatives(complement_square, content=1)
        )
        residues = _degenerate_radical_residues(lattice, k, content)
        for point in primitive:
            coordinates = point.to_vector()
            complement_labels = tuple(complement.module_generating_set())
            for residue in residues:
                coefficients = (
                    *(content * int(coordinates(label)) for label in complement_labels),
                    *residue,
                )
                reduced_point = reduced.linear_combination({label: reduced.base_ring()(value)
                    for label, value in zip(labels, coefficients, strict=True) if value})
                classes.append(category.object(category.base_object().hom([
                    reduction.isometry(reduced_point)
                ], codomain=lattice)))
    if acting_group is not None:
        classes = [split for marked in classes for split in lattice.O().split_orbit(marked, acting_group)]
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
    if int(lattice.rank()) == 0:
        return finite_ordered_set(())
    if square == 0 and not lattice.is_nondegenerate() and any(
        lattice.b(x, y) != 0 for x in lattice.module_generators()
        for y in lattice.module_generators()
    ):
        reduction, n, _k, complement = _degenerate_quotient_and_radical(lattice)
        if not complement.is_definite() and complement.is_isotropic():
            raise ValueError("mixed degenerate isotropic fixed-content orbit sets may be infinite; use iter_pointed_content_classes")
        # An anisotropic nondegenerate quotient forces every null vector
        # into the radical, where GL(rad L) is transitive at fixed content.
        radical_generator = reduction.reduced.module_generator(tuple(reduction.reduced.module_generating_set())[n])
        vector = reduction.isometry(reduction.reduced.scalar_multiple(
            reduction.reduced.base_ring()(content), radical_generator
        ))
        marked = category.object(source.hom([vector], codomain=lattice))
        if acting_group is None:
            return finite_ordered_set((marked,))
        return finite_ordered_set(tuple(lattice.O().split_orbit(marked, acting_group)))
    if not lattice.is_nondegenerate() and square:
        classes = _degenerate_nonzero_pointed_classes(category, lattice, square, acting_group=acting_group)
        return finite_ordered_set(tuple(marked for marked in classes if int(_marked_vector(marked).content()) == content))
    if square % (content * content):
        return finite_ordered_set(())
    primitive_square = square // (content * content)
    primitive_source = lattice.lattice_category()([[primitive_square]])
    primitive_category = CosliceCategory(lattice.lattice_category(), primitive_source)
    if square == 0 and all(lattice.b(x, y) == 0 for x in lattice.module_generators()
                           for y in lattice.module_generators()):
        # GL_n(ZZ) is transitive on primitive nonzero vectors of a
        # zero-form lattice. The indefinite Witt-index dispatcher is not
        # defined for this degenerate form.
        primitive_classes = tuple(primitive_marking_classes(primitive_category, primitive_source, lattice))
    else:
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

    if not isinstance(category, CosliceCategory) or int(category.base_object().rank()) != 1:
        raise ValueError("nonzero-norm pointed orbits require a rank-one coslice category")
    source = category.base_object()
    if int(source.q(next(iter(source.module_generators())))) != int(square):
        raise ValueError("the pointed source form must have the requested square")
    if square == 0:
        raise ValueError("isotropic contents are unbounded; request pointed_content_classes with a prescribed positive content")
    if int(lattice.rank()) == 0:
        return finite_ordered_set(())
    if not lattice.is_nondegenerate() and square:
        return _degenerate_nonzero_pointed_classes(category, lattice, square, acting_group=acting_group)
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
    if int(lattice.rank()) == 0:
        return
    if not lattice.is_nondegenerate():
        yield from _degenerate_isotropic_pointed_iterator(category, lattice, acting_group=acting_group)
        return
    if not lattice.is_isotropic():
        return
    for content in count(1):
        yield from pointed_content_classes(category, lattice, 0, content, acting_group=acting_group)




def iter_pointed_content_classes(category: Category, lattice: Parent, square: int, content: int, *, acting_group=None):
    r"""Enumerate exact pointed content even when its orbit set is infinite.

    Mixed degenerate zero-norm lattices may have infinitely many primitive
    orbits at content one. Filter the exhaustive fair stream by full-vector
    content rather than asserting a finite classification.
    """
    if not isinstance(category, CosliceCategory) or int(category.base_object().rank()) != 1:
        raise ValueError("content strata require a rank-one pointed coslice category")
    source = category.base_object()
    if int(source.q(next(iter(source.module_generators())))) != int(square):
        raise ValueError("the pointed source form must have the requested square")
    if content < 1:
        raise ValueError("nonzero pointed-vector content must be positive")
    if square != 0 or lattice.is_nondegenerate() or all(
        lattice.b(x, y) == 0 for x in lattice.module_generators()
        for y in lattice.module_generators()
    ):
        yield from pointed_content_classes(category, lattice, square, content, acting_group=acting_group)
        return
    _reduction, _n, _k, complement = _degenerate_quotient_and_radical(lattice)
    if complement.is_definite() or not complement.is_isotropic():
        yield from pointed_content_classes(category, lattice, square, content, acting_group=acting_group)
        return
    for marked in iter_pointed_isotropic_classes(category, lattice, acting_group=acting_group):
        if int(_marked_vector(marked).content()) == content:
            yield marked


def iter_primitive_isotropic_line_classes(lattice: Parent, *, acting_group=None):
    r"""Fair exact moving primitive-null-line orbit stream.

    Every primitive null vector occurs in the coordinate-height walk.
    A moving rank-one isotropic sublattice is identified with the pair
    ``{v,-v}``; exact pointed transport decides equivalence of those
    generators.  The stream is exhaustive without a finite-class claim.
    Supported finite-character arithmetic subgroups split each moving orbit.
    """
    if int(lattice.rank()) == 0:
        return
    if lattice.is_nondegenerate() and (lattice.is_definite() or int(lattice.witt_index()) == 0):
        return
    seen = []
    for vector in RECURSION.exact_vector_orbit_iterator(lattice, 0, primitive_only=True):
        if any(
            RECURSION.pointed_isometry(vector, previous) is not None
            or RECURSION.pointed_isometry(vector, -previous) is not None
            for previous in seen
        ):
            continue
        seen.append(vector)
        source = lattice.lattice_category()([[0]])
        marked = ARROWS.object(source.hom([vector], codomain=lattice))
        if acting_group is None:
            yield marked
        else:
            yield from lattice.O().split_orbit(marked, acting_group)


def iter_isotropic_line_content_classes(lattice: Parent, content: int, *, acting_group=None):
    r"""Exact exhaustive moving-null-line stream at fixed saturation index.

    Scaling by positive ``content`` is an equivariant bijection between
    primitive isotropic lines and rank-one isotropic sublattices with that
    saturation index. No finite orbit count is asserted.
    """
    if content < 1:
        raise ValueError("rank-one sublattice index in saturation must be positive")
    source = lattice.lattice_category()([[0]])
    generator = next(iter(source.module_generators()))
    for marking in iter_primitive_isotropic_line_classes(lattice, acting_group=acting_group):
        yield ARROWS.object(source.hom([
            lattice.scalar_multiple(lattice.base_ring()(content), marking.arrow()(generator))
        ], codomain=lattice))


def iter_all_isotropic_line_classes(lattice: Parent, *, acting_group=None):
    r"""Fair moving-null-line quotient, including all positive contents.

    Walk complete finite coordinate-height shells and decide exact full
    O(L)-equivalence for every new primitive null line. Split each new
    full-group class by the supplied subgroup, then emit its content-c
    multiple at height ``discovery_height+c-1``. Each shell is finite,
    so every primitive orbit and every positive content is reached, provided
    the exact pointed-equivalence backend terminates on each comparison.
    Unlike an iterator that waits for the next *new* primitive class, this
    remains productive when the primitive orbit quotient is finite but an
    unbounded search cannot certify that no further classes exist.
    """
    from itertools import count
    from sage_indefinite_port.indefinite.recursive import _coordinate_box, _coordinate_element

    if int(lattice.rank()) == 0:
        return
    if lattice.is_nondegenerate() and int(lattice.witt_index()) == 0:
        return
    full_representatives = []
    representatives = []
    source = lattice.lattice_category()([[0]])
    generator = next(iter(source.module_generators()))
    for height in count(1):
        for coordinates in _coordinate_box(lattice, height):
            point = _coordinate_element(lattice, coordinates)
            if point.q() != lattice.base_ring().zero() or int(point.content()) != 1:
                continue
            if any(RECURSION.pointed_isometry(point, prior) is not None or
                   RECURSION.pointed_isometry(point, -prior) is not None
                   for prior in full_representatives):
                continue
            full_representatives.append(point)
            marked = ARROWS.object(source.hom([point], codomain=lattice))
            split = (marked,) if acting_group is None else lattice.O().split_orbit(marked, acting_group)
            representatives.extend((height, orbit) for orbit in split)
        for discovery_height, marked in representatives:
            content = height - discovery_height + 1
            point = marked.arrow()(generator)
            yield ARROWS.object(source.hom([
                lattice.scalar_multiple(lattice.base_ring()(content), point)
            ], codomain=lattice))


def canonical_moving_isotropic_line(vector, *, acting_group=None):
    r"""Canonical moving isotropic line and a certified ambient transporter.

    The moving line generated by v is the union of the two pointed
    O(L)-orbits of v and -v. Their independently exact height-lexicographic
    minima therefore have a minimum equal to the canonical moving line.
    The returned isometry takes its canonical generator to v or -v.
    With an explicit arithmetic subgroup, search the same finite height
    segment using the exact marked transporter in that subgroup; this
    retains setwise source-sign freedom without enlarging the action.
    """
    lattice = vector.parent()
    if vector == lattice.zero() or vector.q() != lattice.base_ring().zero():
        raise ValueError("a moving isotropic line requires a nonzero null generator")
    if acting_group is not None:
        from sage_indefinite_port.indefinite.recursive import _coordinate_box, _coordinate_element, _coordinate_height

        source = lattice.lattice_category()([[0]])
        target = ARROWS.object(source.hom([vector], codomain=lattice))
        for height in range(_coordinate_height(lattice, vector) + 1):
            for coordinates in _coordinate_box(lattice, height):
                candidate = _coordinate_element(lattice, coordinates)
                if candidate.q() != lattice.base_ring().zero() or candidate.content() != vector.content():
                    continue
                marked = ARROWS.object(source.hom([candidate], codomain=lattice))
                witness = ARROWS.Core().Mor(marked, target).an_element_in(acting_group)
                if witness is not None:
                    return marked, witness.forward().right()
        raise ArithmeticError("the supplied subgroup did not identify an isotropic line with itself")
    forward = RECURSION.marked_orbit_canonical(vector)
    backward = RECURSION.marked_orbit_canonical(-vector)
    labels = tuple(lattice.module_generating_set())

    def key(result):
        coordinates = tuple(int(result.representative.to_vector()(label)) for label in labels)
        return max(abs(entry) for entry in coordinates), coordinates

    selected = min((forward, backward), key=key)
    source = lattice.lattice_category()([[0]])
    marked = ARROWS.object(source.hom([selected.representative], codomain=lattice))
    return marked, selected.transporter


def nonprimitive_isotropic_line_classes(lattice: Parent, content: int, *, acting_group=None) -> Parent:
    r"""Moving rank-one isotropic sublattice orbits of prescribed index in saturation.

    A rank-one sublattice of index ``content`` in its saturation is
    ``content * S`` for a unique primitive isotropic line ``S``.  The
    operation commutes with integral isometries, so the primitive moving
    line classes determine these orbits without a second classification.
    An optional arithmetic subgroup splits the *moving* line orbits and
    retains source sign changes.
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
    if acting_group is not None:
        classes = tuple(split for marking in classes
                        for split in lattice.O().split_orbit(marking, acting_group))
    return finite_ordered_set(tuple(classes))



def nonisotropic_line_classes(lattice: Parent, square: int, *, acting_group=None) -> Parent:
    r"""Moving rank-one sublattice orbits of nonzero generator square.

    First enumerate every pointed vector orbit, including imprimitive
    contents.  Passing from pointed vectors to rank-one sublattices
    identifies exactly the orbits related by sign, as the source
    automorphism group is ``{+1,-1}``.  Decide that identification with
    the existing exact moving-arrow transporter, not with a heuristic
    comparison of norms or divisibilities.  An ``acting_group`` splits
    the resulting *line* orbits, retaining the source sign equivalence.
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
    if acting_group is not None:
        # Split *moving* line orbits, not pointed vector orbits: source
        # negation remains an allowed automorphism even when -id is not
        # contained in the selected arithmetic subgroup.
        kept = [split for moving in kept for split in lattice.O().split_orbit(moving, acting_group)]
    return finite_ordered_set(tuple(kept))



def homothetic_isotropic_sublattice_classes(lattice: Parent, rank: int, multiplier: int, *, acting_group=None) -> Parent:
    r"""Moving isotropic sublattice classes ``J=multiplier*S`` for primitive rank-``rank`` ``S``.

    For a free rank-``rank`` saturation, ``S/J`` has elementary divisors
    ``(multiplier,)*rank``.  The map ``S -> multiplier*S`` is equivariant
    under every integral isometry, so it induces a bijection of orbit sets.
    Other elementary-divisor types require the actual finite-module
    stabilizer action and are not silently conflated with this family.
    An optional subgroup splits moving orbits without fixing source bases.
    """
    if rank < 1 or multiplier < 1:
        raise ValueError("rank and saturation multiplier must be positive")
    if rank > int(lattice.rank()):
        return finite_ordered_set(())
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
    if acting_group is not None:
        classes = [split for marked in classes for split in lattice.O().split_orbit(marked, acting_group)]
    return finite_ordered_set(tuple(classes))




def _split_radical_saturation_stabilizer_generators(saturated: Parent):
    r"""Exact stabilizer of a primitive split-radical graph over all of C.

    After integral adaptation, ``S=C+R'`` for a primitive rank-r direct
    summand ``R'`` of ``R``.  The stabilizer is ``O(C)`` times the parabolic
    ``Stab_GL(R)(R')``, extended by ``Hom(C,R')``.  Complete generators
    of both general linear blocks and all integral shears give the full
    group, not merely its finite image.
    """
    from sage.matrix.constructor import identity_matrix, matrix
    from sage.rings.integer_ring import ZZ as SageZZ
    from sage_extensions.categories.group.groups import OwnedGroups
    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation
    from sage_indefinite_port.indefinite.recursive import (
        _gram_in_label_order, _isometry_from_row_matrix, _conjugated,
    )

    lattice = saturated.inclusion().codomain()
    reduction = indefinite_reduced_presentation(lattice)
    reduced = reduction.reduced
    gram = _gram_in_label_order(reduced)
    n, k = int(gram.rank()), int(gram.nrows() - gram.rank())
    r = int(saturated.inclusion().domain().rank()) - n
    if not n or not (1 <= r <= k):
        raise ValueError("the split-radical stabilizer requires the full nondegenerate quotient plus a radical summand")
    labels = tuple(reduced.module_generating_set())
    inverse = reduction.isometry.inverse()
    actual = matrix(SageZZ, [[int(inverse(saturated.inclusion()(x)).to_vector()(label)) for label in labels]
        for x in saturated.module_generators()]).hermite_form()
    expected = matrix(SageZZ, [
        [int(i == j) for j in range(n + k)] for i in range(n)
    ] + [[int(j == n + i) for j in range(n + k)] for i in range(r)]).hermite_form()
    size = n + k
    # A primitive graph over all of C with one primitive radical line
    # is conjugate to the standard C+Zz_1 by an integral radical change
    # and an unrestricted shear C -> R.  Construct that conjugation,
    # rather than requiring the caller's basis to be adapted.
    if actual.nrows() != n + r or actual[:n, :n] != identity_matrix(SageZZ, n) or any(actual[n + i, j] for i in range(r) for j in range(n)):
        raise ValueError("the saturation must project unimodularly onto the nondegenerate quotient")
    radical_rows = matrix(SageZZ, [list(actual.row(n + i)[n:]) for i in range(r)])
    _, _, d = radical_rows.smith_form()
    if abs(int(d.det())) != 1:
        raise ArithmeticError("the radical summand must have an integral unimodular completion")
    adapted = identity_matrix(SageZZ, size)
    for i in range(k):
        for j in range(k):
            adapted[n + i, n + j] = d[i, j]
    for i in range(n):
        for j in range(k):
            adapted[i, n + j] = -sum(int(actual[i, n + t]) * int(d[t, j]) for t in range(k))
    if (actual * adapted).hermite_form() != expected:
        raise ValueError("the saturation radical summand is not primitive or the graph is not integral")
    transport = _isometry_from_row_matrix(reduced, reduced, adapted)
    blocks = []
    complement = reduced.lattice_category()([[int(gram[i, j]) for j in range(n)] for i in range(n)])
    clabels = tuple(complement.module_generating_set())
    for generator in RECURSION.orthogonal_group_generators(complement):
        block = identity_matrix(SageZZ, size)
        for i, label in enumerate(clabels):
            row = generator(complement.module_generator(label)).to_vector()
            for j, target_label in enumerate(clabels):
                block[i, j] = int(row(target_label))
        blocks.append(block)
    def embedded_gl(offset, dimension):
        if not dimension:
            return
        group = OwnedGroups().GL(dimension, ZZ)
        for generator in group.group_generators():
            action = group._to_engine(generator).matrix()
            block = identity_matrix(SageZZ, size)
            for i in range(dimension):
                for j in range(dimension):
                    block[offset + i, offset + j] = int(action[i, j])
            blocks.append(block)

    embedded_gl(n, r)
    embedded_gl(n + r, k - r)
    for i in range(r, k):
        for j in range(r):
            shear = identity_matrix(SageZZ, size)
            shear[n + i, n + j] = 1
            blocks.append(shear)
    for i in range(n):
        for j in range(r):
            shear = identity_matrix(SageZZ, size)
            shear[i, n + j] = 1
            blocks.append(shear)
    reduced_generators = tuple(
        transport.inverse() * _isometry_from_row_matrix(reduced, reduced, block) * transport
        for block in blocks
    )
    generators = _conjugated(reduction.isometry, reduced_generators)
    if any(not saturated.inclusion().is_in_image(g(saturated.inclusion()(x)))
           for g in generators for x in saturated.module_generators()):
        raise ArithmeticError("a claimed split-radical stabilizer generator does not preserve the saturation")
    return generators


def nonprimitive_sublattice_classes(saturated: Parent | Morphism, divisors: tuple[int, ...], *, stabilizer_group=None) -> tuple[Morphism, ...]:
    r"""Finite-index sublattice classes under the actual ambient stabilizer image.

    The port supplies the integral action on the primitive saturation.  The
    generic finite-index submodule orbit construction belongs to the module
    owner and returns represented subobjects, not matrix surrogates.  When
    ``stabilizer_group`` is supplied, it is used without enlargement to the
    full ambient or abstract general linear group.
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
    if stabilizer_group is not None:
        # The supplied group is already the actual arithmetic stabilizer
        # of S (or a subgroup thereof).  Never enlarge it to O(S).
        ambient_generators = tuple(stabilizer_group.group_generators())
        if any(not inclusion.is_in_image(g(point)) for g in ambient_generators for point in basis):
            raise ValueError("the supplied group does not preserve the saturation")
        module = inclusion.domain()
        actions = tuple(module.Isom(module)(
            lambda label, generator=g: inclusion.lift(generator(inclusion(module.module_generator(label))))
        ) for g in ambient_generators)
    elif isotropic:
        if all(lattice.b(x, y) == 0 for x in lattice.module_generators() for y in lattice.module_generators()):
            linear = lattice.general_linear_group()
            preserving = (linear.group_generators()
                          if int(inclusion.domain().rank()) == int(lattice.rank())
                          else linear.flag_stabilizer((saturated,)).group_generators())
            ambient_generators = _zero_form_orthogonal_generators(lattice, preserving)
        elif not lattice.is_nondegenerate() and all(
            RECURSION.vector_lies_in_radical(point) for point in basis
        ):
            # The radical is preserved by every integral isometry and
            # carries the full integral GL-action on its primitive summands.
            # The nondegenerate perpendicular-lattice induction is not
            # defined: L -> S* vanishes identically when S lies in rad(L).
            ambient_generators = _degenerate_radical_flag_generators(
                lattice, (int(inclusion.domain().rank()),), basis
            )
        elif not lattice.is_nondegenerate():
            # A mixed degenerate maximal null saturation meets rad(L)
            # nontrivially. Witt extension for nondegenerate ambient spaces
            # cannot extend arbitrary GL(S) actions across that intersection.
            # The moving-arrow stabilizer owns the integral radical/glue
            # constraints and supplies precisely the permitted action on S.
            ambient_generators = tuple(automorphism_group_generators(
                ARROWS, ARROWS.object(inclusion)
            ))
        else:
            ambient_generators = INDUCTION.plane_automorphism_generators(saturated)
        module = inclusion.domain()
        actions = tuple(module.Isom(module)(
            lambda label, generator=g: inclusion.lift(generator(inclusion(module.module_generator(label))))
        ) for g in ambient_generators)
    else:
        if not inclusion.domain().is_nondegenerate():
            ambient_generators = _split_radical_saturation_stabilizer_generators(saturated)
            module = inclusion.domain()
            actions = tuple(module.Isom(module)(
                lambda label, generator=g: inclusion.lift(generator(inclusion(module.module_generator(label))))
            ) for g in ambient_generators)
            submodules = module.finite_index_submodule_orbits(divisors, actions)
            return tuple(lattice.subobject_on(tuple(inclusion(subobject(point))
                for point in subobject.domain().module_generators())).inclusion() for subobject in submodules)
        # The actual stabilizer image can be smaller than O(S) when S and
        # S^perp have nontrivial glue.  Its generators come from the full
        # finite discriminant-image preimage, not a generator filter.
        ambient_generators = _nondegenerate_stabilizer_generators(inclusion, based=False)
        module = inclusion.domain()
        actions = tuple(module.Isom(module)(
            lambda label, generator=g: inclusion.lift(generator(inclusion(module.module_generator(label))))
        ) for g in ambient_generators)
    submodules = inclusion.domain().finite_index_submodule_orbits(divisors, actions)
    return tuple(lattice.subobject_on(tuple(inclusion(subobject(point))
        for point in subobject.domain().module_generators())).inclusion() for subobject in submodules)



def nonprimitive_isotropic_sublattice_classes(lattice: Parent, rank: int, divisors: tuple[int, ...], *, acting_group=None) -> Parent:
    r"""Moving totally isotropic sublattices with prescribed saturation Smith type.

    Enumerate primitive saturated sublattice orbits, split them under the
    *actual* requested arithmetic group, and classify finite-index
    submodules under the stabilizer of each saturation.  Thus the finite
    quotient step does not replace the ambient arithmetic group by GL(S).
    """
    if rank < 1 or len(divisors) != rank:
        raise ValueError("positive isotropic rank and matching elementary-divisor tuple required")
    if any(d < 1 for d in divisors) or any(right % left for left, right in pairwise(divisors)):
        raise ValueError("elementary divisors must be positive and ordered by divisibility")
    if acting_group is lattice.O():
        acting_group = None
    if rank > int(lattice.rank()):
        return finite_ordered_set(())
    labels = tuple(lattice.module_generating_set())
    if all(lattice.b(x, y) == 0 for x in lattice.module_generators() for y in lattice.module_generators()):
        # A zero-form lattice has O(L)=GL(L), transitive on primitive
        # rank-r direct summands. Its saturation stabilizer induces the
        # *full* GL_r(ZZ) on the summand. This avoids calling the
        # nondegenerate Witt-index recursion on a degenerate lattice.
        saturated = lattice.primitive_isotropic_subobject(*(lattice.module_generator(label) for label in labels[:rank]))
        inclusion = saturated.inclusion()
        module = inclusion.domain()
        if acting_group is None:
            submodules = module.finite_index_submodule_orbits(divisors, module.general_linear_group().group_generators())
            return finite_ordered_set(tuple(ARROWS.object(lattice.subobject_on(tuple(
                inclusion(subobject(point)) for point in subobject.domain().module_generators()
            )).inclusion()) for subobject in submodules))
        source = lattice.lattice_category()([[0] * rank for _ in range(rank)])
        standard = ARROWS.object(source.hom(tuple(
            inclusion(point) for point in module.module_generators()
        ), codomain=lattice))
        # A full-rank primitive sublattice is the whole lattice: its
        # moving orbit is a singleton under every arithmetic subgroup.
        primitive = (standard,) if rank == int(lattice.rank()) else lattice.O().split_orbit(standard, acting_group)
        classes = []
        for marked in primitive:
            arrow = marked.arrow()
            saturation = lattice.subobject_on(tuple(arrow(point) for point in arrow.domain().module_generators()))
            stabilizer = acting_group if rank == int(lattice.rank()) else acting_group.intersection(lattice.O().stabilizer(saturation))
            classes.extend(ARROWS.object(embedding) for embedding in nonprimitive_sublattice_classes(
                saturation, divisors, stabilizer_group=stabilizer
            ))
        return finite_ordered_set(tuple(classes))
    source = lattice.lattice_category()([[0] * rank for _ in range(rank)])
    primitive = primitive_marking_classes(ARROWS, source, lattice)
    group = lattice.O() if acting_group is None else acting_group
    if acting_group is not None:
        primitive = tuple(split for marked in primitive for split in lattice.O().split_orbit(marked, group))
    classes = []
    for marked in primitive:
        arrow = marked.arrow()
        saturated = lattice.subobject_on(tuple(arrow(point) for point in arrow.domain().module_generators()))
        if acting_group is None:
            submodules = nonprimitive_sublattice_classes(saturated, divisors)
        else:
            stabilizer = group.intersection(lattice.O().stabilizer(saturated))
            submodules = nonprimitive_sublattice_classes(saturated, divisors, stabilizer_group=stabilizer)
        classes.extend(ARROWS.object(embedding) for embedding in submodules)
    return finite_ordered_set(tuple(classes))


def iter_all_isotropic_sublattice_classes(lattice: Parent, rank: int, *, acting_group=None):
    r"""Enumerate all moving isotropic rank-r sublattices by Smith type.

    The elementary divisors of a finite-index sublattice of its primitive
    saturation satisfy d_1 | ... | d_r. Every such tuple occurs in the
    finite shell with d_r equal to its maximum. The fixed-Smith-type
    classification supplies the exact stabilizer orbits, so these disjoint
    shells exhaust the unrestricted moving quotient without an index bound,
    in the domain where the fixed-Smith-type classifier is complete and
    terminating. No finite cardinality is claimed for the unrestricted set.
    """
    from itertools import count, product

    if rank < 1:
        raise ValueError("isotropic sublattice rank must be positive")
    if rank > int(lattice.rank()):
        return
    positive, negative = (int(index) for index in lattice.signature_pair())
    radical_rank = int(lattice.rank()) - positive - negative
    if rank > radical_rank + min(positive, negative):
        # Every totally isotropic rational subspace has dimension at most
        # dim(radical) + min(p,q), even when the rational Witt index of the
        # nondegenerate quotient is smaller than its real Witt index.
        return
    if not lattice.is_nondegenerate():
        # Projection modulo the radical sends every rational totally
        # isotropic subspace to one in the nondegenerate quotient. Its
        # dimension is therefore at most dim(rad L) + witt(L/rad L).
        # The real signature bound above is not sufficient: an indefinite
        # quotient can be rationally anisotropic, in which case an empty
        # rank-r stream must terminate rather than search every Smith shell.
        _reduction, _n, k, complement = _degenerate_quotient_and_radical(lattice)
        quotient_witt = int(complement.witt_index()) if complement is not None and complement.is_isotropic() else 0
        if rank > k + quotient_witt:
            return
    if rank == 1:
        yield from iter_all_isotropic_line_classes(lattice, acting_group=acting_group)
        return
    if lattice.is_nondegenerate() and rank > int(lattice.witt_index()):
        return
    for maximum in count(1):
        for divisors in product(range(1, maximum + 1), repeat=rank - 1):
            smith = (*divisors, maximum)
            if any(right % left for left, right in pairwise(smith)):
                continue
            yield from nonprimitive_isotropic_sublattice_classes(
                lattice, rank, smith, acting_group=acting_group
            )



def nonprimitive_nondegenerate_sublattice_classes(lattice: Parent, saturation_type: Parent, divisors: tuple[int, ...], *, acting_group=None) -> Parent:
    r"""Moving finite-index sublattices with a prescribed primitive nondegenerate saturation type.

    First classify primitive images of ``saturation_type`` up to source
    isometry, then take the actual arithmetic stabilizer action on the
    finite-index submodules.  No genus-isometry or full O(S) shortcut is
    used in place of the embedded saturation's ambient stabilizer.
    """
    if not saturation_type.is_nondegenerate():
        raise ValueError("use isotropic or radical marking classification for degenerate saturation")
    if len(divisors) != int(saturation_type.rank()):
        raise ValueError("elementary divisors must match the saturation rank")
    if any(d < 1 for d in divisors) or any(right % left for left, right in pairwise(divisors)):
        raise ValueError("elementary divisors must be positive and ordered by divisibility")
    primitive = primitive_marking_classes(ARROWS, saturation_type, lattice)
    group = lattice.O() if acting_group is None else acting_group
    if acting_group is not None:
        primitive = tuple(split for marked in primitive for split in lattice.O().split_orbit(marked, group))
    classes = []
    for marked in primitive:
        arrow = marked.arrow()
        saturated = lattice.subobject_on(tuple(arrow(x) for x in arrow.domain().module_generators()))
        stabilizer = None if acting_group is None else group.intersection(lattice.O().stabilizer(saturated))
        classes.extend(ARROWS.object(embedding) for embedding in
            nonprimitive_sublattice_classes(saturated, divisors, stabilizer_group=stabilizer))
    return finite_ordered_set(tuple(classes))


def iter_all_nondegenerate_sublattice_index_classes(lattice: Parent, saturation_type: Parent, *, acting_group=None):
    r"""Unbounded Smith-type quotient for a fixed primitive saturation type.

    Every finite-index sublattice has a unique saturation and a unique
    elementary-divisor tuple d_1 | ... | d_r. Enumerate these tuples by
    maximum divisor and retain the actual arithmetic saturation stabilizer
    for each stratum. Completeness is relative to the fixed-saturation
    primitive embedding classifier, as in the finite-index owner.
    """
    from itertools import count, product

    if not saturation_type.is_nondegenerate():
        raise ValueError("a nondegenerate saturation type is required")
    rank = int(saturation_type.rank())
    if rank < 1:
        raise ValueError("a positive-rank saturation type is required")
    if rank > int(lattice.rank()):
        return
    for maximum in count(1):
        for divisors in product(range(1, maximum + 1), repeat=rank - 1):
            smith = (*divisors, maximum)
            if any(right % left for left, right in pairwise(smith)):
                continue
            yield from nonprimitive_nondegenerate_sublattice_classes(
                lattice, saturation_type, smith, acting_group=acting_group
            )


def primitive_marking_classes(category: Category, source: Parent | Morphism, lattice: Parent, *, acting_group=None) -> Parent:
    r"""One marked lattice for each isomorphism class of primitive markings of ``lattice`` of shape ``source``.

    ``category`` is the category of the marking: the coslice under ``<n>`` for pointed
    lattices, ``Ar(Lattices)`` for sublattices, ``Fun([t], Lattices)`` for finite flags. ``source``
    is the shape: the lattice ``<n>``, a totally isotropic lattice of rank ``k`` (zero form),
    or the chain of such lattices for a flag. A flag source may be a shorter
    diagram of its isotropic terms or a full marked diagram in ``category``;
    in the latter case its terminal arrow into the ambient lattice does not
    contribute another isotropic rank. ``acting_group`` splits
    the full-group orbit set using the actual supplied arithmetic subgroup
    without changing the source marking category.
    """
    match category:
        case CosliceCategory():
            if int(source.rank()) == 0:
                classes = (category.object(source.hom([], codomain=lattice)),)
            elif int(source.rank()) > int(lattice.rank()):
                classes = ()
            elif (
                all(lattice.b(x, y) == 0 for x in lattice.module_generators() for y in lattice.module_generators())
                and all(source.b(x, y) == 0 for x in source.module_generators() for y in source.module_generators())
            ):
                # GL_n(ZZ) is transitive on ordered primitive r-frames:
                # every such frame extends to a unimodular basis. Unlike
                # the moving arrow, this based coslice fixes the source.
                labels = tuple(lattice.module_generating_set())
                classes = (category.object(source.hom(tuple(
                    lattice.module_generator(label) for label in labels[:int(source.rank())]
                ), codomain=lattice)),)
            elif int(source.rank()) == 1:
                square = int(source.q(next(iter(source.module_generators()))))
                # The nondegenerate Witt-index dispatcher is not a
                # classifier for degenerate ambient lattices.  At rank one
                # primitive markings are precisely the content-one stratum;
                # the radical/glue classifier handles its finite cases and
                # explicitly refuses infinite mixed-null strata.
                classes = (tuple(pointed_content_classes(category, lattice, square, 1))
                           if not lattice.is_nondegenerate()
                           else _pointed_classes(category, lattice, square))
            else:
                source_basis = tuple(source.module_generators())
                if source.is_nondegenerate():
                    # The owned primitive-embedding classification='emb'
                    # fixes the source marking pointwise.  Classification
                    # 'sub' would silently quotient by O(source) instead.
                    classes = tuple(category.object(embedding) for embedding in
                        source.Emb(lattice).primitive_embedding_class_representatives(classification="emb"))
                elif not lattice.is_nondegenerate() and all(
                    source.b(x, y) == 0 for x in source.module_generators()
                    for y in source.module_generators()
                ):
                    reduction, n, k, complement = _degenerate_quotient_and_radical(lattice)
                    reduced = reduction.reduced
                    quotient_witt = int(complement.witt_index()) if complement is not None and complement.is_isotropic() else 0
                    if int(source.rank()) > k + quotient_witt:
                        classes = ()
                    elif quotient_witt:
                        raise ValueError(
                            "mixed degenerate primitive isotropic based-marking orbit sets may be infinite; "
                            "this finite-class interface cannot represent them"
                        )
                    else:
                        basis = tuple(reduction.isometry(reduced.module_generator(label))
                                      for label in tuple(reduced.module_generating_set())[n:n + int(source.rank())])
                        classes = (category.object(source.hom(basis, codomain=lattice)),)
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
            elif rank > int(lattice.rank()):
                classes = ()
            elif (
                rank <= int(lattice.rank())
                and all(lattice.b(x, y) == 0 for x in lattice.module_generators() for y in lattice.module_generators())
                and all(source.b(x, y) == 0 for x in source.module_generators() for y in source.module_generators())
            ):
                # Primitive rank-r direct summands of a zero-form lattice
                # form one GL_n(ZZ)-orbit, by completion of a primitive
                # frame to an integral basis.
                labels = tuple(lattice.module_generating_set())
                classes = (ARROWS.object(source.hom(tuple(
                    lattice.module_generator(label) for label in labels[:rank]
                ), codomain=lattice)),)
            elif rank == 1 and not lattice.is_nondegenerate() and all(
                source.b(x, y) == 0 for x in source.module_generators()
                for y in source.module_generators()
            ):
                # The nondegenerate Witt-index flag induction is undefined
                # for a lattice with radical.  Classify primitive pointed
                # null vectors by exact content, then forget their signs to
                # obtain moving lines.  In the anisotropic-quotient case
                # there is exactly one radical line; an isotropic quotient
                # can have infinitely many and the finite API refuses it.
                pointed_category = CosliceCategory(lattice.lattice_category(), source)
                try:
                    pointed = pointed_content_classes(pointed_category, lattice, 0, 1)
                except ValueError as exc:
                    if "may be infinite" not in str(exc):
                        raise
                    raise ValueError(
                        "mixed degenerate primitive null-line orbit sets may be infinite; "
                        "use iter_primitive_isotropic_line_classes"
                    ) from exc
                classes = tuple(ARROWS.object(lattice.primitive_isotropic_subobject(
                    _marked_vector(marking)
                ).inclusion()) for marking in pointed)
            elif not lattice.is_nondegenerate() and all(
                source.b(x, y) == 0 for x in source.module_generators()
                for y in source.module_generators()
            ):
                reduction, n, k, complement = _degenerate_quotient_and_radical(lattice)
                reduced = reduction.reduced
                quotient_witt = int(complement.witt_index()) if complement is not None and complement.is_isotropic() else 0
                if rank > k + quotient_witt:
                    classes = ()
                elif quotient_witt == 1 and rank == k + 1:
                    # A maximal totally isotropic rational subspace of
                    # C_Q + R_Q contains R_Q when C has Witt index one.
                    # Primitivity then forces the whole integral radical
                    # into the sublattice.  Modulo R its image is exactly a
                    # primitive isotropic line in C; the radical shears
                    # remove all choices of lift.  Thus moving maximal
                    # null sublattices are classified by the quotient's
                    # primitive moving null lines, with no infinite
                    # mixed-content ambiguity.
                    line_source = lattice.lattice_category()([[0]])
                    quotient_lines = primitive_marking_classes(ARROWS, line_source, complement)
                    reduced_labels = tuple(reduced.module_generating_set())
                    complement_labels = tuple(complement.module_generating_set())
                    radical_basis = tuple(reduction.isometry(reduced.module_generator(label))
                                          for label in reduced_labels[n:])
                    classes = []
                    for line in quotient_lines:
                        arrow = line.arrow()
                        point = arrow(next(iter(arrow.domain().module_generators())))
                        coordinates = point.to_vector()
                        lift = reduced.linear_combination({
                            reduced_labels[i]: reduced.base_ring()(int(coordinates(label)))
                            for i, label in enumerate(complement_labels)
                        })
                        basis = (reduction.isometry(lift), *radical_basis)
                        classes.append(ARROWS.object(lattice.primitive_isotropic_subobject(*basis).inclusion()))
                    classes = tuple(classes)
                elif quotient_witt:
                    raise ValueError(
                        "mixed degenerate primitive isotropic sublattice orbit sets may be infinite; "
                        "use iter_all_isotropic_sublattice_classes"
                    )
                # The anisotropic quotient forces every isotropic subspace
                # into rad(L).  Since rad(L) is saturated, its primitive
                # rank-r sublattices form one GL(rad L)-orbit for r <= k.
                else:
                    basis = tuple(reduction.isometry(reduced.module_generator(label))
                                  for label in tuple(reduced.module_generating_set())[n:n + rank])
                    classes = (ARROWS.object(lattice.primitive_isotropic_subobject(*basis).inclusion()),)
            elif source.is_nondegenerate():
                # Moving sublattices allow the full source isometry group.
                # The owned 'sub' classification quotients by both O(source)
                # and O(target), unlike the based 'emb' classification.
                classes = tuple(ARROWS.object(embedding) for embedding in
                    source.Emb(lattice).primitive_embedding_class_representatives(classification="sub"))
            else:
                source_basis = tuple(source.module_generators())
                if any(source.b(left, right) != 0 for left in source_basis for right in source_basis):
                    raise ValueError("mixed degenerate Gram sublattices require their radical/glue classification")
                classes = tuple(ARROWS.object(lattice.primitive_isotropic_subobject(*basis).inclusion()) for basis in INDUCTION.isotropic_flag_bases(lattice, (rank,)))
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            if isinstance(source, Morphism):
                shape_arrows = (source,)
            else:
                shape_arrows = _flag_arrows(source)
            if not shape_arrows:
                raise ValueError("an isotropic flag shape requires at least one inclusion")
            full_diagram_shape = (
                not isinstance(source, Morphism)
                and len(shape_arrows) == int(category.domain_category().object_set().cardinality()) - 1
                and (
                    shape_arrows[-1].codomain() is lattice
                    or not all(
                        shape_arrows[-1].codomain().b(left, right) == 0
                        for left in shape_arrows[-1].codomain().module_generators()
                        for right in shape_arrows[-1].codomain().module_generators()
                    )
                )
            )
            last_shape_lattice = shape_arrows[-1].domain() if full_diagram_shape else shape_arrows[-1].codomain()
            if any(
                arrow.domain().b(left, right) != 0
                for arrow in shape_arrows
                for left in arrow.domain().module_generators()
                for right in arrow.domain().module_generators()
            ) or any(
                last_shape_lattice.b(left, right) != 0
                for left in last_shape_lattice.module_generators()
                for right in last_shape_lattice.module_generators()
            ):
                raise ValueError("nonisotropic flag shapes require a general formed-diagram classification")
            if full_diagram_shape:
                # A full marked diagram includes its final arrow into the
                # ambient lattice; this arrow is not another flag rank.
                dimensions = tuple(int(arrow.domain().rank()) for arrow in shape_arrows)
            else:
                # A shorter source diagram is precisely the source flag
                # shape, with its last codomain the final isotropic term.
                dimensions = (int(shape_arrows[0].domain().rank()), *(int(arrow.codomain().rank()) for arrow in shape_arrows))
            if len(dimensions) != int(category.domain_category().object_set().cardinality()) - 1:
                raise ValueError("the isotropic source shape must specify exactly the target flag's ranks")
            if dimensions[0] < 1 or any(left >= right for left, right in pairwise(dimensions)):
                raise ValueError("isotropic flag ranks must be positive and strictly increasing")
            if dimensions[-1] > int(lattice.rank()):
                classes = ()
            elif all(lattice.b(x, y) == 0 for x in lattice.module_generators() for y in lattice.module_generators()):
                # GL_n(ZZ) acts transitively on primitive integral flags
                # with any fixed strictly increasing rank sequence.
                if dimensions[0] < 1 or any(left >= right for left, right in pairwise(dimensions)):
                    raise ValueError("isotropic flag ranks must be positive and strictly increasing")
                basis = tuple(lattice.module_generators())[:dimensions[-1]]
                classes = (category(_flag_marking_from_basis(lattice, basis, dimensions)),)
            elif not lattice.is_nondegenerate():
                reduction, n, k, complement = _degenerate_quotient_and_radical(lattice)
                reduced = reduction.reduced
                quotient_witt = int(complement.witt_index()) if complement is not None and complement.is_isotropic() else 0
                if dimensions[0] < 1 or any(left >= right for left, right in pairwise(dimensions)):
                    raise ValueError("isotropic flag ranks must be positive and strictly increasing")
                if dimensions[-1] > k + quotient_witt:
                    classes = ()
                elif quotient_witt:
                    raise ValueError(
                        "mixed degenerate primitive isotropic flag orbit sets may be infinite; "
                        "this finite-class interface cannot represent them"
                    )
                else:
                    basis = tuple(reduction.isometry(reduced.module_generator(label))
                                  for label in tuple(reduced.module_generating_set())[n:n + dimensions[-1]])
                    classes = (category(_flag_marking_from_basis(lattice, basis, dimensions)),)
            else:
                classes = tuple(category(_flag_marking_from_basis(lattice, basis, dimensions)) for basis in INDUCTION.isotropic_flag_bases(lattice, dimensions))
        case _:
            raise ValueError(f"{category} is not a category of marked lattices this port classifies: a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t]")
    if acting_group is not None and acting_group is not lattice.O():
        classes = tuple(split for marked in classes
                        for split in lattice.O().split_orbit(marked, acting_group))
    return finite_ordered_set(classes)


def _flag_marking_from_basis(lattice: Lattices.ParentMethods, basis: tuple[Lattices.ElementMethods, ...], dimensions: tuple[int, ...]) -> tuple[Morphism, ...]:
    r"""The chain of primitive isotropic prefix sublattices of an adapted basis."""
    if not dimensions or dimensions[0] < 1 or dimensions[-1] != len(basis) or any(
        left >= right for left, right in zip(dimensions, dimensions[1:])
    ):
        raise ValueError("flag ranks must be positive, strictly increasing, and end at the basis rank")
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
    "iter_pointed_content_classes",
    "iter_isotropic_line_content_classes",
    "iter_all_isotropic_line_classes",
    "canonical_moving_isotropic_line",
    "nonprimitive_isotropic_line_classes",
    "nonisotropic_line_classes",
    "homothetic_isotropic_sublattice_classes",
    "nonprimitive_sublattice_classes",
    "nonprimitive_isotropic_sublattice_classes",
    "iter_all_isotropic_sublattice_classes",
    "nonprimitive_nondegenerate_sublattice_classes",
    "iter_all_nondegenerate_sublattice_index_classes",
]
