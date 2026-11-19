r"""Isomorphisms and automorphisms of marked indefinite lattices.

A marked lattice is a diagram of lattice morphisms ending at a lattice L
(references/port-design.md, "Marked lattices"):

- a plain lattice L, an object of ``Lattices``;
- a pointed lattice ``<n> -> L``, an object of the coslice ``<n>/Lattices``;
- a primitive totally isotropic sublattice ``I -> L`` of any attainable rank, an object of ``Ar(Lattices)``;
- a partial or complete isotropic flag ``I_1 -> ... -> I_t -> L``, an object of ``Fun([t], Lattices)``;
- such a flag with an orientation of each graded piece, an object of ``OrientedIsotropicFlags(Lattices, t)``.

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
from sage_extensions.categories.isotropic_orbits import (
    FlagType,
    _graded_orientation_signs,
    _marked_isotropic_category,
    _marked_isotropic_object,
    _orientation_group_elements,
    _oriented_isotropic_stabilizer_generators,
)
from sage_extensions.categories.lattice_morphisms import LatticeIsometryMethods
from sage_extensions.categories.lattices import OrientedIsotropicFlags
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


def _zero_form_orthogonal_generators(ambient, linear_generators):
    r"""Carry integral GL generators into the owned O(L) homset for q=0."""
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ
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
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ

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
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ

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
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ

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
    from sage.matrix.constructor import matrix
    from sage.rings.rational_field import QQ as SageQQ

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
    from sage.matrix.constructor import matrix
    from sage.rings.rational_field import QQ as SageQQ
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
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ
    from sage.rings.rational_field import QQ as SageQQ

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
    from sage.matrix.constructor import matrix
    from sage.matrix.special import identity_matrix
    from sage.rings.integer_ring import ZZ as SageZZ
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
    # The terminal arrow's domain and the domain of its reconstructed
    # primitive image need not be the *same parent*. In particular a
    # subgroup-translated flag has original constituent arrows in a
    # transported top lattice, while _isotropic_sublattice builds its
    # ambient image anew. Translate each intermediate term through the
    # actual terminal embedding and its unique primitive-image lift.
    # This is the required parent-preserving identification, not equality
    # of two presentations of the same abstract zero-form lattice.
    top_arrow = arrows[-1]
    top_inclusion = space.inclusion()
    top_domain = top_inclusion.domain()
    terms = []
    composite = None
    for arrow in reversed(arrows[:-1]):
        composite = arrow if composite is None else composite * arrow
        original = _isotropic_sublattice(composite)
        original_inclusion = original.inclusion()
        lifted = tuple(
            top_inclusion.lift(top_arrow(original_inclusion(point)))
            for point in original_inclusion.domain().module_generators()
        )
        terms.append(top_domain.primitive_isotropic_subobject(*lifted))
    terms.reverse()
    dimensions = tuple(int(term.inclusion().domain().rank()) for term in terms) + (int(space.inclusion().domain().rank()),)
    if any(left >= right for left, right in zip(dimensions, dimensions[1:])):
        raise ValueError("an isotropic flag must have strictly increasing ranks")
    basis, images = INDUCTION.flag_action_data(space, tuple(terms))
    return dimensions, basis, images


def _oriented_flag_action(marked: Parent) -> tuple[
    tuple[int, ...], tuple[Lattices.ElementMethods, ...],
    tuple[tuple[Lattices.ElementMethods, ...], ...],
]:
    r"""The full orientation-preserving rational Levi image of a marked flag.

    CombinedAlgorithms.h:682-717 constructs the complete integral block
    parabolic P in GL(I). The orientation character of each successive
    graded factor gives a homomorphism P -> {+1,-1}^t. Its *kernel* is
    the required oriented Levi group. Filtering P's original generators
    by their signs would generally lose the kernel: a product of two
    disallowed generators may preserve every graded orientation.

    Enumerate the finite sign image and choose an exact word transversal
    T_s in P. Schreier's lemma gives all generators
    T_s g (T_(s*sign(g)))^-1 of the orientation kernel, not a
    selected subset of P. Express these products on the adapted integral
    basis of I and pass their actual images to the same rational
    stabilizer/lattice-intersection algorithm as for unoriented flags.
    """
    from collections import deque

    from sage.matrix.constructor import identity_matrix, matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    oriented = marked.functor().oriented_flag()
    dimensions, basis, images = _flag_action(marked)
    top = _isotropic_sublattice(_flag_arrows(marked)[-1]).inclusion()
    labels = tuple(top.domain().module_generating_set())
    rank = len(basis)

    def rows(points):
        return matrix(SageZZ, [
            [int(top.lift(point).to_vector()(label)) for label in labels]
            for point in points
        ])

    initial_rows = rows(basis)
    if abs(int(initial_rows.det())) != 1:
        raise ArithmeticError("an adapted oriented isotropic frame must be unimodular in its top sublattice")
    inverse = initial_rows.inverse().change_ring(SageZZ)
    matrices = tuple((rows(family) * inverse).change_ring(SageZZ) for family in images)
    if any(abs(int(transformation.det())) != 1 for transformation in matrices):
        raise ArithmeticError("the full flag Levi generators must be integral unimodular")

    block_starts = (0, *dimensions[:-1])

    def signs(transformation):
        result = tuple(int(transformation[left:right, left:right].det())
                       for left, right in zip(block_starts, dimensions, strict=True))
        if any(sign not in (-1, 1) for sign in result):
            raise ArithmeticError("a flag Levi generator has a nonunimodular graded block")
        return result

    generator_signs = tuple(signs(generator) for generator in matrices)
    identity_signs = (1,) * len(dimensions)
    transversal = {identity_signs: identity_matrix(SageZZ, rank)}
    queue = deque((identity_signs,))
    while queue:
        current = queue.popleft()
        representative = transversal[current]
        for generator, signature in zip(matrices, generator_signs, strict=True):
            next_signs = tuple(a*b for a,b in zip(current, signature, strict=True))
            if next_signs not in transversal:
                transversal[next_signs] = representative * generator
                queue.append(next_signs)

    kernel_matrices = []
    for current, representative in transversal.items():
        for generator, signature in zip(matrices, generator_signs, strict=True):
            next_signs = tuple(a*b for a,b in zip(current, signature, strict=True))
            schreier = (
                representative * generator
                * transversal[next_signs].inverse().change_ring(SageZZ)
            )
            if signs(schreier) != identity_signs:
                raise ArithmeticError("a Schreier generator fails to preserve the graded orientations")
            if schreier != identity_matrix(SageZZ, rank):
                kernel_matrices.append(schreier)

    ambient = oriented.lattice()
    ring = ambient.base_ring()

    def image_family(transformation):
        return tuple(
            sum((ambient.scalar_multiple(ring(int(value)), basis[j])
                 for j, value in enumerate(transformation.row(i)) if value),
                ambient.zero())
            for i in range(rank)
        )

    return dimensions, basis, tuple(image_family(generator) for generator in kernel_matrices)


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
            if (int(source_arrow.domain().rank()) == int(ambient.rank())
                    and int(target_arrow.domain().rank()) == int(target_arrow.codomain().rank())
                    and int(target_arrow.codomain().rank()) == int(ambient.rank())):
                # A full-rank based frame determines the rational ambient
                # isometry uniquely. Its integral existence is exactly the
                # unimodularity of the change-of-frame matrix; there is no
                # need to invoke general lattice isometry or glue induction.
                from sage.matrix.constructor import matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.recursive import _gram_in_label_order, _isometry_from_row_matrix

                left_labels = tuple(ambient.module_generating_set())
                right_labels = tuple(target_arrow.codomain().module_generating_set())
                left = matrix(SageZZ, [[int(source_arrow(point).to_vector()(label))
                    for label in left_labels] for point in source_arrow.domain().module_generators()])
                right = matrix(SageZZ, [[int(target_arrow(point).to_vector()(label))
                    for label in right_labels] for point in target_arrow.domain().module_generators()])
                if left.det() == 0 or right.det() == 0:
                    return None
                action = left.inverse() * right
                if any(not entry.is_integral() for entry in action.list()) or abs(action.det()) != 1:
                    return None
                action = action.change_ring(SageZZ)
                target_ambient = target_arrow.codomain()
                left_gram = _gram_in_label_order(ambient)
                right_gram = _gram_in_label_order(target_ambient)
                if action * right_gram * action.transpose() != left_gram:
                    return None
                return _isometry_from_row_matrix(ambient, target_ambient, action)
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
            nonhomothetic = source_saturation is None or target_saturation is None
            if nonhomothetic and (source_arrow.codomain().is_nondegenerate()
                                  or target_arrow.codomain().is_nondegenerate()):
                from sage.matrix.constructor import matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage.rings.rational_field import QQ as SageQQ

                saturations = (source_arrow.saturation(), target_arrow.saturation())
                inclusions = tuple(saturation.inclusion() for saturation in saturations)
                frames = tuple(tuple(inclusion(point) for point in inclusion.domain().module_generators())
                               for inclusion in inclusions)
                coordinates = tuple(matrix(SageZZ, [[int(inclusion.lift(arrow(point)).to_vector()(label))
                    for label in inclusion.domain().module_generating_set()]
                    for point in arrow.domain().module_generators()])
                    for arrow, inclusion in zip((source_arrow, target_arrow), inclusions, strict=True))
                if coordinates[0].nrows() != coordinates[0].ncols() or coordinates[0].det() == 0:
                    return None
                action = coordinates[0].change_ring(SageQQ).inverse() * coordinates[1]
                if any(entry.denominator() != 1 for entry in action.list()) or abs(action.det()) != 1:
                    return None
                target_images = tuple(sum((int(action[i, j]) * point
                    for j, point in enumerate(frames[1])), target_arrow.codomain().zero())
                    for i in range(len(frames[0])))
                witness = RECURSION.isotropic_subspace_isometry(frames[0], target_images, ())
                if witness is not None and all(witness(source_arrow(point)) == target_arrow(point)
                    for point in source_arrow.domain().module_generators()):
                    return witness
                return None
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
            if nonhomothetic:
                # A nonhomothetic finite-index marking has no common
                # scalar to cancel. Its exact based rows, rather than
                # individually primitive parts, determine the transporter.
                source_primitive, target_primitive = source_images, target_images
            if (not source_arrow.codomain().is_nondegenerate()
                    and not target_arrow.codomain().is_nondegenerate()
                    and len(source_primitive) >= 2):
                from sage.matrix.constructor import identity_matrix, matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage.rings.rational_field import QQ as SageQQ
                from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int
                from sage_indefinite_port.indefinite.recursive import _gram_in_label_order, _isometry_from_row_matrix

                source_ambient = source_arrow.codomain()
                target_ambient = target_arrow.codomain()
                ambient_witness = RECURSION.isometry(source_ambient, target_ambient)
                if ambient_witness is None:
                    return None
                reduction = indefinite_reduced_presentation(target_ambient)
                reduced = reduction.reduced
                gram = _gram_in_label_order(reduced)
                n, size = int(gram.rank()), int(gram.nrows())
                k = size - n
                labels = tuple(reduced.module_generating_set())
                frames = tuple(matrix(SageZZ, [[int((~reduction.isometry)(point).to_vector()(label))
                    for label in labels] for point in points]) for points in (
                        tuple(ambient_witness(point) for point in source_primitive), target_primitive))
                x0, x1 = (frame[:, :n] for frame in frames)
                if n and k and int(x0.rank()) == int(x1.rank()) > 0:
                    # Projection onto L/rad(L) is equivariant for every
                    # ambient isometry. Its integral elementary divisors
                    # are therefore a necessary invariant even for a
                    # *based* marking with equal generator divisibilities.
                    if tuple(abs(int(d)) for d in x0.elementary_divisors()) != tuple(
                        abs(int(d)) for d in x1.elementary_divisors()
                    ):
                        return None
                    smith, left, right = x0.smith_form()
                    s = int(x0.rank())
                    kernel = left[s:, :]
                    if kernel * x1 != 0:
                        return None
                    complement = reduced.lattice_category()([
                        [int(gram[i, j]) for j in range(n)] for i in range(n)])
                    c_labels = tuple(complement.module_generating_set())
                    projected = tuple(complement.subobject_on(tuple(
                        complement.linear_combination({
                            label: complement.base_ring()(int(row[j]))
                            for j, label in enumerate(c_labels) if row[j]
                        }) for row in rows.rows() if any(row))).saturation()
                        for rows in (x0, x1))
                    inclusions = tuple(term.inclusion() for term in projected)
                    saturated_bases = tuple(tuple(inclusion(point)
                        for point in inclusion.domain().module_generators())
                        for inclusion in inclusions)
                    coordinates = tuple(matrix(SageZZ, [[int(inclusion.lift(
                        complement.linear_combination({
                            label: complement.base_ring()(int(row[j]))
                            for j, label in enumerate(c_labels) if row[j]
                        })).to_vector()(label))
                        for label in inclusion.domain().module_generating_set()]
                        for row in rows.rows()]) for rows, inclusion in zip(
                            (x0, x1), inclusions, strict=True))
                    # The based map on the finite-index projections
                    # extends to their saturations precisely when its
                    # matrix in saturated bases is unimodular integral.
                    saturated_action = coordinates[0].change_ring(SageQQ).solve_right(
                        coordinates[1].change_ring(SageQQ))
                    if any(entry.denominator() != 1 for entry in saturated_action.list()) or abs(saturated_action.det()) != 1:
                        return None
                    target_projected = tuple(sum((int(saturated_action[i, j]) * point
                        for j, point in enumerate(saturated_bases[1])), complement.zero())
                        for i in range(s))
                    projected_witness = RECURSION.isotropic_subspace_isometry(
                        saturated_bases[0], target_projected, ())
                    if projected_witness is not None:
                        action = matrix(SageZZ, [[int(projected_witness(point).to_vector()(label))
                            for label in c_labels] for point in complement.module_generators()])
                        if x0 * action == x1:
                            r0, r1 = (kernel * frame[:, n:] for frame in frames)
                            try:
                                # The marked radical kernel need not be
                                # primitive. First solve its exact based
                                # right-GL(k,Z) equivalence by column HNF,
                                # then complete its primitive saturation.
                                h0, h0_change = r0.transpose().hermite_form(transformation=True)
                                h1, h1_change = r1.transpose().hermite_form(transformation=True)
                                if h0 != h1:
                                    return None
                                base_radical_action = (h0_change.transpose()
                                                       * h1_change.transpose().inverse())
                                if r0 * base_radical_action != r1:
                                    raise ArithmeticError("radical Hermite transporter does not map the based kernel")
                                primitive_kernel = r0.row_module().saturation().basis_matrix()
                                r0_frame = primitive_kernel.stack(
                                    _source_subspace_completion_int(primitive_kernel, k))
                                r1_frame = r0_frame * base_radical_action
                            except ValueError:
                                pass
                            else:
                                from math import lcm

                                t = int(r0.nrows())
                                modulus = lcm(*(abs(int(smith[i, i])) for i in range(s)))
                                standard_generators = []
                                for i in range(t, k):
                                    for j in range(k):
                                        if i == j:
                                            continue
                                        elementary = identity_matrix(SageZZ, k)
                                        elementary[i, j] = 1
                                        standard_generators.append(elementary)
                                for i in range(t, k - 1):
                                    swap = identity_matrix(SageZZ, k)
                                    swap.swap_rows(i, i + 1)
                                    standard_generators.append(swap)
                                if t < k:
                                    sign = identity_matrix(SageZZ, k)
                                    sign[t, t] = -1
                                    standard_generators.append(sign)
                                # The pointwise stabilizer of the radical
                                # kernel is Hom(Z^(k-t),Z^t) ⋊ GL_(k-t)(Z).
                                # Its finite reduction modulo the Smith
                                # exponent controls every shear congruence.
                                frontier = [identity_matrix(SageZZ, k)]
                                seen = set()
                                while frontier:
                                    word = frontier.pop()
                                    key = tuple(int(entry) % modulus for entry in word.list())
                                    if key in seen:
                                        continue
                                    seen.add(key)
                                    radical_action = (r0_frame.inverse().change_ring(SageZZ)
                                                      * word * r1_frame)
                                    discrepancy = frames[1][:, n:] - frames[0][:, n:] * radical_action
                                    corrected = left * discrepancy
                                    if (all(corrected[i, j] == 0 for i in range(s, len(source_primitive))
                                            for j in range(k))
                                        and all(corrected[i, j] % smith[i, i] == 0
                                                for i in range(s) for j in range(k))):
                                        solution = matrix(SageZZ, n, k, 0)
                                        for i in range(s):
                                            for j in range(k):
                                                solution[i, j] = corrected[i, j] // smith[i, i]
                                        shear = right * solution
                                        block = identity_matrix(SageZZ, size)
                                        block[:n, :n] = action
                                        block[:n, n:] = shear
                                        block[n:, n:] = radical_action
                                        witness = (reduction.isometry
                                                   * _isometry_from_row_matrix(reduced, reduced, block)
                                                   * (~reduction.isometry) * ambient_witness)
                                        if all(witness(a) == b for a, b in zip(
                                            source_primitive, target_primitive, strict=True)):
                                            return witness
                                    frontier.extend(generator * word for generator in standard_generators)
            return (None if nonhomothetic else RECURSION.isotropic_subspace_isometry(
                source_primitive, target_primitive, ()))
        case OrientedIsotropicFlags():
            # Any isometry g carrying the underlying flags, then a correction h in
            # Stab(F') with the graded signs of g: h g has every sign +1. The
            # orientation-preserving isometries are the coset Stab+(F') h g, so
            # no h means the sign vector of g lies outside the Levi sign image
            # of Stab(F'), and then no oriented isometry exists.
            source_flag, target_flag = source.functor().oriented_flag().flag(), target.functor().oriented_flag().flag()
            length = int(source_flag.flag_length())
            if int(target_flag.flag_length()) != length:
                return None
            flags = _marked_isotropic_category(source_flag.lattice(), length)
            witness = isometry(flags, _marked_isotropic_object(source_flag, flag=True), _marked_isotropic_object(target_flag, flag=True))
            if witness is None:
                return None
            _generators, _signs, transversal = _orientation_group_elements(target_flag.lattice().O(), target_flag)
            correction = transversal.get(_graded_orientation_signs(witness, source_flag, target_flag))
            return None if correction is None else correction * witness
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
                from sage.matrix.constructor import matrix
                from sage.matrix.special import identity_matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

                source_labels = tuple(ambient.module_generating_set())
                target_labels = tuple(target_arrow.codomain().module_generating_set())
                left = tuple(source_arrow(point) for point in source_arrow.domain().module_generators())
                right = tuple(target_arrow(point) for point in target_arrow.domain().module_generators())
                left_rows = matrix(SageZZ, [[int(point.to_vector()(label)) for label in source_labels] for point in left])
                right_rows = matrix(SageZZ, [[int(point.to_vector()(label)) for label in target_labels] for point in right])
                if left_rows.rank() != right_rows.rank():
                    return None
                left_smith, _left_change, left_basis = left_rows.smith_form(transformation=True)
                right_smith, _right_change, right_basis = right_rows.smith_form(transformation=True)
                if left_smith != right_smith:
                    return None
                dimension = int(left_rows.ncols())
                # Arrow-category isometries may act by an arbitrary
                # isomorphism on the source subobject; unlike coslice
                # isometries, they need not preserve its ordered frame.
                # Smith transforms identify the two sublattices with the
                # same diagonal module, so the identity there suffices.
                rows = left_basis * right_basis.inverse()
                if abs(rows.det()) != 1:
                    return None
                witness = _isometry_from_row_matrix(ambient, target_arrow.codomain(), rows)
                target_subobject = target_arrow.codomain().subobject_on(right)
                if not all(target_subobject.inclusion().is_in_image(witness(point)) for point in left):
                    raise ArithmeticError("Smith arrow transporter failed its image certificate")
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
                    # CombinedAlgorithms.h:1401-1404: unequal isotropic k-plane invariants
                    # (INDEF_FORM_Invariant_IsotropicKplane, :2424-2428) decide inequivalence.
                    # The pinned invariant uses a nondegenerate radical
                    # pairing L -> I*. It is undefined on a degenerate
                    # ambient lattice; its exact transporter branch below
                    # handles such lattices without this prefilter.
                    if source.arrow().codomain().is_nondegenerate() and target.arrow().codomain().is_nondegenerate():
                        if INDUCTION.flag_invariant(source.arrow().codomain(), INDUCTION.plane_vectors(source_sublattice), (source_rank,)) != INDUCTION.flag_invariant(
                            target.arrow().codomain(), INDUCTION.plane_vectors(target_sublattice), (target_rank,)
                        ):
                            return None
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
            # CombinedAlgorithms.h:1401-1404 uses an invariant of the
            # *nondegenerate* ambient form. It constructs an invertible
            # rational pairing with the null subspace; that pairing has
            # non-full rank on a radical ambient. Degenerate flags instead
            # use their characteristic-radical and projected Smith data
            # below. Invoking the nondegenerate prefilter here wrongly
            # rejects even two equivalent flags in zero-form lattices.
            if ambient.is_nondegenerate() and target_ambient.is_nondegenerate():
                if INDUCTION.flag_invariant(ambient, source_basis, source_dimensions) != INDUCTION.flag_invariant(
                    target_ambient, target_basis, target_dimensions
                ):
                    return None
            if ambient is target_ambient and len(source_basis) == len(target_basis):
                # Extend the prescribed matching of coordinate generators
                # to a permutation of the full coordinate basis. This
                # includes overlapping frames and cyclic permutations, not
                # only disjoint transpositions. Verify the full Gram form
                # before accepting the resulting integral isometry.
                from sage.matrix.constructor import matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.recursive import _gram_in_label_order, _isometry_from_row_matrix

                # The Gram matrix and the row-matrix isometry constructor
                # both use coordinate-label order, which need not coincide
                # with an arbitrary module-generator iteration order.
                generators = tuple(ambient.module_generator(label)
                                   for label in ambient.module_generating_set())
                def signed_coordinate(point):
                    for index, generator in enumerate(generators):
                        if point == generator:
                            return index, 1
                        if point == -generator:
                            return index, -1
                    return None

                coordinates = tuple(signed_coordinate(point) for point in (*source_basis, *target_basis))
                if all(coordinate is not None for coordinate in coordinates):
                    s = len(source_basis)
                    fixed = {coordinates[i][0]: (coordinates[s + i][0],
                              coordinates[i][1] * coordinates[s + i][1]) for i in range(s)}
                    gram = _gram_in_label_order(ambient)
                    # In a hyperbolic direct sum the Gram form determines
                    # each coordinate's unique dual. Extend the prescribed
                    # isotropic-coordinate action to these duals as well.
                    for i, (j, sign) in tuple(fixed.items()):
                        source_duals = [p for p in range(len(generators)) if gram[i, p]]
                        target_duals = [p for p in range(len(generators)) if gram[j, p]]
                        if (len(source_duals) == len(target_duals) == 1
                                and abs(gram[i, source_duals[0]]) == abs(gram[j, target_duals[0]])):
                            dual_sign = sign * int(gram[i, source_duals[0]] / gram[j, target_duals[0]])
                            fixed.setdefault(source_duals[0], (target_duals[0], dual_sign))
                    if len({entry[0] for entry in fixed.values()}) != len(fixed):
                        fixed = {}
                    unused_sources = [i for i in range(len(generators)) if i not in fixed]
                    unused_targets = [i for i in range(len(generators)) if i not in
                                      {entry[0] for entry in fixed.values()}]
                    fixed.update((i, (j, 1)) for i, j in zip(unused_sources, unused_targets, strict=True))
                    permutation = matrix(SageZZ, len(generators), len(generators),
                                         lambda i, j: fixed[i][1] if fixed[i][0] == j else 0)
                    if permutation * gram * permutation.transpose() == gram:
                        witness = _isometry_from_row_matrix(ambient, ambient, permutation)
                        if all(witness(point) == target for point, target in zip(
                                source_basis, target_basis, strict=True)):
                            return witness
            if not ambient.is_nondegenerate() and not target_ambient.is_nondegenerate():
                # Each stage's intersection with the characteristic
                # radical is invariant under ambient isometries. Test
                # every stage, not just the terminal isotropic subspace.
                from sage.matrix.constructor import matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.recursive import _gram_in_label_order
                intersection_profiles = []
                projection_smith_profiles = []
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis)):
                    radical = owner.radical()
                    radical_points = tuple(radical.inclusion()(point) for point in radical.module_generators())
                    intersection_profiles.append(tuple(
                        rank + int(radical.rank()) - int(owner.subobject_on(
                            (*basis[:rank], *radical_points)
                        ).rank()) for rank in source_dimensions
                    ))
                    reduction = indefinite_reduced_presentation(owner)
                    reduced = reduction.reduced
                    labels = tuple(reduced.module_generating_set())[:int(_gram_in_label_order(reduced).rank())]
                    rows = tuple(tuple(int((~reduction.isometry)(point).to_vector()(label))
                                       for label in labels) for point in basis)
                    projection_smith_profiles.append(tuple(tuple(
                        abs(int(entry)) for entry in matrix(SageZZ, rows[:rank]).elementary_divisors()
                        if entry
                    ) for rank in source_dimensions))
                if intersection_profiles[0] != intersection_profiles[1]:
                    return None
                if projection_smith_profiles[0] != projection_smith_profiles[1]:
                    return None
            source_radical_flag = not ambient.is_nondegenerate() and all(
                RECURSION.vector_lies_in_radical(point) for point in source_basis
            )
            target_radical_flag = not target_ambient.is_nondegenerate() and all(
                RECURSION.vector_lies_in_radical(point) for point in target_basis
            )
            if source_radical_flag != target_radical_flag:
                return None
            if source_radical_flag:
                # GL(rad L) is transitive on primitive flags with the
                # same rank sequence. Complete both adapted flag bases
                # to unimodular bases of the radical and lift the exact
                # transporter with identity on the nondegenerate quotient.
                from sage.matrix.constructor import identity_matrix, matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int
                from sage_indefinite_port.indefinite.recursive import _gram_in_label_order, _isometry_from_row_matrix

                phi = RECURSION.isometry(ambient, target_ambient)
                if phi is None:
                    return None
                reduction = indefinite_reduced_presentation(target_ambient)
                reduced = reduction.reduced
                gram = _gram_in_label_order(reduced)
                n, size = int(gram.rank()), int(gram.nrows())
                labels = tuple(reduced.module_generating_set())[n:]
                completed = []
                for points in (tuple(phi(point) for point in source_basis), target_basis):
                    rows = matrix(SageZZ, [[int((~reduction.isometry)(point).to_vector()(label))
                                            for label in labels] for point in points])
                    completed.append(rows.stack(_source_subspace_completion_int(rows, size - n)))
                action = completed[0].inverse().change_ring(SageZZ) * completed[1]
                block = identity_matrix(SageZZ, size)
                for i in range(size - n):
                    for j in range(size - n):
                        block[n + i, n + j] = action[i, j]
                witness = reduction.isometry * _isometry_from_row_matrix(reduced, reduced, block) * (~reduction.isometry) * phi
                target_terms = tuple(target_ambient.primitive_isotropic_subobject(*target_basis[:rank])
                                     for rank in target_dimensions)
                if not all(term.inclusion().is_in_image(witness(point))
                           for rank, term in zip(source_dimensions, target_terms, strict=True)
                           for point in source_basis[:rank]):
                    raise ArithmeticError("the lifted radical-flag transporter does not preserve every term")
                return witness
            if len(source_dimensions) >= 2 and all(
                not owner.is_nondegenerate()
                and any(all(RECURSION.vector_lies_in_radical(point) for point in basis[:rank])
                        and rank == source_dimensions[-1] + int(owner.radical().rank())
                        - int(owner.subobject_on((*basis, *(
                    owner.radical().inclusion()(point)
                    for point in owner.radical().module_generators()
                ))).rank()) for rank in source_dimensions[:-1])
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis))
            ):
                # A stage is exactly J∩rad(L) on both sides.
                # Any transporter of the top isotropic subspaces must
                # carry this characteristic intersection to its target.
                return INDUCTION.plane_isometry(
                    ambient.primitive_isotropic_subobject(*source_basis),
                    target_ambient.primitive_isotropic_subobject(*target_basis),
                    flag_dimensions=source_dimensions,
                    source_flag_basis=source_basis,
                    target_flag_basis=target_basis,
                )
            if len(source_dimensions) >= 2 and all(
                not owner.is_nondegenerate()
                and (t := source_dimensions[-1] + int(owner.radical().rank())
                     - int(owner.subobject_on((*basis, *(
                         owner.radical().inclusion()(point)
                         for point in owner.radical().module_generators()
                     ))).rank())) > 0
                and t < source_dimensions[-1] and t not in source_dimensions
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis))
            ):
                # The characteristic radical intersection may be implicit
                # in the supplied flag. Insert it only when it is nested
                # with the given stages, adapting within each interval.
                adapted_pairs = []
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis)):
                    t = (source_dimensions[-1] + int(owner.radical().rank())
                         - int(owner.subobject_on((*basis, *(
                             owner.radical().inclusion()(point)
                             for point in owner.radical().module_generators()
                         ))).rank()))
                    previous = 0
                    ordered = []
                    for rank in source_dimensions:
                        ordered.extend(sorted(basis[previous:rank], key=lambda point:
                                              not RECURSION.vector_lies_in_radical(point)))
                        previous = rank
                    if not all(RECURSION.vector_lies_in_radical(point) for point in ordered[:t]):
                        break
                    adapted_pairs.append((t, tuple(ordered)))
                if len(adapted_pairs) == 2 and adapted_pairs[0][0] == adapted_pairs[1][0]:
                    t = adapted_pairs[0][0]
                    augmented = tuple(sorted((*source_dimensions, t)))
                    source_adapted, target_adapted = (pair[1] for pair in adapted_pairs)
                    witness = INDUCTION.plane_isometry(
                        ambient.primitive_isotropic_subobject(*source_adapted),
                        target_ambient.primitive_isotropic_subobject(*target_adapted),
                        flag_dimensions=augmented,
                        source_flag_basis=source_adapted,
                        target_flag_basis=target_adapted,
                    )
                    if witness is not None:
                        return witness
            if len(source_dimensions) >= 2 and all(
                not owner.is_nondegenerate()
                and (t := source_dimensions[-1] + int(owner.radical().rank())
                     - int(owner.subobject_on((*basis, *(
                         owner.radical().inclusion()(point)
                         for point in owner.radical().module_generators()
                     ))).rank())) > 0
                and t < source_dimensions[-1]
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis))
            ):
                # Transport the characteristic intersection before testing
                # stages transverse to it, including projected first lines.
                radical_pairs = []
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis)):
                    t = (source_dimensions[-1] + int(owner.radical().rank())
                         - int(owner.subobject_on((*basis, *(
                             owner.radical().inclusion()(point)
                             for point in owner.radical().module_generators()
                         ))).rank()))
                    ordered = tuple(sorted(basis, key=lambda point:
                                           not RECURSION.vector_lies_in_radical(point)))
                    if not all(RECURSION.vector_lies_in_radical(point) for point in ordered[:t]):
                        break
                    radical_pairs.append((t, ordered))
                if len(radical_pairs) == 2 and radical_pairs[0][0] == radical_pairs[1][0]:
                    t = radical_pairs[0][0]
                    target_terms = tuple(target_ambient.primitive_isotropic_subobject(*target_basis[:rank])
                                         for rank in source_dimensions)
                    for dimensions in (tuple(sorted(set((*source_dimensions, t)))),
                                       (t, source_dimensions[-1])):
                        witness = INDUCTION.plane_isometry(
                            ambient.primitive_isotropic_subobject(*radical_pairs[0][1]),
                            target_ambient.primitive_isotropic_subobject(*radical_pairs[1][1]),
                            flag_dimensions=dimensions,
                            source_flag_basis=radical_pairs[0][1],
                            target_flag_basis=radical_pairs[1][1],
                        )
                        if witness is not None and all(
                            term.inclusion().is_in_image(witness(point))
                            for rank, term in zip(source_dimensions, target_terms, strict=True)
                            for point in source_basis[:rank]
                        ):
                            return witness
            if all(
                not owner.is_nondegenerate()
                and source_dimensions[0] >= int(owner.radical().rank())
                and source_dimensions[-1] > int(owner.radical().rank())
                and all(owner.primitive_isotropic_subobject(*basis[:source_dimensions[0]]).inclusion().is_in_image(
                    owner.radical().inclusion()(point)) for point in owner.radical().module_generators())
                for owner, basis in ((ambient, source_basis), (target_ambient, target_basis))
            ):
                return INDUCTION.plane_isometry(
                    ambient.primitive_isotropic_subobject(*source_basis),
                    target_ambient.primitive_isotropic_subobject(*target_basis),
                    flag_dimensions=source_dimensions,
                    source_flag_basis=source_basis,
                    target_flag_basis=target_basis,
                )
            if (len(source_dimensions) >= 2 and source_dimensions[0] == 1
                    and not ambient.is_nondegenerate() and not target_ambient.is_nondegenerate()
                    and not RECURSION.vector_lies_in_radical(source_basis[0])
                    and not RECURSION.vector_lies_in_radical(target_basis[0])):
                # First align the terminal plane. If the projected first
                # lines agree, a primitive radical-complement shear can
                # then align their actual graph lifts while fixing J∩R.
                from sage.matrix.constructor import identity_matrix, matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int
                from sage_indefinite_port.indefinite.recursive import _gram_in_label_order, _isometry_from_row_matrix

                radical_tail = all(RECURSION.vector_lies_in_radical(point)
                                   for point in (*source_basis[1:], *target_basis[1:]))
                source_top = (*source_basis[1:], source_basis[0]) if radical_tail else source_basis
                target_top = (*target_basis[1:], target_basis[0]) if radical_tail else target_basis
                if not radical_tail:
                    # The adapted basis of the requested flag need not
                    # display its characteristic radical intersection:
                    # e.g. (2u+t+z,w,2u+t) contains z as a difference.
                    # Extract the integral kernel of projection on each
                    # stage and adapt its nested radical flag exactly.
                    radical_tops = []
                    for owner, basis in ((ambient, source_basis), (target_ambient, target_basis)):
                        reduction = indefinite_reduced_presentation(owner)
                        reduced = reduction.reduced
                        labels = tuple(reduced.module_generating_set())
                        n = int(_gram_in_label_order(reduced).rank())
                        rows = matrix(SageZZ, [[int((~reduction.isometry)(point).to_vector()(label))
                                                for label in labels] for point in basis])
                        projection_rank = int(rows[:, :n].rank())
                        if not projection_rank:
                            break
                        radical_terms = []
                        for rank in source_dimensions:
                            prefix = rows[:rank, :]
                            _diagonal, left, _right = prefix[:, :n].smith_form()
                            kernel = (left * prefix)[int(prefix[:, :n].rank()):, :]
                            points = tuple(reduction.isometry(reduced.linear_combination({
                                label: reduced.base_ring()(int(row[j]))
                                for j, label in enumerate(labels) if row[j]
                            })) for row in kernel.rows())
                            if points and (not radical_terms or len(points) > int(radical_terms[-1].rank())):
                                radical_terms.append(owner.primitive_isotropic_subobject(*points))
                        if not radical_terms or int(radical_terms[-1].rank()) != len(basis) - projection_rank:
                            break
                        radical_top = radical_terms[-1]
                        domain = radical_top.inclusion().domain()
                        nested = tuple(domain.primitive_isotropic_subobject(*(
                            radical_top.inclusion().lift(term.inclusion()(point))
                            for point in term.module_generators()
                        )) for term in radical_terms)
                        radical_basis, _ = INDUCTION.flag_action_data(radical_top, nested)
                        projected_basis = []
                        projected_rows = matrix(SageZZ, 0, n)
                        for point, row in zip(basis, rows.rows(), strict=True):
                            extension = projected_rows.stack(matrix(SageZZ, [list(row[:n])]))
                            if int(extension.rank()) > int(projected_rows.rank()):
                                projected_basis.append(point)
                                projected_rows = extension
                        reordered = (*radical_basis, *projected_basis)
                        original_rows = matrix(SageZZ, [[int(point.to_vector()(label))
                            for label in owner.module_generating_set()] for point in basis])
                        reordered_rows = matrix(SageZZ, [[int(point.to_vector()(label))
                            for label in owner.module_generating_set()] for point in reordered])
                        if original_rows.hermite_form() != reordered_rows.hermite_form():
                            break
                        reordered_top = owner.primitive_isotropic_subobject(*reordered)
                        if not all(reordered_top.inclusion().is_in_image(point) for point in basis):
                            break
                        radical_tops.append(reordered)
                    if len(radical_tops) == 2:
                        source_top, target_top = radical_tops
                        radical_tail = True
                top_witness = INDUCTION.plane_isometry(
                    ambient.primitive_isotropic_subobject(*source_top),
                    target_ambient.primitive_isotropic_subobject(*target_top),
                    flag_dimensions=tuple(range(1, len(source_basis) + 1)) if radical_tail else (),
                    source_flag_basis=source_top if radical_tail else (),
                    target_flag_basis=target_top if radical_tail else (),
                )
                if top_witness is not None:
                    target_line = target_ambient.primitive_isotropic_subobject(target_basis[0])
                    target_stages = tuple(target_ambient.primitive_isotropic_subobject(*target_basis[:rank])
                                          for rank in source_dimensions)
                    if all(stage.inclusion().is_in_image(top_witness(point))
                           for rank, stage in zip(source_dimensions, target_stages, strict=True)
                           for point in source_basis[:rank]):
                        return top_witness
                    reduction = indefinite_reduced_presentation(target_ambient)
                    reduced = reduction.reduced
                    gram = _gram_in_label_order(reduced)
                    n, size = int(gram.rank()), int(gram.nrows())
                    labels = tuple(reduced.module_generating_set())
                    rows = matrix(SageZZ, [[int((~reduction.isometry)(point).to_vector()(label))
                                            for label in labels] for point in target_basis])
                    if int(rows[:, :n].rank()) == 1:
                        _smith, left, _right = rows[:, :n].smith_form()
                        radical_rows = (left * rows)[1:, n:]
                        if int(radical_rows.rank()) == int(radical_rows.nrows()):
                            first = matrix(SageZZ, [[int((~reduction.isometry)(point).to_vector()(label))
                                                     for label in labels]
                                                    for point in (top_witness(source_basis[0]), target_basis[0])])
                            for sign in (1, -1):
                                if first[0, :n] != sign * first[1, :n]:
                                    continue
                                source_r = first[0, n:]
                                target_r = sign * first[1, n:]
                                try:
                                    source_frame = radical_rows.stack(source_r)
                                    target_frame = radical_rows.stack(target_r)
                                    source_frame = source_frame.stack(_source_subspace_completion_int(source_frame, size - n))
                                    target_frame = target_frame.stack(_source_subspace_completion_int(target_frame, size - n))
                                except ValueError:
                                    continue
                                radical_action = source_frame.inverse().change_ring(SageZZ) * target_frame
                                block = identity_matrix(SageZZ, size)
                                block[n:, n:] = radical_action
                                correction = (reduction.isometry
                                              * _isometry_from_row_matrix(reduced, reduced, block)
                                              * (~reduction.isometry))
                                candidate = correction * top_witness
                                if (target_line.inclusion().is_in_image(candidate(source_basis[0]))
                                        and all(target_ambient.primitive_isotropic_subobject(*target_basis[:rank])
                                                .inclusion().is_in_image(candidate(point))
                                                for rank in source_dimensions for point in source_basis[:rank])):
                                    return candidate
            if (not ambient.is_nondegenerate() and not target_ambient.is_nondegenerate()
                    and len(source_dimensions) >= 2
                    and all(all(RECURSION.vector_lies_in_radical(point)
                                or all(int((~reduction.isometry)(point).to_vector()(label)) == 0
                                       for label in tuple(reduction.reduced.module_generating_set())[
                                           int(reduction.reduced.rank()) - int(owner.radical().rank()):])
                                for point in basis)
                            for owner, basis, reduction in (
                                (ambient, source_basis, indefinite_reduced_presentation(ambient)),
                                (target_ambient, target_basis, indefinite_reduced_presentation(target_ambient))))):
                # The projected and radical constituent flags are split;
                # the integral shear may therefore align their graph-free
                # realizations after the quotient and radical transports.
                return INDUCTION.plane_isometry(
                    ambient.primitive_isotropic_subobject(*source_basis),
                    target_ambient.primitive_isotropic_subobject(*target_basis),
                    flag_dimensions=source_dimensions,
                    source_flag_basis=source_basis,
                    target_flag_basis=target_basis,
                )
            if (int(ambient.rank()) == int(target_ambient.rank())
                    and all(ambient.b(x, y) == 0 for x in ambient.module_generators() for y in ambient.module_generators())
                    and all(target_ambient.b(x, y) == 0 for x in target_ambient.module_generators()
                            for y in target_ambient.module_generators())):
                # An adapted integral basis of each primitive top flag
                # gives an explicit GL_n(ZZ) transporter. No unbounded
                # unimodular enumeration is necessary.
                from sage.matrix.constructor import matrix
                from sage.rings.integer_ring import ZZ as SageZZ
                from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int
                from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

                source_labels = tuple(ambient.module_generating_set())
                target_labels = tuple(target_ambient.module_generating_set())
                source_rows = matrix(SageZZ, [[int(point.to_vector()(label))
                    for label in source_labels] for point in source_basis])
                target_rows = matrix(SageZZ, [[int(point.to_vector()(label))
                    for label in target_labels] for point in target_basis])
                source_frame = source_rows.stack(_source_subspace_completion_int(
                    source_rows, int(ambient.rank())))
                target_frame = target_rows.stack(_source_subspace_completion_int(
                    target_rows, int(target_ambient.rank())))
                rows = source_frame.inverse() * target_frame
                witness = _isometry_from_row_matrix(ambient, target_ambient, rows)
                target_terms = tuple(target_ambient.primitive_isotropic_subobject(*target_basis[:rank])
                                     for rank in source_dimensions)
                if all(term.inclusion().is_in_image(witness(point))
                       for rank, term in zip(source_dimensions, target_terms, strict=True)
                       for point in source_basis[:rank]):
                    return witness
                raise ArithmeticError("integral adapted flag completion failed its transport certificate")
            try:
                return RECURSION.isotropic_subspace_isometry(source_basis, target_basis, target_images)
            except ValueError as exc:
                if "partial isometry does not identify the ambient radical intersections" in str(exc):
                    # An isometry must preserve the marking's intersection
                    # with the characteristic radical. This obstruction is
                    # checked by the exact partial-extension owner.
                    return None
                raise
    raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices), Fun([t], Lattices) for finite [t] or oriented isotropic flags")


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
        if int(source.rank()) == 1 and not ambient.is_nondegenerate():
            # The rank-one primitive mark can project to an *imprimitive*
            # vector in L/rad(L). For example, 5e+r in <2> + <0>
            # is primitive in L but has quotient projection 5e.
            # degenerate_marking_generators assumes that the projected
            # span is primitive and consequently rejects this valid mark.
            # The exact pointed algorithm already solves the integral
            # radical congruence for arbitrary projection content.
            point = arrow(next(iter(source.module_generators())))
            if point == ambient.zero():
                return finite_ordered_set(tuple(RECURSION.orthogonal_group_generators(ambient)))
            primitive = point.primitive_part()
            fixed = tuple(RECURSION.pointed_automorphism_generators(primitive))
            if isinstance(category, CosliceCategory):
                return finite_ordered_set(fixed)
            reversal = RECURSION.pointed_isometry(primitive, -primitive)
            return finite_ordered_set((*fixed, *((reversal,) if reversal is not None else ())))
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
                ambient = marked.arrow().codomain()
                generators = (_zero_form_orthogonal_generators(ambient, ambient.general_linear_group().group_generators())
                              if all(ambient.b(x, y) == 0 for x in ambient.module_generators()
                                     for y in ambient.module_generators())
                              else RECURSION.orthogonal_group_generators(ambient))
            elif rank == 1:
                generators = RECURSION.pointed_automorphism_generators(_marked_vector(marked).primitive_part())
            elif not _is_totally_isotropic(marked.arrow()) and _is_nondegenerate_primitive(marked.arrow()):
                from sage_extensions.categories.lattice_embeddings import _nondegenerate_marking_stabilizer_generators

                generators = _nondegenerate_marking_stabilizer_generators(marked.arrow(), pointwise=True)
            else:
                arrow = marked.arrow()
                if all(
                    arrow.codomain().b(x, y) == 0
                    for x in arrow.codomain().module_generators()
                    for y in arrow.codomain().module_generators()
                ):
                    from sage.matrix.special import identity_matrix
                    from sage.rings.integer_ring import ZZ as SageZZ
                    from sage.matrix.constructor import matrix
                    from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int
                    from sage_indefinite_port.indefinite.recursive import _isometry_from_row_matrix

                    ambient = arrow.codomain()
                    n, k = int(ambient.rank()), int(arrow.domain().rank())
                    labels = tuple(ambient.module_generating_set())
                    # Pointwise fixation of an integral frame is equivalent
                    # to pointwise fixation of its primitive saturation:
                    # if m*x is fixed, torsion-freeness forces x fixed.
                    frame_arrow = arrow if arrow.is_primitive() else arrow.saturation().inclusion()
                    fixed = tuple(frame_arrow(point) for point in frame_arrow.domain().module_generators())
                    frame = tuple(tuple(int(point.to_vector()(label)) for label in labels) for point in fixed)
                    fixed_rows = matrix(SageZZ, frame)
                    completion = fixed_rows.stack(_source_subspace_completion_int(fixed_rows, n))
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
                elif int(arrow.domain().rank()) == int(arrow.codomain().rank()):
                    # An injective full-rank marking spans the ambient
                    # rational vector space. Any isometry fixing its frame
                    # pointwise is therefore the identity, regardless of
                    # degeneracy or finite index of the integral embedding.
                    return finite_ordered_set((arrow.codomain().O().identity(),))
                elif arrow.domain().is_nondegenerate():
                    generators = (_nondegenerate_stabilizer_generators(arrow, based=True)
                                  if arrow.is_primitive() else _nondegenerate_stabilizer_generators(arrow.saturation().inclusion(), based=True))
                    return finite_ordered_set(tuple(generators))
                saturation = _isotropic_sublattice(arrow) if arrow.is_primitive() else _homothetic_isotropic_saturation(arrow)
                mixed_projection = False
                if not arrow.codomain().is_nondegenerate() and int(arrow.domain().rank()) >= 2:
                    from sage.matrix.constructor import matrix
                    from sage.rings.integer_ring import ZZ as SageZZ

                    reduction, n, _k, _quotient = _degenerate_quotient_and_radical(arrow.codomain())
                    labels = tuple(reduction.reduced.module_generating_set())[:n]
                    projected_rows = matrix(SageZZ, [[int((~reduction.isometry)(arrow(point))
                        .to_vector()(label)) for label in labels]
                        for point in arrow.domain().module_generators()])
                    mixed_projection = int(projected_rows.rank()) > 0
                if saturation is None or (not arrow.codomain().is_nondegenerate()
                                          and mixed_projection):
                    if arrow.codomain().is_nondegenerate():
                        saturation = arrow.saturation()
                    else:
                        from sage.matrix.constructor import identity_matrix, matrix
                        from sage.rings.integer_ring import ZZ as SageZZ
                        from sage.rings.rational_field import QQ as SageQQ
                        from sage_extensions.all import Modules
                        from sage_extensions.categories.modules.orbit_spans import OrbitSpanConstructions
                        from sage_extensions.categories.modules.pure.modules import ModulesOverIntegralDomains
                        from sage_extensions.categories.rational_integral_stabilizers import IntegralStructureAction
                        from sage_indefinite_port.indefinite.indefinite_reduction import _source_subspace_completion_int
                        from sage_indefinite_port.indefinite.recursive import _conjugated, _gram_in_label_order, _isometry_from_row_matrix, _integral_generator_set

                        ambient = arrow.codomain()
                        saturation = arrow.saturation()
                        reduction = indefinite_reduced_presentation(ambient)
                        reduced = reduction.reduced
                        gram = _gram_in_label_order(reduced)
                        n, size = int(gram.rank()), int(gram.nrows())
                        k = size - n
                        labels = tuple(reduced.module_generating_set())
                        rows = matrix(SageZZ, [[int((~reduction.isometry)(saturation.inclusion()(point))
                            .to_vector()(label)) for label in labels]
                            for point in saturation.module_generators()])
                        projected = rows[:, :n]
                        smith, left, _right = projected.smith_form()
                        if n and k and int(projected.rank()) == 1:
                            adapted = left * rows
                            xrow, yrow = adapted.row(0)[:n], adapted.row(0)[n:]
                            radical_rows = adapted[1:, n:]
                            t = int(radical_rows.nrows())
                            content = abs(int(smith[0, 0]))
                            primitive = matrix(SageZZ, [[int(entry) // content for entry in xrow]])
                            cframe = primitive.stack(_source_subspace_completion_int(primitive, n))
                            rframe = radical_rows.stack(_source_subspace_completion_int(radical_rows, k))
                            complement = reduced.lattice_category()([
                                [int(gram[i, j]) for j in range(n)] for i in range(n)])
                            c_labels = tuple(complement.module_generating_set())
                            line = complement.linear_combination({label: complement.base_ring()(int(primitive[0, j]))
                                for j, label in enumerate(c_labels) if primitive[0, j]})
                            ring = SageZZ if content == 1 else SageQQ
                            normalizer = identity_matrix(ring, size)
                            shear = matrix(ring, n, k, 0)
                            shear[0, :] = -yrow / content
                            normalizer[:n, n:] = cframe.inverse() * shear
                            blocks = []
                            for action in RECURSION.pointwise_stabilizer_generators(line):
                                block = identity_matrix(SageZZ, size)
                                for i, point in enumerate(complement.module_generators()):
                                    for j, label in enumerate(c_labels):
                                        block[i, j] = int(action(point).to_vector()(label))
                                blocks.append(block)
                            for i in range(t, k):
                                for j in range(t):
                                    block = identity_matrix(SageZZ, size)
                                    standard = identity_matrix(SageZZ, k)
                                    standard[i, j] = 1
                                    block[n:, n:] = rframe.inverse() * standard * rframe
                                    blocks.append(block)
                            if t < k:
                                radical_complement = reduced.lattice_category()([[0] * (k - t) for _ in range(k - t)])
                                r_labels = tuple(radical_complement.module_generating_set())
                                for action in radical_complement.general_linear_group().group_generators():
                                    standard = identity_matrix(SageZZ, k)
                                    for i, point in enumerate(radical_complement.module_generators()):
                                        for j, label in enumerate(r_labels):
                                            standard[t + i, t + j] = int(action(point).to_vector()(label))
                                    block = identity_matrix(SageZZ, size)
                                    block[n:, n:] = rframe.inverse() * standard * rframe
                                    blocks.append(block)
                            for i in range(1, n):
                                for j in range(k):
                                    standard = matrix(SageZZ, n, k, 0)
                                    standard[i, j] = 1
                                    block = identity_matrix(SageZZ, size)
                                    block[:n, n:] = cframe.inverse() * standard
                                    blocks.append(block)
                            if content > 1:
                                rational = reduced.base_change(reduced.base_ring().fraction_field_map())
                                inclusion_r = ModulesOverIntegralDomains(reduced.base_ring())(reduced).generic_fibre_map()
                                conjugates = []
                                for block in blocks:
                                    conjugated = normalizer * block * normalizer.inverse()
                                    conjugates.append(rational.Isom(rational)(tuple(
                                        rational.linear_combination({label: rational.base_ring()(conjugated[i, j])
                                            for j, label in enumerate(labels) if conjugated[i, j]})
                                        for i in range(size))))
                                action = IntegralStructureAction(rational.O().subgroup(tuple(conjugates)), inclusion_r)
                                bound = content ** 2
                                containing = Modules(reduced.base_ring()).Mono(reduced, inclusion_r.codomain())(
                                    lambda label: inclusion_r.codomain().wrap(
                                        rational.scalar_multiple(rational.base_ring().one() / rational.base_ring()(bound),
                                            inclusion_r(reduced.module_generator(label)).underlying_element())))
                                span = OrbitSpanConstructions(reduced.base_ring())(
                                    action.ambient_action_functor(), inclusion_r, containing=containing)
                                generators = _conjugated(reduction.isometry, tuple(_integral_generator_set(
                                    reduced, action.with_stable_containing_bound(span.inclusion()))))
                            else:
                                generators = _conjugated(reduction.isometry, tuple(
                                    _isometry_from_row_matrix(reduced, reduced,
                                        normalizer * block * normalizer.inverse())
                                    for block in blocks))
                            generators = generators or (ambient.O().identity(),)
                            if not all(generator(arrow(point)) == arrow(point)
                                       for generator in generators
                                       for point in arrow.domain().module_generators()):
                                raise ArithmeticError("the mixed based pointwise stabilizer failed its frame certificate")
                            return finite_ordered_set(tuple(generators))
                        s = int(projected.rank())
                        if n and k and s > 1:
                            # The rational projected plane splits off its graph.
                            # Intersect the resulting rational parabolic with
                            # the integral lattice when projection has index.
                            # Its pointwise Levi factor is the
                            # actual isotropic stabilizer of the quotient,
                            # not an assumed full orthogonal group.
                            adapted = left * rows
                            xrows, yrows = adapted[:s, :n], adapted[:s, n:]
                            radical_rows = adapted[s:, n:]
                            t = int(radical_rows.nrows())
                            divisors = tuple(abs(int(smith[i, i])) for i in range(s))
                            primitive_rows = matrix(SageZZ, [
                                [int(xrows[i, j]) // divisors[i] for j in range(n)]
                                for i in range(s)])
                            cframe = primitive_rows.stack(_source_subspace_completion_int(primitive_rows, n))
                            rframe = radical_rows.stack(_source_subspace_completion_int(radical_rows, k))
                            complement = reduced.lattice_category()([
                                [int(gram[i, j]) for j in range(n)] for i in range(n)])
                            c_labels = tuple(complement.module_generating_set())
                            quotient_points = tuple(complement.linear_combination({
                                label: complement.base_ring()(int(primitive_rows[i, j]))
                                for j, label in enumerate(c_labels) if primitive_rows[i, j]
                            }) for i in range(s))
                            ring = SageZZ if all(d == 1 for d in divisors) else SageQQ
                            normalizer = identity_matrix(ring, size)
                            graph = matrix(ring, n, k, 0)
                            for i in range(s):
                                graph[i, :] = -yrows.row(i) / divisors[i]
                            normalizer[:n, n:] = cframe.inverse() * graph
                            blocks = []
                            for action in INDUCTION._stabilizer_generators(quotient_points, ()):
                                block = identity_matrix(SageZZ, size)
                                for i, point in enumerate(complement.module_generators()):
                                    for j, label in enumerate(c_labels):
                                        block[i, j] = int(action(point).to_vector()(label))
                                blocks.append(block)
                            for i in range(t, k):
                                for j in range(t):
                                    block = identity_matrix(SageZZ, size)
                                    standard = identity_matrix(SageZZ, k)
                                    standard[i, j] = 1
                                    block[n:, n:] = rframe.inverse() * standard * rframe
                                    blocks.append(block)
                            if t < k:
                                radical_complement = reduced.lattice_category()([[0] * (k - t) for _ in range(k - t)])
                                r_labels = tuple(radical_complement.module_generating_set())
                                for action in radical_complement.general_linear_group().group_generators():
                                    standard = identity_matrix(SageZZ, k)
                                    for i, point in enumerate(radical_complement.module_generators()):
                                        for j, label in enumerate(r_labels):
                                            standard[t + i, t + j] = int(action(point).to_vector()(label))
                                    block = identity_matrix(SageZZ, size)
                                    block[n:, n:] = rframe.inverse() * standard * rframe
                                    blocks.append(block)
                            for i in range(s, n):
                                for j in range(k):
                                    standard = matrix(SageZZ, n, k, 0)
                                    standard[i, j] = 1
                                    block = identity_matrix(SageZZ, size)
                                    block[:n, n:] = cframe.inverse() * standard
                                    blocks.append(block)
                            if ring is SageQQ:
                                rational = reduced.base_change(reduced.base_ring().fraction_field_map())
                                inclusion_r = ModulesOverIntegralDomains(reduced.base_ring())(reduced).generic_fibre_map()
                                conjugates = []
                                for block in blocks:
                                    conjugated = normalizer * block * normalizer.inverse()
                                    conjugates.append(rational.Isom(rational)(tuple(
                                        rational.linear_combination({label: rational.base_ring()(conjugated[i, j])
                                            for j, label in enumerate(labels) if conjugated[i, j]})
                                        for i in range(size))))
                                action = IntegralStructureAction(rational.O().subgroup(tuple(conjugates)), inclusion_r)
                                from math import lcm
                                bound = lcm(*divisors) ** 2
                                containing = Modules(reduced.base_ring()).Mono(reduced, inclusion_r.codomain())(
                                    lambda label: inclusion_r.codomain().wrap(
                                        rational.scalar_multiple(rational.base_ring().one() / rational.base_ring()(bound),
                                            inclusion_r(reduced.module_generator(label)).underlying_element())))
                                span = OrbitSpanConstructions(reduced.base_ring())(
                                    action.ambient_action_functor(), inclusion_r, containing=containing)
                                generators = _conjugated(reduction.isometry, tuple(_integral_generator_set(
                                    reduced, action.with_stable_containing_bound(span.inclusion()))))
                            else:
                                generators = _conjugated(reduction.isometry, tuple(
                                    _isometry_from_row_matrix(reduced, reduced,
                                        normalizer * block * normalizer.inverse()) for block in blocks))
                            generators = generators or (ambient.O().identity(),)
                            if not all(generator(arrow(point)) == arrow(point)
                                       for generator in generators
                                       for point in arrow.domain().module_generators()):
                                raise ArithmeticError("the primitive projected plane stabilizer failed its frame certificate")
                            return finite_ordered_set(tuple(generators))
                        raise ValueError("based isotropic stabilizers require a nonzero quotient projection")
                # Pointwise fixation of the marked finite-index frame is
                # equivalent to pointwise fixation of its primitive
                # saturation, since the ambient module is torsion-free.
                inclusion = saturation.inclusion()
                basis = tuple(inclusion(point) for point in inclusion.domain().module_generators())
                generators = INDUCTION._stabilizer_generators(basis, ())
        case OrientedIsotropicFlags():
            oriented = marked.functor().oriented_flag()
            generators = _oriented_isotropic_stabilizer_generators(oriented.lattice().O(), oriented)
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
                # Full-radical maximal isotropic subspaces use the common
                # exact parabolic owner below. Filtering a generating set
                # by whether each generator preserves J does not generate
                # Stab(J), and quotient lifts must retain the entire
                # Hom(C,R) kernel. INDUCTION.plane_automorphism_generators
                # constructs all three blocks with integral certificates.
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
                    from sage.matrix.constructor import matrix
                    from sage.rings.integer_ring import ZZ as SageZZ

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
            elif not ambient.is_nondegenerate() and len(dimensions) >= 2 and any(
                all(RECURSION.vector_lies_in_radical(point) for point in basis[:rank])
                and rank == dimensions[-1] + int(ambient.radical().rank())
                - int(ambient.subobject_on((*basis, *(
                    ambient.radical().inclusion()(point)
                    for point in ambient.radical().module_generators()
                ))).rank())
                for rank in dimensions[:-1]
            ):
                # One stage T is precisely J∩rad(L), characteristic in J.
                # Preserve both its preceding radical flag and its
                # subsequent quotient flag in the integral parabolic.
                generators = INDUCTION.plane_automorphism_generators(
                    ambient.primitive_isotropic_subobject(*basis),
                    flag_dimensions=dimensions if len(dimensions) > 2 else (),
                    flag_basis=basis if len(dimensions) > 2 else (),
                )
            elif (not ambient.is_nondegenerate() and len(dimensions) >= 2
                  and (radical_rank := dimensions[-1] + int(ambient.radical().rank())
                       - int(ambient.subobject_on((*basis, *(
                           ambient.radical().inclusion()(point)
                           for point in ambient.radical().module_generators()
                       ))).rank())) > 0
                  and radical_rank < dimensions[-1]
                  and radical_rank not in dimensions):
                # J∩rad(L) is characteristic in J. An omitted stage at
                # exactly its rank imposes no extra restriction: insert
                # that stage for the graph-parabolic construction rather
                # than falling back to nondegenerate Witt extension.
                augmented = tuple(sorted((*dimensions, radical_rank)))
                # Within an interval of the supplied flag, the adapted
                # basis is not ordered by the characteristic intersection.
                # Reorder only within intervals, retaining every given
                # stage while exposing J intersect rad(L) as a prefix.
                previous = 0
                reordered = []
                for rank in dimensions:
                    reordered.extend(sorted(basis[previous:rank], key=lambda point:
                                            not RECURSION.vector_lies_in_radical(point)))
                    previous = rank
                adapted = tuple(reordered)
                if not all(RECURSION.vector_lies_in_radical(point)
                           for point in adapted[:radical_rank]):
                    generators = INDUCTION.plane_automorphism_generators(
                        ambient.primitive_isotropic_subobject(*basis),
                        flag_dimensions=dimensions, flag_basis=basis,
                    )
                else:
                    generators = INDUCTION.plane_automorphism_generators(
                        ambient.primitive_isotropic_subobject(*adapted),
                        flag_dimensions=augmented, flag_basis=adapted,
                    )
            elif (not ambient.is_nondegenerate()
                  and dimensions[0] >= int(ambient.radical().rank())
                  and dimensions[-1] > int(ambient.radical().rank())
                  and all(ambient.primitive_isotropic_subobject(*basis[:dimensions[0]]).inclusion().is_in_image(point)
                          for point in (ambient.radical().inclusion()(x)
                                        for x in ambient.radical().module_generators()))):
                # Each flag term contains the characteristic radical R.
                # Quotienting by R gives a primitive isotropic flag in C;
                # the full stabilizer is its inverse image under the
                # projection O(C+R) -> O(C), with unrestricted GL(R) and
                # Hom(C,R) blocks. This includes the radical-first flag.
                generators = INDUCTION.plane_automorphism_generators(
                    ambient.primitive_isotropic_subobject(*basis),
                    flag_dimensions=dimensions, flag_basis=basis,
                )
            elif (not ambient.is_nondegenerate() and len(dimensions) == 2
                  and dimensions[0] == 1
                  and not RECURSION.vector_lies_in_radical(basis[0])
                  and dimensions[-1] + int(ambient.radical().rank())
                  - int(ambient.subobject_on((*basis, *(
                      ambient.radical().inclusion()(point)
                      for point in ambient.radical().module_generators()
                  ))).rank()) == dimensions[-1] - 1):
                # The rank-one characteristic intersection is transverse
                # to the marked first line; its rank coincides with the
                # marked stage but the subobjects are not equal.
                generators = INDUCTION.plane_automorphism_generators(
                    ambient.primitive_isotropic_subobject(*basis),
                    flag_dimensions=dimensions, flag_basis=basis,
                )
            elif (not ambient.is_nondegenerate() and len(dimensions) > 2
                  and dimensions[0] == 1
                  and not RECURSION.vector_lies_in_radical(basis[0])
                  and all(RECURSION.vector_lies_in_radical(point) for point in basis[1:])):
                # All later stages extend the projected first line by a
                # flag inside the characteristic radical intersection.
                generators = INDUCTION.plane_automorphism_generators(
                    ambient.primitive_isotropic_subobject(*basis),
                    flag_dimensions=dimensions, flag_basis=basis,
                )
            elif (not ambient.is_nondegenerate() and len(dimensions) >= 2
                  and (split_reduction := indefinite_reduced_presentation(ambient)) is not None
                  and all(RECURSION.vector_lies_in_radical(point)
                          or all(int((~split_reduction.isometry)(point)
                                     .to_vector()(label)) == 0
                                 for label in tuple(split_reduction.reduced.module_generating_set())[
                                     int(split_reduction.reduced.rank()) - int(ambient.radical().rank()):])
                          for point in basis)):
                # A split flag has independent projected and radical
                # constituent flags, but the shear kernel must respect
                # the radical stage present when each projected line enters.
                generators = INDUCTION.plane_automorphism_generators(
                    ambient.primitive_isotropic_subobject(*basis),
                    flag_dimensions=dimensions, flag_basis=basis,
                )
            elif (not ambient.is_nondegenerate() and len(dimensions) > 2
                  and dimensions[-1] + int(ambient.radical().rank())
                  - int(ambient.subobject_on((*basis, *(
                      ambient.radical().inclusion()(point)
                      for point in ambient.radical().module_generators()
                  ))).rank()) > 0):
                # A nonprimitive projected graph with interleaved stages
                # requires the rational graph parabolic and its integral
                # intersection, not the nondegenerate Witt-extension path.
                generators = INDUCTION.plane_automorphism_generators(
                    ambient.primitive_isotropic_subobject(*basis),
                    flag_dimensions=dimensions, flag_basis=basis,
                )
            else:
                generators = INDUCTION._stabilizer_generators(basis, images)
        case _:
            raise ValueError(f"{category} is not a category of marked lattices this port decides: Lattices, a coslice <n>/Lattices, Ar(Lattices), Fun([t], Lattices) for finite [t] or oriented isotropic flags")
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
                    for orbit in MarkedCellOrbitAlgorithm(
                        perfect_complex(reduction.reduced, NONNEGATIVE_SQUARES, retain_facet_orbits=False),
                        signed_square,
                    ).global_orbits()
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
    from sage.matrix.constructor import matrix
    from sage.matrix.special import identity_matrix
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
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        # In odd rank, -id is a determinant-minus-one stabilizer of every
        # moving line, so SO and O have identical moving-line orbits.
        acting_group = None
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

    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        acting_group = None
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
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        acting_group = None
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
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        acting_group = None
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
    if rank == 0 and multiplier >= 1:
        return nonprimitive_isotropic_sublattice_classes(lattice, 0, (), acting_group=acting_group)
    if rank < 1 or multiplier < 1:
        raise ValueError("rank must be nonnegative and saturation multiplier positive")
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        acting_group = None
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
    r"""Complete integral stabilizer of a primitive degenerate saturation.

    The invariant of a degenerate marked submodule S is not its rank
    minus rank(L/rad(L)): the projection of S to L/rad(L) can have
    proper rank, or meet the radical in any proper primitive submodule.
    The maintained degenerate_marking_generators takes the *actual*
    isotropic/radical intersections and projected integral glue, and
    computes the full setwise stabilizer by the integral parabolic
    action, including its congruence inverse image. Reuse that
    complete owner for every properly degenerate marking, rather
    than guessing a full-quotient block from its rank.
    """
    from sage_indefinite_port.indefinite.degenerate_markings import degenerate_marking_generators
    from sage_indefinite_port.indefinite.recursive import _isometry_from_column_matrix

    inclusion = saturated.inclusion()
    lattice = inclusion.codomain()
    if not inclusion.is_primitive():
        raise ValueError("a saturation stabilizer requires a primitive sublattice")
    if lattice.is_nondegenerate():
        raise ValueError("the degenerate saturation stabilizer requires a nonzero ambient radical")
    basis = tuple(inclusion(point) for point in saturated.module_generators())
    matrices = degenerate_marking_generators(lattice, basis, pointwise=False)
    generators = tuple(
        _isometry_from_column_matrix(lattice, lattice, action) for action in matrices
    ) or (lattice.O().identity(),)
    if any(not inclusion.is_in_image(generator(point))
           for generator in generators for point in basis):
        raise ArithmeticError(
            "the integral degenerate stabilizer does not preserve its marked saturation"
        )
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
    if int(inclusion.domain().rank()) == 0:
        if divisors != ():
            raise ValueError("the zero saturation has the unique empty Smith type")
        return (inclusion,)
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
    if rank == 0 and divisors == ():
        # The zero sublattice is primitive, has the unique empty Smith
        # tuple, and is fixed by every subgroup of O(L).
        source = lattice.lattice_category()(0)
        return primitive_marking_classes(ARROWS, source, lattice)
    if rank < 1 or len(divisors) != rank:
        raise ValueError("nonnegative isotropic rank and matching elementary-divisor tuple required")
    if any(d < 1 for d in divisors) or any(right % left for left, right in pairwise(divisors)):
        raise ValueError("elementary divisors must be positive and ordered by divisibility")
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        # In odd ambient rank, central -id_L has determinant -1 and fixes
        # every moving sublattice setwise. Hence every O(L)-orbit on
        # sublattices is already one SO(L)-orbit: for g of determinant -1,
        # (-id_L)g has determinant +1 and acts identically on sublattices.
        # This is an equality of the *moving* orbit sets, not of the
        # stabilizer groups or the based-vector orbit sets.
        acting_group = None
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

    Rank zero is the unique zero sublattice, emitted once.
    The elementary divisors of a finite-index sublattice of its primitive
    saturation satisfy d_1 | ... | d_r. Every such tuple occurs in the
    finite shell with d_r equal to its maximum. The fixed-Smith-type
    classification supplies the exact stabilizer orbits, so these disjoint
    shells exhaust the unrestricted moving quotient without an index bound,
    in the domain where the fixed-Smith-type classifier is complete and
    terminating. No finite cardinality is claimed for the unrestricted set.
    """
    from itertools import count, product

    if rank == 0:
        yield from nonprimitive_isotropic_sublattice_classes(lattice, 0, (), acting_group=acting_group)
        return
    if rank < 0:
        raise ValueError("isotropic sublattice rank must be nonnegative")
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        acting_group = None
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
    if int(saturation_type.rank()) == 0:
        # The zero lattice is nondegenerate, primitive in every L, and
        # its empty Smith type has exactly one moving class.
        return primitive_marking_classes(ARROWS, saturation_type, lattice)
    if any(d < 1 for d in divisors) or any(right % left for left, right in pairwise(divisors)):
        raise ValueError("elementary divisors must be positive and ordered by divisibility")
    if acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        # In odd rank, -id_L is central, has determinant -1, and fixes
        # every moving sublattice. Hence the O(L) and SO(L) orbit
        # partitions agree, also for nondegenerate saturations and their
        # finite-index sublattices.
        acting_group = None
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

    The zero saturation contributes its unique zero sublattice once.
    Every positive-rank finite-index sublattice has a unique saturation and a unique
    elementary-divisor tuple d_1 | ... | d_r. Enumerate these tuples by
    maximum divisor and retain the actual arithmetic saturation stabilizer
    for each stratum. Completeness is relative to the fixed-saturation
    primitive embedding classifier, as in the finite-index owner.
    """
    from itertools import count, product

    if not saturation_type.is_nondegenerate():
        raise ValueError("a nondegenerate saturation type is required")
    rank = int(saturation_type.rank())
    if rank == 0:
        yield from nonprimitive_nondegenerate_sublattice_classes(
            lattice, saturation_type, (), acting_group=acting_group
        )
        return
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
    lattices, ``Ar(Lattices)`` for sublattices, ``Fun([t], Lattices)`` for finite flags,
    ``OrientedIsotropicFlags(Lattices, t)`` for flags with oriented graded pieces. ``source``
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
        case OrientedIsotropicFlags():
            # A full oriented marked object has the ambient lattice as its
            # underlying parent. Reading source.rank() would therefore
            # add the ambient rank as a fictitious isotropic flag rank.
            # Use the actual oriented flag's successive term dimensions.
            dimensions = (
                source.functor().oriented_flag().flag().flag_type().dimensions
                if source in category
                else _isotropic_flag_shape_dimensions(
                    source, lattice, int(category.domain_category().object_set().cardinality()) - 1
                )
            )
            classes = tuple(category(oriented) for oriented in lattice.O().oriented_isotropic_orbit_representatives(FlagType(dimensions)))
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            dimensions = _isotropic_flag_shape_dimensions(
                source, lattice, int(category.domain_category().object_set().cardinality()) - 1
            )
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
            raise ValueError(f"{category} is not a category of marked lattices this port classifies: a coslice <n>/Lattices, Ar(Lattices) or Fun([t], Lattices) for finite [t] or oriented isotropic flags")
    if (category is ARROWS or isinstance(category, CosliceCategory)) and classes and all(
        int(marked.arrow().domain().rank()) == 0 for marked in classes
    ):
        # Every isometry fixes the unique zero embedding (and the unique
        # rank-zero based marking). Splitting this singleton through a
        # finite-character quotient would require computing O(0), which
        # is unnecessary and is not supported by the definite engine.
        acting_group = None
    if (category is ARROWS or isinstance(category, _FiniteSequenceFunctorCategory)) and acting_group is not None and int(lattice.rank()) % 2 and acting_group is lattice.SO():
        # Central -id_L has determinant -1 in odd rank and preserves
        # every moving sublattice and every term of an unoriented flag.
        # Consequently O(L) and SO(L) have identical moving-arrow and
        # unoriented-flag orbit sets, unlike based or oriented markings.
        acting_group = None
    if acting_group is not None and acting_group is not lattice.O():
        classes = tuple(split for marked in classes
                        for split in lattice.O().split_orbit(marked, acting_group))
    return finite_ordered_set(classes)


def _isotropic_flag_shape_dimensions(source: Parent | Morphism, lattice: Parent, expected_terms: int) -> tuple[int, ...]:
    r"""Ranks of the isotropic terms of a source shape, excluding the ambient.

    A shorter diagram encodes only its null terms, while a full marked
    diagram has one final ambient leg. These are different even if the
    ambient Gram form vanishes or belongs to another isometric parent.
    Both oriented and unoriented flag classifiers consume the same shape.
    """
    if isinstance(source, Parent):
        shape_lattices = (source,)
        dimensions = (int(source.rank()),)
    else:
        shape_arrows = (source,) if isinstance(source, Morphism) else _flag_arrows(source)
        if not shape_arrows:
            raise ValueError("an isotropic flag shape requires at least one inclusion")
        terminal = shape_arrows[-1].codomain()
        # A full marked diagram has the same number of arrows as the
        # requested isotropic terms, with its last codomain representing
        # the ambient lattice. It can be a distinct isometric presentation
        # or a pure zero form, so neither parent identity nor nonzero Gram
        # detects this leg. Its ambient rank must equal the target rank.
        # Otherwise, a diagram of the same arrow count may specify one
        # additional isotropic source term, which must not be discarded.
        full_diagram_shape = (
            not isinstance(source, Morphism)
            and len(shape_arrows) == expected_terms
            and int(terminal.rank()) == int(lattice.rank())
        )
        shape_lattices = tuple(arrow.domain() for arrow in shape_arrows) + (
            (() if full_diagram_shape else (terminal,))
        )
        dimensions = (
            tuple(int(arrow.domain().rank()) for arrow in shape_arrows)
            if full_diagram_shape
            else (int(shape_arrows[0].domain().rank()),
                  *(int(arrow.codomain().rank()) for arrow in shape_arrows))
        )
    if any(
        term.b(left, right) != 0
        for term in shape_lattices
        for left in term.module_generators()
        for right in term.module_generators()
    ):
        raise ValueError("nonisotropic flag shapes require a general formed-diagram classification")
    if len(dimensions) != expected_terms:
        raise ValueError("the isotropic source shape must specify exactly the target flag's ranks")
    if dimensions[0] < 1 or any(left >= right for left, right in pairwise(dimensions)):
        raise ValueError("isotropic flag ranks must be positive and strictly increasing")
    return dimensions

def iter_primitive_isotropic_flag_classes(lattice: Parent, dimensions: tuple[int, ...], *, acting_group=None):
    r"""Exhaustive height-shell stream of primitive moving isotropic flags.

    Unlike the finite flag classifier, this permits infinitely many orbits
    when the ambient form is degenerate. Each finite shell is exhausted
    before the next; exact public flag transporters eliminate duplicates.
    Completeness and productivity are conditional on termination of those
    exact transporter decisions for the supplied lattice. Supported
    finite-image arithmetic subgroups split each full-group orbit.
    """
    from itertools import count, product
    from sage_indefinite_port.indefinite.recursive import _coordinate_box, _coordinate_element

    if not dimensions or dimensions[0] < 1 or any(a >= b for a, b in pairwise(dimensions)):
        raise ValueError("isotropic flag ranks must be positive and strictly increasing")
    rank = int(lattice.rank())
    top = dimensions[-1]
    if top > rank:
        return
    if lattice.is_nondegenerate() and top > int(lattice.witt_index()):
        return
    reduction = None
    quotient_witt = None
    if not lattice.is_nondegenerate():
        reduction, nondegenerate_rank, radical_rank, quotient = _degenerate_quotient_and_radical(lattice)
        quotient_witt = int(quotient.witt_index()) if quotient is not None and quotient.is_isotropic() else 0
        if top > radical_rank + quotient_witt:
            return
    if acting_group is lattice.O() or (
        acting_group is not None and rank % 2 and acting_group is lattice.SO()
    ):
        # Central -id preserves every unoriented flag in odd rank.
        acting_group = None
    category = Cat().Mor(FiniteOrdinalCategory(len(dimensions) + 1), LATTICES)
    if lattice.is_nondegenerate():
        # The source-backed isotropic induction already certifies the
        # complete finite orbit set in its nondegenerate algorithmic
        # domain. Do not continue searching height shells indefinitely
        # after its final representative.
        for basis in INDUCTION.isotropic_flag_bases(lattice, dimensions):
            arrows = _flag_marking_from_basis(lattice, basis, dimensions)
            marked = category(arrows[0] if len(arrows) == 1 else arrows)
            if acting_group is None:
                yield marked
            else:
                yield from lattice.O().split_orbit(marked, acting_group)
        return
    if all(lattice.b(x, y) == lattice.base_ring().zero()
           for x in lattice.module_generators() for y in lattice.module_generators()):
        # Every primitive flag in a zero-form lattice is a primitive
        # module flag, and GL(L) is transitive on flags of fixed ranks.
        # The complete full-group quotient is finite, so terminate it.
        basis = tuple(lattice.module_generators())[:top]
        arrows = _flag_marking_from_basis(lattice, basis, dimensions)
        marked = category(arrows[0] if len(arrows) == 1 else arrows)
        if acting_group is None:
            yield marked
        else:
            yield from lattice.O().split_orbit(marked, acting_group)
        return
    if reduction is not None and quotient_witt == 0:
        # A rationally anisotropic nondegenerate quotient forces every
        # totally isotropic subspace into rad(L). The radical is saturated,
        # and its full GL action extends to the split reduced lattice.
        # Thus each primitive flag shape has exactly one full O(L)-orbit.
        reduced = reduction.reduced
        labels = tuple(reduced.module_generating_set())[nondegenerate_rank:]
        basis = tuple(reduction.isometry(reduced.module_generator(label))
                      for label in labels[:top])
        arrows = _flag_marking_from_basis(lattice, basis, dimensions)
        marked = category(arrows[0] if len(arrows) == 1 else arrows)
        if acting_group is None:
            yield marked
        else:
            yield from lattice.O().split_orbit(marked, acting_group)
        return
    seen = []
    zero = lattice.base_ring().zero()
    for height in count(1):
        # The maximum of the coordinate heights must equal this shell.
        # Enumerating bounded tuples afresh would duplicate earlier shells.
        vectors = tuple(_coordinate_element(lattice, coordinates)
                        for bound in range(1, height + 1)
                        for coordinates in _coordinate_box(lattice, bound))
        for basis in product(vectors, repeat=top):
            if max(max((abs(int(point.to_vector()(label)))
                        for label in lattice.module_generating_set()), default=0)
                       for point in basis) != height:
                continue
            if any(lattice.b(x, y) != zero for x in basis for y in basis):
                continue
            if any(not lattice.subobject_on(basis[:d]).inclusion().is_primitive()
                   or int(lattice.subobject_on(basis[:d]).rank()) != d
                   for d in dimensions):
                continue
            arrows = _flag_marking_from_basis(lattice, basis, dimensions)
            marked = category(arrows[0] if len(arrows) == 1 else arrows)
            if any(isometry(category, marked, previous) is not None for previous in seen):
                continue
            seen.append(marked)
            if acting_group is None:
                yield marked
            else:
                yield from lattice.O().split_orbit(marked, acting_group)


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


def isotropic_stabilizer_right_cosets(category: Category, marked: Parent) -> Parent:
    r"""A transversal of ``G/H`` for a primitive isotropic marking.

    CombinedAlgorithms.h:1565-1609 INDEF_FORM_RightCosets_IsotropicKstuff_Reduced,
    specialized at :2461-2471. ``G`` is the rational group extending the
    source's stabilizer of ``I^perp`` and ``H = G`` meets ``GL(L)``.
    """
    match category:
        case CosliceCategory():
            # Based markings use the pointwise rational stabilizer:
            # its permitted GL(I) image is the identity, not the
            # setwise GL(I) stabilizer of an unbased sublattice. This
            # is the trivial-Gamma case of the pinned isotropic
            # stabilizer construction (IndefiniteForms.g:1710-1723).
            arrow = marked.arrow()
            ambient = arrow.codomain()
            if not arrow.is_primitive():
                raise ValueError("based isotropic right cosets require a primitive integral marking")
            if int(arrow.domain().rank()) == 0:
                return finite_ordered_set((ambient.O().identity(),))
            basis = tuple(arrow(point) for point in arrow.domain().module_generators())
            if any(ambient.b(first, second) != ambient.base_ring().zero()
                   for first in basis for second in basis):
                raise ValueError("a based isotropic right-coset marking must have zero Gram form")
            images = ()
        case _ if category is ARROWS:
            sublattice = _isotropic_sublattice(marked.arrow())
            if int(sublattice.inclusion().domain().rank()) == 0:
                # For I=0, I-perp=L. The rational extension group already
                # preserves the entire integral lattice, so G=G_L=O(L).
                # There is one right coset, represented by the identity;
                # the rank-zero general linear engine is not needed.
                ambient = marked.arrow().codomain()
                return finite_ordered_set((ambient.O().identity(),))
            basis, images = INDUCTION.plane_vectors(sublattice), INDUCTION.sublattice_action_images(sublattice)
        case OrientedIsotropicFlags():
            _dimensions, basis, images = _oriented_flag_action(marked)
        case _FiniteSequenceFunctorCategory() if category.codomain_category() is LATTICES:
            _dimensions, basis, images = _flag_action(marked)
        case _:
            raise ValueError(f"{category} is not a based isotropic embedding, Ar(Lattices), or Fun([t], Lattices)")
    cosets = RECURSION.isotropic_stabilizer_action(basis, images).right_cosets()
    representative = cosets.representative_map()
    # The selected orbit is G/G_L, with G_L acting on the RIGHT.
    # Inverting its section produces left-coset representatives instead;
    # for a nonnormal stabilizer they can duplicate right cosets.
    return finite_ordered_set(tuple(representative(coset) for coset in cosets))


__all__ = [
    "automorphism_group_generators",
    "isometry",
    "isotropic_stabilizer_right_cosets",
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
    "iter_primitive_isotropic_flag_classes",
    "nonprimitive_nondegenerate_sublattice_classes",
    "iter_all_nondegenerate_sublattice_index_classes",
]
