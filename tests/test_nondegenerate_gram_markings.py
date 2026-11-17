"""Primitive fixed-source Gram embeddings retain the based source action."""
from tests.acceptance.consumer import lattice


def test_public_primitive_hyperbolic_plane_embedding_classes_in_two_u():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    source = Lattices(ZZ)("U")
    ambient = source + source
    category = Lattices(ZZ).CosliceUnder(source)
    classes = primitive_marking_classes(category, source, ambient)
    assert classes.cardinality() == cardinal(1)
    for marked in classes:
        arrow = marked.arrow()
        assert arrow.domain() is source and arrow.codomain() is ambient
        assert arrow.is_primitive()
        assert all(ambient.b(arrow(x), arrow(y)) == source.b(x, y)
                   for x in source.module_generators() for y in source.module_generators())


def test_split_nondegenerate_based_and_moving_transporters():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    source = Lattices(ZZ)("U")
    ambient = source + source
    e, f, u, v = tuple(ambient.module_generators())
    first_arrow = source.hom([e, f], codomain=ambient)
    second_arrow = source.hom([u, v], codomain=ambient)
    based = Lattices(ZZ).CosliceUnder(source)
    witness = isometry(based, based(first_arrow), based(second_arrow))
    assert witness is not None
    assert witness(e) == u and witness(f) == v
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f, u, v) for y in (e, f, u, v))
    moving = isometry(ARROWS, ARROWS(first_arrow), ARROWS(second_arrow))
    assert moving is not None
    assert all(second_arrow.is_in_image(moving(x)) for x in (e, f))


def test_split_nondegenerate_transporter_uses_unimodular_noncoordinate_frame():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    source = Lattices(ZZ)("U")
    ambient = source + source
    e, f, u, v = tuple(ambient.module_generators())
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e + u, f], codomain=ambient))
    second = category(source.hom([u, v], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(e + u) == u and witness(f) == v
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f, u, v) for y in (e, f, u, v))


def test_glued_nondegenerate_marking_has_identity_without_fictitious_split():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, 2]])
    category = Lattices(ZZ).CosliceUnder(source)
    marked = category(source.hom([e + f, u + v], codomain=ambient))
    witness = isometry(category, marked, marked)
    assert witness is not None
    assert all(witness(x) == x for x in (e, f, u, v))


def test_split_nondegenerate_based_and_moving_stabilizer_generators():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    source = Lattices(ZZ)("U")
    ambient = source + source
    e, f, u, v = tuple(ambient.module_generators())
    arrow = source.hom([e, f], codomain=ambient)
    based = Lattices(ZZ).CosliceUnder(source)
    fixed = tuple(automorphism_group_generators(based, based(arrow)))
    moving = tuple(automorphism_group_generators(ARROWS, ARROWS(arrow)))
    assert all(g(e) == e and g(f) == f for g in fixed)
    assert all(arrow.is_in_image(g(e)) and arrow.is_in_image(g(f)) for g in moving)
    assert ambient.O().subgroup(fixed).cardinality() == cardinal(4)
    assert ambient.O().subgroup(moving).cardinality() == cardinal(16)


def test_split_nondegenerate_noncoordinate_based_stabilizer_is_full_complement_group():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import automorphism_group_generators

    source = Lattices(ZZ)("U")
    ambient = source + source
    e, f, u, v = tuple(ambient.module_generators())
    arrow = source.hom([e + u, f], codomain=ambient)
    category = Lattices(ZZ).CosliceUnder(source)
    family = tuple(automorphism_group_generators(category, category(arrow)))
    assert all(g(e + u) == e + u and g(f) == f for g in family)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y) for g in family for x in (e, f, u, v) for y in (e, f, u, v))
    assert ambient.O().subgroup(family).cardinality() == cardinal(4)


def test_glued_rank_one_hyperbolic_line_has_exact_setwise_and_based_stabilizers():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U")
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2]])
    arrow = source.hom([e + f], codomain=ambient)
    based = Lattices(ZZ).CosliceUnder(source)
    fixed = tuple(automorphism_group_generators(based, based(arrow)))
    moving = tuple(automorphism_group_generators(ARROWS, ARROWS(arrow)))
    assert all(g(e + f) == e + f for g in fixed)
    assert all(g(e + f) in (e + f, -e - f) for g in moving)
    assert ambient.O().subgroup(fixed).cardinality() == cardinal(2)
    assert ambient.O().subgroup(moving).cardinality() == cardinal(4)


def test_public_unbased_primitive_nondegenerate_embedding_classes():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, primitive_marking_classes

    source = Lattices(ZZ)("U")
    ambient = source + source
    classes = primitive_marking_classes(ARROWS, source, ambient)
    assert classes.cardinality() == cardinal(1)
    assert all(marking.arrow().is_primitive() and marking.arrow().codomain() is ambient for marking in classes)


def test_glued_rank_two_definite_source_and_complement_have_exact_stabilizers():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, 2]])
    arrow = source.hom([e + f, u + v], codomain=ambient)
    category = Lattices(ZZ).CosliceUnder(source)
    fixed = tuple(automorphism_group_generators(category, category(arrow)))
    moving = tuple(automorphism_group_generators(ARROWS, ARROWS(arrow)))
    assert all(g(e + f) == e + f and g(u + v) == u + v for g in fixed)
    assert all(arrow.is_in_image(g(e + f)) and arrow.is_in_image(g(u + v)) for g in moving)
    assert ambient.O().subgroup(fixed).cardinality() == cardinal(4)
    assert ambient.O().subgroup(moving).cardinality() == cardinal(32)


def test_glued_rank_two_based_transporter_exchanges_positive_generators():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, 2]])
    first_arrow = source.hom([e + f, u + v], codomain=ambient)
    second_arrow = source.hom([u + v, e + f], codomain=ambient)
    category = Lattices(ZZ).CosliceUnder(source)
    witness = isometry(category, category(first_arrow), category(second_arrow))
    assert witness is not None
    assert witness(e + f) == u + v and witness(u + v) == e + f
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f, u, v) for y in (e, f, u, v))
    moving = isometry(ARROWS, ARROWS(first_arrow), ARROWS(second_arrow))
    assert moving is not None
    assert all(second_arrow.is_in_image(moving(x)) for x in (e + f, u + v))


def test_glue_extension_rejects_integral_nonunimodular_rational_isometry():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.lattice_embeddings import _finite_overlattice_isometry_extension

    source_lattice = Lattices(ZZ)([[2, 0], [0, 2]])
    target_lattice = Lattices(ZZ)([[1, 0], [0, 1]])
    u, v = tuple(target_lattice.module_generators())
    source_inclusion = source_lattice.subobject_on(source_lattice.module_generators()).inclusion()
    target_inclusion = target_lattice.subobject_on([u + v, u - v]).inclusion()
    source_module, target_module = source_inclusion.domain(), target_inclusion.domain()
    isometry = source_module.Isom(target_module)(list(target_module.module_generators()))
    assert _finite_overlattice_isometry_extension(source_inclusion, target_inclusion, isometry) is None


def test_odd_index_three_glue_rejects_independent_factor_signs():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)([[1, 0], [0, 2]])
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[3]])
    arrow = source.hom([e + f], codomain=ambient)
    saturated = ambient.subobject_on([e + f])
    complement = ambient.orthogonal_complement(saturated)
    combined = ambient.subobject_on([e + f, *(complement.inclusion()(x) for x in complement.module_generators())])
    assert combined.inclusion().cokernel().cardinality() == cardinal(3)
    based = Lattices(ZZ).CosliceUnder(source)
    fixed = tuple(automorphism_group_generators(based, based(arrow)))
    moving = tuple(automorphism_group_generators(ARROWS, ARROWS(arrow)))
    assert ambient.O().subgroup(fixed).cardinality() == cardinal(1)
    assert ambient.O().subgroup(moving).cardinality() == cardinal(2)
    assert all(arrow.is_in_image(g(e + f)) for g in moving)


def test_odd_index_three_glue_based_transporter_has_integral_inverse():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = Lattices(ZZ)([[1, 0], [0, 2]])
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[3]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e + f], codomain=ambient))
    second = category(source.hom([-e - f], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(e + f) == -e - f
    assert (~witness)(-e - f) == e + f
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f) for y in (e, f))


def test_glued_positive_line_in_two_u_has_full_indefinite_complement_kernel():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import _nondegenerate_stabilizer_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2]])
    arrow = source.hom([e + f], codomain=ambient)
    generators = tuple(_nondegenerate_stabilizer_generators(arrow, based=True))
    assert all(g(e + f) == e + f for g in generators)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in generators for x in (e, f, u, v) for y in (e, f, u, v))
    assert any(g(u) != u or g(v) != v for g in generators)


def test_glued_positive_line_in_two_u_setwise_stabilizer_includes_negation():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import _nondegenerate_stabilizer_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2]])
    arrow = source.hom([e + f], codomain=ambient)
    generators = tuple(_nondegenerate_stabilizer_generators(arrow, based=False))
    assert all(g(e + f) in (e + f, -e - f) for g in generators)
    assert any(g(e + f) == -e - f for g in generators)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in generators for x in (e, f, u, v) for y in (e, f, u, v))


def test_glued_indefinite_rank_two_source_and_complement_use_full_preimage():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, -2]])
    arrow = source.hom([e + f, u - v], codomain=ambient)
    category = Lattices(ZZ).CosliceUnder(source)
    fixed = tuple(automorphism_group_generators(category, category(arrow)))
    moving = tuple(automorphism_group_generators(ARROWS, ARROWS(arrow)))
    assert all(g(e + f) == e + f and g(u - v) == u - v for g in fixed)
    assert all(arrow.is_in_image(g(e + f)) and arrow.is_in_image(g(u - v)) for g in moving)
    assert ambient.O().subgroup(fixed).cardinality() == cardinal(4)
    assert ambient.O().subgroup(moving).cardinality() == cardinal(16)


def test_generic_finite_overlattice_isometry_extension_requires_both_integral_directions():
    from sage_extensions.all import Lattices, ZZ
    from sage_extensions.categories.lattice_embeddings import _finite_overlattice_isometry_extension

    plane = Lattices(ZZ)("U")
    e, f = tuple(plane.module_generators())
    scaled = plane.twist(2)
    e2, f2 = tuple(scaled.module_generators())
    inclusion = scaled.Emb(plane)([2 * e, f])
    negation = scaled.O()(lambda label: -scaled.module_generator(label))
    extension = _finite_overlattice_isometry_extension(inclusion, inclusion, negation)
    assert extension is not None and extension(e) == -e and extension(f) == -f
    swap = scaled.O()({0: f2, 1: e2})
    assert _finite_overlattice_isometry_extension(inclusion, inclusion, swap) is None


def test_glued_indefinite_rank_two_based_transporter_uses_finite_glue_orbit():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, -2]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e + f, u - v], codomain=ambient))
    second = category(source.hom([e + f, v - u], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(e + f) == e + f and witness(u - v) == v - u
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, u, v) for y in (e, f, u, v))


def test_odd_glue_transporter_requires_nonidentity_finite_orbit_correction():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import _glued_nondegenerate_embedding_isometry

    ambient = Lattices(ZZ)([[1, 0], [0, 2]])
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[3]])
    first = source.hom([e + f], codomain=ambient)
    second = source.hom([-e - f], codomain=ambient)
    witness = _glued_nondegenerate_embedding_isometry(first, second, based=True)
    assert witness is not None
    assert witness(e + f) == -e - f
    assert (~witness)(-e - f) == e + f
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in (e, f) for y in (e, f))


def test_finite_glue_orbit_transporter_rejects_incompatible_overlattice_indices():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import _glued_nondegenerate_embedding_isometry

    source_lattice = Lattices(ZZ)([[2, 0], [0, 2]])
    target_lattice = Lattices(ZZ)([[1, 0], [0, 1]])
    source = Lattices(ZZ)([[2]])
    e, _f = tuple(source_lattice.module_generators())
    u, v = tuple(target_lattice.module_generators())
    first = source.hom([e], codomain=source_lattice)
    second = source.hom([u + v], codomain=target_lattice)
    assert _glued_nondegenerate_embedding_isometry(first, second, based=True) is None


def test_glued_pell_binary_factors_return_infinite_stabilizer_generators():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U").twist(2)
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, -4]])
    arrow = source.hom([e + f, u - v], codomain=ambient)
    category = Lattices(ZZ).CosliceUnder(source)
    fixed = tuple(automorphism_group_generators(category, category(arrow)))
    assert all(g(e + f) == e + f and g(u - v) == u - v for g in fixed)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in fixed for x in (e, f, u, v) for y in (e, f, u, v))
    assert any(g(e - f) not in (e - f, -e + f) for g in fixed)
    labels = tuple(ambient.module_generating_set())
    def complement_trace(g):
        return int(g(e - f).to_vector()(labels[0]) + g(u + v).to_vector()(labels[2]))
    assert any(abs(complement_trace(g * h)) > 2 for g in fixed for h in fixed)


def test_glued_infinite_pell_complement_transporter_requires_glue_correction():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import isometry

    ambient = Lattices(ZZ)([[1, 0], [0, 2]]) + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[3, 0], [0, 2]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([e + f, u + v], codomain=ambient))
    second = category(source.hom([-e - f, u + v], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(e + f) == -e - f and witness(u + v) == u + v
    assert (~witness)(-e - f) == e + f
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, u, v) for y in (e, f, u, v))


def test_public_nondegenerate_embedding_transporter_rejects_distinct_glue_indices():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    source_lattice = Lattices(ZZ)([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    target_lattice = Lattices(ZZ)([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    e, _f, g = tuple(source_lattice.module_generators())
    u, v, w = tuple(target_lattice.module_generators())
    source = Lattices(ZZ)([[2, 0], [0, 2]])
    first = ARROWS(source.hom([e, g], codomain=source_lattice))
    second = ARROWS(source.hom([u + v, w], codomain=target_lattice))
    assert isometry(ARROWS, first, second) is None


def test_public_nondegenerate_glued_transporter_between_distinct_ambient_presentations():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    source_lattice = Lattices(ZZ)([[1, 0, 0], [0, 2, 0], [0, 0, 2]])
    target_lattice = Lattices(ZZ)([[2, 0, 0], [0, 1, 0], [0, 0, 2]])
    e, f, g = tuple(source_lattice.module_generators())
    u, v, w = tuple(target_lattice.module_generators())
    source = Lattices(ZZ)([[3, 0], [0, 2]])
    first = ARROWS(source.hom([e + f, g], codomain=source_lattice))
    second = ARROWS(source.hom([u + v, w], codomain=target_lattice))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert second.arrow().is_in_image(witness(e + f))
    assert second.arrow().is_in_image(witness(g))
    assert all(target_lattice.b(witness(x), witness(y)) == source_lattice.b(x, y)
               for x in (e, f, g) for y in (e, f, g))
    assert all((~witness)(witness(x)) == x for x in (e, f, g))


def test_nonprimitive_nondegenerate_rank_two_moving_transport_and_stabilizer():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 2], [2, 0]])
    first = ARROWS(source.hom([2 * e, f], codomain=ambient))
    second = ARROWS(source.hom([e, 2 * f], codomain=ambient))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert all(second.arrow().is_in_image(witness(first.arrow()(x))) for x in source.module_generators())
    stabilizer = tuple(automorphism_group_generators(ARROWS, first))
    assert all(first.arrow().is_in_image(g(2 * e)) and first.arrow().is_in_image(g(f)) for g in stabilizer)
    assert ambient.O().subgroup(stabilizer).cardinality() == cardinal(8)


def test_nonprimitive_nondegenerate_based_markings_distinguish_source_coordinates():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import isometry, automorphism_group_generators

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 2], [2, 0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([2 * e, f], codomain=ambient))
    second = category(source.hom([2 * f, e], codomain=ambient))
    third = category(source.hom([e, 2 * f], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None
    assert witness(2 * e) == 2 * f and witness(f) == e
    assert isometry(category, first, third) is None
    fixed = tuple(automorphism_group_generators(category, first))
    assert all(g(2 * e) == 2 * e and g(f) == f for g in fixed)
    assert ambient.O().subgroup(fixed).cardinality() == cardinal(4)


def test_public_nonprimitive_nondegenerate_classes_from_primitive_saturation():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import nonprimitive_nondegenerate_sublattice_classes

    source = Lattices(ZZ)("U")
    ambient = source + source
    classes = nonprimitive_nondegenerate_sublattice_classes(ambient, source, (1, 2))
    assert classes.cardinality() == cardinal(2)
    assert all(marked.arrow().codomain() is ambient for marked in classes)


def test_unrestricted_nondegenerate_saturation_index_iterator_retains_two_index_two_orbits():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import (
        ARROWS, isometry, iter_all_nondegenerate_sublattice_index_classes,
    )

    source = Lattices(ZZ)("U")
    ambient = source + source
    classes = tuple(islice(iter_all_nondegenerate_sublattice_index_classes(ambient, source), 3))
    assert len(classes) == 3
    indices = []
    for marked in classes:
        arrow = marked.arrow()
        saturation = arrow.saturation()
        embedded = saturation.inclusion().domain().subobject_on(tuple(
            saturation.inclusion().lift(arrow(x)) for x in arrow.domain().module_generators()
        )).inclusion()
        indices.append(int(embedded.cokernel().cardinality()))
    assert indices == [1, 2, 2]
    assert isometry(ARROWS, classes[1], classes[2]) is None


def test_unrestricted_rank_one_nondegenerate_index_iterator_preserves_saturation():
    from itertools import islice
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_all_nondegenerate_sublattice_index_classes

    ambient = Lattices(ZZ)("U")
    source = Lattices(ZZ)([[2]])
    classes = tuple(islice(iter_all_nondegenerate_sublattice_index_classes(ambient, source), 3))
    assert len(classes) == 3
    for content, marked in enumerate(classes, 1):
        arrow = marked.arrow()
        vector = arrow(next(iter(arrow.domain().module_generators())))
        assert int(vector.content()) == content
        assert int(vector.q()) == 2 * content * content
