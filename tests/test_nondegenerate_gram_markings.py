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
