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
