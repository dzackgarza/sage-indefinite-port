"""Odd full-group vector classes transported through even doubling, not stable kernels."""
from tests.acceptance.consumer import lattice


def test_odd_binary_split_full_group_orbits_transfer_back_to_original_form():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import odd_primitive_pointed_classes_via_doubling

    ambient = lattice([[1, 0], [0, -1]])
    source = Lattices(ZZ)([[1]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = odd_primitive_pointed_classes_via_doubling(category, ambient, 1)
    assert classes.cardinality() == cardinal(1)
    for marked in classes:
        point = marked.arrow()(source.module_generators()[0])
        assert point.q() == ZZ(1)
        assert point.content() == ZZ(1)
    doubled, _ = ambient.doubled_form()
    assert ambient.discriminant_group().cardinality() == cardinal(1)
    assert doubled.discriminant_group().cardinality() == cardinal(4)


def test_public_odd_split_primitive_markings_use_even_double():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[1, 0], [0, -1]])
    source = Lattices(ZZ)([[1]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = primitive_marking_classes(category, source, ambient)
    assert classes.cardinality() == cardinal(1)
    assert all(marked.arrow()(source.module_generators()[0]).q() == ZZ(1) for marked in classes)


def test_odd_split_doubling_preserves_zero_and_negative_square():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[1, 0], [0, -1]])
    for square in (0, -1):
        source = Lattices(ZZ)([[square]])
        category = Lattices(ZZ).CosliceUnder(source)
        classes = primitive_marking_classes(category, source, ambient)
        assert classes.cardinality() == cardinal(1)
        assert all(marked.arrow()(source.module_generators()[0]).q() == ZZ(square) for marked in classes)


def test_odd_split_nonprimitive_norm_four_orbit_preserves_original_square():
    from sage_extensions.all import Lattices, ZZ, cardinal
    from sage_indefinite_port.indefinite.marked_lattices import pointed_nonzero_norm_classes

    ambient = lattice([[1, 0], [0, -1]])
    source = Lattices(ZZ)([[4]])
    category = Lattices(ZZ).CosliceUnder(source)
    classes = pointed_nonzero_norm_classes(category, ambient, 4)
    assert classes.cardinality() == cardinal(1)
    marked = tuple(classes)[0]
    point = marked.arrow()(source.module_generators()[0])
    assert point.q() == ZZ(4)
    assert point.content() == ZZ(2)
