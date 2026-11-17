"""Exact low-rational-Witt-index isometry witnesses through public lattice Mor."""
from tests.acceptance.consumer import lattice


def test_low_witt_nonprimitive_isotropic_line_stream_preserves_saturation_index():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_isotropic_line_content_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[2, 0], [0, -4]])
    for content in (2, 3):
        marking = next(iter_isotropic_line_content_classes(ambient, content))
        arrow = marking.arrow()
        vector = arrow(next(iter(arrow.domain().module_generators())))
        assert vector.q() == ZZ(0)
        assert int(vector.content()) == content


def test_rank_four_witt_zero_and_one_changed_basis_witnesses():
    from sage.all import ZZ as SageZZ, matrix
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram in (
        [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]],
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, -4]],
    ):
        original = lattice(gram)
        change = matrix(SageZZ, [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]])
        changed = change * matrix(SageZZ, gram) * change.transpose()
        target = lattice([[int(value) for value in row] for row in changed.rows()])
        witness = RECURSION.isometry(original, target)
        assert witness is not None
        assert all(target.b(witness(x), witness(y)) == original.b(x, y)
                   for x in original.module_generators() for y in original.module_generators())
        assert all((~witness)(witness(x)) == x for x in original.module_generators())


def test_low_witt_nonisometric_same_signature_rejected_exactly():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    first = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    second = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -2, 0], [0, 0, 0, -18]])
    assert RECURSION.isometry(first, second) is None


def test_public_isometry_mor_constructs_low_witt_rank_four_witness():
    from sage.all import ZZ as SageZZ, matrix

    gram = matrix(SageZZ, [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    source = lattice([[int(value) for value in row] for row in gram.rows()])
    change = matrix(SageZZ, [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    transformed = change * gram * change.transpose()
    target = lattice([[int(value) for value in row] for row in transformed.rows()])
    witness = source.Isom(target).an_element()
    assert all(target.b(witness(x), witness(y)) == source.b(x, y)
               for x in source.module_generators() for y in source.module_generators())


def test_low_witt_exact_orthogonal_iterator_yields_distinct_integral_maps():
    from itertools import islice
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries

    ambient = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    family = tuple(islice(iter_exact_integral_orthogonal_isometries(ambient), 3))
    assert len(set(family)) == 3
    assert all(all(ambient.b(g(x), g(y)) == ambient.b(x, y)
                   for x in ambient.module_generators() for y in ambient.module_generators())
               for g in family)


def test_witt_zero_pointed_vector_transporter_swaps_positive_coordinates():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    e, f, u, v = tuple(ambient.module_generators())
    witness = RECURSION.pointed_isometry(e, f)
    assert witness is not None and witness(e) == f
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, u, v) for y in (e, f, u, v))


def test_low_witt_public_lazy_fixed_norm_orbit_iterator_returns_exact_vector():
    from sage_extensions.all import Lattices, ZZ

    ambient = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    representative = next(ambient.O().exact_vector_orbit_iterator(2, primitive_only=True))
    assert representative.q() == ZZ(2)
    assert int(representative.content()) == 1


def test_witt_one_primitive_isotropic_vector_lazy_stream():
    from sage_extensions.all import Lattices, ZZ

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[2, 0], [0, -4]])
    representative = next(ambient.O().exact_vector_orbit_iterator(0, primitive_only=True))
    assert representative.q() == ZZ(0)
    assert int(representative.content()) == 1


def test_low_witt_one_primitive_null_pointed_transporter():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[2, 0], [0, -4]])
    e, f, u, v = tuple(ambient.module_generators())
    witness = RECURSION.pointed_isometry(e, f)
    assert witness is not None and witness(e) == f
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, u, v) for y in (e, f, u, v))


def test_low_witt_one_public_lazy_primitive_null_line_orbit_stream():
    from sage_indefinite_port.indefinite.marked_lattices import iter_primitive_isotropic_line_classes
    from sage_extensions.all import Lattices, ZZ

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[2, 0], [0, -4]])
    marked = next(iter_primitive_isotropic_line_classes(ambient))
    arrow = marked.arrow()
    vector = arrow(next(iter(arrow.domain().module_generators())))
    assert vector.q() == ZZ(0)
    assert int(vector.content()) == 1
    assert arrow.is_primitive()


def test_low_witt_one_primitive_null_line_moving_transporter_and_inverse():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[2, 0], [0, -4]])
    e, f, u, v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0]])
    first = ARROWS(source.hom([e], codomain=ambient))
    second = ARROWS(source.hom([f], codomain=ambient))
    witness = isometry(ARROWS, first, second)
    assert witness is not None
    assert witness(e) in (f, -f)
    assert (~witness)(witness(e)) == e
