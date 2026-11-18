"""Exact low-rational-Witt-index isometry witnesses through public lattice Mor."""
from tests.acceptance.consumer import lattice


def test_finite_primitive_classes_reject_impossible_norm_without_eichler_model(monkeypatch):
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    def forbidden(*args, **kwargs):
        raise AssertionError("an impossible norm must not request a rational Eichler envelope")

    monkeypatch.setattr(RECURSION, "_model", forbidden)
    assert RECURSION.primitive_vector_classes(ambient, 1) == ()


def test_finite_primitive_classes_reject_prime_power_obstruction_without_eichler_model(monkeypatch):
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    # q=2(x^2+y^2)-6(z^2+w^2).  Its value ideal is (2),
    # so norm 18 passes the ideal test.  Modulo 9, q=0 forces
    # x=y=0 modulo 3 and then z=w=0 modulo 3: no primitive
    # residue exists.
    ambient = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    def forbidden(*args, **kwargs):
        raise AssertionError("a certified empty primitive fibre must not request an Eichler envelope")

    monkeypatch.setattr(RECURSION, "_model", forbidden)
    assert RECURSION.primitive_vector_classes(ambient, 18) == ()


def test_finite_primitive_classes_reject_mod_five_obstruction_without_eichler_model(monkeypatch):
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    # Modulo 5, q=x^2, and 2 is not a square.  The integral
    # value ideal is (1), so it does not detect this obstruction.
    ambient = lattice([[1, 0, 0, 0], [0, 5, 0, 0], [0, 0, -5, 0], [0, 0, 0, -5]])
    def forbidden(*args, **kwargs):
        raise AssertionError("mod-five nonrepresentability must precede Eichler reduction")

    monkeypatch.setattr(RECURSION, "_model", forbidden)
    assert RECURSION.primitive_vector_classes(ambient, 2) == ()
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2)) == ()


def test_exact_vector_orbit_iterator_terminates_on_certified_empty_and_singleton_null_fibres():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    zero = Lattices(ZZ)(0)
    anisotropic = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    for ambient in (zero, anisotropic):
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0)) == (ambient.zero(),)
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True)) == ()
    assert tuple(RECURSION.exact_vector_orbit_iterator(zero, 2)) == ()


def test_exact_vector_orbit_iterator_terminates_on_norm_ideal_obstruction():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    hyperbolic = lattice([[0, 1], [1, 0]])
    even_indefinite = lattice([[2, 0], [0, -6]])
    for ambient in (hyperbolic, even_indefinite):
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 1)) == ()
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 1, primitive_only=True)) == ()


def test_exact_vector_orbit_iterator_terminates_on_modular_obstruction():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    # q(x,y)=2*x^2-6*y^2 is 0 or 2 modulo 3, never 4=1.
    ambient = lattice([[2, 0], [0, -6]])
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 4)) == ()
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 4, primitive_only=True)) == ()
    # In U, q(x,y)=2xy; every primitive vector has an odd coordinate,
    # so zero modulo 2 is possible and must not be rejected.
    hyperbolic = lattice([[0, 1], [1, 0]])
    assert next(RECURSION.exact_vector_orbit_iterator(hyperbolic, 0, primitive_only=True)).content() == 1


def test_exact_vector_orbit_iterator_detects_prime_power_and_primitive_obstructions():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    sums_of_two_squares = lattice([[1, 0], [0, 1]])
    # 3 is representable modulo 2 but not modulo 4.
    assert tuple(RECURSION.exact_vector_orbit_iterator(sums_of_two_squares, 3)) == ()
    # 4 has integral representations, but every such representation
    # is divisible by 2; the primitive residue obstruction is mod 4.
    assert tuple(RECURSION.exact_vector_orbit_iterator(sums_of_two_squares, 4, primitive_only=True)) == ()
    assert next(RECURSION.exact_vector_orbit_iterator(sums_of_two_squares, 4)).q() == 4


def test_exact_vector_orbit_iterator_exhausts_definite_norm_fibres():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram in ([[1, 0], [0, 1]], [[-1, 0], [0, -1]]):
        ambient = lattice(gram)
        sign = gram[0][0]
        assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, sign))) == 1
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, -sign)) == ()
        assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0))) == 1
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True)) == ()


def test_exact_vector_orbit_iterator_exhausts_primitive_zero_form_orbits():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for rank in (1, 2, 4):
        ambient = lattice([[0] * rank for _ in range(rank)])
        classes = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))
        assert len(classes) == 1
        assert int(classes[0].content()) == 1
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2, primitive_only=True)) == ()


def test_exact_vector_orbit_iterator_zero_form_content_stream_is_complete_and_ordered():
    from itertools import islice
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for rank in (1, 2, 4):
        ambient = lattice([[0] * rank for _ in range(rank)])
        classes = tuple(islice(RECURSION.exact_vector_orbit_iterator(ambient, 0), 6))
        assert [int(point.content()) for point in classes] == list(range(6))
        assert all(point.q() == 0 for point in classes)


def test_zero_rank_exact_integral_orthogonal_stream_has_one_identity():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries

    ambient = Lattices(ZZ)(0)
    stream = iter_exact_integral_orthogonal_isometries(ambient)
    identity = next(stream)
    assert identity.domain() is ambient and identity.codomain() is ambient
    assert identity == ambient.O().identity()
    assert tuple(stream) == ()


def test_exact_integral_orthogonal_stream_terminates_on_definite_lattices():
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries

    for gram in ([[1, 0], [0, 1]], [[-1, 0], [0, -1]], [[2, 0], [0, 6]]):
        ambient = lattice(gram)
        elements = tuple(iter_exact_integral_orthogonal_isometries(ambient))
        expected = 8 if abs(gram[0][0]) == abs(gram[1][1]) else 4
        assert len(elements) == expected
        assert len({tuple(element(next(iter(ambient.module_generators()))).to_vector()(label)
                         for label in ambient.module_generating_set()) for element in elements}) >= 2
        assert all(element.domain() is ambient and element.codomain() is ambient for element in elements)


def test_exact_integral_orthogonal_stream_exhausts_rank_one_degenerate_group():
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries

    ambient = lattice([[0]])
    point = next(iter(ambient.module_generators()))
    elements = tuple(iter_exact_integral_orthogonal_isometries(ambient))
    assert len(elements) == 2
    assert {element(point) for element in elements} == {point, -point}
    assert all(element.domain() is ambient and element.codomain() is ambient for element in elements)


def test_exact_integral_orthogonal_stream_exhausts_split_binary_indefinite_groups():
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries

    for gram in ([[0, 1], [1, 0]], [[0, 3], [3, 0]], [[2, 0], [0, -2]]):
        ambient = lattice(gram)
        elements = tuple(iter_exact_integral_orthogonal_isometries(ambient))
        assert len(elements) == 4
        assert all(element.domain() is ambient and element.codomain() is ambient for element in elements)


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


def test_low_witt_isotropic_line_stream_explicit_full_group():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import iter_isotropic_line_content_classes

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[2, 0], [0, -4]])
    marked = next(iter_isotropic_line_content_classes(ambient, 2, acting_group=ambient.O()))
    arrow = marked.arrow()
    vector = arrow(next(iter(arrow.domain().module_generators())))
    assert vector.q() == ZZ(0)
    assert int(vector.content()) == 2


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


def test_quaternion_norm_one_unit_gives_explicit_integral_witt_zero_isometry():
    """For B=(-1,3), left multiplication by 2+j has reduced norm one."""
    from sage.all import ZZ as SageZZ, matrix

    ambient = lattice([[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, -6]])
    # Coordinates are (1,i,j,ij); the row matrix is the transpose of
    # the usual column-action matrix for left multiplication by 2+j.
    action = matrix(SageZZ, [[2, 0, 1, 0], [0, 2, 0, -1],
                             [3, 0, 2, 0], [0, -3, 0, 2]])
    gram = matrix(SageZZ, [[2, 0, 0, 0], [0, 2, 0, 0],
                            [0, 0, -6, 0], [0, 0, 0, -6]])
    assert abs(action.det()) == 1
    assert action * gram * action.transpose() == gram
    basis = tuple(ambient.module_generators())
    witness = ambient.Aut()({label: sum((int(action[i, j]) * basis[j] for j in range(4)), ambient.zero())
                             for i, label in enumerate(ambient.module_generating_set())})
    assert witness in ambient.O()
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y) for x in basis for y in basis)
    assert all((~witness)(witness(x)) == x for x in basis)


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
