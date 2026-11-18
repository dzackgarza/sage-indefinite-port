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


def test_definite_norm_orbits_use_generator_action_not_pointed_isometry(monkeypatch):
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[1, 0], [0, 2]])
    def forbidden(*args, **kwargs):
        raise AssertionError("finite definite orbits must use exact generator actions")

    monkeypatch.setattr(RECURSION, "pointed_isometry", forbidden)
    # The four solutions (+/-1,+/-1) of x^2+2y^2=3 form one
    # O(L)-orbit.  The norm-2 fibre is another single orbit.
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 3))) == 1
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2))) == 1
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 3, primitive_only=True))) == 1


def test_definite_norm_orbits_respect_explicit_generator_labels():
    from sage_extensions.all import Lattices, ZZ
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = Lattices(ZZ)([[1, 0], [0, 2]], module_generators=("right", "left"))
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 3))
    assert len(representatives) == 1
    assert int(representatives[0].q()) == 3


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


def test_exact_integral_orthogonal_stream_zero_form_general_linear(monkeypatch):
    from itertools import islice
    from sage_indefinite_port.indefinite import recursive

    for rank in (2, 3):
        ambient = lattice([[0] * rank for _ in range(rank)])
        def forbidden(*args, **kwargs):
            raise AssertionError("zero-form group must use the maintained GL generators")
        monkeypatch.setattr(recursive, "_fair_unimodular_matrices", forbidden)
        elements = tuple(islice(recursive.iter_exact_integral_orthogonal_isometries(ambient), 12))
        matrices = tuple(recursive._column_matrix(g) for g in elements)
        assert len({tuple(g.list()) for g in matrices}) == 12
        assert all(abs(g.det()) == 1 for g in matrices)
        assert all(g.domain() is ambient and g.codomain() is ambient for g in elements)


def test_exact_integral_orthogonal_stream_higher_rank_radical_generators(monkeypatch):
    from itertools import islice
    from sage_indefinite_port.indefinite import recursive

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    def forbidden(*args, **kwargs):
        raise AssertionError("degenerate enumeration must not traverse ambient GL_n")
    monkeypatch.setattr(recursive, "_fair_unimodular_matrices", forbidden)
    family = tuple(islice(recursive.iter_exact_integral_orthogonal_isometries(ambient), 12))
    matrices = tuple(recursive._column_matrix(g) for g in family)
    assert len({tuple(g.list()) for g in matrices}) == 12
    assert all(abs(g.det()) == 1 for g in matrices)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in family for x in ambient.module_generators() for y in ambient.module_generators())
    from sage.all import ZZ as SageZZ, matrix
    original = matrix(SageZZ, [[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    change = matrix(SageZZ, [[1, 0, 0], [1, 1, 0], [0, 1, 1]])
    changed = change * original * change.transpose()
    twisted = lattice([[int(value) for value in row] for row in changed.rows()])
    transported = tuple(islice(recursive.iter_exact_integral_orthogonal_isometries(twisted), 12))
    assert len({tuple(recursive._column_matrix(g).list()) for g in transported}) == 12
    assert all(twisted.b(g(x), g(y)) == twisted.b(x, y)
               for g in transported for x in twisted.module_generators() for y in twisted.module_generators())
    mixed = lattice([[2, 0, 0, 0], [0, -6, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    mixed_family = tuple(islice(recursive.iter_exact_integral_orthogonal_isometries(mixed), 12))
    assert len({tuple(recursive._column_matrix(g).list()) for g in mixed_family}) == 12
    assert all(mixed.b(g(x), g(y)) == mixed.b(x, y)
               for g in mixed_family for x in mixed.module_generators() for y in mixed.module_generators())


def test_exact_integral_orthogonal_stream_exhausts_split_binary_indefinite_groups():
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries

    for gram in ([[0, 1], [1, 0]], [[0, 3], [3, 0]], [[2, 0], [0, -2]]):
        ambient = lattice(gram)
        elements = tuple(iter_exact_integral_orthogonal_isometries(ambient))
        assert len(elements) == 4
        assert all(element.domain() is ambient and element.codomain() is ambient for element in elements)


def test_exact_vector_orbit_iterator_split_binary_finite_and_content_classes():
    from itertools import islice
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram in ([[0, 1], [1, 0]], [[2, 0], [0, -2]]):
        ambient = lattice(gram)
        primitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))
        assert len(primitive) == 1
        classes = tuple(islice(RECURSION.exact_vector_orbit_iterator(ambient, 0), 5))
        assert [int(point.content()) for point in classes] == [0, 1, 2, 3, 4]
        assert all(point.q() == 0 for point in classes)
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 1)) == ()
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2, primitive_only=True))


def test_nonsplit_binary_norm_orbit_iterator_terminates_with_exact_pell_classes():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[1, 0], [0, -2]])
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 3)) == ()
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True)) == ()
    classes = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8))
    assert classes
    assert all(int(point.q()) == 8 for point in classes)
    assert all(int(point.content()) == 2 for point in classes)
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8, primitive_only=True)) == ()
    negative = tuple(RECURSION.exact_vector_orbit_iterator(ambient, -1))
    assert negative and all(int(point.q()) == -1 for point in negative)
    primitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 1, primitive_only=True))
    assert primitive and all(int(point.q()) == 1 and int(point.content()) == 1 for point in primitive)
    from sage.all import ZZ as SageZZ, matrix
    gram = matrix(SageZZ, [[1, 0], [0, -2]])
    change = matrix(SageZZ, [[1, 1], [0, 1]])
    transformed = change * gram * change.transpose()
    changed = lattice([[int(value) for value in row] for row in transformed.rows()])
    changed_classes = tuple(RECURSION.exact_vector_orbit_iterator(changed, 8))
    assert len(changed_classes) == len(classes)
    assert all(int(point.q()) == 8 and int(point.content()) == 2 for point in changed_classes)
    even = lattice([[2, 0], [0, -6]])
    assert tuple(RECURSION.exact_vector_orbit_iterator(even, 1)) == ()
    represented = tuple(RECURSION.exact_vector_orbit_iterator(even, 2))
    assert represented and all(int(point.q()) == 2 for point in represented)
    assert tuple(RECURSION.exact_vector_orbit_iterator(even, 0, primitive_only=True)) == ()


def test_exact_vector_orbit_iterator_binary_rank_one_radical_classes():
    from itertools import islice
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram in ([[2, 0], [0, 0]], [[0, 0], [0, -6]]):
        ambient = lattice(gram)
        null_classes = tuple(islice(RECURSION.exact_vector_orbit_iterator(ambient, 0), 5))
        assert [int(point.content()) for point in null_classes] == [0, 1, 2, 3, 4]
        assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))) == 1
    ambient = lattice([[2, 0], [0, 0]])
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8))) == 2
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8, primitive_only=True))) == 1
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2))) == 1


def test_exact_vector_orbit_iterator_rank_one_quotient_higher_radical():
    from itertools import islice
    from sage.all import ZZ as SageZZ, matrix
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for rank in (3, 4):
        gram = [[2 if i == j == 0 else 0 for j in range(rank)] for i in range(rank)]
        original = matrix(SageZZ, gram)
        change = matrix.identity(SageZZ, rank)
        change[1, 0] = 1
        changed = change * original * change.transpose()
        ambient = lattice([[int(value) for value in row] for row in changed.rows()])
        nulls = tuple(islice(RECURSION.exact_vector_orbit_iterator(ambient, 0), 5))
        assert [int(v.content()) for v in nulls] == [0, 1, 2, 3, 4]
        assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))) == 1
        # x=6, q=72; gcd(6,y) ranges over 1,2,3,6.
        representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 72))
        assert len(representatives) == 4
        assert all(int(v.q()) == 72 for v in representatives)
        assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 72, primitive_only=True))) == 1
        negative = lattice([[-int(value) for value in row] for row in changed.rows()])
        negative_representatives = tuple(RECURSION.exact_vector_orbit_iterator(negative, -72))
        assert len(negative_representatives) == 4
        assert all(int(v.q()) == -72 for v in negative_representatives)
        assert tuple(RECURSION.exact_vector_orbit_iterator(negative, 72)) == ()


def test_rank_one_quotient_radical_orbit_representatives_are_pairwise_inequivalent():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 72))
    assert len(representatives) == 4
    for i, first in enumerate(representatives):
        for second in representatives[i + 1:]:
            assert RECURSION.pointed_isometry(first, second) is None
    labels = tuple(ambient.module_generating_set())
    point = ambient.linear_combination({labels[0]: 6, labels[1]: 4, labels[2]: 2})
    assert any(RECURSION.pointed_isometry(point, representative) is not None
               for representative in representatives)


def test_rank_one_quotient_large_prime_norm_uses_exact_divisor_classes():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    prime = 10007
    square = 2 * prime * prime
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, square))
    assert len(representatives) == 2
    assert {int(point.content()) for point in representatives} == {1, prime}
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, square, primitive_only=True))) == 1


def test_anisotropic_quotient_degenerate_null_orbits_are_radical_contents():
    from itertools import islice
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram in (
        [[1, 0, 0], [0, -2, 0], [0, 0, 0]],
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
        [[2, 0, 0, 0], [0, 2, 0, 0], [0, 0, -6, 0], [0, 0, 0, 0]],
        [[2, 0, 0, 0, 0], [0, 2, 0, 0, 0], [0, 0, -6, 0, 0],
         [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]],
    ):
        ambient = lattice(gram)
        primitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))
        assert len(primitive) == 1
        assert int(primitive[0].content()) == 1
        nulls = tuple(islice(RECURSION.exact_vector_orbit_iterator(ambient, 0), 6))
        assert [int(v.content()) for v in nulls] == list(range(6))
        assert all(int(v.q()) == 0 for v in nulls)
    from sage.all import ZZ as SageZZ, matrix
    original = matrix(SageZZ, [[2, 0, 0, 0], [0, 2, 0, 0],
                               [0, 0, -6, 0], [0, 0, 0, 0]])
    change = matrix(SageZZ, [[1, 0, 0, 0], [0, 1, 0, 0],
                             [0, 0, 1, 0], [1, 1, 0, 1]])
    transformed = change * original * change.transpose()
    changed = lattice([[int(value) for value in row] for row in transformed.rows()])
    nulls = tuple(islice(RECURSION.exact_vector_orbit_iterator(changed, 0), 5))
    assert [int(v.content()) for v in nulls] == list(range(5))
    assert all(int(v.q()) == 0 for v in nulls)


def test_definite_quotient_degenerate_nonzero_norm_orbits_are_complete():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram, expected in (
        ([[1, 0, 0], [0, 1, 0], [0, 0, 0]], 2),
        ([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], 2),
    ):
        ambient = lattice(gram)
        representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8))
        assert len(representatives) == expected
        assert {int(v.content()) for v in representatives} == {1, 2}
        assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8, primitive_only=True))) == 1
        assert all(int(v.q()) == 8 for v in representatives)
    ambient = lattice([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 25))) == 3
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 25, primitive_only=True))) == 2
    assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, -1)) == ()
    negative = lattice([[-int(value) for value in row] for row in
                        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]])
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(negative, -25))) == 3
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(negative, -25, primitive_only=True))) == 2
    assert tuple(RECURSION.exact_vector_orbit_iterator(negative, 25)) == ()


def test_definite_quotient_degenerate_norm_orbits_under_basis_change():
    from sage.all import ZZ as SageZZ, matrix
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    original = matrix(SageZZ, [[1, 0, 0, 0], [0, 1, 0, 0],
                               [0, 0, 0, 0], [0, 0, 0, 0]])
    change = matrix(SageZZ, [[1, 0, 0, 0], [0, 1, 0, 0],
                             [1, 1, 1, 0], [0, 1, 0, 1]])
    transformed = change * original * change.transpose()
    ambient = lattice([[int(value) for value in row] for row in transformed.rows()])
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 25))
    assert len(representatives) == 3
    for i, first in enumerate(representatives):
        for second in representatives[i + 1:]:
            assert RECURSION.pointed_isometry(first, second) is None
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 25, primitive_only=True))) == 2


def test_indefinite_quotient_degenerate_norm_orbits_are_finite_and_exact():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    for gram in (
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    ):
        ambient = lattice(gram)
        representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8))
        assert representatives
        assert all(int(point.q()) == 8 for point in representatives)
        for i, first in enumerate(representatives):
            for second in representatives[i + 1:]:
                assert RECURSION.pointed_isometry(first, second) is None
        primitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8, primitive_only=True))
        assert primitive
        assert all(int(point.content()) == 1 for point in primitive)
        assert tuple(RECURSION.exact_vector_orbit_iterator(ambient, 1)) == ()
    nonsplit = lattice([[2, 0, 0], [0, -6, 0], [0, 0, 0]])
    classes = tuple(RECURSION.exact_vector_orbit_iterator(nonsplit, 8))
    assert classes and all(int(point.q()) == 8 for point in classes)
    assert tuple(RECURSION.exact_vector_orbit_iterator(nonsplit, 1)) == ()
    for i, first in enumerate(classes):
        for second in classes[i + 1:]:
            assert RECURSION.pointed_isometry(first, second) is None
    lorentzian = lattice([[0, 1, 0, 0], [1, 0, 0, 0],
                         [0, 0, -2, 0], [0, 0, 0, 0]])
    lorentzian_classes = tuple(RECURSION.exact_vector_orbit_iterator(lorentzian, 2))
    assert lorentzian_classes and all(int(point.q()) == 2 for point in lorentzian_classes)
    for i, first in enumerate(lorentzian_classes):
        for second in lorentzian_classes[i + 1:]:
            assert RECURSION.pointed_isometry(first, second) is None


def test_definite_quotient_radical_residue_orbit_transporters():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[1, 0, 0, 0], [0, 1, 0, 0],
                       [0, 0, 0, 0], [0, 0, 0, 0]])
    labels = tuple(ambient.module_generating_set())
    def point(x, y, r, s):
        return ambient.linear_combination(dict(zip(labels, (x, y, r, s), strict=True)))

    first = point(5, 0, 1, 0)
    second = point(5, 0, 2, 0)
    third = point(5, 0, 0, 5)
    witness = RECURSION.pointed_isometry(first, second)
    assert witness is not None and witness(first) == second
    assert RECURSION.pointed_isometry(first, third) is None
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 25))
    assert len(representatives) == 3
    for i, representative in enumerate(representatives):
        for other in representatives[i + 1:]:
            assert RECURSION.pointed_isometry(representative, other) is None
    # The nonprimitive class is represented by a radical residue of
    # content five; a shear removes multiples of the quotient content.
    nonprimitive = point(5, 0, 5, 0)
    assert any(RECURSION.pointed_isometry(nonprimitive, representative) is not None
               for representative in representatives)


def test_exact_integral_orthogonal_stream_anisotropic_binary_uses_exact_generators(monkeypatch):
    from itertools import islice
    from sage_indefinite_port.indefinite import recursive

    ambient = lattice([[2, 0], [0, -6]])
    generators = (ambient.O().identity(),)
    monkeypatch.setattr(recursive, "lorentzian_orthogonal_group_generators", lambda _: generators)
    def forbidden(*args, **kwargs):
        raise AssertionError("anisotropic binary enumeration must not traverse GL_2")
    monkeypatch.setattr(recursive, "_fair_unimodular_matrices", forbidden)
    assert tuple(recursive.iter_exact_integral_orthogonal_isometries(ambient)) == generators


def test_exact_integral_orthogonal_stream_includes_inverse_pell_period(monkeypatch):
    from itertools import islice
    from sage.all import ZZ as SageZZ, matrix
    from sage_indefinite_port.indefinite import recursive

    ambient = lattice([[2, 0], [0, -6]])
    period = recursive._isometry_from_row_matrix(
        ambient, ambient, matrix(SageZZ, [[2, 1], [3, 2]])
    )
    monkeypatch.setattr(recursive, "lorentzian_orthogonal_group_generators", lambda _: (period,))
    family = tuple(islice(recursive.iter_exact_integral_orthogonal_isometries(ambient), 3))
    matrices = {tuple(recursive._column_matrix(g).list()) for g in family}
    forward = recursive._column_matrix(period)
    assert matrices == {
        tuple(matrix(SageZZ, [[1, 0], [0, 1]]).list()),
        tuple(forward.list()),
        tuple(forward.inverse().list()),
    }


def test_exact_integral_orthogonal_stream_anisotropic_binary_actual_edgewalk():
    from itertools import islice
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries, _column_matrix

    ambient = lattice([[2, 0], [0, -6]])
    family = tuple(islice(iter_exact_integral_orthogonal_isometries(ambient), 8))
    assert len({_column_matrix(g).__repr__() for g in family}) == 8
    assert all(g.domain() is ambient and g.codomain() is ambient for g in family)
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in family for x in ambient.module_generators() for y in ambient.module_generators())


def test_exact_integral_orthogonal_stream_lorentzian_higher_rank_generators(monkeypatch):
    from itertools import islice
    from sage.all import ZZ as SageZZ, matrix
    from sage_indefinite_port.indefinite import recursive

    for gram in ([[1, 0, 0], [0, -1, 0], [0, 0, -1]],
                 [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]):
        ambient = lattice(gram)
        identity = ambient.O().identity()
        monkeypatch.setattr(recursive.EichlerRecursion, "orthogonal_group_generators", lambda self, _: (identity,))
        def forbidden(*args, **kwargs):
            raise AssertionError("Lorentzian enumeration must not traverse GL_n")
        monkeypatch.setattr(recursive, "_fair_unimodular_matrices", forbidden)
        family = tuple(recursive.iter_exact_integral_orthogonal_isometries(ambient))
        assert len(family) == 1
        assert tuple(recursive._column_matrix(family[0]).list()) == tuple(matrix.identity(SageZZ, 3).list())


def test_exact_integral_orthogonal_stream_lorentzian_rank_three_real_generators():
    from itertools import islice
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries, _column_matrix

    ambient = lattice([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
    family = tuple(islice(iter_exact_integral_orthogonal_isometries(ambient), 6))
    assert len({tuple(_column_matrix(g).list()) for g in family}) == 6
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in family for x in ambient.module_generators() for y in ambient.module_generators())
    opposite = lattice([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    opposite_family = tuple(islice(iter_exact_integral_orthogonal_isometries(opposite), 6))
    assert len({tuple(_column_matrix(g).list()) for g in opposite_family}) == 6
    assert all(opposite.b(g(x), g(y)) == opposite.b(x, y)
               for g in opposite_family for x in opposite.module_generators() for y in opposite.module_generators())


def test_exact_integral_orthogonal_stream_higher_witt_index_real_generators():
    from itertools import islice
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries, _column_matrix

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    family = tuple(islice(iter_exact_integral_orthogonal_isometries(ambient), 8))
    assert len({tuple(_column_matrix(g).list()) for g in family}) == 8
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for g in family for x in ambient.module_generators() for y in ambient.module_generators())


def test_lorentzian_fixed_norm_iterator_uses_complete_perfect_cell_orbits():
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    ambient = lattice([[0, 1, 0], [1, 0, 0], [0, 0, -2]])
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2))
    assert representatives
    assert all(int(point.q()) == 2 for point in representatives)
    for i, first in enumerate(representatives):
        for second in representatives[i + 1:]:
            assert RECURSION.pointed_isometry(first, second) is None
    primitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2, primitive_only=True))
    assert len(primitive) == len(representatives)
    nonprimitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8))
    assert nonprimitive and all(int(point.q()) == 8 for point in nonprimitive)
    assert all(int(point.content()) in (1, 2) for point in nonprimitive)
    isotropic = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))
    assert isotropic and all(int(point.q()) == 0 and int(point.content()) == 1 for point in isotropic)
    opposite = lattice([[0, -1, 0], [-1, 0, 0], [0, 0, 2]])
    opposite_representatives = tuple(RECURSION.exact_vector_orbit_iterator(opposite, -2))
    assert opposite_representatives
    assert all(int(point.q()) == -2 for point in opposite_representatives)
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(opposite, 0, primitive_only=True))) == len(isotropic)


def test_higher_witt_fixed_norm_iterator_uses_complete_eichler_cover():
    from itertools import islice
    from sage_indefinite_port.indefinite.marked_lattices import RECURSION

    plane = lattice([[0, 1], [1, 0]])
    ambient = plane + plane
    representatives = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2))
    assert representatives and all(int(point.q()) == 2 for point in representatives)
    for i, first in enumerate(representatives):
        for second in representatives[i + 1:]:
            assert RECURSION.pointed_isometry(first, second) is None
    assert len(tuple(RECURSION.exact_vector_orbit_iterator(ambient, 2, primitive_only=True))) == len(representatives)
    nonprimitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 8))
    assert nonprimitive and all(int(point.q()) == 8 for point in nonprimitive)
    assert all(int(point.content()) in (1, 2) for point in nonprimitive)
    null_primitive = tuple(RECURSION.exact_vector_orbit_iterator(ambient, 0, primitive_only=True))
    assert null_primitive and all(int(point.q()) == 0 and int(point.content()) == 1 for point in null_primitive)
    nulls = tuple(islice(RECURSION.exact_vector_orbit_iterator(ambient, 0), 5))
    assert [int(point.content()) for point in nulls] == [0, 1, 2, 3, 4]


def test_exact_integral_orthogonal_stream_binary_radical_parametrization():
    from itertools import islice
    from sage_indefinite_port.indefinite.recursive import iter_exact_integral_orthogonal_isometries, _column_matrix

    for gram in ([[2, 0], [0, 0]], [[0, 0], [0, -6]]):
        ambient = lattice(gram)
        elements = tuple(islice(iter_exact_integral_orthogonal_isometries(ambient), 20))
        matrices = [tuple(_column_matrix(element).list()) for element in elements]
        assert len(set(matrices)) == 20
        assert all(element.domain() is ambient and element.codomain() is ambient for element in elements)
        assert all(_column_matrix(element).det() in (-1, 1) for element in elements)


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
