"""Exact low-rational-Witt-index isometry witnesses through public lattice Mor."""
from tests.acceptance.consumer import lattice


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
