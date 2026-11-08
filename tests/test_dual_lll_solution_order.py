"""CVP witnesses follow the pinned outer-to-inner recursive coordinates."""

from sage.all import QQ, ZZ, matrix


def test_solutions_follow_recursive_coordinates_with_exact_witnesses():
    from sage_indefinite_port.indefinite.lorentzian_cells import _dual_lll_ordered_solutions

    transform = matrix(QQ, [[1, 2, 0], [0, 1, 1], [0, 0, 1]])
    large = ZZ(2) ** 100
    # The recursive columns are (3,2,2), (11,5,0), (8,3,1),
    # (large,0,-1), and (large+1,0,-1). Coordinate zero varies fastest.
    solutions = (
        ((ZZ(3), ZZ(0), ZZ(2)), QQ(7) / 3),
        ((ZZ(1), ZZ(5), ZZ(0)), QQ(11) / 5),
        ((ZZ(4), ZZ(2), ZZ(1)), QQ(13) / 7),
        ((large - 1, ZZ(1), ZZ(-1)), QQ(17) / 11),
        ((large - 2, ZZ(1), ZZ(-1)), QQ(19) / 13),
    )
    expected = (solutions[4], solutions[3], solutions[1], solutions[2], solutions[0])
    assert _dual_lll_ordered_solutions(solutions, transform) == expected
    assert all(actual is original for actual, original in zip(
        _dual_lll_ordered_solutions(solutions, transform), expected, strict=True
    ))


def test_empty_and_rank_zero_solution_families():
    from sage_indefinite_port.indefinite.lorentzian_cells import _dual_lll_ordered_solutions

    assert _dual_lll_ordered_solutions((), matrix(ZZ, 3, 3)) == ()
    solutions = (((), QQ(0)), ((), QQ(1)))
    assert _dual_lll_ordered_solutions(solutions, matrix(ZZ, 0, 0)) == solutions
