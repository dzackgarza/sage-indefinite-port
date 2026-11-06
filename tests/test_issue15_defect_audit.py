"""#15: the migration-inventory defect rows, checked against the current code.

`references/migration-inventory.md` records six correctness defects in the Lorentzian
cell and edgewalk code with an inference- or read-basis citation, and #15 requires each
to be *checked* rather than assumed resolved because the code moved. Each test below is
the reproducer for one row: it is written against the behaviour the defect would produce,
so it fails if that behaviour returns.

Rows not represented here are recorded in the issue rather than guessed at:

- the traversal signature (`perfect_domain_traversal` taking Gram rows) is checked by the
  frozen acceptance module, which passes a lattice and reads the owned result;
- the rank-2 `U` special case no longer keys on the basis: the branch that depended on both
  basis vectors being isotropic is gone, and `U` itself is covered by
  `tests/acceptance/test_reflective.sage`'s rank-two cases.
"""

import pytest

from sage.all import QQ, ZZ as SageZZ
from sage_extensions.all import ZZ, Lattices

from sage_indefinite_port.indefinite import edgewalk
from sage_indefinite_port.indefinite import lorentzian_cells


def test_isotropic_mode_on_an_anisotropic_lattice_refuses_instead_of_searching_forever() -> None:
    """Row 1: `_search_initial_vectors` doubled a shell bound with no isotropic vector to find.

    The pinned isotropic search has no pair of isotropic lines to locate, so the correct
    answer is a typed refusal naming the lattice, not an unbounded increase of the bound.
    A regression here does not fail an assertion: it does not return at all, which is why
    the reproducer runs the search directly and requires the refusal.
    """
    anisotropic = Lattices(ZZ)([[2, 0], [0, -6]])
    assert anisotropic.is_isotropic() is False

    with pytest.raises(ValueError, match="anisotropic"):
        lorentzian_cells._search_initial_vectors(anisotropic, (0, 1), lorentzian_cells.ZERO_SQUARES)


def test_a_rational_coordinate_is_refused_rather_than_truncated() -> None:
    """Row 6: `int(entry)` silently truncated a rational coordinate to its floor.

    Every edgewalk scalar conversion goes through `_integer`; a lattice vector whose
    coordinate is not an integer is a defect in the caller, and the floor of that
    coordinate is not a value the algorithm may proceed with.
    """
    with pytest.raises(ArithmeticError, match="not integral"):
        edgewalk._integer(QQ(3) / QQ(2))
    with pytest.raises(ArithmeticError, match="not integral"):
        edgewalk._integer(QQ(-1) / QQ(3))

    # The integral path is unchanged, including a negative coordinate and a Sage integer.
    assert edgewalk._integer(SageZZ(5)) == 5
    assert edgewalk._integer(SageZZ(-7)) == -7
    assert edgewalk._integer(QQ(4) / QQ(2)) == 2


def test_the_rank_two_branch_reports_every_wall_it_found() -> None:
    """Row 5: the rank-2 branch reported reflectivity with one wall for signature (1,1).

    The anisotropic binary cycle terminates with the two walls bounding the chamber, so a
    reflective rank-two lattice reports both roots and both vertices. Reporting a single
    wall would claim a fundamental domain bounded on one side only.
    """
    for gram in ([[1, 0], [0, -2]], [[2, 1], [1, -2]]):
        lattice_ = Lattices(ZZ)(gram)
        report = edgewalk.edgewalk_fundamental_domain(lattice_)
        assert report.is_reflective() is True, gram
        assert len(tuple(report.simple_roots())) == 2, gram
        assert len(tuple(report.vertices())) == 2, gram


def test_a_split_rank_two_lattice_is_not_reported_reflective_with_a_wall() -> None:
    """Row 5's other half: a split (1,1) lattice whose reflection group is infinite.

    `U` has infinitely many walls, so the edgewalk must terminate with a non-reflective
    verdict and no fundamental polygon rather than with whatever partial chamber it
    reached when its stopping rule fired.
    """
    hyperbolic_plane = Lattices(ZZ)([[0, 1], [1, 0]])
    report = edgewalk.edgewalk_fundamental_domain(hyperbolic_plane)

    assert report.is_reflective() is False
    assert len(tuple(report.simple_roots())) == 0
    assert len(tuple(report.vertices())) == 0
