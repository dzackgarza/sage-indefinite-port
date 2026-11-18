"""The isotropic flag invariant at the public flag-equivalence entry (#18).

Source: CombinedAlgorithms.h:1361-1395 (the invariant), :1401-1404 (Equivalence_IsotropicKstuff
returns no equivalence when the invariants differ), :2429-2433 (INDEF_FORM_Invariant_IsotropicKflag).

L = 2U + <-2> + <-6> with basis e1, f1, e2, f2, v, w and l = 2e2 + 2f2 + v + w
(references/dawes-L2-boundary.md). The flags e1 < E1 = <e1, e2> and e1 < E2 = <e1, l>
are not equivalent: E1^perp/E1 = <-2> + <-6> has determinant 12, and E2^perp/E2 has
determinant 3, because every x.l is even while x.e1 = 1 for x = f1. The invariant records
E^perp/E, so it separates them. An isometry g of L carries e1 < E1 to an equivalent flag,
on which the invariant is equal.
"""

from tests.acceptance.consumer import are_isomorphic, basis, flag_category, lattice
from sage_indefinite_port.indefinite.marked_lattices import INDUCTION

GRAM = [
    [0, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, -2, 0],
    [0, 0, 0, 0, 0, -6],
]


def _flag(lattice_, line, other):
    plane = lattice_.primitive_isotropic_subobject(line, other)
    inclusion = plane.inclusion()
    line_in_plane = inclusion.domain().primitive_isotropic_subobject(inclusion.lift(line))
    return flag_category()((line_in_plane.inclusion(), inclusion))


def _invariant(lattice_, line, other):
    return INDUCTION.flag_invariant(lattice_, (line, other), (1, 2))


def test_invariant_is_equal_on_equivalent_flags():
    lattice_ = lattice(GRAM)
    e1, f1, e2, f2, v, w = basis(lattice_)
    g = lattice_.eichler_transvection(f1, v)
    moved_line, moved_other = g(e1), g(e2)
    assert moved_line.q() == 0 and moved_other.q() == 0 and moved_line.b(moved_other) == 0
    assert moved_line != e1

    assert _invariant(lattice_, e1, e2) == _invariant(lattice_, moved_line, moved_other)
    assert are_isomorphic(flag_category(), _flag(lattice_, e1, e2), _flag(lattice_, moved_line, moved_other))


def test_invariant_separates_inequivalent_flags():
    lattice_ = lattice(GRAM)
    e1, f1, e2, f2, v, w = basis(lattice_)
    null = e2 + e2 + f2 + f2 + v + w
    assert null.q() == 0 and null.b(e1) == 0

    assert _invariant(lattice_, e1, e2) != _invariant(lattice_, e1, null)
    assert not are_isomorphic(flag_category(), _flag(lattice_, e1, e2), _flag(lattice_, e1, null))
