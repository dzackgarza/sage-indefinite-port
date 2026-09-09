"""Integral-structure actions for rational isometry groups.

These are source specimens for TASK-03-1.  They are committed unexecuted while
the research checkout remains before its terminal verification phase.
"""

from dzack_research.preamble.all import FreeModule, Lattices, Modules, QQ, ZZ, module_embedding

from sage_indefinite_port.groups.integral_structures import (
    IntegralStructureAction,
    RationalMatrixGroup,
)


def _standard_hyperbolic_lattice():
    plane = Lattices(QQ)("U")
    restriction = Modules(QQ).restriction_of_scalars(
        ZZ.Mor(QQ)(lambda element: QQ(element))
    )
    space = restriction(plane)
    e, f = plane.module_generators()
    standard = module_embedding(
        FreeModule(ZZ, 2),
        space,
        {0: space.wrap(e), 1: space.wrap(f)},
    )
    return plane, space, standard, e, f


def test_an_integral_group_keeps_the_selected_lattice() -> None:
    plane, _space, standard, e, f = _standard_hyperbolic_lattice()
    swap = plane.Aut()({0: f, 1: e})
    action = IntegralStructureAction(RationalMatrixGroup(plane, (swap,)), standard)

    assert action.invariant_overlattice() is standard
    assert action.quotient_exponent() == 1


def test_a_rational_involution_closes_to_the_smallest_invariant_overlattice() -> None:
    plane, space, standard, e, f = _standard_hyperbolic_lattice()
    involution = plane.Aut()(
        {
            0: plane.scalar_multiple(QQ(2), f),
            1: plane.scalar_multiple(QQ(1) / 2, e),
        }
    )
    assert involution * involution == plane.Aut().one()

    action = IntegralStructureAction(
        RationalMatrixGroup(plane, (involution,)),
        standard,
    )
    invariant = action.invariant_overlattice()
    half_e = space.wrap(plane.scalar_multiple(QQ(1) / 2, e))

    assert invariant.codomain() is space
    assert invariant.is_in_image(space.wrap(e))
    assert invariant.is_in_image(space.wrap(f))
    assert invariant.is_in_image(half_e)
    assert action.quotient_exponent() == 2

    for basis_vector in invariant.domain().module_generators():
        embedded = invariant(basis_vector).underlying_element()
        moved = space.wrap(involution(embedded))
        assert invariant.is_in_image(moved)
