"""Independently sourced nonisometry oracle, including genus equality.

Pinned source: references/vendor/Hecke.jl@e2ab5716/test/QuadForm/Quad/
ZLatticeAutIso.jl, lines 732--736 ("Example from Conway Sloane Chapter 15
p.393"). The two displayed positive-definite ternary Gram matrices have
equal genera but are not integrally isometric. The negative verdict is also
proved below by the distinct numbers of vectors of squared norm two.

The *indefinite* Conway--Sloane same-genus counterexample is already tested
in tests/acceptance/test_isometry.sage; it is intentionally not duplicated.
"""

from tests.acceptance.consumer import MARKED_ISOMETRY, lattice, require


def test_hecke_conway_sloane_ternary_same_genus_but_nonisometric():
    r"""Hecke.jl@e2ab5716, ZLatticeAutIso.jl:732--736.

    For the first lattice, Q=2(x²+xy+y²)+18z². Q=2 forces
    z=0 and x²+xy+y²=1, with six solutions
    (±1,0), (0,±1), (1,-1), (-1,1).

    For the second, Q=6(x²+xy+y²)+2z². The binary form is
    nonnegative integral and vanishes only at (0,0), so Q=2
    forces x=y=0, z=±1: exactly two vectors. An integral
    isometry bijects norm-two vectors, proving nonisometry
    independently of the port implementation. Equality of the
    genera is the separately pinned Hecke source assertion.
    """
    require(MARKED_ISOMETRY)
    first = lattice([[2, 1, 0], [1, 2, 0], [0, 0, 18]])
    second = lattice([[6, 3, 0], [3, 6, 0], [0, 0, 2]])

    assert first.genus() == second.genus()
    assert not first.is_isometric(second)
    assert first.Isom(second).is_empty()
