"""The hyperbolic-plane canonicalization and transporter regression, deliberately kept separate and slow.

This is the case that was removed from the rapid file only for runtime: on `U = [[0,1],[1,0]]` with
`q(x, y) = 2xy`, transporting `(2, 1)` onto its orbit's minimum `(-2, -1)` does **not** take the already-integral
branch of `pointed_isometry`; it reaches `_lattice_action(target).transporter(image)`, the correction locus recorded
as correct-but-costly. Measured at 57 s for the transport alone.

An "already canonical, so the identity path" test does **not** cover this obligation: the identity path returns
before any of the complement, extension or correction machinery runs, so a real transport between two distinct
members of the orbit is the thing that has to be exercised. Run this file when `pointed_isometry`, the correction
locus or the canonicalizer is touched -- not on every iteration.
"""
import pytest

from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import element, lattice

_RECURSION = EichlerRecursion()
_LATTICE = None


def _hyperbolic():
    global _LATTICE
    if _LATTICE is None:
        _LATTICE = lattice([[0, 1], [1, 0]])
    return _LATTICE


def _coordinates(lattice_, vector) -> tuple[int, ...]:
    values = vector.to_vector()
    return tuple(int(values(label)) for label in lattice_.module_generating_set())


@pytest.mark.slow
def test_hyperbolic_canonicalization_transports_a_distinct_member() -> None:
    r"""`can(2, 1) = (-2, -1)` on `U`, and the transporter is a real isometry carrying one onto the other.

    Square 4 of `U` is the single orbit `{(+-1,+-2), (+-2,+-1)}`, whose least member in the declared
    height-then-lex order is `(-2, -1)`. Reaching it from `(2, 1)` is a genuine transport, so this asserts the
    correction branch produces an isometry that satisfies `transporter(representative) == input` and preserves the
    pairing -- not merely that some vector was returned.
    """
    L = _hyperbolic()
    input_vector = element(L, (2, 1))

    canonical = _RECURSION.marked_orbit_canonical(input_vector)

    assert _coordinates(L, canonical.representative) == (-2, -1)
    assert canonical.height == 2
    assert canonical.transporter(canonical.representative) == input_vector, (
        "the transporter must carry the representative onto the input"
    )
    # it is an isometry, not a coordinate relabelling: norms and pairings are preserved
    for probe in ((-2, -1), (1, 2), (3, -5), (0, 1)):
        other = element(L, probe)
        assert int(canonical.transporter(other).q()) == int(other.q())
        assert int(L.b(canonical.transporter(other), canonical.transporter(input_vector))) == int(
            L.b(other, input_vector)
        )


@pytest.mark.slow
def test_hyperbolic_enumeration_emits_the_derived_minima() -> None:
    r"""The two square slices whose minima are derivable by hand: square 4 gives `(-2,-1)`, square 2 gives `(-1,-1)`."""
    L = _hyperbolic()

    square_four = list(_RECURSION.marked_orbit_representatives(L, 4, limit=1))
    assert [_coordinates(L, vector) for vector in square_four] == [(-2, -1)]

    square_two = list(_RECURSION.marked_orbit_representatives(L, 2, limit=1))
    assert [_coordinates(L, vector) for vector in square_two] == [(-1, -1)]
