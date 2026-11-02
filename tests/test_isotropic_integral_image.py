"""Regression: `_integral_image` must read rational-space coordinates without `int()`.

The boundary-building acceptance tests (`test_boundaries.sage`, Dawes family) failed with

    TypeError: int() argument must be a string, a bytes-like object or a real number,
    not 'OwnedJoinCategory.element_class'

at `isotropic_induction.py::IsotropicInduction._integral_image`. The value being converted is `moved`'s
coordinate, and `moved` lives in the *rational* space, so its coordinates are base-ring elements rather
than Python integers. `int()` is the wrong conversion for them; the port's exact-rational helper reads them
and keeps the integrality invariant instead of truncating.

These two cases pin the contract at the owning API and need no lattice bigger than the hyperbolic plane.
"""
import pytest

from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
from tests.acceptance.consumer import lattice


def _hyperbolic_plane():
    lattice_ = lattice([[0, 1], [1, 0]])
    anchor = next(iter(lattice_.module_generators()))
    return lattice_, anchor


def test_integral_image_reads_rational_space_coordinates() -> None:
    """An identity of the rational space carries an integral vector to itself, not to a TypeError."""
    lattice_, vector = _hyperbolic_plane()
    image = IsotropicInduction._integral_image(lattice_, lambda point: point, vector)

    assert image == vector
    for label in lattice_.module_generating_set():
        coordinate = image.to_vector()(label)
        # The container type is a base-ring element; the point is that it converts exactly.
        assert int(coordinate) == (1 if label == lattice_.module_generating_set()[0] else 0)


def test_integral_image_rejects_a_non_integral_image_as_an_arithmetic_defect() -> None:
    """A genuinely non-integral image is an arithmetic defect, reported as such rather than truncated."""
    lattice_, vector = _hyperbolic_plane()
    rational = lattice_.base_change(lattice_.base_ring().fraction_field_map())
    half = rational.base_ring().one() / rational.base_ring()(2)

    with pytest.raises(ArithmeticError):
        IsotropicInduction._integral_image(lattice_, lambda point: rational.scalar_multiple(half, point), vector)
