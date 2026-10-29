"""Coordinate invariants for Allcock's anisotropic binary cycle."""

from sage_indefinite_port.indefinite.edgewalk_rank2 import _OrientedBinaryCycle
from tests.acceptance.consumer import lattice


def test_binary_cycle_form_uses_the_coordinate_label_order() -> None:
    """The form, coordinate extraction, and reconstruction share one ordered basis."""
    lattice_ = lattice([[1, 1], [1, -1]])
    cycle = _OrientedBinaryCycle(lattice_)
    first, second = tuple(lattice_.module_generator(label) for label in cycle.labels)

    assert cycle.form == (int(first.q()), int(lattice_.b(first, second)), int(second.q()))
    for coordinates in ((1, 0), (0, 1), (2, -3), (-5, 2)):
        point = cycle.point(coordinates)
        assert cycle.coordinates(point) == coordinates
        assert cycle.pairing(coordinates, coordinates) == int(point.q())

    result = cycle.cycle(1, cycle.point((1, 0)))
    assert result is not None
    period, representatives = result
    assert representatives
    assert all(0 < int(point.q()) <= 1 and int(point.content()) == 1 for point in representatives)
    basis = tuple(lattice_.module_generator(label) for label in cycle.labels)
    assert period.determinant() == 1
    assert all(lattice_.b(period(left), period(right)) == lattice_.b(left, right) for left in basis for right in basis)
