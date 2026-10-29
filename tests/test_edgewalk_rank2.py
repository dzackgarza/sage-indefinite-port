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
    basis = tuple(lattice_.module_generators())
    assert period.determinant() == 1
    assert all(lattice_.b(period(left), period(right)) == lattice_.b(left, right) for left in basis for right in basis)


def test_edgewalk_period_action_matrix_uses_coordinate_label_order() -> None:
    """The period action rows align with the coordinate frame of the sublattice."""
    from dzack_research.preamble.all import ZZ, Lattices
    from sage.all import QQ, matrix

    from sage_indefinite_port.indefinite.edgewalk import (
        _isometry_row_matrix_in_coordinate_order,
    )

    lattice_ = Lattices(ZZ)([[2, 1], [1, 2]], module_generators=("right", "left"))
    labels = tuple(lattice_.module_generating_set())
    row_matrix = matrix(ZZ, [[-1, 1], [0, 1]])
    images = {
        label: lattice_.linear_combination(
            {target_label: row_matrix[row, column] for column, target_label in enumerate(labels) if row_matrix[row, column]}
        )
        for row, label in enumerate(labels)
    }
    isometry = lattice_.O()(images)

    assert _isometry_row_matrix_in_coordinate_order(lattice_, isometry) == row_matrix.change_ring(QQ)
