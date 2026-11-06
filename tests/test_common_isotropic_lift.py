"""A common rational framing controls integrality of perpendicular lifts."""
from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
from sage_indefinite_port.indefinite.edgewalk_scalars import _rational
from tests.acceptance.consumer import lattice


def test_common_frame_lift_preserves_its_integral_structure():
    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(4)] for i in range(4)])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    points = tuple(rational.module_generator(label) for label in labels)
    field = rational.base_ring()
    scales = (field(2), field.one() / field(2), field.one(), field.one())
    frame = rational.Isom(rational)({label: rational.scalar_multiple(scale, point) for label, scale, point in zip(labels, scales, points)})
    basis = (frame(points[0]), frame(points[2]))
    lifted = IsotropicInduction.lift_invariant_frame(frame, basis, (basis[1], basis[0]))
    assert all(lifted(point) == image for point, image in zip(basis, (basis[1], basis[0])))
    conjugate = ~frame * lifted * frame
    for point in points:
        assert all(_rational(conjugate(point).to_vector()(label)).denominator() == 1 for label in labels)
        assert all(_rational((~conjugate)(point).to_vector()(label)).denominator() == 1 for label in labels)
