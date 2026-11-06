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


def test_invariant_module_frame_retains_the_module_and_pairing():
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[0, 1], [1, 0]])
    points = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    action = EichlerRecursion().isotropic_stabilizer_action((points[0],), ((-points[0],),))
    inclusion = action.invariant_module()
    frame = IsotropicInduction.invariant_module_frame(ambient, inclusion)
    space = frame.domain()
    module = inclusion.domain()
    module_points = tuple(module.module_generator(label) for label in module.module_generating_set())
    framed_points = tuple(space.module_generator(label) for label in space.module_generating_set())
    assert tuple(frame(point) for point in framed_points) == tuple(inclusion(point).underlying_element() for point in module_points)
    for left in framed_points:
        for right in framed_points:
            assert space.b(left, right) == frame.codomain().b(frame(left), frame(right))
