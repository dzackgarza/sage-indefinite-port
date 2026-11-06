"""Integral stabilizer generators retain actions invisible modulo two."""
from sage_indefinite_port.indefinite.recursive import _integral_generator_set
from tests.acceptance.consumer import lattice


def test_generating_set_keeps_sign_and_unipotent_kernel_actions():
    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    e, f, x, y = (rational.module_generator(label) for label in labels)
    ring = rational.base_ring()
    identity = rational.Isom(rational).identity()
    sign = rational.Isom(rational)({label: -rational.module_generator(label) for label in labels})
    shear = rational.Isom(rational)(dict(zip(labels, (e, f + rational.scalar_multiple(ring(2), x), x, y - rational.scalar_multiple(ring(2), e)), strict=True)))
    reduced = _integral_generator_set(ambient, (identity, sign, shear, sign, shear, identity))
    assert len(reduced) == 2
    integral_basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    assert any(all(generator(point) == -point for point in integral_basis) for generator in reduced)
    assert any(generator(integral_basis[1]) == integral_basis[1] + ambient.scalar_multiple(ambient.base_ring()(2), integral_basis[2]) for generator in reduced)
