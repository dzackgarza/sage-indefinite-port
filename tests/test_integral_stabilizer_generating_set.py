"""Integral stabilizer generators retain actions invisible modulo two."""
import pytest

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
    expected = integral_basis[1] + ambient.scalar_multiple(ambient.base_ring()(2), integral_basis[2])
    assert any(
        action(integral_basis[1]) == expected
        for generator in reduced
        for action in (generator, ~generator)
    )


@pytest.mark.parametrize("powers", ((2, 3, 5), (2, 3, 5, -2, -3, -5)))
def test_unipotent_powers_generate_the_primitive_integral_action(powers):
    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    e, f, x, y = (rational.module_generator(label) for label in labels)
    ring = rational.base_ring()

    def shear(coefficient):
        return rational.Isom(rational)(dict(zip(
            labels,
            (e, f + rational.scalar_multiple(ring(coefficient), x),
             x, y - rational.scalar_multiple(ring(coefficient), e)),
            strict=True,
        )))

    reduced = _integral_generator_set(ambient, tuple(shear(power) for power in powers))
    assert len(reduced) == 1
    integral_labels = tuple(ambient.module_generating_set())
    integral_basis = tuple(ambient.module_generator(label) for label in integral_labels)
    image = reduced[0](integral_basis[1])
    increment = ambient.scalar_multiple(ambient.base_ring().one(), integral_basis[2])
    assert image in (integral_basis[1] + increment, integral_basis[1] - increment)


def test_equal_norm_block_permutations_reduce_to_generators_of_symmetric_group():
    from sage.all import MatrixGroup, ZZ, matrix

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    basis = tuple(rational.module_generator(label) for label in labels)

    def block_permutation(permutation):
        return rational.Isom(rational)({
            label: basis[2 * permutation[index // 2] + index % 2]
            for index, label in enumerate(labels)
        })

    first = block_permutation((1, 0, 2))
    second = block_permutation((0, 2, 1))
    reduced = _integral_generator_set(ambient, (first, second, first * second))
    integral_labels = tuple(ambient.module_generating_set())
    integral_basis = tuple(ambient.module_generator(label) for label in integral_labels)
    actions = [matrix(ZZ, [
        [int(generator(point).to_vector()(label)) for label in integral_labels]
        for point in integral_basis
    ]) for generator in reduced]
    assert len(reduced) == 2
    assert MatrixGroup(actions).order() == 6
    assert all(sum(abs(entry) for entry in action.list()) == 6 for action in actions)


@pytest.mark.parametrize("powers", ((2, 3), (2, 3, -2, -3)))
def test_large_unipotent_powers_generate_their_full_integral_cyclic_subgroup(powers):
    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    e, f, x, y = tuple(rational.module_generators())
    ring = rational.base_ring()
    step = 2**80

    def shear(coefficient):
        return rational.Isom(rational)(dict(zip(
            labels,
            (e, f + rational.scalar_multiple(ring(coefficient), x),
             x, y - rational.scalar_multiple(ring(coefficient), e)),
            strict=True,
        )))

    reduced = _integral_generator_set(ambient, tuple(shear(power * step) for power in powers))
    e, f, x, y = tuple(ambient.module_generators())
    increment_x = ambient.scalar_multiple(ambient.base_ring()(step), x)
    increment_e = ambient.scalar_multiple(ambient.base_ring()(step), e)
    actions = tuple(tuple(generator(point) for point in (e, f, x, y)) for generator in reduced)

    assert len(actions) == 1
    assert actions[0] in (
        (e, f + increment_x, x, y - increment_e),
        (e, f - increment_x, x, y + increment_e),
    )


def test_large_noncommuting_actions_retain_their_exact_integral_generating_family():
    ambient = lattice([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    rational = ambient.base_change(ambient.base_ring().fraction_field_map())
    labels = tuple(rational.module_generating_set())
    e, f, x, y = tuple(rational.module_generators())
    ring = rational.base_ring()
    step = 2**40

    def multiple(coefficient, point):
        return rational.scalar_multiple(ring(coefficient), point)

    first = rational.Isom(rational)(dict(zip(
        labels, (e, f - multiple(step, y), x + multiple(step, e), y), strict=True,
    )))
    second = rational.Isom(rational)(dict(zip(
        labels, (e + multiple(step, x), f, x, y - multiple(step, f)), strict=True,
    )))
    assert (first * second)(e) == multiple(1 + step**2, e) + multiple(step, x)

    reduced = _integral_generator_set(ambient, (first, second))
    integral_basis = tuple(ambient.module_generators())
    e, f, x, y = integral_basis
    ke, kf, kx, ky = (
        ambient.scalar_multiple(ambient.base_ring()(step), point)
        for point in integral_basis
    )
    first_images = (e, f - ky, x + ke, y)
    second_images = (e + kx, f, x, y - kf)
    inverse_first_images = (e, f + ky, x - ke, y)
    inverse_second_images = (e - kx, f, x, y + kf)
    for expected in (first_images, second_images):
        assert any(
            tuple(action(point) for point in integral_basis) == expected
            for generator in reduced for action in (generator, ~generator)
        )
    assert all(
        tuple(generator(point) for point in integral_basis) in (
            first_images, second_images, inverse_first_images, inverse_second_images,
        )
        for generator in reduced
    )
