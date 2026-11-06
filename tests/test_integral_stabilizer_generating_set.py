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
    expected = integral_basis[1] + ambient.scalar_multiple(ambient.base_ring()(2), integral_basis[2])
    assert any(
        action(integral_basis[1]) == expected
        for generator in reduced
        for action in (generator, ~generator)
    )


def test_unipotent_powers_generate_the_primitive_integral_action():
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

    reduced = _integral_generator_set(ambient, (shear(2), shear(3), shear(5)))
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
