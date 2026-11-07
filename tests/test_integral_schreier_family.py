"""Factor realization retains every integral stabilizer action and inverse."""

import pytest

from tests.acceptance.consumer import lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion


@pytest.mark.parametrize(
    "gram",
    (
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        [[0, 4, 0], [4, 0, 0], [0, 0, -8]],
    ),
    ids=("two_u", "scaled_hyperbolic_line"),
)
def test_factored_integral_family_equals_every_source_action(gram):
    from sage_indefinite_port.indefinite.recursive import _integral_stabilizer_family

    ambient = lattice(gram)
    integral_basis = tuple(ambient.module_generators())
    line = integral_basis[0]
    action = EichlerRecursion().isotropic_stabilizer_action((line,), ((-line,),))
    original = tuple(action.finite_representation().lattice_stabilizer().group_generators())
    realized = tuple(_integral_stabilizer_family(action))
    rational = action.lattice_inclusion().codomain().module_over_extension()
    rational_basis = tuple(rational.module_generators())
    inclusion = action.lattice_inclusion()

    assert len(realized) == len(original)
    for source, target in zip(original, realized, strict=True):
        assert all(source(point) == target(point) for point in rational_basis)
        assert all((~source)(point) == (~target)(point) for point in rational_basis)
        assert all(
            inclusion.is_in_image(inclusion.codomain().wrap(mapping(inclusion(point).underlying_element())))
            for mapping in (target, ~target)
            for point in integral_basis
        )
