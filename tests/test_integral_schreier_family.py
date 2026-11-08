"""Factor realization retains every integral stabilizer action and inverse."""

import pytest

from tests.acceptance.consumer import lattice


@pytest.mark.parametrize(
    "gram",
    (
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        [[0, 4, 0], [4, 0, 0], [0, 0, -8]],
    ),
    ids=("two_u", "scaled_hyperbolic_line"),
)
def test_factored_integral_family_equals_every_source_action(gram):
    from sage_extensions.all import IntegralStructureAction
    from sage_extensions.categories.modules.pure.modules import ModulesOverIntegralDomains
    from sage_indefinite_port.indefinite.recursive import _integral_stabilizer_family

    ambient = lattice(gram)
    integral_basis = tuple(ambient.module_generators())
    inclusion = ModulesOverIntegralDomains(ambient.base_ring())(ambient).generic_fibre_map()
    space = inclusion.codomain()
    rational = space.module_over_extension()
    rational_basis = tuple(rational.module_generators())
    if len(rational_basis) == 4:
        e, f, u, v = rational_basis
        images = (
            (e, f, u, v), (e, f, -u, -v),
            (e, f, v, u), (e, f, -v, -u),
            (-e, -f, u, v), (e, f - v, u + e, v), (e, f - u, u, v + e),
        )
        bound = inclusion
    else:
        e, f, w = rational_basis
        field = rational.base_ring()
        extra = rational.scalar_multiple(field(1) / field(4), e) + rational.scalar_multiple(field(1) / field(2), w)
        images = ((e, f, -w), (-e, -f, w), (e, f + extra, e + w))
        bound = space.subobject_on((*tuple(inclusion(point) for point in integral_basis), space.wrap(extra))).inclusion()
    group_generators = tuple(rational.O()(family) for family in images)
    assert all(
        rational.b(generator(left), generator(right)) == rational.b(left, right)
        for generator in group_generators for left in rational_basis for right in rational_basis
    )
    action = IntegralStructureAction(rational.O().subgroup(group_generators), inclusion).with_stable_containing_bound(bound)
    original = tuple(action.finite_representation().lattice_stabilizer().group_generators())
    realized = tuple(_integral_stabilizer_family(action))
    from sage_indefinite_port.indefinite.recursive import _integral_stabilizer_action_columns

    columns = tuple(_integral_stabilizer_action_columns(ambient, action))
    labels = tuple(ambient.module_generating_set())
    assert len(columns) == len(original)
    for source, action_columns in zip(original, columns, strict=True):
        inverse_columns = action_columns.inverse()
        for index, point in enumerate(integral_basis):
            for mapping, coefficients in ((source, action_columns), (~source, inverse_columns)):
                image = ambient.linear_combination({
                    label: ambient.base_ring()(int(coefficients[row, index]))
                    for row, label in enumerate(labels)
                })
                assert inclusion(image).underlying_element() == mapping(inclusion(point).underlying_element())

    assert len(realized) == len(original)
    for source, target in zip(original, realized, strict=True):
        assert all(source(point) == target(point) for point in rational_basis)
        assert all((~source)(point) == (~target)(point) for point in rational_basis)
        assert all(
            inclusion.is_in_image(inclusion.codomain().wrap(mapping(inclusion(point).underlying_element())))
            for mapping in (target, ~target)
            for point in integral_basis
        )


def test_factored_family_uses_the_noncommuting_lattice_action():
    """S3 permutes three distinct axis indices; central sign remains in the stabilizer."""
    from sage_extensions.all import IntegralStructureAction, Lattices, Modules, QQ, ZZ
    from sage_indefinite_port.indefinite.recursive import _integral_stabilizer_family

    rational = Lattices(QQ)(3)
    e, f, w = tuple(rational.module_generators())
    space = Modules(QQ).restriction_of_scalars(ZZ.fraction_field_map())(rational)
    selected = space.subobject_on((space(e), space(2 * f), space(3 * w)))
    bound = space.subobject_on(tuple(space(point) for point in (e, f, w)))
    group = rational.O().subgroup(tuple(
        rational.O()(images) for images in ((f, w, e), (f, e, w), (-e, -f, -w))
    ))
    action = IntegralStructureAction(group, selected.inclusion()).with_stable_containing_bound(bound.inclusion())
    original = tuple(action.stabilizer().group_generators())
    realized = tuple(_integral_stabilizer_family(action))

    assert action.right_cosets().cardinality() == 6
    assert len(original) == len(realized) == 18
    assert all(
        source(point) == target(point) and (~source)(point) == (~target)(point)
        for source, target in zip(original, realized, strict=True)
        for point in (e, f, w)
    )
    assert any(all(generator(point) == -point for point in (e, f, w)) for generator in realized)

    from sage_indefinite_port.indefinite.recursive import _integral_stabilizer_action_columns

    columns = tuple(_integral_stabilizer_action_columns(selected, action))
    labels = tuple(selected.module_generating_set())
    basis = tuple(selected.module_generators())
    assert len(columns) == len(original)
    for source, coefficients in zip(original, columns, strict=True):
        for mapping, matrix_columns in ((source, coefficients), (~source, coefficients.inverse())):
            for index, point in enumerate(basis):
                image = selected.linear_combination({
                    label: ZZ(int(matrix_columns[row, index]))
                    for row, label in enumerate(labels)
                })
                assert selected.inclusion()(image).underlying_element() == mapping(
                    selected.inclusion()(point).underlying_element()
                )
