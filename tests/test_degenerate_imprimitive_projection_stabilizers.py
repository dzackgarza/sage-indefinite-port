"""Complete rank-one degenerate-ambient stabilizers with imprimitive quotient projection."""
import pytest

from tests.acceptance.consumer import lattice
from sage_extensions.all import Lattices, ZZ
from sage_indefinite_port.indefinite.marked_lattices import ARROWS, automorphism_group_generators


@pytest.mark.parametrize('multiplier,pointwise_order,moving_order', [
    (1, 2, 4), (2, 2, 4), (5, 1, 2),
])
def test_primitive_mixed_radical_marking_uses_full_integral_congruence(multiplier, pointwise_order, moving_order):
    # O(<2> + <0>) consists of e -> s e + b r, r -> d r,
    # where s,d are independent signs and b an arbitrary integer.
    # For v=m e+r, the based stabilizer has s=1 and m b+d=1;
    # the moving stabilizer has s=t and m b+d=t for t=+-1.
    # Since d=+-1, these conditions give exactly the stated orders.
    ambient = lattice([[2, 0], [0, 0]])
    e, r = tuple(ambient.module_generators())
    v = multiplier * e + r
    source = Lattices(ZZ)([[2 * multiplier * multiplier]])
    coslice = Lattices(ZZ).CosliceUnder(source)
    based = coslice(source.hom([v], codomain=ambient))
    moving = ARROWS(ambient.subobject_on((v,)).inclusion())
    assert int(v.content()) == 1
    assert int(v.to_vector()(tuple(ambient.module_generating_set())[0])) == multiplier

    def matrix_key(g):
        labels = tuple(ambient.module_generating_set())
        return tuple(tuple(int(g(p).to_vector()(label)) for label in labels) for p in (e, r))

    for category, marked, expected, signs in (
        (coslice, based, pointwise_order, (1,)),
        (ARROWS, moving, moving_order, (1, -1)),
    ):
        generators = tuple(automorphism_group_generators(category, marked))
        assert all(g(v) in tuple(sign * v for sign in signs) for g in generators)
        seen = {matrix_key(ambient.O().identity()): ambient.O().identity()}
        frontier = [ambient.O().identity()]
        while frontier:
            current = frontier.pop()
            for generator in generators:
                moved = generator * current
                key = matrix_key(moved)
                if key not in seen:
                    seen[key] = moved
                    frontier.append(moved)
                assert len(seen) <= expected, 'invalid generator enlarges the exact stabilizer'
        assert len(seen) == expected


def test_isotropic_imprimitive_projection_preserves_entire_shear_kernel():
    # L=U + <0>, v=5e+r. Every pointwise stabilizer fixes e and r,
    # and sends f to f+t*r. Conversely all such t are allowed, so its
    # shear parameter ideal must be the entire ZZ, not a proper subgroup.
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    v = 5 * e + r
    source = Lattices(ZZ)([[0]])
    coslice = Lattices(ZZ).CosliceUnder(source)
    mark = coslice(source.hom([v], codomain=ambient))
    generators = tuple(automorphism_group_generators(coslice, mark))
    radical_label = tuple(ambient.module_generating_set())[-1]
    params = []
    for generator in generators:
        assert generator(v) == v
        assert generator(e) == e
        assert generator(r) == r
        parameter = generator(f).to_vector()(radical_label)
        assert generator(f) == f + parameter * r
        params.append(parameter)
    assert ambient.base_ring().ideal(params) == ambient.base_ring().ideal(1)
