"""Exact full-projection degenerate marked stabilizers with nonprimitive quotient image."""
import pytest
from sage.matrix.constructor import identity_matrix, matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.modules.free_module_element import vector
from sage_extensions.all import Lattices, ZZ
from sage_indefinite_port.indefinite.degenerate_markings import degenerate_marking_generators
from sage_indefinite_port.indefinite.marked_lattices import ARROWS


def _finite_closure(generators, expected):
    identity = identity_matrix(SageZZ, 3)
    key = lambda g: tuple(map(int, g.list()))
    seen = {key(identity): identity}
    queue = [identity]
    while queue:
        old = queue.pop()
        for generator in generators:
            moved = generator * old
            if key(moved) not in seen:
                seen[key(moved)] = moved
                queue.append(moved)
            assert len(seen) <= expected, "the proposed stabilizer has unexpected additional elements"
    assert len(seen) == expected
    return tuple(seen.values())


@pytest.mark.parametrize("m, based_order, moving_order", [(2, 2, 4), (5, 1, 2)])
def test_nonprimitive_full_projection_has_complete_glue_stabilizer(m, based_order, moving_order):
    # L=U+R and S=<e,mf+r> are primitive, but S -> L/R has index m.
    # O(U) is {+/-id, +/-swap}. For m>1, only +/-id preserves
    # <e,mf>. Writing g(r)=d r and g(f)=t f+b r, integrality forces
    # m*b+d=t. Pointwise S: t=1; moving S: t=+-1.
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    s = m*f + r
    for fixed, expected in ((True, based_order), (False, moving_order)):
        matrices = tuple(degenerate_marking_generators(ambient, (e, s), fixed))
        group = _finite_closure(matrices, expected)
        lattice_rows = matrix(SageZZ, [[1, 0, 0], [0, m, 1]]).row_module()
        for g in group:
            assert abs(g.det()) == 1
            gram = matrix(SageZZ, [[int(ambient.b(x, y)) for y in (e, f, r)] for x in (e, f, r)])
            assert g.transpose() * gram * g == gram
            images = tuple(g * vector(SageZZ, row) for row in ((1, 0, 0), (0, m, 1)))
            if fixed:
                assert images == (vector(SageZZ, [1, 0, 0]), vector(SageZZ, [0, m, 1]))
            else:
                assert all(image in lattice_rows for image in images)


def test_public_primitive_full_projection_stabilizer_does_not_require_primitive_quotient():
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    inclusion = ambient.subobject_on((e, 5*f + r)).inclusion()
    assert inclusion.is_primitive()
    assert len(tuple(ambient.O().pointwise_stabilizer(inclusion).group_generators())) == 0 or all(
        g(e) == e and g(5*f+r) == 5*f+r
        for g in ambient.O().pointwise_stabilizer(inclusion).group_generators()
    )
    assert all(
        inclusion.is_in_image(g(e)) and inclusion.is_in_image(g(5*f+r))
        for g in ambient.O().setwise_stabilizer(inclusion).group_generators()
    )


@pytest.mark.parametrize("pointwise,expected_image", [(True, 10), (False, 20)])
def test_rank_two_radical_full_projection_integral_gl_image(pointwise, expected_image):
    # For S=<e,5f+r1> in U+0^2, the radical action D in GL_2(Z)
    # must send r1 to +r1 mod 5 when fixing S pointwise, or to
    # +/-r1 mod 5 setwise. As det(D)=+/-1, the respective images
    # in GL_2(F5) have 5*2=10 and 2*5*2=20 elements.
    from sage.rings.finite_rings.finite_field_constructor import GF

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, r1, r2 = tuple(ambient.module_generators())
    marked = (e, 5*f+r1)
    basis = tuple(ambient.module_generators())
    grams = matrix(SageZZ, [[int(ambient.b(x, y)) for y in basis] for x in basis])
    generators = tuple(degenerate_marking_generators(ambient, marked, pointwise))
    for action in generators:
        assert abs(action.det()) == 1 and action.transpose()*grams*action == grams
        moved = tuple(action*vector(SageZZ, [1,0,0,0]) if i == 0 else
                      action*vector(SageZZ, [0,5,1,0]) for i in range(2))
        if pointwise:
            assert moved == (vector(SageZZ, [1,0,0,0]), vector(SageZZ, [0,5,1,0]))
        else:
            span = matrix(SageZZ, [[1,0,0,0],[0,5,1,0]]).row_module()
            assert all(v in span for v in moved)

    field = GF(5)
    residue = tuple(matrix(field, [[g[i,j] for j in (2,3)] for i in (2,3)])
                    for g in generators)
    one = identity_matrix(field, 2)
    key = lambda m: tuple(int(c) for c in m.list())
    seen = {key(one)}
    frontier = [one]
    while frontier:
        element = frontier.pop()
        for g in residue:
            moved = g*element
            k = key(moved)
            if k not in seen:
                seen.add(k)
                frontier.append(moved)
            assert len(seen) <= expected_image
    assert len(seen) == expected_image
