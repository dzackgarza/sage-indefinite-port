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


@pytest.mark.parametrize("pointwise,expected_image", [(True, 1), (False, 20)])
def test_full_projection_with_characteristic_radical_intersection(pointwise, expected_image):
    # S=<e,5f+r1,r2> in U+<0>^2 is primitive. Its quotient
    # projection <e,5f> has index five, and W=S cap rad(L)=<r2>.
    # Pointwise S fixes r2. A shear r1 -> r1+a*r2 forces
    # f -> f-(a/5)*r2, so integrality requires 5|a; every a/5
    # is realized. Setwise S additionally permits arbitrary shifts
    # of the marked vector by W, so the radical image is the
    # 2*2*5-element triangular subgroup modulo five.
    from sage.rings.finite_rings.finite_field_constructor import GF

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0, 0], [0, 0]])
    e, f, r1, r2 = tuple(ambient.module_generators())
    marked = (e, 5*f+r1, r2)
    sub = ambient.subobject_on(marked)
    inclusion = sub.inclusion()
    assert inclusion.is_primitive()
    assert int(inclusion.domain().rank()) == 3
    basis = tuple(ambient.module_generators())
    gram = matrix(SageZZ, [[int(ambient.b(x, y)) for y in basis] for x in basis])
    generators = tuple(degenerate_marking_generators(ambient, marked, pointwise))
    for g in generators:
        assert abs(g.det()) == 1 and g.transpose()*gram*g == gram
        images = tuple(g*vector(SageZZ, list(x)) for x in (
            vector(SageZZ, [1, 0, 0, 0]),
            vector(SageZZ, [0, 5, 1, 0]),
            vector(SageZZ, [0, 0, 0, 1]),
        ))
        if pointwise:
            assert images == tuple(vector(SageZZ, coords) for coords in (
                (1,0,0,0), (0,5,1,0), (0,0,0,1)
            ))
        else:
            lattice_span = matrix(SageZZ, [[1,0,0,0],[0,5,1,0],[0,0,0,1]]).row_module()
            assert all(item in lattice_span for item in images)
    field = GF(5)
    residues = tuple(matrix(field, [[g[i,j] for j in (2,3)] for i in (2,3)]) for g in generators)
    one = identity_matrix(field,2)
    key = lambda m: tuple(int(c) for c in m.list())
    seen = {key(one)}
    frontier = [one]
    while frontier:
        value = frontier.pop()
        for generator in residues:
            moved = generator*value
            item = key(moved)
            if item not in seen:
                seen.add(item)
                frontier.append(moved)
            assert len(seen) <= expected_image
    assert len(seen) == expected_image
    if pointwise:
        assert all(int(g[3,2]) % 5 == 0 for g in generators)
        assert SageZZ.ideal([int(g[3,2]) // 5 for g in generators]) == SageZZ.ideal(1)


@pytest.mark.parametrize("pointwise,expected_residue_order", [(True, 6), (False, 12)])
def test_nonprimitive_proper_nondegenerate_projection_retains_perpendicular_shears(
    pointwise, expected_residue_order
):
    # C=<2>+<4>, R=<0>, and S=<5e+r>. The projected subspace is
    # the nondegenerate line Qe, but its integral image is 5Ze.
    # Any ambient isometry preserving S has
    # e -> s*e, r -> s*r, f -> t*f+b*r (s,t=+-1, b in Z).
    # For pointwise S, s=1; for moving S, both s are permitted.
    # Thus the reduction modulo 3 has exactly 2*3=6 and 2*2*3=12
    # elements respectively, and the shear parameter ideal is Z.
    from sage.rings.finite_rings.finite_field_constructor import GF

    ambient = Lattices(ZZ)([[2, 0, 0], [0, 4, 0], [0, 0, 0]])
    e, f, r = tuple(ambient.module_generators())
    v = 5*e+r
    s = ambient.subobject_on((v,))
    assert s.inclusion().is_primitive()
    generators = tuple(degenerate_marking_generators(ambient, (v,), pointwise))
    gram = matrix(SageZZ, [[2,0,0],[0,4,0],[0,0,0]])
    for g in generators:
        assert abs(g.det()) == 1 and g.transpose()*gram*g == gram
        assert g*vector(SageZZ,[5,0,1]) in (
            (vector(SageZZ,[5,0,1]),) if pointwise
            else (vector(SageZZ,[5,0,1]),vector(SageZZ,[-5,0,-1]))
        )
    assert SageZZ.ideal([int(g[2,1]) for g in generators]) == SageZZ.ideal(1)
    field = GF(3)
    reduced = tuple(g.change_ring(field) for g in generators)
    one = identity_matrix(field,3)
    key = lambda m: tuple(int(c) for c in m.list())
    seen = {key(one)}
    frontier = [one]
    while frontier:
        element = frontier.pop()
        for generator in reduced:
            moved = generator*element
            item = key(moved)
            if item not in seen:
                seen.add(item)
                frontier.append(moved)
            assert len(seen) <= expected_residue_order
    assert len(seen) == expected_residue_order


@pytest.mark.parametrize("pointwise,expected_mod_three_order", [(True, 6), (False, 12)])
def test_nondegenerate_projection_stabilizer_survives_integral_basis_change(
    pointwise, expected_mod_three_order
):
    # Change the C-basis e,f to a=e+f,b=f. The Gram becomes
    # [[6,4],[4,4]], and the marked vector is 5(a-b)+r.
    # The subgroup is conjugate by an integral basis change, hence
    # retains its exact image order modulo three.
    from sage.rings.finite_rings.finite_field_constructor import GF

    ambient = Lattices(ZZ)([[6,4,0],[4,4,0],[0,0,0]])
    a,b,r = tuple(ambient.module_generators())
    marked = 5*(a-b)+r
    family = tuple(degenerate_marking_generators(ambient, (marked,), pointwise))
    gram = matrix(SageZZ, [[6,4,0],[4,4,0],[0,0,0]])
    original_mark = vector(SageZZ,[5,-5,1])
    for g in family:
        assert abs(g.det()) == 1
        assert g.transpose()*gram*g == gram
        assert g*original_mark in ((original_mark,) if pointwise else (original_mark,-original_mark))

    field = GF(3)
    generators = tuple(g.change_ring(field) for g in family)
    one = identity_matrix(field,3)
    key = lambda g: tuple(map(int,g.list()))
    reached = {key(one)}
    pending = [one]
    while pending:
        element = pending.pop()
        for generator in generators:
            moved = generator*element
            value = key(moved)
            if value not in reached:
                reached.add(value)
                pending.append(moved)
            assert len(reached) <= expected_mod_three_order
    assert len(reached) == expected_mod_three_order
