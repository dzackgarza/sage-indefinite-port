"""Setwise transport under a primitive degenerate rank-two marking."""

from sage_extensions.all import Lattices, ZZ
from sage_indefinite_port.indefinite.pairing_slices import setwise_transporter


def test_glued_mixed_degenerate_setwise_marking_allows_hyperbolic_pair_swap():
    r"""L=U + <0> and S=<e+f,r> are primitive, with Gram
    diag(2,0). Its orthogonal complement is <e-f>, and
    S + S^perp has index two in L. The integral isometry
    exchanging e and f and fixing r preserves S and sends
    e to f. A singular Gram matrix for S must not prevent
    that exact, glue-compatible transporter.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    g = setwise_transporter(ambient, (e + f, r), e, f)
    assert g is not None
    assert g(e) == f
    assert all(ambient.b(g(x), g(y)) == ambient.b(x, y)
               for x in (e, f, r) for y in (e, f, r))
    sub = ambient.subobject_on((e + f, r)).inclusion()
    assert all(sub.is_in_image(g(x)) for x in (e + f, r))


def test_glued_mixed_setwise_transporter_decides_radical_mod_two_obstruction():
    r"""For S=<e+f,r> in U+<0>, a moving vector x=2e+r
    must retain an odd radical coordinate modulo two under every
    stabilizer: the shear from 2e is even and GL_1(ZZ)={+-1}.
    Hence 2f has no transporter, while 2f+3r does.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    fixed = (e + f, r)
    assert setwise_transporter(ambient, fixed, 2*e + r, 2*f) is None
    g = setwise_transporter(ambient, fixed, 2*e + r, 2*f + 3*r)
    assert g is not None and g(2*e + r) == 2*f + 3*r
    inclusion = ambient.subobject_on(fixed).inclusion()
    assert all(inclusion.is_in_image(g(point)) for point in fixed)


def test_glued_mixed_setwise_transporter_pure_radical_GL_action():
    r"""At zero quotient projection, the full-radical marked
    stabilizer acts on R by GL(R). For a rank-one radical the
    two primitive vectors r and -r are equivalent, but r and
    2r are not: their contents differ.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    fixed = (e + f, r)
    g = setwise_transporter(ambient, fixed, r, -r)
    assert g is not None and g(r) == -r
    assert setwise_transporter(ambient, fixed, r, 2*r) is None


def test_glued_mixed_setwise_transporter_uses_actual_gl_two_mod_five():
    r"""For L=U+<0>^2 and S=<e+f,r1,r2>, the full
    radical GL_2(ZZ) image modulo 5 sends (1,0) to (2,0),
    although rank-one GL_1(ZZ) would not. A lift is
    D=[[2,5],[5,13]] with determinant 1. A projected
    5e shear then carries r1 D to the target radical row.
    The zero residue (0,0) cannot reach (0,1).
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0,0],[0,0]])
    e, f, r1, r2 = tuple(ambient.module_generators())
    fixed = (e+f,r1,r2)
    g = setwise_transporter(ambient, fixed, 5*e+r1, 5*f+2*r1)
    assert g is not None and g(5*e+r1) == 5*f+2*r1
    inclusion = ambient.subobject_on(fixed).inclusion()
    assert all(inclusion.is_in_image(g(point)) for point in fixed)
    assert setwise_transporter(ambient, fixed, 5*e+5*r1, 5*f+r2) is None


def test_full_radical_proper_projection_stabilizer_generators_are_complete():
    r"""The full stabilizer of S=<e+f,r> in U+<0> is

        Hom_Z(U,Zr) semidirect (O(U) x GL_1(Z)).

    Every O(U) element preserves the unoriented line <e+f>,
    giving a quotient of order 4*2=8. The kernel is Z^2 of
    independent shears e -> e+ar, f -> f+br, so its integral
    lattice is generated iff the displayed shear parameters
    generate Z^2.
    """
    from collections import deque
    from math import gcd

    from sage.matrix.constructor import identity_matrix, matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    ambient = Lattices(ZZ)("U") + Lattices(ZZ)([[0]])
    e, f, r = tuple(ambient.module_generators())
    inclusion = ambient.subobject_on((e+f,r)).inclusion()
    group = ambient.O().setwise_stabilizer(inclusion)
    generators = tuple(group.group_generators())
    assert generators
    labels = tuple(ambient.module_generating_set())

    def row(g, point):
        coordinates = g(point).to_vector()
        return tuple(int(coordinates(label)) for label in labels)

    shears = []
    finite = []
    for g in generators:
        assert all(inclusion.is_in_image(g(point)) for point in (e+f,r))
        A = matrix(SageZZ, [row(g,e)[:2], row(g,f)[:2]])
        d = row(g,r)[-1]
        assert abs(A.det()) == abs(d) == 1
        finite.append((A, d))
        if A == identity_matrix(SageZZ,2) and d == 1:
            shears.append((row(g,e)[-1], row(g,f)[-1]))
    assert shears
    assert gcd(*(a*b2-a2*b for a,b in shears for a2,b2 in shears)) == 1

    def key(A,d):
        return tuple(int(entry) for entry in A.list()) + (int(d),)

    identity = identity_matrix(SageZZ,2)
    reached = {key(identity,1)}
    frontier = deque(((identity,1),))
    while frontier:
        A,d = frontier.popleft()
        for generator, sign in finite:
            moved_A,moved_d = A*generator, d*sign
            label = key(moved_A,moved_d)
            if label not in reached:
                reached.add(label)
                frontier.append((moved_A,moved_d))
            assert len(reached) <= 8
    assert len(reached) == 8
