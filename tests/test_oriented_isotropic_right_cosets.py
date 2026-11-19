"""Exact rational/integral right cosets of oriented primitive isotropic flags."""

from sage_extensions.all import Lattices, ZZ
from sage_extensions.categories.isotropic_orbits import FlagType, OrientedIsotropicFlag
from sage_extensions.categories.lattices import OrientedIsotropicFlags

from sage_indefinite_port.indefinite.marked_lattices import isotropic_stabilizer_right_cosets


def test_oriented_complete_flag_in_two_u_has_one_integral_right_coset():
    r"""For the complete flag <e> < <e,u> of 2U, its rational
    orientation-preserving Levi action on the plane consists of integral
    unitriangular matrices (block determinants +1). The pointwise
    unipotent radical has integral parameter, so every rational right
    coset in the finite lattice orbit has an integral representative;
    the quotient consists of one coset.
    """
    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U")
    e, _f, u, _v = tuple(ambient.module_generators())
    category = OrientedIsotropicFlags(lattices, 2)
    marked = category(OrientedIsotropicFlag(
        ambient, (e, u), flag_type=FlagType((1, 2)),
    ))
    representatives = tuple(isotropic_stabilizer_right_cosets(category, marked))
    assert len(representatives) == 1
    witness = representatives[0]
    rational = witness.domain()
    assert witness.codomain() is rational
    assert all(rational.b(witness(x), witness(y)) == rational.b(x, y)
               for x in rational.module_generators() for y in rational.module_generators())


def test_complete_oriented_flag_rational_levi_is_full_integral_unitriangular_group():
    r"""For the oriented complete flag e < (e,u) in 2U, the Levi
    action on its ordered null basis is precisely
    { [[1,0],[n,1]] : n in ZZ }. Both graded determinant signs are +1.
    Its generator images must contain shears with gcd of n equal to 1;
    merely filtering out negative-sign generators is insufficient.
    """
    from math import gcd

    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    from sage_indefinite_port.indefinite.marked_lattices import _oriented_flag_action

    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U")
    e, _f, u, _v = tuple(ambient.module_generators())
    category = OrientedIsotropicFlags(lattices, 2)
    marked = category(OrientedIsotropicFlag(
        ambient, (e, u), flag_type=FlagType((1, 2)),
    ))
    dimensions, basis, images = _oriented_flag_action(marked)
    assert dimensions == (1, 2) and len(basis) == 2
    top = ambient.primitive_isotropic_subobject(*basis).inclusion()
    labels = tuple(top.domain().module_generating_set())

    def rows(points):
        return matrix(SageZZ, [
            [int(top.lift(point).to_vector()(label)) for label in labels]
            for point in points
        ])

    frame_inverse = rows(basis).inverse()
    matrices = tuple(rows(images_of_basis) * frame_inverse for images_of_basis in images)
    assert matrices
    assert all(g[0, 0] == g[1, 1] == 1 and g[0, 1] == 0 and g.det() == 1
               for g in matrices)
    assert gcd(*(int(g[1, 0]) for g in matrices)) == 1


def test_oriented_rank_two_plane_levi_has_complete_sl_two_mod_two_image():
    r"""An oriented rank-two null plane has Levi SL_2(ZZ), not
    GL_2(ZZ). Its mod-two reduction is all GL_2(F_2), of order six.
    The complete integral kernel generators are obtained by Schreier
    reduction along the determinant-sign character of GL_2(ZZ).
    """
    from collections import deque

    from sage.all import GF
    from sage.matrix.constructor import identity_matrix, matrix
    from sage.rings.integer_ring import ZZ as SageZZ

    from sage_indefinite_port.indefinite.marked_lattices import _oriented_flag_action

    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U")
    e, _f, u, _v = tuple(ambient.module_generators())
    category = OrientedIsotropicFlags(lattices, 1)
    marked = category(OrientedIsotropicFlag(ambient, (e, u)))
    dimensions, basis, images = _oriented_flag_action(marked)
    assert dimensions == (2,)
    top = ambient.primitive_isotropic_subobject(*basis).inclusion()
    labels = tuple(top.domain().module_generating_set())

    def rows(points):
        return matrix(SageZZ, [
            [int(top.lift(point).to_vector()(label)) for label in labels]
            for point in points
        ])

    inverse = rows(basis).inverse()
    generators = tuple((rows(family) * inverse).change_ring(SageZZ)
                       for family in images)
    assert generators and all(generator.det() == 1 for generator in generators)
    field = GF(2)
    images_mod_two = tuple(generator.change_ring(field) for generator in generators)
    identity = identity_matrix(field, 2)
    key = lambda element: tuple(int(entry) for entry in element.list())
    reached = {key(identity)}
    pending = deque((identity,))
    while pending:
        current = pending.popleft()
        for generator in images_mod_two:
            moved = current * generator
            if key(moved) not in reached:
                reached.add(key(moved))
                pending.append(moved)
            assert len(reached) <= 6
    assert len(reached) == 6
