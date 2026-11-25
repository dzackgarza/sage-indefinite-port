"""Exact marked-transport oracles from complete binary integral isometry groups.

For U with Gram [[0,1],[1,0]], a row isometry [[a,b],[c,d]]
satisfies ac=bd=0, ad+bc=1. Thus O(U)={I,-I,S,-S}, where
S exchanges the basis vectors: either a,d are units and b=c=0,
or b,c are units and a=d=0.

For D=diag(2,-2), an integral row isometry [[a,b],[c,d]] obeys
a^2-b^2=1, c^2-d^2=-1, ac-bd=0. Factoring
(a-b)(a+b)=1 and (d-c)(d+c)=1 forces b=c=0 and
a,d in {1,-1}. Thus O(D) consists of the four diagonal signs.

In particular any integral invertible matrix preserves the gcd of the
coordinates of each vector: both its matrix and inverse are integral.
The tests below use these complete independently proved groups, not
values produced by the port, to certify positive and negative cases.
"""

from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from sage_indefinite_port.indefinite.marked_lattices import isometry
from sage_extensions.all import Lattices, ZZ
from tests.acceptance.consumer import lattice


def test_u_exact_pointed_transport_and_content_obstruction():
    plane = lattice(((0, 1), (1, 0)))
    e, f = tuple(plane.module_generators())
    recursion = EichlerRecursion()
    # S maps (2,4) to (4,2), and both have square 16.
    source, target = 2 * e + 4 * f, 4 * e + 2 * f
    witness = recursion.pointed_isometry(source, target)
    assert witness is not None and witness(source) == target
    # Same norm 16, but contents 2 and 1; none of the four isometries works.
    assert recursion.pointed_isometry(source, e + 8 * f) is None
    # Both vectors below are primitive of square 12. Nevertheless,
    # O(U) can only change both signs or exchange the coordinates.
    assert recursion.pointed_isometry(e + 6 * f, 2 * e + 3 * f) is None


def test_diagonal_binary_exact_pointed_positive_and_negative():
    diagonal = lattice(((2, 0), (0, -2)))
    e, f = tuple(diagonal.module_generators())
    recursion = EichlerRecursion()
    # The complete group consists of independent coordinate signs.
    x, y = 3 * e + 2 * f, -3 * e + 2 * f
    witness = recursion.pointed_isometry(x, y)
    assert witness is not None and witness(x) == y
    assert recursion.pointed_isometry(x, 3 * e - 2 * f) is not None
    # Norm and content are both necessary, but neither alone suffices
    # for an arbitrary lattice. Here the four signs are exhaustive.
    assert recursion.pointed_isometry(2 * e, e) is None
    assert recursion.pointed_isometry(2 * e, -2 * f) is None


def test_u_based_rank_one_marking_requires_exact_generator_transport():
    ambient = Lattices(ZZ)([[0, 1], [1, 0]])
    e, f = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0]])
    category = Lattices(ZZ).CosliceUnder(source)
    first = category(source.hom([2 * e], codomain=ambient))
    second = category(source.hom([2 * f], codomain=ambient))
    witness = isometry(category, first, second)
    assert witness is not None and witness(2 * e) == 2 * f
    # Content is an invariant of every ambient integral automorphism.
    third = category(source.hom([3 * f], codomain=ambient))
    assert isometry(category, first, third) is None
