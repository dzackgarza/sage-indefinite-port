r"""A negative Hermitian transporter over ZZ[2i], and a parabolic meet a joint finite-image preimage.

Negative transporter. Let O = ZZ[2i], with ZZ-basis 1, 2i and units +-1, and
pi = 1 + 2i, so pi pi^bar = 5. On O^3 with h(x, y) = x G y^* and

    G = [[0, 5, pi], [5, 0, pi^bar], [pi^bar, pi, 0]],

b1 and b2 are null and, as basis vectors, span primitive free lines. The
pairing ideals are h(b1, L) = (5, pi) = pi O and h(b2, L) = (5, pi^bar) = pi^bar O.
They differ: (pi^bar / pi) = (-3 - 4i)/5 is not in ZZ[i]. An isometry g with
g(O b1) = O b2 sends b1 to +-b2, so h(b2, L) = h(g b1, g L) = h(b1, L), which is
false; no transporter exists, setwise or based.

The quotients are isometric. h(x, b1) = 5 x2 + pi^bar x3, so b1^perp = O b1 + O u1
with u1 = (0, 1, -pi), as 5 / pi^bar = pi. Then h(u1, u1) = -(pi^2 + pi^bar^2) = 6.
Likewise h(x, b2) = 5 x1 + pi x3 gives u2 = (1, 0, -pi^bar) with h(u2, u2) = 6.
Both quotients are the free rank-one lattice <6>, and det G = -30 = -6 * N(pi).
So the search reaches the lifting step: each candidate (a, u) in {+-1}^2 is an
isometry of the quotients, and none lifts.

Joint preimage. Over ZZ[i] on ZZ[i]^3 with J = antidiag(1,1,1), the setwise
stabilizer P of ZZ[i] e1 is l(a,u) n(Z,T) with a, u in mu_4, Z in (1+i)ZZ[i]
and T in iZZ. Then g e2 = (-a Z^bar, u, 0) and det g = a^2 u. The map
rho(g) = (g mod 2, det g) is a homomorphism to a finite group, and S is the
subgroup of pairs fixing e2 modulo 2 with determinant one. On P,
(a, u, Z mod 2) is a homomorphism onto mu_4 x mu_4 x ZZ/2, since units act
trivially on (1+i)ZZ[i] / 2ZZ[i]. The conditions are u = +-1, Z in 2ZZ[i] and
a^2 u = 1, so u = a^-2 and the meet is the graph of order 4: index 32 / 4 = 8.
Fixing e1 forces a = 1, so u = 1 and Z in 2ZZ[i]: index 8 in mu_4 x ZZ/2.
"""

import pytest

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.number_field.number_field import CyclotomicField as SageCyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice, HermitianWittFrame


def _assert_unitary_integral(lattice, g):
    form = lattice.gram_matrix().transpose()
    assert g.conjugate_transpose() * form * g == form
    for entry in (g.list() + (~g).list()):
        lattice._coordinates_in_order(entry)


def _two_i_pairing_lattice():
    field = SageCyclotomicField(int(4))
    i = field.gen()
    pi = field.one() + int(2) * i
    bar = pi.conjugate()
    gram = sage_matrix(field, int(3), int(3), [int(0), int(5), pi, int(5), int(0), bar, bar, pi, int(0)])
    return field, HermitianLattice((field.one(), int(2) * i), gram)


def test_two_i_lines_with_isometric_quotients_have_no_transporter():
    field, lattice = _two_i_pairing_lattice()
    e = lambda k: [int(1) if j == k else int(0) for j in range(int(3))]
    first, second = lattice.witt_frame([e(int(0))]), lattice.witt_frame([e(int(1))])
    assert len(first.quotient_isometries_to(second)) == int(2)
    assert lattice.isotropic_transporter([e(int(0))], [e(int(1))]) is None
    assert lattice.isotropic_transporter([e(int(0))], [e(int(1))], based=True) is None
    _assert_unitary_integral(lattice, lattice.isotropic_transporter([e(int(0))], [e(int(0))]))


def _gaussian_line_frame():
    field = SageCyclotomicField(int(4))
    gram = sage_matrix(field, int(3), int(3), [int(0), int(0), int(1), int(0), int(1), int(0), int(1), int(0), int(0)])
    lattice = HermitianLattice((field.one(), field.gen()), gram)
    return field, lattice, HermitianWittFrame(lattice, gram.parent().identity_matrix(), sage_matrix(field, int(1), int(1), [int(1)]), int(1))


@pytest.mark.parametrize("pointwise", (False, True), ids=("setwise", "pointwise"))
def test_gaussian_line_stabilizer_meets_a_joint_residue_and_determinant_preimage(pointwise):
    field, lattice, frame = _gaussian_line_frame()
    identity = frame.residue_coordinates(sage_matrix.identity(field, int(3)), int(2))
    # Entry (row, 1) occupies the two order coordinates starting at 2 * (3 * row + 1).
    column = [position for row in range(int(3)) for position in (int(2) * (int(3) * row + int(1)), int(2) * (int(3) * row + int(1)) + int(1))]
    image = lambda g: (frame.residue_coordinates(g, int(2)), g.det())
    in_subgroup = lambda pair: pair[int(1)] == field.one() and all(pair[int(0)][p] == identity[p] for p in column)
    generators, index = frame.preimage_stabilizer_generators(image, in_subgroup, pointwise)
    assert index == int(8)
    e1 = sage_matrix(field, int(3), int(1), [int(1), int(0), int(0)])
    for g in generators:
        _assert_unitary_integral(lattice, g)
        assert in_subgroup(image(g))
        moved = g * e1
        assert moved[int(1), int(0)] == moved[int(2), int(0)] == int(0)
        if pointwise:
            assert moved == e1
