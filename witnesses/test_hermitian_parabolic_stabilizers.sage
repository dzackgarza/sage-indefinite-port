r"""Hermitian parabolic witnesses whose stabilizer computations take minutes; small cases run locally."""


from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as _sage_identity
from sage.rings.number_field.number_field import CyclotomicField as SageCyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice, HermitianWittFrame
from sage.rings.integer_ring import ZZ as SZZ
from sage.rings.rational_field import QQ as SQQ
from sage.rings.integer_ring import ZZ as SZZ
from sage.rings.rational_field import QQ as SQQ


def _gaussian():
    field = SageCyclotomicField(4)
    i = field.gen()
    one, zero = field.one(), field.zero()
    gram = sage_matrix(field, 3, 3, [zero, zero, one, zero, one, zero, one, zero, zero])
    return field, i, HermitianLattice((one, i), gram), sage_matrix(field, 1, 1, [one])


def _frame(field, lattice, middle, first, last):
    half = sage_matrix(field, 3, 3, [first, field.zero(), field.zero(), field.zero(), field.one(), field.zero(), field.zero(), field.zero(), last])
    return HermitianWittFrame(lattice, half, middle, int(1))


def _parameters(field, Z, T):
    return sage_matrix(field, 1, 1, [Z]), sage_matrix(field, 1, 1, [T])


def _assert_unitary_integral(lattice, g):
    form = lattice.gram_matrix().transpose()
    assert g.conjugate_transpose() * form * g == form
    assert all(entry.is_integral() for entry in g.list())
    assert all(entry.is_integral() for entry in (~g).list())


def _vector(field, *entries):
    return sage_matrix(field, len(entries), int(1), [field(int(x)) if isinstance(x, int) else x for x in entries])


def _split_rank_two():
    field = SageCyclotomicField(int(4))
    zero, one = field.zero(), field.one()
    J = sage_matrix(field, int(5), int(5), [zero, zero, zero, one, zero, zero, zero, zero, zero, one, zero, zero, one, zero, zero, one, zero, zero, zero, zero, zero, one, zero, zero, zero])
    lattice = HermitianLattice((one, field.gen()), J)
    return field, lattice, HermitianWittFrame(lattice, J.parent().identity_matrix(), sage_matrix(field, int(1), int(1), [one]), int(2))


def _glued_eisenstein():
    field = SageCyclotomicField(int(3))
    omega = field.gen()
    theta = omega - omega.conjugate()
    middle = field(int(9))
    zero, one = field.zero(), field.one()
    witt = sage_matrix(field, int(3), int(3), [zero, zero, one, zero, middle, zero, one, zero, zero])
    basis = sage_matrix(field, int(3), int(3), [one, theta**int(-1), zero, zero, theta**int(-1), zero, zero, zero, theta])
    gram = basis.conjugate_transpose() * witt * basis
    lattice = HermitianLattice((one, omega), gram.transpose())
    return field, theta, lattice, HermitianWittFrame(lattice, ~basis, sage_matrix(field, int(1), int(1), [middle]), int(1))


def test_flag_stabilizers_based_and_setwise_on_the_split_frame():
    r"""Split Z[i]^5, flag <e1> < <e1,e2>. Pointwise: A = 1, and all four of U(1) lift.

    Setwise: A is upper triangular with unit diagonal and some A has a nonzero (1,2) entry.
    Every unipotent n(Z,T) fixes e1 and e2, so it lies in both.
    """
    field, lattice, frame = _split_rank_two()
    one, zero = field.one(), field.zero()
    e1, e2 = _vector(field, one, zero, zero, zero, zero), _vector(field, zero, one, zero, zero, zero)
    assert len(frame.pointwise_levi_image()) == int(4)
    pointwise = frame.stabilizer_generators(True, (int(1), int(1)))
    for g in pointwise:
        _assert_unitary_integral(lattice, g)
        assert g * e1 == e1 and g * e2 == e2
    assert {frame.induced_action(g)[1][0, 0] for g in pointwise} == {field.gen()**int(k) for k in range(int(4))}
    units = [field.gen()**int(k) for k in range(int(4))]
    setwise = frame.stabilizer_generators(False, (int(1), int(1)))
    actions = [frame.induced_action(g)[0] for g in setwise]
    assert all(A[1, 0] == zero and A[0, 0] in units and A[1, 1] in units for A in actions)
    assert any(A[0, 1] != zero for A in actions)
    for g in setwise:
        _assert_unitary_integral(lattice, g)


def test_setwise_line_stabilizer_of_the_glued_lattice_induces_the_hand_levi_image():
    r"""A pair lifts iff a = u mod theta, so the generators induce exactly 18 of the 36 unit pairs."""
    field, theta, lattice, frame = _glued_eisenstein()
    generators = frame.stabilizer_generators(False, (int(1),))
    pairs = set(tuple(a[0, 0] for a in frame.induced_action(g)) for g in generators)
    assert all(((a - u) / theta).is_integral() for a, u in pairs)
    known, frontier = set([(field.one(), field.one())]), [(field.one(), field.one())]
    while frontier:
        a, u = frontier.pop()
        for g, h in pairs:
            image = (a * g, u * h)
            if image not in known:
                known.add(image)
                frontier.append(image)
    assert len(known) == int(18)


def test_act_eisenstein_cusp_has_the_full_split_levi_image():
    r"""E^{4,1} = H + I_3 over Z[omega], with rho in the hyperbolic plane H.

    H is unimodular, so it splits off and every pair (a, u) in mu_6 x U(I_3) lifts.
    U(I_3) is the monomial group mu_6 wr S_3, of order 6^3 * 6 = 1296, so the setwise Levi
    image of the line Z[omega] rho has order 6 * 1296 = 7776.
    """
    field = SageCyclotomicField(int(3))
    zero, one = field.zero(), field.one()
    size = int(5)
    witt = sage_matrix(field, size, size, 0)
    witt[int(0), int(4)] = witt[int(4), int(0)] = one
    for k in range(int(1), int(4)):
        witt[k, k] = one
    lattice = HermitianLattice((one, field.gen()), witt)
    middle = _sage_identity(field, int(3))
    frame = HermitianWittFrame(lattice, witt.parent().identity_matrix(), middle, int(1))
    assert len(frame.levi_candidates()) == int(7776)
    assert len(frame.levi_image()) == int(7776)
