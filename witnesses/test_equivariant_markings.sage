"""Long witnesses moved from tests/test_equivariant_markings.sage."""
from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.integer_ring import ZZ as SageZZ
from sage_indefinite_port.indefinite.equivariant_markings import CentralizerTupleOrbit, _right_orbit_restriction, EquivariantMarking, _endomorphism_unit_generators, invisible_source_automorphism_generators, source_automorphism_kernel_generators


def _m(rows):
    return sage_matrix(SageZZ, [[int(c) for c in row] for row in rows])


UUTWO = [[0, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 2]]


JORDAN5 = _m([[0, -1, 0, 1, 0], [-1, 4, 0, -7, -4], [1, -7, 1, 6, 2], [0, 0, 0, 1, 0], [0, 2, 0, -1, -1]])


KER3 = _m([[1, 0, -2], [1, 0, 0], [0, 1, 0], [0, 0, 0], [0, 0, 1]])


KER3_ACTION = _m([[3, 0, -2], [-6, 1, 0], [2, 0, -1]])


def test_unipotent_transporter_beyond_the_linear_step():
    r"""F = JORDAN5 is unipotent with log F = x a single Jordan block of size 5, so A = Q[x] and A / rad A = Q.
    X = ker x^3 (columns (1,1,0,0,0), e3, (-2,0,0,0,1), primitive) is F-stable with singular Gram matrix, and
    x^2 X != 0, so exp(n x) X = X + n x X + n^2 x^2 X / 2 + ... is not linear in n. F X and F^-1 X lie in the
    orbit of X, with transporters in C(F)."""
    assert JORDAN5.transpose() * _m(UUTWO) * JORDAN5 == _m(UUTWO) and JORDAN5 * KER3 == KER3 * KER3_ACTION
    nilpotent = JORDAN5 - sage_identity(SageZZ, int(5))
    assert not (nilpotent ** 2 * KER3).is_zero()
    marking = EquivariantMarking(UUTWO, KER3, [JORDAN5], [KER3_ACTION])
    assert marking._is_unipotent_type()
    for target in (JORDAN5 * KER3, JORDAN5.inverse().change_ring(SageZZ) * KER3):
        g = marking.transporter(target)
        assert g is not None and g * KER3 == target and g * JORDAN5 == JORDAN5 * g
