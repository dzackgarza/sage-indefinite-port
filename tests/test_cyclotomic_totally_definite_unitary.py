"""Compact and noncompact cyclotomic unitary centralizers on indefinite trace lattices.

The scalar form (zeta_5+zeta_5^-1)(x_1 * conjugate(y_1) +
x_2 * conjugate(y_2)) is definite at each real place of Q(sqrt5),
but its two signatures are opposite. Its rational trace lattice has
signature (4,4), so ordinary definite-lattice automorphism computation
does not apply. The integral unitary group nevertheless has order 200.
"""

import pytest

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.groups.matrix_gps.finitely_generated import MatrixGroup
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.number_field.number_field import CyclotomicField

from sage_extensions.all import Lattices, ZZ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant


def _mixed_signature_cyclotomic():
    field = CyclotomicField(5)
    zeta = field.gen()
    scale = zeta + zeta**(-1)
    basis = tuple(zeta**i for i in range(4))
    base_gram = matrix(SageZZ, [
        [(scale * x * y.conjugate()).trace() for y in basis] for x in basis
    ])
    multiplication = matrix(SageZZ, [list(zeta**(i+1)) for i in range(4)]).transpose()
    gram = matrix(SageZZ, 8, 8)
    action = matrix(SageZZ, 8, 8)
    for pos in (0, 4):
        gram.set_block(pos, pos, base_gram)
        action.set_block(pos, pos, multiplication)
    lattice = Lattices(ZZ)([list(row) for row in gram.rows()])
    frame = tuple(lattice.module_generators())
    automorphism = lattice.O()(tuple(
        sum((int(action[i, j]) * frame[i] for i in range(8) if action[i, j]), lattice.zero())
        for j in range(8)
    ))
    return lattice, automorphism, gram, action


def test_mixed_archimedean_signs_have_the_complete_finite_unitary_centralizer():
    """U_2(Z[zeta_5],(zeta+zeta^-1)I) is mu_10 wr C_2, of order 200."""
    lattice, rotation, gram, action = _mixed_signature_cyclotomic()
    assert tuple(int(index) for index in lattice.signature_pair()) == (4, 4)
    assert action.transpose() * gram * action == gram
    assert action**5 == identity_matrix(SageZZ, 8)
    commutant = IsometryCommutant(gram, (action,))
    assert commutant.is_compact() is True
    exact = commutant.compact_centralizer_elements()
    assert len(exact) == 200

    group = lattice.O()
    generators = tuple(group.centralizer(rotation).group_generators())
    assert generators
    assert all(g.domain() is lattice and g.codomain() is lattice for g in generators)
    assert all(g * rotation == rotation * g for g in generators)
    labels = tuple(lattice.module_generating_set())
    matrices = [matrix(SageZZ, [
        [int(g(v).to_vector()(label)) for v in lattice.module_generators()] for label in labels
    ]) for g in generators]
    assert MatrixGroup([M.change_ring(GF(3)) for M in matrices]).order() == 200


def test_public_quartic_cm_hermitian_unitary_group_has_200_integral_elements():
    r"""The indefinite trace of (zeta_5+zeta_5^-1)I_2 has compact U(h) of order 200.

    The exact unitary generators are now available directly from the
    Hermitian lattice, without supplying a root-of-unity scalar operator
    or pretending the rational trace lattice is definite.
    """
    from sage_extensions.categories.hermitian_lattices import HermitianLattice

    field = CyclotomicField(5)
    zeta = field.gen()
    scale = zeta + zeta**(-1)
    hermitian = HermitianLattice(
        tuple(zeta**k for k in range(4)),
        matrix(field, 2, 2, [scale, 0, 0, scale]),
    )
    lattice = hermitian.trace_lattice()
    assert tuple(int(i) for i in lattice.signature_pair()) == (4, 4)
    generators = tuple(hermitian.unitary_group_generators())
    assert generators and all(hermitian.is_order_linear(g) for g in generators)
    labels = tuple(lattice.module_generating_set())
    matrices = [matrix(SageZZ, [
        [int(g(v).to_vector()(label)) for v in lattice.module_generators()]
        for label in labels
    ]) for g in generators]
    assert MatrixGroup([M.change_ring(GF(3)) for M in matrices]).order() == 200


def test_indefinite_eisenstein_binary_unitary_has_complete_public_generators():
    r"""The order-three centralizer on A2+A2(-1) is the full infinite U(1,1;O_E).

    The arithmetic commutant model gives seven integral generators,
    including a nonidentity unipotent of infinite order. Their actual
    mod-3 image is a group of order 108. The public O(L)-centralizer must
    expose the complete integer generators, not the finite mod-3 image
    or an unconditional refusal of an indefinite Hermitian component.
    """
    lattice = Lattices(ZZ)([[2, -1], [-1, 2]]) + Lattices(ZZ)([[-2, 1], [1, -2]])
    labels = tuple(lattice.module_generating_set())
    a1, a2, b1, b2 = lattice.module_generators()
    rotation = lattice.O()({
        labels[0]: a2, labels[1]: -a1 - a2,
        labels[2]: b2, labels[3]: -b1 - b2,
    })
    generators = tuple(lattice.O().centralizer(rotation).group_generators())
    assert generators
    assert all(g.domain() is lattice and g.codomain() is lattice for g in generators)
    assert all(g * rotation == rotation * g for g in generators)

    matrices = [matrix(SageZZ, [
        [int(g(v).to_vector()(label)) for v in lattice.module_generators()]
        for label in labels
    ]) for g in generators]
    assert MatrixGroup([M.change_ring(GF(3)) for M in matrices]).order() == 108
    one = identity_matrix(SageZZ, 4)
    # A nontrivial unipotent element has infinite order over ZZ:
    # (1+N)^k = 1 + kN for N^2=0 and k != 0.
    assert any(M != one and (M - one)**2 == 0 for M in matrices)


@pytest.mark.parametrize("middle_rank", (1, 2))
def test_eisenstein_hyperbolic_unitary_centralizer_generators(middle_rank):
    r"""U(H+I_m; Z[omega]) has complete public generators for m=1,2.

    The port's source-faithful Hermitian unitary algorithm constructs the
    full group through null-vector reduction and a finite monomial Levi
    action. Its trace-lattice orthogonal centralizer must use this route
    for the rank-six and rank-eight cases, not only the hand-normalized
    binary hyperbolic or rank-ten ACT routes.
    """
    from sage_indefinite_port.indefinite.eisenstein_hyperbolic_unitary import (
        OMEGA, lattice as hermitian_lattice,
    )

    hermitian = hermitian_lattice(middle_rank)
    lattice = hermitian.trace_lattice()
    scalar = hermitian.order_action(OMEGA)
    rotation = lattice.O()(tuple(scalar(x) for x in lattice.module_generators()))
    generators = tuple(lattice.O().centralizer(rotation).group_generators())
    assert generators
    assert all(g.domain() is lattice and g.codomain() is lattice for g in generators)
    assert all(g * rotation == rotation * g for g in generators)

    labels = tuple(lattice.module_generating_set())
    matrices = [matrix(SageZZ, [
        [int(g(v).to_vector()(label)) for v in lattice.module_generators()]
        for label in labels
    ]) for g in generators]
    one = identity_matrix(SageZZ, 2 * (middle_rank + 2))
    assert any(M != one and (M - one)**3 == 0 for M in matrices)
