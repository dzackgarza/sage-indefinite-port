"""Compact cyclotomic unitary centralizers on indefinite trace lattices.

The scalar form (zeta_5+zeta_5^-1)(x_1 * conjugate(y_1) +
x_2 * conjugate(y_2)) is definite at each real place of Q(sqrt5),
but its two signatures are opposite. Its rational trace lattice has
signature (4,4), so ordinary definite-lattice automorphism computation
does not apply. The integral unitary group nevertheless has order 200.
"""

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
