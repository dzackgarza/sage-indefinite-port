r"""Commutant-group witnesses that need long runs (#54)."""

from sage_extensions.all import *
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import QuadraticField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.commutant_groups import centralizer_generators, equivariant_transporter

def _eichler():
    gram = sage_matrix(SageZZ, [[int(0), int(1), int(0)], [int(1), int(0), int(0)], [int(0), int(0), int(2)]])
    E = sage_matrix(SageZZ, [[int(1), int(-1), int(-2)], [int(0), int(1), int(0)], [int(0), int(1), int(1)]])
    return gram, E


def _is_transporter(t, gram, target_gram, sources, targets):
    return (
        t.det() in (int(1), int(-1))
        and t.transpose() * target_gram * t == gram
        and all(t * F == T * t for F, T in zip(sources, targets))
    )


def test_levi_transporter_to_the_inverse_on_one_summand():
    gram, E = _eichler()
    big = gram.block_sum(gram)
    F = E.block_sum(E)
    commutant = IsometryCommutant(big, [F])
    target = E.block_sum(E.inverse().change_ring(SageZZ))
    t = equivariant_transporter(commutant, big, [target])
    assert t is not None and _is_transporter(t, big, big, [F], [target])
    assert equivariant_transporter(commutant, int(3) * big, [F]) is None
