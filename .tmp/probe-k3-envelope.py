"""Does the K3 covering go through the Eichler envelope -- i.e. through BOUND's finite representation?

primitive_covering_vectors (eichler.py:980-995) takes the envelope branch when self.envelope is not None, and
that branch immediately asks envelope.right_coset_transversal(), which asks the lattice stabilizer's finite
representation. If K3 has an envelope, the A_VECTORS K3 cost and BOUND's finite-representation obligation are
the same mathematical locus. Cheap to check: the model builds in ~9 s.
"""
import time

from sage.all import ZZ, CartanMatrix, block_diagonal_matrix, matrix

from sage_indefinite_port.invariants import AttackProfile
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import lattice

hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
e8 = -matrix(ZZ, CartanMatrix(["E", 8]))
gram = block_diagonal_matrix([hyperbolic, hyperbolic, hyperbolic, e8, e8])
L = lattice([list(row) for row in gram.rows()])
profile = AttackProfile.from_lattice(L)
model = EichlerRecursion()._model(profile.signed_view)
print(f"model built: envelope is None? {model.envelope is None}", flush=True)
if model.envelope is not None:
    envelope = model.envelope
    print(f"  similarity: {type(envelope.similarity).__name__}", flush=True)
    started = time.monotonic()
    exponent = envelope.similarity.cokernel().exponent()
    print(f"  exponent = {exponent} in {time.monotonic() - started:.2f}s", flush=True)
