"""How long does one _gram_in_label_order call take on the rank-22 K3 rational envelope?"""
import time

from sage.all import ZZ, CartanMatrix, block_diagonal_matrix, matrix

from sage_indefinite_port.invariants import AttackProfile
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from sage_indefinite_port.indefinite.groupoid_cache import _gram_in_label_order
from tests.acceptance.consumer import lattice

hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
e8 = -matrix(ZZ, CartanMatrix(["E", 8]))
gram = block_diagonal_matrix([hyperbolic, hyperbolic, hyperbolic, e8, e8])
L = lattice([list(row) for row in gram.rows()])
model = EichlerRecursion()._model(AttackProfile.from_lattice(L).signed_view)
envelope = model.envelope
rational_envelope = envelope.rational_isometry.codomain()
print(f"envelope rank {rational_envelope.module_rank()}, base ring {rational_envelope.base_ring()}", flush=True)

for attempt in range(2):
    started = time.monotonic()
    gram_matrix = _gram_in_label_order(rational_envelope)
    print(f"call {attempt + 1}: {time.monotonic() - started:.2f}s  shape {gram_matrix.nrows()}x{gram_matrix.ncols()}",
          flush=True)
