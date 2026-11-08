"""Does the Gram memo make repeated calls cheap, and is the value unchanged?"""
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
rational_envelope = model.envelope.rational_isometry.codomain()

timings = []
values = []
for attempt in range(3):
    started = time.monotonic()
    values.append(_gram_in_label_order(rational_envelope))
    timings.append(time.monotonic() - started)
print("timings:", [f"{t:.2f}s" for t in timings], flush=True)
print("all equal:", all(v == values[0] for v in values), flush=True)
print("same object:", all(v is values[0] for v in values), flush=True)

# an integral lattice too, to exercise the ZZ branch
small = lattice([[2, 0], [0, 0]])
a = _gram_in_label_order(small)
b = _gram_in_label_order(small)
print("integral lattice: equal", a == b, "| same object", a is b, "| value", a.rows(), flush=True)
