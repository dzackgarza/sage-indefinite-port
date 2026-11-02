"""Where does the K3 vector-class work actually spend its time? One stage at a time, printed as it starts.

The call graph says: lattice -> AttackProfile.from_lattice (signature_pair) -> _model -> primitive_covering_vectors
-> pointed_isometry dedup. That is a call graph, not an attribution. This prints BEFORE each stage, so a cap
still leaves the attribution, and each stage is bounded by the outer timeout.
"""
import itertools
import time

START = time.monotonic()
def stage(label):
    print(f"[t={time.monotonic() - START:7.1f}s] {label}", flush=True)

stage("importing")
from sage.all import ZZ, CartanMatrix, block_diagonal_matrix, matrix
from sage_indefinite_port.invariants import AttackProfile
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import lattice

stage("building the Gram")
hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
e8 = -matrix(ZZ, CartanMatrix(["E", 8]))
gram = block_diagonal_matrix([hyperbolic, hyperbolic, hyperbolic, e8, e8])

stage("constructing the lattice")
L = lattice([list(row) for row in gram.rows()])
stage(f"lattice constructed: rank {L.rank()}")

stage("signature_pair()")
signature = L.signature_pair()
stage(f"signature_pair = {signature}")

stage("AttackProfile.from_lattice")
profile = AttackProfile.from_lattice(L)
stage(f"profile: positive_index {profile.positive_index}, sign {profile.sign}")

stage("EichlerRecursion()")
recursion = EichlerRecursion()
stage("building the model on the signed view")
model = recursion._model(profile.signed_view)
stage("model built")

for square in (-2, 2):
    stage(f"cover enumeration for square {square} (first 200 vectors)")
    count = 0
    for _ in itertools.islice(model.primitive_covering_vectors(profile.sign * square), 200):
        count += 1
    stage(f"square {square}: {count} cover vectors in the first 200-vector slice")
