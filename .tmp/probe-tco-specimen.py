"""Current-pin reproduction of the T_Co discriminant-surjectivity handoff.

Specimen copied verbatim from the research artifact at
dzackgarza/research@a145cea4 computations/reports/coble_tco_surjectivity_specimen.py:
the Gram is <2> + U(2) + E8(2), and `Lattices(ZZ)` is given [list(row) for row in GRAM.rows()],
not a Sage matrix. Both same-session controls are kept. Nothing is repaired here.

Old-pin record (provider c546077d, Sage 11.0.beta5): raised TypeError at the engine's
recursive.py:154 `int(coordinates(label))` with an OwnedJoinCategory.element_class.
"""
import hashlib
import inspect
import os
import sys
import traceback

from sage.all import CartanMatrix, ZZ, block_diagonal_matrix, matrix
from sage.version import version as sage_version

import sage_extensions
import sage_indefinite_port
from sage_extensions.all import Lattices

sys.stdout.reconfigure(line_buffering=True)

provider_root = inspect.getfile(sage_extensions).rsplit("/", 1)[0]
port_root = inspect.getfile(sage_indefinite_port).rsplit("/", 1)[0]
print("provenance")
print("  sage-extensions module :", inspect.getfile(sage_extensions))
print("  sage-extensions pin    :", os.environ.get("PROVIDER_PIN", "unset"))
print("  port module           :", inspect.getfile(sage_indefinite_port))
print("  port engine revision  :", os.environ.get("ENGINE_REV", "unset"))
print("  Sage                  :", sage_version)
for rel in ("categories/rational_integral_stabilizers.py", "categories/lattices.py"):
    path = f"{provider_root}/{rel}"
    print(f"  {rel:44s} sha256 {hashlib.sha256(open(path, 'rb').read()).hexdigest()[:16]}")
print()

U = matrix(ZZ, [[0, 1], [1, 0]])
E8 = -CartanMatrix(["E", 8])
GRAM = block_diagonal_matrix(matrix(ZZ, [[2]]), 2 * U, 2 * E8)

L = Lattices(ZZ)([list(row) for row in GRAM.rows()])
print("construction: Lattices(ZZ)(Gram of <2> + U(2) + E8(2))")
print("  rank %d, determinant %s, signature %s, even %s"
      % (L.rank(), L.determinant(), L.signature_pair(), L.is_even()))
print("call: L.discriminant_representation_is_surjective()")
status = "raised"
try:
    value = L.discriminant_representation_is_surjective()
    status = f"returned {value}"
    print("  returned", value)
except Exception:
    print("  raised; full traceback follows")
    print("-" * 78)
    traceback.print_exc(file=sys.stdout)
    print("-" * 78)

print()
print("control specimens, same call, same session:")
for name, gram in (
    ("S_Co = <2> + <-2>^10", block_diagonal_matrix(matrix(ZZ, [[2]]), -2 * matrix.identity(ZZ, 10))),
    ("M = S_En + <delta> = U(2) + E8(2) + <-2>",
     block_diagonal_matrix(matrix(ZZ, [[-2]]), 2 * U, 2 * E8)),
):
    control = Lattices(ZZ)([list(row) for row in gram.rows()])
    try:
        print("  %-42s -> %s" % (name, control.discriminant_representation_is_surjective()))
    except Exception as exc:
        print("  %-42s -> raised %s" % (name, type(exc).__name__))
print()
print("RESULT:", status)
