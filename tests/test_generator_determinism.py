"""Two processes return the same automorphism-group generating family (issue #72).

The indefinite reduction, the Eichler envelope's hyperbolic-pair search and
the Lorentzian perturbations draw random matrices, as the pinned source does.
Each process starts with a fresh interpreter, so an unseeded draw gives a
different family in each process.
"""

import subprocess
import sys

# U(2) + U + <2> does not split U + U, so the recursion reduces it with the
# randomized indefinite reduction and builds the Eichler envelope.
CHILD = """
from sage.all import ZZ, matrix
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import lattice, rows_of

gram = matrix(ZZ, [[0, 2, 0, 0, 0], [2, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 2]])
generators = EichlerRecursion().orthogonal_group_generators(lattice(gram.rows()))
print([rows_of(generator) for generator in generators])
"""


def _family_in_fresh_process() -> str:
    completed = subprocess.run([sys.executable, "-c", CHILD], capture_output=True, text=True, check=True)
    return completed.stdout


def test_two_processes_return_the_same_generating_family():
    first = _family_in_fresh_process()
    second = _family_in_fresh_process()
    assert first
    assert first == second
