"""Is GAP's matrix-group recognition affordable at degree 10 over F_5?

Recognition would identify the image as a classical group, which fixes its order from the classical formula
and avoids building a BSGS on F_p^10. A small control (GL_2(5)) is run first, so that a failure on the real
group is a cost statement rather than an API mistake.
"""
import time

from sage.all import GF, MatrixGroup, matrix, libgap

from tests.acceptance.consumer import lattice, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases


def phase(label, fn):
    started = time.monotonic()
    try:
        value = fn()
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:200]} "
              f"after {time.monotonic() - started:.2f}s", flush=True)
        return None
    print(f"[{label}] {value} in {time.monotonic() - started:.2f}s", flush=True)
    return value


print("available recognition entry points:", flush=True)
for name in ("IsomorphismTypeOfMatrixGroup", "IsomorphismTypeOfMatrixGroupGeneric",
             "SizeOfMatrixGroup", "NiceMonomorphism", "IsomorphismPermGroup"):
    print(f"  {name}: {bool(libgap.eval(f'IsBound({name})'))}", flush=True)

print("--- control: GL_2(5) ---", flush=True)
control = MatrixGroup([matrix(GF(5), [[0, 1], [1, 0]]), matrix(GF(5), [[1, 1], [0, 1]])])
phase("control IsomorphismPermGroup", lambda: libgap.Size(libgap.IsomorphismPermGroup(control.gap()).Range()))
phase("control IsomorphismTypeOfMatrixGroup",
      lambda: libgap.IsomorphismTypeOfMatrixGroup(control.gap()))

print("--- the real case: O(L) mod 5, degree 10 ---", flush=True)
case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup2")
lattice_ = lattice(case["gram"])
computed_rows = [rows_of(isometry) for isometry in lattice_.O().group_generators()]
group = MatrixGroup([matrix(GF(5), rows) for rows in computed_rows])
print(f"  {len(computed_rows)} generators built", flush=True)
phase("p=5 IsomorphismTypeOfMatrixGroup",
      lambda: libgap.IsomorphismTypeOfMatrixGroup(group.gap()))
