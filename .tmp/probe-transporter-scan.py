"""Find a pointed pair whose rational transport is NON-integral, so the correction branch runs.

`pointed_isometry` (recursive.py:795-805) takes that branch when
`rational_transport.is_integral_on(source_lattice, target_lattice)` is False, and only then calls
`self._lattice_action(target).transporter(image)` -- the provider's finite-representation transporter that the
BOUND note says fails on a represented finite torsion module. This scans small lattices for such a pair,
bounded, and reports exactly which pair enters the branch and what the transporter does.
"""
import itertools
import time

from sage_extensions.categories.rational_integral_stabilizers import IntegralStructureAction
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import element, lattice

events = []
original = IntegralStructureAction.transporter

def counted(self, *args, **kwargs):
    try:
        result = original(self, *args, **kwargs)
        events.append(("ok", result is not None))
        return result
    except Exception as error:
        events.append(("raised", type(error).__name__, str(error)[:140]))
        raise

IntegralStructureAction.transporter = counted
RECURSION = EichlerRecursion()

GRAMMARS = {
    "U + <-2>": [[0, 1, 0], [1, 0, 0], [0, 0, -2]],
    "U + <2>": [[0, 1, 0], [1, 0, 0], [0, 0, 2]],
    "<2> + <2> + <-2>": [[2, 0, 0], [0, 2, 0], [0, 0, -2]],
    "U(2) + <2>": [[0, 2, 0], [2, 0, 0], [0, 0, 2]],
}
BOX = 2

started_all = time.monotonic()
for name, gram in GRAMMARS.items():
    L = lattice(gram)
    rank = int(L.rank())
    vectors = []
    for coordinates in itertools.product(range(-BOX, BOX + 1), repeat=rank):
        if not any(coordinates):
            continue
        vector = element(L, list(coordinates))
        if int(vector.content()) != 1:
            continue
        vectors.append(vector)
    for source, target in itertools.permutations(vectors, 2):
        if int(source.q()) != int(target.q()):
            continue
        events.clear()
        try:
            witness = RECURSION.pointed_isometry(source, target)
        except Exception as error:
            if events:
                print(f"[{name}] {tuple(int(c) for c in source.to_vector())} -> "
                      f"{tuple(int(c) for c in target.to_vector())} RAISED {type(error).__name__}: {str(error)[:110]}",
                      flush=True)
                print(f"    transporter events: {events}", flush=True)
                raise SystemExit(0)
            continue
        if events:
            print(f"[{name}] q={int(source.q())} {tuple(int(c) for c in source.to_vector())} -> "
                  f"{tuple(int(c) for c in target.to_vector())} | witness {witness is not None} | "
                  f"transporter events {events}", flush=True)
            raise SystemExit(0)
    if time.monotonic() - started_all > 240:
        print("time budget reached without hitting the branch", flush=True)
        break
print(f"scan finished; no transporter branch entered "
      f"({time.monotonic() - started_all:.1f}s)", flush=True)
