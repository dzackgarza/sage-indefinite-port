"""One pointed-E10 Isom unit from the q>=6 plateau, profiled.

Specimen is the exact one the acceptance node builds: root r_-1 and the recorded orbit
representative at index 3 (q = 6), which measured 56.12 s uninstrumented.

Attribution uses cProfile because targeted wrappers have already missed a hot path once
in this session (the OrbitCoverModel guess). cProfile's overhead distorts absolute
times but not relative attribution; the uninstrumented figure is recorded alongside.
faulthandler dumps every 5 s give the in-flight picture for the slow section.
"""
import cProfile
import faulthandler
import pstats
import time

from sage.all import ZZ, identity_matrix, matrix, vector

from tests.acceptance.consumer import (
    element,
    isomorphisms,
    lattice,
    pointed,
    pointed_category,
    underlying_isometry,
)
from tests.fixtures.oracle_fixtures import load_e10_fundamental_domain

_E10 = load_e10_fundamental_domain()
gram = _E10["gram"]
lattice_ = lattice(gram)


def _reflection_rows(gram, root):
    form = matrix(ZZ, gram)
    r = vector(ZZ, root)
    return [list(bv + (bv * form * r) * r) for bv in identity_matrix(ZZ, len(root)).rows()]


root = next(r for r in _E10["simple_roots"] if str(r["label"]) == "-1")
reflection = matrix(ZZ, _reflection_rows(gram, root["vector"]))
INDEX = 3
representative = _E10["orbit_representatives"][INDEX]
source = vector(ZZ, representative["vector"])
target = source * reflection
source_vector = element(lattice_, source)
norm = int(source_vector.q())
assert norm == 6, f"specimen changed: representative {INDEX} has q={norm}, expected 6"
print(f"specimen: root r_-1, representative index {INDEX}, vector {representative['vector']}, q={norm}")
print("uninstrumented duration for this exact specimen was 56.12 s")

faulthandler.dump_traceback_later(30, repeat=True)
profiler = cProfile.Profile()
started = time.monotonic()
profiler.enable()
isometry = underlying_isometry(
    isomorphisms(
        pointed_category(source_vector.q()),
        pointed(lattice_, source_vector),
        pointed(lattice_, element(lattice_, target)),
    ).an_element()
)
profiler.disable()
elapsed = time.monotonic() - started
faulthandler.cancel_dump_traceback_later()
print(f"profiled unit elapsed: {elapsed:.2f}s")

profiler.dump_stats("/logs/pointed-q6.prof")
stats = pstats.Stats(profiler)
print("\n=== top 22 by CUMULATIVE ===")
stats.sort_stats("cumulative").print_stats(22)
print("\n=== top 16 by TOTTIME (this function only) ===")
stats.sort_stats("tottime").print_stats(16)

print("\n=== eager-representation check: repr/str entries in the profile ===")
hits = 0
for (filename, line, name), (_cc, nc, tt, ct, callers) in stats.stats.items():
    if "repr" in name or "_repr_" in name or name in ("__str__", "format"):
        print(f"   {name:28s} {filename.split('/')[-1]}:{line} calls={nc} cum={ct:.2f}s")
        hits += 1
print(f"   ({hits} repr/str entries)")
