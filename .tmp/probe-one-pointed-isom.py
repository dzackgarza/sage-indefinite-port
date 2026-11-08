"""Time a single pointed E10 Isom: the unit that the r_* node repeats 61 times."""
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
print("root r_-1 vector:", root["vector"], "reflection shape:", reflection.nrows(), "x", reflection.ncols())

representative = _E10["orbit_representatives"][0]
source = vector(ZZ, representative["vector"])
target = source * reflection
source_vector = element(lattice_, source)
print("first representative:", representative["vector"], "-> q =", source_vector.q())
print("orbit representatives per node:", len(_E10["orbit_representatives"]))

started = time.monotonic()
isometry = underlying_isometry(
    isomorphisms(
        pointed_category(source_vector.q()),
        pointed(lattice_, source_vector),
        pointed(lattice_, element(lattice_, target)),
    ).an_element()
)
elapsed = time.monotonic() - started
print(f"ONE pointed Isom: {elapsed:.2f}s")
print(f"61 of these in one r_* node: ~{elapsed * 61:.0f}s")

print("\n=== full node loop: all 61 representatives for r_-1 ===")
total = 0.0
slowest = []
for index, representative in enumerate(_E10["orbit_representatives"]):
    src = vector(ZZ, representative["vector"])
    tgt = src * reflection
    sv = element(lattice_, src)
    t0 = time.monotonic()
    iso = underlying_isometry(
        isomorphisms(
            pointed_category(sv.q()),
            pointed(lattice_, sv),
            pointed(lattice_, element(lattice_, tgt)),
        ).an_element()
    )
    dt = time.monotonic() - t0
    total += dt
    slowest.append((dt, index, representative["vector"], int(sv.q())))
    if dt > 20:
        print(f"  rep {index:2d} q={int(sv.q()):4d} {dt:7.2f}s  <-- slow", flush=True)
slowest.sort(reverse=True)
print(f"TOTAL for 61: {total:.2f}s")
print("slowest five:", [(round(d, 2), i, int(q)) for d, i, v, q in slowest[:5]])
