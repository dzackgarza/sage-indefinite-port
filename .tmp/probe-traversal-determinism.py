"""Is `O(L).group_generators()` for LORENTZ_OrthogonalGroup0 the same in every process?

Two runs of the same acceptance case on the same pins returned 16 and 17 generators and 1 and 2 traversal
expansions. This counts the decomposition itself, per process, under a chosen PYTHONHASHSEED.
"""
import os

import sage_indefinite_port.indefinite.lorentzian_cells as lc
from tests.acceptance.consumer import lattice, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

SEED = os.environ.get("PYTHONHASHSEED", "unset")
print(f"[determinism] PYTHONHASHSEED={SEED} pid={os.getpid()}", flush=True)

COMPLEXES = []
original = lc._cached_perfect_complex


def capturing(*args, **kwargs):
    result = original(*args, **kwargs)
    if result not in COMPLEXES:
        COMPLEXES.append(result)
    return result


lc._cached_perfect_complex = capturing

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup0")
lattice_ = lattice(case["gram"])
generators = lattice_.O().group_generators()
print(f"[determinism] generators: {len(generators)}", flush=True)

for position, complex_ in enumerate(COMPLEXES):
    cells = complex_.quotient_cells()
    adjacencies = complex_.adjacencies()
    transportered = sum(1 for adjacency in adjacencies if adjacency.transporter is not None)
    sizes = tuple(len(cell.rows) for cell in cells)
    print(
        f"[determinism] complex {position}: cells={len(cells)} sizes={sizes} "
        f"adjacencies={len(adjacencies)} with_transporter={transportered} "
        f"hash_seed_seen={SEED}",
        flush=True,
    )

rows = sorted(
    tuple(tuple(int(entry) for entry in row) for row in rows_of(generator))
    for generator in generators
)
import hashlib
digest = hashlib.sha256(repr(rows).encode()).hexdigest()[:16]
print(f"[determinism] generator-set digest: {digest}", flush=True)

# Persist the rows so two seeds' generating sets can be compared as groups later.
out = f"/work/.tmp/generators-seed{os.environ.get('LORENTZ_SEED', 'x')}.json"
with open(out, "w") as handle:
    import json
    json.dump({"seed": os.environ.get("LORENTZ_SEED"), "count": len(rows), "rows": rows,
               "cells": [len(c.rows) for c in COMPLEXES[0].quotient_cells()] if COMPLEXES else []}, handle)
print(f"[determinism] wrote {out}", flush=True)
