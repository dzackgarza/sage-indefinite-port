"""Do two hash-seed-dependent generating sets span the SAME group?

The generator set for LORENTZ_OrthogonalGroup0 varies with PYTHONHASHSEED (16 vs 17 generators). The
acceptance question is whether the group is affected, or only the choice of generators. Compared at the
two primes where an image order is affordable: p=2 and p=3.
"""
import glob
import json

from sage.all import GF, MatrixGroup, matrix

files = sorted(glob.glob("/work/.tmp/generators-seed*.json"))
payloads = [json.load(open(path)) for path in files]
for payload in payloads:
    print(f"[compare] seed {payload['seed']}: {payload['count']} generators, cells {payload['cells']}", flush=True)

for prime in (2, 3):
    groups = [
        MatrixGroup([matrix(GF(prime), [list(row) for row in rows]) for rows in payload["rows"]])
        for payload in payloads
    ]
    for left in range(len(groups)):
        for right in range(left + 1, len(groups)):
            forward = all(g in groups[right] for g in groups[left].gens())
            backward = all(g in groups[left] for g in groups[right].gens())
            print(
                f"[compare] p={prime} seed {payloads[left]['seed']} vs {payloads[right]['seed']}: "
                f"left<=right {forward}, right<=left {backward}, equal {forward and backward}",
                flush=True,
            )
