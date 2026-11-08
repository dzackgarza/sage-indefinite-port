"""Can the LORENTZ mod-p image be certified structurally, without a BSGS on F_p^10?

GAP's classical group object carries its order by construction:
  libgap.Size(libgap.GO(+1, 10, 5)) -> 54461311079175000000000000000000   in 0.0 s
(after the closed form |O^+_{2n}(q)| = 2 q^{n(n-1)} (q^n - 1) prod (q^{2i} - 1), n = 5).

If membership of generators in GO is also cheap, then "the image is the full orthogonal group" becomes a
bounded certificate and the harness's order computation is not the only route. This probe measures exactly
that, cheapest phases first, so a hang still leaves the earlier numbers.
"""
import time

from sage.all import GF, MatrixGroup, matrix, libgap

from tests.acceptance.consumer import lattice, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup2")
lattice_ = lattice(case["gram"])
started = time.monotonic()
computed_rows = [rows_of(isometry) for isometry in lattice_.O().group_generators()]
print(f"computed {len(computed_rows)} generators in {time.monotonic() - started:.1f}s", flush=True)


def phase(label, fn):
    started = time.monotonic()
    try:
        value = fn()
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:150]} "
              f"after {time.monotonic() - started:.1f}s", flush=True)
        return None
    print(f"[{label}] {value} in {time.monotonic() - started:.2f}s", flush=True)
    return value


print("--- orders of the classical groups, both primes, no lattice involved ---", flush=True)
for prime in (5, 7):
    phase(f"p={prime} Size(GO(+1,10,q))", lambda q=prime: libgap.Size(libgap.GO(1, 10, q)))

print("--- membership against the classical group, p=5 ---", flush=True)
prime = 5
group = MatrixGroup([matrix(GF(prime), rows) for rows in computed_rows])
go = libgap.GO(1, 10, prime)
phase("p=5 IsSubgroup(GO, <computed>)", lambda: libgap.IsSubgroup(go, group.gap()))
phase("p=5 IsSubgroup(<computed>, GO)", lambda: libgap.IsSubgroup(group.gap(), go))
phase("p=5 every generator in GO",
      lambda: all(bool(libgap.IsContained(go, libgap(g))) for g in group.gens()[:2]))
phase("p=5 Size(<computed>)", lambda: libgap.Size(group.gap()))
