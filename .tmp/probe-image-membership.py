"""Do the computed and recorded generators span the SAME subgroup mod p, at p=5 and p=7?

The acceptance assertion compares two independently computed group ORDERS, which at p=5 and p=7 means
enumerating F_p^10. Containment is the question the criterion's wording ("the same finite images") actually
asks, and it is a membership problem.

Sage's `ParentLibGAP.is_subgroup()` means "was constructed as a subgroup of a bigger group" and is NOT a
containment test (sage/groups/libgap_wrapper.pyx:154), so two candidate mechanisms are exercised on a tiny
fixture with a known answer and only the one that passes is used on the real question.
"""
import json
import os
import time

from sage.all import GF, ZZ, MatrixGroup, libgap, matrix

from tests.acceptance.consumer import lattice, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases


def group(rows, prime):
    return MatrixGroup([matrix(GF(prime), row) for row in rows])


def contains(outer, inner):
    """Whether inner <= outer, by whichever mechanism the fixture validated."""
    if MECHANISM == "membership":
        return all(generator in outer for generator in inner.gens())
    return bool(libgap.IsSubgroup(outer.gap(), inner.gap()))


MECHANISM = None


def fixture_selftest() -> str:
    a = group([[[2, 0], [0, 1]]], 5)
    bigger = group([[[2, 0], [0, 1]], [[1, 0], [0, 2]]], 5)
    outside = group([[[2, 0], [0, 1]], [[0, 1], [1, 0]]], 5)
    for name, fn in (
        ("membership", lambda outer, inner: all(g in outer for g in inner.gens())),
        ("libgap.IsSubgroup", lambda outer, inner: bool(libgap.IsSubgroup(outer.gap(), inner.gap()))),
    ):
        try:
            got = (fn(bigger, a), fn(a, bigger), fn(bigger, outside))
        except Exception as error:
            print(f"[fixture] {name}: FAILED {type(error).__name__}: {error}", flush=True)
            continue
        ok = got == (True, False, False)
        print(f"[fixture] {name}: {got} (expect (True, False, False)) {'VALIDATED' if ok else 'WRONG'}", flush=True)
        if ok:
            return name
    raise RuntimeError("no containment mechanism passed the fixture")


MECHANISM = fixture_selftest()
print(f"[fixture] using {MECHANISM}", flush=True)

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == os.environ.get("LORENTZ_CASE", "LORENTZ_OrthogonalGroup0"))
gram = case["gram"]

lattice_ = lattice(gram)
started = time.monotonic()
computed = lattice_.O().group_generators()
computed_rows = [rows_of(isometry) for isometry in computed]
print(f"[rows] {len(computed_rows)} computed generators in {time.monotonic() - started:.1f}s "
      f"| {len(case['generators'])} recorded", flush=True)

report = {}
for prime in (5, 7):
    started = time.monotonic()
    computed_group = group(computed_rows, prime)
    recorded_group = group(case["generators"], prime)
    built = time.monotonic() - started
    started = time.monotonic()
    forward = bool(contains(recorded_group, computed_group))
    first = time.monotonic() - started
    started = time.monotonic()
    backward = bool(contains(computed_group, recorded_group))
    second = time.monotonic() - started
    report[prime] = {
        "build_s": round(built, 2),
        "computed_le_recorded": forward,
        "recorded_le_computed": backward,
        "equal_subgroups": forward and backward,
        "membership_s": [round(first, 2), round(second, 2)],
        "generator_counts": [len(computed_rows), len(case["generators"])],
    }
    print(f"[membership p={prime}] {report[prime]}", flush=True)

print("[result] " + json.dumps(report), flush=True)
