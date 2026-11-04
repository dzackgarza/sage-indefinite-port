"""Attribute the cost of one A_GROUPS Lorentzian node, without changing what it computes.

Runs the exact acceptance path -- `lattice(case["gram"]).O().group_generators()`, then the
original `assert_same_finite_images` -- and wraps (in this process only) the perfect-complex
traversal and the exact-transporter search, so a wall cap still leaves a named locus and a
split between enumeration and transporter work. Nothing in src/ is edited.
"""
import atexit
import faulthandler
import os
import sys
import time
import traceback

import sage_indefinite_port.indefinite.recursive as recursive
from sage_indefinite_port.indefinite import lorentzian_cells as lc
from tests.acceptance.consumer import MARKED_AUTOMORPHISMS, assert_same_finite_images, lattice, require
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

STACK_EVERY = int(os.environ.get("STACK_EVERY", "120"))
faulthandler.dump_traceback_later(STACK_EVERY, repeat=True)

TOTALS = {}
STAGE_OUT = os.environ.get("STAGE_OUT", "/work/.tmp/lorentzian-stages.json")
_LAST_WRITE = [0.0]


def snapshot(force=False):
    """Persist the totals, so a hard kill still leaves the attribution on disk."""
    import json

    now = time.monotonic()
    if not force and now - _LAST_WRITE[0] < 1.0:
        return
    _LAST_WRITE[0] = now
    try:
        with open(STAGE_OUT, "w") as handle:
            json.dump({"totals": TOTALS, "elapsed": now}, handle, indent=1)
    except OSError:
        pass


def bump(key, seconds, rows=None, **extra):
    entry = TOTALS.setdefault(key, {"n": 0, "total": 0.0, "max": 0.0, "rows_max": 0, "last": 0.0})
    entry["n"] += 1
    entry["total"] += seconds
    entry["last"] = seconds
    if seconds > entry["max"]:
        entry["max"] = seconds
    if rows is not None and rows > entry["rows_max"]:
        entry["rows_max"] = rows


def wrap(cls, name, label, rows_of_args=None):
    original = getattr(cls, name)

    def wrapper(self, *args, **kwargs):
        started = time.monotonic()
        rows = rows_of_args(args, kwargs) if rows_of_args else None
        try:
            return original(self, *args, **kwargs)
        finally:
            elapsed = time.monotonic() - started
            bump(label, elapsed, rows)
            snapshot()
            entry = TOTALS[label]
            if label == "expand_traversal_once":
                # Print early and often, then thin out: a traversal can have many cells,
                # and every line must still be cheap relative to the work it describes.
                if entry["n"] <= 200 or entry["n"] % 25 == 0:
                    print(
                        f"[stage] {label} #{entry['n']} {elapsed:8.2f}s rows={rows} "
                        f"cum={entry['total']:8.2f}s",
                        flush=True,
                    )
            elif elapsed >= 5.0:
                print(
                    f"[stage] {label} #{entry['n']} {elapsed:8.2f}s rows={rows} "
                    f"all_images={kwargs.get('all_images')} "
                    f"generators_only={kwargs.get('generators_only')}",
                    flush=True,
                )

    setattr(cls, name, wrapper)


def _source_rows(args, kwargs):
    try:
        return len(args[0].rows)
    except Exception:
        return None


import tests.acceptance.consumer as consumer


def own_rss_mib():
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024
    except OSError:
        pass
    return 0.0


def wrap_assertion(name):
    """Time the assertion's finite-image comparisons, including the prime they reduce by."""
    original = getattr(consumer, name)

    def wrapper(*args, **kwargs):
        tag = name if name != "reduction_image_order" else f"reduction_image_order[p={args[-1]}]"
        before = own_rss_mib()
        print(f"[assert] {tag} ... starting at rss {before:.0f} MiB", flush=True)
        started = time.monotonic()
        result = original(*args, **kwargs)
        elapsed = time.monotonic() - started
        bump(tag, elapsed)
        snapshot()
        print(
            f"[assert] {tag} {elapsed:8.2f}s rss {before:.0f}->{own_rss_mib():.0f} MiB -> {result}",
            flush=True,
        )
        return result

    setattr(consumer, name, wrapper)


wrap_assertion("discriminant_image_order")
wrap_assertion("reduction_image_order")

wrap(lc.LorentzianPerfectComplex, "_expand_traversal_once", "expand_traversal_once", _source_rows)
wrap(lc.LorentzianPerfectLocalBackend, "cell_transporter", "cell_transporter", _source_rows)
wrap(lc.LorentzianPerfectLocalBackend, "cell_stabilizer", "cell_stabilizer", _source_rows)


def report():
    snapshot(force=True)
    print("\n================ STAGE TOTALS ================", flush=True)
    for label in sorted(TOTALS):
        entry = TOTALS[label]
        print(
            f"[stage] {label:26s} n={entry['n']:6d} total={entry['total']:9.2f}s "
            f"max={entry['max']:8.2f}s rows_max={entry['rows_max']}",
            flush=True,
        )


atexit.register(report)

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == os.environ.get("LORENTZ_CASE", "LORENTZ_OrthogonalGroup0"))
require(MARKED_AUTOMORPHISMS)
gram = case["gram"]
print(f"=== case {case['id']} | rank {len(gram)} | recorded generators {len(case['generators'])}", flush=True)

# The reduction is what the public path actually traverses; time it and name its target.
original_reduction = recursive.indefinite_reduced_presentation


def timed_reduction(target):
    started = time.monotonic()
    result = original_reduction(target)
    elapsed = time.monotonic() - started
    bump("indefinite_reduced_presentation", elapsed)
    if elapsed >= 1.0:
        print(
            f"[stage] indefinite_reduced_presentation {elapsed:8.2f}s "
            f"-> rank {result.reduced.rank()} signature {result.reduced.signature_pair()}",
            flush=True,
        )
    return result


recursive.indefinite_reduced_presentation = timed_reduction

lattice_ = lattice(gram)
started = time.monotonic()
computed = lattice_.O().group_generators()
print(f"\n=== group_generators returned {len(computed)} generators in {time.monotonic() - started:.1f}s", flush=True)

started = time.monotonic()
assert_same_finite_images(gram, computed, case["generators"])
print(f"=== ORIGINAL ASSERTION PASSED in {time.monotonic() - started:.1f}s", flush=True)
