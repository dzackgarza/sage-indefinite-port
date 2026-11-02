"""LORENTZ image-order obligation: is the p=5 computation capped, or is the assertion wrong?

The acceptance assertion compares, for each prime in comparison_primes(gram) = {2,3,5,7},

    reduction_image_order(rows, p) = MatrixGroup([matrix(GF(p), rows)]).order()

for the computed and the recorded generating set. On LORENTZ_OrthogonalGroup0 the p=2 and p=3 values
agreed on both sides and p=5 did not finish inside the cap, so the record is CAPPED, not a demonstrated
wrong value. This probe decides the smallest unresolved case (LORENTZ_OrthogonalGroup2, 14 recorded
generators) and asks whether the maintained provider semantic

    O(L) --base change--> O(L (x) F_p)   (a finite group)

reproduces the same numbers where the harness is affordable, and whether it terminates at p=5.

Nothing is weakened: the harness function is still called and its value compared, and comparison primes are
unchanged.
"""
import os
import time

from sage.all import GF, ZZ
from sage_extensions.all import PrimeFields
from sage_extensions.all import ZZ as OwnedZZ

from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import (
    comparison_primes,
    lattice,
    reduction_image_order,
    rows_of,
)
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

CASE = os.environ.get("LORENTZ_CASE", "LORENTZ_OrthogonalGroup2")
case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == CASE)
gram = case["gram"]
print(f"=== {CASE}: rank {len(gram)}, determinant {int(__import__('sage').all.matrix(ZZ, gram).det())}, "
      f"{len(case['generators'])} recorded generators", flush=True)
print(f"    comparison primes: {comparison_primes(gram)}", flush=True)

lattice_ = lattice(gram)
started = time.monotonic()
computed = lattice_.O().group_generators()
computed_rows = [rows_of(isometry) for isometry in computed]
print(f"    computed {len(computed_rows)} generators in {time.monotonic() - started:.1f}s", flush=True)


def ring_map_to(prime):
    """ZZ -> F_p as a ring map the provider's base_change can consume (provider rings, own session)."""
    field = PrimeFields()(prime)
    try:
        return OwnedZZ.hom(field)
    except Exception as error:
        print(f"    OwnedZZ.hom(PrimeFields()({prime})) unavailable: {type(error).__name__}: {error}",
              flush=True)
        raise


def provider_image_order(prime):
    """Order of the image of O(L) in O(L (x) F_p), through the provider's own base-change semantic."""
    change = lattice_.orthogonal_group_base_change(ring_map_to(prime))
    target = lattice_.base_change(ring_map_to(prime)).Aut()
    images = tuple(change(isometry) for isometry in computed)
    return int(target.subgroup_on(images).cardinality())


for prime in (2, 3):
    started = time.monotonic()
    harness = reduction_image_order(computed_rows, prime)
    harness_seconds = time.monotonic() - started
    started = time.monotonic()
    try:
        provider = provider_image_order(prime)
    except Exception as error:
        print(f"[p={prime}] provider route failed: {type(error).__name__}: {error}", flush=True)
        provider = None
    provider_seconds = time.monotonic() - started
    status = "agree" if provider == harness else "DISAGREE"
    print(f"[p={prime}] harness {harness} ({harness_seconds:.2f}s) | provider {provider} "
          f"({provider_seconds:.2f}s) -> {status}", flush=True)
    if prime == 2 or prime == 3:
        print(f"[p={prime}] recorded-side harness "
              f"{reduction_image_order(case['generators'], prime)}", flush=True)

print(f"--- now p=5 on {CASE} (unchanged prime, no cap change here beyond the harness guard) ---", flush=True)
started = time.monotonic()
try:
    value = provider_image_order(5)
    print(f"[p=5] provider route returned {value} in {time.monotonic() - started:.1f}s", flush=True)
except Exception as error:
    print(f"[p=5] provider route raised {type(error).__name__}: {error} "
          f"after {time.monotonic() - started:.1f}s", flush=True)
