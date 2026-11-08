"""Verify the admission fix at 42d7361c WITHOUT re-running the engine BSGS at p=5.

p=2 is cheap (1024 points) so the provider route is checked end-to-end against the harness number there.
At p=5 only admission is checked: finiteness is asked (cheap), and the placement it should induce is
inspected. cardinality() is deliberately NOT called at p=5 -- the engine cost is unchanged by this fix.
"""
import time

from sage_extensions.all import PrimeFields
from sage_extensions.categories.group.groups import OwnedFiniteGroups
from tests.acceptance.consumer import lattice, reduction_image_order, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup2")
lattice_ = lattice(case["gram"])
computed = lattice_.O().group_generators()
computed_rows = [rows_of(isometry) for isometry in computed]
print(f"computed {len(computed)} generators", flush=True)

for prime in (2, 5):
    ring_map = PrimeFields()(prime).quotient_map()
    target = lattice_.base_change(ring_map).Aut()
    images = tuple(target(lattice_.orthogonal_group_base_change(ring_map)(g)) for g in computed)
    subgroup = target.subgroup(images)
    started = time.monotonic()
    finite = subgroup.is_finite()
    asked = time.monotonic() - started
    placed = subgroup in OwnedFiniteGroups()
    print(f"[p={prime}] is_finite() {finite} in {asked:.2f}s | in OwnedFiniteGroups() now {placed}",
          flush=True)
    if prime == 2:
        started = time.monotonic()
        value = int(subgroup.cardinality())
        print(f"[p={prime}] cardinality {value} in {time.monotonic() - started:.2f}s | "
              f"harness {reduction_image_order(computed_rows, 2)} | "
              f"agree {value == reduction_image_order(computed_rows, 2)}", flush=True)
    else:
        print(f"[p={prime}] cardinality NOT called: admission fixed, engine cost unchanged", flush=True)
