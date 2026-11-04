"""Find the supported subgroup/order surface of O(L (x) F_p), then ask it for the p=5 image order."""
import time
import traceback

from sage_extensions.all import PrimeFields
from tests.acceptance.consumer import lattice, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup2")
lattice_ = lattice(case["gram"])
computed = lattice_.O().group_generators()
print(f"computed {len(computed)} generators", flush=True)

ring_map = PrimeFields()(2).quotient_map()
change = lattice_.orthogonal_group_base_change(ring_map)
target = lattice_.base_change(ring_map).Aut()
print("type(target):", type(target).__module__ + "." + type(target).__name__, flush=True)
for name in ("subgroup", "subgroup_on", "cardinality", "order", "group_generators", "Aut", "domain"):
    print(f"  has {name}: {hasattr(target, name)}", flush=True)
images = tuple(target(change(g)) for g in computed)
print("first image type:", type(images[0]).__module__ + "." + type(images[0]).__name__, flush=True)
print("image parent is target:", all(image.parent() is target for image in images), flush=True)

for label, fn in (
    ("target.subgroup(images)", lambda: target.subgroup(images)),
):
    started = time.monotonic()
    try:
        value = fn()
        print(f"[{label}] {value} in {time.monotonic() - started:.2f}s", flush=True)
        subgroup = value
        for order_label, order_fn in (("cardinality", lambda: subgroup.cardinality()),
                                      ("order", lambda: subgroup.order())):
            if not hasattr(subgroup, order_label):
                print(f"  has {order_label}: False", flush=True)
                continue
            started = time.monotonic()
            try:
                print(f"  [{order_label}] {order_fn()} in {time.monotonic() - started:.2f}s", flush=True)
            except Exception:
                print(f"  [{order_label}] raised after {time.monotonic() - started:.2f}s", flush=True)
                traceback.print_exc()
    except Exception:
        print(f"[{label}] raised after {time.monotonic() - started:.2f}s", flush=True)
        traceback.print_exc()
