"""The maintained provider route to the mod-p image order, invoked through the canonical quotient map.

Earlier I called `OwnedZZ.hom(PrimeFields()(5))` and got "the identity morphism exists only on Mor(R,R)". That
was the wrong invocation, not evidence that the route is unsupported: a quotient ring Q = R/I carries the
retained canonical map `Q.quotient_map() : R -> Q`, and that is what `orthogonal_group_base_change` should be
given.

The provider route is validated at p=2 and p=3 against the harness numbers already established (those are an
API check, not a conclusion about p=5), and only then applied at p=5.
"""
import time

from sage.all import GF, ZZ, matrix

from sage_extensions.all import PrimeFields

from tests.acceptance.consumer import lattice, reduction_image_order, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup2")
lattice_ = lattice(case["gram"])
started = time.monotonic()
computed = lattice_.O().group_generators()
computed_rows = [rows_of(isometry) for isometry in computed]
print(f"computed {len(computed_rows)} generators in {time.monotonic() - started:.1f}s", flush=True)


def phase(label, fn):
    started = time.monotonic()
    try:
        value = fn()
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:220]} "
              f"after {time.monotonic() - started:.1f}s", flush=True)
        return None
    print(f"[{label}] {value} in {time.monotonic() - started:.2f}s", flush=True)
    return value


def quotient_map(prime):
    field = PrimeFields()(prime)
    return field.quotient_map()


for prime in (2, 3):
    field = PrimeFields()(prime)
    ring_map = quotient_map(prime)
    print(f"--- p={prime}: map {ring_map.domain()} -> {ring_map.codomain()}; "
          f"lattice base ring {lattice_.base_ring()}; same object "
          f"{ring_map.domain() is lattice_.base_ring()}", flush=True)
    change = phase(f"p={prime} orthogonal_group_base_change", lambda m=ring_map: lattice_.orthogonal_group_base_change(m))
    if change is None:
        continue
    target = phase(f"p={prime} target.Aut()", lambda m=ring_map: lattice_.base_change(m).Aut())
    if target is None:
        continue
    images = phase(f"p={prime} image of generators", lambda c=change, t=target: tuple(t(c(g)) for g in computed))
    if images is None:
        continue
    subgroup = phase(f"p={prime} subgroup_on", lambda t=target, i=images: t.subgroup_on(i))
    if subgroup is None:
        continue
    provider = phase(f"p={prime} cardinality", lambda s=subgroup: int(s.cardinality()))
    harness = reduction_image_order(computed_rows, prime)
    print(f"[p={prime}] provider {provider} | harness {harness} | agree {provider == harness}", flush=True)

print("--- p=5: the obligation itself ---", flush=True)
ring_map = quotient_map(5)
print(f"    map {ring_map.domain()} -> {ring_map.codomain()}", flush=True)
change = phase("p=5 orthogonal_group_base_change", lambda: lattice_.orthogonal_group_base_change(ring_map))
if change is not None:
    target = phase("p=5 target.Aut()", lambda: lattice_.base_change(ring_map).Aut())
    if target is not None:
        images = phase("p=5 image of generators", lambda: tuple(target(change(g)) for g in computed))
        if images is not None:
            subgroup = phase("p=5 subgroup_on", lambda: target.subgroup_on(images))
            if subgroup is not None:
                value = phase("p=5 cardinality", lambda: int(subgroup.cardinality()))
                print(f"[p=5] provider image order = {value}", flush=True)
