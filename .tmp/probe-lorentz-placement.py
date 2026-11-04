"""Is O(L (x) F_p) placed finite, and is that the missing decision?"""
from sage_extensions.all import PrimeFields
from sage_extensions.categories.group.groups import OwnedFiniteGroups
from tests.acceptance.consumer import lattice, rows_of
from tests.fixtures.oracle_fixtures import load_lorentzian_stabilizer_cases

case = next(c for c in load_lorentzian_stabilizer_cases() if c["id"] == "LORENTZ_OrthogonalGroup2")
lattice_ = lattice(case["gram"])
computed = lattice_.O().group_generators()

for prime in (2, 5):
    ring_map = PrimeFields()(prime).quotient_map()
    changed = lattice_.base_change(ring_map)
    target = changed.Aut()
    images = tuple(target(lattice_.orthogonal_group_base_change(ring_map)(g)) for g in computed)
    subgroup = target.subgroup(images)
    print(f"--- p={prime}: {changed}", flush=True)
    print(f"    base ring finite: {changed.base_ring().is_finite()}", flush=True)
    print(f"    O(V) in OwnedFiniteGroups(): {target in OwnedFiniteGroups()}", flush=True)
    print(f"    O(V).is_finite(): {target.is_finite()}", flush=True)
    print(f"    subgroup in OwnedFiniteGroups(): {subgroup in OwnedFiniteGroups()}", flush=True)
    print(f"    subgroup.is_finite(): {subgroup.is_finite()}", flush=True)
    print(f"    subgroup.has_selected_group_resolution(): {subgroup.has_selected_group_resolution()}",
          flush=True)
    print(f"    O(V).has_selected_group_resolution(): {target.has_selected_group_resolution()}", flush=True)
