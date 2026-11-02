"""Build the G-set of O(L) acting on a lattice (degenerate specimen) through GObjects."""
from sage_extensions.all import GObjects, Modules, ZZ as OwnedZZ
from tests.acceptance.consumer import element, lattice

L = lattice([[2, 0], [0, 0]])
v = element(L, [1, 0])
radical_vector = element(L, [0, 1])
group = L.O()
print("group object:", type(group).__name__, flush=True)

candidates = [
    ("GObjects(L.O(), Modules(ZZ))", lambda: GObjects(group, Modules(OwnedZZ))),
    ("GObjects(L.O(), Modules(ZZ)(L))", lambda: GObjects(group, Modules(OwnedZZ)(L))),
]
for label, fn in candidates:
    try:
        action = fn()
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:120]}", flush=True)
        continue
    print(f"[{label}] -> {type(action).__name__} | methods: "
          f"{[m for m in ('orbits', 'orbit_of', 'stabilizer', 'in_same_orbit', 'action_on') if hasattr(action, m)]}",
          flush=True)
    try:
        orbits = action.orbits()
        print(f"    orbits() -> {type(orbits).__name__}", flush=True)
        orbit = orbits.orbit_of(v)
        points = orbit.points()
        print(f"    orbit_of((1,0)).points() -> {type(points).__name__}", flush=True)
        print(f"    v in points: {v in points} | (3,0) in points: {element(L, [3, 0]) in points}",
              flush=True)
        print(f"    is_finite: {points.is_finite() if hasattr(points, 'is_finite') else 'n/a'}", flush=True)
    except Exception as error:
        print(f"    orbits raised {type(error).__name__}: {str(error)[:150]}", flush=True)
