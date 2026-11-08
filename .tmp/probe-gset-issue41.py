"""The supported semantic action for O(L) on a degenerate lattice, through GObjects(...).on_set(...)."""
from sage_extensions.all import GObjects, Sets
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import element, lattice

RECURSION = EichlerRecursion()
L = lattice([[2, 0], [0, 0]])
v = element(L, [1, 0])
radical = element(L, [0, 1])
other = element(L, [3, 0])


def orbit_relation(left, right):
    """O(L)-equivalence of two vectors, decided by the port's pointed transporter."""
    return RECURSION.pointed_isometry(left, right) is not None


action = GObjects(L.O(), Sets()).on_set(L, lambda g, x: g(x), orbit_relation=orbit_relation)
print("action:", type(action).__name__, flush=True)
print("  surface:", [m for m in ("orbits", "orbit_of", "stabilizer", "in_same_orbit") if hasattr(action, m)],
      flush=True)

for label, point in (("(1,0)", v), ("(0,1) radical, c=0", radical)):
    try:
        orbit = action.orbits().orbit_of(point)
        points = orbit.points()
        print(f"[{label}] orbit {type(orbit).__name__} | points {type(points).__name__}", flush=True)
        print(f"    point in points: {point in points}", flush=True)
        print(f"    (3,0) in points: {other in points}", flush=True)
        finite = points.is_finite() if hasattr(points, "is_finite") else "n/a"
        print(f"    is_finite: {finite}", flush=True)
        if label == "(1,0)":
            print(f"    (1,1) in points (shear image): {element(L, [1, 1]) in points}", flush=True)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:150]}", flush=True)
