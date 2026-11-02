"""Does the native lazy orbit of a marked vector exist for a degenerate lattice?"""
from tests.acceptance.consumer import element, lattice

L = lattice([[2, 0], [0, 0]])
v = element(L, [1, 0])
radical_vector = element(L, [0, 1])
print("L degenerate:", not L.is_nondegenerate(), flush=True)

for label, vector in (("(1,0) content 1", v), ("(0,1) in the radical, c=0", radical_vector)):
    try:
        orbit = L.O().orbit(vector)
        print(f"[{label}] orbit object: {type(orbit).__name__}", flush=True)
        points = orbit.points()
        print(f"    points(): {type(points).__name__}", flush=True)
        print(f"    representative: {orbit.representative()}", flush=True)
        print(f"    is_finite: {getattr(points, 'is_finite', lambda: 'n/a')()}", flush=True)
        print(f"    vector in points: {vector in points}", flush=True)
        print(f"    unrelated (3,0) in points: {element(L, [3, 0]) in points}", flush=True)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:150]}", flush=True)
