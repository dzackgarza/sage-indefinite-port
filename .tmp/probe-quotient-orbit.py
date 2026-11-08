"""The quotient L[n]/O(L) as a lazy orbit set, and the maintained Hom/evaluation route."""
from sage_extensions.all import GObjects, Modules, Sets, ZZ as OwnedZZ
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
from tests.acceptance.consumer import element, lattice

RECURSION = EichlerRecursion()
L = lattice([[2, 0], [0, 0]])


def relation(left, right):
    return RECURSION.pointed_isometry(left, right) is not None


print("--- the square-n slice L[n] as a set ---", flush=True)
for label, fn in (
    ("Sets().condition_set(L, pred)", lambda: Sets().condition_set(L, lambda v: int(v.q()) == 2)),
    ("L.condition_set(pred)", lambda: L.condition_set(lambda v: int(v.q()) == 2)),
):
    try:
        subset = fn()
        print(f"  {label} -> {type(subset).__name__}; (1,0) in: {element(L, [1, 0]) in subset}; "
              f"(0,1) in: {element(L, [0, 1]) in subset}", flush=True)
    except Exception as error:
        print(f"  {label} raised {type(error).__name__}: {str(error)[:110]}", flush=True)

print("--- the quotient L[2]/O(L) ---", flush=True)
try:
    slice_set = Sets().condition_set(L, lambda v: int(v.q()) == 2)
    action = GObjects(L.O(), Sets()).on_set(slice_set, lambda g, x: g(x), orbit_relation=relation)
    quotient = action.orbits()
    print("  quotient:", type(quotient).__name__, flush=True)
    print("  surface:", [m for m in ("orbit_of", "orbits", "is_finite", "cardinality", "__iter__")
                         if hasattr(quotient, m)], flush=True)
    orbit = quotient.orbit_of(element(L, [1, 0]))
    print("  orbit_of((1,0)):", type(orbit).__name__, "| (1,1) in it:",
          element(L, [1, 1]) in orbit.points(), flush=True)
except Exception as error:
    print(f"  raised {type(error).__name__}: {str(error)[:200]}", flush=True)

print("--- maintained Hom + evaluation ---", flush=True)
try:
    ring = OwnedZZ
    free_one = Modules(ring).an_object()
    print("  Modules(ZZ).an_object():", free_one, flush=True)
    hom = Modules(ring).Mor(free_one, free_one)
    print("  Modules(ZZ).Mor(F,F):", type(hom).__name__, "| evaluation:",
          hasattr(hom, "evaluation"), flush=True)
except Exception as error:
    print(f"  Hom route raised {type(error).__name__}: {str(error)[:160]}", flush=True)
