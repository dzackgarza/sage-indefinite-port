"""Issue #41 baseline: what the current public pointed operations do on L = <2> + 0."""
import traceback

from tests.acceptance.consumer import (
    automorphisms, element, isomorphisms, lattice, pointed, pointed_category, rows_of,
)
from sage_indefinite_port.indefinite.recursive import EichlerRecursion

L = lattice([[2, 0], [0, 0]])
print("L:", L, "| rank", L.rank(), "| is_nondegenerate", L.is_nondegenerate(),
      "| is_definite", L.is_definite(), flush=True)
v, w = element(L, [1, 0]), element(L, [1, 1])
u1, u6, u2 = element(L, [5, 1]), element(L, [5, 6]), element(L, [5, 2])
for name, vector in (("(1,0)", v), ("(1,1)", w), ("(5,1)", u1), ("(5,6)", u6), ("(5,2)", u2)):
    print(f"  {name}: q={int(vector.q())} content={int(vector.content())} div={int(vector.div())}",
          flush=True)

recursion = EichlerRecursion()


def attempt(label, fn):
    try:
        value = fn()
        print(f"[{label}] -> {value}", flush=True)
        return value
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:180]}", flush=True)
        traceback.print_exc()
        return None


print("--- unmarked route (must keep working) ---", flush=True)
attempt("L.O().group_generators()", lambda: tuple(L.O().group_generators()))

print("--- pointed route (the issue's obligation) ---", flush=True)
witness = attempt("recursion.pointed_isometry((1,0) -> (1,1))", lambda: recursion.pointed_isometry(v, w))
attempt("recursion.pointed_isometry((5,1) -> (5,6))", lambda: recursion.pointed_isometry(u1, u6))
attempt("recursion.pointed_isometry((5,1) -> (5,2))", lambda: recursion.pointed_isometry(u1, u2))
attempt("recursion.pointed_automorphism_generators((1,0))",
        lambda: tuple(recursion.pointed_automorphism_generators(v)))

print("--- public pointed consumers ---", flush=True)
category = pointed_category(2)
source = pointed(L, v)
target = pointed(L, w)
attempt("isomorphisms(Isom((L,(1,0)),(L,(1,1)))).is_empty()",
        lambda: isomorphisms(category, source, target).is_empty())
attempt("automorphisms(Aut(L,(1,0))).group_generators()",
        lambda: tuple(automorphisms(category, source).group_generators()))
