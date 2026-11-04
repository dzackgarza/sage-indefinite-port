"""The radical not in the final Gram coordinates: is the marked route presentational?"""
from sage.all import ZZ, matrix
from tests.acceptance.consumer import element, lattice, rows_of
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix

recursion = EichlerRecursion()

# Permuted presentation: coordinates (r, x) with q(r, x) = 2 x^2.  Radical is the FIRST coordinate.
L = lattice([[0, 0], [0, 2]])
a, b = element(L, [0, 1]), element(L, [1, 1])
print("L:", L, "| is_nondegenerate", L.is_nondegenerate(), flush=True)
print("  a = (r,x) = (0,1): q =", int(a.q()), "content =", int(a.content()), flush=True)
print("  b = (r,x) = (1,1): q =", int(b.q()), "content =", int(b.content()), flush=True)

try:
    witness = recursion.pointed_isometry(a, b)
    print("[radical-first (0,1)->(1,1)] ->", witness, flush=True)
    if witness is not None:
        image = witness(a)
        print("   acts correctly:", image == b,
              "| image", tuple(int(image.to_vector()(l)) for l in L.module_generating_set()),
              "| column matrix", _column_matrix(witness).rows(), flush=True)
except Exception as error:
    print(f"[radical-first (0,1)->(1,1)] raised {type(error).__name__}: {str(error)[:150]}", flush=True)

try:
    generators = tuple(recursion.pointed_automorphism_generators(a))
    print(f"[radical-first stabilizer of (0,1)] {len(generators)} generators; all fix it:",
          all(g(a) == a for g in generators), flush=True)
    if generators:
        from sage.all import MatrixGroup
        print("   order:", MatrixGroup([matrix(ZZ, rows_of(g)) for g in generators]).order(), flush=True)
except Exception as error:
    print(f"[radical-first stabilizer] raised {type(error).__name__}: {str(error)[:150]}", flush=True)

try:
    print("[radical-first O(L).group_generators()]",
          len(tuple(L.O().group_generators())), "generators", flush=True)
except Exception as error:
    print(f"[radical-first unmarked] raised {type(error).__name__}: {str(error)[:150]}", flush=True)
