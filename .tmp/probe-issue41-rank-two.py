"""The two #41 requirements the retained rank-one examples do not touch.

(1) rank-two radical: a refusal that is NOT about primitivity of a column -- r = 0 mod c but s != 0 mod c.
(2) nontrivial shear kernel: it needs rank(C) >= 2, since for rank(C) = 1 every B is determined by B(x0)
    and c*B(x0) = 0 forces B = 0.
"""
from sage.all import MatrixGroup, ZZ, matrix
from tests.acceptance.consumer import element, lattice, rows_of
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix

recursion = EichlerRecursion()

R2 = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]          # C = <2>, R = Z^2
L = lattice(R2)


def show(label, a, b, expected):
    source, target = element(L, a), element(L, b)
    try:
        witness = recursion.pointed_isometry(source, target)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:120]}", flush=True)
        return
    if witness is None:
        print(f"[{label}] None | expected {expected} | {'OK' if expected is False else 'WRONG'}", flush=True)
    else:
        image = witness(source)
        cols = _column_matrix(witness)
        print(f"[{label}] witness -> {tuple(int(image.to_vector()(l)) for l in L.module_generating_set())}"
              f" | D-block {[cols[1, 1], cols[2, 1], cols[1, 2], cols[2, 2]]}"
              f" | expected {expected} | {'OK' if (image == target and expected) else 'WRONG'}", flush=True)


print("--- rank-two radical ---", flush=True)
show("transport (5,1,0)->(5,6,5)", [5, 1, 0], [5, 6, 5], True)
show("refusal  (5,5,0)->(5,0,1)", [5, 5, 0], [5, 0, 1], False)

print("--- nontrivial shear kernel: C of rank 2, radical Z ---", flush=True)
L2 = lattice([[2, 0, 0], [0, 2, 0], [0, 0, 0]])
v = element(L2, [1, 0, 0])
try:
    generators = tuple(recursion.pointed_automorphism_generators(v))
    print(f"  {len(generators)} generators; all fix v: {all(g(v) == v for g in generators)}", flush=True)
    for index, g in enumerate(generators):
        print(f"    g{index}: column matrix {_column_matrix(g).rows()}", flush=True)
    forms = [matrix(ZZ, rows_of(g)) for g in generators]
    if forms:
        group = MatrixGroup(forms)
        print(f"  group order: {group.order() if group.is_finite() else 'infinite'}", flush=True)
except Exception as error:
    print(f"  raised {type(error).__name__}: {str(error)[:160]}", flush=True)
