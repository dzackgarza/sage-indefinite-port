"""Does the split-radical path now produce the pointwise stabilizer of a radical-marked vector?"""
from sage.all import ZZ, matrix
from tests.acceptance.consumer import element, lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix

RECURSION = EichlerRecursion()
specimens = (
    ("C=<2>, R=Z   (0,1)", [[2, 0], [0, 0]], [0, 1]),
    ("C=<2>+<2>, R=Z   (0,0,1)", [[2, 0, 0], [0, 2, 0], [0, 0, 0]], [0, 0, 1]),
    ("C=<2>, R=Z^2  (0,1,0)", [[2, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 1, 0]),
    ("C=<2>, R=Z^2  (0,1,2) non-primitive", [[2, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 1, 2]),
)
for label, gram, coords in specimens:
    L = lattice(gram)
    v = element(L, coords)
    try:
        generators = tuple(RECURSION.pointed_automorphism_generators(v))
        fixes = all(g(v) == v for g in generators)
        print(f"[{label}] {len(generators)} generators | all fix v: {fixes}", flush=True)
        for g in generators[:6]:
            print(f"    {_column_matrix(g).rows()}", flush=True)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:130]}", flush=True)
