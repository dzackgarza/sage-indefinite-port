"""Where does the current degenerate pointed code stop being right?

The four retained <2>+0 checks already pass. This battery pushes the same construction onto cases where the
issue's congruence criterion s - D r in cR is genuinely exercised: rank-2 radical, x=0, larger C.
"""
from sage.all import ZZ, matrix, MatrixGroup
from tests.acceptance.consumer import element, lattice, rows_of
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix

recursion = EichlerRecursion()


def check(label, gram, source_coords, target_coords, expected):
    L = lattice(gram)
    source = element(L, source_coords)
    target = element(L, target_coords)
    try:
        witness = recursion.pointed_isometry(source, target)
    except Exception as error:
        print(f"[{label}] raised {type(error).__name__}: {str(error)[:110]} (expected {expected})", flush=True)
        return
    if witness is None:
        verdict = "None"
        ok = expected is False
    else:
        image = witness(source)
        ok = (image == target) and expected is True
        verdict = f"witness -> {tuple(int(image.to_vector()(l)) for l in L.module_generating_set())}"
    print(f"[{label}] {verdict} | expected {expected} | {'OK' if ok else 'WRONG'}", flush=True)


# A. the retained examples
check("A1 <2>+0 (1,0)->(1,1)", [[2, 0], [0, 0]], [1, 0], [1, 1], True)
check("A2 <2>+0 (5,1)->(5,6)", [[2, 0], [0, 0]], [5, 1], [5, 6], True)
check("A3 <2>+0 (5,1)-/->(5,2)", [[2, 0], [0, 0]], [5, 1], [5, 2], False)

# B. rank-2 radical: c = 5, r = (1,0). d11 = 1 mod 5, d21 = 0 mod 5 is reachable, so positive.
check("B1 rank2 rad (5,1,0)->(5,6,5)", [[2, 0, 0], [0, 0, 0], [0, 0, 0]], [5, 1, 0], [5, 6, 5], True)
# C. same but first radical column would be (2,0): not primitive, so no unimodular D.
check("B2 rank2 rad (5,1,0)-/->(5,2,0)", [[2, 0, 0], [0, 0, 0], [0, 0, 0]], [5, 1, 0], [5, 2, 0], False)

# D. x = 0 in the radical: the condition is exactly Dr = s.
check("C1 rank2 rad (0,1,0)->(0,0,1)", [[2, 0, 0], [0, 0, 0], [0, 0, 0]], [0, 1, 0], [0, 0, 1], True)
check("C2 rank1 rad (0,1)->(0,-1)", [[2, 0], [0, 0]], [0, 1], [0, -1], True)

# E. larger C, content 1: any radical coordinate is reachable.
check("D1 C=<2>+<2> (1,0,0)->(1,0,7)", [[2, 0, 0], [0, 2, 0], [0, 0, 0]], [1, 0, 0], [1, 0, 7], True)
check("D2 C=<2>+<2> (1,0,0)-/->(0,1,0)", [[2, 0, 0], [0, 2, 0], [0, 0, 0]], [1, 0, 0], [0, 1, 0], True)
