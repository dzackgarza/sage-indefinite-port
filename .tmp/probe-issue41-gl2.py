"""Does the port genuinely use a non-identity integral D in GL_2(Z)?

(5,1,0) -> (5,6,5) is satisfiable with D = I because s - r = (5,5) is already in 5R.
(5,1,0) -> (5,2,0) is NOT: s - r = (1,0) is not in 5R, so D r must differ from r modulo 5.
"""
from sage.all import ZZ
from tests.acceptance.consumer import element, lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _column_matrix

recursion = EichlerRecursion()
L = lattice([[2, 0, 0], [0, 0, 0], [0, 0, 0]])

for label, a, b in (("D=I admissible (5,1,0)->(5,6,5)", [5, 1, 0], [5, 6, 5]),
                    ("D!=I required   (5,1,0)->(5,2,0)", [5, 1, 0], [5, 2, 0])):
    source, target = element(L, a), element(L, b)
    witness = recursion.pointed_isometry(source, target)
    if witness is None:
        print(f"[{label}] None", flush=True)
        continue
    columns = _column_matrix(witness)
    # block order: C = coordinate 0, R = coordinates 1,2
    A = columns[0, 0]
    D = columns.submatrix(1, 1, 2, 2)
    B = columns.submatrix(1, 0, 2, 1)   # rows R, column C
    image = witness(source)
    D_times_r = (D * __import__("sage.all", fromlist=["vector"]).vector(ZZ, [a[1], a[2]]))
    s = __import__("sage.all", fromlist=["vector"]).vector(ZZ, [b[1], b[2]])
    print(f"[{label}] acts {image == target} | A={A} | D={D.rows()} | D r = {tuple(D_times_r)}"
          f" | s - D r = {tuple(s - D_times_r)} | B={B.rows()}", flush=True)
