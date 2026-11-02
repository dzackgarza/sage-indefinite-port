"""Which Sage integer-linear-algebra idioms give a unimodular completion and the kernel basis?"""
from sage.all import ZZ, matrix, vector

x0 = vector(ZZ, [2, 5])
A = matrix(ZZ, [x0])
print("A:", A, flush=True)
for name, fn in (
    ("A.kernel().basis_matrix()", lambda: A.kernel().basis_matrix()),
    ("A.right_kernel_matrix()", lambda: A.right_kernel_matrix()),
    ("A.smith_form()[2] (V)", lambda: matrix(ZZ, [x0]).smith_form()[2]),
):
    try:
        value = fn()
        print(f"  {name} -> {value} shape {value.nrows()}x{value.ncols()}", flush=True)
    except Exception as error:
        print(f"  {name} raised {type(error).__name__}: {str(error)[:90]}", flush=True)

K = A.kernel().basis_matrix()
U = A.stack(K)
print("completion:", U.rows(), "det:", U.det(), flush=True)

# a rank-3 C with radical rank 2
x0b = vector(ZZ, [3, 4, 5])
Ab = matrix(ZZ, [x0b])
Kb = Ab.kernel().basis_matrix()
Ub = Ab.stack(Kb)
print("rank-3 completion det:", Ub.det(), "rows", Ub.rows(), flush=True)
