# Invariants (not timings) of C(A, sigma) for the anisotropic skew-Hermitian fixture of test_commutant_groups.sage.
from sage.all import QuaternionAlgebra, QQ, ZZ, matrix, RealField
from sage_extensions.algebras.isometry_commutant import IsometryCommutant
from sage_indefinite_port.indefinite.clifford_commutant import CliffordModel
from sage_indefinite_port.indefinite.fuchsian_domain import fuchsian_places
A = QuaternionAlgebra(QQ, -1, 3); ONE, I, J, K = A.basis(); B = (ONE, I, J, K)
Z0 = 0 * ONE
PAIRS = tuple((b, Z0) for b in B) + tuple((Z0, b) for b in B)
QUADS = tuple((x, (Z0, Z0)) for x in PAIRS) + tuple(((Z0, Z0), x) for x in PAIRS)
co = lambda q: [QQ(c) for c in list(q)]
pc = lambda x: co(x[0]) + co(x[1])
qc = lambda X: pc(X[0]) + pc(X[1])
qm = lambda mv: matrix(ZZ, [qc(mv(X)) for X in QUADS]).transpose()
add = lambda x, y: (x[0] + y[0], x[1] + y[1])
right = lambda d: qm(lambda X: tuple((x[0] * d, x[1] * d) for x in X))
src = [right(I), right(2 + J), qm(lambda X: (add(X[0], X[1]), X[1])), qm(lambda X: (X[0], add(X[1], X[0])))]
al, be = I, 3 * I + J
def form(X, Y):
    h = lambda x, y: y[0].conjugate() * al * x[0] + y[1].conjugate() * be * x[1]
    return (h(X[0], Y[1]) - h(X[1], Y[0])).reduced_trace()
gram = matrix(ZZ, [[form(X, Y) for Y in QUADS] for X in QUADS])
m = CliffordModel(IsometryCommutant(gram, src))
C = m.algebra
print("delta", m.delta, "field disc", m.field.discriminant(), "class number", m.field.class_number())
print("ramified", C.ramified_places(inf=True))
print("fuchsian places", fuchsian_places(C, RealField(64)) is not None)
