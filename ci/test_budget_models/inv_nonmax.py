# Invariant: order of the GF(2) image of the U(1,1;Z[i]) generators conjugated into the Z[2i]^2 frame (as in the test).
from sage.all import matrix, QQ, ZZ, GF, libgap, block_matrix, identity_matrix
from sage_indefinite_port.indefinite.binary_hyperbolic_unitary import generators
from sage_extensions.algebras.order_hermitian_modules import OrderHermitianModule

def tm(g):
    rows = []
    for i in range(2):
        re, im = [], []
        for j in range(2):
            a, b = g[i, j].vector()
            re += [a, -b]; im += [b, a]
        rows += [re, im]
    return matrix(ZZ, rows)

half = QQ(1) / 2
ion = matrix(QQ, 2, 2, [0, -2, half, 0]); z = matrix(QQ, 2, 2)
iact = block_matrix(QQ, [[ion, z], [z, ion]], subdivide=False)
gram = matrix(ZZ, 4, 4, [0, 0, 2, 0, 0, 0, 0, 8, 2, 0, 0, 0, 0, 8, 0, 0])
module = OrderHermitianModule(gram, [2 * iact], [identity_matrix(QQ, 4), iact])
sat = matrix(QQ, module.saturation_basis)
change = sat.inverse() * matrix.diagonal(QQ, [1, half, 1, half])
moves = [change * tm(g).change_ring(QQ) * change.inverse() for g in generators(1).values()]
print("generators", len(moves), "all integral", all(m.denominator() == 1 for m in moves))
print("GF(2) image order", libgap.Group([libgap(matrix(ZZ, m).change_ring(GF(2))) for m in moves]).Size())
