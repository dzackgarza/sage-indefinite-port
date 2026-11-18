# Invariant: orders of the mod-2 images of the Picard generators, by GAP on Z-coordinate matrices over GF(2).
from sage.all import matrix, GF, QQ, ZZ, libgap
from sage_indefinite_port.indefinite import gaussian_picard as gp, euclidean_picard as ep

def zrep(g, basis):
    change = matrix(QQ, [list(b.vector()) for b in basis]).inverse()
    def mult(x):  # Z-matrix of multiplication by x on basis
        return matrix(QQ, [list((x * b).vector() * change) for b in basis]).transpose()
    n = g.nrows()
    return matrix(QQ, [[mult(g[i, j])[r, c] for j in range(n) for c in range(len(basis))] for i in range(n) for r in range(len(basis))])

F = gp.FIELD
gb = [F(1), gp.I]
G = [zrep(g, gb).change_ring(ZZ).change_ring(GF(2)) for g in gp.GENERATORS.values()]
print("gaussian mod 2 image", libgap.Group([libgap(m) for m in G]).Size())
M = zrep(gp.GENERATORS["M"], gb).change_ring(ZZ).change_ring(GF(2))
print("cyclic <M> mod 2", libgap.Group([libgap(M)]).Size(), "gens", len(G))
b7 = ep.order_basis(7)
Z = [zrep(g, b7).change_ring(ZZ).change_ring(GF(2)) for g in ep.generators(7).values()]
print("zhao d=7 mod 2 image", libgap.Group([libgap(m) for m in Z]).Size(), "gens", len(Z))
