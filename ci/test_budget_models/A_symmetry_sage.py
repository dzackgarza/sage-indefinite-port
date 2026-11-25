"""Integral symmetry group Sym_L(P) of the Coxeter chamber of each root_systems_103 case."""
import json, sys
from sage.all import Graph, matrix, QQ, ZZ, PermutationGroup
cases = json.load(open(sys.argv[1]))
out = {}
for c in cases:
    G = matrix(ZZ, c["gram"]); R = [vector for vector in c["roots"]]
    M = matrix(ZZ, R); F = M * G * M.transpose()
    k = len(R)
    g = Graph([list(range(k)), [(i, j, F[i, j]) for i in range(k) for j in range(i + 1, k) if F[i, j]]], format="vertices_and_edges")
    norms = sorted({F[i, i] for i in range(k)})
    A = g.automorphism_group(partition=[[i for i in range(k) if F[i, i] == x] for x in norms], edge_labels=True)
    piv = M.pivot_rows(); B = M.matrix_from_rows(piv); Binv = B.change_ring(QQ).inverse()
    integral = []
    for s in A:
        img = matrix(ZZ, [R[s(i + 1) - 1] if A.domain()[0] == 1 else R[s(i)] for i in piv])
        T = Binv * img
        if all(x in ZZ for x in T.list()):
            integral.append(s)
    H = PermutationGroup(integral, domain=A.domain()) if integral else None
    gens = [[ (s(i + 1) - 1) if A.domain()[0] == 1 else s(i) for i in range(k)] for s in (H.gens() if H else [])]
    out[c["id"]] = {"diagram_order": int(A.order()), "integral_order": int(H.order()) if H else 1, "gens": gens}
json.dump(out, open(sys.argv[2], "w"))
print(len(out), max(v["diagram_order"] for v in out.values()))
