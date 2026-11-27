r"""Prescribed starting vertices of Allcock's edgewalk (#31, T8).

`edgewalk.h:2100-2115` takes the starting vertex from the caller (`FileVertex`), builds its
simple cone with `get_simple_cone` (`:1848-2049`) and normalises the stored generator with
`RemoveFractionVector` (`:2107`, `:2143`); `:1985-1999` raises `NonReflectivityException`
with a reason and `StandardEdgewalkAnalysis` turns it into the report (`:2175-2183`). A
vertex is a *ray*: scaling the prescribed vector is not part of it.

Record 1 of `allcock_rank3_reflective.json` is the published lattice
`G = [[4,0,0],[0,0,-1],[0,-1,-2]]` with simple roots `(-1,-2,0), (1,0,0), (0,2,-1)`
(`references/vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/RK3_all`, line 1).
`e_2 = (0,1,0)` is a null primitive vector of it, so it is an admissible ideal corner;
`2*e_2` is the same ray. Both must produce the recorded chamber and a reflective verdict.

Implementer diagnostic: asserts the source rule above; it is not acceptance evidence under
CONTRIBUTING.md (which requires an independent test author).
"""

from sage.all import QQ, ZZ, matrix, vector
from sage_extensions.all import ZZ as PreambleZZ
from sage_extensions.all import HyperbolicLattices

from tests.acceptance.consumer import coordinates, lattice
from sage_indefinite_port.indefinite.edgewalk import edgewalk_fundamental_domain

GRAM = [[4, 0, 0], [0, 0, -1], [0, -1, -2]]
RECORDED = [(-1, -2, 0), (1, 0, 0), (0, 2, -1)]

form = matrix(ZZ, GRAM)
negative = sum(1 for value in matrix(QQ, GRAM).eigenvalues() if value < 0)
oriented = form if negative == 1 else -form
L = HyperbolicLattices(PreambleZZ)(lattice([list(row) for row in oriented.rows()]))
labels = list(L.module_generating_set())


def ray(coordinates_):
    return L.linear_combination({label: L.base_ring()(entry)
                                 for label, entry in zip(labels, coordinates_) if entry})


def root_gram(roots):
    return [[int(u * form * v) for v in (vector(ZZ, r) for r in roots)]
            for u in (vector(ZZ, r) for r in roots)]


def sign_equivalent(left, right):
    """Same chamber up to a reindexing and a choice of wall normal per root."""
    from itertools import permutations, product
    size = len(left)
    for order in permutations(range(size)):
        for signs in product((1, -1), repeat=size):
            if all(left[order[i]][order[j]] * signs[i] * signs[j] == right[i][j]
                   for i in range(size) for j in range(size)):
                return True
    return False


for coordinates_ in [(0, 1, 0), (0, 2, 0), (0, -3, 0)]:
    vertex = ray(coordinates_)
    assert vertex.q() == 0, (coordinates_, vertex.q())
    report = edgewalk_fundamental_domain(L, initial_vertex=vertex)
    assert report.is_reflective(), (coordinates_, report.reason_non_reflective())
    roots = [coordinates(root) for root in report.simple_roots()]
    assert len(roots) == 3, (coordinates_, roots)
    assert sign_equivalent(root_gram(roots), root_gram(RECORDED)), (coordinates_, roots)
    print("prescribed null ray", coordinates_, "->", roots)

print("all prescribed scalings of the ideal corner reproduce the recorded chamber")
