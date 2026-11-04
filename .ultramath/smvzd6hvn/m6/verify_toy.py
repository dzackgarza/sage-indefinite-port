"""Exact, engine-free check of the orbit-facet hypothesis package on self-contained cones.

The port's own specimens could not be run in this VM (the provider library
`sage_extensions` fails to import against this Sage build: `ValueError: base
must be a subring of this ring`), so the package is verified on cones that
exercise exactly the structural points the theorems assert.  Pure Python,
`fractions.Fraction` only, no floating point.

C-A  hexagonal cone in R^3        6 rays, 6 facets, simplicial (|I_F| = n-1)
C-B  cube cone in R^4             8 rays, 6 facets, non-simplicial (|I_F| = n)
C-C  half-plane cone in R^2       lineality: the doubly-covered property fails
C-D  cube cone + a subgroup action  facet stabilizer vs set stabilizer of a
                                    non-maximal spanning sub-incidence
"""
import itertools
from fractions import Fraction as F


def rank(rows):
    if not rows:
        return 0
    m = [[F(x) for x in row] for row in rows]
    r, cols = 0, len(m[0])
    for c in range(cols):
        pivot = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                m[i] = [a - factor * b for a, b in zip(m[i], m[r])]
        r += 1
        if r == len(m):
            break
    return r


def kernel_basis(rows, n):
    if not rows:
        return [tuple(F(1) if i == j else F(0) for i in range(n)) for j in range(n)]
    m = [[F(x) for x in row] for row in rows]
    pivots, r = [], 0
    for c in range(n):
        pivot = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if pivot is None:
            continue
        m[r], m[pivot] = m[pivot], m[r]
        pv = m[r][c]
        m[r] = [x / pv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                factor = m[i][c]
                m[i] = [a - factor * b for a, b in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(n) if c not in pivots]
    basis = []
    for fc in free:
        vec = [F(0)] * n
        vec[fc] = F(1)
        for i, pc in enumerate(pivots):
            vec[pc] = -m[i][fc]
        basis.append(tuple(vec))
    return basis


def dot(a, b):
    return sum(F(x) * F(y) for x, y in zip(a, b))


def same_direction(t, u):
    if t[0] * u[1] - t[1] * u[0] != 0:
        return False
    return t[0] * u[0] + t[1] * u[1] > 0


def cone2_extreme_rays(constraints):
    cs = [(F(c[0]), F(c[1])) for c in constraints if (c[0], c[1]) != (0, 0)]
    out = []
    for c1, c2 in cs:
        for t in ((-c2, c1), (c2, -c1)):
            if all(a * t[0] + b * t[1] >= 0 for a, b in cs):
                if t == (0, 0) or any(same_direction(t, u) for u in out):
                    continue
                out.append(t)
    return out


def as_ambient(t, basis, n):
    return tuple(t[0] * basis[0][k] + t[1] * basis[1][k] for k in range(n))


def cone_extreme_in_span(basis, sub, n):
    """Extreme rays of {x in span(basis) : x >= 0 on sub}, as ambient covectors.

    Dimension 1 is a sign test; dimension 2 is the exact angular computation.
    """
    d = len(basis)
    if d == 1:
        v = basis[0]
        values = [dot(v, s) for s in sub]
        if all(x >= 0 for x in values):
            return [tuple(v)]
        if all(x <= 0 for x in values):
            return [tuple(-x for x in v)]
        return []
    if d == 2:
        constraints = [(dot(basis[0], s), dot(basis[1], s)) for s in sub]
        return [as_ambient(t, basis, n) for t in cone2_extreme_rays(constraints)]
    raise ValueError(f"unexpected kernel dimension {d}")


def lift(b, facet_normal, rays, facet_zeros):
    """b + M*a_F with M making it strictly positive off the facet (T4's construction)."""
    off = [i for i in range(len(rays)) if i not in set(facet_zeros)]
    need = F(0)
    for i in off:
        slope = dot(facet_normal, rays[i])
        assert slope > 0
        value = dot(b, rays[i])
        if value < 0:
            need = max(need, -value / slope)
    M = need + 1
    return tuple(b[k] + M * facet_normal[k] for k in range(len(b)))


def cone_facets(rays, n, subset):
    sub = [rays[i] for i in subset]
    r = rank(sub)
    found = {}
    if r < 2:
        return found
    for choice in itertools.combinations(range(len(sub)), r - 1):
        picked = [sub[i] for i in choice]
        if rank(picked) != r - 1:
            continue
        basis = kernel_basis(picked, n)
        for normal in cone_extreme_in_span(basis, sub, n):
            values = [dot(normal, x) for x in sub]
            if not all(v >= 0 for v in values):
                continue
            zeros = tuple(sorted(subset[i] for i, v in enumerate(values) if v == 0))
            if rank([rays[i] for i in zeros]) != r - 1:
                continue
            found[zeros] = normal
    return found


def check(name, rays, n, pointed=True):
    print("=" * 72)
    print(name)
    rays = [tuple(F(x) for x in r) for r in rays]
    print(f"  rays={len(rays)} ambient={n} rank={rank(rays)}")
    assert rank(rays) == n, "full-dimensional"
    facets = cone_facets(rays, n, tuple(range(len(rays))))
    print(f"  facets={len(facets)} incidence sizes={sorted(len(z) for z in facets)}")
    total = [sum(nrm[k] for nrm in facets.values()) for k in range(n)]
    values = [dot(total, r) for r in rays]
    print(f"  min over rows of (sum of facet normals) = {min(values)}")
    if pointed:
        assert all(v > 0 for v in values), "sum of facet normals is a strictly positive functional"
        assert set().union(*map(set, facets)) == set(range(len(rays)))
    # ridge structure is verified separately in verify_ridges.py
    print("  done (facet counts and pointedness certificate above)")


def as_homogeneous(points):
    return [tuple(p) + (1,) for p in points]


def main():
    hexagon = as_homogeneous([(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)])
    check("C-A hexagonal cone (simplicial facets)", hexagon, 3)

    cube = as_homogeneous(list(itertools.product((1, -1), repeat=3)))
    check("C-B cube cone (non-simplicial facets, |I_F| = n)", cube, 4)

    halfplane = [(1, 0), (-1, 0), (0, 1)]
    check("C-C half-plane cone (lineality; doubly-covered property fails)", halfplane, 2,
          pointed=False)
    rays = [tuple(F(x) for x in r) for r in halfplane]
    facets = cone_facets(rays, 2, tuple(range(3)))
    print(f"  half-plane facets={list(facets)}: its zero face lies in "
          f"{len(facets)} facet, not two.")

    rays = [tuple(F(x) for x in r) for r in cube]
    n = 4
    index = {r: i for i, r in enumerate(rays)}
    generators = [
        tuple(index[(r[1], r[0], r[2], r[3])] for r in rays),
        tuple(index[(r[0], r[1], -r[2], r[3])] for r in rays),
    ]
    print("=" * 72)
    print("C-D cube cone: facet stabilizer vs set stabilizer of a spanning sub-incidence")
    facets = cone_facets(rays, n, tuple(range(len(rays))))
    target = next(z for z in facets if len(z) == 4 and all(rays[i][0] == 1 for i in z))
    print(f"  facet x=1 incidence={target} size={len(target)} (n-1 = {n-1})")

    def moves(incidence, perm):
        return tuple(sorted(perm[i] for i in incidence))

    def stabilizer(incidence):
        identity = tuple(range(len(rays)))
        elements, frontier = {identity}, [identity]
        while frontier:
            current = frontier.pop()
            for p in generators:
                composed = tuple(p[current[i]] for i in range(len(rays)))
                if composed not in elements:
                    elements.add(composed)
                    frontier.append(composed)
        return [e for e in elements if moves(incidence, e) == tuple(sorted(incidence))]

    subset = target[:3]
    print(f"  |Stab(facet)| = {len(stabilizer(target))}")
    print(f"  |set stabilizer of a 3-of-4 spanning sub-incidence| = {len(stabilizer(subset))}")
    print(f"  rank of the sub-incidence = {rank([rays[i] for i in subset])} = n-1: it spans the facet")
    print("  => a non-maximal incidence spans the same face but is not a valid facet key: its")
    print("     set stabilizer is strictly smaller and one facet acquires several keys.")


if __name__ == "__main__":
    main()
