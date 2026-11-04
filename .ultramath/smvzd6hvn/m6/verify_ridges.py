"""Exact check of the ridge half of the hypothesis package: ridges of a facet are
the facets of the cone spanned by that facet's incidence *inside* its own span,
every ridge lies in exactly two facets, and each ridge's ambient normal cone is
2-dimensional with exactly two extreme rays = those two facets.

Pure Python, Fractions only.
"""
import itertools
from fractions import Fraction as F


def rref(rows, n):
    """Reduced row-echelon basis of the row space, with pivot columns."""
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
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    return m[:r], pivots


def rank(rows):
    return len(rref(rows, len(rows[0]))[0]) if rows else 0


def dot(a, b):
    return sum(F(x) * F(y) for x, y in zip(a, b))


def kernel_of(v, k):
    """A basis of the annihilator of a single nonzero vector in QQ^k (dim k-1)."""
    piv = next(i for i in range(k) if v[i] != 0)
    basis = []
    for j in range(k):
        if j == piv:
            continue
        w = [F(0)] * k
        w[piv] = -v[j] / v[piv]
        w[j] = F(1)
        basis.append(tuple(w))
    return basis


def same_dir(t, u):
    return t[0] * u[1] - t[1] * u[0] == 0 and t[0] * u[0] + t[1] * u[1] > 0


def cone2(constraints):
    cs = [(F(a), F(b)) for a, b in constraints if (a, b) != (0, 0)]
    out = []
    for c1, c2 in cs:
        for t in ((-c2, c1), (c2, -c1)):
            if all(a * t[0] + b * t[1] >= 0 for a, b in cs):
                if t != (0, 0) and not any(same_dir(t, u) for u in out):
                    out.append(t)
    return out


def facets_of_coordinate_cone(coords, k):
    """Facet normals of cone(coords) in QQ^k (pointed, full-dimensional), with zero sets."""
    found = {}
    if k < 2:
        return found
    for choice in itertools.combinations(range(len(coords)), k - 1):
        picked = [coords[i] for i in choice]
        if rank(picked) != k - 1:
            continue
        ker = kernel_of_multi(picked, k)
        assert len(ker) == 1
        nrm = ker[0]
        vals = [dot(nrm, x) for x in coords]
        if all(v >= 0 for v in vals):
            pass
        elif all(v <= 0 for v in vals):
            nrm = tuple(-x for x in nrm)
            vals = [-v for v in vals]
        else:
            continue
        zeros = tuple(sorted(i for i, v in enumerate(vals) if v == 0))
        if rank([coords[i] for i in zeros]) != k - 1:
            continue
        found[zeros] = nrm
    return found


def main():
    for name, homog in (
        ("hexagonal cone", [(1, 0, 1), (1, 1, 1), (0, 1, 1), (-1, 0, 1), (-1, -1, 1), (0, -1, 1)]),
        ("cube cone", [(x, y, z, 1) for x in (1, -1) for y in (1, -1) for z in (1, -1)]),
    ):
        rays = [tuple(F(x) for x in r) for r in homog]
        n = len(rays[0])
        print("=" * 70)
        print(f"{name}: rays={len(rays)} ambient={n} rank={rank(rays)}")

        # all facets in the ambient space, by the same (n-1)-subset method
        facets = facets_of_coordinate_cone(rays, n)
        print(f"  facets={len(facets)} incidence sizes={sorted(len(z) for z in facets)}")

        owner = {}
        for zeros, _normal in facets.items():
            sub = [rays[i] for i in zeros]
            # change coordinates into span(F): RREF basis, read coordinates at pivots
            base, pivots = rref(sub, n)
            k = len(base)
            assert k == n - 1
            coords = [tuple(r[p] for p in pivots) for r in sub]
            for j, back in enumerate(coords):  # sanity: coordinates reconstruct the ray
                rebuilt = [sum(back[i] * base[i][c] for i in range(k)) for c in range(n)]
                assert tuple(rebuilt) == sub[j]
            ridge_facets = facets_of_coordinate_cone(coords, k)
            assert ridge_facets, "the facet cone has facets"
            for local_zeros, _ln in ridge_facets.items():
                ridge = tuple(sorted(zeros[i] for i in local_zeros))
                assert rank([rays[i] for i in ridge]) == n - 2
                assert set(ridge) <= set(zeros)
                # ambient 2-dimensional normal cone at the ridge
                kernel = kernel_of_multi([rays[i] for i in ridge], n)
                assert len(kernel) == 2
                constraints = [(dot(kernel[0], r), dot(kernel[1], r)) for r in rays]
                extremes = cone2(constraints)
                assert len(extremes) == 2, f"{ridge}: {len(extremes)} extreme rays"
                normals = []
                for t in extremes:
                    nrm = tuple(t[0] * kernel[0][c] + t[1] * kernel[1][c] for c in range(n))
                    vals = [dot(nrm, r) for r in rays]
                    assert all(v >= 0 for v in vals)
                    normals.append(tuple(sorted(i for i, v in enumerate(vals) if v == 0)))
                owner.setdefault(ridge, set()).add(zeros)
                assert len(set(normals)) == 2
                owner.setdefault(("pair", ridge), set()).update(normals)
        pairs = {z: o for z, o in owner.items() if not isinstance(z[0], str)}
        hist = sorted({len(o) for o in pairs.values()})
        print(f"  distinct ridges={len(pairs)} multiplicity histogram={hist}")
        assert hist == [2], "every ridge lies in exactly two facets"
        checks = 0
        for ridge, owning in pairs.items():
            both = owner[("pair", ridge)]
            assert both == owning, (ridge, both, owning)
            checks += 1
        print(f"  on {checks} ridges: the two extreme rays of the ridge normal cone are")
        print("  exactly the two facets through that ridge")
        print("  OK")


def kernel_of_multi(rows, n):
    """Basis of {x : rows . x = 0} in QQ^n."""
    base, pivots = rref(rows, n)
    free = [c for c in range(n) if c not in pivots]
    out = []
    for fc in free:
        vec = [F(0)] * n
        vec[fc] = F(1)
        for i, pc in enumerate(pivots):
            vec[pc] = -base[i][fc]
        out.append(tuple(vec))
    return out


if __name__ == "__main__":
    main()
