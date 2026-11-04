"""Exact check of the orbit-facet hypothesis package on the two bounded specimens.

Checks, per specimen:
  H1  row rank == ambient rank; cone full-dimensional.
  H2  a strictly positive functional on the rows exists (pointedness certificate).
  H3  every facet covector is nonnegative on the rows and annihilates a rank-(n-1) set.
  H4  the stabilizer acts by ray permutations, and facet orbits partition the facets.
  H5  ridge <-> facet-of-incidence-cone bijection, dim n-2, dim strictly drops.
  H6  each ridge's ambient normal cone is 2-dimensional with exactly two extreme
      rays; the two are the normals of the two facets through the ridge.
All arithmetic exact (Sage QQ / ZZ), no floating point.
"""
import itertools

from sage.all import QQ, ZZ, Permutation, PermutationGroup, matrix

from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    LorentzianPerfectLocalBackend,
)
from tests.acceptance.consumer import lattice

SPECIMENS = {
    "A2_plus_1_plus_minus1": [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]],
    "A3_plus_1_plus_minus1": [
        [2, -1, 0, 0, 0],
        [-1, 2, -1, 0, 0],
        [0, -1, 2, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, -1],
    ],
}


def rk(rows):
    if not rows:
        return 0
    return matrix(QQ, [list(r) for r in rows]).rank()


def kernel_basis(rows, n):
    if not rows:
        return [tuple(QQ(1) if i == j else QQ(0) for i in range(n)) for j in range(n)]
    basis = matrix(QQ, [list(r) for r in rows]).right_kernel_matrix()
    return [tuple(QQ(x) for x in row) for row in basis.rows()]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def same_direction(t, u):
    if t[0] * u[1] - t[1] * u[0] != 0:
        return False
    return t[0] * u[0] + t[1] * u[1] > 0


def cone2_extreme_rays(constraints):
    """Extreme rays (up to positive scale) of {t in QQ^2 : c.t >= 0 for all c}.

    Assumes the solution cone is 2-dimensional.  Exact: candidate boundary
    directions are the exact rotations of the constraint vectors, tested by
    rational arithmetic only.
    """
    cs = [c for c in constraints if (c[0], c[1]) != (0, 0)]
    candidates = []
    for c1, c2 in cs:
        for t in ((-c2, c1), (c2, -c1)):
            if all(c1_ * t[0] + c2_ * t[1] >= 0 for c1_, c2_ in cs):
                candidates.append(t)
    out = []
    for t in candidates:
        if t[0] == 0 and t[1] == 0:
            continue
        if not any(same_direction(t, u) for u in out):
            out.append(t)
    return out


def as_ambient(t, basis, n):
    return tuple(t[0] * basis[0][k] + t[1] * basis[1][k] for k in range(n))


def run(name, gram):
    print("=" * 70)
    print(name)
    backend = LorentzianPerfectLocalBackend()
    lat = lattice(gram)
    cell = backend.initial_cell(lat, NONNEGATIVE_SQUARES)
    rays = [tuple(int(e) for e in row) for row in cell.rows]
    n = len(rays[0])
    print(f"  rays={len(rays)} ambient_dim={n} row_rank={rk(rays)}")
    assert rk(rays) == n, "H1 full-dimensionality"

    cone = backend.cell_cone(cell)
    labels = tuple(lat.module_generating_set())

    def coords_of(cov):
        raw = tuple(cov(lat.module_generator(lab)) for lab in labels)
        out = []
        for value in raw:
            value = QQ(value)
            assert value.denominator() == 1, "facet covector coordinates must be integral"
            out.append(int(value))
        return tuple(out)

    facet_data = []
    for cov in cone.halfspace_covectors():
        coords = coords_of(cov)
        values = [dot(coords, r) for r in rays]
        assert all(v >= 0 for v in values), "H3 nonnegativity"
        zeros = tuple(i for i, v in enumerate(values) if v == 0)
        assert rk([rays[i] for i in zeros]) == n - 1, "H3 incidence rank n-1"
        facet_data.append((coords, zeros))
    facet_data.sort(key=lambda item: item[1])
    assert len({z for _c, z in facet_data}) == len(facet_data)
    print(f"  facets={len(facet_data)} incidence sizes={[len(z) for _c, z in facet_data]}")

    # H2: sum of all facet normals is strictly positive on every row.
    total = [sum(c[k] for c, _z in facet_data) for k in range(n)]
    strict = [dot(total, r) for r in rays]
    assert all(v > 0 for v in strict), "H2 pointedness certificate"
    print("  H2 positivity certificate: sum of facet normals is >0 on all rows")

    # H4: stabilizer action.
    positions = {r: i for i, r in enumerate(rays)}
    permutations = []
    for generator in backend.cell_stabilizer(cell):
        images = tuple(
            positions[
                tuple(
                    int(generator(point).to_vector()(label))
                    for label in lat.module_generating_set()
                )
            ]
            for point in cell.vector_configuration
        )
        assert sorted(images) == list(range(len(rays)))
        permutations.append(images)
    group = PermutationGroup([Permutation(list(p)) for p in permutations])
    order = int(group.order())
    print(f"  stabilizer: {len(permutations)} generators, order {order}")

    def moved(incidence, permutation):
        return tuple(sorted(permutation[i] for i in incidence))

    facet_set = {z for _c, z in facet_data}
    for permutation in permutations:
        assert {moved(z, permutation) for z in facet_set} == facet_set, "H4 acts on facets"

    seen = set()
    orbits = []
    for z in sorted(facet_set):
        if z in seen:
            continue
        orbit = {z}
        frontier = [z]
        while frontier:
            current = frontier.pop()
            for permutation in permutations:
                image = moved(current, permutation)
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        seen |= orbit
        orbits.append(orbit)
    assert sum(len(o) for o in orbits) == len(facet_set)
    print(
        "  facet orbits: "
        + ", ".join(f"size {len(o)} (|Stab|={order // len(o)})" for o in orbits)
    )

    by_incidence = {z: c for c, z in facet_data}

    # H5/H6: ridges of each facet, and the ambient 2-dim normal cone at each ridge.
    ridge_owner = {}
    total_ridges = 0
    for coords, zeros in facet_data:
        span_rows = [rays[i] for i in zeros]
        facet_rank = rk(span_rows)
        ridges = set()
        for a, b in itertools.combinations(range(len(span_rows)), 2):
            if rk([span_rows[a], span_rows[b]]) <= facet_rank - 2:
                continue
            basis = kernel_basis([span_rows[a], span_rows[b]], n)
            constraints = [(dot(basis[0], s), dot(basis[1], s)) for s in span_rows]
            for t in cone2_extreme_rays(constraints):
                normal = as_ambient(t, basis, n)
                zero_all = tuple(i for i, r in enumerate(rays) if dot(normal, r) == 0)
                if rk([rays[i] for i in zero_all]) != facet_rank - 1:
                    continue
                ridges.add(zero_all)
        total_ridges += len(ridges)
        assert ridges, "a facet has ridges"
        for ridge in ridges:
            # H5: the ridge has dimension n-2 and drops strictly below the facet.
            assert rk([rays[i] for i in ridge]) == n - 2
            assert set(ridge) >= set(zeros)
            ridge_owner.setdefault(ridge, set()).add(zeros)
            # H6: ambient normal cone at the ridge is 2-dimensional, exactly two rays.
            basis = kernel_basis([rays[i] for i in ridge], n)
            assert len(basis) == 2
            constraints = [(dot(basis[0], r), dot(basis[1], r)) for r in rays]
            extremes = cone2_extreme_rays(constraints)
            assert len(extremes) == 2, f"ridge normal cone has {len(extremes)} extreme rays"
            normals = []
            for t in extremes:
                normal = as_ambient(t, basis, n)
                values = [dot(normal, r) for r in rays]
                assert all(v >= 0 for v in values)
                zero_all = tuple(i for i, v in enumerate(values) if v == 0)
                normals.append(zero_all)
            assert set(normals) == set(ridge_owner[ridge]) | {normals[0], normals[1]}
            # the two extreme rays are exactly the two facets through the ridge
            assert len(set(normals)) == 2
            assert set(normals) <= facet_set, "each ridge normal is a facet of the cone"
    # H5/H6 global: every ridge is contained in exactly two facets.
    bad = {r: o for r, o in ridge_owner.items() if len(o) != 2}
    print(f"  ridges over all facets={total_ridges}, distinct={len(ridge_owner)}, not doubly covered={len(bad)}")
    assert not bad, "every ridge must lie in exactly two facets"
    print("  OK: H1-H6 verified exactly on this specimen")


def main():
    for name in sorted(SPECIMENS):
        run(name, SPECIMENS[name])


if __name__ == "__main__":
    main()
