"""Validation and cost measurements for code/primitives.py.

Two independent engines are used as ground truth:

  * analytic -- the cone over the n-cube, whose facets are known in closed form
    (x_0 -+ x_j >= 0, 2n of them, each tight on 2**(n-1) rays), at the specimen's
    dimension n = 10;
  * SageMath's PPL-backed facet enumeration on random small cones, i.e. the same
    "generic engine" role Normaliz plays in the port.

For every facet the script recomputes the normal from the incidence set alone
(P3), re-certifies it (P1) and decides adjacency twice, by pairings (P2) and by
the reference's joint-rank criterion.  Then it times the primitives at the sizes
the A7 specimen reaches (49 861 rays, rank 10, incidences up to thousands).

Run with plain python3 in this VM (sage.all is importable there).
"""

from __future__ import annotations

import json
import random
import time

import numpy as np
from sage.all import QQ, ZZ, Polyhedron, matrix, vector

from primitives import (
    negative_support,
    supports_are_disjoint,
    adjacent_from_zero_sets,
    zero_support,
    adjacent_rank_test,
    exact_minor_is_nonzero,
    is_facet_certificate,
    pairings,
    pairings_exact,
    primitive_normal_from_incidence,
    rank_lower_certificate,
)


def cube_cone(dimension: int):
    """Rays (1, v) for v in {-1,1}^dimension and the closed-form facet normals."""
    rays = []
    for mask in range(2**dimension):
        rays.append((1,) + tuple(1 if (mask >> bit) & 1 else -1 for bit in range(dimension)))
    normals = []
    for j in range(1, dimension + 1):
        plus = [0] * (dimension + 1)
        plus[0], plus[j] = 1, -1
        minus = [0] * (dimension + 1)
        minus[0], minus[j] = 1, 1
        normals.extend([tuple(plus), tuple(minus)])
    return rays, normals


def incidence_of(normal, rays):
    return frozenset(index for index, row in enumerate(rays)
                     if sum(a * b for a, b in zip(normal, row, strict=True)) == 0)


def check_facet_family(name, rays, normals):
    """(P1)+(P3) on a known facet family: recover every normal from its incidence."""
    recovered, certified, primary = 0, 0, 0
    for normal in normals:
        incidence = incidence_of(normal, rays)
        rebuilt = primitive_normal_from_incidence([rays[i] for i in incidence], len(normal))
        if rebuilt in (tuple(normal), tuple(-entry for entry in normal)):
            recovered += 1
        if is_facet_certificate(normal, incidence, rays):
            certified += 1
        if all(sum(a * b for a, b in zip(rebuilt, row, strict=True)) >= 0 for row in rays):
            primary += 1
    return {"case": name, "facets": len(normals), "normal_recovered_up_to_sign": recovered,
            "certificate_accepted": certified, "rebuilt_primitive_nonnegative": primary}


def adjacency_agreement(name, rays, normals):
    """(P2): the bitset joint-rank decision, and how much of J it had to read."""
    pairs = adjacent_pairs = 0
    examined = 0
    supports = [zero_support(normal, rays) for normal in normals]
    for i in range(len(normals)):
        for j in range(i + 1, len(normals)):
            pairs += 1
            adjacent, rows = adjacent_from_zero_sets(supports[i], supports[j], rays, len(normals[0]))
            adjacent_pairs += adjacent
            examined = max(examined, rows)
    return {"case": name, "facet_pairs": pairs, "adjacent": adjacent_pairs,
            "max_joint_rows_examined": examined}


def sage_facets(rays):
    """Facet normals of the cone over `rays`, from Sage's exact H-representation.

    Sage writes each facet as A x + b >= 0 with rational A; a cone through the
    origin has b = 0, so the primitive integer version of A is the facet normal,
    already oriented nonnegatively on every ray.
    """
    from math import gcd, lcm

    polyhedron = Polyhedron(rays=[vector(ZZ, row) for row in rays], base_ring=QQ)
    normals = []
    for inequality in polyhedron.inequalities():
        if QQ(inequality.b()) != 0:
            continue
        coefficients = [QQ(value) for value in inequality.A()]
        denominator = 1
        for value in coefficients:
            denominator = lcm(denominator, value.denominator())
        integral = [int(value * denominator) for value in coefficients]
        divisor = 0
        for entry in integral:
            divisor = gcd(divisor, abs(entry))
        if divisor == 0:
            continue
        primitive = tuple(entry // divisor for entry in integral)
        if any(sum(a * b for a, b in zip(primitive, row, strict=True)) < 0 for row in rays):
            primitive = tuple(-entry for entry in primitive)
        if primitive not in normals:
            normals.append(primitive)
    return normals


def random_cone(dimension, count, seed, entry=3):
    """A random cone over a point configuration, i.e. the specimen's own shape.

    Rays are (1, x) with x a small integer vector, so they lie in the affine
    hyperplane x_0 = 1 and the cone they generate is pointed and full-dimensional
    -- exactly the structure `LorentzianPerfectCell.__post_init__` certifies for a
    perfect cell (augmented rank n).
    """
    generator = random.Random(seed)
    rays = []
    seen = set()
    while len(rays) < count:
        point = (1,) + tuple(generator.randint(-entry, entry) for _ in range(dimension - 1))
        if point not in seen:
            seen.add(point)
            rays.append(point)
    return rays, sage_facets(rays)


def time_certified_rank(dimension, sizes, seeds=3):
    """Cost of one exact rank certificate at incidence sizes up to the specimen's."""
    rows_all = []
    generator = random.Random(7)
    for _ in range(max(sizes)):
        rows_all.append(tuple(generator.randint(-15, 15) for _ in range(dimension)))
    results = []
    for size in sizes:
        subset = rows_all[:size]
        best = min(_timed_certificate(subset, dimension) for _ in range(seeds))
        sage_qq = min(_timed_sage_rank(subset, QQ) for _ in range(seeds))
        sage_zz = min(_timed_sage_rank(subset, ZZ) for _ in range(seeds))
        results.append({"rows": size, "certificate_s": best, "sage_QQ_rank_s": sage_qq,
                        "sage_ZZ_rank_s": sage_zz})
    return results


def _timed_certificate(rows, dimension):
    start = time.perf_counter()
    indices, _rank_mod_p = rank_lower_certificate(rows, dimension)
    elapsed = time.perf_counter() - start
    assert len(indices) == dimension
    return elapsed


def _timed_sage_rank(rows, ring):
    start = time.perf_counter()
    matrix(ring, [list(row) for row in rows]).rank()
    return time.perf_counter() - start


def time_pairing_tests(dimension, ray_count, repeats=5):
    """Cost of one (P2) adjacency decision at the A7 ray count."""
    generator = random.Random(11)
    rays = [tuple(generator.randint(-15, 15) for _ in range(dimension)) for _ in range(ray_count)]
    first = tuple(generator.randint(0, 1) for _ in range(dimension))
    second = tuple(generator.randint(0, 1) for _ in range(dimension))
    best = min(_timed_pairing(first, second, rays) for _ in range(repeats))
    sage_best = min(_timed_sage_pairing(first, rays) for _ in range(repeats))
    return {"rays": ray_count, "dimension": dimension,
            "pairing_adjacency_s": best, "sage_ZZ_single_pairing_s": sage_best}


def _timed_pairing(first, second, rays):
    start = time.perf_counter()
    pairings(first, rays); pairings(second, rays)
    return time.perf_counter() - start


def _timed_sage_pairing(normal, rays):
    start = time.perf_counter()
    matrix(ZZ, [list(row) for row in rays]) * vector(ZZ, list(normal))
    return time.perf_counter() - start



def time_adjacency_tests(dimension, ray_count, orbit_counts=(50, 200, 600)):
    """(P2) cost of the exact bitset adjacency test at the A7 ray count."""
    generator = random.Random(19)
    rays = [tuple(generator.randint(-15, 15) for _ in range(dimension)) for _ in range(ray_count)]
    results = []
    for orbits in orbit_counts:
        normals = [tuple(generator.choice((0, 0, 0, 0, 0, 1, -1)) for _ in range(dimension))
                   for _ in range(orbits)]
        start = time.perf_counter()
        supports = [zero_support(normal, rays) for normal in normals]
        build = time.perf_counter() - start
        start = time.perf_counter()
        pairs = adjacent = 0
        for i in range(orbits):
            for j in range(i + 1, orbits):
                pairs += 1
                adjacent += adjacent_from_zero_sets(supports[i], supports[j], rays, dimension)[0]
        elapsed = time.perf_counter() - start
        results.append({"case": "adjacency_cost", "facet_orbits": orbits, "rays": ray_count,
                        "zero_support_build_s": round(build, 4),
                        "pairs": pairs, "adjacent": adjacent,
                        "exact_tests_total_s": round(elapsed, 4),
                        "exact_test_us": round(1e6 * elapsed / max(pairs, 1), 2)})
    return results

def main():
    report = []
    # 1. dimension 10, the specimen's dimension, with 256-ray incidences.
    rays, normals = cube_cone(9)
    report.append(check_facet_family("cube_cone_dim10", rays, normals))
    report.append(adjacency_agreement("cube_cone_dim10", rays, normals))
    # 2. the documented small specimen shape: rank 4, 27 rays, 6 facets.
    rays4, normals4 = cube_cone(3)
    report.append(check_facet_family("cube_cone_dim4", rays4, normals4))
    report.append(adjacency_agreement("cube_cone_dim4", rays4, normals4))
    # 3. random degenerate cones through Sage's own facet enumeration.
    for seed, (dimension, count) in enumerate(((4, 20), (5, 30), (6, 40))):
        rays_r, normals_r = random_cone(dimension, count, seed)
        report.append(check_facet_family(f"random_dim{dimension}_{count}rays", rays_r, normals_r))
        if len(normals_r) <= 80:
            report.append(adjacency_agreement(f"random_dim{dimension}_{count}rays", rays_r, normals_r))
    # 4. exactness cross-check of the rank certificate against Sage's rank over Q.
    generator = random.Random(3)
    mismatches = 0
    for _ in range(120):
        dimension = generator.randint(2, 10)
        count = generator.randint(dimension, dimension + 6)
        rows = [tuple(generator.randint(-15, 15) for _ in range(dimension)) for _ in range(count)]
        indices, rank_mod_p = rank_lower_certificate(rows, dimension, full_rank=True)
        sage_rank = matrix(QQ, [list(row) for row in rows]).rank()
        certified = len(indices) if indices else None
        if certified is not None and certified != sage_rank:
            mismatches += 1
        if rank_mod_p != sage_rank:
            mismatches += 1
    report.append({"case": "rank_certificate_vs_sage_QQ", "trials": 120, "mismatches": mismatches})
    # 5. cost measurements.
    report.extend({"case": "certified_rank_cost", **row}
                  for row in time_certified_rank(10, (64, 256, 1024, 4096, 16384)))
    report.append({"case": "pairing_cost", **time_pairing_tests(10, 49861)})
    report.extend(time_bitset_closure(10, 49861))
    report.extend(time_adjacency_tests(10, 49861))
    print(json.dumps(report, indent=1))




def time_bitset_closure(dimension, ray_count, orbit_counts=(50, 200, 600)):
    """(P2) cost model of the orbit-graph closure at the A7 ray count.

    One pairing pass and one packed negative support per facet orbit, then one
    big-integer AND per facet-orbit pair -- no rank test anywhere.
    """
    generator = random.Random(17)
    rays = [tuple(generator.randint(-15, 15) for _ in range(dimension)) for _ in range(ray_count)]
    results = []
    for orbits in orbit_counts:
        # Sparse-ish, small normals: the shape of real facet covectors, whose
        # negative supports are a small fraction of the 49 861 rays.
        normals = [tuple(generator.choice((0, 0, 0, 1, -1)) for _ in range(dimension))
                   for _ in range(orbits)]
        start = time.perf_counter()
        supports = [negative_support(normal, rays) for normal in normals]
        build = time.perf_counter() - start
        start = time.perf_counter()
        pairs = 0
        adjacent = 0
        for i in range(orbits):
            left = supports[i]
            for j in range(i + 1, orbits):
                pairs += 1
                adjacent += supports_are_disjoint(left, supports[j])
        test = time.perf_counter() - start
        results.append({
            "case": "bitset_closure", "facet_orbits": orbits, "rays": ray_count,
            "support_build_s": round(build, 4), "pairs": pairs,
            "pair_tests_total_s": round(test, 4),
            "pair_test_us": round(1e6 * test / max(pairs, 1), 2),
            "adjacent_pairs": adjacent,
        })
    return results


if __name__ == "__main__":
    main()
