"""Per-primitive costs of the orbit-wise route at the target specimen's scale.

The target cell cone has 49,861 rays in dimension 10.  The synthetic
configuration used here has the same shape (rank 10, signature (9,1), rays in
the affine hyperplane ``x_10 = 1``, entries of size at most 3) and 44,640 rays;
it is *not* the A7 specimen and nothing here re-runs that specimen.

Measured here:

  * the seed linear program (one per facet-orbit expansion) with ``m``
    constraints in ``n`` variables, exact rational arithmetic through PPL;
  * the ridge-pencil operator's two real costs: the ``2 x n`` by ``n x m``
    exact integer product and the ``O(m)`` angular extreme scan;
  * the setwise-stabiliser (Schreier) work per facet orbit, which is
    ``|orbit| * |generators|`` packed permutation applications.
"""

from __future__ import annotations

import random
import time

from sage.all import QQ, ZZ, matrix, vector

from orbit_addm import primitive, seed_facet
from validate_specimens import shell_points


def timed(label, function, repeats=3):
    best = None
    value = None
    for _ in range(repeats):
        started = time.perf_counter()
        value = function()
        elapsed = time.perf_counter() - started
        best = elapsed if best is None else min(best, elapsed)
    print(f"{label:<58} {best * 1000:>9.1f} ms")
    return best, value


def ray_matrix(k, square):
    points = shell_points(k, square)
    return [tuple(list(point) + [1]) for point in points]


def measure_seed(rows, dimension):
    height = tuple([0] * (dimension - 1) + [1])

    def seed():
        return seed_facet(rows, height, dimension)

    elapsed, value = timed(f"seed LP: m={len(rows):>6} n={dimension} (exact PPL, {dimension} LPs)", seed)
    face = value
    print(f"{'':<58}   seed incidence size {len(face)}")
    return elapsed


def measure_pencil_costs(m, dimension, repeats=3):
    """The two real costs of one ridge pencil, as a function of the ray count."""
    rng = random.Random(20240607)
    rows = [tuple(rng.randint(-15, 15) for _ in range(dimension)) for _ in range(m)]
    ray_transpose = matrix(ZZ, rows).transpose()
    left = tuple(rng.randint(-3, 3) for _ in range(dimension))
    right = tuple(rng.randint(-3, 3) for _ in range(dimension))
    basis = matrix(ZZ, 2, dimension, [*left, *right])

    def product():
        return basis * ray_transpose

    elapsed_product, projection_matrix = timed(
        f"ridge pencil: 2x{dimension} . {dimension}x{m} exact integer product", product, repeats
    )

    def scan():
        best_max = None
        best_min = None
        for index in range(m):
            first = int(projection_matrix[0, index])
            second = int(projection_matrix[1, index])
            if first == 0 and second == 0:
                continue
            if best_max is None:
                best_max = best_min = (first, second)
                continue
            if best_max[0] * second - best_max[1] * first > 0:
                best_max = (first, second)
            if best_min[0] * second - best_min[1] * first < 0:
                best_min = (first, second)
        return best_max, best_min

    elapsed_scan, _ = timed(f"ridge pencil: O(m) angular extreme scan (m={m})", scan, repeats)

    def feasibility():
        a, b = 1, 1
        feasible = all(
            a * int(projection_matrix[0, index]) + b * int(projection_matrix[1, index]) >= 0
            for index in range(m)
        )
        return feasible

    elapsed_feasible, _ = timed(
        f"ridge pencil: two O(m) exact feasibility passes (m={m})", feasibility, repeats
    )
    return elapsed_product, elapsed_scan, elapsed_feasible


def measure_schreier(point_count, generator_count, orbit_size):
    """One facet orbit's setwise-stabiliser work: |orbit| * |gens| bit applications."""
    rng = random.Random(11)
    generators = []
    base = list(range(point_count))
    for _ in range(generator_count):
        shuffled = base[:]
        rng.shuffle(shuffled)
        generators.append(shuffled)
    face = set(rng.sample(range(point_count), max(2, point_count // 20)))

    def closure():
        members = {frozenset(face): ()}
        frontier = [frozenset(face)]
        while frontier and len(members) < orbit_size:
            current = frontier.pop()
            for generator in generators:
                image = frozenset(generator[position] for position in current)
                if image not in members:
                    members[image] = ()
                    frontier.append(image)
        return len(members)

    elapsed, count = timed(
        f"orbit closure: {generator_count} generators on {point_count} points, {orbit_size} members",
        closure,
        2,
    )
    print(f"{'':<58}   closed {count} members")
    return elapsed


def main():
    print("=== A7-scale (rank 10) primitives ===")
    rows = ray_matrix(9, 10)
    print(f"configuration: {len(rows)} rays, rank {matrix(QQ, rows).rank()}, dimension 10")
    measure_seed(rows, 10)
    print()
    for m in (2034, 7392, 22608, 44640):
        print(f"--- m = {m} ---")
        measure_pencil_costs(m, 10)
        print()
    print("=== facet-orbit bookkeeping at 49,861 points ===")
    measure_schreier(49861, 17, 5000)


if __name__ == "__main__":
    main()
