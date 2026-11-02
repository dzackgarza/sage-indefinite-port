"""Cost scaling of the orbit-wise route against the complete dual description.

Each specimen is a perfect-cell-shaped configuration: the lattice vectors of
norm ``square - 1`` in the affine hyperplane ``x_d = 1`` for the lattice of Gram
``I_{d-1} (+) <-1>``, with the signed coordinate permutations as stabiliser
generators.  ``rays`` is the input size, ``facets`` what the *complete* dual
description has to produce, ``orbits`` what the orbit-wise route produces.

The complete-dual-description column is PPL's exact H-enumeration of the same
ray matrix (``Polyhedron(rays=...)``); it is skipped for specimens where it is
known to be out of reach, and the count is then reported as unavailable.
"""

from __future__ import annotations

import resource
import sys
import time

from orbit_addm import Stats, facet_orbits, facet_normal, incidence, pair
from validate_specimens import (
    full_facet_incidences,
    generators_of,
    shell_points,
)


def peak_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def run(k, square, *, reference: bool):
    points = shell_points(k, square)
    rows = tuple(tuple(list(point) + [1]) for point in points)
    dimension = k + 1
    height = tuple([0] * k + [1])
    generators, _actions = generators_of(rows, k)
    before = peak_mib()
    stats = Stats()
    started = time.perf_counter()
    orbits = facet_orbits(rows, generators, dimension, height=height, stats=stats)
    elapsed = time.perf_counter() - started
    after = peak_mib()

    members = set()
    for record in orbits:
        assert not (members & set(record["members"]))
        members |= set(record["members"])
        normal = record["normal"]
        assert all(pair(normal, row) >= 0 for row in rows)
        assert incidence(rows, normal) == record["face"]
        assert facet_normal(rows, record["face"], dimension) == normal

    facets = None
    reference_elapsed = None
    if reference:
        ref_started = time.perf_counter()
        facets = len(full_facet_incidences(rows))
        reference_elapsed = time.perf_counter() - ref_started
        assert facets == len(members), (facets, len(members))
    sizes = sorted(len(record["members"]) for record in orbits)
    print(
        f"k={k:>2} shell={square:>3} rays={len(rows):>6} dim={dimension:>2} "
        f"facets={('?' if facets is None else facets):>7} orbits={len(orbits):>3} "
        f"max_orbit={max(sizes):>6} nodes={stats.nodes:>5} seeds={stats.seeds:>5} "
        f"pencils={stats.pencil_calls:>5} "
        f"adm={elapsed:>8.2f}s full={('n/a' if reference_elapsed is None else format(reference_elapsed, '.2f') + 's'):>9} "
        f"peak={after:>7.0f}MiB (+{after - before:.0f})",
        flush=True,
    )
    return {
        "rays": len(rows),
        "dimension": dimension,
        "facets": facets,
        "orbits": len(orbits),
        "orbit_sizes": sizes,
        "nodes": stats.nodes,
        "seeds": stats.seeds,
        "pencils": stats.pencil_calls,
        "elapsed": elapsed,
        "reference_elapsed": reference_elapsed,
        "peak_mib": after,
    }


SERIES = [
    (2, 2, True),
    (3, 2, True),
    (4, 4, True),
    (5, 4, True),
    (6, 4, True),
    (7, 4, True),
    (8, 4, True),
    (9, 4, True),
]


def main():
    results = []
    for k, square, reference in SERIES:
        if reference is None:
            reference = "--reference" in sys.argv
        results.append(run(k, square, reference=reference))
    return results


if __name__ == "__main__":
    main()
