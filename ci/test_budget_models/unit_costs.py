"""Measure the unit costs of the primitives the budget models price, and write unit_costs.json.

Each entry times one primitive operation at the size the models use, in one warm Sage session: the median
of repeated runs, in seconds per call. It times primitives only, never a test, so the measured slowness of
a test is not absorbed into its budget. Run in the Sage container from the repository root:

    sage-docker-python ci/test_budget_models/unit_costs.py
"""

import json
import platform
import statistics
import time
from pathlib import Path

import sage.all  # noqa: F401
from sage.all import QQ, ZZ, libgap, matrix, pari, random_matrix, vector
from sage.combinat.root_system.cartan_matrix import CartanMatrix

OUT = Path(__file__).with_name("unit_costs.json")


def median_seconds(call, calls_per_run, runs=7):
    """Median over runs of the time per call, after one warm-up run."""
    for _ in range(calls_per_run):
        call()
    samples = []
    for _ in range(runs):
        start = time.perf_counter()
        for _ in range(calls_per_run):
            call()
        samples.append((time.perf_counter() - start) / calls_per_run)
    return statistics.median(samples)


def integer_matrices(n):
    a = random_matrix(ZZ, n, n, x=-5, y=6)
    b = random_matrix(ZZ, n, n, x=-5, y=6)
    return a, b


def main():
    costs = {}
    for n in (3, 4, 8, 10, 12, 18):
        a, b = integer_matrices(n)
        costs[f"zz_matrix_product_{n}"] = median_seconds(lambda: a * b, 2000)
        q = a.change_ring(QQ)
        # Sage caches a matrix's rank, so each call ranks a fresh copy; the copy's cost is subtracted.
        copy_cost = median_seconds(lambda: q.__copy__(), 500)
        costs[f"qq_rank_{n}"] = median_seconds(lambda: q.__copy__().rank(), 500) - copy_cost
    v, w = vector(ZZ, range(12)), vector(ZZ, range(1, 13))
    costs["zz_dot_12"] = median_seconds(lambda: v * w, 20000)

    e8 = 2 * matrix(ZZ, CartanMatrix(["E", 8]))
    costs["pari_qfauto_e8"] = median_seconds(lambda: pari(e8).qfauto(), 5)
    costs["pari_qfisom_e8"] = median_seconds(lambda: pari(e8).qfisom(e8), 5)
    a4 = 2 * matrix(ZZ, CartanMatrix(["A", 4]))
    costs["pari_qfauto_a4"] = median_seconds(lambda: pari(a4).qfauto(), 50)
    points = int(pari(e8).qfminim(8)[0])
    costs["pari_qfminim_e8_point"] = median_seconds(lambda: pari(e8).qfminim(8), 5) / points
    basis = integer_matrices(10)[0] + 20
    gram10 = basis * basis.transpose()
    costs["pari_lllgram_10"] = median_seconds(lambda: pari(gram10).qflllgram(), 200)

    costs["libgap_call"] = median_seconds(lambda: libgap.IdentityMat(2), 20000)
    group = libgap.SymmetricGroup(527)
    generators = list(group.GeneratorsOfGroup())
    costs["gap_permutation_image_527"] = median_seconds(lambda: libgap.OnPoints(17, generators[0]), 20000)
    o = libgap.SO(1, 10, 2)
    costs["gap_order_matrix_group_so_plus_10_2"] = median_seconds(lambda: libgap.Group(o.GeneratorsOfGroup()).Size(), 3)

    from sage_extensions.all import Lattices

    gram = [[0, 1, 0], [1, 0, 0], [0, 0, -2]]
    costs["lattice_construction_rank_3"] = median_seconds(lambda: Lattices(ZZ)(gram), 50)
    lattice = Lattices(ZZ)(gram)
    x, y = lattice.module_generators()[0], lattice.module_generators()[1]
    costs["lattice_element_add"] = median_seconds(lambda: x + y, 5000)
    costs["lattice_element_pairing"] = median_seconds(lambda: lattice.b(x, y), 5000)

    from PyNormaliz import Cone

    rays = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]]
    costs["normaliz_cone_6_rays"] = median_seconds(lambda: Cone(cone=rays).SupportHyperplanes(), 50)

    OUT.write_text(
        json.dumps(
            {
                "source": "ci/test_budget_models/unit_costs.py",
                "machine": platform.platform(),
                "seconds_per_call": {k: float(v) for k, v in sorted(costs.items())},
            },
            indent=1,
        )
        + "\n"
    )
    for key, value in sorted(costs.items()):
        print(f"{value:.3e}\t{key}")


if __name__ == "__main__":
    main()
