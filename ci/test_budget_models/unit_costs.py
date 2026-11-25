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
    # A primitive already measured keeps its value; a run measures only the keys that are missing.
    costs = json.loads(OUT.read_text())["seconds_per_call"]
    for n in (3, 4, 8, 10, 12, 18):
        a, b = integer_matrices(n)
        if f"zz_matrix_product_{n}" not in costs:
            costs[f"zz_matrix_product_{n}"] = median_seconds(lambda: a * b, 2000)
        q = a.change_ring(QQ)
        # Sage caches a matrix's rank, so each call ranks a fresh copy; the copy's cost is subtracted.
        copy_cost = median_seconds(lambda: q.__copy__(), 500)
        if f"qq_rank_{n}" not in costs:
            costs[f"qq_rank_{n}"] = median_seconds(lambda: q.__copy__().rank(), 500) - copy_cost
    v, w = vector(ZZ, range(12)), vector(ZZ, range(1, 13))
    if "zz_dot_12" not in costs:
        costs["zz_dot_12"] = median_seconds(lambda: v * w, 20000)

    e8 = 2 * matrix(ZZ, CartanMatrix(["E", 8]))
    if "pari_qfauto_e8" not in costs:
        costs["pari_qfauto_e8"] = median_seconds(lambda: pari(e8).qfauto(), 5)
    if "pari_qfisom_e8" not in costs:
        costs["pari_qfisom_e8"] = median_seconds(lambda: pari(e8).qfisom(e8), 5)
    a4 = 2 * matrix(ZZ, CartanMatrix(["A", 4]))
    if "pari_qfauto_a4" not in costs:
        costs["pari_qfauto_a4"] = median_seconds(lambda: pari(a4).qfauto(), 50)
    points = int(pari(e8).qfminim(8)[0])
    if "pari_qfminim_e8_point" not in costs:
        costs["pari_qfminim_e8_point"] = median_seconds(lambda: pari(e8).qfminim(8), 5) / points
    basis = integer_matrices(10)[0] + 20
    gram10 = basis * basis.transpose()
    if "pari_lllgram_10" not in costs:
        costs["pari_lllgram_10"] = median_seconds(lambda: pari(gram10).qflllgram(), 200)

    rank3 = matrix(ZZ, [[2, 1, 0], [1, 4, 1], [0, 1, 6]])
    if "pari_qfauto_rank3" not in costs:
        costs["pari_qfauto_rank3"] = median_seconds(lambda: pari(rank3).qfauto(), 500)
    if "pari_qfisom_rank3" not in costs:
        costs["pari_qfisom_rank3"] = median_seconds(lambda: pari(rank3).qfisom(rank3), 500)
    for n in (3, 18):
        snf = integer_matrices(n)[0]
        if f"pari_matsnf_{n}" not in costs:
            costs[f"pari_matsnf_{n}"] = median_seconds(lambda: pari(snf).matsnf(), 200)
    ternary = matrix(ZZ, [[1, 0, 0], [0, 1, 0], [0, 0, -13]])
    if "pari_qfsolve_3" not in costs:
        costs["pari_qfsolve_3"] = median_seconds(lambda: pari(ternary).qfsolve(), 200)

    if "libgap_call" not in costs:
        costs["libgap_call"] = median_seconds(lambda: libgap.IdentityMat(2), 20000)
    group = libgap.SymmetricGroup(527)
    generators = list(group.GeneratorsOfGroup())
    if "gap_permutation_image_527" not in costs:
        costs["gap_permutation_image_527"] = median_seconds(lambda: libgap.OnPoints(17, generators[0]), 20000)
    o = libgap.SO(1, 10, 2)
    if "gap_order_matrix_group_so_plus_10_2" not in costs:
        costs["gap_order_matrix_group_so_plus_10_2"] = median_seconds(lambda: libgap.Group(o.GeneratorsOfGroup()).Size(), 3)

    from sage_extensions.all import Lattices

    gram = [[0, 1, 0], [1, 0, 0], [0, 0, -2]]
    if "lattice_construction_rank_3" not in costs:
        costs["lattice_construction_rank_3"] = median_seconds(lambda: Lattices(ZZ)(gram), 50)
    lattice = Lattices(ZZ)(gram)
    x, y = lattice.module_generators()[0], lattice.module_generators()[1]
    if "lattice_element_add" not in costs:
        costs["lattice_element_add"] = median_seconds(lambda: x + y, 5000)
    if "lattice_element_pairing" not in costs:
        costs["lattice_element_pairing"] = median_seconds(lambda: lattice.b(x, y), 5000)

    from PyNormaliz import Cone

    rays = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]]
    if "normaliz_cone_6_rays" not in costs:
        costs["normaliz_cone_6_rays"] = median_seconds(lambda: Cone(cone=rays).SupportHyperplanes(), 50)

    # Scope G (hermitian-and-arithmetic-groups.md, "Unit costs of scope G").
    from sage.all import QuaternionAlgebra, QuadraticField, CyclotomicField, ComplexField, RealIntervalField, Polyhedron, identity_matrix

    h = QuaternionAlgebra(QQ, -1, 3)
    p1, p2 = h([1, 2, -1, 3]), h([2, -1, 1, 1])
    k5 = QuadraticField(5, "s")
    b5 = QuaternionAlgebra(k5, 11, 19)
    q1, q2 = b5([1, k5.gen(), 2, -1]), b5([k5.gen(), 1, -1, 3])
    k3 = QuadraticField(3, "s")
    m1, m2 = matrix(k3, 2, 2, [1, k3.gen(), 2, -1]), matrix(k3, 2, 2, [k3.gen(), 1, -1, 3])
    gi = CyclotomicField(4)
    g1, g2 = matrix(gi, 2, 2, [1, gi.gen(), 2, -1]), matrix(gi, 2, 2, [gi.gen(), 1, -1, 3])
    km5 = QuadraticField(-5, "w")
    lam, mu = km5([1, 1]), km5(3)
    d4 = matrix(ZZ, 4, 4, [2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 6, 0, 0, 0, 0, 6])
    g8 = matrix(ZZ, 8, 8, lambda i, j: 4 if i == j else (1 if abs(i - j) == 1 else 0))
    h8 = random_matrix(ZZ, 8, 12, x=-20, y=21)
    q16a, q16b = random_matrix(QQ, 16, 16, num_bound=5, den_bound=3), random_matrix(QQ, 16, 16, num_bound=5, den_bound=3)
    rif, c256 = RealIntervalField(128), ComplexField(256)
    r1, r2, z1, z2 = rif(1, 2), rif(3, 4), c256(1, 2), c256(3, 4)
    rays6 = [list(r) for r in identity_matrix(QQ, 6).rows()] + [[1, 1, -1, 0, 0, 1], [1, -1, 1, 1, 0, 0], [0, 1, 1, -1, 1, 0], [1, 0, 0, 1, 1, -1], [0, 0, 1, 1, -1, 1], [1, 1, 0, 0, 1, 1]]
    scope_g = {
        "qq_quaternion_product": (lambda: p1 * p2, 20000),
        "qq_quaternion_reduced_norm": (lambda: p1.reduced_norm(), 20000),
        "nf_quaternion_product_sqrt5": (lambda: q1 * q2, 5000),
        "nf_matrix_product_2_sqrt3": (lambda: m1 * m2, 5000),
        "nf_matrix_product_2_gaussian": (lambda: g1 * g2, 5000),
        "nf_ideal_equality_sqrt_minus5": (lambda: km5.ideal([lam, mu]) == km5.ideal(1), 200),
        "pari_qfminim_rank4_call": (lambda: pari(d4).qfminim(12, None, 2), 2000),
        "pari_qfminim_rank8_call": (lambda: pari(g8).qfminim(8, None, 2), 500),
        "pari_bnfinit_deg2": (lambda: pari("bnfinit(x^2-3)"), 50),
        "pari_bnfinit_deg4": (lambda: pari("bnfinit(x^4-16*x^2+4)"), 20),
        "zz_hermite_8x12": (lambda: h8.__copy__().hermite_form(), 500),
        "qq_matrix_product_16": (lambda: q16a * q16b, 2000),
        "rif128_mult": (lambda: r1 * r2, 50000),
        "complex256_mult": (lambda: z1 * z2, 50000),
        "ppl_cone_dim6_12_rays": (lambda: Polyhedron(rays=rays6, base_ring=QQ, backend="ppl").inequalities(), 50),
    }
    for key, (call, calls) in scope_g.items():
        if key not in costs:
            costs[key] = median_seconds(call, calls)
    for n in (8, 16):
        # The commutant equations F X - X F = 0 for one isometry F of size n (IsometryCommutant._integral_basis).
        F = identity_matrix(ZZ, n) + matrix(ZZ, n, n, lambda i, j: 1 if j == i + 1 else 0)
        rows = []
        for i in range(n):
            for j in range(n):
                row = [0] * (n * n)
                for k in range(n):
                    row[k * n + j] -= F[i, k]
                    row[i * n + k] += F[k, j]
                rows.append(row)
        equations = matrix(ZZ, rows)
        if f"zz_right_kernel_{n * n}" not in costs:
            costs[f"zz_right_kernel_{n * n}"] = median_seconds(lambda: equations.__copy__().right_kernel(), 3)

    # Scope G (hermitian-and-arithmetic-groups.md, "Unit costs of scope G").
    from sage.all import CyclotomicField, Polyhedron, QuadraticField, QuaternionAlgebra, identity_matrix

    h = QuaternionAlgebra(QQ, -1, 3)
    p1, p2 = h([1, 2, -1, 3]), h([2, -1, 1, 1])
    if "qq_quaternion_product" not in costs:
        costs["qq_quaternion_product"] = median_seconds(lambda: p1 * p2, 20000)
    k5 = QuadraticField(5, "s")
    b5 = QuaternionAlgebra(k5, 11, 19)
    q1, q2 = b5([1, k5.gen(), 2, -1]), b5([k5.gen(), 1, -1, 3])
    if "nf_quaternion_product_sqrt5" not in costs:
        costs["nf_quaternion_product_sqrt5"] = median_seconds(lambda: q1 * q2, 5000)
    k3 = QuadraticField(3, "s")
    m1, m2 = matrix(k3, 2, 2, [1, k3.gen(), 2, -1]), matrix(k3, 2, 2, [k3.gen(), 1, -1, 3])
    if "nf_matrix_product_2_sqrt3" not in costs:
        costs["nf_matrix_product_2_sqrt3"] = median_seconds(lambda: m1 * m2, 5000)
    gi = CyclotomicField(4)
    g1, g2 = matrix(gi, 2, 2, [1, gi.gen(), 2, -1]), matrix(gi, 2, 2, [gi.gen(), 1, -1, 3])
    if "nf_matrix_product_2_gaussian" not in costs:
        costs["nf_matrix_product_2_gaussian"] = median_seconds(lambda: g1 * g2, 5000)
    d4 = matrix(ZZ, 4, 4, [2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 6, 0, 0, 0, 0, 6])
    if "pari_qfminim_rank4_call" not in costs:
        costs["pari_qfminim_rank4_call"] = median_seconds(lambda: pari(d4).qfminim(12, None, 2), 2000)
    g8 = matrix(ZZ, 8, 8, lambda i, j: 4 if i == j else (1 if abs(i - j) == 1 else 0))
    if "pari_qfminim_rank8_call" not in costs:
        costs["pari_qfminim_rank8_call"] = median_seconds(lambda: pari(g8).qfminim(8, None, 2), 500)
    if "pari_bnfinit_deg2" not in costs:
        costs["pari_bnfinit_deg2"] = median_seconds(lambda: pari("bnfinit(x^2-3)"), 50)
    if "pari_bnfinit_deg4" not in costs:
        costs["pari_bnfinit_deg4"] = median_seconds(lambda: pari("bnfinit(x^4-16*x^2+4)"), 20)
    h8 = random_matrix(ZZ, 8, 12, x=-20, y=21)
    if "zz_hermite_8x12" not in costs:
        costs["zz_hermite_8x12"] = median_seconds(lambda: h8.__copy__().hermite_form(), 500)
    for n in (8, 16):
        f = random_matrix(ZZ, 2 * n * n, n * n, x=-1, y=2)
        if f"zz_right_kernel_{n * n}" not in costs:
            costs[f"zz_right_kernel_{n * n}"] = median_seconds(lambda: f.__copy__().right_kernel(), 3)
    q16a, q16b = random_matrix(QQ, 16, 16, num_bound=5, den_bound=3), random_matrix(QQ, 16, 16, num_bound=5, den_bound=3)
    if "qq_matrix_product_16" not in costs:
        costs["qq_matrix_product_16"] = median_seconds(lambda: q16a * q16b, 2000)
    rays6 = [list(r) for r in identity_matrix(QQ, 6).rows()] + [[1, 1, -1, 0, 0, 1], [1, -1, 1, 1, 0, 0], [0, 1, 1, -1, 1, 0], [1, 0, 0, 1, 1, -1], [0, 0, 1, 1, -1, 1], [1, 1, 0, 0, 1, 1]]
    if "ppl_cone_dim6_12_rays" not in costs:
        costs["ppl_cone_dim6_12_rays"] = median_seconds(lambda: Polyhedron(rays=rays6, base_ring=QQ, backend="ppl").inequalities(), 50)

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
