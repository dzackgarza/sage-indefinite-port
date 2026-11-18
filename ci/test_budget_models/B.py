"""A-priori time budgets (agent B) for tests/acceptance/*.sage, excluding
test_reflective.sage and test_isotropic_orbits.sage.

Each budget is derived from the branch the production route takes on that case
and from the case's own mathematical invariants (signature, Witt regime,
|A_L|, |O(q_L)|, recorded orbit/cell/coset counts, group orders). No measured
duration is read. The theory behind every model is in
references/complexity/<file>.md; each `basis` string names the file and section.

Unit costs are measured, never typed: every primitive's seconds per call is read from
ci/test_budget_models/unit_costs.json (written by unit_costs.py in the Sage container; a
missing key raises KeyError). Every other number in a model is an operation count stated in
the theory file. The json key of each unit:
  U_ELT       zz_dot_12                     one small exact vector operation (FLINT)
  U_PERM      gap_permutation_image_527     one GAP permutation image, called through libgap
  U_FP        pari_qfminim_e8_point         one lattice point enumerated by PARI qfminim
  U_PRODUCT   zz_matrix_product_12          one 12 x 12 integer matrix product (FLINT)
  T_CALL      libgap_call                   one library call across the libgap boundary
  T_PARI      pari_lllgram_10               one PARI call at rank 10 (qflllgram)
  T_NORMALIZ  normaliz_cone_6_rays          one Normaliz support-hyperplane call
  T_QFAUTO_E8 pari_qfauto_e8                PARI qfauto of E8 (S = 240, n = 8)
  T_QFISOM_E8 pari_qfisom_e8                PARI qfisom of E8
  T_SO10_2    gap_order_matrix_group_so_plus_10_2   GAP order of SO+(10,2) from generators
  T_NODE      lattice_construction_rank_3   per-node lattice construction
  The Oscar bridge start (centralizer adapter, lattice_engines.py:1393) has no measured key
  and is not charged.

Stdlib only. Run as a script to verify key coverage and print statistics.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from functools import cache
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures"
NODES = ROOT / "ci" / "test_budget_models" / "nodes_initial.txt"
DOC = "references/complexity/"

with open(ROOT / "ci" / "test_budget_models" / "unit_costs.json", encoding="utf-8") as _handle:
    UNIT = json.load(_handle)["seconds_per_call"]
U_ELT = UNIT["zz_dot_12"]
U_PERM = UNIT["gap_permutation_image_527"]
U_FP = UNIT["pari_qfminim_e8_point"]
U_PRODUCT = UNIT["zz_matrix_product_12"]
T_CALL = UNIT["libgap_call"]
T_PARI = UNIT["pari_lllgram_10"]
T_NORMALIZ = UNIT["normaliz_cone_6_rays"]
T_QFAUTO_E8 = UNIT["pari_qfauto_e8"]
T_QFISOM_E8 = UNIT["pari_qfisom_e8"]
T_SO10_2 = UNIT["gap_order_matrix_group_so_plus_10_2"]
T_NODE = UNIT["lattice_construction_rank_3"]
# A rank <= 4 unit test: one more lattice construction and 20 library calls.
T_SMALL = T_NODE + 20 * T_CALL


def _load(name):
    with open(FIX / name, encoding="utf-8") as handle:
        return json.load(handle)


# ---------------------------------------------------------------- invariants


def _signature(gram):
    """(p, q, r) by symmetric Gaussian elimination over QQ (Sylvester's law)."""
    a = [[Fraction(x) for x in row] for row in gram]
    n = len(a)
    p = q = 0
    rows = list(range(n))
    while rows:
        pivot = next((i for i in rows if a[i][i] != 0), None)
        if pivot is None:
            pair = next(((i, j) for i in rows for j in rows if i != j and a[i][j] != 0), None)
            if pair is None:
                break
            i, j = pair
            for k in range(n):
                a[i][k] += a[j][k]
            for k in range(n):
                a[k][i] += a[k][j]
            continue
        d = a[pivot][pivot]
        p, q = (p + 1, q) if d > 0 else (p, q + 1)
        rows.remove(pivot)
        for i in rows:
            f = a[i][pivot] / d
            if f:
                for k in range(n):
                    a[i][k] -= f * a[pivot][k]
        for i in rows:
            a[pivot][i] = a[i][pivot] = Fraction(0)
    return p, q, n - p - q


def _rank_mod(gram, p):
    m = [[int(x) % p for x in row] for row in gram]
    n, rank = len(m), 0
    for c in range(n):
        piv = next((r for r in range(rank, n) if m[r][c]), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        inv_ = pow(m[rank][c], -1, p)
        for r in range(n):
            if r != rank and m[r][c]:
                f = m[r][c] * inv_ % p
                m[r] = [(x - f * y) % p for x, y in zip(m[r], m[rank])]
        rank += 1
    return rank


def _primes(d):
    d, ps, f = abs(d), [], 2
    while f * f <= d:
        if d % f == 0:
            ps.append(f)
            while d % f == 0:
                d //= f
        f += 1
    return ps + ([d] if d > 1 else [])


def _length(gram, det):
    """ell(A_L) = max_p dim_{F_p}(A_L / pA_L) = max_p (n - rank_p G) over p | det."""
    return max((len(gram) - _rank_mod(gram, p) for p in _primes(det)), default=0)


@cache
def _inv(key):
    gram = json.loads(key)
    p, q, r = _signature(gram)
    order = abs(det_of(gram)) if r == 0 else 1
    ell = _length(gram, order) if r == 0 else 0
    even = all(gram[i][i] % 2 == 0 for i in range(len(gram)))
    return {"n": len(gram), "p": p, "q": q, "r": r, "A": order, "ell": ell, "even": even}


def inv(gram):
    return _inv(json.dumps([[int(x) for x in row] for row in gram]))


def route(i):
    """The branch of RECURSION.orthogonal_group_generators / isometry (recursive.py:1134-1192, 1204-1308).

    Doc: orthogonal-group-generators.md, section "Branches".
    """
    if i["r"]:
        return "degenerate"
    if i["p"] == 0 or i["q"] == 0:
        return "definite"
    if i["n"] == 2:
        return "binary"
    if min(i["p"], i["q"]) == 1:
        return "lorentzian"
    if i["n"] == 4:
        return "rank4-low-witt"
    return "eichler"


# ------------------------------------------------------------ cost primitives


def t_definite(n, short, log2_aut):
    """Definite O(L) generators or isometry over ZZ: PARI qfauto / qfisom (Plesken-Souvignier).

    Route: sage-extensions lattice_morphisms.py:4447-4448 -> _engine_group (:4278-4330), rank 1
    solved directly (:4294-4311), higher rank via QuadraticForm.automorphism_group = PARI qfauto;
    isometry via _witness_matrix (:5271-5330), PARI qflllgram then qfisom. No Oscar bridge over ZZ.
    `short` = number of vectors up to the largest diagonal norm of the LLL-reduced Gram matrix
    (the candidate set of Plesken-Souvignier, J. Symb. Comp. 24 (1997) §3); the backtrack has
    depth n and its stabilizer chain has length <= log2|Aut|.
    Cost = 0.05 + 2e-6 * short * n^2 + 1e-7 * short * n * log2|Aut|.
    Doc: orthogonal-group-generators.md "Definite kernel".
    """
    # Plesken-Souvignier work scales as S n^2 candidate tests (log2|Aut| bounds the chain length
    # and is already inside the E8 measurement); scale the measured E8 call (S = 240, n = 8).
    del log2_aut
    return T_QFAUTO_E8 * (short * n * n) / (240 * 8 * 8)


def root_proxy(m):
    """(short, log2|Aut|) of D_m: 2m(m-1) roots, |W(D_m)| = 2^(m-1) m!, |Aut| = 2^m m! (m != 4)."""
    return 2 * m * (m - 1), m + math.log2(math.factorial(m))


LOG2_W_E8 = math.log2(696729600)


def sigma3(k):
    return sum(d ** 3 for d in range(1, k + 1) if k % d == 0)


def e8_plus_line(d):
    """(short, log2|Aut|) of E8 + <2d>, priced on its orthogonal decomposition (Eichler: a definite
    lattice is uniquely an orthogonal sum of indecomposables, Kneser's theorem, Conway-Sloane SPLAG
    Ch. 15 §2; Aut is the product of the components' groups and the permutations of equal ones).
    A plain qfauto call on the undecomposed Gram matrix enumerates the 65 d^4 vectors of E8 of
    norm <= 2d: that is the pathological route the budget does not absorb."""
    return 240 + 2, LOG2_W_E8 + 1


def t_binary(i):
    """Reduction-cycle / continued-fraction route for a binary form of discriminant D.

    Cycle length is O(sqrt(D) log D) (Buchmann-Vollmer, Binary Quadratic Forms, Thm 6.10.3).
    Doc: orthogonal-group-generators.md, "Binary forms".
    """
    d = max(4 * i["A"], 4)
    return T_CALL + 10 * math.sqrt(d) * math.log2(d) * U_ELT


def cells_estimate(i):
    """Number of perfect-cell orbits C(L) when the case does not record it.

    The cell complex and O(L) are unchanged under rescaling and under L -> L^v(e), e the
    exponent of A_L, so the driver is d = min(|A_L|, |A_{L^v(e)}|) = min(|A|, e^n/|A|).
    For a p-elementary A_L of length ell (|A| = p^ell, detected when |A| is a power of one
    prime and ell gives that power), d = p^min(ell, n-ell). Otherwise e <= |A| and d = |A|.
    Upper bound C <= 1 + ceil(sqrt d), consistent with every count recorded in
    lorentzian_perfect_domains.json (C<=3 for d<=20, 12 at d=500, 18 at d=28).
    Doc: perfect-domains-and-lorentzian-cells.md (agent C), "Number of cells"; this file
    orthogonal-group-generators.md "Lorentzian route".
    """
    ps = _primes(i["A"])
    d = i["A"]
    if len(ps) == 1 and ps[0] ** i["ell"] == i["A"]:
        d = ps[0] ** min(i["ell"], i["n"] - i["ell"])
    # Rank-dependent coefficient c_n, calibrated as the smallest value covering every recorded
    # count of that rank: rank 3 (dets 1..64, C <= 2) gives c_3 = 1/8, rank 4 (dets 3..20,
    # C <= 2) gives c_4 = 1/4, rank 5 (d = 9, C = 5) gives c_5 = 1.5, rank 6 (d = 28, C = 18)
    # gives c_6 = 3.25. No count is recorded above rank 6; ranks >= 7 use c_6.
    c = {3: 0.125, 4: 0.25, 5: 1.5}.get(i["n"], 3.25)
    return 1 + math.ceil(c * math.sqrt(d))


def t_cell(n, A):
    """One perfect cell: shell enumeration, cone facets, facet orbits, flips.

    V ~ n(n+1)/2 configuration vectors, ~V/2 facet orbits, each flip one
    Fincke-Pohst search in rank n-1 with volume ~ sqrt|A| (LORENTZ_FindPositiveVectors,
    polyhedral_common src_lorentzian/lorentzian_perfect.h:397).
    Doc: perfect-domains-and-lorentzian-cells.md (agent C).
    """
    # One flip (flip_across, lorentzian_cells.py:1809-1848 -> _kernel_flipping) enumerates the
    # first nonempty shell of _PositiveVectorEnumerator (:5387-5600): close vectors of m w in the
    # rank n-1 lattice K = t-perp, multipliers m = 1, 2, ... until nonempty, one PARI CVP call
    # (close_vector_coordinates) per multiplier. The returned shell is the new cell's
    # configuration, ~V vectors; the Fincke-Pohst tree visits ~ n per returned vector. Two
    # multipliers per flip are assumed. A is not a driver of one flip.
    del A
    v = n * (n + 1) / 2
    flip = 2 * T_PARI + n * v * U_FP
    return T_NORMALIZ + (v / 2) * flip + v * v * n * U_ELT


def t_traverse(i, cells=None):
    """Perfect-domain traversal of a Lorentzian lattice with C cell orbits; C^2 transporter tests."""
    c = cells if cells is not None else cells_estimate(i)
    n, A = i["n"], i["A"]
    v = n * (n + 1) / 2
    return c * t_cell(n, A) + c * c * (T_CALL + v * v * U_ELT)


def t_eichler_generators(i, depth=0):
    """Eichler approximate-model generation step (recursive.py:1194-1202; ApproximateModels.h:407-526, 719-755).

    K <= |A_L| covering orbits (Eichler criterion: orbits of primitive vectors of fixed
    norm in a lattice containing U+U are determined by v/div(v) in A_L, Gritsenko-Hulek-
    Sankaran, Compositio 143 (2007) Lemma 3.5); each candidate is a pointed isometry
    reducing to the rank n-1 complement.
    Doc: orthogonal-group-generators.md, "Eichler route".
    """
    n = i["n"]
    k = min(i["A"], 64)
    model = 10 * n ** 3 * U_ELT
    pointed = T_CALL + 10 * n ** 3 * U_ELT
    # The splitting vector is positive in the signed view (positive index min(p,q)),
    # so v-perp has positive index min(p,q)-1 (Lorentzian after one step when
    # min(p,q)=2) and |A_{v-perp}| <= |q(v)| |A_L| = 2|A_L| for a norm-2 vector of U.
    sub = {"n": n - 1, "p": min(i["p"], i["q"]) - 1, "q": max(i["p"], i["q"]),
           "r": 0, "A": 2 * i["A"], "ell": i["ell"] + 1, "even": i["even"]}
    return model + k * pointed + t_group(sub, depth + 1)


def t_group(i, depth=0):
    """Generators of O(L) for the case's branch."""
    match route(i):
        case "degenerate":
            return T_CALL + 10 * i["n"] ** 3 * U_ELT
        case "definite":
            return t_definite(i["n"], *root_proxy(i["n"]))
        case "binary":
            return t_binary(i)
        case "lorentzian":
            return t_traverse(i)
        case "rank4-low-witt":
            return 500 * T_CALL  # quaternion / Clifford order computations in PARI
        case _:
            return t_eichler_generators(i, depth) if depth < 20 else 1.0


def needs_julia(i):
    """Only the centralizer adapter reaches the Oscar bridge (lattice_engines.py:1393)."""
    return False


def t_oq(A, log2_o):
    """Build O(q_L) as a permutation group on |A| points and its order (Brandhorst-Veniani algorithm).

    Doc: discriminant-groups-and-congruence-orbits.md (agent F); orthogonal-group-generators.md "Finite images".
    """
    gens = max(2, math.ceil(log2_o))
    return 10 * T_CALL + A * gens * U_PERM + 1e-7 * A * A + A * U_ELT


def log2_oq_bound(i):
    """log2 |O(q_L)| <= ell * log2|A| (O(q) embeds in Aut(A), |Aut(A)| <= |A|^ell)."""
    return max(1.0, i["ell"] * math.log2(max(i["A"], 2)))


def t_image(i, gens, log2_o=None):
    """Image of `gens` isometries in O(q_L). A discriminant morphism is fixed by the images of the
    ell generators of A_L (ell * n element operations); its permutation on |A| points is |A| U_PERM."""
    lo = log2_o if log2_o is not None else log2_oq_bound(i)
    return t_oq(i["A"], lo) + gens * i["ell"] * i["n"] * U_ELT + gens * i["A"] * U_PERM


def comparison_primes(i, det):
    ps = {3, 5, 7}
    d = abs(2 * det)
    f = 2
    while f * f <= d:
        while d % f == 0:
            ps.add(f)
            d //= f
        f += 1
    if d > 1:
        ps.add(d)
    return sorted(ps)


def t_reduction(n, p, gens):
    """Order of a subgroup of O(n, F_p) by constructive recognition (GAP recog), polynomial in n, log p.

    Doc: orthogonal-group-generators.md, "Finite images" (the orbit-on-vectors method is
    pathological: (p^n-1)/(p-1) points).
    """
    # Measured: the order of SO+(10,2) from its generators. Matrix-group order work scales as
    # n^3 log p per generator (one stabilizer-chain level per basis vector, n^2 log p per image).
    return T_SO10_2 * (n / 10) ** 3 * math.log2(p) * max(gens, 1) / 2


def det_of(gram):
    a = [[Fraction(x) for x in row] for row in gram]
    n = len(a)
    d = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if a[r][c] != 0), None)
        if piv is None:
            return 0
        if piv != c:
            a[c], a[piv] = a[piv], a[c]
            d = -d
        d *= a[c][c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            for k in range(c, n):
                a[r][k] -= f * a[c][k]
    return int(d)


def t_same_finite_images(gram, computed_gens, recorded_gens):
    i = inv(gram)
    det = det_of(gram)
    t = t_image(i, computed_gens) + t_image(i, recorded_gens)
    for p in comparison_primes(i, det):
        t += t_reduction(i["n"], p, computed_gens) + t_reduction(i["n"], p, recorded_gens)
    return t


def budget(model_seconds, julia=False):
    del julia
    return T_NODE + model_seconds


def _fmt(i):
    return f"n={i['n']} sig=({i['p']},{i['q']}) |A|={i['A']}"


# ------------------------------------------------------------ test modules

ORTH = "tests/acceptance/test_orthogonal_group.sage::"
PERF = "tests/acceptance/test_perfect_domains.sage::"
ISO = "tests/acceptance/test_isometry.sage::"
BND = "tests/acceptance/test_boundaries.sage::"
VEC = "tests/acceptance/test_vector_orbits.sage::"
DBL = "tests/acceptance/test_double_cosets.sage::"
CEN = "tests/acceptance/test_centralizers.sage::"
OG = DOC + "orthogonal-group-generators.md"
IT = DOC + "isometry-testing.md"
BC = DOC + "boundary-components.md"
VO = DOC + "vector-orbits.md"
DC = DOC + "double-cosets.md"
CZ = DOC + "centralizers.md"
PD = DOC + "perfect-domains-and-lorentzian-cells.md"


def orthogonal_group():
    out = {}
    m = _load("mertens_generators.json")
    g = m["gram"]
    i = inv(g)
    gens = len(m["omega_generators"]) + 1
    t = t_group(i) + t_same_finite_images(g, 4, gens)
    out[ORTH + "test_mertens_omega_and_minus_one_generate_o_l"] = (budget(t, needs_julia(i)), f"{OG}#branches route={route(i)} {_fmt(i)} + finite images")
    b = _load("binary_form_automorphisms.json")
    for c in b["explicit_generators"]:
        i = inv(c["gram"])
        t = t_group(i) + t_same_finite_images(c["gram"], 2, len(c["automorphism_group_generators"]))
        out[f"{ORTH}test_binary_form_generators_match_hecke[{c['id']}]"] = (budget(t), f"{OG}#binary-forms D={c['discriminant']} + finite images")
    for c in b["improper_automorphisms"]:
        i = inv(c["gram"])
        out[f"{ORTH}test_ambiguous_binary_form_has_an_improper_automorphism[{c['id']}]"] = (budget(t_group(i)), f"{OG}#binary-forms D={c['discriminant']}")
    for c in _load("lorentzian_stabilizers_cases.json"):
        i = inv(c["gram"])
        cells = max(cells_estimate(i), math.ceil(c["num_generators"] / 3))
        t = t_traverse(i, cells) + t_same_finite_images(c["gram"], c["num_generators"], c["num_generators"])
        out[f"{ORTH}test_lorentzian_generators_match_polyhedral_common[{c['id']}]"] = (budget(t), f"{OG}#lorentzian-route cells~{cells} {_fmt(i)} gens={c['num_generators']}")
    for name, n in (("contains_radical_shear_and_gl_action", 2), ("respects_permuted_radical_coordinates", 2),
                    ("generates_rank_two_radical_action", 3), ("generates_permuted_rank_two_radical_action", 3),
                    ("embeds_non_diagonal_complement", 3)):
        # The nonsingular part is <2> (rank one, solved directly) except for A2 in the
        # non-diagonal case, whose automorphism group is PARI qfauto (6 roots, |Aut| = 12).
        bridge = name == "embeds_non_diagonal_complement"
        out[f"{ORTH}test_degenerate_orthogonal_group_{name}"] = (budget(T_CALL + 10 * n ** 3 * U_ELT + (t_definite(2, 6, math.log2(12)) if bridge else 0.0)), f"{OG}#degenerate ambient rank={n}, nonsingular part {'A2 via qfauto' if bridge else '<2>'}")
    out[ORTH + "test_zero_form_orthogonal_group_is_the_full_general_linear_group"] = (budget(T_SMALL), f"{OG}#degenerate zero form: explicit GL_2(Z) generators")
    out[ORTH + "test_zero_form_orthogonal_group_uses_permuted_coordinate_labels"] = (budget(T_SMALL), f"{OG}#degenerate zero form: explicit GL_2(Z) generators")
    out[ORTH + "test_orthogonal_group_of_the_hyperbolic_plane_has_order_four"] = (budget(t_binary({"A": 1}) + T_SMALL), f"{OG}#binary-forms U, finite O(U)")
    oracle = _load("oscar_lattice_oracles.json")["discriminant_images"]
    for c in oracle:
        if "O_qL_order" not in c:
            continue
        i = inv(c["gram"])
        lo = math.log2(c["O_qL_order"])
        t = t_group(i) * 1.5 + 3 * t_image(i, 4, lo)
        out[f"{ORTH}test_discriminant_images_match_oscar[{c['id']}]"] = (budget(t, needs_julia(i)), f"{OG}#finite-images route={route(i)} {_fmt(i)} |O(q)|={c['O_qL_order']}")
    for c in oracle:
        i = inv(c["gram"])
        if c["id"] == "A7_plus_diag_1_1_-1":
            out[ORTH + "test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q"] = (budget(t_group(i) + 2 * t_image(i, 8)), f"{OG}#lorentzian-route {_fmt(i)}")
        if c["id"] == "U_plus_-3":
            out[ORTH + "test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q"] = (budget(1.5 * t_group(i) + 2 * t_image(i, 4)), f"{OG}#lorentzian-route {_fmt(i)}")
    return out


def perfect_domains():
    out = {}
    cases = _load("lorentzian_perfect_domains.json")
    for c in cases:
        i = inv(c["gram"])
        t = t_traverse(i, c["total_count"]) + t_image(i, 2 * c["total_count"])
        out[f"{PERF}test_perfect_domain_orbit_counts_match_polyhedral_common[{c['id']}]"] = (budget(t), f"{PD}#traversal total mode C={c['total_count']} {_fmt(i)} + generation check")
        t = t_traverse(i, c["isotropic_count"])
        out[f"{PERF}test_isotropic_mode_counts_match_polyhedral_common[{c['id']}]"] = (budget(t), f"{PD}#traversal isotropic mode C={c['isotropic_count']} {_fmt(i)}")
    r1 = next(c for c in cases if c["id"] == "Result_Enumeration_1")
    i1 = inv(r1["gram"])
    trav = t_traverse(i1, r1["total_count"])
    out[PERF + "test_lazy_neighbor_cones_retain_exact_transformed_rays"] = (budget(2 * trav), f"{PD}#traversal C=2 n=3 traversal + adjacency materialization")
    out[PERF + "test_traversal_retains_only_adjacency_facet_orbits"] = (budget(trav), f"{PD}#traversal C=2 n=3 traversal")
    out[PERF + "test_facet_orbit_permutation_view_materializes_exact_isometries"] = (budget(t_cell(3, 8)), f"{PD}#traversal one cell n=3 facet orbits")
    rows = 4000 * 6
    out[PERF + "test_large_cell_cone_uses_exact_normaliz_h_presentation"] = (budget(t_cell(3, 8) + 10 * T_CALL + rows * 3 * 10 * U_ELT), f"{PD}#traversal Normaliz on {rows} rows, incidence check rows*facets*n")
    out[PERF + "test_large_pairing_graph_fallback_returns_exact_transporter"] = (budget(T_SMALL + 257 * 257 * U_ELT), f"{PD}#traversal 257-row transporter, n=2")
    for name in ("non_diagonal_positive_shell_uses_the_pinned_cvp_center", "non_diagonal_positive_shell_respects_explicit_coordinate_labels",
                 "positive_shell_survives_a_permuted_gram_basis", "non_diagonal_isotropic_shell_uses_the_exact_sphere_center",
                 "positive_shell_scales_rational_direction_pairing_bound", "non_diagonal_isotropic_shell_respects_explicit_coordinate_labels",
                 "isotropic_shell_survives_a_permuted_gram_basis"):
        out[f"{PERF}test_{name}"] = (budget(T_SMALL), f"{PD}#traversal single first-shell CVP in rank 3")
    return out


def t_isometry(g1, g2):
    """is_isometric + Isom().an_element() on the case's branch. Doc: isometry-testing.md."""
    i = inv(g1)
    genus = T_CALL * i["n"] * (1 + i["ell"])
    match route(i):
        case "lorentzian":
            c = cells_estimate(i)
            return genus + t_traverse(i, c) + t_cell(i["n"], i["A"]) + c * (T_CALL + i["n"] ** 4 * U_ELT), f"{IT}#lorentzian C~{c} {_fmt(i)}", False
        case "eichler":
            k = min(i["A"], 64)
            return genus + 10 * T_CALL * i["ell"] + 10 * i["n"] ** 3 * U_ELT + k * 2 * (T_CALL + 10 * i["n"] ** 3 * U_ELT), f"{IT}#eichler-covering K<={k} {_fmt(i)}", False
        case "binary":
            return genus + 2 * t_binary(i), f"{IT}#binary {_fmt(i)}", False
        case "rank4-low-witt":
            return genus + 500 * T_CALL, f"{IT}#low-witt {_fmt(i)}", False
        case "definite":
            return genus + t_definite(i["n"], *root_proxy(i["n"])), f"{IT}#definite qfisom {_fmt(i)}", False
    return genus + 10 * T_CALL, f"{IT}#degenerate {_fmt(i)}", False


def isometry():
    out = {}
    for c in _load("indefinite_isometry_pairs.json"):
        t, basis, j = t_isometry(c["gram1"], c["gram2"])
        out[f"{ISO}test_isometry_verdicts_match_certified_pairs[{c['id']}]"] = (budget(t, j), basis)
    # Only 2 distinct mat1 lattices occur among the 145 cases. The source traversal is memoized
    # per reduced lattice (recursive.py groupoid cache), so it is charged to the first case of each
    # mat1 in collection (fixture) order; every case pays the target reduction, its initial cell,
    # and C cell-transporter tests.
    seen = set()
    for c in _load("lorentzian_equivalence_145.json"):
        i = inv(c["mat1"])
        key = json.dumps(c["mat1"])
        cells = cells_estimate(i)
        genus = T_CALL * i["n"] * (1 + i["ell"])
        t = 2 * genus + t_cell(i["n"], i["A"]) + cells * (T_CALL + i["n"] ** 4 * U_ELT)
        first = key not in seen
        if first:
            t += t_traverse(i, cells)
            seen.add(key)
        out[f"{ISO}test_lorentzian_pairs_are_isometric[{c['id']}]"] = (
            budget(t), f"{IT}#lorentzian C~{cells} {_fmt(i)}" + (" +source traversal (first of its mat1)" if first else " source traversal shared"))
    for c in _load("ci_indefinite_comp.json"):
        t, basis, j = t_isometry(c["gram"], c["gram"])
        out[f"{ISO}test_change_of_basis_gives_an_isometric_lattice[{c['id']}]"] = (budget(t, j), basis)
    cs = _load("conway_sloane_cases.json")["spinor_genus_pair_determinant_minus_128"]
    i = inv(cs["form_a"]["gram"])
    out[ISO + "test_conway_sloane_same_genus_pair_is_not_isometric"] = (budget(3 * (T_CALL * i["n"] * (1 + i["ell"])) + 10 * T_CALL), f"{IT}#spinor-genus decision {_fmt(i)}, two genus symbols")
    for name, t, why in (
        ("row_matrix_morphism_roundtrip_uses_coordinate_labels", T_SMALL, "rank-2 morphism construction"),
        ("column_matrix_morphism_roundtrip_uses_coordinate_labels", 2 * T_SMALL, "rank-2 morphism constructions"),
        ("equal_gram_groupoid_identification_preserves_coordinate_labels", T_SMALL, "coordinate identification"),
        ("degenerate_isometry_witness_respects_complement_coordinates", 10 * T_CALL + t_definite(2, 6, math.log2(12)), "degenerate rank 3, A2 kernel by qfisom"),
        ("eichler_generator_families_preserve_nonlexical_complement", 10 * 6 ** 3 * U_ELT + 6 ** 4 * 30 * U_ELT, "Eichler model U+U+<2>+<4> n=6, generator checks n^2 per generator"),
        ("two_u_full_generating_data_accepts_raw_sage_integer_norm", 10 * 5 ** 3 * U_ELT + T_SMALL, "Eichler model U+U+<2> n=5, one covering family"),
        ("odd_eichler_envelope_uses_declared_coordinate_order", 4 * T_SMALL, "four odd envelopes n<=6"),
        ("large_cell_transporter_invariant_rejection_returns_none", 2 * 513 * 2 * U_ELT * 10 + T_SMALL, "513-row invariant rejection, n=2"),
    ):
        out[ISO + "test_" + name] = (budget(t), f"{IT}#unit {why}")
    return out


# ------------------------------------------------------------ boundaries


def _cartan(kind, m):
    if kind == "D":
        c = [[2 if i == j else (-1 if abs(i - j) == 1 and max(i, j) < m - 1 else 0) for j in range(m)] for i in range(m)]
        c[m - 3][m - 1] = c[m - 1][m - 3] = -1
        return c
    e8 = [[2, 0, -1, 0, 0, 0, 0, 0], [0, 2, 0, -1, 0, 0, 0, 0], [-1, 0, 2, -1, 0, 0, 0, 0], [0, -1, -1, 2, -1, 0, 0, 0],
          [0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0], [0, 0, 0, 0, 0, -1, 2, -1], [0, 0, 0, 0, 0, 0, -1, 2]]
    return e8


def _block(*blocks):
    n = sum(len(b) for b in blocks)
    g = [[0] * n for _ in range(n)]
    o = 0
    for b in blocks:
        for r, row in enumerate(b):
            for c, x in enumerate(row):
                g[o + r][o + c] = x
        o += len(b)
    return g


def _neg(b):
    return [[-x for x in row] for row in b]


U = [[0, 1], [1, 0]]




def word_length(degree, log2_order):
    """Length w of a GAP word-section word (orthogonal_quotients.py:176-191, PreImagesRepresentative
    of a free-group epimorphism through the stabilizer chain). Base length b = ceil(log2|Q| / log2 deg);
    each level contributes at most one Schreier-tree path of depth <= 2 ceil(log2 deg) (GAP's shallow
    trees), and each strong generator is a word of length <= 2.5 in the family on average."""
    b = math.ceil(log2_order / math.log2(degree))
    return math.ceil(b * 2 * math.ceil(math.log2(degree)) * 2.5)


def t_lift(i, degree, log2_order):
    """One section lift (splitting_isometries_from_stabilizer_image, orthogonal_quotients.py:1600-1630):
    PreImagesRepresentative in GAP, O(w + degree); _evaluate_letters, w products of n x n isometries;
    _checked_lift, one character image (ell discriminant generators x n coordinates); one apex action."""
    w = word_length(degree, log2_order)
    return T_CALL + (w + degree) * U_PERM + w * U_PRODUCT + i["ell"] * i["n"] * U_ELT + 2 * i["n"] * U_ELT


def t_split(i, count, orbit_points, gens=12, degree=None, log2_order=None):
    """Split O(L)-classes into Gamma-orbits (split_orbit, lattice_morphisms.py:4638; DoubleCosets in the
    permutation model, orthogonal_quotients.py:1600-1630): GAP double cosets ~ orbit_points * gens
    permutation images, then exactly one section lift per double coset = per output representative
    (`count`). Doc: boundary-components.md "Splitting"."""
    deg = degree if degree is not None else max(i["A"] - 1, 2)
    lo = log2_order if log2_order is not None else log2_oq_bound(i)
    return orbit_points * gens * U_PERM + count * t_lift(i, deg, lo)


def t_classes(i, k, complement=None):
    """O(L)-classes of primitive isotropic rank-k sublattices / flags (isotropic_induction.py:471) and
    their stabilizers. Lines: stabilizer through v-perp/v (rank n-2, Lorentzian for signature (2,n)).
    Planes: classes of E-perp/E, a definite rank n-4 genus (Kneser neighbours, class number h).
    Doc: boundary-components.md "Isotropic classes".
    """
    n = i["n"]
    if k == 1:
        sub = {"n": n - 2, "p": 1, "q": n - 3, "r": 0, "A": i["A"], "ell": i["ell"], "even": i["even"]}
        return T_CALL + t_traverse(sub)
    m = n - 4
    h = 4
    short, log2_aut = complement if complement is not None else root_proxy(m)
    return h * h * t_definite(m, short, log2_aut) + 20 * h * (T_CALL + m * m * U_ELT * 10) + h * T_CALL


def t_counts(i, counts, orbit_points):
    lines, planes, flags = counts
    t = t_classes(i, 1) + t_classes(i, 2) + t_classes(i, 2) + 10 * T_CALL
    t += t_split(i, lines, orbit_points[0]) + t_split(i, planes, orbit_points[1]) + t_split(i, flags, orbit_points[2])
    return t


def t_splits_enriques(i, counts, orbit_points):
    """Per-case Enriques splits: GAP double cosets on the shared configuration sets of A_N, one
    0.042 s section lift per output representative."""
    lo = math.log2(2 ** 21 * 3 ** 5 * 5 ** 2 * 7 * 17 * 31)
    return sum(t_split(i, c, p, degree=527, log2_order=lo) for c, p in zip(counts, orbit_points))


def isotropic_points_bound(A, k):
    """|G/G_m| <= number of k-step isotropic configurations in A_L: <= A^k."""
    return min(A ** k, 10 ** 7)


def boundaries():
    out = {}
    d = _load("dawes_buildings.json")
    groups = {"dawes_2U_A2_stable": 1, "dawes_2U2_A2_stable": 1, "dawes_2U_minus2_minus6": 2, "dawes_L2_gamma2": 1}
    for c in d["buildings"]:
        i = inv(c["gram"])
        b = c["building"]
        pts = tuple(isotropic_points_bound(i["A"], k) for k in (1, 2, 2))
        t = groups[c["id"]] * (t_counts(i, (b["points"], b["curves"], b["edges"]), pts) + t_oq(i["A"], log2_oq_bound(i)))
        out[f"{BND}test_dawes_buildings_have_the_published_counts[{c['id']}]"] = (budget(t), f"{BC}#dawes groups={groups[c['id']]} {_fmt(i)} counts={b['points']},{b['curves']},{b['edges']}")
    for gph in d["incidence_graphs"]:
        c = next(x for x in d["buildings"] if x["id"] == gph["building_id"])
        i = inv(c["gram"])
        b = c["building"]
        pts = tuple(isotropic_points_bound(i["A"], k) for k in (1, 2, 2))
        t = t_counts(i, (b["points"], b["curves"], b["edges"]), pts) + T_CALL
        out[f"{BND}test_dawes_building_incidence_graph_is_the_published_figure[{gph['building_id']}]"] = (budget(t), f"{BC}#dawes cusps+incidence {_fmt(i)}")
    c = next(x for x in d["buildings"] if x["id"] == "dawes_2U_A2_stable")
    i = inv(c["gram"])
    out[BND + "test_dawes_last_index_in_the_chain_is_the_discriminant_image_of_o_plus"] = (budget(1.5 * t_group(i) + t_image(i, 8)), f"{BC}#dawes O+ generators {_fmt(i)}")
    dt = _load("dtower_boundaries.json")

    def lam(N):
        m = N - 2
        dd = [[-4]] if m == 1 else ([[-2, 0], [0, -2]] if m == 2 else _neg(_cartan("D", m)))
        return _block(U, U, dd)

    for c in dt["type_counts"]:
        i = inv(lam(c["N"]))
        pts = (isotropic_points_bound(i["A"], 1), isotropic_points_bound(i["A"], 2), 0)
        t = t_classes(i, 1) + t_classes(i, 2) + t_split(i, c["type_iii_components"], pts[0]) + t_split(i, c["type_ii_components"], pts[1]) + t_oq(i["A"], log2_oq_bound(i)) + t_group(i)
        out[f"{BND}test_dtower_boundary_component_counts[N_{c['N']}]"] = (budget(t), f"{BC}#dtower Gamma_xi {_fmt(i)} III={c['type_iii_components']} II={c['type_ii_components']}")
    for p in dt["boundary_pictures"]:
        i = inv(p["gram"])
        pts = tuple(isotropic_points_bound(i["A"], k) for k in (1, 2, 2))
        t = t_counts(i, (p["type_iii_points"], p["type_ii_curves"], len(p["incidences"])), pts) + t_group(i)
        out[f"{BND}test_dtower_boundary_pictures[{p['id']}]"] = (budget(t), f"{BC}#dtower {_fmt(i)}")
    e8n = _neg(_cartan("E", 8))
    k3 = _load("k3_modular_strata.json")
    for k, name, curves in ((1, "test_degree_two_k3_type_ii_components_and_root_types", k3["degree_two_polarized_k3"]["baily_borel_boundary"]["type_ii_curves"]),
                            (2, "test_degree_four_k3_has_nine_type_ii_components", k3["degree_four_polarized_k3"]["baily_borel_boundary"]["type_ii_curves"])):
        i = inv(_block(U, U, e8n, e8n, [[-2 * k]]))
        h = curves
        m = i["n"] - 4
        # E-perp/E = <-2k> + 2E8(-1), priced on its decomposition (see e8_plus_line):
        # 2*240 + 2 short vectors, |Aut| = 2 |W(E8)|^2 * 2.
        short = 482
        t = h * h * t_definite(m, short, 2 + 2 * LOG2_W_E8) + 20 * h * T_CALL + h * T_CALL + 10 * T_CALL
        out[BND + name] = (budget(t), f"{BC}#k3-type-ii genus of E-perp/E rank {m}, h={h}, {_fmt(i)}")
    # Enriques: N = U + U(2) + E8(-2), sig (2,10), A_N = (Z/2)^10, O(q_N) = O^+(10,2) x? (order recorded).
    n_gram = _block(U, [[0, 2], [2, 0]], [[2 * x for x in row] for row in e8n])
    i = inv(n_gram)
    order_o_plus = 2 ** 21 * 3 ** 5 * 5 ** 2 * 7 * 17 * 31
    lo = math.log2(order_o_plus)
    # One-time setup cached by _enriques_transport: primitive extension of the
    # Enriques involution and an isometry E10(2) -> M (Lorentzian rank 10, |A|=2^10).
    m_i = {"n": 10, "p": 1, "q": 9, "r": 0, "A": 1024, "ell": 10, "even": True}
    setup = 500 * T_CALL + t_traverse(m_i, 1) + t_cell(10, 1024)
    # Shared per session (charged to case_1, first in collection order): the transport setup,
    # O(N) generators, the O(N)-classes of lines, planes and flags, and the permutation model.
    shared = setup + t_group(i) + t_classes(i, 1) + 2 * t_classes(i, 2, (240, LOG2_W_E8)) + 10 * T_CALL + t_oq(i["A"], lo)
    # Per case: Stab(h) in O(M) through the definite rank-9 complement, its image in O(q_M),
    # the glue transport and preimage generators (one section lift per image generator, <= 20),
    # and the image of Gamma_h in O(q_N) from <= 40 generators.
    gamma_h = t_definite(9, *e8_plus_line(1)) + t_image(m_i, 10, lo) + 20 * t_lift(i, 527, lo) + t_image(i, 40, lo)
    rows = _load("enriques_87_polarizations.json")
    # Orbit sizes |O(L)/O(L)_m| on A_N: singular points of O^+(10,2): (2^5-1)(2^4+1)=527;
    # totally singular lines 527*135/3=23715; point-line flags 527*135=71145.
    pts = (527, 23715, 71145)
    for idx, c in enumerate(rows):
        counts = (c["line_orbits"], c["plane_orbits"], c["flag_orbits"])
        t = gamma_h + t_splits_enriques(i, counts, pts) + (shared if idx == 0 else 0.0)
        out[f"{BND}test_enriques_modular_groups_match_dutour_sikiric_hulek[case_{c['case']}]"] = (
            budget(t), f"{BC}#enriques N sig(2,10) |A|=1024 |Gamma image|={c['group_order']} orbits={counts}" + (" +shared setup" if idx == 0 else ""))
    ci = _load("enriques_classical_indices.json")
    for c in ci["consistent_with_table"] + ci["inconsistent_with_table"]:
        out[f"{BND}test_classical_polarization_indices[degree_{c['degree']}]"] = (budget(gamma_h + t_image(i, 40, lo)), f"{BC}#enriques Gamma_h image only, case {c['case']}")
    return out


# ------------------------------------------------------------ vector orbits


def t_lorentz_vector_classes(i, square, orbits):
    """Positive-cone vector classes by MarkedCellOrbitAlgorithm over the perfect complex
    (marked_lattices.py:600-606). Enumeration volume of norm-s vectors in a cell with bounded pairing:
    N(s) ~ 10 * (s/2)^{(n-1)/2} points (Fincke-Pohst ellipsoid volume, Cohen GTM 138 Alg. 2.7.5).
    Doc: vector-orbits.md "Lorentzian cone".
    """
    n = i["n"]
    pts = 10 * (max(square, 2) / 2) ** ((n - 1) / 2)
    return t_traverse(i, 1) + pts * U_FP + orbits * (T_CALL + n ** 3 * U_ELT + t_definite(n - 1, *e8_plus_line(max(square, 2) // 2)) * 0.2)


def t_embedding_classes(i):
    """Negative squares: native Nikulin primitive-embedding classes (lattice_morphisms.py:2670).
    In a Lorentzian L of signature (1, n-1) the complement of a negative vector is Lorentzian of
    rank n-1, an indefinite genus of one class per spinor genus; the cost is its genus symbol and
    the glue double cosets in O(q) of order <= 2|s|. Doc: vector-orbits.md "Embedding branch"."""
    m = i["n"] - 1
    return 100 * T_CALL + m * T_CALL


def vector_orbits():
    out = {}
    e10 = _load("e10_fundamental_domain.json")
    i = inv(e10["gram"])
    for c in _load("enriques_polarization_orbits.json"):
        s = c["two_d"]
        out[f"{VEC}test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_{s}]"] = (
            budget(t_lorentz_vector_classes(i, s, c["primitive_vector_orbits"])), f"{VO}#lorentzian-cone E10 s={s} orbits={c['primitive_vector_orbits']}")
    reps = e10["orbit_representatives"]
    for deg in sorted({r["degree"] for r in reps}):
        k = sum(r["degree"] == deg for r in reps)
        t = t_lorentz_vector_classes(i, deg, k) + k * k * (t_definite(9, *e8_plus_line(deg // 2)) + T_CALL)
        out[f"{VEC}test_published_representatives_are_one_per_orbit[{deg}]"] = (budget(t), f"{VO}#lorentzian-cone E10 s={deg} k={k}, k^2 pointed isometry tests")
    orbits = _load("e10_vector_orbits.json")["orbits"]
    for s in sorted({o["h_squared"] for o in orbits}):
        k = sum(o["h_squared"] == s for o in orbits)
        if s == 0:
            t = t_traverse(i, 1) + 10 * T_CALL * k
            basis = f"{VO}#isotropic E10 null classes k={k}"
        elif s < 0:
            t = t_embedding_classes(i)
            basis = f"{VO}#embedding-branch E10 s={s} k={k}"
        else:
            t = t_lorentz_vector_classes(i, s, k)
            basis = f"{VO}#lorentzian-cone E10 s={s} k={k}"
        out[f"{VEC}test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[{s}]"] = (budget(t), basis)
    a = _load("allcock_i_2_10_orbits.json")
    ia = inv(a["gram"])
    out[VEC + "test_i_2_10_is_transitive_on_norm_minus_one_vectors"] = (budget(t_eichler_generators(ia) + 100 * T_CALL), f"{VO}#eichler odd unimodular {_fmt(ia)}")
    out[VEC + "test_anisotropic_binary_vector_classes_use_primitive_embeddings"] = (budget(4 * t_binary({"A": 6}) + 8 * T_CALL), f"{VO}#binary anisotropic disc 24")
    e8n = _neg(_cartan("E", 8))
    ik = inv(_block(U, U, U, e8n, e8n))
    for s in _load("k3_modular_strata.json")["k3_unimodular_lattice"]["tested_represented_norms"]:
        t = 10 * 22 ** 3 * U_ELT + T_CALL + 100 * T_CALL
        out[f"{VEC}test_k3_lattice_has_one_primitive_orbit_per_norm[{s}]"] = (budget(t), f"{VO}#eichler unimodular {_fmt(ik)} K=1 s={s}")
    return out


def double_cosets():
    out = {}
    cases = _load("double_coset_cases.json")
    for c in cases:
        deg = c["degree"]
        index = c["big_group_order"] // c["small_group_order"]
        gg, gh = len(c["group_g"]["generators"]), len(c["group_h"]["generators"])
        ss = 2 * (T_CALL + deg * deg * (gg + gh) * U_PERM * 50)
        orbit = index * (gh + gg) * (U_PERM * deg + U_ELT * 4)
        t = T_CALL + ss + orbit + index * U_ELT * 20
        out[f"{DBL}test_double_cosets_partition_the_rational_group[{c['id']}]"] = (budget(t), f"{DC} degree={deg} |G|={c['big_group_order']} |H|={c['small_group_order']} [G:H]={index}")
    out[DBL + "test_the_recorded_case_families_are_all_present"] = (budget(0.0), f"{DC} fixture assertions only")
    return out


def centralizers():
    out = {}
    for c in _load("isometry_centralizers.json")["centralizer_cases"]:
        i = inv(c["gram"])
        o = c.get("centralizer_image_order", 2 ** log2_oq_bound(i))
        t = 100 * T_CALL + 10 * T_CALL * i["n"] ** 2 + t_oq(i["A"], math.log2(max(o, 2))) + c["isometry_order"] * i["n"] ** 3 * U_ELT
        out[f"{CEN}test_centralizer_images_match_oscar[{c['id']}]"] = (budget(t, julia=True), f"{CZ} Oscar hermitian Miranda-Morrison {_fmt(i)} ord(f)={c['isometry_order']} image={o}")
    return out


# ------------------------------------------------------------ added scope (missing_B.txt)

MISSING = ROOT / "ci" / "test_budget_models" / "missing_B.txt"
VS = "tests/acceptance/test_vector_stabilizers.sage::"
HILBERT = {5: (1, 5), 2: (1, 8), 10: (2, 40)}  # d: (class number h(Q(sqrt d)), |disc K|)


def t_hilbert(d):
    """Witt-index-one route U + (O_K, Tr) for K = Q(sqrt d): Hilbert cusps = h(K); reduction bound
    C^2 = |delta| Disc/2 gives O(h sqrt(Disc)) candidates per norm (witt_one_fixed_norm_classes.py).
    Doc: vector-orbits.md "Witt index one"."""
    h, disc = HILBERT[d]
    return 100 * T_CALL + 100 * T_CALL * h * math.sqrt(disc)


def t_small_isotropic(rank, k=1):
    """Isotropic induction on rank <= 4 (isotropic_induction.py:744-750): e-perp/e has rank rank-2k,
    binary or zero; Witt extension and containing bound are O(rank^3) element operations.
    Doc: orthogonal-group-generators.md "Isotropic stabilizers"."""
    return 10 * T_CALL + k * (t_binary({"A": 16}) + rank ** 3 * 200 * U_ELT)


def vector_stabilizers():
    out = {}
    e10 = _load("e10_fundamental_domain.json")
    i = inv(e10["gram"])
    orbits = _load("e10_vector_orbits.json")["orbits"]
    for s in sorted({o["h_squared"] for o in orbits}):
        k = sum(o["h_squared"] == s for o in orbits)
        if s > 0:
            classes, stab = t_lorentz_vector_classes(i, s, k), t_definite(9, *e8_plus_line(s // 2))
        elif s < 0:
            classes = t_embedding_classes(i)
            stab = t_traverse({"n": 9, "p": 1, "q": 8, "r": 0, "A": abs(s), "ell": 1, "even": True})
        else:
            classes, stab = t_traverse(i, 1) + 10 * T_CALL * k, t_definite(8, 240, LOG2_W_E8)
        mod2 = 10 * T_CALL + 20 * 1023 * U_PERM
        out[f"{VS}test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[{s}]"] = (
            budget(classes + k * (stab + mod2)), f"{VO}#stabilizers E10 s={s} k={k}, v-perp {'definite' if s > 0 else 'Lorentzian' if s < 0 else 'E8 via e-perp/e'}, mod-2 image on 1023 vectors")
    reps = len(e10["orbit_representatives"])
    for root in e10["simple_roots"]:
        out[f"{VS}test_equivalence_witness_carries_a_representative_to_its_reflection[r_{root['label']}]"] = (
            budget(sum(t_definite(9, *e8_plus_line(r["degree"] // 2)) + 10 * T_CALL for r in e10["orbit_representatives"])), f"{IT}#pointed-isometries {reps} pointed isometries in E10, qfisom on E8+<2d> complements")
    small = {
        "isotropic_pointed_automorphism_generators_preserve_marking_and_form": (t_small_isotropic(4), "2U line and pointed stabilizers"),
        "divisibility_two_isotropic_pointed_isometry_has_exact_witness": (t_small_isotropic(4), "U(2)+U(2) div 2"),
        "odd_lattice_isotropic_pointed_isometry_uses_general_divisibility_two_route": (t_small_isotropic(2), "odd binary div 2"),
        "general_isotropic_pointed_isometry_carries_negative_marking": (t_small_isotropic(4), "U(2)+U(2)"),
        "isotropic_pointed_isometry_transports_between_changed_basis_parents": (t_small_isotropic(4), "U(2)+U(2), two parents"),
        "split_divisibility_one_pointed_isometry_carries_actual_negative_marking": (t_small_isotropic(2), "U"),
        "isotropic_pointed_isometry_rejects_different_divisibility": (T_SMALL, "divisibility invariant rejects"),
        "isotropic_plane_isometry_carries_changed_basis_plane": (t_small_isotropic(4, 2), "U(4)+U(4) plane"),
        "isotropic_plane_isometry_transports_between_changed_basis_parents": (t_small_isotropic(4, 2), "U(4)+U(4) plane"),
        "isotropic_flag_isometry_preserves_line_inside_plane": (t_small_isotropic(4, 3), "U(4)+U(4) flag"),
        "isotropic_flag_isometry_transports_between_changed_basis_parents": (t_small_isotropic(4, 3), "U(4)+U(4) flag"),
        "isotropic_plane_and_flag_automorphism_generators_preserve_their_marks": (t_small_isotropic(4, 3) + T_SMALL, "plane/flag stabilizers, F_3 images 48 and 12"),
        "divisibility_four_isotropic_stabilizer_constructs_exact_containing_bound": (t_small_isotropic(3), "U(4)+<-8> cokernel 4"),
        "divisibility_four_isotropic_plane_stabilizer_constructs_exact_containing_bound": (t_small_isotropic(4, 2) * 2, "U(4)+U(4) two planes, cokernel 16"),
        "vector_stabilizer_has_oscars_order": (t_binary({"A": 36}) + T_SMALL, "rank-3 Lorentzian pointed, v-perp binary; finite order 2"),
    }
    for name, (t, why) in small.items():
        out[f"{VS}test_{name}"] = (budget(t), f"{OG}#isotropic-stabilizers {why}")
    return out


def other_scope(keys):
    """Remaining ids from missing_B.txt, classified by file."""
    out = {}
    for key in keys:
        if key in out or key.startswith(VS):
            continue
        name = key.split("::")[1]
        if "witt_one_fixed_norm" in key:
            d = 10 if "sqrt_ten" in name else 2 if ("sqrt_two" in name or "u_plus_two_minus_four" in name) else 5
            t = t_hilbert(d) * (3 if key.startswith("witnesses/") else 1) + (t_hilbert(10) if "recursion_uses" in name else 0.0)
            out[key] = (budget(t), f"{VO}#witt-index-one K=Q(sqrt {d}) h={HILBERT[d][0]} Disc={HILBERT[d][1]}")
        elif "test_fixed_norm_classes" in key:
            out[key] = (budget(100 * T_CALL + 20 * T_CALL * (2 if "norm_two" in name or "eight" in name else 1)), f"{VO}#quaternion-norm-lattices (-1,3), normalizers in box 3")
        elif "pairing_slices" in key:
            t = 6 * T_SMALL if "centralizer" in name else 4 * T_SMALL
            out[key] = (budget(t), f"{IT}#pairing-slices rank<=4 slices, finite pairing shells")
        elif "eichler_covering_classes" in key:
            i = inv(_block(U, U, [[-2 * int(a == b) for b in range(4)] for a in range(4)]) if "four_negative" in name else _block(U, U))
            k = min(i["A"], 64)
            out[key] = (budget(10 * i["n"] ** 3 * U_ELT + k * (T_CALL + 10 * i["n"] ** 3 * U_ELT)), f"{OG}#eichler-route covering K<={k} {_fmt(i)}")
        elif "spinor_genus_negative" in key or "negative_isometry_lorentzian" in key:
            out[key] = (budget(2 * (T_CALL * 4 * 3) + 10 * T_CALL), f"{IT}#decision spinor-genus negative, rank<=4")
        elif "canonical_hyperbolic" in key:
            out[key] = (budget(2 * t_binary({"A": 1}) + T_SMALL), f"{IT}#binary U canonicalization and correction transporter")
        elif "gram_memoization" in key:
            out[key] = (budget(0.0), "memoized Gram access, no search")
        elif "prescribed_gram" in key:
            out[key] = (budget(2 * (10 * 4 ** 3 * U_ELT + T_SMALL)), f"{VO}#eichler-branch embedding orbits in U+U, box 1")
        else:
            raise AssertionError(f"unclassified {key}")
    return out


def missing_nodes():
    with open(MISSING, encoding="utf-8") as handle:
        return {line for line in handle.read().splitlines() if line.strip()}


def pending_constants():
    """Constants requested by agents A and F for the Enriques lattice N = U + U(2) + E8(-2)."""
    e8n = _neg(_cartan("E", 8))
    i = inv(_block(U, [[0, 2], [2, 0]], [[2 * x for x in row] for row in e8n]))
    return {
        "line_class_seconds": t_classes(i, 1),
        "plane_class_seconds": t_classes(i, 2),
        "split_orbit_lift_seconds": t_lift(i, 527, math.log2(2 ** 21 * 3 ** 5 * 5 ** 2 * 7 * 17 * 31)),
        "word_length": word_length(527, math.log2(2 ** 21 * 3 ** 5 * 5 ** 2 * 7 * 17 * 31)),
    }


def estimates() -> dict[str, tuple[float, str]]:
    result = {}
    for part in (orthogonal_group, perfect_domains, isometry, boundaries, vector_orbits, double_cosets, centralizers, vector_stabilizers):
        result.update(part())
    result.update(other_scope(sorted(missing_nodes())))
    return dict(result)


def scope_nodes():
    with open(NODES, encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    return {
        line for line in lines
        if line.startswith("tests/acceptance/")
        and "test_reflective.sage" not in line
        and "test_isotropic_orbits.sage" not in line
    }


if __name__ == "__main__":
    table = estimates()
    scope = scope_nodes() | missing_nodes()
    missing, extra = scope - table.keys(), table.keys() - scope
    print("count", len(table), "scope", len(scope))
    for key in sorted(missing):
        print("MISSING", key)
    for key in sorted(extra):
        print("EXTRA", key)
    functions = {}
    for key, (seconds, _basis) in table.items():
        functions.setdefault(key.split("[")[0], []).append(seconds)
    for name, values in sorted(functions.items()):
        print(f"{name.split('tests/', 1)[-1]}  n={len(values)} min={min(values):.1f} median={median(values):.1f} max={max(values):.1f}")
    print("pending constants (seconds):", {k: round(v, 2) for k, v in pending_constants().items()})
    assert not missing and not extra
