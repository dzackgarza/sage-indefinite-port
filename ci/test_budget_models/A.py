"""A-priori time budgets for tests/acceptance/test_reflective.sage and test_isotropic_orbits.sage.

Derived by inspection only; no measured duration is used. Each case is first classified
into the branch the production code takes, then the work of that branch is estimated from
the mathematical invariants that drive it.

Anchor: polyhedral_common@1592b246 .github/workflows/ci_20_reflective.yml runs, in one job
with timeout-minutes 120 that also installs GAP and builds the C++, the edgewalk
reflectivity test (simple roots, isotropy) on all 8821 ListReflect lattices
(CI_tests/20_Reflective/TestReflectivity.g FullTest). ci_09A_isotropic.yml decides isotropy
and finds a vector for all 8821 forms under the same 120 min limit. So the C++ route needs
< 120*60/8821 = 0.82 s per lattice on average. Port factor 10: the heavy primitives
(Smith form, qfauto/qfisom, short vectors, LLL, dual description) are the same PARI/FLINT/
cdd C code; only the walk bookkeeping is interpreted Python, about 10x the C++ loop.
Family consistency: the sum of the ListReflect budgets must stay below 8821*0.82*10 s.

Unit costs, derived from the operation at the fixtures' sizes (rank <= 18, entries < 2^25),
in one shared Sage session (no per-case start-up):
  OP   = 2e-5 s   one exact small linear-algebra or element operation: a FLINT/PARI call on
                  a rank <= 18 integer vector/matrix is ~1-5 us; the Python call and parent
                  coercion around it is ~10-20 us.
  STAB = 5e-3 s   stabilizer/transporter of a coloured family of <= n+1 vectors at rank 3:
                  PARI qfauto/qfisom on a 3-dimensional definite form with a handful of
                  short vectors is < 1 ms (Plesken-Souvignier backtrack over O(10) vectors)
                  plus wrapping; scaled by n^2/9 and log2|O(q_L)| above rank 3.
  LP   = 1e-6 s   one halfspace-ray incidence in cdd's double description (C arithmetic).
  BASE = 0.05 s   lattice construction, Smith form and root norms, isotropy decision and
                  qfsolve, in the shared session.
"""

from __future__ import annotations

import json
import math
import itertools
import statistics
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures"
NODES = ROOT / "ci" / "test_budget_models" / "nodes_initial.txt"
REFL = "tests/acceptance/test_reflective.sage::"
ISO = "tests/acceptance/test_isotropic_orbits.sage::"

OP, STAB, LP, BASE = 2e-5, 5e-3, 1e-6, 0.05


def _load(name):
    with open(FIX / f"{name}.json", encoding="utf-8") as handle:
        return json.load(handle)


def smith_invariants(gram: list[list[int]]) -> list[int]:
    """Elementary divisors of the Gram matrix = invariant factors of A_L = L^v/L.

    Reference: Cohen, A Course in Computational Algebraic Number Theory, Alg. 2.4.14
    (Smith normal form by elementary row/column operations).
    """
    a = [row[:] for row in gram]
    n = len(a)
    out = []
    for t in range(n):
        while True:
            entries = [(abs(a[i][j]), i, j) for i in range(t, n) for j in range(t, n) if a[i][j]]
            if not entries:
                return out + [0] * (n - t)
            _, i, j = min(entries)
            a[t], a[i] = a[i], a[t]
            for row in a:
                row[t], row[j] = row[j], row[t]
            p = a[t][t]
            for i in range(t + 1, n):
                q = a[i][t] // p
                for j in range(t, n):
                    a[i][j] -= q * a[t][j]
            for j in range(t + 1, n):
                q = a[t][j] // p
                for i in range(t, n):
                    a[i][j] -= q * a[i][t]
            if any(a[i][t] for i in range(t + 1, n)) or any(a[t][j] for j in range(t + 1, n)):
                continue
            bad = next((i for i in range(t + 1, n) if any(a[i][j] % p for j in range(t + 1, n))), None)
            if bad is None:
                out.append(abs(p))
                break
            for j in range(t, n):
                a[t][j] += a[bad][j]
    return out


def factor(m: int) -> dict[int, int]:
    f, p = {}, 2
    while p * p <= m:
        while m % p == 0:
            f[p] = f.get(p, 0) + 1
            m //= p
        p += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def discriminant_invariants(gram):
    """|A_L|, exponent, p-ranks, root-norm count, log2 of an upper bound for |O(q_L)|.

    Root norms: possible_root_lengths = divisors d of 2*exp(A_L) (a root r has
    r^2 | 2 b(r, L), so r^2 | 2*exp); edgewalk.py:262, upstream edgewalk.h:55-67.
    |O(q_L)| <= prod_p |GL_{m_p}(Z/p^{e_p})|-type bound; log2 bound = sum_p m_p^2 log2(p^{e_p}).
    """
    inv = [d for d in smith_invariants(gram) if d != 1]
    order = math.prod(inv) if inv else 1
    exponent = max(inv) if inv else 1
    ranks: dict[int, int] = {}
    for d in inv:
        for p in factor(d):
            ranks[p] = ranks.get(p, 0) + 1
    norms = sum(1 for d in range(1, 2 * exponent + 1) if (2 * exponent) % d == 0)
    fe = factor(exponent)
    log_oq = sum(m * m * math.log2(p ** fe.get(p, 1)) for p, m in ranks.items())
    return order, exponent, ranks, norms, log_oq


ISOTROPIC = {json.dumps(c["gram"]): c["has_isotropic"] for c in _load("isotropic_cases_8821")}


def _det(m) -> Fraction:
    m = [[Fraction(x) for x in row] for row in m]
    n, d = len(m), Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if m[r][i]), None)
        if p is None:
            return Fraction(0)
        if p != i:
            m[i], m[p] = m[p], m[i]
            d = -d
        d *= m[i][i]
        for r in range(i + 1, n):
            f = m[r][i] / m[i][i]
            for c in range(i, n):
                m[r][c] -= f * m[i][c]
    return d


def _kernel(m) -> list[Fraction]:
    """One nonzero kernel vector of a singular square rational matrix (Gaussian elimination)."""
    n = len(m)
    a = [[Fraction(x) for x in row] for row in m]
    pivots, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, n) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        a[r] = [x / a[r][c] for x in a[r]]
        for i in range(n):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [x - f * y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
    free = next(c for c in range(n) if c not in pivots)
    v = [Fraction(0)] * n
    v[free] = Fraction(1)
    for row, c in enumerate(pivots):
        v[c] = -a[row][free]
    return v


def coxeter_vertices(gram, roots):
    """Vertices of the Coxeter polytope P bounded by the simple roots' walls.

    P is acute-angled, so (Vinberg, Hyperbolic reflection groups, Russian Math. Surveys 40
    (1985), Thm. 3.1 and Prop. 4.2):
      finite vertex  <-> elliptic subset of n-1 walls (Gram positive definite);
      ideal vertex   <-> parabolic subdiagram of rank n-2; every (n-1)-subset of it with
                         PSD rank-(n-2) Gram has its kernel combination on the same null ray.
    Returns a list of (roots at the vertex, invariant class key).
    """
    n = len(gram)
    form = [[sum(r[i] * gram[i][j] * s[j] for i in range(n) for j in range(n)) for s in roots] for r in roots]
    sign = 1 if form[0][0] > 0 else -1
    form = [[sign * x for x in row] for row in form]
    vertices: dict[tuple, set[int]] = {}
    for subset in itertools.combinations(range(len(roots)), n - 1):
        sub = [[form[i][j] for j in subset] for i in subset]
        if any(sub[i][j] * sub[i][j] > sub[i][i] * sub[j][j] for i in range(n - 1) for j in range(i)):
            continue
        minors = [_det([row[:k] for row in sub[:k]]) for k in range(1, n)]
        if all(x > 0 for x in minors):
            vertices[("f", subset)] = set(subset)
            continue
        if minors[-1] != 0:
            continue
        if any(_det([[sub[i][j] for j in idx] for i in idx]) < 0
               for k in range(1, n - 1) for idx in itertools.combinations(range(n - 1), k)):
            continue
        c = _kernel(sub)
        ray = [sum(c[t] * roots[i][x] for t, i in enumerate(subset)) for x in range(n)]
        g = Fraction(0)
        for x in ray:
            g = Fraction(math.gcd(g.numerator * x.denominator, x.numerator * g.denominator), g.denominator * x.denominator) if g else abs(x)
        key = tuple(x / g for x in ray)
        key = key if next(x for x in key if x) > 0 else tuple(-x for x in key)
        vertices.setdefault(("i", key), set()).update(subset)
    out = []
    for (kind, _), walls in vertices.items():
        idx = sorted(walls)
        sub = [[form[i][j] for j in idx] for i in idx]
        out.append((idx, (kind, len(idx), tuple(sorted(tuple(sorted(row)) for row in sub)))))
    return out


# Sym_L(P): diagram automorphisms of the recorded root system that are integral on L,
# computed once in Sage by A_symmetry_sage.py (an invariant of the case, not a timing).
with open(Path(__file__).with_name("A_symmetry.json"), encoding="utf-8") as _h:
    _SYM_BY_ID = json.load(_h)
SYMMETRY = {json.dumps(c["gram"]): _SYM_BY_ID[c["id"]] for c in _load("root_systems_103")}

ROOTS: dict[str, list] = {}
for _c in _load("allcock_rank3_reflective"):
    ROOTS.setdefault(json.dumps(_c["gram"]), _c["simple_roots"])
for _c in _load("root_systems_103"):
    ROOTS.setdefault(json.dumps(_c["gram"]), _c["roots"])


def edgewalk_cost(gram, walls: int, isotropic: bool, cusps: int | None = None, roots=None) -> tuple[float, str]:
    """Allcock edgewalk on a reflective Lorentzian lattice (edgewalk.py:1509-1618, edgewalk.h:1697-1792).

    Theory: references/complexity/vinberg-and-edgewalk.md, "Edgewalk branches" and "Main walk".
    Branch (edgewalk.py:710-716, edgewalk.h:1880-2055): cusp start if L is isotropic, cocompact
    start otherwise. Vertices: when the simple roots are recorded (Allcock table, root_systems_103,
    OSCAR, E10) the Coxeter polytope's vertices are enumerated exactly (coxeter_vertices); the
    walk visits at most that many (one per Sym(P)-orbit). Otherwise V <= R vertex orbits and
    two representatives per invariant class (stated assumption, see theory file).
    Transporter tests (edgewalk.py:1584-1588): each edge's neighbour is compared with the known
    configurations in order, stopping at the first match. A configuration with a different
    invariant (number of roots, root Gram) is rejected by its Gram (10 OP); a full test (STAB)
    happens only inside the neighbour's invariant class.
    """
    n = len(gram)
    order, exponent, ranks, norms, log_oq = discriminant_invariants(gram)
    roots = roots if roots is not None else ROOTS.get(json.dumps(gram))
    if roots is not None and walls:
        verts = coxeter_vertices(gram, roots)
        sym = SYMMETRY.get(json.dumps(gram), {"gens": [], "integral_order": 1})
        index = {tuple(w): i for i, (w, _) in enumerate(verts)}
        parent = list(range(len(verts)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for perm in sym["gens"]:
            for i, (w, _) in enumerate(verts):
                parent[find(i)] = find(index[tuple(sorted(perm[x] for x in w))])
        verts = [verts[i] for i in range(len(verts)) if find(i) == i]
        classes: dict[tuple, int] = {}
        for _, key in verts:
            classes[key[1:]] = classes.get(key[1:], 0) + 1
        v = len(verts)
        degree = [len(w) if key[0] == "f" else len(w) for w, key in verts]
        e = sum(degree)
        full_tests = sum(d * classes[key[1:]] for d, (_, key) in zip(degree, verts))
        ideal = sum(1 for _, key in verts if key[0] == "i")
        source = f"exact Coxeter vertex orbits V={v} under |Sym_L(P)|={sym['integral_order']} (ideal {ideal}) classes={len(classes)}"
    else:
        v = walls
        e = v * (n - 1)
        full_tests = e * 2
        ideal = cusps if cusps is not None else (1 if isotropic else 0)
        source = f"no recorded roots: V<=R={v}, 2 reps/class"
    search = math.log2(2 * exponent * max(order, 2))
    stab = STAB * max(1.0, n * n / 9) * (1 + log_oq / 64)
    t_vertices = v * stab + full_tests * stab + e * v * 10 * OP
    t_edges = e * norms * (30 * OP + search * 2 * OP * n)
    t_cert = v * (50 * OP + LP * walls * v)
    if isotropic:
        c = max(1, ideal)
        t_start = (1 + c) * norms * 40 * OP * (n - 2) ** 2
        branch = f"cusp(isotropic) cusps={c}"
    else:
        t_start = norms * (40 * OP * n * n + 1e-5 * math.sqrt(order))
        branch = "cocompact(anisotropic)"
    t = (BASE + t_start + t_vertices + t_edges + t_cert)
    basis = (f"[vinberg-and-edgewalk.md#edgewalk-branches] edgewalk {branch}: n={n} |A_L|={order} exp={exponent} "
             f"p-ranks={ranks} log2|O(q)|<={log_oq:.0f} norms={norms} walls={walls} {source} E={e} full-transporter-tests={full_tests}")
    return t, basis


def root_compare_cost(walls: int) -> float:
    """Graph.canonical_label on a walls-vertex labelled graph (test_reflective.sage:34-41)."""
    return 2 * (0.01 + 1e-4 * walls * walls)


REFLECTIVE_WALLS = {json.dumps(c["gram"]): c["num_simple_roots"] for c in _load("reflective_forms_8821")}


def isotropic_line_cost(case) -> tuple[float, str]:
    """isotropic_classes(L, 1) (consumer.py:223) -> isotropic_induction.py:170-200.

    Theory: references/complexity/isotropic-vector-search.md, sections "Branch selection"
    and "Lorentzian branch". Every fixture form is hyperbolic (positive index 1).
      anisotropic (rank 3/4, Hasse-Minkowski fails at some p | 2|A_L|; rank >= 5 never,
        Meyer): lattice.is_isotropic() only -> local Hilbert symbols.
      isotropic: indefinite reduction (CombinedAlgorithms.h:2483-2492) then orbits of null
        vertices of the Lorentzian perfect tessellation (lorentzian_perfect.h:1625-1644),
        whose cost is owned by perfect-domains-and-lorentzian-cells.md. Quotient cells C <=
        vertices of the Coxeter polyhedron (exact from recorded roots, else <= R); each cell
        has <= n adjacencies and is identified against <= 2 representatives of its class
        (lorentzian_cells.py:3514-3560).
    """
    gram = case["gram"]
    n = len(gram)
    order, exponent, ranks, norms, log_oq = discriminant_invariants(gram)
    primes = len(factor(2 * order))
    local = BASE + primes * 20 * OP
    if not case["has_isotropic"]:
        return local, f"[isotropic-vector-search.md#branch-selection] anisotropic (Hasse-Minkowski, n={n}): primes|2|A_L|={primes} |A_L|={order}"
    walls = REFLECTIVE_WALLS[json.dumps(gram)]
    roots = ROOTS.get(json.dumps(gram))
    if roots is not None:
        verts = coxeter_vertices(gram, roots)
        cusps = sum(1 for _, key in verts if key[0] == "i")
        cells = len(verts)
        source = f"exact Coxeter vertices={cells} cusps={cusps}"
    else:
        cells = walls
        cusps = max(1, walls // 4)
        source = f"no recorded roots: cells<=R={walls}"
    stab = STAB * max(1.0, n * n / 9)
    # lorentzian_cells.py:3514-3560: one local stabilizer per quotient cell, one transport of
    # marks per adjacency (n facets per cell), identification within the cell's invariant
    # class (<= 2 representatives) owned by perfect-domains-and-lorentzian-cells.md.
    t = (local + 10 * OP * n + cells * (stab + 40 * OP * n) + cells * n * (2 * stab + 10 * OP))
    return t, (f"[isotropic-vector-search.md#lorentzian-branch] isotropic Lorentzian perfect-cell orbits: n={n} |A_L|={order} "
               f"exp={exponent} reflective walls R={walls} {source}; quotient cells<=vertices, n adjacencies/cell, 2 reps/class")


def enriques_costs() -> dict[str, tuple[float, str]]:
    """N = U + U(2) + E8(-2): rank 12, A_N = (Z/2)^10, |O(q_N)| = 46998591897600 (fixture).

    Line orbits under O(N): 2 (primitive isotropic e with div 1 or 2); plane orbits: 2.
    Splitting O(N)-orbits along a subgroup (lattice_morphisms.py:4638 split_orbit) costs one
    coset/double-coset step per resulting orbit; for the stable group the index is
    |O(q_N)| and the split orbits number 528 / 24242 / 72199 (fixture counts), each needing
    a stabilizer image in O(q_N) acting on 2^10 points (GAP Schreier-Sims, ~5 ms each with
    wrapping; per-representative cost 0.042 s from B). O(q_N) order: Schreier-Sims on 1023 points, ~2 s.
    """
    q_order = 2.0
    image = 10.0  # generators of O(N), discriminant images, order of generated group
    # Constants from orthogonal-group-generators.md (owner B), sections "Isotropic
    # stabilizers" and "Lorentzian route": one line class = one Lorentzian traversal of
    # e^perp/e (rank 10, C <= 2) = 79.6 s; one plane class = 39.0 s; 
    line_class, plane_class = 79.6, 39.0
    # boundary-components.md (owner B), section "Splitting": one section lift per double
    # coset via GAP DoubleCosets, word length 300 at degree 527 = 0.042 s per output representative.
    split_rep = 0.042
    line = 2 * line_class    # 2 O(N)-classes of lines
    plane = 2 * plane_class  # 2 O(N)-classes of planes
    split_o_plus = 4 * split_rep  # O^+ split outputs 2 + 2 representatives
    out = {}
    out[ISO + "test_o_plus_n_has_two_line_and_two_plane_orbits"] = (
        (BASE + line + plane + split_o_plus),
        "[isotropic-vector-search.md#eichler-branch] Enriques N rank12 A=(Z/2)^10: line+plane O(N)-classes then O^+ index-2 split of 4 orbits")
    flags = 2 * plane_class  # one plane-class traversal per line class (flags through a line)
    stable = (528 + 24242 + 72199) * split_rep
    out[ISO + "test_stable_o_plus_n_orbits_match_dutour_sikiric_hulek"] = (
        (BASE + line + plane + flags + image + stable),
        "[isotropic-vector-search.md#eichler-branch] Enriques N: classes of lines/planes/flags, split along stable group of index |O(q)|~4.7e13 "
        "into 528+24242+72199 orbits x 0.042 s section lift (B: boundary-components.md#splitting)")
    out[ISO + "test_o_n_maps_onto_o_of_the_discriminant_form"] = (
        (BASE + q_order + image),
        "[isotropic-vector-search.md#eichler-branch] Enriques N: |O(q_N)| by Schreier-Sims on 1023 points; image of O(N) generators")
    for rank_id in ("rank0", "rank1"):
        classes = line if rank_id == "rank0" else plane
        out[ISO + f"test_published_representatives_lie_in_distinct_o_plus_orbits[{rank_id}]"] = (
            (BASE + classes + split_o_plus + 2 * 2 * 10.0),
            f"[isotropic-vector-search.md#eichler-branch] Enriques N {rank_id}: classes (B: 79.6/39.0 s) + O^+ split (0.042 s/rep) + 2x2 marked isometry tests (10 s each)")
    for rep, idx in (("I_1_1", 1), ("I_1_2", 527), ("I_2_1", 527), ("I_2_2", 23715)):
        t = BASE + image + (line_class if rep.startswith("I_1") else plane_class) + 2 + math.log2(idx + 1) * 0.5
        out[ISO + f"test_stabilizer_images_have_dutour_sikiric_hulek_indices[{rep}]"] = (
            t,
            f"[isotropic-vector-search.md#eichler-branch] Enriques N {rep}: setwise stabilizer in O^+(N) (parabolic, rank12) + image subgroup order; index {idx}")
    return out


def i_2_10_costs() -> dict[str, tuple[float, str]]:
    """I_{2,10} odd unimodular: A_L trivial, isotropic. Two vector orbits, two plane orbits.

    Unimodular branch: orbits classified by isotropic_reduction genus (parity) without
    discriminant candidates; each class: reduction (rank 10 / 8) + 1 isometry test.
    """
    return {
        ISO + "test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank0-primitive_isotropic_vector_orbits]": (
            (BASE + 2 * 3.0 + 1 * 5.0), "[isotropic-vector-search.md#eichler-branch] I_{2,10} unimodular isotropic: 2 line classes, reductions rank 10"),
        ISO + "test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank1-isotropic_plane_orbits]": (
            (BASE + 2 * 8.0 + 1 * 8.0), "[isotropic-vector-search.md#eichler-branch] I_{2,10} unimodular isotropic: 2 plane classes, reductions rank 8"),
    }


def estimates() -> dict[str, tuple[float, str]]:
    out: dict[str, tuple[float, str]] = {}
    for case in _load("reflective_forms_8821"):
        t, b = edgewalk_cost(case["gram"], case["num_simple_roots"], ISOTROPIC[json.dumps(case["gram"])])
        out[REFL + f"test_simple_root_counts_match_polyhedral_common[{case['id']}]"] = (t, b)
    for case in _load("root_systems_103"):
        t, b = edgewalk_cost(case["gram"], case["num_roots"], ISOTROPIC[json.dumps(case["gram"])])
        out[REFL + f"test_root_systems_match_polyhedral_common[{case['id']}]"] = (t + root_compare_cost(case["num_roots"]), b)
    for case in _load("allcock_rank3_reflective"):
        t, b = edgewalk_cost(case["gram"], case["num_simple_roots"], ISOTROPIC[json.dumps(case["gram"])])
        key = f"{case['id']}_w{case['weyl_group_id']}_l{case['lattice_id']}"
        out[REFL + f"test_root_systems_match_allcock[{key}]"] = (t + root_compare_cost(case["num_simple_roots"]), b)
    # OSCAR Vinberg cases: each Gram contains <1> + <-1> or U, hence isotropic (rank>=5 by Meyer).
    vinberg = {c["id"]: c for c in _load("oscar_lattice_oracles")["vinberg"]}
    for cid, case in vinberg.items():
        t, b = edgewalk_cost(case["gram"], case["num_simple_roots"], True, case.get("cusps"))
        out[REFL + f"test_simple_root_counts_match_oscar_vinberg[{cid}]"] = (t, b)
    t, b = edgewalk_cost(vinberg["vinberg_scharlau_U_A2"]["gram"], 4, True)
    out[REFL + "test_vinberg_coxeter_diagram_has_oscars_edges"] = (t, b)
    case = vinberg["vinberg_diag_1_-1_-3"]
    t, b = edgewalk_cost(case["gram"], case["num_simple_roots"], True, roots=case["simple_roots_of_test_chamber"])
    out[REFL + "test_vinberg_roots_match_oscars_chamber_up_to_isometry"] = (t + root_compare_cost(case["num_simple_roots"]), b)
    e10 = _load("e10_fundamental_domain")
    t, b = edgewalk_cost(e10["gram"], len(e10["simple_roots"]), True, 1, roots=[r["vector"] for r in e10["simple_roots"]])
    out[REFL + "test_e10_simple_roots_match_dutour_sikiric_hulek"] = (t + root_compare_cost(10), "E10=U+E8(-1) " + b)
    t, b = edgewalk_cost([[0, 1, 0], [1, 0, 0], [0, 0, 46]], 0, True, 1)
    # Nonreflective branch: three vertex orbits, then the finiteness test fails
    # (edgewalk.py:1538-1541): V=3, no certificate success.
    t46 = (BASE + 3 * STAB + 3 * 2 * 8 * (30 * OP + 2 * OP * 3 * 6))
    out[REFL + "test_u_plus_46_nonreflective_certificate_terminates_after_three_vertices"] = (
        t46, "[vinberg-and-edgewalk.md#worst-cases-and-nontermination] nonreflective U+<46>: |A_L|=46 norms=8, 3 vertex orbits then infinite-group stop")
    out[REFL + "test_promised_go_left_preserves_left_supplement"] = (
        (BASE + 50 * OP), "[vinberg-and-edgewalk.md#binary-branch] rank-2 promised_step: constant number of binary-form ops")
    out[REFL + "test_edgewalk_period_action_matrix_uses_coordinate_label_order"] = (
        (2 * BASE + 100 * OP), "[vinberg-and-edgewalk.md#binary-branch] two rank-2 lattices, two isometries, coordinate matrices")
    out[REFL + "test_binary_cycle_coordinate_order_and_period_certificate"] = (
        (2 * BASE + 2 * 200 * OP), "[vinberg-and-edgewalk.md#binary-branch] two rank-2 binary cycles (bound<=2): reduction cycle length <= 4")
    out[REFL + "test_binary_cycle_period_respects_permuted_coordinate_basis"] = (
        (2 * BASE + 2 * 200 * OP), "[vinberg-and-edgewalk.md#binary-branch] two rank-2 binary cycles bound 2")
    allcock = _load("allcock_rank3_reflective")
    sample = [edgewalk_cost(c["gram"], c["num_simple_roots"], ISOTROPIC[json.dumps(c["gram"])])[0]
              for c in allcock if c["gram"] in ([[0, -1, 0], [-1, 0, 0], [0, 0, 1]], [[-3, 0, 0], [0, 1, 0], [0, 0, 1]])]
    rs0 = _load("root_systems_103")[0]
    extra = [edgewalk_cost(rs0["gram"], rs0["num_roots"], True)[0],
             edgewalk_cost([[0, 1, 0], [1, 0, 0], [0, 0, 2]], 3, True)[0] * 2,
             t46,
             edgewalk_cost([[-1, 0, 0], [0, 3, 1], [0, 1, 7]], 6, False)[0] * 2,
             (BASE + 2 * STAB) * 2]
    out[REFL + "test_edgewalk_branch_representative_sample"] = (
        sum(sample) + sum(extra), "sum of branch sample: 2 Allcock + polyhedral root system + U+<2> (x2) + U+<46> + compact hexagon (x2) + 2 binary")
    out.update(enriques_costs())
    out.update(i_2_10_costs())
    for case in _load("isotropic_cases_8821"):
        out[ISO + f"test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[{case['id']}]"] = isotropic_line_cost(case)
    return out


if __name__ == "__main__":
    est = estimates()
    scope = [line.strip() for line in NODES.read_text().splitlines()
             if line.startswith((REFL, ISO))]
    assert set(scope) == set(est), (sorted(set(scope) - set(est))[:5], sorted(set(est) - set(scope))[:5])
    print("count", len(est))
    by = {}
    for key, (t, _) in est.items():
        by.setdefault(key.split("[")[0], []).append(t)
    for fn, ts in sorted(by.items()):
        print(f"{fn.split('::')[1]}: n={len(ts)} min={min(ts):.1f} median={statistics.median(ts):.1f} max={max(ts):.1f}")
