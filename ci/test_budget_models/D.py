"""A-priori time budgets for scope D (facet orbits, Lorentzian cells, edgewalk, Eichler, Schreier).

Derived by inspection only: no test was run and no measured duration was read.
Each budget is: (branch taken) x (work fixed by that branch's mathematical invariants) x unit costs,
with no multiplicative factor. Paths cite src/sage_indefinite_port/indefinite/<file>:<line>.

Branch rules used throughout:
  * facet_orbits.py:461-471 (POLY_Heuristics StandardHeuristicSplitting): direct Normaliz dual description
    iff r - dim < 20, or (not order > 5000 and not (r >= 40 and order >= 100)); otherwise adjacency
    decomposition, recursing on each orbit representative's facet modulo its stabilizer (:512-548).
  * cell_cone (lorentzian_cells.py:1061): r < 10000 -> owned RationalPolyhedralCones.from_rays.
  * recursive.py:1134-1185: definite -> O().group_generators; indefinite -> reduction, then
    positive index 1 -> Lorentzian perfect complex; Witt index >= 2 -> Eichler _generation_step.
  * edgewalk.py:1530-1612: early_termination stops at the first infinite-order chamber isometry
    (finiteness test); reflective lattices stop at the root_chamber_certificate.
"""

from __future__ import annotations

import statistics


# Unit costs (runner: 4 vCPU GitHub ubuntu-latest, Sage 10.x).
C_OWN = 2e-3  # one owned-category lattice element op (Python dispatch in sage_extensions)
C_GAP = 2e-3  # one libGAP call on a permutation group of degree <= 64
C_RANK = 3e-4  # exact FLINT rank/kernel of a <= 64 x 8 integral matrix
C_NMZ = 0.05  # one Normaliz dual description, <= 64 rays, dim <= 8
C_PPL = 0.1  # owned from_rays cone construction (PPL) <= 30 rays, dim <= 5
C_PY = 1e-6  # one pure-Python set/tuple step
C_TRANS3 = 0.5  # one O(L).transporter of a vertex configuration, hyperbolic rank 3-4
SAGE_BOOT = 8.0  # one real fresh `python -c` child importing sage.all + sage_extensions + port (only test_generator_determinism spawns these)
C_LAT = 0.1  # one owned Lattices(ZZ)(gram) / HyperbolicLattices wrap with its first O() parent
C_ISO = 0.05  # one owned isometry construction from images, rank <= 7
C_VSTAB = 0.5  # one O(L).stabilizer(vertex configuration), hyperbolic rank 3-4
C_EDGE = 0.3  # one edge_step (wall-crossing enumeration) in rank 3-4
C_TINY = 0.01  # <= 10 exact arithmetic ops on <= 4x4 data
SPECIMEN = {
    # Exact invariants of the initial NONNEGATIVE_SQUARES cell, computed in the Sage container
    # (facet-orbits.md "Exact invariants of the test specimens"): rank n, rays r, facets F, facet
    # incidence size k, stabilizer order G, ray-action generators g, quotient cells c, facet orbits O.
    # A2: r-n=23 >= 20, |G|=48 <= 5000, r < 40 -> direct (facet_orbits.py:470).
    # A3: r=89 >= 40, |G|=384 >= 100, level 1 -> adjacency decomposition (:471); the single facet
    #     orbit's cone has 11 rays in dimension 4, r-d=7 < 20 -> direct at level 2.
    "A2_plus_1_plus_minus1": dict(n=4, r=27, F=6, k=9, G=48, g=4, c=1, O=1, split=False),
    "A3_plus_1_plus_minus1": dict(n=5, r=89, F=24, k=11, G=384, g=5, c=1, O=1, split=True),
}

R0 = dict(n=3, r=4, F=4, k=3, G=8, g=2, c=1, O=1, split=False)  # Result_Enumeration_0, exact (test docstring)


def _initial_cell(s) -> float:
    """lorentzian_cells.py:1243-1277: positive direction, initial vector search, then at most n-1
    _kernel_flipping steps (:6053), each a fixed-distance LDL enumeration over ~100*n candidates
    at ~n*C_OWN per candidate."""
    return 100 * s["n"] * C_OWN + (s["n"] - 1) * 100 * s["n"] * s["n"] * C_OWN * 0.1


def _cell_cone(s) -> float:
    """lorentzian_cells.py:1061-1064: r < 10000 branch, owned from_rays plus F owned covectors."""
    return C_PPL + s["F"] * s["n"] * C_OWN


def _stabilizer(s) -> float:
    """cell_transporter(all_images, generators_only) (:1956): basis-image backtrack over n levels with r
    candidates, pruned by the growing permutation subgroup; ~log2|G| generators, each certified by
    n*n owned pairings; plus r*r integral pairing invariants."""
    import math

    gens = s["g"]
    return s["n"] * s["r"] * s["r"] * C_PY * 50 + gens * (s["n"] ** 2 * C_OWN + s["r"] * s["n"] * C_OWN)


def _ray_perms(s) -> float:
    """facet_orbits.ray_permutations_from_stabilizer: each generator applied to r owned rays."""
    import math

    return s["g"] * s["r"] * s["n"] * C_OWN


def _direct_walk(s) -> float:
    """facet_orbits.py:472-492 direct branch: one Normaliz call, F orient/certify (rank of incidence
    rows), plus PermutationGroup order and orbit closure under G."""
    if s["split"]:
        # seed facet + one orbit rep: GAP facet stabilizer, recursive direct call on k rays in dim n-1,
        # one ridge rotation per ridge orbit, Balinski test.
        return 2 * C_RANK * s["r"] + 3 * C_GAP + C_NMZ + s["n"] * (2 * C_RANK + C_GAP) + s["F"] * s["g"] * s["k"] * C_PY
    return C_NMZ + s["F"] * 2 * C_RANK + 5 * C_GAP + s["F"] * s["g"] * s["k"] * C_PY


def _cone_facets_full(s) -> float:
    """facet_orbits.py:261-300 ridge closure: F facets x (n-1..F) ridges x pencil rank tests."""
    return s["F"] * s["F"] * 3 * C_RANK


def _flip(s) -> float:
    """flip_across (:1809): one pinned nullspace + one _kernel_flipping enumeration."""
    return 100 * s["n"] * s["n"] * C_OWN * 0.1


def _base(s) -> float:
    return C_LAT + _initial_cell(s)


THEORY = "references/complexity/facet-orbits.md"
SECTION = {
    "tests/test_facet_orbits.py": "#exact-invariants-of-the-test-specimens, #direct-branch, #adjacency-decomposition, #ridge-closure-cone_facets, #group-steps",
    "tests/test_facet_certificates.py": "#ridge-closure-cone_facets",
    "tests/test_facet_orbit_contracts.py": "#group-steps, #direct-branch",
    "tests/test_facet_traversal_seam.py": "#direct-branch, #group-steps; traversal: perfect-domains-and-lorentzian-cells.md (C)",
    "tests/test_cell_faces_facet_orbits.py": "#group-steps (flags of a cell)",
    "tests/test_lorentzian_flipping.py": "perfect-domains-and-lorentzian-cells.md (C)",
    "tests/test_edgewalk_options.sage": "vertex-orbit counts depend on vinberg-and-edgewalk.md (A), edgewalk branches",
    "tests/test_edgewalk_rank2.py": "vinberg-and-edgewalk.md (A)",
    "tests/test_generator_determinism.py": "covering-vector count depends on orthogonal-group-generators.md (B), Eichler generation step",
    "tests/test_integral_schreier_family.py": "#group-steps (Schreier generators)",
    "tests/test_dual_lll_solution_order.py": "#lll-ordered-solutions",
}
_CURRENT = [""]


def b(seconds: float, text: str) -> tuple[float, str]:
    return round(seconds, 3), f"{text}; theory {THEORY} {SECTION[_CURRENT[0]]}"


def facet_orbits(name: str, param: str) -> tuple[float, str]:
    """tests/test_facet_orbits.py."""
    if name == "test_square_edge_stabilizer_has_the_exact_induced_action":
        return b(15 * C_GAP, f"[{param}] D4 subgroup of order <=8 on 4 points; _stabilizer_on_facet GAP stabilizer; group construction")
    if name in ("test_same_ray_contract_is_proportionality_and_positive_orientation",
                "test_a_rational_coordinate_is_refused_rather_than_truncated"):
        return b(C_TINY, "pure integer predicates on <=4-vectors; no lattice")
    if name == "test_the_walk_refuses_a_cone_it_cannot_certify":
        return b(C_TINY + 2 * C_RANK, "3 rays dim 2: dimension==2 branch (:278), non-pointed pencil raises after one rank test")
    if name == "test_cube_cone_facets_and_orbits_match_the_closed_form":
        f = 2 * (6 + 8)  # n=4,5: 6 and 8 facets
        return b(f * 8 * 3 * C_RANK, "ridge closure on 3-cube and 4-cube cones (8,16 rays; 6,8 facets), certify each")
    if name == "test_six_cube_facet_records_retain_the_complete_signed_permutation_action":
        r, n = 64, 7
        perms = 3 * r * n * C_OWN * 2  # stabilizer image on rows + test's own ray_permutations
        walk = 2 * C_GAP + 12 * C_NMZ + 12 * 32 * C_PY * 50  # B6 order 46080>5000 -> split; facets are 5-cubes (32 rays, r-d=26, B5 3840, r<40 -> direct)
        return b(C_LAT + 3 * C_ISO + perms + walk, "64 rays dim 7, B6 |G|=46080 -> adjacency-decomposition split (:466); one facet orbit of 12; 3 owned generators on 64 rays")
    if name == "test_null_quadrant_facets_retain_the_swap_transport":
        return b(C_LAT + C_ISO + C_NMZ, "U, 2 rays dim 2: dimension==2 branch; one swap generator; one orbit of 2")
    if name == "test_signed_cube_facet_representatives_have_the_full_orbit_size":
        sign, d = param.split("-")
        d = int(d)
        r = 2 ** d
        if d < 6:
            return b(10 * C_GAP + C_NMZ + 2 * d * 2 * C_RANK + 10 * C_GAP,
                     f"[{param}] {r} rays dim {d+1}, r-d<20 -> direct Normaliz; {2*d} facets; retain_words=False RepresentativeAction per pair")
        order = 46080 if sign == "True" else 720
        rule = ">5000 split" if sign == "True" else "r>=40 & order>=100 split"
        orbits = 1 if sign == "True" else 2
        sub = orbits * (C_NMZ + 10 * C_GAP)
        return b(10 * C_GAP + sub + 2 * d * 4 * C_GAP + 2 * d * 2 * C_RANK,
                 f"[{param}] 64 rays dim 7, |G|={order} ({rule}); {orbits} orbit(s) recurse to 5-cube facet (32 rays, direct); GAP RepresentativeAction/Stabilizer per candidate")
    if name == "test_second_hypersimplex_has_the_complete_coordinate_facet_orbits":
        moving, words = param.split("-")
        order = 5040 if moving == "7" else 40320
        orbits = 4 if moving == "7" else 2
        sub = orbits * (C_NMZ + 30 * C_RANK)
        gap = 0 if words == "True" else 16 * orbits * C_GAP
        return b(10 * C_GAP + sub + gap + 16 * 28 * C_PY,
                 f"[{param}] Delta(8,2): 28 rays dim 8, r-d=20, |S_{moving}|={order}>5000 -> split; {orbits} orbits; facets Delta(7,2) (21 rays) and simplex (7 rays) -> direct")
    s = SPECIMEN[param]
    common = _base(s)
    if name == "test_walk_reproduces_the_engine_facet_family_exactly":
        t = common + _cone_facets_full(s) + _cell_cone(s) + s["F"] * s["n"] * 4 * C_OWN
        return b(t, f"[{param}] initial cell (n={s['n']}, r~{s['r']}); full ridge closure F~{s['F']}; owned cell_cone; covector coords")
    if name == "test_orbits_partition_the_facets_and_words_transport":
        t = common + _stabilizer(s) + _ray_perms(s) + _cell_cone(s) + _direct_walk(s)
        return b(t, f"[{param}] stabilizer |G|~{s['G']} transporter; walk branch split={s['split']}; words over F~{s['F']}")
    if name == "test_pruned_walk_finds_every_orbit_exactly_once":
        t = common + _stabilizer(s) + _ray_perms(s) + _cone_facets_full(s) + _cell_cone(s) + _direct_walk(s)
        return b(t, f"[{param}] stabilizer + full closure + engine cone + direct walk")
    if name == "test_orbit_stabilizer_and_transport_witnesses_are_exact":
        t = common + _stabilizer(s) + 2 * _ray_perms(s) + _direct_walk(s) + _cone_facets_full(s) + s["G"] * s["F"] * s["r"] * C_PY * 5
        return b(t, f"[{param}] stabilizer (cached second call); Python closure of |G|~{s['G']}; per-orbit stabilizer count")
    if name == "test_walk_permutations_match_the_port_action_convention":
        t = common + _stabilizer(s) + 2 * _ray_perms(s) + _cell_cone(s)
        return b(t, f"[{param}] stabilizer + _facet_mask_data(use_owned_cone) masks")
    raise KeyError(name)


def facet_certificates(name: str, param: str) -> tuple[float, str]:
    if param:
        s = SPECIMEN[param]
        common = _base(s)
        if name == "test_every_engine_facet_normal_is_recovered_from_its_incidence":
            return b(common + _cell_cone(s) + s["F"] * (3 * C_RANK + s["n"] * 4 * C_OWN), f"[{param}] initial cell, owned cone F~{s['F']}, cofactor normal per facet")
        if name == "test_adjacency_certificate_matches_the_ridge_rank_criterion":
            pairs = s["F"] * (s["F"] - 1) // 2
            return b(common + _cell_cone(s) + pairs * 2 * C_RANK, f"[{param}] initial cell, F~{s['F']}, {pairs} facet pairs x 2 ranks")
        if name == "test_pointedness_is_certified_by_the_constructor_hypothesis":
            return b(common + 2 * C_RANK, f"[{param}] initial cell + one augmented kernel")
    if name == "test_cube_cone_adjacency_matches_the_closed_form":
        return b((15 + 45) * 2 * C_RANK, "closed-form cube facets n=4,6; 15+45 pairs, one ridge rank each")
    return b(C_TINY, "a handful of exact ranks/determinants on <=3x3 matrices; no lattice")


def facet_orbit_contracts(name: str, param: str) -> tuple[float, str]:
    s = SPECIMEN[param]
    common = _base(s)
    if name.startswith("test_facet_covectors"):
        return b(common + _cell_cone(s) + s["F"] * (s["n"] * 4 * C_OWN + C_RANK), f"[{param}] owned cone F~{s['F']}, rank per facet")
    if name.startswith("test_stabilizer_permutes"):
        return b(common + _stabilizer(s) + _ray_perms(s) * 2 + _cell_cone(s), f"[{param}] stabilizer |G|~{s['G']}, owned images of r~{s['r']} rays per generator")
    return b(common + _cell_cone(s) + 100 * s["n"] * C_OWN, f"[{param}] owned cone + _positive_direction short-vector search")


def facet_traversal_seam(name: str, param: str) -> tuple[float, str]:
    s = SPECIMEN[param]
    common = _base(s) + _stabilizer(s) + _ray_perms(s)
    if name == "test_walk_orbits_match_the_engine_route_at_the_seam":
        return b(common + _cell_cone(s) + 2 * _direct_walk(s), f"[{param}] engine (_facet_mask_data, Normaliz-free owned route) and walk orbits")
    if name == "test_walk_words_transport_and_neighbours_agree":
        return b(common + _cell_cone(s) + 2 * _direct_walk(s) + 2 * s["O"] * _flip(s), f"[{param}] both routes + flip_across on O~{s['O']} reps per route")
    # complete traversal twice (cache cleared): per route, c quotient cells; each: stabilizer, walk,
    # O flips, and a cell_transporter of each neighbour against up to c known cells; then
    # Poincare generator assembly.
    per_route = s["c"] * (_stabilizer(s) + _cell_cone(s) + _direct_walk(s) + s["O"] * (_flip(s) + s["c"] * _stabilizer(s)))
    return b(_base(s) + 2 * per_route, f"[{param}] two full perfect_complex traversals, c~{s['c']} cells, O~{s['O']} facet orbits/cell, c*O transporters")


def cell_faces(name: str, param: str) -> tuple[float, str]:
    # Exact (test docstring, Result_Enumeration_0): 4 rays, 4 facets, |G|=8, 1 orbit, 1 flag.
    # Module fixture (Result_Enumeration_0, rank 3, 4 rays, 4 facets, |G|=8) charged to the first test.
    fixture = _base(R0) + _stabilizer(R0) + _cell_cone(R0)
    own = {
        "test_facet_set_and_incidences": fixture,
        "test_cell_group_order": C_GAP,
        "test_facet_orbits_under_the_cell_group": 4 * C_GAP,
        "test_flag_count_under_the_cell_group": 3 * 2 * C_GAP + 3 * C_RANK,
        "test_facet_stabilizer_is_a_subgroup_of_the_cell_group": 3 * C_GAP,
    }[name]
    return b(own, "rank-3 cell, 4 rays/4 facets, |G|=8, one facet orbit, one flag; fixture cost on first test")


def lorentzian_flipping(name, param):
    return b(C_LAT + 2 * 100 * 4 * C_OWN * 0.1, "<1>+<-1>, two _kernel_flipping calls on a 1-ray face; slice x=1 has 3 vectors")


def edgewalk_options(name, param):
    hyp = C_LAT
    walk_reflective = 2 * (C_VSTAB + 3 * C_EDGE) + C_NMZ  # <=2 vertices, 3 walls, one chamber certificate
    if name == "test_full_walk_on_u_plus_46_closes_every_edge":
        V, E = 4, 3  # nonreflective; early_termination=False walks until vertex quotient closes (>=3 vertices)
        walk = V * (C_VSTAB + E * (C_EDGE + V * C_TRANS3))
        check = V * E * (C_EDGE + V * C_TRANS3)
        return b(hyp + walk + check, f"U+<46> not reflective, full walk (early_termination False, edgewalk.py:1542): V~{V} vertex orbits, E~{E} edges/vertex, V*E*V transporters, recheck in test")
    if name == "test_k3_norms_drop_the_norm_four_wall_of_u_plus_4":
        return b(hyp + 2 * walk_reflective, "U+<4> reflective (all and K3 norms): certificate after <=2 vertices each, two walks")
    if name == "test_prescribed_initial_vertex_on_u_plus_2":
        V = 2
        return b(hyp + (1 + 2 * V) * walk_reflective, f"U+<2> reflective, chamber ~3 walls, V~{V} vertices; 1 + 2V reflective walks stopping at certificate")
    if name == "test_u_plus_46_reports_the_finiteness_reason":
        return b(hyp + 2 * C_VSTAB + 2 * C_TRANS3, "U+<46> early termination: stops at first infinite-order isometry (edgewalk.py:1584/1608), <=2 vertices")
    if name == "test_edgewalk_isometry_of_u_plus_2_under_a_basis_change":
        return b(2 * hyp + 2 * walk_reflective + 2 * C_TRANS3 * 4, "edgewalk_isometry: two reflective walks of U+<2> (early_termination False, :119) + vertex transporters")
    return b(hyp + C_VSTAB + C_NMZ + 4 * 2 * C_EDGE, "U+<2>+<2> K3: isotropic cusp vertex with 4 roots (2 affine A1); 4 facets x 2 omitted walls edge_step")


def edgewalk_rank2(name, param):
    if name == "test_binary_cycle_form_uses_the_coordinate_label_order":
        return b(C_LAT + 20 * 10 * C_OWN, "anisotropic binary form disc 2^2+1... cycle of reduced forms, period of length O(log); <=10 owned ops/step")
    return b(C_LAT + 10 * C_OWN, "rank-2 lattice construction and coordinate frame only")


def generator_determinism(name, param):
    # U(2)+U+<2>: rank 5, signature (3,2), positive index 3, Witt index 2 -> Eichler _generation_step
    # (recursive.py:1180). A_L = (Z/2)^3, |A_L|=8; covering vectors ~ |O(q_L)| orbits (<=6), one
    # pointed isometry + envelope stabilizer each, recursion into rank-4 U(2)+<..> complement.
    compute = 8 * (C_TRANS3 + 2 * C_ISO) + 2.0  # 8 covering vectors x pointed isometry; reduction + envelope
    return b(2 * (SAGE_BOOT + compute), "two fresh Sage processes; U(2)+U+<2> Eichler generation step, |A_L|=8, ~8 covering vectors")


def integral_schreier(name, param):
    if name == "test_eichler_envelope_retains_the_supplied_lattice_action":
        return b(C_LAT + 2.0 + 6 * 4 * 2 * 4 * C_OWN, "U(2)+U envelope, |A_L|=4: stabilizer of the lattice in finite rep, Schreier gens <=6, 2*4 owned images each")
    if param == "two_u":
        return b(C_LAT + 7 * C_ISO + 7 * 4 * 2 * 4 * C_OWN, "U+U unimodular, bound = inclusion: index 1, one coset, 7 Schreier generators")
    if param == "scaled_hyperbolic_line":
        return b(C_LAT + 3 * C_ISO + 4 * 3 * 3 * 2 * 3 * C_OWN, "U(4)+<-8>, bound index 8 (1/4 e+1/2 w), <=4 cosets x 3 generators")
    return b(C_LAT + 3 * C_ISO + 18 * 3 * 4 * C_OWN, "S3 x sign on <1>^3 with index-6 sublattice: 6 cosets, 18 Schreier generators")


def dual_lll(name, param):
    return b(C_TINY, "5 solution tuples sorted by a 3x3 rational transform; pure ordering")


FILES = {
    "tests/test_facet_orbits.py": facet_orbits,
    "tests/test_facet_certificates.py": facet_certificates,
    "tests/test_facet_orbit_contracts.py": facet_orbit_contracts,
    "tests/test_facet_traversal_seam.py": facet_traversal_seam,
    "tests/test_cell_faces_facet_orbits.py": cell_faces,
    "tests/test_lorentzian_flipping.py": lorentzian_flipping,
    "tests/test_edgewalk_options.sage": edgewalk_options,
    "tests/test_edgewalk_rank2.py": edgewalk_rank2,
    "tests/test_generator_determinism.py": generator_determinism,
    "tests/test_integral_schreier_family.py": integral_schreier,
    "tests/test_dual_lll_solution_order.py": dual_lll,
}


def _scope_nodes() -> list[str]:
    from pathlib import Path

    nodes = Path(__file__).resolve().parent / "nodes_initial.txt"
    return [line.strip() for line in nodes.read_text().splitlines() if line.split("::")[0] in FILES]


def estimates() -> dict[str, tuple[float, str]]:
    result = {}
    for node in _scope_nodes():
        path, rest = node.split("::")
        name, _, param = rest.partition("[")
        _CURRENT[0] = path
        result[node] = FILES[path](name, param.rstrip("]"))
    return result


if __name__ == "__main__":
    table = estimates()
    assert set(table) == set(_scope_nodes()) and len(table) == 75, len(table)
    print("count", len(table))
    groups: dict[str, list[float]] = {}
    for node, (seconds, _basis) in table.items():
        groups.setdefault(node.split("[")[0], []).append(seconds)
    for function, values in sorted(groups.items()):
        print(f"{min(values):8.2f} {statistics.median(values):8.2f} {max(values):8.2f}  {function}")
