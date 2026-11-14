# Suite run-time triage

The suite estimate (about 236 h, nearly all in one fire below) is the problem. A fix is worth doing only by the hours it removes from that estimate: a saving of seconds counts only when it compounds over enough cases to remove hours. The suite is not run whole while its estimate exceeds 6 hours. Work proceeds in this order:

1. Read the code of each family to find its time sinks (the theory files in this directory).
2. Test the estimates with small probes: a few cases per stratum, short timeouts. Reuse existing timings; never measure a node again at the same revision.
3. Total per family = measured mean per case × case count. Rank by total hours.
4. Fix the largest totals first: one large inefficiency, or a small one compounded over thousands of cases. After each fix, probe the same nodes again and update this file.
5. Run the suite only when the probe-backed estimate is under 6 hours.

## What establishes correctness

Correctness of a family is established by a branch-representative sample: unchanged cases that together exercise every branch of the pinned algorithm (AGENTS.md, #33 amendment, for A_REFLECTIVE, A_GROUPS and A_ISOTROPIC). Thousands of cases on one branch add no evidence that the sample lacks. The full corpus is nonblocking long-tail evidence. Priorities:

1. The minimal subset of tests that establishes the correctness of every capability (the branch-representative sample) passes, and runs fast.
2. The long tail is still an extreme problem (about 236 h). Its fires are fixed next, by hours removed, until its estimate is under 6 hours; only then does it run.

## Measured families (port 1cbb692 era, sage-extensions before 731fad7)

Source: `.tmp/probe/sample_report.jsonl` and `measure_report*.jsonl` (local container, 1200 s timeout). A stratum with timeouts is a lower bound.

| Test function | Rank | Cases | Measured | Mean s | At 1200 s timeout | Total h (lower bound if timeouts) |
|---|---|---|---|---|---|---|
| test_simple_root_counts_match_polyhedral_common | 3 | 8602 | 20 | 49.0 | 0 | 117.1 |
| test_root_systems_match_allcock | - | 8595 | 34 | 33.1 | 0 | 79.1 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 3 | 8602 | 8 | 16.6 | 0 | 39.7 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 6 | 84 | 4 | 61.5 | 0 | 1.4 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 9 | 4 | 4 | 777.1 | 1 | 0.9 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 15 | 2 | 2 | 1200.0 | 2 | 0.7 |
| test_simple_root_counts_match_polyhedral_common | 4 | 56 | 18 | 37.8 | 0 | 0.6 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 10 | 3 | 2 | 661.4 | 1 | 0.6 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 18 | 1 | 1 | 1200.0 | 1 | 0.3 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 16 | 1 | 1 | 1200.0 | 1 | 0.3 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 17 | 1 | 1 | 1200.0 | 1 | 0.3 |
| test_degree_two_k3_type_ii_components_and_root_types | - | 1 | 1 | 1200.0 | 1 | 0.3 |
| test_simple_root_counts_match_polyhedral_common | 6 | 84 | 19 | 13.2 | 0 | 0.3 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 5 | 47 | 4 | 22.4 | 0 | 0.3 |
| test_simple_root_counts_match_polyhedral_common | 5 | 47 | 15 | 6.9 | 0 | 0.1 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 4 | 56 | 5 | 5.4 | 0 | 0.1 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 8 | 4 | 4 | 42.3 | 0 | 0.0 |
| test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors | 7 | 4 | 4 | 26.3 | 0 | 0.0 |
| test_u_plus_46_reports_the_finiteness_reason | - | 1 | 1 | 11.2 | 0 | 0.0 |
| test_k3_norms_drop_the_norm_four_wall_of_u_plus_4 | - | 1 | 1 | 8.0 | 0 | 0.0 |

Sum of these rows: 242 h

## Fires, ranked

| Rank | Fire | Hours | Lever |
| --- | --- | --- | --- |
| 1 | ListReflect and Allcock rank-3 cases: `allcock_edgewalk` runs on each of 8,595 Gram matrices twice, once per test function (the Allcock Gram matrices are a subset of ListReflect's); IsotropicCases runs `isotropic_classes` on the same 8,818 matrices | ≈ 236 | Run each edgewalk once per session: the two functions construct separate lattice objects from one Gram matrix, so the per-object cache never hits. Then remove the per-call construction overhead in sage-extensions (sage-extensions#226) |
| 2 | IsotropicCases rank ≥ 9 | ≥ 3 | Isotropic search at high rank: compare with the pinned source's route |
| 3 | `test_degree_two_k3_type_ii_components_and_root_types` | ≥ 0.3 | Not yet read |

Unmeasured: the 145 Lorentzian isometry pairs, the 87 Enriques modular groups (every case errors before it computes), and the single-case tests.

## Fix log, fire 1 (case allcock_rk3_30_w183_l30, local container, unprofiled)

| sage-extensions | Per-case time | Change |
| --- | --- | --- |
| eac20b9 | 113 s (profiled) | baseline |
| 6d8ed19 | 50-57 s | presented-extension preimage solver; indexed-family ranking cached |
| c62175f | 45 s | Mor_C(A, B) arrow class built once per C, not per endpoint pair |
| 7f9b6cb | 31 s (27 s warm) | port: cusp roots by ResolveLattEquation on the cached rows of C_N, not an affine fibre per call |
| 7f9b6cb | 25 s | port: edge orientation once per edge; binary plane lattices shared by Gram |
| b77866a | 23.5 s (session with 4 earlier cases; 29.8 s before) | divisible sublattice of a ZZ lattice as one integral kernel |
| 05536b5 | 23.9 s (mean of six cases in one session: 27.0 s at b77866a) | lattice pairing as one engine product; port: split-plane null basis once per plane lattice, primitive plane vectors by integer gcd |
| a15ec92 | 24.1 s and 18.9 s (w183, w217; 23.9 s and 20.1 s at 05536b5 in the same container) | covectors and matrix-Mor sums, negatives, scalar multiples as engine matrices; to_covector as v^T G |
| 1481da9 | 19.7 s and 16.1 s (w183, w217) | category graph signature, depth and provider order key retained per category, not re-derived per class linearization |
| be530ea | 17.5 s and 17.0 s (w183, w217) | Witt extension raised from its engine matrix; the coisotropic solution admitted to its locus by construction; locus no longer a cached_method keyed by printing the inclusions |

At 24 s per case, the two edgewalk families still total about 2 x 8,600 x 24 s, roughly 115 h; the reference cost is under 1 s per lattice. The step counts match the source (one C_N per norm, one edge plane per edge and norm, one get_next per Dynkin extension, one ResolveLattEquation per parabolic extension); the remaining cost is the construction of owned objects, 2-100 ms each.

## Correctness subset (priority 1): sample nodes that failed in Acceptance 43505668003

Local container, no budget plugin, 300 s cutoff, sage-extensions changing during the run (14 November 2026). 46 of the 61 nodes reported; the job stopped before the other 15.

| Outcome | Nodes |
| --- | --- |
| Passed | 10 (8 s to 290 s; budgets 0.1 s to 8.9 s) |
| At the 300 s cutoff | 32 (budgets 0.07 s to 23 s; test_stable_o_plus_n budget 315 s) |
| biproduct_factor AttributeError (fixed at sage-extensions 3261c97, loaded before the fix) | 4 |

Families at the cutoff: vector orbits (published representatives, 11 nodes), equivalence witnesses (r_0, r_1, r_7, r_8), E10 stabilizer indices, d-tower components and pictures, Dawes buildings and incidence, isotropic stabilizer images, K3 type II components, O+(N) orbits, I_{2,10} plane classes. Every one exceeds its derived estimate by a factor of 100 or more, so each needs a diagnosis: the implementation does more work than the algorithm, or the model misses a term.

## Correctness subset diagnoses

One representative node per family, cProfile under a 180 s or 240 s cutoff, port at 1aa0b69, sage-extensions as loaded on 14 November 2026. Ranked by nodes at the cutoff times the factor over the estimate. SE is sage-extensions `src/sage_extensions/categories/`; P is `src/sage_indefinite_port/indefinite/`.

| Family | Nodes | Before | Class | Change or sink | After |
| --- | --- | --- | --- | --- | --- |
| test_published_representatives_are_one_per_orbit (and test_e10_stabilizer_indices, same route) | 11 + 3; profiled [2] and stabilizer [2] | > 240 s, budget 0.05 s | c | P/lorentzian_cells.py `_CellFaces.flag_simplices` (204 s of 240) orbits facet subsets through SE FiniteGSets: 62208 calls of SE/group/groups.py:7562 `_permutation_point_action`, 167 s, 2.7 ms per point image (owned FiniteOrderedSet ranking and element construction per point); SE/group/g_sets.py:1295 `_orbits` with 766 `in_same_orbit` decisions, 122 s. The pinned cell orbit work is a GAP permutation action at microseconds per point. | unchanged |
| test_equivalence_witness_carries_a_representative_to_its_reflection | 4 (r_0 profiled) | > 180 s, budget below 23 s | c | SE/lattices.py:6958 `discriminant_module`, 6 calls, 149 s, inside SE/lattice_morphisms.py:5155 `_isometry_decision`; 480 calls of SE/lattices.py:6877 `rational_pairing`, 87 s (0.18 s each), from discriminant_modules.py:951 `_descended_quadratic`. Recomputed per isometry decision of the same lattices. | unchanged |
| test_dtower_boundary_component_counts and pictures, test_dawes_buildings and incidence graph | 7 + 4 (N_3, dawes_2U_A2_stable profiled) | > 180 s | c | SE/lattice_morphisms.py:4651 `split_orbit` -> SE/isotropic_parabolics.py:3243 `quotient_reflection_generators` -> SE/lattices.py:9871 `shortest_vectors` -> SE/lattice_engines.py:1265 `definite_minimum`, 1 call, 91 s (N_3) and 82 s (Dawes); 81 s of it in lattice_engines.py:1257 `_bridge`. | unchanged |
| test_stabilizer_images_have_dutour_sikiric_hulek_indices | 2 (I_1_1 profiled) | > 180 s | c | SE/group/groups.py:5835 `group_generators` of a predicate subgroup: SE/group/predicate_subgroups.py:1119 `schreier_generating_morphism` -> `_preimage_coset_space` -> `_free_cover`, 110 s, 1 call. | unchanged |
| test_degree_two/degree_four K3 type II, test_o_plus_n two line and plane orbits, test_i_2_10 rank1, distinct O+ orbits rank0, stable O+(N) | 6 (degree_two, O+(N), I_2_10 rank1 profiled) | > 180 s | c, with one port term | SE/rational_integral_stabilizers.py:861 `right_cosets` -> `finite_representation` -> SE/group/g_sets.py:325 `from_finite_orbits`, 135 s (O+(N)), 92 s (I_2_10); SE `split_covering_vectors` 58 s and lattice_engines.py:645 `_rational_witt_index` 30 s (K3 degree two). Port term: P/eichler.py:211 `_source_hyperbolic_pair` (GetHyperbolicPlane, ApproximateModels.h:996-1180), 64 s: 175 `_source_find_isotropic` calls, 45 `_indefinite_lll` runs with 5464 Gram-Schmidt passes. Replacing the scalar Gram-Schmidt by Sage's compiled LDL^T factorization gave identical output and no measurable gain (107 s against 110 s on 300 random forms), so the cost is the LLL iteration count the source also performs, not a translation departure; not committed. The source-level count is not yet in the model of boundary-components.md. | unchanged |

## Correctness sinks fixed at sage-extensions 686a72a

| Node | Before | After | Change |
| --- | --- | --- | --- |
| d-tower N_3 | over 300 s | 77 s, passes | definite minimum by PARI qfminim, not the Julia bridge |
| Dawes dawes_2U_A2_stable | over 300 s | 144 s, passes | same |
| equivalence witness r_0 | over 300 s | over 300 s | discriminant module removed from the decision; the remaining sink is PARI qfisom on a rank-9 definite complement inside the pointed-isometry recursion |
| published representatives [4] | over 300 s | over 300 s | GAP orbits; still over the cutoff, not yet profiled again |

Edgewalk cases at this revision: w183 20.3 s, w217 15.9 s, ListReflect_2917 121.8 s, branch sample 28.0 s; all pass.
