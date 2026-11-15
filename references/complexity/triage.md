# Suite run-time triage

The suite estimate (about 236 h, nearly all in one fire below) is the problem. A fix is worth doing only by the hours it removes from that estimate: a saving of seconds counts only when it compounds over enough cases to remove hours. The suite is not run whole while its estimate exceeds 6 hours. Work proceeds in this order:

1. Read the code of each family to find its time sinks (the theory files in this directory).
2. Test the estimates with small probes: a few cases per stratum, short timeouts. Reuse existing timings; never measure a node again at the same revision.
3. Total per family = measured mean per case × case count. Rank by total hours.
4. Price each sink by the method in README.md ('How to price a step'), as a function of its input's invariants. A negligible step that shows up in a profile is a defect. A hard step is priced by the analysis of its algorithm, never by a stated number of minutes.
5. Probe each sink directly: read the code to name the operation that costs the time, then time that one operation on its real input in a few lines that finish in seconds. The probe confirms the sink before a fix and its removal after; a whole test runs once, as the final check.
6. Fix the largest totals first: one large inefficiency, or a small one compounded over thousands of cases. After each fix, probe the same nodes again and update this file.
7. Run the suite only when the probe-backed estimate is under 6 hours.

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

## Correctness sinks fixed after 686a72a

| sage-extensions | Node | Probe: operation, input, before -> after | Node before -> after | Change |
| --- | --- | --- | --- | --- |
| bbb6608 | test_published_representatives_are_one_per_orbit[4] | finite G-set stabilizer generators, S_12 on its 220 three-subsets: 8.53 s -> 0.96 s (440 -> 4 generators); 1000 owned permutation point images: 0.100 s -> 0.002 s | over 300 s -> 131 s, passes | stabilizer as GAP PreImage(rho, Stabilizer(rho(G), x)); integer points keyed by engine integers. Remaining: pointed_isometry 61 s (rational-to-integral transporter 22 s, _lattice_action 13 s, definite isometry 9 s), cell construction 28 s |
| 9b29712 | test_equivalence_witness_carries_a_representative_to_its_reflection[r_0] | qfauto of the E10 complement E8 + <2d>: 26.5 s (d=7), 62.5 s (d=8), 148 s (d=9) -> 0.02 s at every d; definite isometry of two bases of E8 + <14>: over the 60 s cap -> 0.07 s | over 300 s -> 227 s, passes | definite Aut and isometry on the decomposition into indecomposables (Kneser's minimal-vector construction, unique by Eichler): S = 242 instead of S(d) = 1 + sum_{k<=d} 240 sigma_3(k). Remaining, about 25 s per representative of degree 4 and 6: the integral correction of the rational transport (rational_integral_stabilizers.py transporter: invariant_module 9.6 s through owned module preimages, finite residue orbit 14 s) |
| da0fa12 | test_equivalence_witness_carries_a_representative_to_its_reflection[r_0] | IntegralStructureAction.transporter on the r_0 pointed isometry: 11.4 s -> 3.3-5.1 s per call (invariant module 5.7 s -> 0; 90 integral preimages 7.3 s -> 1.3 s) | 227 s -> 209 s, passes | residues M/dM over the supplied stable bound M = N/d instead of the least invariant module (any G-stable M with dM <= L <= M carries the same orbit); integer preimages on one engine Smith form D = U A V per arrow. Remaining in the transporter: owned ring-element construction (63k raisings per call), orbit-word evaluation 2.1 s |
| 0a4460d | test_o_plus_n_has_two_line_and_two_plane_orbits | _FiniteIntegralRepresentation orbit sweep (right_cosets -> from_finite_orbits): 66.5 s -> 13.3 s (restriction 26.9 s -> under 1 s; isometry hash 22.9 s -> under 0.1 s) | over 180 s -> over 770 s cap, a later step holds the GIL (not yet located) | residue restriction as B^-1 A B on the engine; isometry equality and hash on the retained engine matrix; tensor.from_matrix in one engine read |
| 5d408ca | test_o_plus_n_has_two_line_and_two_plane_orbits | sampled stacks every 60 s: held in the kernel spanning family's eager name (printing a morphism computes endpoint signatures) for over 12 min -> line classes finish, plane classes (isotropic_induction.lift_invariant_frame) running at the 240 s cap | over 770 s -> over 240 s | indexed family name evaluated only on printing |

Edgewalk gates at bbb6608: w183 18.9 s, w217 17.0 s, ListReflect_2917 99.7 s, branch sample 28.6 s; all pass.

Edgewalk gates at 9b29712: w183 15.1 s, w217 14.9 s, ListReflect_2917 115.4 s, branch sample 36.0 s; all pass. published_representatives[4] at the same revision: 62.9 s.

Edgewalk gates at da0fa12: w183, w217_l26, ListReflect_2917 and branch sample pass; 191 s for the four in one session.

Sample nodes at 5d408ca (local container, wall time with Sage start-up, one node per session): test_e10_stabilizer_indices [2] 53 s, [4] 58 s, [8] 57 s; equivalence witnesses r_0 209 s (at da0fa12), r_1 124 s, r_7 135 s, r_8 107 s; all pass. Edgewalk gates at 5d408ca: the four pass, 194 s in one session.

## Sinks fixed after 5d408ca

| sage-extensions | Node | Probe: operation, input, before -> after | Node before -> after | Change |
| --- | --- | --- | --- | --- |
| 874f431 (3239c91, 6179cee); port 8f95b7d | test_o_plus_n_has_two_line_and_two_plane_orbits | lift_invariant_frame, plane classes: 1.7 s -> 0.8-1.1 s per call (rank-n rational linear algebra; price negligible). Residue stabilizer of the lattice in Sym(X): GAP PreImage(rho, Stabilizer(rho(Sym(X)), x)) on the full symmetric group held the GIL for over 10 min (the in-process cap never fired) -> 0.03 s at \|X\| = 2000 | over 19 min with the GIL held -> over the 950 s cap with the GIL released; the plane classes now reach representative descent and the second line class | Sym(X)_x named as Sym(X - x) by a transposition and a cycle; isometry composition at distinct endpoints by framed matrix product; F E F^-1 retained as the composite; sparse engine images skip coefficient re-conversion; denominator ideal from the framed matrix. Remaining in the node: `_extend_perpendicular_isometry` 2-5 s per call and `isotropic_stabilizer_action` 44 s per line class (the I_1_1 sinks), Lorentzian quotient cells for O(K) |
| a106d23 (9f3fa4e) | test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1] | `_extend_perpendicular_isometry` (rank-n rational linear algebra and one Smith form; price negligible), six calls under the profiler: 12.0 s -> 3.7-7.6 s (0.6-1.3 s per call); the node's first `isotropic_stabilizer_action` profiled at over the 150 s cap (107 s of it in the action: extensions 53.7 s, quotient O(K) generators 37.3 s) | over 650 s -> over 650 s | elements of QQ raised from the engine retain the engine value, and a lowered element memoizes it; exact tensor equality on engine components; zero and one retained per localization. Remaining: O(K) of the rank-n Lorentzian quotient through `quotient_cells` (31 s, `_facet_orbit_records` 17 s, facet registration 14.5 s), the orbit span, and per-entry lowering of matrices built from owned images |
| 5fabcdc (0c7f730); port 162fd5a | test_degree_two_k3_type_ii_components_and_root_types | stack samples every 5 s over 400 s (inclusive): O(K) of the Lorentzian quotient through `quotient_cells` 265 s, of which facet registration 225 s; `split_covering_vectors` 30 s; `_rational_witt_index` printed a kernel family's morphism (endpoint signatures) at construction. `_cell_facet_orbits_from_permutations` on the first quotient: 162 s -> 18 s | over 1200 s -> over 1300 s (cap; the four edgewalk gates pass in the same session) | the non-word facet walk records each decided incidence and lists every orbit of at most 20000 facets, so later members are found by lookup instead of GAP RepresentativeAction(OnSets) (SetFunctionality.g:174-201); every indexed-family name is a closure evaluated on printing. Remaining: the later stages of the node, not yet sampled after this fix |
| e4e901e (dc88f00) | test_simple_root_counts_match_polyhedral_common[ListReflect_2917] | `edgewalk_fundamental_domain` under the profiler: 161 s -> 22.8 s. The definite square shell of the vertex root enumeration (Fincke-Pohst; price: the shell size at the root norms) ran through the Julia bridge, 103 s with 7.5 s of bridge start-up | about 100 s -> passes; the four edgewalk gates take 84 s in one session (194 s at 5d408ca) | square shells by PARI qfminim with the exact-square filter; the Julia adapter and the coordinate-box branch deleted. Remaining in the profile: edge planes 8.3 s (`divisible_sublattice` 4.0 s), `Mor.Of` and object construction 3.4 s, `_new_sparse_free_module` 1.7 s |

w217 passes at 6179cee (one session with the O+(N) node) and at the I_1_1 revision. The four edgewalk gates pass at 0c7f730 with port 162fd5a (one session).

## Open sinks after e4e901e

Stack samples every 5 s (inclusive time), K3 type II degree two at sage-extensions e4e901e with port 11e0650: over 900 s, of which 780 s in O(K) of the Lorentzian quotient (`quotient_cells` traversal). The facet walk is now 100 s of it.

| Node | State at e4e901e | Sink | Class |
| --- | --- | --- | --- |
| test_degree_two_k3_type_ii_components_and_root_types, and the O(K) step of I_1_1 and O+(N) | over 900 s | P/lorentzian_cells.py `_traversal`: `cell_transporter` 380 s (SE lattice_morphisms.py `transporter` -> coloured_vector_families.py `_transporter`: owned morphism evaluation and owned ring equality per candidate), `flip_across` 230 s (`first_shell` 160 s, `_kernel_flipping` 145 s), `cell_stabilizer_permutations` 100 s. Price: one Lorentzian cell transporter is a coloured-vector-family isometry search on the cell's rays; the step counts are the source's | c: per-element owned construction in the transporter search; needs elements and morphisms that carry engine data (item-4 redesign) |
| test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1] | over 650 s at a106d23; not re-measured after the facet and square-shell fixes | the quotient O(K) above (37 s per action at a106d23), `_extend_perpendicular_isometry` 0.6-1.3 s per call, the orbit span | c |
| test_o_plus_n_has_two_line_and_two_plane_orbits | over 950 s at 6179cee with the GIL released; not re-measured after the facet and square-shell fixes | the quotient O(K) above per line class; `_integral_generator_set` in representative descent | c |
| test_equivalence_witness_carries_a_representative_to_its_reflection[r_0] | transporter 3.3-5.1 s per call (5d408ca) | owned ring-element raising and orbit-word evaluation | c |
| Item-4 redesign | not done | indexed families and `binary_candidates`; Mor_C(A, B) per endpoint pair and Category objects built in loops (`Mor.Of` 3.4 s in ListReflect_2917); engine-data constructors for `subobject_on`, `kernel`, `divisible_sublattice` (4.0 s in ListReflect_2917), `_new_sparse_free_module` (1.7 s), `_lattice`, `finite_ordered_set`, cone. Done in this pass: isometry composition at distinct endpoints is the framed matrix product; raised sparse elements and raised rationals carry their engine values; family names are lazy | blocker: a cross-cutting change of the element and morphism constructors in sage-extensions; each constructor needs its own probe and the edgewalk gates |
