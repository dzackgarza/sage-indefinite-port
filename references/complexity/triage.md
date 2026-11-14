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

At 24 s per case, the two edgewalk families still total about 2 x 8,600 x 24 s, roughly 115 h; the reference cost is under 1 s per lattice. The step counts match the source (one C_N per norm, one edge plane per edge and norm, one get_next per Dynkin extension, one ResolveLattEquation per parabolic extension); the remaining cost is the construction of owned objects, 2-100 ms each.

## Correctness subset (priority 1): sample nodes that failed in Acceptance 43505668003

Local container, no budget plugin, 300 s cutoff, sage-extensions changing during the run (14 November 2026). 46 of the 61 nodes reported; the job stopped before the other 15.

| Outcome | Nodes |
| --- | --- |
| Passed | 10 (8 s to 290 s; budgets 0.1 s to 8.9 s) |
| At the 300 s cutoff | 32 (budgets 0.07 s to 23 s; test_stable_o_plus_n budget 315 s) |
| biproduct_factor AttributeError (fixed at sage-extensions 3261c97, loaded before the fix) | 4 |

Families at the cutoff: vector orbits (published representatives, 11 nodes), equivalence witnesses (r_0, r_1, r_7, r_8), E10 stabilizer indices, d-tower components and pictures, Dawes buildings and incidence, isotropic stabilizer images, K3 type II components, O+(N) orbits, I_{2,10} plane classes. Every one exceeds its derived estimate by a factor of 100 or more, so each needs a diagnosis: the implementation does more work than the algorithm, or the model misses a term.
