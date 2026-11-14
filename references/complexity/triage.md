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
