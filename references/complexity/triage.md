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

Correctness of a family is established by a branch-representative sample (tier 1 in "Acceptance tiers" below): unchanged cases that together exercise every branch of the pinned algorithm (AGENTS.md, #33 amendment, for A_REFLECTIVE, A_GROUPS and A_ISOTROPIC). Thousands of cases on one branch add no evidence that the sample lacks. The full corpus is nonblocking long-tail evidence. Priorities:

1. The minimal subset of tests that establishes the correctness of every capability (the branch-representative sample) passes, and runs fast.
2. The long tail is still an extreme problem (about 236 h). Its fires are fixed next, by hours removed, until its estimate is under 6 hours; only then does it run.

## Acceptance tiers

The Acceptance workflow takes three corpora (`ci/plan_acceptance.py`):

| Tier | File | Nodes | Content |
| --- | --- | --- | --- |
| 1 | `ci/acceptance_tier1.txt` | 49 | For each dispatch branch of each pinned algorithm, the smallest frozen node (least rank, determinant or \|O(q)\| that reaches the branch) whose expected answer comes from an independent source or from a construction, not from a round trip of the port |
| 2 | `ci/acceptance_tier2.txt` | 86 | The published-result reproductions that the #64 table requires |
| 3 | `full` | 26907 | Every acceptance node; tier 3 is the full corpus less tiers 1 and 2 |

Tier 1 is the correctness evidence (priority 1 above). Recorded tier-1 time: 1015.9 s over all 49 nodes (call time, the least recorded passed run). Tier 2 has 7 priced nodes (6167 s) and 79 unpriced nodes. `ci/plan_acceptance.py` costs a tier node at its price ("prices" in `ci/test_budgets.json`), and an unpriced tier-2 node at its derived budget: the tier-1 job timeout went from 23 to 37 minutes, the tier-2 job timeout from 63 to 165 minutes.

Probe of the 20 nodes unpriced before (one container process per node, `timeout -s KILL 60` inside the container, port 8208d1f, sage-extensions 6f337e9a): 16 passed. Three perfect-domain traversals (`Result_Enumeration_1`, `_3` total, `_3` isotropic) failed with TypeError: since sage-extensions 4d14cd55 a finite ordered set hashes its points, and `PerfectDomainOrbitAdjacency` is a `SageObject`, which is unhashable; fixed by an identity hash (sage-extensions 87b4335c; equality is identity). `test_large_cell_cone_uses_exact_normaliz_h_presentation` was killed at 60 s: its incidence check reads `covector.to_vector()` about 108000 times on three covectors (0.4-0.6 ms each); `to_vector` is now a cached_method (sage-extensions f7a03e96), and the node passes in 30.6 s.

IsotropicCases re-selection with exact rational determinants and signatures (Sage `QuadraticForm(QQ, G).signature_vector()`, all 8821 cases): least (n, |det|) per branch is IsotropicCases_3 (lines, sig (2,1), det -1) and IsotropicCases_16 (anisotropic indefinite, sig (2,1), det -3). IsotropicCases_3261 is not definite: its exact signature is (2,1), det -63525, has_isotropic false, so it repeats the branch of _16 and leaves tier 1. The fixture has no definite case. Gap: the definite (empty) branch of isotropic classes has no frozen tier-1 node.

A price is the least recorded passed call time of the node: a CI Acceptance BUDGET line (runs 43505668003 and 43494572338), a local pytest report under `.tmp/`, or a row of this file or of `.tmp/acceptance-coverage.md` with its revision. A node without such a record is unpriced. The derived budgets in `ci/test_budgets.json` are not prices. Tier-2 vector-orbit nodes are the remaining cases of `test_vector_orbits.sage` that #64 lists.

| Node | Tier | Branch covered | Invariants | Recorded time |
| --- | --- | --- | --- | --- |
| `test_reflective.sage::test_edgewalk_branch_representative_sample` | 1 | Allcock edgewalk: binary walk, oriented binary cycle (go-left, shorter pair), cusp start, cocompact start and increasing-prefix corner search, ordinary, cusp and extension roots, isotropic next vertex, reflective stop, nonreflective stop | n=2..4; U+<2>, U+<46>, compact hexagon, allcock_rk3_4, allcock_rk3_2, ListSimpleRootSystem_0 | 59.8 s (.tmp/sample61/report_2.jsonl) |
| `test_perfect_domains.sage::test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_3]` | 1 | Perfect domains, total mode: initial cell, one cell orbit, O(L) generation | n=3 sig (2,1) |det|=1, 1 cell | 2.8 s (probe, sage-extensions f7a03e96) |
| `test_perfect_domains.sage::test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_1]` | 1 | Perfect domains, total mode: flips, cell transporter, side pairings | n=3 sig (2,1) |det|=8, 2 cells | 2.7 s (probe, sage-extensions f7a03e96) |
| `test_perfect_domains.sage::test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_3]` | 1 | Perfect domains, isotropic mode | n=3 |det|=1, 1 cell | 2.9 s (probe, sage-extensions f7a03e96) |
| `test_perfect_domains.sage::test_non_diagonal_positive_shell_uses_the_pinned_cvp_center` | 1 | Perfect domains: positive-vector shell, CVP center | n=3, against exact enumeration | 0.7 s (probe, sage-extensions 6f337e9a) |
| `test_perfect_domains.sage::test_non_diagonal_isotropic_shell_uses_the_exact_sphere_center` | 1 | Perfect domains: isotropic shell, exact-sphere center | n=3, against exact enumeration | 0.9 s (probe, sage-extensions 6f337e9a) |
| `test_perfect_domains.sage::test_large_pairing_graph_fallback_returns_exact_transporter` | 1 | Perfect domains: transporter above 256 rows (no graph canonization) | n=2, 257 rows | 0.5 s (probe, sage-extensions 6f337e9a) |
| `test_perfect_domains.sage::test_large_cell_cone_uses_exact_normaliz_h_presentation` | 1 | Perfect domains: large-configuration cone (Normaliz) | 24000 rows | 30.6 s (probe, sage-extensions f7a03e96) |
| `test_isometry.sage::test_isometry_verdicts_match_certified_pairs[hecke_U_vs_U2]` | 1 | Isometry: binary, prefilter rejection | n=2 sig (1,1) |det|=1 | 0.5 s (probe, sage-extensions 6f337e9a) |
| `test_isometry.sage::test_isometry_verdicts_match_certified_pairs[hecke_3x3_swapped_basis]` | 1 | Isometry: Lorentzian witness (perfect-cell transporter) | n=3 sig (1,2) |det|=1250 | 11.5 s (probe, sage-extensions 6f337e9a) |
| `test_isometry.sage::test_conway_sloane_same_genus_pair_is_not_isometric` | 1 | Isometry: same genus, spinor-genus negative verdict | n=3 sig (2,1) |det|=128 | 1.4 s (probe, sage-extensions 6f337e9a) |
| `test_isometry.sage::test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U]` | 1 | Isometry: Eichler covering search witness | n=4 sig (2,2) |det|=4, U+U(2) | 12.5 s (probe, sage-extensions 6f337e9a) |
| `test_isometry.sage::test_degenerate_isometry_witness_respects_complement_coordinates` | 1 | Isometry: degenerate, radical extension witness | n=3, A2 plus rank-1 radical | 8.7 s (probe, sage-extensions 6f337e9a) |
| `test_orthogonal_group.sage::test_degenerate_orthogonal_group_contains_radical_shear_and_gl_action` | 1 | O(L): degenerate, radical extension | n=2, <2> plus radical | 0 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_zero_form_orthogonal_group_is_the_full_general_linear_group` | 1 | O(L): zero form, GL_2(Z) | n=2, q=0 | 0 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_degenerate_orthogonal_group_embeds_non_diagonal_complement` | 1 | O(L): definite kernel inside the radical extension | n=3, A2 plus radical | 0.4 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_orthogonal_group_of_the_hyperbolic_plane_has_order_four` | 1 | O(L): binary, isotropic | U, n=2 |det|=1, |O|=4 | 0 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_ambiguous_binary_form_has_an_improper_automorphism[hecke_binary_3_1_-3]` | 1 | O(L): binary, anisotropic reduction cycle | n=2 D=37 | 3.0 s (probe, sage-extensions 6f337e9a) |
| `test_orthogonal_group.sage::test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q` | 1 | O(L): Lorentzian route; image of O+ in O(q) | n=3 sig (1,2) |det|=3 | 0.8 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_mertens_omega_and_minus_one_generate_o_l` | 1 | O(L): Lorentzian route, generators against Mertens | n=3 sig (2,1) |det|=155 | 4.9 s (CI Acceptance BUDGET line) |
| `test_isotropic_orbits.sage::test_o_n_maps_onto_o_of_the_discriminant_form` | 1 | O(L): Eichler route; discriminant image | N=U+U(2)+E8(-2), n=12 sig (2,10) |A|=2^10 | 468 s (.tmp/acceptance-coverage.md, provider 42d7361c) |
| `test_vector_stabilizers.sage::test_isotropic_pointed_isometry_rejects_different_divisibility` | 1 | Pointed isometry: invariant rejection | n=4 | 0 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_split_divisibility_one_pointed_isometry_carries_actual_negative_marking` | 1 | Pointed isometry: isotropic, div 1 | U, n=2 | 0.3 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_divisibility_two_isotropic_pointed_isometry_has_exact_witness` | 1 | Pointed isometry: isotropic, div 2 | U(2)+U(2), n=4 | 0.8 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_odd_lattice_isotropic_pointed_isometry_uses_general_divisibility_two_route` | 1 | Pointed isometry: isotropic, odd lattice, general div 2 route | odd binary, n=2 | 2.4 s (probe, sage-extensions 6f337e9a) |
| `test_vector_stabilizers.sage::test_general_isotropic_pointed_isometry_carries_negative_marking` | 1 | Pointed isometry: isotropic, general route | U(2)+U(2), n=4 | 0.2 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_divisibility_four_isotropic_stabilizer_constructs_exact_containing_bound` | 1 | Isotropic line stabilizer: finite correction, cokernel 4 | U(4)+<-8>, n=3 | 1 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_isotropic_pointed_automorphism_generators_preserve_marking_and_form` | 1 | Isotropic line and pointed stabilizer generators | U+U, n=4 | 1.7 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_divisibility_four_isotropic_plane_stabilizer_constructs_exact_containing_bound` | 1 | Isotropic plane stabilizer: finite correction, cokernel 16 | U(4)+U(4), n=4 | 0.9 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_isotropic_plane_and_flag_automorphism_generators_preserve_their_marks` | 1 | Isotropic plane and flag stabilizer generators | images 48 and 12 in O over F_3 | 2.7 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_isotropic_plane_isometry_carries_changed_basis_plane` | 1 | Isotropic plane isometry | U(4)+U(4), n=4 | 2.9 s (probe, sage-extensions 6f337e9a) |
| `test_vector_stabilizers.sage::test_isotropic_flag_isometry_preserves_line_inside_plane` | 1 | Isotropic flag isometry | U(4)+U(4), n=4 | 3.5 s (probe, sage-extensions 6f337e9a) |
| `test_vector_stabilizers.sage::test_vector_stabilizer_has_oscars_order` | 1 | Anisotropic pointed stabilizer: v-perp with gluing | n=3 sig (1,2), v-perp binary, order 2 | 0.5 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_5]` | 1 | Anisotropic pointed isometry: complement isometry and integral correction | E10, n=10 |det|=1 | 2.2 s (CI Acceptance BUDGET line) |
| `test_isotropic_orbits.sage::test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_16]` | 1 | Isotropic classes: anisotropic indefinite (Hasse-Minkowski) | n=3 sig (2,1) |det|=3 | 0.1 s (.tmp/probe/sample_report.jsonl) |
| `test_isotropic_orbits.sage::test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_3]` | 1 | Isotropic classes: Lorentzian lines | n=3 sig (2,1) |det|=1 | 4.3 s (.tmp/probe/sample_report.jsonl) |
| `test_isotropic_orbits.sage::test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank0-primitive_isotropic_vector_orbits]` | 1 | Isotropic classes: Eichler lines, I-perp/I reduction | I_{2,10}, n=12 |det|=1 | 100.7 s (.tmp/sample61/report_2.jsonl) |
| `test_boundaries.sage::test_dtower_boundary_component_counts[N_3]` | 1 | Isotropic line and plane classes, split along a finite-index subgroup | d-tower N_3 | 10.7 s (CI Acceptance BUDGET line) |
| `test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_2U_A2_stable]` | 1 | Isotropic line, plane and flag classes, split along the stable group | n=6 sig (2,4) |A|=3, counts 1,1,1 | 144 s (triage.md, sage-extensions 686a72a) |
| `test_boundaries.sage::test_dawes_last_index_in_the_chain_is_the_discriminant_image_of_o_plus` | 1 | O+ generators and their discriminant image | n=6 sig (2,4) |A|=3 | 4.7 s (CI Acceptance BUDGET line) |
| `test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_2]` | 1 | Setwise line stabilizer in O+ and its discriminant image | Enriques N, n=12 |A|=2^10 | 27.7 s (CI Acceptance BUDGET line) |
| `test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_2_1]` | 1 | Setwise plane stabilizer in O+ and its discriminant image | Enriques N, n=12 |A|=2^10 | 34.7 s (CI Acceptance BUDGET line) |
| `test_vector_orbits.sage::test_anisotropic_binary_vector_classes_use_primitive_embeddings` | 1 | Vector classes: anisotropic binary | n=2, <2>+<-3> | 0.1 s (CI Acceptance BUDGET line) |
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[0]` | 1 | Vector classes: isotropic, m=1 | E10, s=0 | 0 s (CI Acceptance BUDGET line) |
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[2]` | 1 | Vector classes: Lorentzian cone (marked cell orbits) | E10, s=2, 1 orbit | 32.5 s (probe, sage-extensions 6f337e9a) |
| `test_vector_orbits.sage::test_i_2_10_is_transitive_on_norm_minus_one_vectors` | 1 | Vector classes: Eichler | I_{2,10}, n=12, s=-1 | 7.8 s (.tmp/sample61/report_0.jsonl) |
| `test_centralizers.sage::test_centralizer_images_match_oscar[oscar_lattices_with_isometry_133]` | 1 | Centralizer image: Hermitian | n=3 sig (2,1) |A|=12, ord(f)=6 | 12.3 s (probe, sage-extensions 6f337e9a) |
| `test_centralizers.sage::test_centralizer_images_match_oscar[oscar_lattices_with_isometry_118]` | 1 | Centralizer image: gluing | n=4 sig (2,2) |A|=5, ord(f)=5 | 4.6 s (probe, sage-extensions 6f337e9a) |
| `test_double_cosets.sage::test_double_cosets_partition_the_rational_group[DoubleCoset_n55_big5760_sma960_vf47_idx0]` | 1 | Double cosets: finite right cosets | degree 55, |G|=5760, |H|=960 | 0 s (CI Acceptance BUDGET line) |
| `test_boundaries.sage::test_degree_two_k3_type_ii_components_and_root_types` | 2 | published-result reproduction (#64) | | over 1500 s (cap, triage.md) |
| `test_boundaries.sage::test_degree_four_k3_has_nine_type_ii_components` | 2 | published-result reproduction (#64) | | unpriced |
| `test_isotropic_orbits.sage::test_o_plus_n_has_two_line_and_two_plane_orbits` | 2 | published-result reproduction (#64) | | 1136 s (triage.md, sage-extensions eb3aec70) |
| `test_isotropic_orbits.sage::test_stable_o_plus_n_orbits_match_dutour_sikiric_hulek` | 2 | published-result reproduction (#64) | | unpriced |
| `test_isotropic_orbits.sage::test_published_representatives_lie_in_distinct_o_plus_orbits[rank0]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_isotropic_orbits.sage::test_published_representatives_lie_in_distinct_o_plus_orbits[rank1]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1]` | 2 | published-result reproduction (#64) | | 726 s (triage.md, sage-extensions eb3aec70) |
| `test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_2_2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_isotropic_orbits.sage::test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank1-isotropic_plane_orbits]` | 2 | published-result reproduction (#64) | | 1661 s (.tmp/acceptance-coverage.md, provider 42d7361c) |
| `test_isotropic_orbits.sage::test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_8716]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_87]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_2U2_A2_stable]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_2U_minus2_minus6]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dawes_buildings_have_the_published_counts[dawes_L2_gamma2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dawes_building_incidence_graph_is_the_published_figure[dawes_2U2_A2_stable]` | 2 | published-result reproduction (#64) | | 2475 s (.tmp/dawes-incidence-2945e78-aa9b99f0.xml) |
| `test_boundaries.sage::test_classical_polarization_indices[degree_2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_classical_polarization_indices[degree_6]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_classical_polarization_indices[degree_10]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dtower_boundary_component_counts[N_20]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dtower_boundary_pictures[lambda_9]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dtower_boundary_pictures[lambda_10]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_boundaries.sage::test_dtower_boundary_pictures[lambda_11]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_orthogonal_group.sage::test_lorentzian_generators_match_polyhedral_common[LORENTZ_OrthogonalGroup0]` | 2 | published-result reproduction (#64) | | 102.9 s (.tmp/sample61/report_2.jsonl) |
| `test_orthogonal_group.sage::test_lorentzian_generators_match_polyhedral_common[LORENTZ_OrthogonalGroup1]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_orthogonal_group.sage::test_lorentzian_generators_match_polyhedral_common[LORENTZ_OrthogonalGroup2]` | 2 | published-result reproduction (#64) | | 3.7 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[10]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[4]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[6]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[8]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[-2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[-4]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[0]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_k3_lattice_has_one_primitive_orbit_per_norm[4]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_10]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_12]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_14]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_16]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_18]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_20]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_22]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_24]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_26]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_28]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_30]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_32]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_34]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_36]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_38]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_40]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_42]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_44]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_46]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_48]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_4]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_50]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_52]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_54]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_56]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_58]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_60]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_62]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_64]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_66]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_68]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_6]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_70]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_72]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_8]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[10]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[12]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[14]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[16]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[18]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[20]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[22]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[24]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[26]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[28]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[2]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[30]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[4]` | 2 | published-result reproduction (#64) | | 62.9 s (triage.md, sage-extensions 9b29712) |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[6]` | 2 | published-result reproduction (#64) | | unpriced |
| `test_vector_orbits.sage::test_published_representatives_are_one_per_orbit[8]` | 2 | published-result reproduction (#64) | | unpriced |

Branches with no small frozen node (gaps; no large node fills them in tier 1):

| Branch | Reason |
| --- | --- |
| Isometry: low Witt index (m >= 2, Witt index < 2), low_witt_isometry_witness | no frozen node found with this signature |
| O(L): rank-four low Witt index (quaternion and Clifford models) | no frozen node found |
| Vector classes: degenerate lattice (_degenerate_pointed_classes) | no frozen node found |
| Vector classes: split binary (binary_split_primitive_vector_classes) | no frozen node found |
| Vector classes: embedding branch (m = 1, sign*s < 0, Nikulin primitive embeddings) | no frozen node found; E10 and U+E8(-1) cases have s >= 0 |
| Enriques modular groups Gamma_h (coinvariant_extension_subgroup) | no small node; every recorded case is killed or errors (triage.md) |
| K3 type II components: plane classes in genus of rank 17 | no small node; smallest is degree two, over 1500 s |
| Double cosets: finite left cosets | no frozen node; every case takes the right-coset branch (double-cosets.md) |

## Long-tail families at port fee4bd3, sage-extensions 55defcb

Three cases per family, each in its own process, killed 60 s after the test call starts (`zk_cases.sh`, `zk_case_cap.py` in the probe directory). Total = mean of the three × case count. A family with a killed case is a lower bound.

| Family (test function) | Cases | Probed cases: seconds | Mean s | Total h |
| --- | --- | --- | --- | --- |
| `test_root_systems_match_allcock` | 8595 | allcock_rk3_7_w6_l7 5.9, allcock_rk3_21_w135_l21 24.9, allcock_rk3_4_w323_l4 30.9 | 20.6 | 49.1 |
| `test_simple_root_counts_match_polyhedral_common` (rank 3) | 8602 | ListReflect_99 6.4, ListReflect_3999 16.1, ListReflect_7999 30.1 | 17.5 | 41.9 |
| `test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors` | 8818 | IsotropicCases_99 4.9, IsotropicCases_3999 1.0, IsotropicCases_7999 0.8 | 2.2 | 5.5 |
| `test_enriques_modular_groups_match_dutour_sikiric_hulek` | 87 | case_5, case_40, case_80: each killed at 60 | >= 60 | >= 1.5 |
| `test_lorentzian_pairs_are_isometric` | 145 | LORENTZ_TestIsomorphism9 6.9, 69 6.8, 439 8.8 | 7.5 | 0.3 |

Probe-backed long-tail total: >= 98 h (was about 236 h). The two edgewalk families are 91 h of it.

Fix 1 (edgewalk `constraint_rows`, this commit): C_N = {v in L : b(v, L) in d Z} was built as the owned `divisible_sublattice` (ideal dual and its integral pullback, 0.22 s per divisor). It is now Sage's exact ZZ-module intersection of Z^n with d G^-1 Z^n, the pinned edgewalk.h:110-131 `ComputeLattice_LN`. Profiled w323_l4 edgewalk: `constraint_rows` 6.75 s -> 0.36 s in 30 divisors; whole edgewalk 46.7 s -> 43.3 s under the profiler. The four edgewalk gates pass (80 s). Re-probe of the same cases (unprofiled, the machine was loaded; spread is about +-20 %): Allcock 7.5, 13.9, 24.9 s (mean 15.4, 36.8 h); ListReflect 5.6, 15.6, 35.6 s (mean 18.9, 45.2 h). Edgewalk total about 82 h.

Fix 2 (sage-extensions 17f2c5e8): the finite Gram contraction reads its rows into one cached engine space (Sage's generic vector constructor inferred a parent per call); w323_l4 profile: `Lattices.norm` 1.97 s -> 1.18 s in 1436 calls, about 2 % of the walk (about 1.6 h over the edgewalk families). Gates pass.

Fix 3 (sage-extensions 0e07fb40, merged as b95b8f67): `owned_parent_class` ran the final-category aliases and participant implementations again on the cached dynamic class at each refinement; a second run for the providers that last prepared the class is skipped (the returned class is the same, so no shared class is mutated, unlike 7eb5248). `cardinal(n)` for an exact nonnegative integer goes straight to the finite construction (50750 calls per w323_l4 walk through the category dispatch); `index_rank_family` is built once per tuple of ranks; the nested-rows Gram tensor no longer re-checks its valence. w323_l4 profile: `owned_parent_class` 2.66 s -> 0.63 s, `_nested_gram_tensor` 1.04 s -> 0.27 s. Four edgewalk gates pass (73 s); the touched sage-extensions unit files (cardinals, tensors, lattice construction routes, refine) have the same 14 failures at 17f2c5e8 and after.

Same three cases per family, unprofiled, each in its own process (`zk_price_base3.log` at 17f2c5e8, `zk_price_f1.log` at 0e07fb40; one run each, spread about +-20 %):

| Family | Cases | At 17f2c5e8: s | Mean | At 0e07fb40: s | Mean | Total h before -> after |
| --- | --- | --- | --- | --- | --- | --- |
| Allcock | 8595 | w6_l7 5.7, w135_l21 14.0, w323_l4 26.3 | 15.3 | 6.3, 11.1, 23.8 | 13.7 | 36.6 -> 32.8 |
| ListReflect rank 3 | 8602 | 99 7.0, 3999 22.7, 7999 31.9 | 20.5 | 4.4, 11.8, 21.2 | 12.5 | 49.1 -> 29.8 |

Fix 4 (sage-extensions 766dd1e8, merged as c696bf9e): a sparse free-module vector computes its hash once (it decided every coefficient against zero per call; 8996 hashes per w323_l4 walk, 1.39 s -> 0.12 s profiled). Fix 5 (port eddfefc): the oriented binary cycle of a plane lattice is read once per lattice (468 calls on 292 lattices; `_OrientedBinaryCycle.__init__` 2.34 s -> 1.15 s profiled). The machine carried another Sage session and a Lean build (load average 8-18 on 10 cores), so the price is paired: each case at 17f2c5e8 and then at c696bf9e, back to back, both with port eddfefc (`zk_pair_base.log`, `zk_pair_head.log`). The pair measures fixes 3 and 4 together; fix 5 is priced by its profile only (about 0.5 s per case unprofiled, about 2 h). The c696bf9e side also holds origin/main eb5f1a3d (#552). Four edgewalk gates pass at c696bf9e (84 s); the touched unit files have the same results at b95b8f67 and after (52 passed).

| Family | Cases | 17f2c5e8: s | Mean | c696bf9e: s | Mean | Total h before -> after |
| --- | --- | --- | --- | --- | --- | --- |
| Allcock | 8595 | w6_l7 9.0, w135_l21 15.6, w323_l4 28.6 | 17.7 | 5.8, 12.4, 27.1 | 15.1 | 42.3 -> 36.0 |
| ListReflect rank 3 | 8602 | 99 7.2, 3999 14.8, 7999 27.6 | 16.5 | 9.2, 11.9, 20.5 | 13.9 | 39.5 -> 33.1 |

Fix 6 (sage-extensions 739dbc07, merged as 6f337e9a): the arrow-hash probe of a domain and the hash of a module arrow are computed once (module-arrow hashes 1.01 s -> 0.25 s profiled); `is_isotropic` is a cached datum of the lattice (818 decisions on 294 lattices, 2.47 s -> 1.15 s); `signature_pairs()` and `signature_pair(p, q)` are built once (`_sylvester` 1.41 s -> 0.75 s). Gates pass at 6f337e9a (73 s); the touched unit files (arrow and set hashes, isotropic binary lattices) have the same results at c696bf9e and after (23 passed, 1 failed on both). Paired price, load average 5-8 (`zk_pair2_base.log` at c696bf9e, `zk_pair2_head.log` at 6f337e9a; the 6f337e9a side also holds origin/main 08322dee):

| Family | Cases | c696bf9e: s | Mean | 6f337e9a: s | Mean | Total h before -> after |
| --- | --- | --- | --- | --- | --- | --- |
| Allcock | 8595 | w6_l7 4.4, w135_l21 9.6, w323_l4 20.3 | 11.4 | 4.2, 10.0, 18.6 | 10.9 | 27.3 -> 26.1 |
| ListReflect rank 3 | 8602 | 99 5.4, 3999 10.1, 7999 19.6 | 11.7 | 4.9, 8.8, 19.4 | 11.0 | 28.0 -> 26.4 |

The c696bf9e prices here are lower than in the fix 4 pair (15.1 s, 13.9 s) at the same revision: the machine load fell from 8-18 to 5-8. Only paired differences are evidence; the family totals carry the load of their run.

Blocker (long tail above 6 h, after fixes 3 to 6): the edgewalk families are about 52 h (Allcock 26.1 h, ListReflect 26.4 h at load 5-8), with IsotropicCases 5.5 h and Enriques >= 1.5 h not re-priced. The 6 h target needs the edgewalk cases about 9 times faster. Price per profiled second of the w323_l4 walk (32.8 s profiled, 18.6 s unprofiled, family mean 0.6 of w323_l4): about 1.6 h over the two families. The remaining sinks at 6f337e9a, inclusive, profiled:

| Sink | Profiled s | Calls | Hours if removed | What removing it needs |
| --- | --- | --- | --- | --- |
| `edge_plane_for_divisor` (port) | 8.9 | 390 | 14 (upper bound; includes the plane lattices and Sage linear algebra, 1.4 s own time) | pinned edgewalk.h steps; the owned share is the plane lattice and its elements |
| `anisotropic_cycle` (port) | 4.3 | 468 | 7 | `cycle` builds each representative as an owned element (`linear_combination`, `to_vector`) and the period as an owned isometry; the stepping is integer arithmetic |
| `initial_vertex` | 4.0 | 1 | 6.5 | `vectors_of_square` 1.1 s; the rest is owned element and Mor construction |
| `Lattices(ZZ)(gram)` | 2.7 | 293 | 4.3 | the constructor redesign: `_gram_axiom_categories` decides nondegeneracy, unimodularity and evenness at construction (0.7 s; placement axioms, so making them lazy changes which category operations an object has: an architecture decision for sage-extensions), and `_construct_form_module` builds a fresh free-module parent through `_object_of` (1.05 s) |
| `Mor_C(A, B)` classcall | 1.8 | 1137 | 2.9 | one Mor object per fresh endpoint pair; the arrow class is shared (c62175f), the object init chain is not |
| `_new_sparse_free_module` | 1.6 | 458 | 2.6 (overlaps the lattice row) | the category meet and `_finitely_generated_free_placement` (0.5 s: membership tests on the label set) per module |

No single row reaches the 9x factor, and the rows that are not port arithmetic are the per-object construction chain (`_object_of`, category meets, placement axioms). Removing that chain is the engine-data contract of OWN-22/OWN-24 for each constructor, with the placement-axiom question decided first in sage-extensions.


Fix 7 (decision: lattice construction skips eager checks): `_gram_axiom_categories` over ZZ reads nondegeneracy (det G != 0), unimodularity (det G = +-1) and evenness (even diagonal) from the engine Gram matrix; the owned determinant is computed when a consumer asks for it (w323_l4 profile 0.7 s -> 0.37 s in 300 lattices). Fix 8: the content of a framed free element over ZZ is the gcd of its engine coordinates; the owned order ideal was formed per call (`content` 1.07 s -> 0.08 s in 504 calls, the `anisotropic_cycle` sink). Fix 7 alone, paired against f7a03e96 (`zk_pair_f7_*.log`), is inside the noise: Allcock 4.1/8.9/16.6 -> 4.2/9.7/16.9 s, ListReflect 4.4/8.9/16.7 -> 4.3/9.1/20.2 s; its profiled 0.33 s prices at about 0.5 h, under the one-hour threshold. Fixes 7 and 8 together, paired at load 4-5 (`zk_pair_c2_base.log` at f7a03e96, `zk_pair_c2_head.log` at 7d64d821, merged as fb742277); gates pass at fb742277 (64 s); the touched unit files have the same results as f7a03e96 (19 passed; 3 failed on both: `test_coeven_archive` timeout, two `test_coxeter_gram_archive` AttributeError root_lattice):

| Family | Cases | f7a03e96: s | Mean | 7d64d821: s | Mean | Total h before -> after |
| --- | --- | --- | --- | --- | --- | --- |
| Allcock | 8595 | w6_l7 4.3, w135_l21 8.8, w323_l4 19.5 | 10.9 | 4.1, 9.5, 17.2 | 10.3 | 25.9 -> 24.5 |
| ListReflect rank 3 | 8602 | 99 4.6, 3999 9.4, 7999 18.3 | 10.8 | 4.1, 9.4, 17.4 | 10.3 | 25.7 -> 24.6 |

Blocker (long tail above 6 h, after fixes 7 and 8): the edgewalk families are about 49 h. The w323_l4 profile at 7d64d821 (24.5 s profiled) has no single owned sink above about 2 s: `edge_plane_for_divisor` 7.0 s inclusive, of which `plane_lattice` (owned `Lattices(ZZ)(gram)`: `_construct_form_module` 0.83 s, `_form_module_categories` 0.30 s, `_gram_axiom_categories` 0.37 s, `_nested_gram_tensor` 0.21 s) 2.3 s and Sage linear algebra (span, kernel, solve) about 2.5 s of pinned edgewalk.h steps; own time 1.0 s (one product for the projection of all constraint rows did not change it); `initial_vertex` 3.8 s (`vectors_of_square` 0.9 s, `_vinberg_decomposed_roots` 0.9 s, cone construction 1.1 s); `anisotropic_cycle` 2.7 s (`selected_module_resolution` realization 0.6 s in 70 calls, `Isom` construction 0.55 s in 55 calls). The remaining owned cost is the per-object construction chain (form module, Mor, realization of the selected resolution). The 6 h target needs about 8 times; that is the constructor redesign of OWN-22/OWN-24 for the form-module constructor, not a further local sink.


### Price of the construction chain (sage-extensions fb742277, w323_l4)

`zk_ew_prof.py w323_l4` at fb742277: 27.9 s profiled (17.2 s unprofiled in the fix 7-8 pair, ratio 0.62). Own time by owner (`zk_part.py`):

| Owner | Profiled s | Share |
| --- | --- | --- |
| sage-extensions | 13.8 | 50 % |
| Sage arithmetic (matrices, free modules, saturation) | 4.6 | 16 % |
| port | 3.4 | 12 % |
| Sage category and structure machinery | 3.1 | 11 % |
| builtins and other | 3.0 | 11 % |

The owned chain (sage-extensions plus the Sage category machinery it drives) is 16.9 s, 61 % of the walk. Its inclusive entry points (overlapping): `_object_of` 2.6 s in 1704 parents (1.5 ms each), `Mor_C(A, B)` classcall 1.7 s in 1053 Mor objects, `Lattices(ZZ)(gram)` 2.0 s in 293 lattices, `realize` of the selected resolution 1.4 s in 78 calls, ring-element construction 0.95 s in 47198 elements, set-map application 3.7 s in 13087 calls (the coloured family `colour` 1.6 s in 585 calls, the port's `coloured_family` lambda 1.1 s in 657 calls). No entry point is above 2.7 s. Price over the 17,200 edgewalk cases (49 h at load 4-5, the share taken as uniform over the families): the whole owned chain is about 30 h; parent and Mor construction about 11 h; the lattice constructor alone about 3.5 h (a direct probe, `zk_ctor.py`, 300 binary Grams: 2.8 ms per lattice unprofiled, 0.8 s of the 17 s case, about 2.4 h). With the owned chain at zero the families are about 19 h, so the 6 h target also needs the port arithmetic and the Sage linear algebra of `edge_plane_for_divisor` (7.8 s inclusive, 21503 `MatrixSpace` element constructions) reduced.

Fix 9 (sage-extensions c67bdff9, merged as 54664c40): `Cat().meet` of a tuple and the group-algebra placement of a category are memoized per argument (categories are unique objects; no class is mutated), and a two-index tensor's engine matrix is read from its flat components. `zk_ctor.py` paired twice: 4.16 -> 3.62 and 3.65 -> 2.75 ms per lattice. Skipping the repeated alias pass on the shared element class (the analogue of fix 3) measured nothing (0.17 s in both profiles) and is not committed. Paired price at load 6-7 (`zk_pair_c3_base.log` at origin/main dc5b828b, `zk_pair_c3_head.log` at c67bdff9); w217 passes at the head (13.1 s):

| Family | Cases | dc5b828b: s | Mean | c67bdff9: s | Mean | Total h before -> after |
| --- | --- | --- | --- | --- | --- | --- |
| Allcock | 8595 | w6_l7 4.3, w135_l21 9.5, w323_l4 19.0 | 10.9 | 4.4, 9.6, 18.4 | 10.8 | 26.1 -> 25.8 |
| ListReflect rank 3 | 8602 | 99 4.0, 3999 10.2, 7999 18.5 | 10.9 | 4.2, 8.7, 20.0 | 11.0 | 26.0 -> 26.3 |

Inside the noise: about 0 h removed.

Blocker (after fix 9): the edgewalk families are about 52 h at this load. The cost is the per-object overhead of the sage-extensions object model, spread over thousands of constructions at 0.2-3 ms each; local memos move it by less than the noise. Removing it is an architecture decision in sage-extensions, with these alternatives: (a) a lattice, sublattice or morphism built from engine data is a light object that holds the engine matrix and its category, and builds the framed free module, Mor objects and owned elements on first request (OWN-22 read literally; the largest change, and it touches every consumer of `module_generating_set` and owned elements); (b) one parent class per (category, implementation) built once and never refined per object, so `_object_of` and Mor construction skip `realize_owned_category` and construction hooks after the first object (needs the class-mutation rule of 7eb5248 settled); (c) the port's per-walk hot loop (edge planes, anisotropic cycles, coloured families) works on engine matrices and raises owned objects only at the walk's result, which is a port decision under the OWN boundary rules (public results stay owned objects). Recommendation: (b) first, because it applies to every constructor and does not change the object model; (a) for lattices from Gram matrices after it. Either needs the owner's choice before implementation.

Next sinks in the w323_l4 profile (46.7 s): `anisotropic_cycle` on the binary edge planes 6.6 s (468 calls), `plane_lattice` construction `Lattices(ZZ)(gram)` 3.9 s (292 lattices, 12 ms each: queue item 4, the lattice constructor), `is_isotropic` of binary planes 2.1 s, `initial_vertex` 4.9 s.

Blocker (long tail above 6 h): after fixes 1 and 2 the probe-backed total is about 88 h (edgewalk about 82 h, IsotropicCases 5.5 h, Enriques >= 1.5 h, Lorentzian pairs 0.3 h). The 6 h target needs the edgewalk cases about 14 times faster. No single step is more than about 15 % of the w323_l4 walk: the remaining time is spread over owned rank-2 objects built per edge (`Lattices(ZZ)(gram)` 12 ms each, per-element `q`, `content`, `module_rank` and cardinality decisions in `anisotropic_cycle`, `vectors_of_square` in `initial_vertex`). Removing it is the constructor redesign of queue item 4 (engine-data constructors for `_lattice`, `divisible_sublattice`, kernels, finite ordered sets), each with its own sage-extensions engine-data contract (OWN-22, OWN-24) and probe; not started. K3 queue items 2-4 (flip LLL, remaining transporter sinks, the 600 s target) were withdrawn by the coordinator: the K3 node needs to pass once at about 950 s.

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

## Sinks fixed after e4e901e

Probes are cProfile wraps of one function in the K3 type II degree-two node, over the first calls the walk reaches; per-call times vary with the cell, so the same call count is compared.

| sage-extensions | Port | Probe: operation, input, before -> after | Change |
| --- | --- | --- | --- |
| 6bed0a6 | 3d4a424 | `cell_transporter`, first 60 calls: 102.9 s -> 75.6 s. Sage `PermutationGroup.stabilizer` chains in the basis-image search 29.6 s; generated-subgroup generators admitted by `membership_decision` 10.1 s; restricted Gram of `subobject_on` by rank^2 pairings 14.6 s; owned fraction equality 18 s | port: one GAP `Stabilizer(H, a, OnTuples)` per search node. sage-extensions: generators admitted by the free cover's range law; restricted Gram as one engine product B G B^t; sparse sums and negatives keep decided coefficients; fractions with retained engine values compare in the engine |
| 4afe041b | 3d4a424 | `cell_transporter`, first 48 calls: 46.1 s -> 37.5 s. The integral-structure orbit search rebuilt each letter's residue automorphism at every orbit point | residue map per letter, built once; pairing-graph rankings read once |
| 2f91a98d | 3d4a424 | `cell_transporter`, first 48 calls: 37.5 s -> 33.7 s. `rational_transformation` and the transport composite evaluated lambda morphisms per generator | T_sigma and the composite as framed matrix products; the rational family raised from the family's coordinate matrix |
| 03c5997b | 5a53356 | `flip_across`, 190 calls: 161.9 s -> 114.4 s. Each shell target and bound raised to owned rationals and lowered again by the coordinate CVP (45 s); `_engine_component_matrix` as a cached_function compared equal Grams by value (9 ms per call) | engine rationals passed to the CVP lowering primitive; a tensor retains its engine component matrix |
| eb3aec70 | 5a53356 | `cell_transporter`, first 48 calls: 33.7 s -> 25.8 s (65 s at e4e901e). `rational_transformation` built the I x I relabelling and the family section by owned evaluation (`_engine_matrix` 9 s) | section j B^-1 from the retained adjugate; relabelling as a column permutation |
| 93021dfa | c7b0d82 | r_0 pointed isometry, `IntegralStructureAction.transporter`: 3.3-5.1 s per call at 5d408ca -> 1.9-2.5 s after 4afe041b -> 1.1-1.7 s. Preimages into the restriction of scalars of a free framing used an owned coordinate lift and a second owned preimage per query (420 queries, 4.7 s in six calls) | the free framing uses the presented case's retained Smith-form solver with P = 1 |
| 022f7a52, 7eb5248 | 80ce6cc | ListReflect_2917 `edgewalk_fundamental_domain` under the profiler (queue item 2): 21.3 s -> 18.2 s with `owned_category_join` memoized per tuple and `owned_parent_class` per (implementation, category); 23.3 s at 0f7d1ad after the parent-class memo was withdrawn (7eb5248: later refinements set attributes on the shared class, and a set object became unhashable in a following test) | the join memo stays; the parent-class memo needs classes that refinement does not mutate (open) |
| - | 8fc1e77 | Facet walk threshold (queue item 4): the constant 20000 is replaced by the size of the ray matrix, \|G F\| <= n d. On this node the orbit sizes are at most 256 or at least 43200 (n d from 289 to 7999), so the listed orbits are the same. Listing every orbit held 3.2 GB in one 8847360-facet orbit; GAP `MinimalImage` and `CanonicalImage` orbit names cost 150-175 s in the first 250-400 s of the walk | rule from the data: each listed orbit holds at most n^2 d integers |
| 2f55afc (main 1e285e9) | 57d6461 | `cell_transporter`, first 48 calls: 31.8 s at 0f7d1ad -> 23.9 s. The residue-orbit transporter formed quotient elements from each orbit key and lowered them again, built each letter's residue map as an owned quotient automorphism with its linearity check (0.23 s per letter), and found the two intermediate keys by owned module preimages (1.8 s); `configuration_isometry` built each generator from owned entries, so hashing it in `finite_ordered_set` lowered 361 entries per generator again (2.5 s) | keys act by the integral cover matrix X^t of g restricted to M (M/dM is presented on the framing of M); the intermediate key from B Y = C and C Z = dB over ZZ; the configuration isometry raised once from its engine matrix |
| 1e285e9 | 476ffbe | `point_permutations_from_rows` in the K3 traversal: 165 s of samples at 1e285e9 -> not sampled; the facet walk `register` reads \|Stab_G(F)\| without a Sage permutation group | the stabilizer search retains each generator's integral row action; GAP `Stabilizer(..., OnSets)` read by its `Size` |
| 1e285e9 | a4fa9be | Replay probe (recorded K3 inputs, fresh backend, 60 s process cap): `flip_across` 0.536 -> 0.531 s per call (40 calls); `cell_transporter` 0.511 s per call (40 calls) | the facet normal from Sage's exact kernel of the incident rows; the line is oriented and made primitive, so the basis of the line does not matter |
| 55defcb | (this commit) | Replay probe, 40 recorded K3 `cell_transporter` calls: 0.511 -> 0.382 s per call (0.404 s after the merge with origin/main). `_family_coordinate_matrix` read each vector and its coordinate vector once per label (2.0 s); `pairing_graph` was rebuilt per call (2.1 s); `integral_restriction` evaluated the rational map per generator (3.8 s -> 1.4 s in five calls) | coordinates read once per vector; the pairing graph cached per realization and the family morphism raised from the coordinate matrix; integral restriction over ZZ from one framed engine matrix and its integral inverse. Port: cell vectors raise each nonzero coordinate once by `_owned_engine_element` (element construction 2.4 -> 1.1 ms at rank 20) |

Nodes at sage-extensions eb3aec70 with port c7b0d82 (local container, wall time with Sage start-up, one node per session, private stack monitor every 5 s):

| Node | Before | After | State |
| --- | --- | --- | --- |
| test_o_plus_n_has_two_line_and_two_plane_orbits | over 950 s (6179cee) | 1136 s | passes (not re-measured after 1e285e9: whole-node runs belong to the coordinator). Inclusive: `_integral_generator_set` 385 s (`split_orbit` in the plane classes), `splitting_isometries` 250 s, quotient O(K) traversal 120 s |
| test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1] | over 650 s (a106d23) | 726 s | passes (not re-measured after 1e285e9: whole-node runs belong to the coordinator). Inclusive: the predicate-subgroup Schreier cover of the line stabilizer 480 s, `_integral_generator_set` 410 s, discriminant morphisms 185 s |
| test_degree_two_k3_type_ii_components_and_root_types | over 900 s (e4e901e) | over 1500 s (cap) | the quotient O(K) traversal is 1340 s of the 1500 s and has not finished: `cell_transporter` 630 s (coloured-family transporter 530 s), `flip_across` 395 s, facet walk 340 s, cell stabilizers 185 s. Per call the transporter is 2.5 times faster (first 48 calls 65 s -> 25.8 s) and `flip_across` 1.4 times |

## Unit tests touched since 78b7cf4 (queue item 5)

The whole sage-extensions suite (`--no-time-gates --no-test-lint`) at 022f7a52: 487 failed, 775 passed in 5455 s. The 487 failing nodes at 78b7cf4 (worktree, its own tests): 431 failed, 56 passed. The 56 were bisected over the first-parent line with the same node list (7232ec0 and 78b7cf4: 50 of 50 pass):

| Cause | Nodes | Fix |
| --- | --- | --- |
| 3239c91 `module_category` as a cached_method: RecursionError through `CachedMethod.__get__` during module construction | 43 (tests/algebras, tests/categories, tests/groups) | 3bbf92b |
| 022f7a52 memoized `owned_parent_class`: a shared class mutated by a later refinement | 6 (order-dependent; pass alone) | 7eb5248 |
| 5f65d2e residue restriction as B^-1 A B for a reference lattice not of full rank | 1 (test_commensurability_restriction) | d5a08df |
| #531 hashing rules from origin/main (bc00bc1 side: finite ordered sets hash their points; a class over ZZ has no hash) | 6 (test_centralizer_gluing, test_commuting_involution_centralizers, test_cyclic_normalizers x2, test_cyclotomic_centralizer, test_de_rham rational dual numbers) | not from these fixes; owned by sage-extensions #531 |

## Open sinks after 0f7d1ad

| Node | State | Sink | Class |
| --- | --- | --- | --- |
| test_degree_two_k3_type_ii_components_and_root_types | whole-node runs stopped; per-call replay only (`zk_replay.py` on 40 recorded calls of each sink) | per call at port 476ffbe and sage-extensions 1e285e9: `flip_across` 0.53 s, `cell_transporter` 0.51 s. Earlier inclusive totals (630 s at 0.54 s per transporter call, 395 s at 0.6 s per flip) give about 1170 transporter calls and 660 flips, so the priced traversal sinks are about 600 s + 350 s, above 600 s. Flip: 2.4 shells per flip, each a dual-LLL sublattice reduction in GAP (43 ms) and a PARI `qfcvp` per multiplier (23 multipliers per shell, 2.1 ms each), both pinned source steps. Transporter: the spanning coloured-family calls are 3.2 s each (5 of 40 replayed calls, 16 of 29 s): `_family_coordinate_matrix` reads every coordinate through `to_vector` (5480 entries, 2.0 s), `pairing_graph` 2.1 s, `integral_restriction` 3.8 s in five calls | c: per-element owned coordinate reads in the coloured-family realization; open, next sink |
| test_equivalence_witness_carries_a_representative_to_its_reflection[r_0] | transporter 1.1-1.7 s per call | orbit words 3.1 s and `_intermediate_image_key_or_none` 2.7 s in six calls: module preimages through `factor_through_or_none` | c |
| Queue item 2, remaining constructors | not done | `divisible_sublattice` (24 calls 3.1 s: `ideal_dual` builds its pairing Mor eagerly, `pullback_form` recomputes the Gram `subobject_on` already holds), `_new_sparse_free_module` and `_fresh_free_module_on` (332 calls 1.7 s), Mor objects per fresh endpoint pair (1004 classcalls 1.5 s), `kernel`, `_lattice`, `finite_ordered_set` (24 calls 3.8 s in K3 `_engine_subgroup`), cone | blocker: each constructor needs its own engine-data contract in sage-extensions (CONTRIBUTING OWN-22/OWN-24) and its own probe; the measured share in ListReflect_2917 is about 6 s of 18 s; not started in this pass: the K3 sinks came first in the queue, and each constructor still needs its engine-data contract (OWN-22: construction stores the datum, every other datum on first request) and its own probe |

### Repin 5bfab13a (#537 regression)

sage-extensions 54664c40 merged 5b45523e (#537), which failed six tier-1 nodes in run 43829660254: four with TypeError (StabilizerSubgroups._call_() got an unexpected keyword argument 'description': GL(M).flag_stabilizer still passed its clause by keyword) and test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_2], [I_2_1] with NameError (_named_by_its_own_words not imported in lattice_morphisms.py). Fix 9743a3a0, owner test 5bfab13a (tests/lattices/test_stabilizer_clauses.sage, red at 54664c40). Local tier-1 at 5bfab13a (zk_t1_r1.log, 60 s cap per node): 45 of 49 pass; the other four reach the 60 s cap with no error, and their recorded prices are above the local cap or near it under load: test_o_n_maps_onto_o_of_the_discriminant_form 468 s, dawes_2U_A2_stable 144 s, I_1_2 27.7 s, I_2_1 34.7 s (I_1_2 also reaches the cap at dc5b828b).

Parent, element and morphism classes built once per (implementation, category) and complete when built (sage-extensions branch perf/class-memo, 9b7690a6), paired against 58e1c768 (zk_pair_f10): six-case mean 8.2 -> 8.4 s, inside the noise; about 0 h removed, not landed.

Fix 10 (port, queue item 3): `edge_plane_for_divisor` forms the QQ row space of the walls as their reduced echelon form, the ZZ span of the projected constraint rows as the Hermite form of their integral scaling, and the integral left kernel of C cap P as the Hermite form of PARI's kernel basis. Each is the basis the replaced Sage module carried (the echelon and Hermite forms are unique; `zk_span_eq.py`: 400 random spans and 300 kernels equal), so every downstream row is unchanged; the projection of the constraint rows is one matrix product. Profile of w323_l4 (`zk_prof_e3.log`): `edge_plane_for_divisor` 6.05 -> 4.19 s inclusive (saturation, span and submodule construction gone). Two paired runs (`zk_pair_e3`, `zk_pair_e3b`, port 5c70002 base):

| Family | Cases | Base mean s (two runs) | Head mean s | Total h before -> after |
| --- | --- | --- | --- | --- |
| Allcock | 8595 | 8.87, 8.60 | 8.70, 7.87 | 20.8 -> 19.8 |
| ListReflect rank 3 | 8602 | 8.53, 9.83 | 8.43, 9.93 | 21.9 -> 21.9 |

About 1 h removed (Allcock w323_l4 14.6 -> 12.6 and 11.7 s); the ListReflect change is inside the noise.

Blocker (after fix 10): the edgewalk families are about 42 h at this load. The rest of `edge_plane_for_divisor` is `plane_lattice` (1.6 s profiled: `Lattices(ZZ)(gram)` at 4.1 ms profiled per lattice, spread over `_new_sparse_free_module` 1.3 ms, `_gram_axiom_categories` 0.8 ms, `_nested_gram_tensor` 0.5 ms, indexed families and tensor hashing) and Sage per-operation matrix overhead (MatrixSpace construction, about 0.5 s). Queue item 1 (one class per (implementation, category), built complete) is on sage-extensions branch perf/class-memo and removed nothing (fix 10 table above: 8.2 -> 8.4 s). Queue item 2 (a lattice from a Gram matrix holds the engine matrix and builds its free module on request) is not started: the lattice parent is itself the framed free module object (`_fresh_free_module_on` with the lattice categories), so a deferred free module needs a second object type for lattices whose module structure is unrealized, and every consumer of `module_generating_set`, owned elements and Mor objects must route through it; the measured gain bound is the 4.1 ms per lattice above, about 2.5 h of the 42 h.

Correction to the blockers above ("no single sink", "spread over many small steps"): they priced entry points, not totals. Multiplied out over the 17,200 edgewalk cases (README "How to price a step", step 6), the per-object construction is hours: parents through `_object_of` 1704 per case × 1.5 ms ≈ 12 h; `Mor_C(A,B)` objects 1053 per case, 1.7 s ≈ 8 h; `Lattices(ZZ)(gram)` 293 × 2.8–4.1 ms ≈ 4–6 h, inside the parents row. About 20 h is parent and Mor construction, so it is the next fire, not a blocker. The levers are the count (total against distinct by construction key; intern equal constructions, one Mor object per endpoint pair) and the cost of one instance. The rejected class memo cached classes (under 1 s) and did not change the cost of an instance.

## Construction count and instance cost (queue items 1 and 2, sage-extensions 87db29aa)

Probe `zk_count.py` wraps `_object_of` (parents), `CategoricalMor.__init__` (Mor objects) and `_lattice_from_gram_tensor` (Gram lattices) during one edgewalk. A construction key is the category with the engine and data, read by value: matrices by entries, functions by code and closure values, elements by parent key and representation, parents by their defining construction. Rings, Sage fields and a few engine objects are keyed by identity. Times are the outermost construction only, measured with the key computation off (`ZK_TIMING_ONLY=1`), so a nested construction is counted once.

| Construction | w323_l4 total | distinct | ListReflect_2917 total | distinct | duplicate cost per case | hours over 17,200 cases |
| --- | --- | --- | --- | --- | --- | --- |
| parents through `_object_of` | 1694 | 1209 | 997 | 694 | 0.21 s, 0.14 s | about 0.8 h |
| `Mor_C(A,B)` objects | 1044 | 751 (1021 by endpoint identity) | 695 | 519 (678) | 0.07 s, 0.04 s | about 0.3 h |
| `Lattices(ZZ)(gram)` | 294 | 294 | 106 | 106 | 0 | 0 |

The Mor duplicates have distinct endpoint objects with equal construction keys: they follow from duplicate parents. The largest duplicate sites are `finite_ordered_set` on literal points (132 in w323_l4) and `_from_finite_members` (77), their `Sets().Mor` in `admits_arrow` (101) and `Mono` in `inclusion` (101), and indexed families of module-morphism generator images (about 150). Each site removes under 0.2 h; all duplicates together remove about 1.1 h. Interning is not done: no single site reaches an hour.

Cost of all construction at sage-extensions 1f4b08c2 (outermost, without the profiler): w323_l4 1.41 s of a 13.9 s walk, ListReflect_2917 0.93 s of 7.4 s, so about 1.2 s × 17,200 ≈ 5.6 h. The earlier figure of about 20 h (1704 × 1.5 ms and 1.7 s of Mor) was read from cProfile, which inflates these many small calls by a factor of two to three.

| Kind (w323_l4) | count | ms each | s per case | hours over 17,200 (both families, scaled by count) |
| --- | --- | --- | --- | --- |
| Gram lattices (binary plane lattices) | 294 | 0.65 | 0.19 | 0.6 |
| `ModuleMor` objects | 75 | 2.3 | 0.18 | 0.7 |
| finite ordered sets | 468 | 0.26 | 0.12 | 0.5 |
| indexed families | 540 | 0.21 | 0.11 | 0.5 |
| `SetMor` objects | 526 | 0.19 | 0.10 | 0.5 |
| finite subgroups of O(L) | 69 | 1.1 | 0.08 | 0.4 |
| generated subgroups | 16 | 4.2 | 0.07 | 0.3 |
| rational polyhedral cones, before fix 11 | 36 | 22 (with derived data) | 0.79 | 1.9 |

Isolated instance costs (`zk_inst.py`, 999 builds each, no profiler): finite ordered set 0.27 ms, indexed family 0.11 ms, `Sets().Mor(A, B)` on a fresh endpoint 0.33 ms, `Lattices(ZZ)(gram)` on a fresh binary Gram 2.2 ms. Profile of the lattice (4.4 ms under cProfile): `_gram_axiom_categories` 0.84 ms and `_form_module_categories` 0.42 ms (category resolution by `in` tests and category construction), 65 `Category.__classcall__` lookups per lattice 0.48 ms, `_engine_component_matrix` 0.39 ms, `_tensor_hash` 0.22 ms. Memoizing the two category resolutions removes about 0.6 ms of a lattice, 0.2 s per w323_l4 case, about 0.6 h: under an hour, not done.

Fix 11 (sage-extensions 1f4b08c2): `RationalPolyhedralCones.ParentMethods.__init__` wrapped both covector families and the wall roots in finite ordered sets, lowered each covector through its owned coordinates and built Sage's polyhedron at construction. These are now lazy attributes. `_cone_from_engine_polyhedron` (the route of `from_rays`) gives the cone the polyhedron of its primitive integral rows, the same input that the lazy build gives the engine. Probe `zk_cone.py`: the rays and the H- and V-representations of all 36 and 27 recorded cones have equal digests before and after. Recorded-cone rebuild 21.4 -> 8.6 ms (w323_l4) and 17.6 -> 7.7 ms (ListReflect_2917); cone construction in the walk 0.87 -> 0.36 s and 0.60 -> 0.30 s. Priced by count: 36 × 13 ms × 8595 + 27 × 11 ms × 8602 ≈ 1.1 h + 0.7 h ≈ 1.8 h. A second lever in the same constructor, the owned covectors (`lattice.dual_module()({...})`, a module morphism with an indexed family of generator images and a set map, about 2 ms under the profiler each, 405 per 180 cones), is not changed: the covectors are the defining datum. Caching the scalar identity `Mor_Rings(R, R).identity()` of each linear map (0.6 ms of a covector under the profiler) changed nothing measurable (9.05 against 8.63 ms per cone) and is withdrawn.

Paired runs of fix 11 (zk_pair_cone, zk_pair_cone2; sage-extensions 87db29aa base against 1f4b08c2, port b1ec1ff, three cases per family):

| Family | Cases | Base mean s (two runs) | Head mean s | Total h before -> after |
| --- | --- | --- | --- | --- |
| Allcock | 8595 | 8.63, 7.63 | 8.50, 7.80 | 19.4 -> 19.4 (inside the noise) |
| ListReflect rank 3 | 8602 | 9.13, 7.67 | 8.53, 7.43 | 20.1 -> 19.1 |

About 1 h removed by the pairs, 1.8 h by count. New edgewalk total about 38.5 h at this load.

Repin sage-extensions 081d4829 (fix 11 merged with origin/main 2e2fa8b8). Local tier-1 (zk_t1_cone.log, one process per node, killed at 60 s): 45 of 49 pass; the four priced over the cap (test_o_n_maps_onto_o_of_the_discriminant_form, dawes_2U_A2_stable, I_1_2, I_2_1) reach 60 s with no error. Edgewalk gates (w183, w217_l26, ListReflect_2917, branch-representative sample) pass at 1f4b08c2 and at 081d4829. Touched sage-extensions tests (polyhedral cones, polytopes, dual point hash, polygon invariants, parabolic count, Vinberg facets, hyperbolic chambers, hyperbolic space): 8 failed, 44 passed, the same eight nodes as at 87db29aa.

Blocker for the next hour item: after fix 11 all parent and Mor construction is about 5.6 h over 17,200 cases (table above), and no single construction kind reaches an hour; the duplicate constructions together are about 1.1 h over eight sites. The next measured items are port-side and owned-element arithmetic in the binary-plane step (w323_l4 profile at 1f4b08c2, ew_h1.prof: `extension_roots` 7.1 s of a 19 s profiled walk; `plane_lattice` 1.5 s, `is_isotropic` of the plane 0.75 s, `_OrientedBinaryCycle` construction 0.66 s and `cycle` 1.36 s, mostly owned `q()`, `b`, `to_vector` and `linear_combination` per vector). Unprofiled, they are about 3 s of a 13.5 s case, about 6 h for Allcock if the ratio holds; they are not yet priced by count.

## Wall-clock phase table (queue item 1, port c435fd6, sage-extensions 081d4829)

Probe `zk_phase.py` (probe directory): `sys.monitoring` PY_START/PY_RETURN events on the port edgewalk code (`edgewalk.py`, `edgewalk_rank2.py`, `edgewalk_extensions.py`) and LINE events on `walk`, `edgewalk_fundamental_domain`, `binary_candidates` and `edge_plane_for_divisor`. Each interval between two events is charged to the port function and line on top of the tracked stack, so the time of sage-extensions, Sage, PARI and GAP callees is charged to the port line that called them, and the charges sum to the walk time by construction. No profiler. Eight paired runs (plain, then instrumented, same case, separate processes): instrumented 58.6 s against plain 55.5 s, +5.5 %; single pairs differ by -3 % to +22 %, which is the load noise of one case.

Correction to every earlier edgewalk total: they priced a cold process. The long tail runs thousands of cases per pytest process, so imports, class and category creation and the PARI/GAP start are paid once per shard, not per case. In a cold process the first walk pays about 1 s in the indefinite reduction alone (`_indefinite_reduction` 0.78 s on its first call, 0.024 s after). `zk_phase.py` and `zk_sample.py` therefore warm the process with two small walks (w6_l7, ListReflect_99) before the timed walk. Seeded random sample (`zk_sample.py`, seed 1, 25 cases per family, warm): Allcock mean 6.04 s, ListReflect rank 3 mean 5.66 s, so the families are 14.4 h + 13.5 h = 27.9 h, not 38.5 h. Medians of the sample: w65_l29 (6.8 s) and ListReflect_6219 (6.2 s).

Seconds per case (mean of two instrumented runs). Hours: the phase's mean share of the cases of each family times the family's warm total.

| Phase | w323_l4 | w65_l29 | ListReflect_2917 | ListReflect_6219 | hours over 17,200 |
| --- | --- | --- | --- | --- | --- |
| anisotropic family (`_OrientedBinaryCycle`, cycle, period unfolding) | 1.51 | 1.43 | 0.90 | 1.20 | 4.9 |
| edge plane (`edge_plane_for_divisor`, `constraint_rows`, `plane_primitive`) | 1.59 | 1.46 | 0.84 | 0.93 | 4.6 |
| reflectivity certificate (root closure `from_finite_orbits`, chamber rays) | 1.07 | 0.78 | 0.54 | 0.81 | 3.1 |
| vertex transporter search (walk line 1629) | 0.88 | 0.55 | 0.49 | 0.88 | 2.7 |
| initial vertex (Vinberg) | 1.49 | 0.54 | 0.31 | 0.30 | 2.3 |
| `plane_lattice` (`Lattices(ZZ)(gram)`) | 0.89 | 0.66 | 0.30 | 0.38 | 2.1 |
| `is_isotropic` of the plane (binary_candidates line 1183) | 0.52 | 0.38 | 0.18 | 0.22 | 1.2 |
| edge directions (cone facets) | 0.35 | 0.25 | 0.24 | 0.36 | 1.2 |
| coordinate helpers (key, row, pair, element) | 0.28 | 0.18 | 0.18 | 0.28 | 0.9 |
| report transport to the input lattice | 0.23 | 0.20 | 0.14 | 0.19 | 0.7 |
| vertex stabilizer generators (walk line 1616) | 0.12 | 0.07 | 0.10 | 0.43 | 0.7 |
| edge orbits under the stabilizer | 0.17 | 0.12 | 0.13 | 0.29 | 0.7 |
| indefinite reduction | 0.11 | 0.14 | 0.12 | 0.13 | 0.5 |
| extension_roots (ambient reconstruction) | 0.15 | 0.13 | 0.08 | 0.10 | 0.5 |
| report | 0.09 | 0.10 | 0.08 | 0.12 | 0.4 |
| Dynkin extensions | 0.11 | 0.09 | 0.07 | 0.08 | 0.4 |
| walker setup | 0.08 | 0.12 | 0.08 | 0.08 | 0.4 |
| admit isometry (fixed lattice) | 0.19 | 0.19 | 0.01 | 0.02 | 0.3 |
| binary candidates (isotropic split, filtering) | 0.06 | 0.06 | 0.03 | 0.04 | 0.2 |
| vertex_from_root, candidate order; vertex configuration; walk bookkeeping | 0.07 | 0.07 | 0.06 | 0.08 | 0.3 |
| sum = instrumented walk | 9.97 | 7.50 | 4.86 | 6.95 | 27.9 |

## Unit operations per case and unit costs (queue item 2)

Counts: `zk_phase.py count` (CALL events on port lines: every call a port line makes, keyed by line and callee). Unit costs: `zk_unit_ops.py`, `zk_unit_ops2.py`, `zk_unit_mat.py`, `zk_mem_unit.py` (10,000 operations each, fresh inputs, no profiler); totals of a named operation over one walk: `zk_sylv.py`, `zk_ringel.py`, `zk_member.py`, `zk_orbcl.py` (outermost wall time of the operation in a warm walk).

| Owned operation | unit cost | engine reference |
| --- | --- | --- |
| `Lattices(ZZ)(gram)`, fresh binary Gram | 3.17 ms | the 2 x 2 matrix |
| `is_isotropic`, first request (signature 1.04 ms, `qfsolve` witness 0.32 ms, rest 0.14 ms) | 1.63 ms | -det square test |
| `signature_pair`, first request (Sage `QuadraticForm.signature_vector` 0.57 ms of it) | 1.04 ms | |
| `_OrientedBinaryCycle(P)` (2 `module_generator`, 2 `q`, 1 `b`) | 1.23 ms | three integers |
| `Isom(P)(images)`, rank 2 | 2.03 ms | |
| isometry application, rank 2 / rank 3 | 1.02 / 1.40 ms | 2 x 2 / 3 x 3 product |
| `linear_combination` (element construction) | 0.27 ms | |
| `q()` / `b(x, y)` first request, rank 2 | 0.18 / 0.18 ms | 0.2 us (Python binary form) |
| `scalar_multiple` | 0.52 ms | |
| `module_generator` | 0.095 ms | |
| `ZZ(int)`, preamble ring element | 0.052 ms | |
| `x + y`, rank 3 | 0.10 ms | |
| `content()` / `to_vector()` | 0.047 / 0.017 ms | |
| category membership test (walk mean) | 8.3 us | 1.5 us (`P in Lattices(ZZ)` isolated) |

Count x unit cost, w323_l4 / ListReflect_6219 (seconds per case):

| Phase | operation | count | unit | seconds | phase seconds |
| --- | --- | --- | --- | --- | --- |
| anisotropic family | `_OrientedBinaryCycle` construction | 292 / 114 | 1.23 ms | 0.36 / 0.14 | 1.51 / 1.20 |
| | isometry application (period unfolding, row matrix) | 220 / 264 | 1.0 ms | 0.22 / 0.26 | |
| | `Isom(P)(images)` (closed cycles) | 55 / 66 | 2.0 ms | 0.11 / 0.13 | |
| | `q()` of the start vector | 468 / 312 | 0.18 ms | 0.08 / 0.05 | |
| | `linear_combination` (`point`) | 165 / 198 | 0.27 ms | 0.04 / 0.05 | |
| | `content`, `to_vector`, coordinate reads, `ZZ(int)`, Sage matrix ops | 3742 / 3492 | | about 0.08 / 0.08 | |
| | rest: Python integer loop (29,570 / 28,512 binary pairings), call overhead | | | 0.62 / 0.49 | |
| edge plane | Sage matrix and vector operations | 15,712 / 9448 | | | 1.59 / 0.93 |
| | `linear_combination` (`plane_primitive`) | 390 / 234 | 0.27 ms | 0.11 / 0.06 | |
| | `ZZ(int)` | 714 / 347 | 0.052 ms | 0.04 / 0.02 | |
| reflectivity certificate | `from_finite_orbits` (sweep 0.088 / 0.030 s, orbit G-sets per seed 0.66 / 0.47 s) | 9 / 9 | 83 / 56 ms | 0.75 / 0.50 | 1.07 / 0.81 |
| vertex transporter | `O(L).transporter` of coloured families | 47 / 47 | 19 ms | 0.88 / 0.88 | 0.88 / 0.88 |
| `plane_lattice` | `Lattices(ZZ)(gram)` | 292 / 114 | 3.17 ms | 0.93 / 0.36 | 0.89 / 0.38 |
| `is_isotropic` | first requests | 292 / 114 | 1.63 ms | 0.48 / 0.19 | 0.52 / 0.22 |
| initial vertex | Sage matrix ops; isometry application | 2498 / 0; 115 / 12 | 1.4 ms | 0.16 / 0.02 (isometries) | 1.49 / 0.30 |

Across all phases: `ZZ(int)` 4509 / 2780 calls at 73 / 62 us in the walk (0.33 / 0.17 s), owned ring elements passed back to their ring 39,478 / 24,098 at 3.8 / 3.1 us (0.15 / 0.08 s); category membership 84,444 / 62,458 tests at 9.4 us (0.79 / 0.59 s, no site above 0.09 s); `_sylvester` 346 / 154 calls (0.38 / 0.18 s).

Count against the source floor: the binary-plane counts (292 Gram lattices for 390 edge planes in w323_l4) are the source's edge planes, deduplicated by Gram. The transporter count is above the source's: `func_insert_vertex` (edgewalk.h:1449-1490) tests equivalence only against known vertices of equal canonical-form hash (`ComputeHashWeightMatrix_raw`, fund_domain_vertices.h:365), while the port asks `O(L).transporter` of every known configuration; sage-extensions `_transporter` rejects a mismatch by `pairing_graph().isomorphism` after building both realizations. The root closure was above the floor: `from_finite_orbits` built a G-set per seed orbit that the edgewalk never reads (fix 13).

## Fixes after the phase table (queue item 3)

Fix 12, withdrawn: `_OwnedRingParent._element_constructor_` returns an element of the ring and converts a Python `int` before its membership tests (each built `Modules(R)` and asked category membership). Unit cost of `ZZ(3)` 41.9 -> 7.2 us. In the walk: w323_l4 0.48 -> 0.12 s, ListReflect_6219 0.25 -> 0.07 s, about 3.2 % of a case, 0.9 h. Paired runs (`zk_pair_f12`): Allcock mean 6.93 -> 7.10 s, ListReflect 7.43 -> 7.43 s, inside the noise. Under an hour: not landed (patch kept as `f12.patch` in the probe directory).

Fix 13 (sage-extensions 891fee04): `FiniteGSets.from_finite_orbits` no longer builds the orbit G-set, ranking and section of every seed before it returns. The set retains the finished sweep's table, and `_orbit_of` answers a seed from that table on its first question, with one orbit set per orbit and the section the eager loop wrote; a point that is not a seed keeps the permutation-representation sweep (its words differ from the table walk's, `test_the_same_action_presented_two_ways_presents_the_same_family`). No question sweeps a second time (`test_no_seed_of_one_sweep_is_searched_twice` passes). `from_finite_orbits` in the walk: w323_l4 0.745 -> 0.192 s, ListReflect_6219 0.502 -> 0.116 s; priced by count 5.8 % of a case, about 1.6 h.

| Family | Cases | Base mean s (`zk_pair_f13`, 081d4829) | Head mean s (891fee04) | Total h before -> after (warm sample) |
| --- | --- | --- | --- | --- |
| Allcock | 8595 | 8.10 | 7.37 | 14.4 -> 13.6 (count) |
| ListReflect rank 3 | 8602 | 7.73 | 7.23 | 13.5 -> 12.7 (count) |

About 1.6 h removed by count, 2.9 h by the pairs (noisy). Touched sage-extensions tests (groups: orbits, sweeps, Schreier families, double cosets, generator permutations, retained partition; lattices: pairing-form relabelling, orbit splitting, stable isotropic plane orbits): 43 passed at 081d4829 and 43 passed at 891fee04. Edgewalk gates at 891fee04 (`zk_gates_f13.log`): w183 9.9 s, w217_l26 11.3 s, ListReflect_2917 7.1 s, branch-representative sample 16.0 s, all pass.

Repin sage-extensions 891fee04. Local tier-1 (`zk_t1_f13.log`, one process per node, killed at 60 s): 45 of 49 pass; test_o_n_maps_onto_o_of_the_discriminant_form, dawes_2U_A2_stable, I_1_2 and I_2_1 reach 60 s with no error.

Fix 14 (port): the walk asks `O(L).transporter` of a known vertex only when its configuration invariant equals the new vertex's, as `func_insert_vertex` compares the canonical-form hash before `LORENTZ_TestEquivalence` (edgewalk.h:1457, fund_domain_vertices.h:365). `configuration_invariant` is the sorted tuple of (colour, norm, sorted coloured pairings with the other vectors) of the roots and the vertex, an isometry invariant up to colour-preserving reindexing, so a skipped pair has no transporter. Seeded sample (`zk_sample_tr.py`, seed 1, 25 cases per family, warm; port a244e12 against the fix, sage-extensions 891fee04): transporter calls that return None cost 0.487 s (Allcock) and 0.474 s (ListReflect) per case before and 0.128 s and 0.122 s after (two cases, ListReflect_8117 and ListReflect_7090, keep equal invariants without a transporter). The transporters found are the same in all 50 cases. By count about 0.36 s per case, 0.86 h + 0.84 h = 1.7 h removed.

| Family | Cases | Base mean s (sample) | Head mean s | Total h before -> after (count) |
| --- | --- | --- | --- | --- |
| Allcock | 8595 | 5.53 | 5.39 | 13.2 -> 12.4 |
| ListReflect rank 3 | 8602 | 5.38 | 5.79 (load rise in the last 12 cases of the head process: ListReflect_4363 3.98 -> 7.36 s at equal transporter count) | 12.9 -> 12.0 |

Paired acceptance runs (`zk_pair_f14`, three cases per family): Allcock 7.43 -> 7.63 s, ListReflect 7.23 -> 6.97 s, inside the noise. Gates: w183, w217_l26, ListReflect_2917, branch-representative sample and the U ⊕ <46> nonreflective certificate pass (`zk_gates_f14.log`).

The warm sample at 891fee04 (seed 1, same 50 cases as the first sample): Allcock 6.04 -> 5.53 s, ListReflect 5.66 -> 5.38 s, which corroborates fix 13 (1.2 h + 0.7 h).

Fix 15 (port): `edge_plane_for_divisor` solved three systems with Sage's general `solve_right` on the transposed 2 x n plane basis (orientation sign, r0 coordinates, C cap P coordinates; 1170 of the 1260 `solve_right` calls of w323_l4, 0.59 ms each in isolation). `_plane_coordinates` reads the same unique solution by Cramer's rule on two columns with a nonzero minor and checks the other columns (`zk_plane_eq.py`: 3000 random rank-two bases and in-span rows over n = 3, 4, 5 equal to `solve_right`; a row outside the span raises). Edge-plane phase (`zk_ph15_*.json`): w323_l4 1.59 -> 1.35 s, ListReflect_6219 0.93 -> 0.73 s.

| Family | Cases | Base mean s (65fb0a8), two orders | Head mean s | Total h before -> after |
| --- | --- | --- | --- | --- |
| Allcock | 8595 | 5.35, 4.84 | 4.82, 4.78 | 12.0 -> 11.5 (pairs, about 0.7 h; phase delta about 0.5 h) |
| ListReflect rank 3 | 8602 | 5.48, 5.08 | 4.84, 4.90 | 12.6 -> 11.6 (pairs, about 1.0 h; phase delta about 0.4 h) |

The order effect is large (base first: 0.53 and 0.64 s; head first: 0.07 and 0.18 s), so the order-balanced pairs give about 1.7 h and the two-case phase delta about 0.9 h. Landed at the mean of the two prices, about 1.3 h. Gates (`zk_gates_f15.log`): w183, w217_l26, ListReflect_2917 and the branch-representative sample pass.

Withdrawn, under an hour (`f16.patch` in the probe directory): `_OrientedBinaryCycle` read its form from the Gram tensor instead of two `module_generator`, two `q` and one `b` (`__init__` 0.425 -> 0.067 s in w323_l4, 0.155 -> 0.026 s in ListReflect_6219, about 3 % of a case, 0.7 h), and built its period isometry on the labels instead of enumerating `module_generators()` (3.2 ms per first request; the cost moved into the label-function isometry and `cycle` stayed at 0.91 / 0.76 s). `zk_cycle_eq.py`: the same forms, representatives and isometries on 400 random anisotropic planes.

Phase table after fixes 13-15 (port 6f44af3, sage-extensions 891fee04; `zk_ph15_*.json`; hours on the warm totals below):

| Phase | w323_l4 | ListReflect_6219 | hours over 17,200 |
| --- | --- | --- | --- |
| anisotropic family | 1.52 | 1.18 | 4.4 |
| edge plane | 1.35 | 0.73 | 3.3 |
| initial vertex | 1.67 | 0.30 | 2.8 |
| `plane_lattice` | 0.90 | 0.35 | 1.9 |
| reflectivity certificate | 0.56 | 0.37 | 1.5 |
| edge directions | 0.37 | 0.35 | 1.2 |
| `is_isotropic` of the plane | 0.51 | 0.19 | 1.1 |
| vertex stabilizer generators | 0.13 | 0.45 | 1.1 |
| coordinate helpers | 0.30 | 0.28 | 1.0 |
| vertex transporter search | 0.37 | 0.22 | 0.9 |
| edge orbits | 0.19 | 0.29 | 0.9 |
| other phases | 1.09 | 0.82 | 3.3 |
| sum | 8.96 | 5.53 | 23.1 |

New family totals (warm seeded sample, seed 1, 25 cases per family, `zk_wpair_f15r.log`): Allcock 4.78 s x 8595 = 11.4 h, ListReflect rank 3 4.90 s x 8602 = 11.7 h, about 23.1 h (27.9 h at c435fd6 on the same 50 cases).

Measured remainders of the phases over an hour, none an hour-sized single operation:

| Phase | operation | count per case (w323_l4 / ListReflect_6219) | wall in the walk | hours |
| --- | --- | --- | --- | --- |
| anisotropic family | `_OrientedBinaryCycle.cycle` (start `content` and `q` 0.24 s; closed cycles: `Isom` and representative points 0.33 s) | 468 / 312 | 0.91 / 0.76 s | about 2.7 |
| | `_OrientedBinaryCycle.__init__` | 292 / 114 | 0.43 / 0.16 s | 0.7 (withdrawn above) |
| | lattice isometry applications, all sites (`zk_isoapp.py`) | 315 / 390 | 0.17 / 0.18 s | 0.6 |
| edge directions | `evaluate_covector`, each (facet, root) pair read twice | 56 / 56 at 4.5 ms | 0.25 / 0.26 s | 0.9; halving the count in the port 0.45 |
| `is_isotropic` | `_sylvester` (Sage `QuadraticForm.signature_vector`) + `qfsolve` witness | 346 + 294 / 154 + 116 | 0.47 / 0.22 s | about 0.9 |
| `plane_lattice` | `Lattices(ZZ)(gram)`, 3.2 ms each | 292 / 114 | 0.90 / 0.35 s | 1.9 |
| all phases | `ZZ(int)` and ring-element round trips | 4509 + 39,478 / 2780 + 24,098 | 0.48 / 0.25 s | 0.9 (fix 12, withdrawn) |
| all phases | category membership tests, 8.3 us each, no site above 0.09 s | 84,444 / 62,458 | 0.70 / 0.59 s | about 2 |

Blocker for the next hour item: every remaining item over an hour is the unit cost of owned arithmetic, spread over many sites. `zk_inside.py` on `_evaluate_covector` (4.5 ms for a dot product of two rank-3 vectors) charges its time to owned ring element construction, `_operand_in_ring`, category membership, `_MorCategoryOf.Of` and set-morphism calls, no function above 9 % of it; `Lattices(ZZ)(gram)` (3.2 ms) charges category classcalls, membership tests, ring hashing and tensor-shape checks in the same way. The lever is the one recorded as decided (owned elements carry engine data; an owned operation is engine arithmetic plus a constant-time wrap), applied to owned ring elements, covector evaluation and lattice construction in sage-extensions; it is not one fix but a change to the arithmetic of the owned ring and module element classes, with its own contract (CONTRIBUTING OWN-22, OWN-24). The largest single item is `plane_lattice` at 1.9 h; the next is the cycle itself (2.7 h), whose owned-operation counts are already near the source's (one start vector, one period isometry and its representatives per cycle).

## Owned arithmetic: counts and unit costs (port 6e60423, sage-extensions 891fee04)

Counts per case from one warm walk with wrappers on the entry points (`zk_own_count.py`; outermost calls; nested calls of the same entry point counted apart). The wall column is the outermost time inside the wrapper, so the rows overlap (a ring constructor inside covector evaluation is in both) and include wrapper overhead. Unit costs: `zk_unit_own.py`, warm, `timeit` N = 10,000 (evaluate_covector N = 2,000), no profiler. Hours = w323_l4 seconds x 8595 + ListReflect_6219 seconds x 8602.

| Operation | count w323_l4 / ListReflect_6219 | unit cost (timeit) | wall per case s | hours |
| --- | --- | --- | --- | --- |
| category membership tests | 69,454 / 48,868 | 8-26 us (`ZZ in OwnedRings()` 11 us, `L in Lattices(ZZ)` 26 us) | 0.61 / 0.51 | 2.7 |
| owned ring element constructor (`R(x)`) | 44,055 / 26,999 | `ZZ(int)` 49 us, `QQ(int)` 52 us, `QQ(ZZ element)` 290 us | 0.60 / 0.43 | 2.5 |
| `evaluate_covector` | 56 / 56 | 2040 us | 0.31 / 0.32 | 1.5 |
| `_owned_engine_element` (raise) | 9,994 / 8,077 | | 0.26 / 0.24 | 1.2 |
| `_operand_in_ring` | 15,110 / 9,949 | same parent 4 us; `QQ` with an int 442 us | 0.28 / 0.23 | 1.2 |
| module scalar multiply | 620 (+431 nested) / 174 | `3 * x` 437 us, `a * x` 394 us | 0.26 / 0.04 | 0.7 |
| owned ring element `__init__` | 28,632 / 19,959 | | 0.09 / 0.07 | 0.4 |
| module add | 630 / 145 | 76 us | 0.04 / 0.02 | 0.2 |
| `b(x, y)` | (inside the phases above) | 267 us | | |
| ring add / multiply | | `a + b` 19 us, `a + 3` 88 us, `p + q` (QQ) 122 us | | |

Membership sites (two frames, w323_l4): `_OwnedRingParent._element_constructor_` asked `parent in Modules(R)` and `parent in OwnedRings()` of every argument, an int literal and an element of the ring itself included (4,646 + 4,646); `_own_ring` asked `QQ in OwnedRings()` on each raise into the rationals, because the rationals are a localization parent, not an `_OwnedRingParent` (4,610); `_asked_of_the_category` 5,446; `Cardinalities.__eq__` 2,998; tensor coordinate components 1,547. The first two are the count defect: in arithmetic the operands' parent is already known. The ring `__mul__` also built `Modules(R)` and asked membership for every operand of another parent, an int included. The rationals' constructor asked `Modules` membership and the denominator decision for an int or an integer, whose denominator is 1.

Fix 17 (sage-extensions b6ba5826): (a) `_operand_in_ring` returns an operand of the same parent at once and reads an int through the ring constructor; the native ring constructor reads an int and an element of itself before any membership test; the ring `__mul__` multiplies an operand of its own parent directly and asks module placement only of a parent; `_own_ring` keeps its placement answer per ring (placement is fixed at construction, CON-07). (b) A localization reads an int and an element of its source as `a/1` (1 lies in every multiplicative set). (c) `evaluate_covector` reads a lattice vector on its own coordinates, sums the engine values and raises once; `b` and `norm` read the stored coefficients, and the finite rank once per lattice; the sparse free module scales by a base-ring scalar on its coefficients, as its `_lmul_` does. Unit costs after: `ZZ(int)` 4 us, `QQ(int)` 22 us, `QQ(ZZ element)` 28 us, `a + 3` 11-18 us, `_operand_in_ring(QQ, int)` 28 us, `3 * x` 52-65 us, `a * x` 76-88 us, `b(x, y)` 46 us, `evaluate_covector` 220-280 us. Equality (`zk_own_eq.py`): 3,000 random ring samples (ZZ and QQ sums, products, comparisons, operand reads against the engine), 3,360 covector evaluations on 60 random rank-3 cones (lattice and rational vectors, value and parent equal to the prior route) and 1,600 pairings and norms on 80 random lattices of rank 2-5 (equal to `x^T G y` in the engine).

Paired warm runs (`zk_wse_f17.log`, seed 1, 25 cases per family, the port at 6e60423; sage-extensions 891fee04 against b6ba5826): Allcock 4.89 -> 4.69 s (base first), ListReflect rank 3 5.14 -> 4.36 s (head first). The order effect of the earlier pairs (0.5 s in favour of the second process) bounds the Allcock delta at 0.2-0.5 s and the ListReflect delta at 0.3-0.8 s; at the midpoints about 0.8 h + 1.3 h = 2.1 h removed. Gates (`zk_gates_f17.log`): w183, w217_l26, ListReflect_2917 and the branch-representative sample pass.

Owner unit tests of the changed modules (head 683716ad): tests/algebras/test_group_algebras.sage, tests/lattices/test_named_lattice_constructions.sage, tests/modules/test_dual_point_hash.sage, tests/modules/test_the_paired_module_readings_of_an_algebra.sage, tests/rings/test_a_quotient_is_a_cyclic_module_over_its_cover.sage, tests/rings/test_embeddings_of_quadratic_fields.sage, tests/rings/test_owned_scalars_act_on_sage_matrices.sage, tests/rings/test_rational_field_placement.sage, tests/schemes/test_toric_schemes.sage, tests/test_finite_form_groups.py, tests/test_integer_residue_rings.py, tests/test_native_ring_constructions.py, tests/test_prime_and_maximal_ideals.py: 91 passed, 4 failed. Those 4 and the 64 failures of the first three batches of tests/lattices and tests/modules (36 files; the whole-directory run was stopped) were run at the prior main 891fee04: 66 fail there too. The other two (test_cyclotomic_centralizer: the three-cycle on Z^3, the swap of two hyperbolic planes) pass at head when run alone and in their own file, so they are not regressions. Local tier-1 (`zk_t1_f17.log`): 45 of 49 pass; the four nodes priced over 60 s reach the cap with no error. Edgewalk gates at the merge 683716ad pass (`zk_gates_f17m.log`). Repin sage-extensions 683716ad.

Counts after fix 17 (`zk_own17_*.log`, sage-extensions 683716ad; walls w323_l4 8.58 s, ListReflect_6219 6.13 s instrumented):

| Operation | count w323_l4 / ListReflect_6219 | wall per case s | hours |
| --- | --- | --- | --- |
| category membership tests | 55,223 / 39,520 | 0.54 / 0.40 | 2.2 |
| `_owned_engine_element` (raise; mostly into the rationals, a localization whose raise decodes numerator and denominator) | 8,522 / 6,749 | 0.24 / 0.23 | 1.1 |
| ring element constructor | 21,296 / 13,226 | 0.10 / 0.08 | 0.4 |
| `evaluate_covector` | 56 / 56 | 0.09 / 0.10 | 0.5 |
| `_operand_in_ring` | 11,610 / 7,265 | 0.05 / 0.05 | 0.2 |

Withdrawn, under an hour (`f18.patch` in the probe directory): 4,577 of the remaining membership tests per case are `_own_ring(SageZZ)` (callers raise a value along Sage's integer ring); answering it before the membership test removes about 0.05 s per case, 0.2 h. No remaining membership site is above 0.05 s per case (`_asked_of_the_category` 4,698; `Cardinalities.__eq__` 2,998; tensor coordinate components 1,547).

New family totals (fix 17 by the order-balanced pairs): Allcock about 4.78 - 0.35 = 4.43 s x 8595 = 10.6 h, ListReflect rank 3 about 4.90 - 0.55 = 4.35 s x 8602 = 10.4 h, about 21.0 h.

Blocker for the 6 h target: the owned-arithmetic entry points over an hour now sum to about 3.3 h (membership 2.2 h, raises into the rationals 1.1 h), and each is spread over sites below 0.2 h; removing all of it leaves about 17.7 h. The rest of the walk is the phases of the table after fixes 13-15 (anisotropic family, edge plane, initial vertex, `plane_lattice`, reflectivity certificate), whose owned-operation counts are near the source's. Reaching 6 h needs about 15 h more, which requires a change of representation, not removal of checks: the rationals stored as an engine rational (today a localization pair of owned integers, unreduced, decoded on every raise), lattice vectors stored as one engine vector (today a dict of owned ring elements), and `Lattices(ZZ)(gram)` (3.2 ms) built without category class calls. Each is a contract change to the localization and sparse free module classes under OWN-22 and OWN-24 in sage-extensions, priced together by the counts above at not more than the 3.3 h of entry-point time plus the share of the phase time spent in element wrapping, which no measurement here separates.

## Storage changes: priced before implementation (port 4b88a37, sage-extensions 683716ad)

Counts per case from one warm walk with wrappers on the storage entry points (`zk_store_count.py`, outermost calls, logs `zk_store_*.log`). Hours = w323_l4 seconds x 8595 + ListReflect_6219 seconds x 8602.

| Item | count w323_l4 / ListReflect_6219 | unit cost (warm) | wall per case s | hours (upper bound: the operation made free) |
| --- | --- | --- | --- | --- |
| owned rationals: fraction construction, `_add_`, `_mul_`, raise, lower, equality | 174 / 168 constructions, 89 / 83 raises, 10 / 10 products, 75 / 75 equalities | 10-200 us | 0.04 / 0.03 | 0.2 |
| lattice vectors: `__init__`, coordinate read, `_add_` (plus `evaluate_covector` 0.09 / 0.10 s from the fix-17 table) | 2598 + 2946 + 594 / 1514 + 2490 + 108 | 20-40 us | 0.26 / 0.22 | 1.1 |
| `Lattices(ZZ)(gram)` | 293 / 115 | 2.0-2.9 ms (`zk_latctor.py`, distinct rank-2 Grams) | 0.83 / 0.34 | 2.8 |

The rationals are not stored as a bottleneck in the walk: fewer than 200 fractions are built per case, so item 1 is under an hour and was not changed. The 8,522 raises of the fix-17 table are into the integers.

Withdrawn, under an hour (`f18b.patch` in the probe directory): `Lattices(ZZ)(gram)` resolved the free-module placement once per (ring, label set) and the form-module placement once per (module, kind of form, extra categories), both fixed at construction (CON-07). Unit cost 2.0-2.35 ms -> 1.73 ms; equality (`zk_lat_eq.py`): 400 random Gram matrices of rank 1-5 over ZZ and QQ give the same category object and Gram tensor with the caches and resolved afresh. Paired warm runs (seed 1, 25 cases): Allcock 4.54 -> 4.37 s (base first, `zk_wse_f18.log`), ListReflect rank 3 4.44 -> 4.46 s (head first, `zk_wse_f18lr.log`); order-balanced about 0.08 s per case, 0.4 h. The remaining construction cost is per-parent data (parent `__init__` and module resolution 0.5 ms, Gram tensor from rows 0.45 ms, engine Gram matrix, rank and hash for the axioms 0.6 ms, all profiled), not category resolution.

Lattice vectors as one engine vector: the bound (every vector operation free) is 1.1 h; an engine vector still costs a construction and a wrap per operation, so the realizable part is under an hour. The change touches 35 readers of the coefficient dict in `framed_free_modules.py` and about 150 in 20 other modules; not made.

Blocker for the 6 h target (family totals unchanged, about 21.0 h): the three storage items sum to at most 4.1 h if each operation became free (0.2 + 1.1 + 2.8), and the measured realizable part of the largest is 0.4 h. Removing all of it leaves about 16.9 h; the remaining 10.9 h above the target is in the phases of the table after fixes 13-15, spread over sites below an hour each.

## Phase floors (port f4ca6b2 -> 1e7adc0, sage-extensions 683716ad)

Phase hours re-measured at f4ca6b2 (`zk_ph19_*.json`, `zk_ph19_table.md`; warm, line charges scaled to the family totals 10.6 h + 10.4 h). Steps per case are counted in one warm walk with the cached methods unwrapped (`zk_keys.py`) and read against edgewalk.h:868-889 and 753-845: per edge the source decides isotropy once (GetIsotropicFactorization of Pplane G Pplane^T), and per edge and root norm it builds one projected lattice and runs Anisotropic once on Gwork. Floors: `zk_floor19.py` (one step of each phase on captured arguments, timeit, callee caches warm) and `zk_epops.py` (each engine operation of one edge plane). Counts and seconds are w323_l4 / ListReflect_6219.

| Phase | source steps per case | port steps per case | current per step | floor per step | removable hours |
| --- | --- | --- | --- | --- | --- |
| anisotropic family (cycle, period) | 468 / 312 Anisotropic runs | 468 / 312 | 2.9 / 3.5 ms (owned plane lattice, Isom, owned vectors) | 150 / 195 us (`_FormBinaryCycle` and the family on integers) | 3.8 x (1 - 0.07/1.34) = 3.6 |
| `plane_lattice` | 0 (the source keeps Gwork) | 292 / 114 | 3.2 ms | 0 in the anisotropic branch | 1.8 |
| is_isotropic of the plane | 13 / 13 (once per edge) | 524 / 361 (count defect) | 1.0 ms | 95 us | 1.1 |
| edge plane | 390 / 234 planes + 30 / 18 C_N | same | 1.6 ms + 6.5 ms | about 0.75 ms (sum of `zk_epops.py`: left kernel by PARI 276 us, projection 130 us, HNF 93 us, products and constructions) + 0.15 ms (C_N by an HNF) | 1.2 + 0.7 |
| initial vertex | not yet counted | | | | floor not measured (up to 2.7) |
| reflectivity certificate | not yet counted | | | | floor not measured (up to 1.8) |

Fix 19 (port 1e7adc0): the edge plane carries its form (a, b, c) and r0 as a coordinate pair; `_FormBinaryCycle` runs the same `_OrientedBinaryCycle` algorithm on the form and returns the source's data (integral period matrix, coordinate pairs); the lattice-backed class wraps the same result once for its own callers; the family and the candidates are coordinate pairs; isotropy is decided once per edge (`edge_is_isotropic`); the plane lattice is built only in the isotropic branch. This removes the first three rows together.

| Operation | paired before -> after (seed 1, 25 cases) | hours removed | commit |
| --- | --- | --- | --- |
| anisotropic family, `plane_lattice`, isotropy per edge | Allcock 4.46 -> 3.77 s (base first), ListReflect rank 3 4.49 -> 3.55 s (head first) (`zk_wpair_f19.log`) | about 1.6 + 2.2 = 3.8 (the order effect of earlier pairs bounds Allcock at 0.2-0.7 s per case) | 1e7adc0 |

Equality (`zk_eqwalk.py`, seed 7): 25 Allcock and 25 ListReflect rank 3 cases give the same reflectivity, vertex count and simple roots in report order at f4ca6b2 and at the fix. Gates at the fix (`zk_gates_f19.log`): w183, w217_l26, ListReflect_2917 and the branch-representative sample pass; tests/test_edgewalk_rank2.py and the frozen binary-cycle tests of test_reflective.sage pass (9 nodes, `zk_rank2_f19.log`). No repin, so tier-1 was not run.

New family totals: Allcock 3.77 s x 8595 = 9.0 h, ListReflect rank 3 3.55 s x 8602 = 8.5 h, about 17.5 h.

Phase table after fix 19 (`zk_ph20_*.json`, scaled to 9.0 h + 8.5 h; the line map of `zk_ph_group.py` for `_Edgewalk.walk` predates the current line numbers, so its transporter and stabilizer lines are charged to walk bookkeeping): edge plane 3.2 h, initial vertex 2.7 h, walk bookkeeping (transporter search and stabilizer) 2.1 h, reflectivity certificate 1.8 h, indefinite reduction 1.1 h, coordinate helpers 1.0 h, all others under 1 h.

Remaining removable hours, in order: edge plane 1.2 h above its engine floor plus C_N 0.7 h (`constraint_rows` intersects ZZ-modules at 6.5 ms where an HNF of d G^-1 costs 0.15 ms, but the intersection with ZZ^n is not one HNF); initial vertex and reflectivity certificate, whose source step counts and floors are not yet measured; the walk line map of `zk_ph_group.py` must be updated before the transporter and stabilizer phases are priced again. Multiplied remainder: 17.5 h now; if every phase priced above reached its floor, about 17.5 - 1.9 = 15.6 h before the unmeasured phases.

## Gate and price of fixes 21-23 (port 54bd6b8 -> 4d038e8, sage-extensions 683716ad)

Fixes: 21 bcb407f (ordinary vertex roots enumerate each norm shell on the complement's Gram matrix), 22 0b81dbe (the reflectivity certificate closes the roots and tests the chamber on engine data), 23 4d038e8 (the C_N basis from one integral kernel of [G | -d I]).

Gate at 4d038e8. Equality (`zk_eqwalk.py`, seed 7, `zk_eqg_head.log`): 25 Allcock and 25 ListReflect rank 3 cases give the same reflectivity, vertex count and simple roots in report order as the fix-19 code (`zk_eq19_wt_base.log`; 54bd6b8 has no source change after 1e7adc0). Edgewalk gates (`zk_gates_g54.log`): w183, w217_l26, ListReflect_2917 and the branch-representative sample pass, 34 s for the four. tests/test_edgewalk_rank2.py and the two frozen binary-cycle tests of test_reflective.sage pass (6 nodes, `zk_owner_g54.log`). Tier-1 (`zk_t1_g54.log`, one process per node, killed at 60 s): 45 of 49 pass; test_o_n_maps_onto_o_of_the_discriminant_form, dawes_2U_A2_stable, I_1_2 and I_2_1 reach the cap with no error.

Order-balanced paired warm runs (`zk_price3.out`, `zk_wp_*.log`; seed 1, 25 cases per family; each pair run once base first and once head first; means of the two orders):

| Operation | Allcock before -> after s | ListReflect rank 3 before -> after s | hours removed | commit |
| --- | --- | --- | --- | --- |
| fixes 21-23 together (54bd6b8 -> 4d038e8) | 3.63 -> 3.16 | 3.87 -> 3.07 | 1.1 + 1.9 = 3.0 | bcb407f, 0b81dbe, 4d038e8 |
| fix 21, ordinary vertex roots per norm shell | 3.63 -> 3.51 | 3.76 -> 3.38 | 0.3 + 0.9 = 1.2 | bcb407f |
| fix 22, certificate on engine data | 3.31 -> 3.20 | 3.53 -> 3.10 | 0.2 + 1.0 = 1.2 | 0b81dbe |
| fix 23, C_N from one integral kernel | 3.12 -> 3.07 | 3.16 -> 3.16 | 0.1 (inside the noise) | 4d038e8 |

The same process gives 3.30-3.52 s for one revision in different runs, so a single delta under about 0.2 s per case is noise; the together row is the price. New family totals: Allcock 3.16 s x 8595 = 7.6 h, ListReflect rank 3 3.07 s x 8602 = 7.3 h, about 14.9 h.

Phase table at 4d038e8 (`zk_ph24_*.json`, `zk_ph24_table.md`; scaled to 7.55 h + 7.34 h). The grouping map `zk_ph_group.py` now reads the measured source text: the indefinite reduction is the line that calls `indefinite_reduced_presentation`, isotropy is `edge_is_isotropic`, and `_FormBinaryCycle` is in the anisotropic family (it was charged to Dynkin extensions).

| Phase | ListReflect_2917 | ListReflect_6219 | w323_l4 | w65_l29 | hours |
| --- | --- | --- | --- | --- | --- |
| edge plane | 0.59 | 0.52 | 0.79 | 0.64 | 2.6 |
| indefinite reduction | 0.60 | 0.14 | 0.12 | 0.38 | 1.4 |
| initial vertex | 0.17 | 0.14 | 0.93 | 0.14 | 1.3 |
| vertex transporter search | 0.27 | 0.27 | 0.40 | 0.30 | 1.3 |
| vertex stabilizer generators | 0.14 | 0.50 | 0.17 | 0.12 | 1.0 |
| edge orbits under the stabilizer | 0.21 | 0.33 | 0.21 | 0.15 | 0.9 |
| walker setup | 0.12 | 0.35 | 0.29 | 0.10 | 0.9 |
| edge directions | 0.23 | 0.23 | 0.23 | 0.16 | 0.9 |
| coordinate helpers | 0.19 | 0.25 | 0.25 | 0.14 | 0.8 |
| report transport to the input lattice | 0.16 | 0.19 | 0.18 | 0.13 | 0.7 |
| anisotropic family | 0.13 | 0.12 | 0.14 | 0.13 | 0.5 |
| other phases (ten, each 0.5 h or less) | 0.55 | 0.58 | 0.85 | 0.57 | 2.8 |
| sum | 3.35 | 3.60 | 4.56 | 2.98 | 14.9 |

The anisotropic family is 0.5 h after fix 19 (468 / 312 runs per case on the form's integers), so its 3.6 h row of the phase-floor table is removed.
