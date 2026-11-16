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
| 1 | `ci/acceptance_tier1.txt` | 50 | For each dispatch branch of each pinned algorithm, the smallest frozen node (least rank, determinant or \|O(q)\| that reaches the branch) whose expected answer comes from an independent source or from a construction, not from a round trip of the port |
| 2 | `ci/acceptance_tier2.txt` | 86 | The published-result reproductions that the #64 table requires |
| 3 | `full` | 26907 | Every acceptance node; tier 3 is the full corpus less tiers 1 and 2 |

Tier 1 is the correctness evidence (priority 1 above). Recorded tier-1 time: 879 s over 30 priced nodes; 20 nodes have no recorded time. Tier 2 has 7 priced nodes (6167 s) and 79 unpriced nodes.

A price is the least recorded passed call time of the node: a CI Acceptance BUDGET line (runs 43505668003 and 43494572338), a local pytest report under `.tmp/`, or a row of this file or of `.tmp/acceptance-coverage.md` with its revision. A node without such a record is unpriced. The derived budgets in `ci/test_budgets.json` are not prices. Tier-2 vector-orbit nodes are the remaining cases of `test_vector_orbits.sage` that #64 lists.

| Node | Tier | Branch covered | Invariants | Recorded time |
| --- | --- | --- | --- | --- |
| `test_reflective.sage::test_edgewalk_branch_representative_sample` | 1 | Allcock edgewalk: binary walk, oriented binary cycle (go-left, shorter pair), cusp start, cocompact start and increasing-prefix corner search, ordinary, cusp and extension roots, isotropic next vertex, reflective stop, nonreflective stop | n=2..4; U+<2>, U+<46>, compact hexagon, allcock_rk3_4, allcock_rk3_2, ListSimpleRootSystem_0 | 59.8 s (.tmp/sample61/report_2.jsonl) |
| `test_perfect_domains.sage::test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_3]` | 1 | Perfect domains, total mode: initial cell, one cell orbit, O(L) generation | n=3 sig (2,1) |det|=1, 1 cell | unpriced |
| `test_perfect_domains.sage::test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_1]` | 1 | Perfect domains, total mode: flips, cell transporter, side pairings | n=3 sig (2,1) |det|=8, 2 cells | unpriced |
| `test_perfect_domains.sage::test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_3]` | 1 | Perfect domains, isotropic mode | n=3 |det|=1, 1 cell | unpriced |
| `test_perfect_domains.sage::test_non_diagonal_positive_shell_uses_the_pinned_cvp_center` | 1 | Perfect domains: positive-vector shell, CVP center | n=3, against exact enumeration | unpriced |
| `test_perfect_domains.sage::test_non_diagonal_isotropic_shell_uses_the_exact_sphere_center` | 1 | Perfect domains: isotropic shell, exact-sphere center | n=3, against exact enumeration | unpriced |
| `test_perfect_domains.sage::test_large_pairing_graph_fallback_returns_exact_transporter` | 1 | Perfect domains: transporter above 256 rows (no graph canonization) | n=2, 257 rows | unpriced |
| `test_perfect_domains.sage::test_large_cell_cone_uses_exact_normaliz_h_presentation` | 1 | Perfect domains: large-configuration cone (Normaliz) | 24000 rows | unpriced |
| `test_isometry.sage::test_isometry_verdicts_match_certified_pairs[hecke_U_vs_U2]` | 1 | Isometry: binary, prefilter rejection | n=2 sig (1,1) |det|=1 | unpriced |
| `test_isometry.sage::test_isometry_verdicts_match_certified_pairs[hecke_3x3_swapped_basis]` | 1 | Isometry: Lorentzian witness (perfect-cell transporter) | n=3 sig (1,2) |det|=1250 | unpriced |
| `test_isometry.sage::test_conway_sloane_same_genus_pair_is_not_isometric` | 1 | Isometry: same genus, spinor-genus negative verdict | n=3 sig (2,1) |det|=128 | unpriced |
| `test_isometry.sage::test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U]` | 1 | Isometry: Eichler covering search witness | n=4 sig (2,2) |det|=4, U+U(2) | unpriced |
| `test_isometry.sage::test_degenerate_isometry_witness_respects_complement_coordinates` | 1 | Isometry: degenerate, radical extension witness | n=3, A2 plus rank-1 radical | unpriced |
| `test_orthogonal_group.sage::test_degenerate_orthogonal_group_contains_radical_shear_and_gl_action` | 1 | O(L): degenerate, radical extension | n=2, <2> plus radical | 0 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_zero_form_orthogonal_group_is_the_full_general_linear_group` | 1 | O(L): zero form, GL_2(Z) | n=2, q=0 | 0 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_degenerate_orthogonal_group_embeds_non_diagonal_complement` | 1 | O(L): definite kernel inside the radical extension | n=3, A2 plus radical | 0.4 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_orthogonal_group_of_the_hyperbolic_plane_has_order_four` | 1 | O(L): binary, isotropic | U, n=2 |det|=1, |O|=4 | 0 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_ambiguous_binary_form_has_an_improper_automorphism[hecke_binary_3_1_-3]` | 1 | O(L): binary, anisotropic reduction cycle | n=2 D=37 | unpriced |
| `test_orthogonal_group.sage::test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q` | 1 | O(L): Lorentzian route; image of O+ in O(q) | n=3 sig (1,2) |det|=3 | 0.8 s (CI Acceptance BUDGET line) |
| `test_orthogonal_group.sage::test_mertens_omega_and_minus_one_generate_o_l` | 1 | O(L): Lorentzian route, generators against Mertens | n=3 sig (2,1) |det|=155 | 4.9 s (CI Acceptance BUDGET line) |
| `test_isotropic_orbits.sage::test_o_n_maps_onto_o_of_the_discriminant_form` | 1 | O(L): Eichler route; discriminant image | N=U+U(2)+E8(-2), n=12 sig (2,10) |A|=2^10 | 468 s (.tmp/acceptance-coverage.md, provider 42d7361c) |
| `test_vector_stabilizers.sage::test_isotropic_pointed_isometry_rejects_different_divisibility` | 1 | Pointed isometry: invariant rejection | n=4 | 0 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_split_divisibility_one_pointed_isometry_carries_actual_negative_marking` | 1 | Pointed isometry: isotropic, div 1 | U, n=2 | 0.3 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_divisibility_two_isotropic_pointed_isometry_has_exact_witness` | 1 | Pointed isometry: isotropic, div 2 | U(2)+U(2), n=4 | 0.8 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_odd_lattice_isotropic_pointed_isometry_uses_general_divisibility_two_route` | 1 | Pointed isometry: isotropic, odd lattice, general div 2 route | odd binary, n=2 | unpriced |
| `test_vector_stabilizers.sage::test_general_isotropic_pointed_isometry_carries_negative_marking` | 1 | Pointed isometry: isotropic, general route | U(2)+U(2), n=4 | 0.2 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_divisibility_four_isotropic_stabilizer_constructs_exact_containing_bound` | 1 | Isotropic line stabilizer: finite correction, cokernel 4 | U(4)+<-8>, n=3 | 1 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_isotropic_pointed_automorphism_generators_preserve_marking_and_form` | 1 | Isotropic line and pointed stabilizer generators | U+U, n=4 | 1.7 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_divisibility_four_isotropic_plane_stabilizer_constructs_exact_containing_bound` | 1 | Isotropic plane stabilizer: finite correction, cokernel 16 | U(4)+U(4), n=4 | 0.9 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_isotropic_plane_and_flag_automorphism_generators_preserve_their_marks` | 1 | Isotropic plane and flag stabilizer generators | images 48 and 12 in O over F_3 | 2.7 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_isotropic_plane_isometry_carries_changed_basis_plane` | 1 | Isotropic plane isometry | U(4)+U(4), n=4 | unpriced |
| `test_vector_stabilizers.sage::test_isotropic_flag_isometry_preserves_line_inside_plane` | 1 | Isotropic flag isometry | U(4)+U(4), n=4 | unpriced |
| `test_vector_stabilizers.sage::test_vector_stabilizer_has_oscars_order` | 1 | Anisotropic pointed stabilizer: v-perp with gluing | n=3 sig (1,2), v-perp binary, order 2 | 0.5 s (CI Acceptance BUDGET line) |
| `test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_5]` | 1 | Anisotropic pointed isometry: complement isometry and integral correction | E10, n=10 |det|=1 | 2.2 s (CI Acceptance BUDGET line) |
| `test_isotropic_orbits.sage::test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_3261]` | 1 | Isotropic classes: definite, empty | n=3 sig (3,0) | unpriced |
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
| `test_vector_orbits.sage::test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[2]` | 1 | Vector classes: Lorentzian cone (marked cell orbits) | E10, s=2, 1 orbit | unpriced |
| `test_vector_orbits.sage::test_i_2_10_is_transitive_on_norm_minus_one_vectors` | 1 | Vector classes: Eichler | I_{2,10}, n=12, s=-1 | 7.8 s (.tmp/sample61/report_0.jsonl) |
| `test_centralizers.sage::test_centralizer_images_match_oscar[oscar_lattices_with_isometry_133]` | 1 | Centralizer image: Hermitian | n=3 sig (2,1) |A|=12, ord(f)=6 | unpriced |
| `test_centralizers.sage::test_centralizer_images_match_oscar[oscar_lattices_with_isometry_118]` | 1 | Centralizer image: gluing | n=4 sig (2,2) |A|=5, ord(f)=5 | unpriced |
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
