"""A-priori time budgets for scope C.
Theory: references/complexity/perfect-domains-and-lorentzian-cells.md (section cited in each basis string).
Operation counts are derived by inspection; unit costs are measured primitives read from
ci/test_budget_models/unit_costs.json (KeyError if a key is absent).

Scope: tests/test_pinned_lorentzian_cells.py, test_perfect_cell_transporter.py,
test_perfect_cell_simultaneous_forms.py, test_large_cell_cone_route.py,
test_target_cell_symmetries.py.

Each budget = branch work from mathematical invariants * unit costs (no multiplicative factor).
Branches cited in src/sage_indefinite_port/indefinite/lorentzian_cells.py (LC):
  cell_transporter LC:1975-2160
    * identity shortcut       LC:2020  (same lattice, same rows)
    * small_graph             LC:2022-2032  n <= _PAIRING_GRAPH_CANONICAL_LABEL_LIMIT (256, LC:910)
                              and no cached stabilizer -> Isom(...).transporter on the coloured
                              pairing graph (AutomorphismPolytope.g:288-339)
    * basis_search            otherwise: invariants (Smith, moment det, row sum) LC:2050-2130,
                              colour refinement if n <= 4096 (LC:2142), then backtracking over
                              target images of a source basis; leaves ~ |Stab(cell)| when
                              all_images, ~1 witness otherwise.
  cell_cone LC:1061: rows < 10_000 -> ray route; >= 10_000 -> Normaliz H route.
  _facet_mask_data LC:1417 / flip_across LC:1874 / _perfect_form_hash_key LC:2876-2910:
    rows > 4096 -> linear-cost large route (no quadratic pairing graph).
Upstream: Lorentzian.g:660-696 (transporter), polyhedral_common fund_domain_vertices.h:146.

Unit costs = operation count x measured primitive (json key in brackets):
  LAT    [lattice_construction_rank_3]           owned lattice construction
  OWN    [lattice_element_pairing]               one owned element/covector evaluation
  DOT    [zz_dot_12]                             one compiled integer dot product / row pass step
  GAPOP  [gap_order_matrix_group_so_plus_10_2]   one GAP group order (upper bound for degree <= 12, order <= 48)
  NMZ    [normaliz_cone_6_rays] + 10 OWN         one Normaliz cone + conversion of <= 10 covectors
  LEAF   2 [zz_matrix_product_4] + 12 DOT        congruence A G A^T + row-set check, N <= 12
  INV    3 [qq_rank_4] + [zz_matrix_product_4]   Smith form, moment adjugate, determinant
  ISOM   256 DOT + 4 OWN + LEAF                  pairing graph on N <= 16, isometry assembly, certificate
  SHELL  300 [pari_qfminim_e8_point] + [pari_lllgram_10]   Fincke-Pohst tree <= 300 nodes after LLL
  FLIP   [qq_rank_4] + SHELL                     facet kernel + one shell
  CALL   [pari_lllgram_10]                       one small exact routine (dim <= 5), bounded by LLL at dim 10
  PERM   [gap_permutation_image_527]             one permutation image of a point
Library start-up is once per session and not charged.
"""
from __future__ import annotations

import json
from pathlib import Path
from statistics import median

_U = json.loads((Path(__file__).with_name("unit_costs.json")).read_text())["seconds_per_call"]
LAT = _U["lattice_construction_rank_3"]
OWN = _U["lattice_element_pairing"]
DOT = _U["zz_dot_12"]
GAPOP = _U["gap_order_matrix_group_so_plus_10_2"]
PERM = _U["gap_permutation_image_527"]
CALL = _U["pari_lllgram_10"]
NMZ = _U["normaliz_cone_6_rays"] + 10 * OWN
LEAF = 2 * _U["zz_matrix_product_4"] + 12 * DOT
INV = 3 * _U["qq_rank_4"] + _U["zz_matrix_product_4"]
ISOM = 256 * DOT + 4 * OWN + LEAF
SHELL = 300 * _U["pari_qfminim_e8_point"] + _U["pari_lllgram_10"]
FLIP = _U["qq_rank_4"] + SHELL
PYROW = DOT

PIN = "tests/test_pinned_lorentzian_cells.py::"
TRA = "tests/test_perfect_cell_transporter.py::"
SIM = "tests/test_perfect_cell_simultaneous_forms.py::"
LCR = "tests/test_large_cell_cone_route.py::"
SYM = "tests/test_target_cell_symmetries.py::"


def _b(seconds: float, basis: str) -> tuple[float, str]:
    return (float(f"{seconds:.3g}"), basis)


def exact_routine(calls: int, lattices: int = 0, extra: float = 0.0, note: str = "") -> tuple[float, str]:
    """Pinned source routine on fixed small exact matrices (dim <= 5): no search branch;
    cost = calls*CALL + lattices*LAT + extra."""
    return _b(calls * CALL + lattices * LAT + extra,
              f"[§Positive-vector shell / §Initial cell] exact-routine branch: {calls} small exact calls, {lattices} owned lattices {note}".strip())


def brute_force_sphere(cases_by_rank: dict[int, int], box: int, calls: int) -> tuple[float, str]:
    """Test-side brute force over (2*box+1)^r integer points per case, priced as r^2+r compiled multiply-adds (DOT each)
    (tests/test_pinned_lorentzian_cells.py:2049-2062, 2201-2212, 2310-2314); port side is the
    pinned Fincke-Pohst/dual-LLL recursion whose output is the (tiny) sphere."""
    points = sum(n * (2 * box + 1) ** r * (r * r + r) for r, n in cases_by_rank.items())
    sec = points * DOT + calls * CALL
    return _b(sec, f"[§Positive-vector shell] CVP sphere branch, cases {cases_by_rank}, box {box}: {points} ops + {calls} calls")


def transporter(n: int, branch: str, stab: int, lattices: int, gap: bool, extra: float = 0.0,
                note: str = "") -> tuple[float, str]:
    """cell_transporter LC:1975 branch work. small_graph: ISOM per call. basis_search:
    INV + leaves (= |Stab| for all_images) * LEAF + n^2 PYROW refinement (n <= 4096)."""
    if branch == "small_graph":
        sec = ISOM
    else:
        sec = INV + stab * LEAF * 2 + n * n * PYROW
    sec += lattices * LAT + (GAPOP if gap else 0.0) + extra
    return _b(sec, f"[§Cell stabilizer and transporter] transporter {branch}, n={n}, |Stab|={stab}, gap={gap} {note}".strip())


def octahedral(case: str) -> tuple[float, str]:
    """3A1+A1(-1), six null rows t+-e_i, |Stab| = 48 (signed permutations).
    graph_limit 0 -> basis_search everywhere; 256 -> small_graph for the witness unless a
    stabilizer is cached (LC:2022-2026), all_images always basis_search. Retained stabilizer
    adds cell_stabilizer (basis_search, 48 leaves) + GAP order. Repeated-row target rejects at
    the Counter/invariant stage (LC:2079)."""
    graph, retained, _shear = case.split("-")
    sec = 4 * LAT + 2 * (INV + 48 * LEAF * 2 + 36 * PYROW)  # complete transporter (x2 calls incl. rejection)
    if graph == "small_graph" and retained == "neither":
        sec += ISOM
    else:
        sec += INV + 48 * LEAF  # witness via basis_search, first certified leaf
    if retained != "neither":
        sec += INV + 48 * LEAF * 2 + GAPOP
    return _b(sec, f"[§Cell stabilizer and transporter] basis_search/small_graph per LC:2022; |Stab|=48, n=6, rank 4, branch {case}")


def large_route(rows: int, normaliz: bool, passes: int, note: str) -> tuple[float, str]:
    """> 4096 or >= 10_000 row large branches (LC:1061, 1417, 1874, 2876): linear passes over the
    rows, each PYROW per row, plus Normaliz on the (deduplicated) 2D cone."""
    sec = LAT + passes * rows * PYROW + (NMZ if normaliz else 0.0)
    return _b(sec, f"[§Large-configuration route / §Cell cone] large-row branch, {rows} rows, {passes} linear passes, normaliz={normaliz}: {note}")


def estimates() -> dict[str, tuple[float, str]]:
    e: dict[str, tuple[float, str]] = {}
    # ---------------- test_target_cell_symmetries
    for g in ("basis_search", "small_graph"):
        for r in ("neither", "source", "target"):
            for s in ("canonical_basis", "sheared_basis"):
                e[f"{SYM}test_octahedral_transport_has_one_witness_and_all_signed_permutations[{g}-{r}-{s}]"] = octahedral(f"{g}-{r}-{s}")
    # ---------------- simultaneous forms: diag(5,-1,-2), 11 rows, |Stab|=4, limits 0 -> basis_search
    e[f"{SIM}test_asymmetric_moment_cell_has_its_complete_integral_transporters"] = transporter(
        11, "basis_search", 4, 2, False, extra=3 * INV + 2 * CALL,
        note="+ positive-vector first shell (anisotropic-direction CVP, 11 vectors) + stabilizer")
    # ---------------- large cell cone route: A2+A1(-2) isotropic cell, 6 rows (< 10_000 ray route)
    e[f"{LCR}test_normaliz_h_presentation_describes_the_same_cone_as_from_rays"] = _b(LAT + 2 * NMZ + 36 * OWN, "[§Cell cone] both cone routes on 6-ray rank-3 cell")
    e[f"{LCR}test_normaliz_route_recovers_exactly_the_configuration_rays"] = _b(LAT + NMZ + 36 * OWN + 36 * DOT, "[§Cell cone] Normaliz H route on 6 rays, 6x6 ratio match")
    e[f"{LCR}test_cell_cone_takes_the_ray_route_below_the_large_cell_threshold"] = _b(LAT + 2 * NMZ, "[§Cell cone] ray route LC:1062 (6 < 10_000), equality check")
    e[f"{LCR}test_positive_ratio_is_sign_sensitive_and_exact"] = _b(7 * DOT, "[§Helpers without search] pure test helper, 7 QQ comparisons")
    # ---------------- perfect cell transporter
    U = "U (isotropic, 2 null rows, |Stab|=2)"
    t = {
        "test_group_generator_complex_does_not_alias_full_domain_cache": _b(LAT + 4 * OWN, "[§Helpers without search] perfect_complex cache lookup only, no traversal"),
        "test_group_only_expansion_discards_facet_orbit_records": _b(2 * (LAT + INV + 4 * OWN), "[§Helpers without search] two expansions with stubbed facets/flip; hash key on 2 rows"),
        "test_orthogonal_group_generators_use_raw_h_rows": _b(LAT + 2 * (SHELL + INV + 6 * LEAF * 2 + FLIP + ISOM + NMZ), "[§Traversal, §Worked tessellations U+A1(-2)] h=1 cell orbit C={e,f,e+f+g}, |G_C|=6 (S3), 1 facet orbit; two traversals (group generators + perfect_complex retain=False), each: initial shell + stabilizer + 1 flip + 1 transporter"),
        "test_owned_row_blocks_release_only_consumed_list_entries": _b(10 * DOT, "[§Helpers without search] block iterator"),
        "test_owned_row_blocks_leave_immutable_sequences_unchanged": _b(10 * DOT, "[§Helpers without search] block iterator"),
        "test_schreier_generator_indices_streams_nested_word_in_action_order": _b(10 * PERM, "[§Helpers without search] word stream"),
        "test_facet_orbit_enumerator_streams_representatives_without_words": _b(6 * PERM, "[§Helpers without search] BFS over 3 facet keys, 1 generator"),
        "test_facet_orbit_materializers_release_consumed_bfs_entries": _b(LAT + INV + 2 * NMZ + 2 * ISOM, f"[§Facets and flips] facet_orbits on {U}"),
        "test_integral_transporter_of_spanning_perfect_cells": transporter(10, "small_graph", 1, 1, False, extra=2 * INV, note="rank-10 simplicial cell; unique transporter"),
        "test_all_integral_transporters_extend_exact_anchor_signatures": transporter(2, "basis_search", 2, 1, False, note=U),
        "test_cell_stabilizer_seeds_certified_gram_automorphisms": transporter(6, "basis_search", 12, 1, False, note="A2+A1(-2) hexagon cell, |Stab|=12"),
        "test_coarse_large_configuration_colours_preserve_exact_stabilizer": transporter(2, "basis_search", 2, 1, False, note=U),
        "test_a3_root_cell_has_its_complete_facet_orbits_and_transport_words": transporter(12, "basis_search", 48, 1, True, extra=NMZ + 14 * 12 * DOT + 14 * 2 * PERM, note="A3+A1(-2), 12 roots, |Stab|=48, 14 facets in 2 orbits"),
        "test_class_sum_pairing_colours_are_isometry_invariants": _b(LAT + 4 * OWN + 16 * DOT, "[§Helpers without search] colour invariants on tiny configuration"),
        "test_cell_stabilizer_seeds_basis_transpositions": transporter(6, "basis_search", 12, 1, False),
        "test_gram_component_sign_flips_are_exact_isometries": _b(LAT + 4 * _U["zz_matrix_product_3"], "[§Helpers without search] sign flips on small Gram"),
        "test_integral_basis_reflection_is_exact_or_rejected": _b(LAT + 2 * _U["zz_matrix_product_3"], "[§Helpers without search] single reflection check"),
        "test_configuration_pairing_diagonal_matches_exact_row_products": _b(8 * DOT, "[§Helpers without search] few rows, exact products"),
        "test_configuration_pairing_diagonal_crosses_temporary_block_boundary": _b(2 * 1025 * DOT, "[§Helpers without search] 1025 rows, test-side Sage vector products"),
        "test_coarse_colours_preserve_changed_configuration_transporters": transporter(2, "basis_search", 2, 2, False),
        "test_configuration_signature_multisets_match_exact_columns": _b(20 * DOT, "[§Helpers without search] small exact columns"),
        "test_packed_pairing_fibres_preserve_exact_position_partition": _b(20 * DOT, "[§Helpers without search] small fibres"),
        "test_empty_stable_configuration_has_empty_permutation_for_each_generator": transporter(2, "basis_search", 2, 1, False),
        "test_permuted_facet_bitsets_preserve_exact_positions": _b(3 * 2000 * PERM, "[§Helpers without search] bitsets up to 2000 positions, one image per position"),
        "test_normaliz_pairing_rows_match_owned_covectors": large_route(10_000, True, 3, "U, 10k duplicated rays -> 2 facets"),
        "test_group_only_facet_keys_do_not_construct_owned_cone": _b(LAT + INV + NMZ, f"raw facet keys on {U}"),
        "test_public_facet_orbit_permutations_keep_owned_cone_lazy": _b(LAT + INV + NMZ, f"facet orbit permutations on {U}"),
        "test_perfect_traversal_keeps_owned_cone_lazy": _b(LAT + SHELL + INV + 2 * LEAF * 2 + 2 * NMZ, "[§Traversal, §Worked tessellations U] h=1 cell {e,f}, |G_C|=2, both facets exit at cusps: no flips"),
        "test_raw_normaliz_facet_route_matches_owned_cone_orbits": _b(LAT + INV + 3 * NMZ, f"both facet routes on {U}"),
        "test_large_facet_keys_use_the_smaller_dense_incidence_side": large_route(10_000, True, 4, "distinct rows, 2 facets, dense incidence"),
        "test_large_perfect_form_hash_skips_quadratic_pairing_graph": large_route(10_000, False, 8, "hash key twice (coarse + refined gcd colours), LC:2876"),
        "test_complex_isometry_uses_cached_coarse_key_for_large_route": _b(LAT + 2 * (INV + 4 * OWN), "[§Helpers without search] two complexes on U, coarse key with stubbed traversal"),
        "test_large_distinct_row_transporter_uses_exact_target_membership": transporter(2, "basis_search", 2, 2, False),
        "test_large_transporter_invariant_rejection_returns_none_not_empty_tuple": _b(LAT + 2 * INV + 2 * 513 * PYROW, "[§Cell stabilizer and transporter] rejection branch LC:2079: rows (i,1) vs (2i,2) have Smith content 1 vs 2, rejected before the quadratic refinement; O(N r) invariants for N=513"),
        "test_large_transporter_matches_joint_gcd_colours_after_basis_change": transporter(8, "basis_search", 8, 2, False),
        "test_large_same_lattice_transporter_skips_counter_identity_probe": transporter(4, "basis_search", 4, 1, False),
        "test_repeated_row_transporter_retains_multiset_certificate": transporter(4, "basis_search", 4, 1, False),
        "test_large_configuration_invariants_include_exact_gcd_colours": _b(LAT + INV, "[§Helpers without search] invariants on U cell"),
        "test_configuration_invariant_cache_is_bounded_and_identity_keyed": _b(LAT + 5 * INV, "[§Helpers without search] cache fill over a few small cells"),
        "test_owned_large_cell_cone_releases_normaliz_engine_cache": large_route(10_000, True, 2, "cell_cone Normaliz route LC:1061"),
        "test_large_cell_owned_cone_bypasses_normaliz_polyhedron_wrapper": large_route(10_000, True, 2, "cell_cone Normaliz route"),
        "test_large_cell_owned_cone_preserves_normaliz_equations": large_route(10_000, True, 2, "cell_cone Normaliz route with equations"),
        "test_large_cell_ray_extraction_bypasses_normaliz_polyhedron_wrapper": large_route(10_000, True, 2, "Normaliz ray extraction"),
        "test_facet_incidence_streaming_preserves_exact_orbits": _b(LAT + 2 * INV + 2 * NMZ + 2 * ISOM, f"[§Facets and flips] facet_orbits twice + permutation orbits on {U}"),
        "test_cell_faces_reuse_owned_halfspace_covectors_exactly": _b(LAT + 2 * NMZ + 36 * OWN, "[§Cell cone] cell faces on small cell"),
        "test_facet_flip_uses_pinned_gap_normals": _b(LAT + 2 * GAPOP + NMZ + FLIP, "[§Facets and flips] flip_across diag(1,-1,-1) 3-ray cell: kernel + positive-vector search; libgap NullspaceMat"),
        "test_rank_three_minimized_halfspaces_pair_against_basis_elements": _b(LAT + 2 * NMZ + 2 * SHELL, "[§Cell cone] rank-3 isotropic cell, 4 facets, both modes"),
    }
    for k, v in t.items():
        e[TRA + k] = v
    # ---------------- pinned lorentzian cells (exact pinned routines on fixed small matrices)
    p = {
        "test_gap_lll_transformation_uses_row_congruence_orientation": exact_routine(10, 0, GAPOP, "(libgap LLL)"),
        "test_dual_lll_candidates_sort_in_reduced_column_coordinates": exact_routine(10),
        "test_pinned_positive_gcd_dot_matches_exact_coefficients": exact_routine(5),
        "test_outside_ray_span_normalizes_each_source_direction_independently": exact_routine(5),
        "test_initial_cell_nullspace_uses_pinned_gap_row_order": exact_routine(5, 0, GAPOP),
        "test_pinned_column_reduction_matches_greedy_rank_increase": exact_routine(3 * 8 * 6, note="(random 8 per size, rank checks)"),
        "test_pinned_positive_diagonal_directions_follow_source_scan_and_l1_norm": exact_routine(10),
        "test_pinned_positive_eigen_directions_preserve_source_row_orientation": exact_routine(10),
        "test_pinned_positive_eigen_directions_rejects_invalid_double_gram": exact_routine(5),
        "test_pinned_approximation_iterator_rejects_nonfinite_float_directions": exact_routine(5),
        "test_pinned_approximation_iterator_advances_tied_steps_in_source_order": exact_routine(10),
        "test_source_positive_candidate_lift_uses_transposed_unimodular_basis": exact_routine(10),
        "test_definite_dual_lll_uses_inverse_of_sage_transform": exact_routine(10),
        "test_simple_indefinite_reduction_certifies_each_l1_descent": exact_routine(20),
        "test_compute_indefinite_reduction_rejects_false_lll_certificate": exact_routine(10),
        "test_nondegenerate_reduction_rejects_false_operation_certificate": exact_routine(10),
        "test_indefinite_reduction_integrates_source_radical_and_completion": exact_routine(30),
        "test_indefinite_gram_schmidt_rejects_invalid_form_or_basis": exact_routine(10),
        "test_indefinite_lll_nearest_integer_uses_source_half_tie_rule": exact_routine(10),
        "test_source_integer_nullspace_reduces_and_completes_radical": exact_routine(20, 0, GAPOP),
        "test_source_subspace_completion_int_returns_unimodular_complements": exact_routine(20),
        "test_source_v2_honors_critical_norm_and_strictness": exact_routine(30, note="(isotropic/anisotropic shortcut branches)"),
        "test_source_integral_vector_allmeth_v2_returns_exact_positive_witness": exact_routine(40),
        "test_source_remove_fraction_vector_returns_primitive_integral_coordinates": exact_routine(5),
        "test_source_integral_reduction_rejects_false_basis_certificates": exact_routine(10),
        "test_source_positive_integral_vector_contract_certifies_exact_norm": exact_routine(10),
        "test_source_integral_vector_allmeth_rejects_invalid_gram_before_reduction": exact_routine(5),
        "test_source_integral_vector_allmeth_selects_isotropic_shortcut_exactly": exact_routine(20),
        "test_integral_positive_vector_consumer_uses_pinned_diagonalization": exact_routine(20),
        "test_pinned_symmetric_diagonalization_returns_exact_congruence_witness": exact_routine(20),
        "test_pinned_diagonalizer_handles_deterministic_rational_symmetric_samples": exact_routine(4 * 10 * 3, note="(10 random per size)"),
        "test_pinned_nullspace_preserves_gap_row_basis_order": exact_routine(4 * 12 * 2 + 10, 0, GAPOP),
        "test_source_random_perturbations_match_pinned_matrices": exact_routine(20),
        "test_source_indefinite_short_positive_vector_integrates_pinned_search_loop": exact_routine(40, note="(scripted search loop)"),
        "test_source_indefinite_short_positive_vector_certifies_search_witnesses": exact_routine(20),
        "test_source_positive_search_best_matches_pinned_termination_state": exact_routine(30),
        "test_source_positive_perturbation_reset_uses_strict_factor_1000": exact_routine(20),
        "test_source_direct_graver_improvement_matches_literal_ordered_passes": exact_routine(30),
        "test_source_direct_graver_improvement_rejects_false_descent_prefix": exact_routine(10),
        "test_source_graver_k_basis_matches_support_and_sign_enumeration": exact_routine(40, note="(dim<=5 Graver bases, <= 242 rows each)"),
        "test_source_random_matrix_perturbation_dispatches_both_pinned_branches": exact_routine(5),
        "test_source_perturbation_matrices_are_unimodular_in_small_dimensions": exact_routine(33 + 384 // 4 + 40, note="(<=432 signed perms + elementary, one determinant each)"),
        "test_pinned_nearest_integer_rounds_half_ties_toward_positive_infinity": exact_routine(2),
        "test_source_negative_vector_respects_radical_and_diagonal_order": exact_routine(20),
        "test_source_negative_vector_certifies_congruence_random_indefinite_forms": exact_routine(4 * 8 * 3),
        "test_source_negative_vector_certifies_rational_congruence_forms": exact_routine(4 * 2 * 4 * 3),
        "test_source_mid_value_matches_pinned_continued_fraction_sequence": exact_routine(1925 // 10, extra=1925 * 5 * DOT, note="(1925 interval pairs, CF levels)"),
        "test_quadratic_descent_prefix_matches_pinned_get_alpha": exact_routine(10),
        "test_flip_upper_bound_uses_rational_isotropic_exit_and_rejects_coincident_exit": exact_routine(10),
        "test_pinned_sublattice_basis_reduction_preserves_row_lattice": exact_routine(20),
        "test_pinned_integer_covector_kernel_matches_exact_contract": exact_routine(20),
        "test_pinned_sublattice_reduction_on_deterministic_full_rank_samples": exact_routine(12 * 4),
        "test_primitive_integral_direction_matches_fraction_multiplier": exact_routine(10),
        "test_ambient_vectors_reconstructs_exact_shell_and_filters_squares": exact_routine(20, 2),
        "test_rational_direction_pairing_bound_scales_linearly": exact_routine(5),
        "test_initial_vector_search_advances_integral_bound_without_rescanning": exact_routine(20, 1),
        "test_isotropic_multiplier_uses_dual_lll_fixed_distance_solver": exact_routine(20, 1, note="(isotropic branch: fixed-distance sphere)"),
        "test_total_multiplier_keeps_bounded_cvp_and_pinned_order": exact_routine(20, 1, note="(anisotropic branch: bounded CVP)"),
        "test_fixed_distance_coordinate_solutions_match_exact_brute_force": brute_force_sphere({1: 2, 2: 3}, 5, 30),
        "test_fixed_distance_dual_lll_frame_reconstructs_exact_sphere": brute_force_sphere({1: 10, 2: 10, 3: 10}, 8, 60),
        "test_fixed_distance_basis_inverse_is_exact_and_two_sided": exact_routine(10),
        "test_bounded_distance_dual_lll_solutions_match_exact_brute_force": brute_force_sphere({1: 1, 2: 1, 3: 1}, 5, 10),
        "test_positive_ldl_factorization_is_exact_and_rejects_nonpositive_pivots": exact_routine(15),
        "test_multiplier_scales_cvp_center_and_squared_distance_bound": exact_routine(20, 1),
        "test_positive_vector_constructor_uses_cpp_semantic_cvp_center": exact_routine(20, 2),
        "test_positive_vector_shell_respects_explicit_coordinate_labels": exact_routine(20, 2),
        "test_positive_vector_shell_scales_pairing_bound_with_rational_direction": exact_routine(20, 2),
        "test_positive_vector_shell_matches_exact_nonzero_center_enumeration": exact_routine(20, 1, note="(first shell: 3 vectors, anisotropic CVP)"),
    }
    for k, v in p.items():
        e[PIN + k] = v
    return e


if __name__ == "__main__":
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[2]
    files = ("test_pinned_lorentzian_cells", "test_perfect_cell_transporter", "test_perfect_cell_simultaneous_forms",
             "test_large_cell_cone_route", "test_target_cell_symmetries")
    nodes = {l.strip() for l in (root / "ci/test_budget_models/nodes_initial.txt").read_text().splitlines()
             if any(f"tests/{f}.py::" in l for f in files)}
    est = estimates()
    assert set(est) == nodes, (sorted(nodes - set(est)), sorted(set(est) - nodes))
    by_fn: dict[str, list[float]] = {}
    for k, (s, _) in est.items():
        by_fn.setdefault(k.split("[")[0], []).append(s)
    print(len(est), "nodes")
    for fn, v in sorted(by_fn.items()):
        print(f"{fn}: n={len(v)} min={min(v)} median={median(v)} max={max(v)}")
