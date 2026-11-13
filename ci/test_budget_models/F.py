"""A-priori time budgets for agent F's test scope, derived from the mathematics.

Theory: references/complexity/discriminant-groups-and-congruence-orbits.md (cited as DG §<section>).
Generator costs of O(L) and of parabolics are taken from the owners' files
orthogonal-group-generators.md (B) and boundary-components.md (B); the constants
below state the bound assumed from each branch. No measured duration is used.
"""

SAFETY = 3.0
C_MOR = 0.03  # DG §Unit costs: one morphism/marking operation
C_LAT = 0.5  # DG §Unit costs: formed lattice with Smith form
PURE = 0.2  # pure-Python fixture/CI check, no Sage computation


def c_ss(degree, gens=20):
    """DG §Unit costs: GAP randomized Schreier-Sims, nearly linear in degree (Seress Thm 4.5.5 / Sec. 4.5)."""
    from math import log2
    return 1 + 1e-7 * gens * degree * max(1.0, log2(degree)) ** 2


def c_tqf(order_a):
    """DG §Unit costs: O(q_L) via Sage fqf_orthogonal, per primary part."""
    return 2 + 0.02 * order_a


# Generator branches of O(L): agent B's model, imported unchanged
# (references/complexity/orthogonal-group-generators.md §Branches, §Eichler route, §Lorentzian route;
# Lorentzian cell counts and per-cell cost from perfect-domains-and-lorentzian-cells.md, via B.t_traverse).
import importlib.util as _ilu
import json
from pathlib import Path as _Path

_spec = _ilu.spec_from_file_location("estimates_B", _Path(__file__).with_name("B.py"))
B = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(B)


def _diag_sum(*blocks):
    n = sum(len(b) for b in blocks)
    g = [[0] * n for _ in range(n)]
    k = 0
    for b in blocks:
        for i, row in enumerate(b):
            for j, x in enumerate(row):
                g[k + i][k + j] = x
        k += len(b)
    return g


def _u(s=1):
    return [[0, s], [s, 0]]


E8 = [[2, -1, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, -1], [0, 0, -1, 2, -1, 0, 0, 0],
      [0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0], [0, 0, 0, 0, 0, -1, 2, 0], [0, 0, -1, 0, 0, 0, 0, 2]]
E8M2 = [[-2 * x for x in row] for row in E8]
G_2U = _diag_sum(_u(), _u())
G_3U = _diag_sum(_u(), _u(), _u())
G_4U = _diag_sum(_u(), _u(), _u(), _u())
G_U6 = _diag_sum(_u(), [[6]])
G_3U_A1_3 = _diag_sum(_u(), _u(), _u(), [[6]])
G_2U_U3 = _diag_sum(_u(), _u(), _u(3))
G_2U_M2_M6 = _diag_sum(_u(), _u(), [[-2]], [[-6]])
G_2U2_A2 = _diag_sum(_u(2), _u(2), [[-2, 1], [1, -2]])
G_U0 = _diag_sum(_u(), [[0]])
G_UUA1 = _diag_sum(_u(), _u(), [[-2]])
G_N = _diag_sum(_u(), _u(2), E8M2)
G_LEVI_E1 = _diag_sum(_u(2), E8M2)  # e1^perp/e1 for e1 in U
G_LEVI_E3 = _diag_sum(_u(), E8M2)  # e3^perp/e3 for e3 in U(2), up to the scale of U(2) (B doc §Eichler route)


def gens(gram):
    """B.t_group: O(L) generators on the branch selected by signature, rank and Witt regime (recursive.py:1134-1188)."""
    return B.t_group(B.inv(gram))


def gbranch(gram):
    i = B.inv(gram)
    return f"O(L) route={B.route(i)} n={i['n']} sig=({i['p']},{i['q']}) |A|={i['A']} [B: orthogonal-group-generators.md]"


def parabolic(rank):
    """Levi + unipotent generators of an isotropic stabilizer (B: boundary-components.md), unimodular Levi."""
    return 2.0 + 0.5 * rank


def finite_split(gram, levi_rank, order_f, orbits, n, m, per_rep=2 * C_MOR):
    """DG §Orbit splitting, finite preimage branch: O(L) gens, image F in GL_n(Z/m), rho(P), double cosets."""
    return (C_LAT + gens(gram) + 10 * n * n * C_MOR + c_ss(min(m**n, order_f))
            + parabolic(levi_rank) + 1e-4 * order_f + orbits * per_rep)


def budget(raw, basis):
    return (round(SAFETY * raw, 1), f"{basis}; raw {raw:.1f}s x safety {SAFETY}")


W = "witnesses/test_congruence_orbit_splitting_across_markings.sage::"
O = "witnesses/test_orbits.sage::"
D = "witnesses/test_degenerate_ambient_consumers.sage::"


_AL1 = [[-239360, -11008, -11520], [-11008, -504, -528], [-11520, -528, -553]]
_AL2 = [[13951232, 6751744, 425216], [6751744, 3267528, 205784], [425216, 205784, 12959]]
_H1 = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]
_H2 = [[2, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]


def _estimates_first():
    e = {}
    # --- witnesses: 2U modulo 3, F = O+_4(3) of order 1152 (DG §Congruence, §Orbit splitting)
    for name, orbits, basis in (
        ("test_based_null_vectors_split_into_thirty_two_gamma_three_orbits", 32, "coslice vector, |F|=1152, |rho(P)|=36"),
        ("test_null_lines_split_into_sixteen_gamma_three_orbits", 16, "arrow line, |rho(P)|=72, + one Gamma(3) transporter"),
        ("test_null_planes_split_into_eight_gamma_three_orbits_with_stabilizers", 8, "arrow plane, |rho(P)|=144, + transporter + Schreier stabilizer index<=144"),
        ("test_null_flags_split_into_thirty_two_gamma_three_orbits", 32, "Fun([2]) flag, Borel |rho(P)|=36"),
        ("test_hyperbolic_summands_split_into_seventy_two_gamma_three_orbits", 72, "nondegenerate summand U, |rho(P)|=16"),
        ("test_based_pairs_of_prescribed_gram_split_into_288_gamma_three_orbits", 288, "coslice under U, |rho(P)|=4"),
    ):
        raw = finite_split(G_2U, 4, 1152, orbits, 4, 3, per_rep=4 * C_MOR)
        if "planes" in name:
            raw += 144 * C_MOR * 4 + parabolic(4)
        e[W + name] = budget(raw, f"DG §Orbit splitting finite-preimage branch: 2U mod 3, {basis}, {orbits} double cosets; " + gbranch(G_2U))
    raw = C_LAT + gens(G_2U) + 2 * (c_ss(81) + parabolic(4) + 1152 * 1e-3)
    e[W + "test_a_based_vector_and_its_negative_under_gamma_three_and_gamma_two"] = budget(raw, "DG §Congruence: two transporter searches through finite images O+_4(3) (1152) and O+_4(2) (72)")
    raw = C_LAT + gens(G_2U) + c_ss(81) + parabolic(4) + 36 * 5 * C_MOR
    e[W + "test_based_vector_stabilizer_in_gamma_three_fixes_the_vector"] = budget(raw, "DG §Congruence Schreier generators: index |rho(Stab e)|=36, ~5 gens each")
    raw = finite_split(G_2U, 4, 1152, 5, 4, 3) + c_ss(81) + parabolic(4)
    e[W + "test_null_vectors_under_the_preimage_of_a_non_normal_stabilizer"] = budget(raw, "DG §Orbit splitting finite preimage, S=rho(P_e) order 36: 5 double cosets; extra subgroup order")
    raw = finite_split(G_3U, 6, 40320, 630, 6, 2, per_rep=8 * C_MOR)
    e[W + "test_complete_flags_of_rank_three_split_into_630_gamma_two_orbits"] = budget(raw, "DG §Orbit splitting: 3U mod 2, F=O+_6(2) order 40320, Borel image 64, 630 double cosets, 3 residue spans per rep")
    raw = C_LAT + gens(G_U6) + 2 * c_ss(216) + 2 * 400 * 2 * C_MOR + 4 * (parabolic(3) + 20 * C_MOR) + c_tqf(6)
    e[W + "test_joint_level_and_discriminant_image_on_u_plus_six"] = budget(raw, "DG §Congruence: U+<6> sig (2,1) Lorentzian gens, |A_L|=6, |O(q)|=2; joint image mod 6 enumerated (<=400 elems, section+2 memberships each); 4 splits")
    # --- witnesses/test_orbits.sage
    split16 = finite_split(G_2U, 4, 1152, 16, 4, 3)
    e[O + "test_gamma_three_splits_null_lines_of_two_hyperbolic_planes_into_sixteen_orbits"] = budget(split16 + 2 * (parabolic(4) + c_ss(81)), "DG §Orbit splitting: 2U mod 3 lines, 16 cosets, + 2 transporters")
    e[O + "test_gamma_two_splits_null_lines_of_two_hyperbolic_planes_into_nine_orbits"] = budget(finite_split(G_2U, 4, 72, 9, 4, 2) + parabolic(4), "DG §Orbit splitting: 2U mod 2, F=O+_4(2) order 72, 9 cosets, + 1 transporter")
    e[O + "test_gamma_three_transporters_between_null_lines_of_two_hyperbolic_planes"] = budget(C_LAT + gens(G_2U) + 2 * (parabolic(4) + c_ss(81)), "DG §Congruence: two Gamma(3) transporters through F order 1152")
    e[O + "test_gamma_three_line_stabilizer_has_schreier_generators_of_index_72"] = budget(C_LAT + gens(G_2U) + 2 * parabolic(4) + c_ss(81) + 2 * 72 * 6 * C_MOR, "DG §Congruence Schreier generators: index 72, ~6 parabolic gens, computed twice")
    raw = C_LAT + gens(G_3U_A1_3) + parabolic(7) + c_tqf(6) + c_ss(6)
    e[O + "test_rank_three_discriminant_image"] = budget(raw, "DG §Discriminant image of stabilizer: 3U+A1(3), |A_L|=6, |O(q)|=2, plane parabolic in rank 7")
    raw = C_LAT + gens(G_2U_U3) + parabolic(6) + c_tqf(9) + c_ss(9)
    e[O + "test_rank_three_null_stabilizer_of_twisted_hyperbolic_sum_fixes_one_isotropic_line"] = budget(raw, "DG §Discriminant image: 2U+U(3), A_L=(Z/3)^2, |O(q)|=4, maximal isotropic parabolic")
    # --- degenerate ambient U+<0> (radical branch, marked_lattices.py:386, recursive.py:1158-1161)
    e[D + "test_degenerate_vector_classes_of_square_two_and_fifty"] = budget(C_LAT + gens(G_U0) + 2 * 10 * C_MOR, "DG §Marking classes: degenerate branch, classes from divisors of 50 (<=4) and (Z/d)^*/+-1")
    e[D + "test_congruence_image_and_kernel_generators"] = budget(C_LAT + gens(G_U0) + c_ss(27) + 72 * 4 * C_MOR, "DG §Congruence: U+<0> mod 3 image order 72, Schreier kernel gens index 72")
    e[D + "test_square_zero_classes_are_a_lazy_set"] = budget(C_LAT + gens(G_U0) + 15 * 5 * C_MOR, "DG §Marking classes: lazy enumeration, 8 classes, d<=5")
    # --- Dawes u2u2a2 incidence witness
    raw = C_LAT + gens(G_2U2_A2) + c_tqf(16) + c_tqf(3) + c_ss(48) + 15 * (parabolic(6) + c_ss(48)) + 18 * 10 * C_MOR
    e["witnesses/test_dawes_u2u2a2_incidence_records.py::test_u2u2a2_stable_building_records"] = budget(raw, "DG §Marking classes + §Orbit splitting char-kernel branch: 2U(2)+A2, |A_L|=48, |O(q)|=144; 9 line + 6 plane cusps, 18 incidences, one parabolic per cusp")
    # --- semisimple Voronoi cubic division (E: clifford-and-hermitian-lattices.md owns the algorithm)
    sv = "witnesses/test_semisimple_voronoi_cubic_division.sage::"
    e[sv + "test_split_model_is_a_homomorphism"] = budget(1.0, "two 9x9 rational products and 3x3 matrices over a cubic field")
    e[sv + "test_units_are_integral_of_reduced_norm_plus_minus_one_and_one_has_infinite_order"] = budget(60.0, "E: semisimple Voronoi over M_3(R)-cone (dim 6), maximal order of degree-3 division algebra disc 2^2*7^2: bounded by 60 s of perfect-form steps")
    # --- tests/ pure Python
    for f in ("test_acceptance_summary_fails_when_a_matrix_shard_has_no_artifact", "test_junit_node_identity_is_not_replaced_by_an_equal_result_count",
              "test_batch_runner_continues_after_malformed_batch_junit", "test_acceptance_summary_accepts_a_present_shard_with_a_junit_report",
              "test_acceptance_summary_rejects_mixed_port_or_provider_revisions", "test_acceptance_summary_rejects_missing_revision_metadata",
              "test_acceptance_summary_rejects_a_failed_family_matrix", "test_acceptance_summary_rejects_an_empty_junit_report",
              "test_acceptance_summary_reports_malformed_junit_without_crashing", "test_acceptance_summary_reports_a_failed_plan_without_valid_matrix_json",
              "test_acceptance_summary_rejects_empty_or_malformed_matrix_includes"):
        e["tests/test_acceptance_ci.py::" + f] = budget(PURE, "no mathematics: CI summary script on tmp JSON/JUnit; subprocess monkeypatched")
    a = "tests/test_allcock_fixture_identity.py::"
    parse = 0.3  # 4.9 MB JSON, 8595 entries
    for f, k in (("test_fixture_keeps_every_case", 1), ("test_id_is_the_ambiguous_lattice_only_label", 1),
                 ("test_retained_source_identity_is_unique", 1), ("test_composed_id_is_injective_and_covers_every_case", 1),
                 ("test_no_two_cases_share_content_so_no_occurrence_counter_is_needed", 3), ("test_every_case_cites_its_source", 1),
                 ("test_mathematical_data_is_preserved_per_case", 1)):
        e[a + f] = budget(parse * k, "no mathematics: parse 4.9 MB fixture (8595 cases); sha256 per case where stated")
    for c in ("allcock_rk3_4-0-4", "allcock_rk3_2-2-2"):
        e[a + f"test_the_recorded_sample_cases_resolve_to_one_lattice[{c}]"] = budget(parse, "fixture lookup only")
    s = "tests/test_allcock_same_genus_nonisometric.py::"
    e[s + "test_fixture_holds_the_two_printed_lattices_with_equal_elementary_divisors"] = budget(parse, "fixture lookup only")
    e[s + "test_same_genus_lattices_with_different_weyl_groups_are_not_isometric"] = budget(
        2 * C_LAT + 2 * B.t_isometry(_AL1, _AL2)[0], "B isometry-testing.md §Branches lorentzian (is_isometric + Isom().is_empty, 2 decisions): " + B.t_isometry(_AL1, _AL2)[1])
    i = "tests/test_issue15_defect_audit.py::"
    e[i + "test_isotropic_mode_on_an_anisotropic_lattice_refuses_instead_of_searching_forever"] = budget(C_LAT + 1, "DG §Non-termination: anisotropy check of binary <2,-6> (disc 12 non-square), immediate refusal")
    e[i + "test_a_rational_coordinate_is_refused_rather_than_truncated"] = budget(PURE, "scalar conversions")
    e[i + "test_the_rank_two_branch_reports_every_wall_it_found"] = budget(2 * (C_LAT + 2), "A vinberg-and-edgewalk: anisotropic binary cycle, 2 walls, disc 8 and 20")
    e[i + "test_a_split_rank_two_lattice_is_not_reported_reflective_with_a_wall"] = budget(C_LAT + 2, "A: split binary U, square-discriminant branch, finite O(U)")
    p = "tests/test_configuration_pairings.py::"
    for c in ("1-3", "1208925819614629174706176-3", "1-1208925819614629174706176"):
        e[p + f"test_coloured_hyperbolic_neighbour_multiplicities_in_every_integral_basis[{c}]"] = budget(PURE, "4x4 pairings of 2x2 integer data, 80-bit integers")
    e[p + "test_hyperbolic_plane_pairings_preserve_integral_change_of_basis"] = budget(1.0, "80x80 pairing columns, 3 passes = 19200 small products")
    e[p + "test_configuration_pairing_blocks_have_a_bounded_exact_working_set"] = budget(PURE, "arithmetic on sizes")
    for f in ("test_isotropic_cell_invariants_preserve_integral_change_of_basis", "test_large_coarse_hash_preserves_gcd_colours_under_integral_isometry",
              "test_nonnegative_square_face_has_the_exact_pairing_multisets"):
        e[p + f] = budget(2 * C_LAT + 0.5, "C perfect-domains: hash key of a cell with 2-3 rows in rank 2")
    d2 = "tests/test_dawes_2U_minus2_minus6_stabilizers.py::"
    e[d2 + "test_discriminant_orthogonal_group_has_order_two"] = budget(C_LAT + c_tqf(4) + c_tqf(3), "DG §Discriminant group: A_L=Z/2+Z/6, |O(q)|=2")
    e[d2 + "test_stable_group_has_index_two_in_o_plus_generated_by_sigma_w"] = budget(C_LAT + gens(G_2U_M2_M6) + c_tqf(12) + 2 * c_ss(12), "DG §Discriminant image: O+ gens via approximate model (sig (2,4), Witt 2), two images in O(q) of order 2")
    e[d2 + "test_line_and_plane_stabilizers_surject_onto_the_discriminant_image"] = budget(C_LAT + gens(G_2U_M2_M6) + c_tqf(12) + 4 * (parabolic(6) + 2 + c_ss(12)), "DG §Discriminant image: 4 setwise parabolics (2 lines, 2 planes) in O+, O+ index-2 character, |A_L|=12")
    l2 = "tests/test_dawes_L2_boundary_witnesses.py::"
    e[l2 + "test_L2_isotropic_discriminant_classes"] = budget(PURE, "12 rational evaluations")
    e[l2 + "test_L2_plane_stabilizer_transvections"] = budget(C_LAT + 4 * 6 * 7 * C_MOR, "4 Eichler transvections checked on 36 basis pairs")
    e["tests/test_dawes_scaled_planes_discriminant_flags.py::test_scaled_planes_discriminant_flags_match_the_published_vertex_colours"] = budget(
        C_LAT + 30 * C_MOR + 6 * 0.2 + 0.5, "DG §Discriminant group: 9 classes in (Z/2)^4 part of |A_L|=48, 6 planes, canonical label of 15-vertex graph")
    m = "tests/test_maximal_split_isotropic_stabilizer.py::"
    e[m + "test_pointwise_maximal_stabilizer_in_four_u_has_the_full_skew_parameter_module"] = budget(C_LAT + gens(G_4U) + parabolic(8) + 6 * 40 * C_MOR, "B boundary-components: pointwise maximal isotropic stabilizer in 4U, unipotent radical rank 6 (skew 4x4)")
    e[m + "test_special_orthogonal_group_has_two_moving_maximal_isotropic_orbits_in_four_u"] = budget(C_LAT + gens(G_4U) + c_ss(4) + 4 * parabolic(8) + 40 * C_MOR, "DG §Orbit splitting char-kernel (det character), |Q|=2, 2 orbits, 3 transporters in 4U")
    t = "tests/test_marked_subgroup_transporter.py::"
    e[t + "test_fixed_null_vector_transporter_respects_special_orthogonal_subgroup"] = budget(C_LAT + gens(_u()) + 4 * 10 * C_MOR, "split binary U, O(U) order 4: finite enumeration")
    for c in ("based-4", "moving-2", "flag-2"):
        e[t + f"test_joint_determinant_and_discriminant_conditions_use_the_actual_image[{c}]"] = budget(C_LAT + gens(_u(3)) + c_tqf(9) + c_ss(9) + 6 * 10 * C_MOR, "DG §Orbit splitting char-kernel: U(3) binary split, A_L=(Z/3)^2, |Q|=4, <=4 cosets")
    h = "tests/test_hyperbolic_odd_isometry_witness.py::"
    e[h + "test_hand_change_of_basis_is_a_witness"] = budget(PURE, "4x4 matrix product")
    e[h + "test_returned_witness_carries_one_gram_matrix_to_the_other"] = budget(2 * C_LAT + B.t_isometry(_H1, _H2)[0], "B isometry-testing.md §Branches: " + B.t_isometry(_H1, _H2)[1])
    e["tests/test_unsupported_quotient_group.py::test_line_stabilizer_with_an_anisotropic_rank_four_quotient_is_unsupported"] = budget(
        2 * C_LAT + 3, "DG §Non-termination: Levi O(K), K sig (2,2) Witt 0 not quaternion model -> UnfinishedCapability (recursive.py:1175); Hasse invariants at 2,7")
    raw = 2 * (C_LAT + gens(G_UUA1) + c_tqf(2) + 2 * 2 + 3 * parabolic(5))
    e["tests/test_orbit_counts_basis_change.py::test_U_U_A1_negative_counts_are_basis_independent"] = budget(raw, "DG §Marking classes Eichler criterion: U+U+A1(-1), |A_L|=2; counts (1,2,1,1,1), two bases")
    lr = "tests/test_lorentzian_reduction_image_orders.py::"
    for c in B._load("lorentzian_stabilizers_cases.json"):
        if not c["id"].startswith("LORENTZ_OrthogonalGroup"):
            continue
        i = B.inv(c["gram"])
        cells = max(B.cells_estimate(i), -(-c["num_generators"] // 3))
        k = c["num_generators"]
        raw = C_LAT + B.t_traverse(i, cells) + sum(2 * B.t_reduction(10, p, k) for p in (2, 3))
        e[lr + f"test_reduction_images_agree_and_are_orthogonal_groups_where_nondegenerate[{c['id']}]"] = budget(
            raw, f"B orthogonal-group-generators.md §Lorentzian route cells~{cells} n=10 |A|={i['A']}; DG §Congruence + B §Finite images: 4 reduction orders (p=2,3, constructive recognition), {k} gens")
    en = "tests/test_enriques_line_stabilizers.py::"
    for c, levi in (("I_1_1", "U(2)+E8(-2), |A|=2^10"), ("I_1_2", "U+E8(-2), |A|=2^8")):
        raw = C_LAT + gens(G_N) + gens(G_LEVI_E1 if c == "I_1_1" else G_LEVI_E3) + parabolic(12) + c_tqf(1024) + c_ss(1024, 40)
        e[en + f"test_line_stabilizer_generators_generate_the_stabilizer[{c}]"] = budget(
            raw, gbranch(G_N) + "; Levi " + gbranch(G_LEVI_E1 if c == "I_1_1" else G_LEVI_E3) + f"; DG §Discriminant image: N sig (2,10), |A_N|=1024, |O(q_N)|=|O+_10(2)|=46998591897600; Levi {levi} Lorentzian; index 527 = isotropic points")
    return e


# ---------------------------------------------------------------- second scope (missing_F.txt)
DM = "DG §Markings on degenerate lattices"
IND = "DG §Markings by isotropic induction"
SH = "DG §Pointwise stabilizer of a null vector and its shell decision"
NZ = "DG §Normalizers of subgroups of O(L)"
EQ = "DG §Equivariant markings"
G_U_U2 = _diag_sum(_u(), _u(2))
G_GLUED = [[0, 2, 0, 0], [2, 0, 1, 0], [0, 1, 0, 2], [0, 0, 2, 0]]
G_U_M2 = _diag_sum(_u(), [[-2]])
G_8716 = [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 10, 5], [0, 0, 0, 0, 5, 10]]
G_DEG1 = [[2, 0], [0, 0]]
G_DEG2 = [[2, 0, 0], [0, 0, 0], [0, 0, 0]]
G_DEG3 = [[2, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
G_PELL = [[1, 0], [0, -2]]
U_ELT = B.U_ELT


C_GAPZ = 1.0  # GAP order/IsFinite of a matrix group over ZZ of degree <= 5 (DG §Unit costs: invariant-form method)
C_MM = 1e-4  # one dense product of matrices of size <= 8 over ZZ, QQ or GF(3) in Sage (DG §Unit costs)
C_OP = 0.5  # one equivariant stabilizer/transporter step after the centralizer (DG §Equivariant markings)


def _radical_rank(gram):
    return sum(1 for row in gram if not any(row))


def _per_call(gram):
    """DG §Markings on degenerate lattices: one stabilizer or transporter = O(t*c + t^2) morphism operations."""
    t = _radical_rank(gram)
    c = len(gram) - t
    return (t * c + t * t + 4) * C_MOR


def deg_ops(ops):
    """Sum the traced public calls of one test body. Each lattice pays C_LAT + its O(L) generators once (B's model)."""
    total, seen, parts = 0.0, set(), []
    for op in ops:
        kind = op[0]
        if kind in ("pi", "stab"):
            _k, g, n = op
            key = json.dumps(g)
            if key not in seen:
                seen.add(key)
                total += C_LAT + gens(g)
            total += n * _per_call(g)
            parts.append(f"{n} {kind} on {g}")
        elif kind == "can":
            _k, g, h, n = op
            key = json.dumps(g)
            if key not in seen:
                seen.add(key)
                total += C_LAT + gens(g)
            total += n * 2 * (2 * h + 1) * _per_call(g)
            parts.append(f"{n} canonicalizations at height {h} (2(2h+1) candidate transporters each)")
        elif kind == "cong":
            _k, t, moduli = op
            total += sum(c_ss(m ** t) for m in moduli)
            parts.append(f"radical congruence images t={t} m={list(moduli)}")
        elif kind == "mgZ":
            total += op[1] * C_GAPZ
            parts.append(f"{op[1]} MatrixGroup orders over ZZ")
        elif kind == "shear":
            total += op[1] * 10 * C_MOR
            parts.append(f"{op[1]} shear kernels")
        elif kind == "closure":
            _k, size, k = op
            total += size * k * C_MM
            parts.append(f"closure of {size} elements x {k} generators")
    return total, "; ".join(parts)


D1, D2, D3 = G_DEG1, G_DEG2, G_DEG3
DP = [[0, 0], [0, 2]]
DPT = {  # tests/test_degenerate_pointed_markings.py, traced call by call
    "test_transports_one_zero_to_one_one_by_an_integral_unimodular_shear": [("pi", D1, 1)],
    "test_transports_five_one_to_five_six_by_the_same_shear": [("pi", D1, 1)],
    "test_equal_norm_content_and_divisibility_do_not_transport_five_one_to_five_two": [("pi", D1, 1)],
    "test_pointwise_stabilizer_of_one_zero_is_the_order_two_reflection_of_the_radical": [("stab", D1, 1), ("mgZ", 1)],
    "test_public_pointed_consumers_reach_the_same_verdicts": [("pi", D1, 4), ("stab", D1, 1), ("mgZ", 1)],
    "test_marking_transports_in_a_presentation_where_the_radical_is_not_last": [("pi", DP, 1)],
    "test_unmarked_radical_automorphism_route_is_preserved": [("stab", D1, 1), ("mgZ", 1)],
    "test_rank_two_radical_transport_needs_a_nonidentity_integral_D": [("pi", D2, 1)],
    "test_rank_two_radical_refuses_when_r_vanishes_modulo_content_but_s_does_not": [("pi", D2, 1)],
    "test_pure_shear_transport_when_d_minus_identity_already_lies_in_the_content_module": [("pi", D2, 1)],
    "test_pointwise_stabilizer_contains_the_nontrivial_shear_kernel": [("stab", [[2, 0, 0], [0, 2, 0], [0, 0, 0]], 1)],
    "test_marked_orbit_is_the_native_lazy_set_not_a_list": [("pi", D1, 3)],
    "test_marked_orbit_of_a_radical_vector_handles_the_c_zero_case": [("pi", D1, 4)],
    "test_marked_orbit_action_exposes_the_native_orbit_surface": [("pi", D1, 2)],
    "test_marked_orbit_quotient_is_the_orbit_set_and_refuses_to_be_listed": [("pi", D1, 3)],
    "test_isotropic_vector_off_the_radical_is_not_classified_as_radical": [("pi", D1, 0), ("pi", [[0, 1, 0], [1, 0, 0], [0, 0, 0]], 0)],
    "test_radical_stabilizer_in_a_permuted_basis": [("stab", DP, 1)],
    "test_radical_stabilizer_in_a_mixed_basis": [("stab", [[2, 0, 0], [0, 0, 0], [0, 0, 2]], 1)],
    "test_radical_stabilizer_of_the_zero_vector_is_the_whole_group": [("stab", D1, 1)],
    "test_radical_stabilizer_of_a_nonprimitive_radical_vector": [("stab", D2, 2)],
    "test_radical_stabilizer_of_a_purely_radical_lattice": [("stab", [[0, 0], [0, 0]], 1)],
    "test_lazy_orbit_representative_enumeration_is_a_stream": [("can", D1, 1, 1), ("pi", D1, 2)],
    "test_lazy_enumeration_yields_distinct_classes_only": [("can", D1, 1, 2), ("pi", D1, 1)],
    "test_the_index_check_detects_a_deficient_sublattice": [("shear", 1)],
    "test_the_congruence_image_is_the_sign_condition_and_not_the_whole_residue_group": [
        ("cong", 1, (2, 3, 4, 5, 6, 7, 8, 11)), ("cong", 2, (2, 3, 4, 5, 6, 7, 8, 11, 5, 6, 13)), ("cong", 3, (2, 3, 4, 5, 6, 7, 8, 11))],
    "test_the_sign_generator_is_the_diagonal_matrix_it_claims_to_be": [("cong", 2, (5,)), ("cong", 3, (5,))],
    "test_the_rank_one_transporter_realizes_only_the_two_signs": [("cong", 1, (5,))],
    "test_the_missing_target_sentinel_is_gap_fail_and_the_transporter_returns_a_proven_none": [("cong", 1, (5,))],
    "test_the_congruence_stabilizer_generators_lift_to_unimodular_matrices": [("cong", 1, (5,)), ("cong", 2, (5,)), ("mgZ", 3)],
    "test_a_proper_image_still_carries_the_rank_two_residue": [("cong", 2, (5,))],
    "test_rank_two_mixed_markings_transport_by_an_explicit_congruence_lift": [("cong", 2, (5,)), ("pi", D2, 1)],
    "test_rank_two_refusals_are_separated_by_content_which_is_a_proved_invariant": [("pi", D2, 3)],
    "test_rank_one_radical_image_is_the_sign_and_the_symmetry_of_c": [("cong", 1, (5,)), ("pi", D1, 1)],
    "test_rank_one_refusal_is_exhaustive_over_the_two_units_and_two_symmetries": [("cong", 1, (5,)), ("pi", D1, 1)],
    "test_the_mixed_marking_residue_obstruction_is_necessary_and_never_blocks_a_transport": [("pi", D1, 2 + 4 + 5 + 2), ("pi", D2, 6)],
    "test_the_pointwise_stabilizer_maps_onto_the_congruence_pairs_and_is_infinite_at_rank_two": [
        ("cong", 1, (5,)), ("cong", 2, (5,)), ("stab", D1, 3), ("stab", D2, 1), ("mgZ", 4)],
    "test_the_orbit_set_stabilizer_consumer_refuses_to_resolve_rather_than_returning_empty": [("stab", D1, 4)],
}


def _dpt_param(fn, name):
    """Parametrized cases of tests/test_degenerate_pointed_markings.py: lattice <2>^c + 0^t from the parameters."""
    params = name.split("[", 1)[1].rstrip("]").split("-")
    c, t = int(params[0]), int(params[1])
    g = [[2 if i == j and i < c else 0 for j in range(c + t)] for i in range(c + t)]
    if fn == "test_shear_kernel_through_the_maintained_hom_evaluation_kernel":
        return [("shear", 1), ("pi", g, 0)]
    if fn == "test_generated_shear_kernel_equals_the_maintained_kernel":
        k = c * t + 3  # generators: O(C) gens, shears, GL_t gens
        return [("stab", g, 1), ("closure", min(k ** 3, 1000), k), ("shear", 1)]
    return [("stab", g, 1)]


MOC = {  # tests/test_marked_orbit_canonical.py, traced call by call; canonicalization of height h scans 2(2h+1) candidates
    "test_canonical_representative_is_the_least_member_of_its_orbit": [("can", D1, 5, 1)],
    "test_canonical_representatives_distinguish_equal_invariant_vectors": [("can", D1, 5, 2), ("pi", D1, 1)],
    "test_canonicalization_is_stable_and_idempotent": [("can", D1, 5, 3), ("pi", D1, 2)],
    "test_enumeration_of_a_finite_slice_is_exact_and_duplicate_free": [("can", D1, 5, 1 + 3 + 3), ("pi", D1, 3)],
    "test_the_stream_is_prefix_stable_and_reports_no_completion": [("can", D1, 5, 2)],
    "test_radical_transport_is_decided_for_every_content_with_an_integral_witness": [("pi", D1, 7)],
    "test_the_zero_vector_is_its_own_radical_class": [("pi", D1, 2), ("can", D1, 2, 1)],
    "test_the_unbounded_content_family_is_a_constructive_infinite_quotient": [("can", D1, 3, 2), ("pi", D1, 4), ("can", D1, 7, 3)],
    "test_radical_transport_at_radical_rank_two_needs_a_nontrivial_unimodular_block": [("pi", D2, 6)],
    "test_radical_transport_at_rank_two_in_a_permuted_basis": [("pi", [[0, 0, 0], [0, 0, 0], [0, 0, 2]], 1)],
    "test_radical_transport_between_distinct_presentations": [("pi", D1, 3), ("pi", DP, 0)],
    "test_radical_transport_between_presentations_at_rank_two": [("pi", D2, 4), ("pi", [[0, 0, 0], [0, 0, 0], [0, 0, 2]], 0)],
    "test_radical_transport_between_presentations_distinguishes_inequivalence_from_undecided": [
        ("pi", D1, 2), ("pi", [[4, 0], [0, 0]], 0), ("pi", [[6, 0], [0, 0]], 1), ("pi", [[0, 0], [0, 6]], 0)],
    "test_radical_transport_in_a_permuted_and_mixed_basis": [("pi", D1, 1), ("pi", DP, 1), ("pi", [[2, 2], [2, 2]], 1)],
    "test_the_quotient_is_unchanged_under_a_permuted_presentation": [("can", DP, 6, 1), ("pi", D1, 2), ("pi", DP, 2)],
    "test_an_undecided_equivalence_is_never_read_as_inequivalence": [("pi", D1, 0)],
    "test_canonicalizer_recognizes_an_already_canonical_vector": [("can", [[0, 1], [1, 0]], 2, 2)],
}
DMS = {  # tests/test_degenerate_markings.sage
    "test_setwise_and_pointwise_stabilizers_of_x": [("stab", D1, 2), ("closure", 4, 3)],
    "test_setwise_stabilizer_of_a_tilted_line": [("stab", D1, 1), ("closure", 4, 3)],
    "test_public_stabilizers_on_a_degenerate_lattice": [("stab", D1, 2)],
    "test_marking_inside_the_radical": [("stab", D2, 2), ("closure", 216 + 108, 6)],
    "test_marking_meeting_the_radical_and_the_quotient": [("stab", D2, 2)],
}


def _quotients(m):
    """(m-j)U for the j-th induction step in mU."""
    return [_diag_sum(*([_u()] * (m - j))) if m - j > 0 else [] for j in range(m + 1)]


def _iso(g):
    return B.t_isometry(g, g)[0] if g else 0.0


def classify(m, k):
    """DG §Markings by isotropic induction: isotropic vector classes (Eichler criterion, |A|=1) then k-1 extension steps,
    each a parabolic in mU and generators of O((m-j)U) (B's model)."""
    q = _quotients(m)
    return C_LAT + c_tqf(1) + sum(parabolic(2 * m) + (gens(q[j]) if q[j] else 0) for j in range(1, k))


def transport(m, k):
    """DG §Markings by isotropic induction: k steps, each a parabolic in mU and one isometry decision in (m-j)U."""
    q = _quotients(m)
    return C_LAT + sum(parabolic(2 * m) + _iso(q[j]) for j in range(k))


GIM = {  # tests/test_general_isotropic_markings.py, traced
    "test_extension_shortening_preserves_the_integral_flag_orbit": (C_LAT + 2 * 10 * C_MOR, "2 shortened extensions in 2U"),
    "test_glued_line_invariant_retains_the_rational_denominator_prime": (C_LAT + 1 * (C_LAT + 10 * C_MOR), "1 flag invariant, k=1"),
    "test_public_pointwise_plane_stabilizer_in_two_u_is_the_integral_skew_shear_group": (C_LAT + parabolic(4), "pointwise parabolic of a maximal plane in 2U"),
    "test_public_pointwise_maximal_stabilizer_in_three_u_has_all_three_skew_parameters": (C_LAT + parabolic(6), "pointwise parabolic of a maximal space in 3U"),
    "test_isotropic_orbit_bases_retain_the_original_lattice_presentation": (C_LAT + c_tqf(1) + 2 * parabolic(4) + gens(_u()), "isotropic_flag_bases (2,) in U+[[2,1],[1,0]]: vector classes, one extension with O(U)"),
    "test_three_u_maximal_isotropic_sublattices_have_a_transporter": (transport(3, 3), "transporter of maximal spaces in 3U, 3 steps"),
    "test_three_u_partial_flags_have_a_transporter": (transport(3, 3), "transporter of (1,3)-flags in 3U, 3 steps"),
    "test_three_u_complete_flags_have_a_transporter": (transport(3, 3), "transporter of complete flags in 3U, 3 steps"),
    "test_complete_isotropic_flags_in_three_u_split_into_two_special_orthogonal_orbits": (C_LAT + parabolic(6) + c_ss(2), "flag stabilizer + character quotient |Q|=2"),
    "test_isotropic_rank_above_witt_index_has_no_markings": (C_LAT + 1.0, "Witt-index check (isotropic_induction.py:481) stops at once"),
    "test_three_u_general_induction_classifies_maximal_isotropic_spaces": (classify(3, 3), "classes of maximal spaces in 3U"),
    "test_complete_flag_marking_retains_each_primitive_term": (C_LAT + 3 * 10 * C_MOR, "3 primitive saturations"),
    "test_partial_flag_above_witt_index_has_no_markings": (C_LAT + 1.0, "Witt-index check stops at once"),
    "test_three_u_general_induction_classifies_complete_and_partial_flags": (2 * classify(3, 3), "two classifications (1,3) and (1,2,3) in 3U"),
    "test_public_partial_two_three_flag_classes_in_three_u": (classify(3, 3), "(2,3)-flag classes in 3U"),
    "test_public_complete_flag_classes_in_three_u": (classify(3, 3), "complete flag classes in 3U"),
    "test_core_partial_flags_have_a_natural_isomorphism_across_presentations": (C_LAT + transport(3, 3), "transporter of (1,3)-flags between two presentations of 3U"),
    "test_zero_isotropic_marking_is_the_unique_empty_embedding": (C_LAT + 2 * (gens(_u()) + 5 * C_MOR), "rank 0: O(U) generators, 2 categories"),
}


def _gim(fn, name):
    if fn == "test_general_flag_invariants_are_preserved_by_hyperbolic_pair_exchange":
        k = {"plane": 1, "partial": 2, "complete": 3}[name.split("[")[1].rstrip("]")]
        return C_LAT + 2 * k * (C_LAT + 10 * C_MOR), f"2 flag invariants with {k} reductions in 3U (no group computation)"
    if fn == "test_complete_flag_transporters_respect_the_requested_special_orthogonal_group":
        return transport(3, 3) + parabolic(6) + c_ss(2), "transporter (3 steps) + stabilizer for the SO character kernel"
    raw, why = GIM[fn]
    return raw, why


def centralizer(n, A, order=1):
    """centralizers.md §Cost, imported formula: 2 + 0.5 n^2 + T_Oq(|A|, log2|image|) + ord(f) n^3 U_ELT (B.t_oq)."""
    import math
    return 2 + 0.5 * n * n + B.t_oq(A, max(1.0, math.log2(max(A, 2)))) + order * n ** 3 * U_ELT


EQM = {  # tests/test_equivariant_markings.sage: (markings as (n, |A|), equivariant ops, helper calls, GAP IsFinite over ZZ)
    "test_unipotent_nondegenerate_stabilizer": ([(4, 1)], 1, 0, 0),
    "test_unipotent_source_in_the_radical": ([(3, 2), (3, 2)], 2, 0, 0),
    "test_unipotent_transporter_is_linear": ([(3, 2)], 2, 0, 0),
    "test_noncommuting_pair_nondegenerate_source": ([(4, 1)], 1, 0, 0),
    "test_noncommuting_pair_degenerate_source": ([(4, 1)], 3, 0, 0),
    "test_noncommuting_pair_same_image_transporter": ([(4, 1)], 1, 0, 0),
    "test_glue_controls_the_sign": ([(2, 1), (2, 1)], 2, 0, 0),
    "test_degenerate_source_with_infinite_orbit_pointwise": ([(4, 1)], 1, 0, 0),
    "test_degenerate_source_transporters_through_the_radical": ([(4, 1)], 3, 0, 0),
    "test_degenerate_source_transporter_needs_the_unipotent_correction": ([(4, 1)], 3, 0, 0),
    "test_degenerate_setwise_stabilizer_moves_the_radical": ([(4, 1)], 1, 0, 0),
    "test_invisible_source_factor_has_a_projection_kernel": ([], 0, 2, 0),
    "test_invisible_source_factor_with_infinite_automorphisms": ([], 0, 2, 1),
    "test_pell_pair_transporter_between_the_copies": ([(4, 4)], 1, 0, 0),
    "test_pell_pair_pointwise_stabilizer_is_the_infinite_centralizer_of_the_complement": ([(4, 4)], 1, 0, 0),
    "test_pell_pair_setwise_stabilizer_moves_the_source": ([(4, 4)], 1, 0, 0),
    "test_degenerate_setwise_stabilizer_with_a_non_scalar_radical_action": ([(4, 1)], 1, 0, 0),
    "test_setwise_stabilizer_of_a_degenerate_span_that_is_not_invariant": ([(4, 1)], 1, 0, 0),
    "test_unit_group_of_an_endomorphism_order_with_radical": ([], 0, 1, 0),
    "test_setwise_stabilizer_of_a_non_split_degenerate_marking": ([(4, 1)], 2, 0, 0),
    "test_line_route_pointwise_stabilizer_is_infinite": ([(4, 1)], 1, 0, 1),
    "test_source_kernel_of_an_infinite_group_with_a_finite_orbit": ([], 0, 2, 0),
    "test_source_kernel_for_the_orthogonal_source_group_is_the_invisible_kernel": ([], 0, 1, 0),
    "test_orthogonal_source_kernel_with_an_infinite_orbit": ([], 0, 1, 1),
    "test_reduced_transporter_finds_a_translate_of_the_marking": ([(4, 1)], 1, 0, 0),
}


def _eqm(markings, ops, helpers, finite_checks):
    raw = sum(centralizer(n, A) for n, A in markings) + ops * C_OP + helpers * C_OP + finite_checks * C_GAPZ
    return raw, f"{len(markings)} markings (n,|A|)={markings}, {ops} stabilizer/transporter steps, {helpers} helpers, {finite_checks} IsFinite"


# A2: its 6 roots are the vectors of norm 2, and |O(A2)| = 12.
A2_ROOTS, A2_LOG2_AUT = 6, 3.585


def _nz(fn):
    """tests/test_finite_normalizers.sage, traced per branch of normalizer_generators (finite_normalizers.py:1317-1343)."""
    iso2 = B.t_isometry(_u(), _u())[0]
    match fn:
        case "test_klein_four_on_two_hyperbolic_planes":
            return centralizer(4, 1) + 6 * 2 * iso2 + 3 * 32 * 2 * C_MM, "finite branch: |H|=4, |Aut(H)|=6 twisted transporters on 2 rank-2 commutant pieces; closures of 16 and 32"
        case "test_klein_four_on_non_isometric_summands":
            return centralizer(4, 4) + 7 * 2 * iso2 + 16 * 2 * C_MM, "finite branch: 6+1 twisted transporters, closure 16"
        case "test_glued_eigenlattices_of_the_hyperbolic_plane":
            return centralizer(2, 1) + 7 * iso2 + 8 * 2 * C_MM, "finite branch on U: 6+1 twisted transporters, closures of 4"
        case "test_order_six_rotation_on_an_indefinite_sum_of_hexagonal_planes":
            return centralizer(4, 9, 6) + 2 * 2 * iso2 + 108 * 2 * C_MM, "finite branch: rotation of order 6, |Aut(Z/6)|=2, closures 36 and 72"
        case "test_pell_isometry_normalizer_contains_the_inversion" | "test_infinite_dihedral_normalizer_is_the_orthogonal_group_of_the_pell_plane":
            return 2 * (centralizer(2, 2) + C_OP), "infinite cyclic branch: commutant Z[f] unit (real quadratic), inversion transporter, twice"
        case "test_unipotent_eichler_normalizer_contains_the_inversion":
            return centralizer(4, 1) + C_OP, "infinite cyclic branch, unipotent f on 2U"
        case "test_free_abelian_normalizer_of_two_pell_units":
            return centralizer(4, 6) + 2 * C_OP + 4 * 4 * C_MM, "abelian branch: torus algebra rank 2, decomposition normalizer"
        case "test_virtually_abelian_normalizer_with_a_nonabelian_centralizer_of_the_kernel":
            return centralizer(4, 6) + 3 * C_OP, "virtually abelian branch, rank 4"
        case "test_congruence_normalizer_of_the_summand_swap_modulo_three":
            return C_LAT + centralizer(4, 1) + 1152 * 6 * C_MM, "congruence branch on 2U: closure of the image mod 3 bounded by |O+_4(3)| = 1152, 6 generators"
        case "test_torus_of_a_pell_unit_is_the_special_orthogonal_group" | "test_torus_of_two_independent_pell_units_is_two_dimensional":
            return 2 * C_OP, "torus algebra: eigenspaces over QQbar of a 2x2 or 4x4 integer matrix"
        case "test_decomposition_normalizer_of_two_rotation_tori_on_a2_plus_a2_negative":
            return 2 * B.t_definite(2, A2_ROOTS, A2_LOG2_AUT) + 144 * 4 * C_MM, "decomposition normalizer, O(A2)xO(A2(-1)), closure 144"
        case "test_decomposition_normalizer_of_two_rotation_tori_on_a2_plus_a2_exchanges_the_planes":
            return 2 * B.t_definite(2, A2_ROOTS, A2_LOG2_AUT) + 288 * 5 * C_MM, "decomposition normalizer with swap, closure 288"
        case "test_decomposition_normalizer_of_the_split_torus_on_two_hyperbolic_planes":
            return 2 * B.t_binary({"A": 1}) + 32 * 4 * C_MM, "decomposition normalizer on 2U, closure 32"
        case "test_decomposition_normalizer_with_a_rank_one_piece":
            return B.t_binary({"A": 1}) + 8 * 3 * C_MM, "decomposition normalizer <2>+U, closure 8"
        case "test_decomposition_normalizer_of_glued_scaled_hyperbolic_planes":
            return 2 * B.t_binary({"A": 4}) + c_tqf(16) + 8 * 3 * C_MM, "glued decomposition, glue in A of order 16, closure 8"
        case "test_unipotent_abelian_normalizer_of_two_eichler_transvections":
            return centralizer(6, 1) + 4 * C_OP + 9 * 4 * C_MM, "abelian unipotent branch on 3U, supplied normalizer of 4 generators"
        case "test_mixed_abelian_normalizer_of_eichler_transvections_and_a_pell_unit":
            return centralizer(8, 2) + 6 * C_OP + 9 * 6 * C_MM, "abelian mixed branch on 3U+<1>+<-2>, 6 supplied generators"
        case "test_congruence_kernel_normalizer_of_a_unipotent_kernel_extended_by_the_plane_swap":
            return centralizer(6, 1) + 5 * C_OP + 2 * 576 * 5 * C_MM, "congruence-kernel branch on 3U: closure mod 3 of the supplied group, at most 2^3 * 9 * 8 = 576 elements, twice"
    raise KeyError(fn)


def _missing_nodes():
    root = _Path(__file__).resolve().parents[2]
    return [l for l in (root / "ci/test_budget_models/missing_F.txt").read_text().splitlines() if l]


def shell(pairs, witnesses):
    """DG §Pointwise shell decision on U+U. Gamma generators once (B's model for 2U plus one parabolic). One decision
    evaluates two invariants of 24 element operations each (4 coordinates, 16 pairing products, 4 gcd steps);
    one witness is at most 4 morphism compositions and a check on 4 basis vectors (4 x 16 element operations)."""
    return C_LAT + gens(G_2U) + parabolic(4) + pairs * 2 * 24 * U_ELT + witnesses * (4 * C_MOR + 64 * U_ELT)


SHELL = {  # tests/test_pointwise_shell_decision.py: (decisions, witnesses); box = 129 isotropic vectors of [-2,2]^4,
    # same-class pairs <= 1150 (pairs with equal content and b(e1, v), computed from the box); sample = 17 vectors
    "test_the_cone_is_decided_exhaustively_on_a_box": (129 * 128 // 2 + 129, 1150),
    "test_the_decision_is_an_equivalence_relation_and_agrees_with_witnesses": (17 * 17 * 3, 17 * 17),
    "test_the_invariant_is_never_finer_than_the_classification": (17 * 17 * 2, 0),
    "test_the_signed_pairing_defect_is_fixed": (17 * 2, 17),
    "test_the_cone_hypotheses_are_checked_not_assumed": (129, 0),
    "test_the_classification_covers_the_whole_cone_and_only_the_cone": (129, 0),
    "test_retained_generators_preserve_content_and_class_on_the_imprimitive_cone": (129 * 8, 0),
    "test_content_splits_the_cone_into_scaled_primitive_orbits": (17 * 6, 17 * 3),
    "test_every_witness_is_a_genuine_gamma_element_of_the_retained_generators_group": (17 * 17, 17 * 17),
    "test_the_class_key_is_invariant_under_the_retained_gamma_generators": (17 * 8, 0),
    "test_opposite_vectors_are_never_equivalent_unless_the_class_says_so": (17 * 2, 17 * 2),
    "test_the_family_v_m_is_separated_by_m_and_reaches_across_by_symmetry": (16, 8),
    "test_e1_and_minus_e1_are_different_orbits": (2, 2),
    "test_the_unipotent_part_moves_a_inside_the_shell": (8, 8),
    "test_the_gamma_set_reaches_the_imprimitive_cone": (40, 40),
    "test_negative_m_classes_behave_like_positive_ones": (16, 8),
    "test_the_proved_generating_set_decomposes_exactly_and_generates_gamma": (16, 16),
    "test_an_isometry_outside_gamma_does_not_decompose": (2, 2),
    "test_the_retained_action_constructor_carries_the_pointwise_shell": (4, 4),
    "test_mixed_symmetry_unipotent_round_trips_across_all_four_symmetries": (16, 16),
    "test_the_reported_mixed_composition_case_is_decided": (4, 4),
    "test_decompose_is_total_exactly_on_the_isometries_fixing_e1": (16, 16),
}
PSS = {  # tests/test_pointwise_stabilizer_shell.py
    "test_the_family_has_the_stated_invariant_and_is_primitive_isotropic": (5, 0),
    "test_pointwise_stabilizer_generators_really_fix_e1": (0, 0),
    "test_line_stabilizer_is_strictly_larger_than_the_pointwise_one": (0, 0),
    "test_the_invariant_proves_distinct_gamma_orbits_for_distinct_m": (6 * 3, 0),
    "test_a_nonidentity_gamma_transport_is_witnessed_and_certified": (4, 4),
    "test_the_certificate_checker_refuses_a_non_gamma_element": (2, 2),
    "test_equal_invariant_is_not_read_as_equivalence": (3, 0),
}


def estimates_missing():
    e = {}
    for node in _missing_nodes():
        path, name = node.split("::", 1)
        fn = name.split("[")[0]
        f = path.rsplit("/", 1)[-1]
        if path.startswith("witnesses/") and f == "test_finite_normalizers.sage":
            raw = C_LAT + gens(G_2U) + centralizer(4, 1) + 2 * 1152 * 6 * C_MM
            basis = f"{NZ}: arithmetic branch on O(2U) generators + Eichler; image closure mod 3 bounded by |O+_4(3)|=1152 with 6 generators, twice; {gbranch(G_2U)}"
        elif path.startswith("witnesses/") and f == "test_equivariant_markings.sage":
            raw, why = _eqm([(5, 2)], 2, 0, 0)
            basis = f"{EQ}: U+U+<2> Jordan transporter; {why}"
        elif f == "test_degenerate_pointed_markings.py":
            raw, why = deg_ops(DPT[fn] if fn in DPT else _dpt_param(fn, name))
            basis = f"{DM}: traced calls: {why}"
        elif f == "test_marked_orbit_canonical.py":
            raw, why = deg_ops(MOC[fn])
            basis = f"{DM}: traced calls: {why}"
        elif f == "test_degenerate_markings.sage":
            raw, why = deg_ops(DMS[fn])
            basis = f"{DM}: traced calls: {why}"
        elif f == "test_pointwise_shell_decision.py":
            pairs, wit = SHELL[fn]
            raw = shell(pairs, wit)
            basis = f"{SH}: U+U, |O(e^perp/e)|=4; {pairs} invariant decisions, {wit} witnesses"
        elif f == "test_pointwise_stabilizer_shell.py":
            pairs, wit = PSS[fn]
            raw = shell(pairs, wit) + (parabolic(4) if "line_stabilizer" in fn else 0)
            basis = f"{SH}: U+U; {pairs} decisions, {wit} witnesses"
        elif f == "test_general_isotropic_markings.py":
            raw, why = _gim(fn, name)
            basis = f"{IND}: traced: {why}"
        elif f == "test_equivariant_markings.sage":
            raw, why = _eqm(*EQM[fn])
            basis = f"{EQ}: traced: {why}"
        elif f == "test_finite_normalizers.sage":
            raw, why = _nz(fn)
            basis = f"{NZ}: traced: {why}"
        elif f == "test_general_isotropic_stabilizers.py":
            match fn:
                case "test_maximal_isotropic_stabilizer_in_three_u_has_full_parabolic_image":
                    raw, why = C_LAT + parabolic(6) + B.t_reduction(6, 3, 12), "maximal parabolic in 3U, image order in GL_6(F_3)"
                case "test_partial_and_complete_flag_stabilizers_have_their_parabolic_images":
                    n, k = (6, 3) if "three_u" in name else (8, 1)
                    raw, why = C_LAT + k * (parabolic(n) + B.t_reduction(n, 3, 12)), f"{k} flag parabolics in rank {n}, image orders in GL_{n}(F_3)"
                case "test_fixed_rank_two_embedding_has_pointwise_parabolic_image":
                    raw, why = C_LAT + parabolic(4) + B.t_reduction(4, 3, 6), "pointwise parabolic in 2U, image order 3"
                case _:
                    raw, why = transport(2, 2), "coslice transporter, 2 steps in 2U"
            basis = f"{IND}: traced: {why}"
        elif f == "test_based_isotropic_orbits.py":
            if fn == "test_bases_of_a_maximal_plane_in_u_plus_two_u_have_three_classes":
                raw, why = C_LAT + c_tqf(4) + 2 * parabolic(4) + gens(_u(2)), "basis classes of a plane in U+U(2): |A|=4, quotient U(2)"
            elif fn.startswith("test_public"):
                raw, why = classify(2, 2), "plane classes in 2U"
            else:
                raw, why = 2 * transport(2, 2), "coslice and arrow transporters, 2 steps each, in U+U(2)"
            basis = f"{IND}: traced: {why}"
        elif f == "test_common_isotropic_lift.py":
            match fn:
                case "test_rank_two_induction_in_three_u_reaches_a_maximal_isotropic_space":
                    raw, why = C_LAT + parabolic(6) + gens(_u()), "one extension step from a plane in 3U, quotient U"
                case "test_invariant_module_frame_retains_the_module_and_pairing":
                    raw, why = C_LAT + gens(_u()) + 10 * C_MOR, "isotropic stabilizer action on U"
                case _:
                    raw, why = C_LAT + 20 * C_MOR + 0.5, "rational frame lift on 2U; GF(3) group order 3"
            basis = f"{IND}: traced: {why}"
        elif f == "test_congruence_isotropic_planes.py":
            if fn.startswith("test_null_planes"):
                raw = finite_split(G_2U, 4, 72, 6, 4, 2) + 15 * (transport(2, 2) + C_MOR)
                why = "split_orbit finite preimage, F=O+_4(2) of order 72, 6 double cosets; 15 pairwise transporters"
            else:
                raw = C_LAT + gens(G_2U) + c_ss(16) + 3 * transport(2, 2)
                why = "image of order 72, 3 plane transporters"
            basis = f"DG §Orbit splitting / {IND}: {why}"
        elif f == "test_flag_invariant_consumer.py":
            plane_quot = [[-2, 0], [0, -6]]
            line_quot = _diag_sum(_u(), [[-2]], [[-6]])
            iso = B.t_isometry(plane_quot, plane_quot)[0] + B.t_isometry(line_quot, line_quot)[0]
            raw = C_LAT + 2 * 2 * (C_LAT + 10 * C_MOR) + 2 * parabolic(6) + iso
            basis = f"{IND}: 2 flag invariants (2 reductions) + flag isomorphism: parabolics and isometry decisions in the quotients {gbranch(line_quot)} and {gbranch(plane_quot)}"
        elif f == "test_isotropic_8716_pointed_stabilizer.py":
            levi = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 10, 5], [0, 0, 5, 10]]
            raw = C_LAT + gens(levi) + parabolic(6) + 96 * 6 * (6 * C_MM) + c_ss(96)
            basis = f"DG §Discriminant image: pointed null stabilizer in {gbranch(G_8716)}; Levi {gbranch(levi)} of order 96 (definite, PARI qfauto); 96 x 6 Schreier words of <= 6 products"
        elif f == "test_isotropic_integral_image.py":
            raw = C_LAT + 4 * U_ELT
            basis = "one coordinate read on U; no group computation"
        elif f == "test_isotropic_quotient_image.py":
            quot = [[0, 2], [2, 0]]
            if "proper_quotient" in fn:
                raw = C_LAT + parabolic(4) + 2 * gens(quot) + 2 * 0.5
                why = "line stabilizer and O(quotient) on the glued lattice; quotient U(2); two GF(3) orders"
            else:
                raw = C_LAT + parabolic(4) + gens(quot) + transport(2, 2)
                why = "one extension step and one flag isometry"
            basis = f"{IND}: traced: {why}"
        elif f == "test_isotropic_radical_actions.py":
            if "induction" in fn:
                raw = C_LAT + parabolic(6) + gens(_u()) + 0.5
            else:
                raw = C_LAT + 10 * C_MOR + 0.5
            basis = f"{IND}: GL_3 radical action images over F_2 (orders 168, 24, 96)"
        elif f == "test_finite_representation_transporter.py":
            raw = C_LAT + gens(G_U_M2) + 2 * parabolic(3) + c_tqf(2)
            basis = f"DG §Congruence: pointed isometry with finite-representation transporter in {gbranch(G_U_M2)}"
        elif f == "test_integral_stabilizer_generating_set.py":
            k = {"test_generating_set_keeps_sign_and_unipotent_kernel_actions": 6,
                 "test_equal_norm_block_permutations_reduce_to_generators_of_symmetric_group": 3,
                 "test_large_noncommuting_actions_retain_their_exact_integral_generating_family": 2}.get(fn)
            if k is None:
                primitive = "unipotent_powers_generate_the_primitive" in fn
                k = (6 if "powers1" in name else 3) if primitive else (4 if "powers1" in name else 2)
            n = 6 if "symmetric_group" in fn else 4
            raw = C_LAT + k * n * n * C_MOR + (C_GAPZ if n == 6 else 0)
            basis = f"{IND}: _integral_generator_set on {k} rational isometries of rank {n}: k n^2 morphism evaluations; no group computation"
        else:
            raise KeyError(node)
        e[node] = budget(raw, basis)
    return e

def estimates():
    e = _estimates_first()
    e.update(estimates_missing())
    return e


if __name__ == "__main__":
    import pathlib, re, statistics
    est = estimates()
    root = pathlib.Path(__file__).resolve().parents[2]
    pat = re.compile(r"^(witnesses/|tests/test_(acceptance_ci|allcock_fixture_identity|allcock_same_genus_nonisometric|issue15_defect_audit|configuration_pairings|dawes_2U_minus2_minus6_stabilizers|dawes_L2_boundary_witnesses|dawes_scaled_planes_discriminant_flags|maximal_split_isotropic_stabilizer|marked_subgroup_transporter|hyperbolic_odd_isometry_witness|unsupported_quotient_group|orbit_counts_basis_change|lorentzian_reduction_image_orders|enriques_line_stabilizers)\.py)")
    scope = {l for l in (root / "ci/test_budget_models/nodes_initial.txt").read_text().split() if pat.match(l)}
    scope |= set(_missing_nodes())
    assert set(est) == scope, (sorted(set(est) - scope), sorted(scope - set(est)))
    print("count", len(est))
    funcs = {}
    for k, (s, _b) in est.items():
        funcs.setdefault(k.split("[")[0], []).append(s)
    for f, v in sorted(funcs.items()):
        print(f"{min(v):8.1f} {statistics.median(v):8.1f} {max(v):8.1f}  {f}")
