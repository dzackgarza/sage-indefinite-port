r"""A-priori time budgets (seconds) for scope G, derived by branch classification.

Scope: tests/test_commutant_groups.sage, witnesses/test_bianchi_groups.sage, tests/test_fuchsian_domain.sage,
tests/test_hermitian_ocean.sage, tests/test_quaternion_voronoi.sage,
witnesses/test_split_quaternary_orthogonal_number_field.sage, witnesses/test_ternary_orthogonal_number_field.sage.

Runner: GitHub ubuntu-latest, 4 vCPU, 16 GB. Every budget is the derived maximum work of the branch. The work counts and the unit
costs are derived in references/complexity/hermitian-and-arithmetic-groups.md; each basis string names the section
in brackets. Rule: a pure computation on one model instance is charged once (a correct implementation computes
the Voronoi tessellation of one order once), so repeated recomputation in the code is exposed, not absorbed.
No budget uses a measured duration.
"""

from __future__ import annotations

import json
from pathlib import Path

COMM = "tests/test_commutant_groups.sage"
OCEAN = "tests/test_hermitian_ocean.sage"
QV = "tests/test_quaternion_voronoi.sage"
WCOMM = "witnesses/test_commutant_groups.sage"
FILES = (
    COMM,
    WCOMM,
    OCEAN,
    QV,
    "tests/test_fuchsian_domain.sage",
    "witnesses/test_bianchi_groups.sage",
    "witnesses/test_split_quaternary_orthogonal_number_field.sage",
    "witnesses/test_ternary_orthogonal_number_field.sage",
)

_UNITS = json.loads((Path(__file__).resolve().parent / "unit_costs.json").read_text())["seconds_per_call"]


def u(key: str) -> float:
    return _UNITS[key]


# Derived units: operation counts (hermitian-and-arithmetic-groups.md, "Unit costs of scope G") times measured primitives.
QP = u("qq_quaternion_product")
M3 = u("nf_matrix_product_2_sqrt3")
Q4 = u("pari_qfminim_rank4_call")
GRAM_Q = 32 * M3  # _gram: 16 pairings, 2 products each
SHORT_Q = GRAM_Q + Q4
MINVEC_Q = 3 * GRAM_Q + 2 * Q4 + 20 * M3  # _minimum + second _short + ~10 values
FLIP_Q = MINVEC_Q + 2 * (SHORT_Q + 20 * M3 + MINVEC_Q)  # <= 2 bound doublings
EQUIV_Q = 2 * MINVEC_Q + 16 * (4 * QP + 3 * M3)  # m = 8 minimal vectors x 2 signs
V_Q = 2 * FLIP_Q + 6 * (FLIP_Q + MINVEC_Q + EQUIV_Q)  # 2 perfecting flips, <= 6 facets
NORM_ONE_Q = 2 * 12 * 2 * QP  # Schreier index 2 on <= 12 generators
UNITARY_G = 3 * u("nf_matrix_product_2_gaussian")  # is_unitary on 2 x 2 over Q(i)
OCEAN_DIAG = 4 * 4 * 4 * 16 * UNITARY_G  # <= 4 cell orbits x 4 sides x 4 reps x 16 candidates
OCEAN_EX = 1 * 4 * 1 * 16 * UNITARY_G
K64 = u("zz_right_kernel_64")  # commutant of rank 8 (also bounds r <= 4)
K256 = u("zz_right_kernel_256")  # commutant of rank 16
TABLE32 = 32 * 32 * u("qq_matrix_product_16")  # multiplication table, dim 32 on Q^16
TABLE16 = 16 * 16 * u("qq_matrix_product_16")
B2, B4 = u("pari_bnfinit_deg2"), u("pari_bnfinit_deg4")
H8 = u("zz_hermite_8x12")
P16 = u("qq_matrix_product_16")
QF8 = u("pari_qfminim_rank8_call")
NFQ = u("nf_quaternion_product_sqrt5")
GRAM_NF = 64 * 2 * 4 * NFQ  # 64 entries, 2 traces of 3-fold products
NEIGHBOUR_NF = 8 * GRAM_NF + 3 * QF8 + 12 * (2 * NFQ + H8)  # _step + transporter on cached minima (m = 12)
PERFECT_NF = 20  # assumption: perfect classes of (11, 19)_K
FACETS_NF = 15  # assumption
VORONOI_NF = PERFECT_NF * (FACETS_NF * NEIGHBOUR_NF + u("ppl_cone_dim6_12_rays")) + 6 * NEIGHBOUR_NF
ACT_NF = 16 * NFQ + H8  # conjugate 8 basis elements, Hermite key
NORMALIZER_NF = 20 * 170 * ACT_NF  # orbit <= 20 (assumption) x 170 generators
SCHREIER_NF = 3400 * (6 * NFQ + 2 * M3)  # conjugation matrices + checks

S_COMM = "Commutant tests: branch classification"
S_QV = "Voronoi over Q: enumeration count"
S_OCEAN = "Ocean tests: cell counts"
S_FUCHS = "Fuchsian domain: enumeration radius"
S_BIANCHI = "Swan covering for Q(sqrt -5)"
S_NF = "Voronoi over Q(sqrt 5) for (11, 19): perfect-class estimate"

# (work seconds, section, derivation)
_WORK: dict[tuple[str, str], tuple[float, str, str]] = {
    # ---- commutant groups ----
    (COMM, "test_unipotent_transporter_to_the_inverse_exists"): (K64 + 2 * H8, S_COMM, "r=3, A/rad=Q: unipotent_transporter, one Hermite form on the coset t_r(+-N)"),
    (COMM, "test_unipotent_transporter_to_the_square_is_absent"): (K64 + 2 * H8, S_COMM, "r=3, radical branch, empty coset after one Hermite form"),
    (COMM, "test_symplectic_quaternion_centralizer_is_the_norm_one_units"): (K64 + V_Q + NORM_ONE_Q, S_COMM, "QuaternionCommutant symplectic: 1 V_Q of Z<1,i,j,k> in (-1,3) + norm-one Schreier (index 2)"),
    (COMM, "test_symplectic_quaternion_transporters"): (K64 + V_Q + 3 * (K64 + Q4 + 2 * H8), S_COMM, "1 V_Q charged once + 3 x (intertwiner, solution lattice, qfminim rank 4 under cell_trace_bound)"),
    (COMM, "test_orthogonal_quaternion_centralizer_is_the_real_torus"): (K64 + B2 + 4 * H8, S_COMM, "r=8 commutant; orthogonal adjoint: unit group of Q(sqrt3), binary representation of -1, Schreier"),
    (COMM, "test_orthogonal_quaternion_transporters"): (K64 + 2 * (K64 + B2 + 4 * H8), S_COMM, "r=8; 2 calls: kernel lattice + binary representation"),
    (COMM, "test_gaussian_hyperbolic_centralizer_contains_the_scalar_units"): (K64 + 8 * P16, S_COMM, "r=4, central field Q(i), hyperbolic model: explicit U(1,1) generators"),
    (COMM, "test_real_quadratic_symplectic_centralizer"): (K64 + B2 + 3 * 8 * 8 * P16, S_COMM, "r=8, central field Q(sqrt2), alternating rank 2: SL2(O_K) model + coset reps of C(Z)"),
    (COMM, "test_gaussian_hyperbolic_transporters"): (K64 + 3 * (K64 + 8 * 8 * P16), S_COMM, "3 calls: hyperbolic model isometry + coset reps"),
    (COMM, "test_anisotropic_gaussian_centralizer_is_infinite"): (K64 + OCEAN_DIAG + 8 * P16, S_COMM, "ocean branch: HermitianOcean(diag(1,-3)) + pullback"),
    (COMM, "test_anisotropic_gaussian_transporters"): (K64 + 2 * (2 * OCEAN_DIAG + 8 * P16), S_COMM, "2 calls: ocean isometry builds two oceans + coset reps"),
    (COMM, "test_real_quadratic_orthogonal_centralizer"): (K64 + B4 + 8 * 8 * P16, S_COMM, "r=8, central field Q(sqrt2), symmetric rank 2: torus K[c]^1 (bnfinit deg 4) + reflection coset"),
    (COMM, "test_vector_orbit_in_the_hyperbolic_plane"): (3 * K64, S_COMM, "r=2, O(U) of order 4: marked transporter of a rank-1 span, 3 calls"),
    (COMM, "test_real_quadratic_orthogonal_transporters"): (K64 + 2 * (K64 + B4 + 8 * 8 * P16), S_COMM, "2 calls: conic point + relative norm equation over Q(sqrt2)"),
    (COMM, "test_product_levi_quotient_transporter"): (K64 + 2 * (K64 + 4 * H8), S_COMM, "r=4, A/rad=QxQ: levi_transporter, rank-1 Levi forms, 2 calls"),
    (COMM, "test_anisotropic_eisenstein_centralizer_is_infinite"): (K64 + V_Q + NORM_ONE_Q, S_COMM, "_BinaryHermitianOrder: 1 V_Q of B cap M2(O_K) (assumed small order) + norm-one Schreier"),
    (COMM, "test_anisotropic_eisenstein_transporters"): (K64 + V_Q + 3 * (K64 + Q4 + 8 * 8 * P16), S_COMM, "1 V_Q charged once + 3 calls x (reduced-norm point qfminim + coset reps)"),
    (COMM, "test_hermitian_quaternion_centralizer_is_infinite"): (K64 + TABLE16 + V_Q + 4 * H8, S_COMM, "r=8, dim 16: M2(D) recovery + qfsolve + 1 V_Q of R + covering certificate"),
    (COMM, "test_golden_orthogonal_centralizer_rank_four"): (K256 + TABLE32 + 32 * 2 * P16, S_COMM, "r=16, dim 32: commutant; H+H split; ~32 Eichler transvections pulled back to 16x16; 2,3 inert: no spinor search"),
    (COMM, "test_split_orthogonal_centralizer_of_minus_one_on_U"): (K64, S_COMM, "r=2, split M2(Q): O(U) of order 4"),
    (COMM, "test_split_orthogonal_centralizer_rank_three_is_infinite"): (K64 + 3 * 8 * P16 + 50 * u('libgap_call'), S_COMM, "split M3(Q): O(I_{2,1}) through Lattices.O() (Lorentzian rank 3, 2-3 Vinberg roots) + GAP IsFinite"),
    (COMM, "test_split_symplectic_centralizer_is_SL2"): (K64, S_COMM, "r=4, split alternating: symplectic_generators"),
    (COMM, "test_two_pell_units_transporters"): (K64 + 2 * B2 + 4 * (K64 + B2), S_COMM, "etale Q(sqrt2) x Q(sqrt3): 4 calls x bnfisnorm deg 2"),
    (COMM, "test_swapped_factor_transporter"): (K64 + B2 + 2 * (K64 + B2), S_COMM, "etale Q(sqrt5)^2 with swapped factors, 2 calls"),
    (COMM, "test_split_orthogonal_transporters"): (3 * K64, S_COMM, "r=2 form isometry, 2 calls"),
    (COMM, "test_split_symplectic_transporters"): (3 * K64, S_COMM, "r=4 alternating isometry, 2 calls"),
    (COMM, "test_definite_hermitian_quaternion_centralizer_is_infinite"): (K64 + TABLE16 + 8 * QF8 + 4 * H8, S_COMM, "r=8, M2(D) definite: units by Fincke-Pohst (8 Lipschitz units) + covering"),
    (COMM, "test_ramified_dyadic_eichler_generators"): (40 * 7 * 7 * M3 + 8 * QF8, S_COMM, "rank 7 over Q(sqrt2): ~20 Eichler transvections + one realizing reflection search"),
    (COMM, "test_skew_hermitian_quaternion_centralizer_is_infinite"): (K256 + TABLE16 + V_Q + 8 * H8, S_COMM, "r=16, dim 16: commutant + M2(D) recovery + qfsolve + 1 V_Q + L(u), n(1), W"),
    (COMM, "test_diagonal_null_vector_beyond_ternary_subforms"): (10 * Q4, S_COMM, "qfsolve on rank 5 over Q(sqrt2)"),
    (COMM, "test_dyadic_spinor_norm_groups_of_binary_two_adic_lattices"): (8 * 8 * M3, S_COMM, "Q_2 square classes (8) and two Beli spinor-norm tables"),
    (COMM, "test_dyadic_spinor_norm_group_after_maximal_norm_splitting"): (8 * 8 * M3, S_COMM, "one splitting + Beli table"),
    (COMM, "test_voronoi_perfect_forms_over_real_quadratic_field"): (4 * (GRAM_NF + QF8), S_COMM, "(-1,-1) over Q(sqrt2) totally definite: Herm = K, cone dim 2, <= 2 classes; 8x8 grams"),
    (COMM, "test_voronoi_sl1_of_lipschitz_order_over_sqrt2"): (4 * (GRAM_NF + QF8) + 4 * 8 * NFQ, S_COMM, "same cone + finite preimage of index <= 4"),
    (COMM, "test_anisotropic_skew_hermitian_quaternion_centralizer"): (K256 + TABLE16 + 4 * (GRAM_NF + QF8) + 3 * QF8 + 20 * ACT_NF, S_COMM, "r=16; Clifford C over Q(sqrt6) totally definite (computed): Herm dim 2 Voronoi + maximal orders + 3 ideal generators + Schreier"),
    (COMM, "test_anisotropic_skew_hermitian_quaternion_transporters"): (K256 + TABLE16 + 4 * (GRAM_NF + QF8) + 3 * QF8 + 20 * ACT_NF + 2 * (K256 + 8 * 8 * P16), S_COMM, "Clifford route charged once + 2 transporter searches"),
    (COMM, "test_maximal_orders_containing_the_lipschitz_order_over_sqrt2"): (40 * (GRAM_NF / 8 + H8), S_COMM, "local maximal order at (sqrt2) + neighbour search"),
    (COMM, "test_witt_one_parabolic_generators_over_sqrt2"): (B2 + 20 * 3 * 3 * M3, S_COMM, "WittOneLattice: unit group of K, translations, O(<2>), Euclidean descent"),
    (COMM, "test_witt_one_parabolic_generators_mixed_sign_rank_one"): (B2 + 20 * 3 * 3 * M3, S_COMM, "same over Q(sqrt5)"),
    (COMM, "test_witt_one_pseudo_basis_over_sqrt10"): (B2 + 30 * 3 * 3 * M3, S_COMM, "same over Q(sqrt10), h=2, pseudo-basis"),
    (COMM, "test_binary_orthogonal_group_of_one_minus_three_over_sqrt5"): (B4 + 20 * 2 * 2 * M3, S_COMM, "torus of Q(sqrt5, sqrt3): bnfinit deg 4 + reflection coset"),
    (COMM, "test_factorwise_centralizer_and_transporter"): (K64 + B2 + 3 * (K64 + 8 * P16), S_COMM, "factorwise M2(Q) x Q(sqrt2): O(U) + Pell factor, centralizer + 2 transporters"),
    (COMM, "test_factorwise_centralizer_with_a_gl_type_pair"): (K64 + B2 + 4 * P16, S_COMM, "r=6: GL pair units of Z[phi] + O(U)"),
    (COMM, "test_factorwise_transporter_with_a_gl_type_pair"): (K64 + 2 * (K64 + B2 + 4 * P16), S_COMM, "2 factorwise transporter calls"),
    # ---- commutant witnesses ----
    (WCOMM, "test_levi_transporter_to_the_inverse_on_one_summand"): (K64 + 2 * (K64 + 4 * H8), S_COMM, "r=6, A/rad = M2(Q): levi_transporter, rank-2 Levi form isometry + radical lift; 3B call stops at the determinant ratio"),
    (WCOMM, "test_nonmaximal_levi_transporter"): (K64 + 2 * (K64 + 3 * 4 * H8), S_COMM, "r=6, Levi order Eichler of level 2: [GL2(Z):Gamma0(2)] = 3 coset lifts + form isometries, 2 calls"),
    (WCOMM, "test_hermitian_quaternion_transporters"): (K64 + TABLE16 + V_Q + 2 * (K64 + Q4 + 4 * H8), S_COMM, "r=8, M2(D): recovery + qfsolve + 1 V_Q charged once + 2 calls (reduced-norm point under cell_trace_bound)"),
    (WCOMM, "test_golden_orthogonal_transporters_rank_four"): (K256 + TABLE32 + 2 * (K256 + 32 * 2 * P16), S_COMM, "r=16, dim 32: commutant + H+H models of both forms + coset BFS over ~32 Eichler generators, 2 calls"),
    (WCOMM, "test_definite_hermitian_quaternion_transporters"): (K64 + TABLE16 + 8 * QF8 + 2 * (K64 + 4 * H8), S_COMM, "r=8, definite M2(D): Fincke-Pohst units + 2 model isometries"),
    (WCOMM, "test_skew_hermitian_quaternion_transporters"): (K256 + TABLE16 + V_Q + 2 * (K256 + Q4 + 8 * H8), S_COMM, "r=16 isotropic skew-Hermitian: commutant + qJ model + 1 V_Q charged once + 2 calls"),
    # ---- quaternion Voronoi over Q ----
    (QV, "test_the_order_has_one_perfect_form_class_and_integral_generators"): (V_Q, S_QV, "1 V_Q (charged once for both calls)"),
    (QV, "test_units_in_a_box_are_words_in_the_generators"): (V_Q + 38000 * QP + 14641 * (4 * QP + u('qq_quaternion_reduced_norm')) + 308 * 4 * 12 * (QP + 3 * M3), S_QV, "V_Q + word ball radius 4 on <=12 moves (~3.8e4 products) + box 11^4 + 308 unit descents"),
    (QV, "test_the_cell_area_equals_the_mass_formula"): (V_Q + 6 * 200 * u('complex256_mult'), S_QV, "1 V_Q + float cell areas"),
    (QV, "test_the_norm_one_extraction_applies_to_the_voronoi_generators"): (V_Q + NORM_ONE_Q, S_QV, "1 V_Q + norm-one Schreier index 2"),
    (QV, "test_the_full_orthogonal_group_generators_of_the_anisotropic_lattice"): (V_Q + NORM_ONE_Q + 2401 * 9 * QP + 16 * 16 * 16 * u('zz_dot_12'), S_QV, "1 V_Q + normalizer box 7^4 = 2401 + ~16 isometries checked on 16 pairs"),
    (QV, "test_public_group_generators_of_the_anisotropic_lattice"): (V_Q + NORM_ONE_Q + 2401 * 9 * QP + 16 * 16 * 16 * u('zz_dot_12') + 2 * 2401 * u('zz_dot_12') + 16 * u('lattice_construction_rank_3'), S_QV, "reduction + Witt index + is_isometric_in_box radius 3 () + same work + 16 wrapped isometries"),
    # ---- ocean ----
    (OCEAN, "test_example_discriminant_and_anisotropy"): (4 * M3, S_OCEAN, "discriminant + Hilbert symbol"),
    (OCEAN, "test_inputs_outside_the_source_hypotheses_are_rejected"): (3 * M3, S_OCEAN, "3 _check_form raises"),
    (OCEAN, "test_example_cell_structure_is_the_printed_one"): (OCEAN_EX + 4 * 16 * UNITARY_G, S_OCEAN, "example fixture + 4 edges x transporters"),
    (OCEAN, "test_example_vertices_are_the_printed_ones"): (OCEAN_EX + 4 * 16 * UNITARY_G, S_OCEAN, "example fixture + 4 vertices x transporters"),
    (OCEAN, "test_example_square_has_the_printed_vertex_labels"): (OCEAN_EX, S_OCEAN, "example fixture + 4 invariants"),
    (OCEAN, "test_vertex_invariants_satisfy_the_edge_relation"): (OCEAN_EX + OCEAN_DIAG, S_OCEAN, "both fixtures + value sums"),
    (OCEAN, "test_example_euler_characteristic"): (OCEAN_EX + 8 * 16 * UNITARY_G, S_OCEAN, "example fixture + orbit sums"),
    (OCEAN, "test_generators_preserve_the_form"): (OCEAN_EX + OCEAN_DIAG + 4 * 12 * UNITARY_G, S_OCEAN, "both fixtures + determinant Schreier (index <= 4)"),
    (OCEAN, "test_scalar_kernel_is_not_trivial"): (OCEAN_DIAG, S_OCEAN, "diagonal fixture"),
    (OCEAN, "test_diagonal_euler_characteristic_is_the_mass_formula_value"): (OCEAN_DIAG + 4 * 4 * 2 * 16 * UNITARY_G, S_OCEAN, "diagonal fixture + orbit sums"),
    (OCEAN, "test_diagonal_form_agrees_with_the_quaternion_route"): (OCEAN_DIAG, S_OCEAN, "diagonal fixture + matrix equality"),
    (OCEAN, "test_quaternion_generators_are_words_in_the_ocean_group"): (OCEAN_DIAG + V_Q + NORM_ONE_Q + 10 * 40 * 4 * UNITARY_G, S_OCEAN, "diagonal fixture + 1 V_Q + ~10 best-first ocean walks"),
    (OCEAN, "test_reduction_of_an_ocean_generator"): (OCEAN_DIAG + 40 * 4 * UNITARY_G, S_OCEAN, "diagonal fixture + one walk"),
    (OCEAN, "test_images_modulo_split_primes_agree[5-2]"): (OCEAN_DIAG + V_Q + NORM_ONE_Q + 2 * 480 * 9 * u('zz_dot_12'), S_OCEAN, "fixture + 1 V_Q + BFS in GL2(F5), 480 elements"),
    (OCEAN, "test_images_modulo_split_primes_agree[13-5]"): (OCEAN_DIAG + V_Q + NORM_ONE_Q + 2 * 26208 * 9 * u('zz_dot_12'), S_OCEAN, "fixture + 1 V_Q + 2 BFS in GL2(F13), <= 26208 elements x ~9 generators"),
    # ---- Fuchsian, Bianchi, number-field orthogonal ----
    ("tests/test_fuchsian_domain.sage", "test_fuchsian_domain_of_the_discriminant_six_maximal_order"): (11000 * (u('pari_qfminim_e8_point') + 4 * QP + 16 * u('complex256_mult')) + 200 * 10 * 8 * u('complex256_mult'), S_FUCHS, "lfun + passes r=1..5: <= ~1.1e4 qfminim vectors + ~200 norm-one reductions"),
    ("witnesses/test_bianchi_groups.sage", "test_swan_generators_of_gl2_over_sqrt_minus_five"): (3000 * 270 * 8 * u('rif128_mult') + 714 * u('nf_ideal_equality_sqrt_minus5') + 60 * 4 * u('libgap_call'), S_BIANCHI, "covering b=1..6 (, box bisection) + 714 ideal tests + ~60 generators + GAP over F3"),
    ("witnesses/test_bianchi_groups.sage", "test_stabilizer_of_the_nonfree_lattice_o_plus_p"): (3000 * 270 * 8 * u('rif128_mult') + 714 * u('nf_ideal_equality_sqrt_minus5') + 6 * 60 * 2 * (H8 + 4 * M3), S_BIANCHI, "Swan generators + Schreier of index <= 6 on ~60 conjugated generators x Hermite keys"),
    ("witnesses/test_ternary_orthogonal_number_field.sage", "test_ternary_orthogonal_group_of_one_minus_eleven_minus_nineteen_over_sqrt5"): (VORONOI_NF + NORMALIZER_NF + SCHREIER_NF, S_NF, "1 Voronoi (P=20, f=15: 300 neighbours) + perfecting + normalizer BFS + ~3400 Schreier conjugations"),
    ("witnesses/test_split_quaternary_orthogonal_number_field.sage", "test_split_quaternary_orthogonal_group_of_the_norm_form_of_eleven_nineteen_over_sqrt5"): (2 * VORONOI_NF + 2 * NORMALIZER_NF + 1400 * 4 * M3 + 28000 * (H8 + 8 * M3) + 28000 * (H8 + 4 * M3), S_NF, "2 Voronoi engines (R1, R2) + 2 normalizer BFS + pair Schreier + lattice-orbit BFS + ~2.8e4 stabilizer checks"),
}


def _scope_nodes() -> list[str]:
    path = Path(__file__).resolve().parent / "nodes.txt"
    return [n for n in path.read_text().splitlines() if n and n.split("::")[0] in FILES]


def estimates() -> dict[str, tuple[float, str]]:
    return {
        f"{path}::{name}": (float(f"{work:.3g}"), f"[{section}] {basis}")
        for (path, name), (work, section, basis) in _WORK.items()
    }


if __name__ == "__main__":
    import statistics

    table = estimates()
    assert set(table) == set(_scope_nodes()), (set(table) ^ set(_scope_nodes()))
    for file in FILES:
        values = [s for n, (s, _) in table.items() if n.split("::")[0] == file]
        print(file, len(values), min(values), statistics.median(values), max(values))
