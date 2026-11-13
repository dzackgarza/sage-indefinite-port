r"""A-priori time budgets (seconds) for scope G, derived by branch classification.

Scope: tests/test_commutant_groups.sage, witnesses/test_bianchi_groups.sage, tests/test_fuchsian_domain.sage,
tests/test_hermitian_ocean.sage, tests/test_quaternion_voronoi.sage,
witnesses/test_split_quaternary_orthogonal_number_field.sage, witnesses/test_ternary_orthogonal_number_field.sage.

Runner: GitHub ubuntu-latest, 4 vCPU, 16 GB. Every budget is (branch work) x 3. The work counts and the unit
costs are derived in references/complexity/hermitian-and-arithmetic-groups.md; each basis string names the section
in brackets. Rule: a pure computation on one model instance is charged once (a correct implementation computes
the Voronoi tessellation of one order once), so repeated recomputation in the code is exposed, not absorbed.
No budget uses a measured duration.
"""

from __future__ import annotations

from pathlib import Path

FACTOR = 3.0
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

V_Q = 3.6  # one QuaternionVoronoi._enumerate of a small-discriminant order over Q [Voronoi over Q: enumeration count]
OCEAN_DIAG = 0.6  # HermitianOcean(diag(1, -3)): <= 4 cell orbits x 4 sides x 4 reps x 16 candidates
OCEAN_EX = 0.1  # HermitianOcean of Example 7.1: 1 cell orbit

S_COMM = "Commutant tests: branch classification"
S_QV = "Voronoi over Q: enumeration count"
S_OCEAN = "Ocean tests: cell counts"
S_FUCHS = "Fuchsian domain: enumeration radius"
S_BIANCHI = "Swan covering for Q(sqrt -5)"
S_NF = "Voronoi over Q(sqrt 5) for (11, 19): perfect-class estimate"

# (work seconds, section, derivation)
_WORK: dict[tuple[str, str], tuple[float, str, str]] = {
    # ---- commutant groups ----
    (COMM, "test_unipotent_transporter_to_the_inverse_exists"): (0.3, S_COMM, "r=3, A/rad=Q: unipotent_transporter, one Hermite form on the coset t_r(+-N)"),
    (COMM, "test_unipotent_transporter_to_the_square_is_absent"): (0.3, S_COMM, "r=3, radical branch, empty coset after one Hermite form"),
    (COMM, "test_symplectic_quaternion_centralizer_is_the_norm_one_units"): (V_Q + 0.3, S_COMM, "QuaternionCommutant symplectic: 1 V_Q of Z<1,i,j,k> in (-1,3) + norm-one Schreier (index 2)"),
    (COMM, "test_symplectic_quaternion_transporters"): (V_Q + 0.9, S_COMM, "1 V_Q charged once + 3 x (intertwiner, solution lattice, qfminim rank 4 under cell_trace_bound) 0.3 s"),
    (COMM, "test_orthogonal_quaternion_centralizer_is_the_real_torus"): (1.0, S_COMM, "r=8 commutant 0.5 s; orthogonal adjoint: unit group of Q(sqrt3), binary representation of -1, Schreier"),
    (COMM, "test_orthogonal_quaternion_transporters"): (1.5, S_COMM, "r=8; 2 calls: kernel lattice + binary representation"),
    (COMM, "test_gaussian_hyperbolic_centralizer_contains_the_scalar_units"): (0.5, S_COMM, "r=4, central field Q(i), hyperbolic model: explicit U(1,1) generators"),
    (COMM, "test_real_quadratic_symplectic_centralizer"): (1.5, S_COMM, "r=8, central field Q(sqrt2), alternating rank 2: SL2(O_K) model + coset reps of C(Z)"),
    (COMM, "test_gaussian_hyperbolic_transporters"): (1.5, S_COMM, "3 calls x 0.5 s: hyperbolic model isometry + coset reps"),
    (COMM, "test_anisotropic_gaussian_centralizer_is_infinite"): (OCEAN_DIAG + 0.3, S_COMM, "ocean branch: HermitianOcean(diag(1,-3)) + pullback"),
    (COMM, "test_anisotropic_gaussian_transporters"): (2 * OCEAN_DIAG + 0.6, S_COMM, "2 calls: ocean isometry builds two oceans + coset reps"),
    (COMM, "test_real_quadratic_orthogonal_centralizer"): (1.5, S_COMM, "r=8, central field Q(sqrt2), symmetric rank 2: torus K[c]^1 (bnfinit deg 4) + reflection coset"),
    (COMM, "test_vector_orbit_in_the_hyperbolic_plane"): (0.3, S_COMM, "r=2, O(U) of order 4: marked transporter of a rank-1 span, 3 calls"),
    (COMM, "test_real_quadratic_orthogonal_transporters"): (2.0, S_COMM, "2 calls: conic point + relative norm equation over Q(sqrt2)"),
    (COMM, "test_product_levi_quotient_transporter"): (0.5, S_COMM, "r=4, A/rad=QxQ: levi_transporter, rank-1 Levi forms, 2 calls"),
    (COMM, "test_anisotropic_eisenstein_centralizer_is_infinite"): (V_Q + 0.5, S_COMM, "_BinaryHermitianOrder: 1 V_Q of B cap M2(O_K) (assumed small order) + norm-one Schreier"),
    (COMM, "test_anisotropic_eisenstein_transporters"): (V_Q + 1.5, S_COMM, "1 V_Q charged once + 3 calls x (reduced-norm point qfminim + coset reps) 0.5 s"),
    (COMM, "test_hermitian_quaternion_centralizer_is_infinite"): (V_Q + 2.5, S_COMM, "r=8, dim 16: M2(D) recovery 1 s + qfsolve + 1 V_Q of R + covering certificate"),
    (COMM, "test_golden_orthogonal_centralizer_rank_four"): (5.0, S_COMM, "r=16, dim 32: commutant 3 s; H+H split; ~32 Eichler transvections pulled back to 16x16; 2,3 inert: no spinor search"),
    (COMM, "test_split_orthogonal_centralizer_of_minus_one_on_U"): (0.2, S_COMM, "r=2, split M2(Q): O(U) of order 4"),
    (COMM, "test_split_orthogonal_centralizer_rank_three_is_infinite"): (3.0, S_COMM, "split M3(Q): O(I_{2,1}) through Lattices.O() (Lorentzian rank 3, 2-3 Vinberg roots) + GAP IsFinite"),
    (COMM, "test_split_symplectic_centralizer_is_SL2"): (0.2, S_COMM, "r=4, split alternating: symplectic_generators"),
    (COMM, "test_two_pell_units_transporters"): (1.2, S_COMM, "etale Q(sqrt2) x Q(sqrt3): 4 calls x bnfisnorm deg 2"),
    (COMM, "test_swapped_factor_transporter"): (0.6, S_COMM, "etale Q(sqrt5)^2 with swapped factors, 2 calls"),
    (COMM, "test_split_orthogonal_transporters"): (0.3, S_COMM, "r=2 form isometry, 2 calls"),
    (COMM, "test_split_symplectic_transporters"): (0.3, S_COMM, "r=4 alternating isometry, 2 calls"),
    (COMM, "test_definite_hermitian_quaternion_centralizer_is_infinite"): (2.5, S_COMM, "r=8, M2(D) definite: units by Fincke-Pohst (8 Lipschitz units) + covering"),
    (COMM, "test_ramified_dyadic_eichler_generators"): (2.0, S_COMM, "rank 7 over Q(sqrt2): ~20 Eichler transvections + one realizing reflection search"),
    (COMM, "test_skew_hermitian_quaternion_centralizer_is_infinite"): (V_Q + 3.5, S_COMM, "r=16, dim 16: commutant 2 s + M2(D) recovery + qfsolve + 1 V_Q + L(u), n(1), W"),
    (COMM, "test_diagonal_null_vector_beyond_ternary_subforms"): (0.5, S_COMM, "qfsolve on rank 5 over Q(sqrt2)"),
    (COMM, "test_dyadic_spinor_norm_groups_of_binary_two_adic_lattices"): (0.3, S_COMM, "Q_2 square classes (8) and two Beli spinor-norm tables"),
    (COMM, "test_dyadic_spinor_norm_group_after_maximal_norm_splitting"): (0.3, S_COMM, "one splitting + Beli table"),
    (COMM, "test_voronoi_perfect_forms_over_real_quadratic_field"): (2.0, S_COMM, "(-1,-1) over Q(sqrt2) totally definite: Herm = K, cone dim 2, <= 2 classes; 8x8 grams"),
    (COMM, "test_voronoi_sl1_of_lipschitz_order_over_sqrt2"): (2.0, S_COMM, "same cone + finite preimage of index <= 4"),
    (COMM, "test_anisotropic_skew_hermitian_quaternion_centralizer"): (8.0, S_COMM, "r=16; Clifford C over Q(sqrt6) totally definite (computed): Herm dim 2 Voronoi + maximal orders + 3 ideal generators + Schreier"),
    (COMM, "test_anisotropic_skew_hermitian_quaternion_transporters"): (11.0, S_COMM, "Clifford route charged once 8 s + 2 transporter searches"),
    (COMM, "test_maximal_orders_containing_the_lipschitz_order_over_sqrt2"): (1.5, S_COMM, "local maximal order at (sqrt2) + neighbour search"),
    (COMM, "test_witt_one_parabolic_generators_over_sqrt2"): (1.0, S_COMM, "WittOneLattice: unit group of K, translations, O(<2>), Euclidean descent"),
    (COMM, "test_witt_one_parabolic_generators_mixed_sign_rank_one"): (1.0, S_COMM, "same over Q(sqrt5)"),
    (COMM, "test_witt_one_pseudo_basis_over_sqrt10"): (1.5, S_COMM, "same over Q(sqrt10), h=2, pseudo-basis"),
    (COMM, "test_binary_orthogonal_group_of_one_minus_three_over_sqrt5"): (1.5, S_COMM, "torus of Q(sqrt5, sqrt3): bnfinit deg 4 + reflection coset"),
    (COMM, "test_factorwise_centralizer_and_transporter"): (1.5, S_COMM, "factorwise M2(Q) x Q(sqrt2): O(U) + Pell factor, centralizer + 2 transporters"),
    (COMM, "test_factorwise_centralizer_with_a_gl_type_pair"): (1.0, S_COMM, "r=6: GL pair units of Z[phi] + O(U)"),
    (COMM, "test_factorwise_transporter_with_a_gl_type_pair"): (1.5, S_COMM, "2 factorwise transporter calls"),
    # ---- commutant witnesses ----
    (WCOMM, "test_levi_transporter_to_the_inverse_on_one_summand"): (2.0, S_COMM, "r=6, A/rad = M2(Q): levi_transporter, rank-2 Levi form isometry + radical lift; 3B call stops at the determinant ratio"),
    (WCOMM, "test_nonmaximal_levi_transporter"): (3.0, S_COMM, "r=6, Levi order Eichler of level 2: [GL2(Z):Gamma0(2)] = 3 coset lifts + form isometries, 2 calls"),
    (WCOMM, "test_hermitian_quaternion_transporters"): (V_Q + 3.5, S_COMM, "r=8, M2(D): recovery + qfsolve + 1 V_Q charged once + 2 calls (reduced-norm point under cell_trace_bound)"),
    (WCOMM, "test_golden_orthogonal_transporters_rank_four"): (6.0, S_COMM, "r=16, dim 32: commutant 3 s + H+H models of both forms + coset BFS over ~32 Eichler generators, 2 calls"),
    (WCOMM, "test_definite_hermitian_quaternion_transporters"): (3.5, S_COMM, "r=8, definite M2(D): Fincke-Pohst units + 2 model isometries"),
    (WCOMM, "test_skew_hermitian_quaternion_transporters"): (V_Q + 4.5, S_COMM, "r=16 isotropic skew-Hermitian: commutant 2 s + qJ model + 1 V_Q charged once + 2 calls"),
    # ---- quaternion Voronoi over Q ----
    (QV, "test_the_order_has_one_perfect_form_class_and_integral_generators"): (V_Q + 0.2, S_QV, "1 V_Q (charged once for both calls)"),
    (QV, "test_units_in_a_box_are_words_in_the_generators"): (14.4, S_QV, "V_Q + word ball radius 4 on <=12 moves (~3.8e4 products) 0.6 s + box 11^4 1.2 s + 308 unit descents x 30 ms"),
    (QV, "test_the_cell_area_equals_the_mass_formula"): (V_Q + 0.2, S_QV, "1 V_Q + float cell areas"),
    (QV, "test_the_norm_one_extraction_applies_to_the_voronoi_generators"): (V_Q + 0.3, S_QV, "1 V_Q + norm-one Schreier index 2"),
    (QV, "test_the_full_orthogonal_group_generators_of_the_anisotropic_lattice"): (V_Q + 1.5, S_QV, "1 V_Q + normalizer box 7^4 = 2401 x 0.4 ms + ~16 isometries checked on 16 pairs"),
    (QV, "test_public_group_generators_of_the_anisotropic_lattice"): (V_Q + 3.5, S_QV, "reduction + Witt index + is_isometric_in_box radius 3 (~1 s) + same work + 16 wrapped isometries"),
    # ---- ocean ----
    (OCEAN, "test_example_discriminant_and_anisotropy"): (0.01, S_OCEAN, "discriminant + Hilbert symbol"),
    (OCEAN, "test_inputs_outside_the_source_hypotheses_are_rejected"): (0.02, S_OCEAN, "3 _check_form raises"),
    (OCEAN, "test_example_cell_structure_is_the_printed_one"): (OCEAN_EX + 0.05, S_OCEAN, "example fixture + 4 edges x transporters"),
    (OCEAN, "test_example_vertices_are_the_printed_ones"): (OCEAN_EX + 0.05, S_OCEAN, "example fixture + 4 vertices x transporters"),
    (OCEAN, "test_example_square_has_the_printed_vertex_labels"): (OCEAN_EX, S_OCEAN, "example fixture + 4 invariants"),
    (OCEAN, "test_vertex_invariants_satisfy_the_edge_relation"): (OCEAN_EX + OCEAN_DIAG, S_OCEAN, "both fixtures + value sums"),
    (OCEAN, "test_example_euler_characteristic"): (OCEAN_EX + 0.1, S_OCEAN, "example fixture + orbit sums"),
    (OCEAN, "test_generators_preserve_the_form"): (OCEAN_EX + OCEAN_DIAG + 0.1, S_OCEAN, "both fixtures + determinant Schreier (index <= 4)"),
    (OCEAN, "test_scalar_kernel_is_not_trivial"): (OCEAN_DIAG, S_OCEAN, "diagonal fixture"),
    (OCEAN, "test_diagonal_euler_characteristic_is_the_mass_formula_value"): (OCEAN_DIAG + 0.2, S_OCEAN, "diagonal fixture + orbit sums"),
    (OCEAN, "test_diagonal_form_agrees_with_the_quaternion_route"): (OCEAN_DIAG, S_OCEAN, "diagonal fixture + matrix equality"),
    (OCEAN, "test_quaternion_generators_are_words_in_the_ocean_group"): (OCEAN_DIAG + V_Q + 1.0, S_OCEAN, "diagonal fixture + 1 V_Q + ~10 best-first ocean walks"),
    (OCEAN, "test_reduction_of_an_ocean_generator"): (OCEAN_DIAG + 0.1, S_OCEAN, "diagonal fixture + one walk"),
    (OCEAN, "test_images_modulo_split_primes_agree[5-2]"): (OCEAN_DIAG + V_Q + 0.1, S_OCEAN, "fixture + 1 V_Q + BFS in GL2(F5), 480 elements"),
    (OCEAN, "test_images_modulo_split_primes_agree[13-5]"): (OCEAN_DIAG + V_Q + 1.4, S_OCEAN, "fixture + 1 V_Q + 2 BFS in GL2(F13), <= 26208 elements x ~9 generators"),
    # ---- Fuchsian, Bianchi, number-field orthogonal ----
    ("tests/test_fuchsian_domain.sage", "test_fuchsian_domain_of_the_discriminant_six_maximal_order"): (5.0, S_FUCHS, "lfun + passes r=1..5: <= ~1.1e4 qfminim vectors x 0.3 ms + ~200 norm-one reductions x 5 ms"),
    ("witnesses/test_bianchi_groups.sage", "test_swan_generators_of_gl2_over_sqrt_minus_five"): (7.0, S_BIANCHI, "covering b=1..6 (~5 s, box bisection) + 714 ideal tests + ~60 generators + GAP over F3"),
    ("witnesses/test_bianchi_groups.sage", "test_stabilizer_of_the_nonfree_lattice_o_plus_p"): (10.0, S_BIANCHI, "Swan generators 7 s + Schreier of index <= 6 on ~60 conjugated generators x Hermite keys"),
    ("witnesses/test_ternary_orthogonal_number_field.sage", "test_ternary_orthogonal_group_of_one_minus_eleven_minus_nineteen_over_sqrt5"): (220.0, S_NF, "1 Voronoi (P=20, f=15: 300 neighbours x 0.6 s) + perfecting 2.4 s + normalizer BFS ~17 s + ~3400 Schreier conjugations ~17 s"),
    ("witnesses/test_split_quaternary_orthogonal_number_field.sage", "test_split_quaternary_orthogonal_group_of_the_norm_form_of_eleven_nineteen_over_sqrt5"): (600.0, S_NF, "2 Voronoi engines (R1, R2) 364 s + 2 normalizer BFS 68 s + pair Schreier 4 s + lattice-orbit BFS ~56 s + ~2.8e4 stabilizer checks ~112 s"),
}


def _scope_nodes() -> list[str]:
    path = Path(__file__).resolve().parent / "nodes.txt"
    return [n for n in path.read_text().splitlines() if n and n.split("::")[0] in FILES]


def estimates() -> dict[str, tuple[float, str]]:
    return {
        f"{path}::{name}": (round(work * FACTOR, 2), f"[{section}] {basis}; work {work:.2f} s x{FACTOR:g}")
        for (path, name), (work, section, basis) in _WORK.items()
    }


if __name__ == "__main__":
    import statistics

    table = estimates()
    assert set(table) == set(_scope_nodes()), (set(table) ^ set(_scope_nodes()))
    for file in FILES:
        values = [s for n, (s, _) in table.items() if n.split("::")[0] == file]
        print(file, len(values), min(values), statistics.median(values), max(values))
