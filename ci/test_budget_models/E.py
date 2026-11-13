r"""A-priori time budgets (seconds) for scope E, derived by branch classification.

Runner: ubuntu-latest, 4 vCPU, 16 GB, Sage 10.x. Unit costs used throughout:
  U1  Sage exact QQ/ZZ small-matrix op (<= 6x6, FLINT): 10 us.
  NF  quadratic number-field element op: 2 us; cubic: 5 us. 3x3 product over a quadratic
      field: 0.1 ms; 5x5: 0.4 ms.
  U2  Sage CliffordAlgebra monomial product: 5 us.
  U3  GAP Size() of a classical matrix group over GF(3): 2 s at degree 6, 0.5 s at degree 5.
  U4  PARI qfisom rank 18 with 18 minimal vectors: 2 s; rank 6: 0.5 s; E8 qfisom/qfauto: 0.05-0.1 s.
  U5  Julia startup + Oscar load where the Nikulin embedding owner calls Oscar: 60 s fixed.
Every budget is the derived maximum work of the branch, with no safety factor and no minimum.
Entries marked "(unit-corrected)" state counts derived at earlier unit costs; their seconds
are rescaled to the units above. Session import cost is not charged to a test. Budgets are independent of
measured durations. Theory: references/complexity/clifford-and-hermitian-lattices.md;
each basis string names its section in brackets.
"""

from __future__ import annotations

from pathlib import Path

SRC = "src/sage_indefinite_port/indefinite/"


def _clifford_even_rank(name: str) -> tuple[float, str]:
    """clifford_even_rank.py:19-98. Branch: EvenCliffordCentre builds CliffordAlgebra(form) (l.25) of
    dim 2^n; the centre is QQ(z) with z^2 = disc (l.45-56). Split branch (disc a square: 3U, diag(2,2,-6,-6))
    computes idempotents (1 +- z/sqrt)/2 and component_dimension by multiplying the idempotent against the
    2^(n-1) even keys (l.66-77). Nonsplit branch (2U+<1,-2>, disc 2) only builds QQ(sqrt 2), no idempotents."""
    if name == "test_three_hyperbolic_planes_split_into_two_degree_four_components":
        # split, n=6: 2 idempotents x 32 even keys x products with 32-term idempotent -> ~2*32*32*32 term pairs
        # (U2) = 65k*20us = 1.3 s, plus volume centrality on 32 keys ~0.1 s.
        return 2.0, "split branch (disc 64 square), n=6, C0 dim 32: ~1.4 s products"
    if name == "test_the_anisotropic_rank_four_form_splits_into_two_quaternion_components":
        return 0.0175, "split branch (disc 144 square), n=4, C0 dim 8: 2*8*8*8=1024 monomial products (U2) 0.02 s + algebra 0.05 s (unit-corrected)"
    if name == "test_a_nonsplit_rank_six_form_has_a_quadratic_discriminant_field":
        return 0.02, "nonsplit branch (disc 2 nonsquare): volume square + QuadraticField(2) nfinit 0.01 s"
    if name == "test_the_nonsplit_rank_six_form_has_an_involution_of_the_second_kind":
        return 0.01, "nonsplit branch: reversal sign (-1)^15 on volume, dimension count: 0.01 s"
    return 0.01, "n=4 reversal sign on volume only (l.83-95): 0.01 s"


def _clifford_o3u(name: str) -> tuple[float, str]:
    """clifford_split_rank_six.py: elementary_generators 12 SL4 elementaries (l.95-101) -> second_compound 6x6
    (l.65-73); hyperbolic_coset_representatives 3 reflections (l.108-135); spinor characters (l.137-147).
    The mod-3 branch builds libgap groups SO^+ kernel (|.|=6065280) and O_6^+(3) (|.|=24261120)."""
    if name == "test_image_modulo_three_has_the_hand_derived_orders":
        return 4.0, "2 GAP Size() on O_6^+(3) images, natural module GF(3)^6 (728 nonzero vectors): 2 s each"
    if name == "test_elementary_images_are_isometries_in_the_spinor_kernel":
        return 0.012, "12 generators x 36 pairs b(f x,f y) + spinor chars: ~1100 U1 = 0.055 s (unit-corrected)"
    return 0.002, "single 6x6 compound / 3 reflection characters: ~200 U1 = 0.01 s (unit-corrected)"


def _clifford_rank_five(name: str) -> tuple[float, str]:
    """clifford_rank_five.py:43-103. CliffordAlgebra of rank 5 (dim 32); scalar = volume square (l.69-75);
    pfaffian gram on 5 pfaffian basis elements (l.77-88); even dims count keys (l.90-101). Same branch for all
    cases; negated case builds two models."""
    if name == "test_negated_forms_are_similar_with_scalars_differing_by_minus_one":
        return 0.25, "two rank-5 models, 25 pfaffian products in dim-32 algebra each (U2) (unit-corrected)"
    return 0.025, "one rank-5 model, <=25 products in dim-32 algebra (U2) + construction: 0.1 s (unit-corrected)"


def _clifford_split_six(name: str) -> tuple[float, str]:
    """clifford_split_rank_six.py:34-93. Fock generators 8x8 (l.34-50); even_part_dimension (l.52-63);
    second_compound 6x6 minors (l.65-73); isometry_of_special_linear (l.90-93)."""
    if name == "test_the_fock_generators_satisfy_the_clifford_relations_of_3u":
        return 0.0008, "36 anticommutators = 72 products of 8x8 QQ matrices (U1): 0.004 s (unit-corrected)"
    if name == "test_the_compound_is_multiplicative_and_preserves_3u_integrally":
        return 0.0016, "3 compounds + 36 pairings on 3U: ~150 U1 = 0.0075 s (unit-corrected)"
    return 0.002, "single compound (36 minors) / dimension count: 0.002 s"


def _clifford_uu2(name: str) -> tuple[float, str]:
    """clifford_rank_five.py:105-154 + clifford_split_rank_six.hyperbolic_coset_representatives. Sp4 generators
    (l.120-140) act on the 5-dim pfaffian complement (l.111-118). Mod-3 branch: libgap Size() of O_5(3)
    images (|.| = 25920, 103680)."""
    if name == "test_image_modulo_three_has_the_hand_derived_orders":
        return 1.0, "2 GAP Size() on O_5(3) images, GF(3)^5 (242 nonzero vectors): 0.5 s each"
    if name == "test_symplectic_images_are_isometries_in_the_spinor_kernel":
        return 0.002, "Sp4 generators x 25 pairings + characters: ~200 U1 = 0.01 s (unit-corrected)"
    return 0.01, "single pfaffian action / characters on rank 5: 0.01 s"


def _binary_unitary(name: str) -> tuple[float, str]:
    """binary_hyperbolic_unitary.py: Q(sqrt -d), d in 1,2,3 are norm-Euclidean (class number 1), so word()
    (l.54-80) terminates by Euclidean division, length O(log height). Public centralizer test goes to
    sage_extensions HermitianLattice.trace_lattice().O().centralizer: trace lattice U+U-type rank 4, signature
    (2,2), Witt index 2 indefinite group; centralizer of the order-3/4 unit isometry = U(1,1) over Z[zeta],
    computed by the indefinite generator route (fixed lattice of zeta is 0, so centralizer equals the unitary
    group): budget allows a full rank-4 indefinite O(L) generator computation."""
    if name == "test_public_centralizer_on_the_trace_of_the_hyperbolic_plane":
        return 40.0, "2 fields (Q(zeta3), Q(zeta4)), class no. 1; rank-4 sig (2,2) centralizer via indefinite route ~20 s each"
    if name == "test_words_recover_long_products":
        return 0.08, "3 fields x 8 words of 18 letters: Euclidean reduction, ~100 steps of 2x2 NF ops each (unit-corrected)"
    return 0.00024, "a few 2x2 NF products / one Euclidean word: ~100 NF ops = 0.005 s (unit-corrected)"


def _eisenstein(name: str) -> tuple[float, str]:
    """eisenstein_hyperbolic_unitary.py: Z[omega] norm-Euclidean, class number 1; word() (l.171-216) reduces by
    nearest-integer steps (l.109-113). trace_generators (l.225-228) build trace lattice of rank 2m+4 <= 8."""
    if name == "test_trace_generators_commute_with_omega":
        return 0.04, "m=1,2: <=~10 generators x rank<=8 points, NF trace lattice (unit-corrected)"
    return 0.004, "m=1,2 Euclidean words on <=4x4 matrices over Q(omega): ~2000 NF ops = 0.1 s (unit-corrected)"


def _quaternion(name: str) -> tuple[float, str]:
    """quaternion_lattices.py: (-1,3)_QQ, ramified at 2 and 3. normalizer_classes_in_box (l.108-117) loops 7^4=2401
    elements, normalizes_order (l.97-100) 4 conjugations, compares against <=12 classes (l.102-106)."""
    if name == "test_the_normalizer_quotient_has_exactly_twelve_classes":
        return 3.0, "box radius 3: 2401 x (4 conj + <=12 compares) ~ 4e4 quaternion ops x 50 us = 2 s"
    return 0.002, "4x4 multiplication / conjugation matrices, <=16 pairings: ~200 U1 = 0.01 s (unit-corrected)"


def _genus_nf(name: str) -> tuple[float, str]:
    """genus_number_field.py: lattice_isometry (l.182-220). Definite branch (sign 1): trace forms over Q(sqrt5)
    (l.128-135) give rank-2n positive definite Z-forms, decided by PARI qfisom (l.210). Indefinite binary
    branch (sign 0): trace-form commutant route (l.153-180). Neighbour test: WittOneLattice + Cusps.frame
    (l.271-297), neighbour_vector (l.353-), then a definite rank-3 isometry test."""
    if name == "test_odd_unimodular_rank_nine_genus_has_two_inequivalent_lattices_over_sqrt5":
        return 4.0, "definite branch: 2 qfisom on rank-18 trace forms, 18 vectors of norm 1, |Aut(I_9 x O_K)| = 2^9 9!: 2 s each"
    if name == "test_sqrt5_neighbour_of_three_squares_is_three_squares":
        return 0.5, "Witt index 1, h(Q(sqrt5))=1: one frame + qfisom rank 6 (I_3 trace form): 0.5 s"
    return 0.01, "indefinite binary branch: commutant linear algebra on 4x4, 2 calls: 0.01 s"


def _voronoi(name: str) -> tuple[float, str]:
    """semisimple_voronoi.py: M_2(Z) with splitting None (l.126-131), field QQ, herm dim 3; one GL_2(Z) class of
    perfect forms (A2), so perfect_forms (l.345-362) explores one form and its facet neighbours (3 facets,
    Polyhedron l.290); stabilizer/transporter candidates (l.304-341) <= |min|*(2|min|)^(k-1) with |min| ~ 12."""
    return 6.0, "QQ branch, 1 perfect class, 3 facets, ~300 candidates per stabilizer, qfminim small: ~5 s"


def _witt_one(name: str) -> tuple[float, str]:
    """witt_one_hilbert.py: D in {2,3,5}; Q(sqrt D) class number 1 and norm-Euclidean, so _reduce_pair (l.140-155)
    terminates. Public tests on U+<2,-4> (det 8 = disc Q(sqrt2)): recursive.py:2280 selects recognized_model
    (l.185-220): null vectors in 7^4 box, inner change search <= 13^4 = 28561 2x2 checks for the first
    hyperbolic pair; h=1 so cusp_links (witt_one_fixed_norm_classes) adds nothing beyond one cusp."""
    if name == "test_public_group_generators_of_the_issue_lattice":
        return 6.0, "recognized Hilbert branch D=2, h=1: box 2401 + <=28561 2x2 checks (~2 s) + pairings"
    if name == "test_one_orbit_of_primitive_isotropic_lines_through_the_public_group":
        return 12.0, "Witt index 1, h=1 -> one cusp; recognition (~2 s) + isotropic orbit dedupe"
    if name == "test_vector_classes_at_witt_index_one_are_decided_exactly":
        return 12.0, "SubgroupVectorOrbits (pairing_slices.py:372-385): 3 vector_equivalence_witness on Witt-1 recognized group"
    if name == "test_null_vectors_reach_e_and_products_factor":
        return 0.064, "2 fields x 10 trials x (10 4x4 products + Euclidean reduction): ~6000 U1 = 0.3 s + bnfinit (unit-corrected)"
    if name == "test_generators_are_integral_isometries":
        return 5.0, "3 fields: unit_group fundamental unit (PARI bnfinit deg 2) + ~7 4x4 checks"
    return 8e-05, "4x4 Gram from 16 NF pairings: 0.002 s (unit-corrected)"


def _voronoi_cubic(name: str) -> tuple[float, str]:
    """semisimple_voronoi.py split branch (l.131-150), division algebra (L/Q, sigma, 2), herm dim 6 over cubic L,
    order rank 9. Per neighbour: ~15 probes x (9x9 Gram over L, ~2500 cubic-field ops = 0.13 s, + qfminim rank 9 0.01 s)
    = 2 s. P ~ 7800 by Prasad covolume scaling from BCNS Sec. 6.4 (431 classes for ramification {2,3}); F <= 30 assumed."""
    if name == "test_split_model_is_a_homomorphism":
        return 0.01, "one 9x9 QQ product + two _split images over L: ~200 ops = 0.01 s"
    return 82000.0, "division branch, ramified {2,7}: P*F*(15 probes x (9x9 Gram over L, 2500 cubic-field ops at 5 us = 12 ms, + qfminim 10 ms)) = 7800*30*0.35 s = 8.2e4 s (P scaled from BCNS 6.4 431 forms by covolume ratio 18; F assumed)"


def _isotropic_reduction(name: str) -> tuple[float, str]:
    """sage_extensions isotropic_reduction on U+E8(-1): unimodular, A_L = 0, |O(q_L)| = 1, so the gluing is trivial
    and the Levi image is all of O(E8) (W(E8), order 696729600, ~2-3 generators from PARI qfauto). Twisted
    specimen U(2)+E8(-2): A_L = (Z/2)^10, gluing check of a reflection in O(q) of order |O_10^+(2)|."""
    if name == "test_isotropic_reduction_recovers_the_e8_leaf":
        return 0.1, "unimodular, one qfisom of E8 (240 roots): 0.1 s"
    if name == "test_levi_generators_lift_to_isometries_fixing_the_line":
        return 0.5, "A_L=0: qfauto E8 0.05 s + one split lift per generator (<=3): 0.5 s"
    return 10.0, "A_L=(Z/2)^10 twisted gluing: one element membership in image in O(q_L)"


def _spacelike(name: str) -> tuple[float, str]:
    """marked_lattices.py:583-600: U+<-2>, positive index 1, sign*square < 0 -> spacelike branch: Nikulin
    primitive_embedding_class_representatives of <-2> (pinned OSCAR embeddings.jl). A_L = Z/2, |O(q_L)| = 1,
    complement genus rank 2 det +-4 sig (1,1), class number small; may load Julia/Oscar (U5)."""
    return 65.0, "spacelike Nikulin branch, A_L=Z/2, |O(q)|=1, tiny genus: ~5 s + Oscar load 60 s"


SECTIONS = {'_clifford_even_rank': 'Even Clifford centre', '_clifford_o3u': 'Rank-five model and transfers', '_clifford_rank_five': 'Rank-five model and transfers', '_clifford_split_six': 'Rank-five model and transfers', '_clifford_uu2': 'Rank-five model and transfers', '_binary_unitary': 'Hermitian lattices', '_eisenstein': 'Hermitian lattices', '_quaternion': 'Quaternion norm lattices', '_genus_nf': 'Lattice isometry over a number field', '_voronoi': 'Semisimple Voronoi', '_witt_one': 'Witt-index-one Hilbert route', '_voronoi_cubic':'Semisimple Voronoi over a cubic division algebra', '_isotropic_reduction': 'isometry-testing.md (owner B)', '_spacelike': 'Spacelike vector classes'}

MODELS = {
    "test_clifford_even_rank.sage": _clifford_even_rank,
    "test_clifford_o3u_transfer.sage": _clifford_o3u,
    "test_clifford_rank_five.sage": _clifford_rank_five,
    "test_clifford_split_rank_six.sage": _clifford_split_six,
    "test_clifford_u_u_two_transfer.sage": _clifford_uu2,
    "test_binary_hyperbolic_unitary.sage": _binary_unitary,
    "test_eisenstein_hyperbolic_unitary.sage": _eisenstein,
    "test_quaternion_lattices.sage": _quaternion,
    "test_genus_number_field.sage": _genus_nf,
    "test_semisimple_voronoi.sage": _voronoi,
    "test_witt_one_hilbert.sage": _witt_one,
    "test_isotropic_reduction.sage": _isotropic_reduction,
    "witnesses/test_semisimple_voronoi_cubic_division.sage": _voronoi_cubic,
    "test_spacelike_vector_classes_index_one.py": _spacelike,
}


# ---- Hermitian and arithmetic-group nodes (theory: references/complexity/hermitian-and-arithmetic-groups.md) ----
# (file, function) -> (seconds, section, derived work)
_T = "tests/"
_W = "witnesses/"

ACT = _T + "test_act_eisenstein.sage"
COMM = _T + "test_commutant_groups.sage"
LMB = _T + "test_li_margulis_box.sage"
LMC = _T + "test_li_margulis_constants.sage"
OCEAN = _T + "test_hermitian_ocean.sage"

TABLE: dict[tuple[str, str], tuple[float, str, str]] = {
    # ACT Eisenstein lattice E^{4,1} over Z[omega] (Allcock-Carlson-Toledo Sec. 7); 5x5 products 0.4 ms.
    (ACT, "test_roots_and_generators"): (0.0112, "Explicit unitary groups", "~15 isometry checks of 5x5 (2 products each) = 0.18 s (unit-corrected)"),
    (ACT, "test_affine_e6_braid_relations"): (0.0268, "Explicit unitary groups", "21 pairs x <=6 products of 5x5 = 126 x 6 ms = 0.4 s (unit-corrected)"),
    (ACT, "test_source_translations"): (0.00268, "Explicit unitary groups", "~10 5x5 products/inverses = 0.04 s (unit-corrected)"),
    (ACT, "test_translation_composition_law"): (0.00134, "Explicit unitary groups", "3 translations, 1 product = 0.02 s (unit-corrected)"),
    (ACT, "test_null_reduction_from_source_shape"): (0.0067, "Explicit unitary groups", "Euclidean descent on height |v_3|, <=10 steps x 5x5 = 0.1 s (unit-corrected)"),
    (ACT, "test_reduction_of_a_group_image_of_rho"): (0.0067, "Explicit unitary groups", "word of length 5 then descent <=10 steps = 0.1 s (unit-corrected)"),
    (ACT, "test_orthogonal_height_one_root_case"): (0.00402, "Explicit unitary groups", "height-1 root branch: one root step = 0.05 s (unit-corrected)"),
    (ACT, "test_unit_multiples_of_rho_reach_rho"): (0.00402, "Explicit unitary groups", "5 units, each a unit-diagonal step = 0.05 s (unit-corrected)"),
    (ACT, "test_eisenstein_primitivity_is_not_integer_primitivity"): (0.00067, "Explicit unitary groups", "one gcd of 5 entries in Z[omega] (unit-corrected)"),
    (ACT, "test_trace_lattice_rank_ten_signature"): (0.00268, "Explicit unitary groups", "rank-10 trace Gram + signature = 0.04 s (unit-corrected)"),
    (ACT, "test_spinor_character_separates_aut_plus"): (0.0335, "Explicit unitary groups", "8 spinor characters, each <=10 reflections on rank 10 = 0.5 s (unit-corrected)"),
    (ACT, "test_spinor_character_is_multiplicative"): (0.335, "Explicit unitary groups", "25 pairs x 3 spinor characters (rank-10 Cartan-Dieudonne, ~0.06 s) = 5 s (unit-corrected)"),
    (ACT, "test_generators_are_trace_isometries_commuting_with_omega"): (0.0067, "Explicit unitary groups", "8 generators x 10 points x 2 maps = 0.1 s (unit-corrected)"),
    (ACT, "test_gamma_theta_contains_triflections_and_not_hexflections"): (0.0067, "Explicit unitary groups", "7 roots x (triflection + R^4) = 0.1 s (unit-corrected)"),
    (ACT, "test_null_transporter_carries_source_vectors"): (0.0067, "Explicit unitary groups", "two descents + inverse = 0.1 s (unit-corrected)"),
    (ACT, "test_short_roots_reduce_to_r7"): (0.0134, "Explicit unitary groups", "7 root descents <=5 steps = 0.2 s (unit-corrected)"),
    (ACT, "test_root_transporter_on_a_group_image"): (0.0067, "Explicit unitary groups", "two root descents = 0.1 s (unit-corrected)"),
    (ACT, "test_nonroot_is_rejected"): (0.00067, "Explicit unitary groups", "norm check raises (unit-corrected)"),
    (ACT, "test_public_centralizer_of_omega_is_generated_by_the_source_group"): (10.0, "Commutant centralizers", "CentralFieldCommutant second kind, Hermitian rank 5 over Z[omega], h=1, coset index 1: model match + rank-10 trace isometries ~10 s"),
    (ACT, "test_cusp_stabilizer_generators_fix_rho_and_decompose_a_conjugate"): (0.0067, "Explicit unitary groups", "Heisenberg word decomposition, <=20 5x5 ops = 0.1 s (unit-corrected)"),
    (ACT, "test_source_reflections_fixing_rho_decompose_in_the_cusp_stabilizer"): (0.0134, "Explicit unitary groups", "3 decompositions + null stabilizer = 0.2 s (unit-corrected)"),
    (ACT, "test_short_root_stabilizer_has_both_spinor_cosets"): (0.067, "Explicit unitary groups", "6 spinor characters rank 10 (0.06 s) + stabilizer = 1 s (unit-corrected)"),
    # Anisotropic binary Hermitian diag(1,-3) over Z[i]: SU = norm-one units of Z<1,i,j,k> in (-1,3).
    (_T + "test_anisotropic_binary_hermitian.sage", "test_the_form_is_anisotropic"): (0.05, "Quaternion Voronoi", "Hilbert symbol (-1,3)_3"),
    (_T + "test_anisotropic_binary_hermitian.sage", "test_generators_are_integral_and_unitary"): (0.5, "Quaternion Voronoi", "Voronoi of Z<1,i,j,k>, 1 perfect class: ~3 s (unit-corrected)"),
    (_T + "test_anisotropic_binary_hermitian.sage", "test_the_special_unitary_group_is_infinite"): (0.5, "Quaternion Voronoi", "same Voronoi route ~3 s (unit-corrected)"),
    (_T + "test_anisotropic_binary_hermitian.sage", "test_scalar_kernel_is_the_fourth_roots_of_unity"): (0.05, "Quaternion Voronoi", "one 2x2 scalar"),
    (_T + "test_definite_principal_ideals.sage", "test_definite_principal_and_nonprincipal_right_ideals_of_discriminant_23"): (0.1, "Quaternion Voronoi", "definite branch: short vectors nrd<=2 in rank 4 (Kirschmer-Voight), 2 qfminim = 0.1 s"),
    (_T + "test_effective_reduction.sage", "test_null_vectors_of_u_plus_u_are_joined"): (0.2, "Commutant centralizers", "box 1: 81 columns, null shells ~9, tuples ~10^3 = 0.2 s"),
    (_T + "test_effective_reduction.sage", "test_different_divisors_have_no_transporter"): (0.1, "Commutant centralizers", "divisor invariant 1 vs 2 separates: 0.1 s"),
    (_T + "test_factorwise_dedekind_transporter.sage", "test_gl_type_transporter_of_rank_two_over_the_golden_order"): (6.0, "Commutant centralizers", "levi_transporter, A/rad = M2(K)xM2(K), K=Q(sqrt5) h=1: Morita reduction + GL pair isomorphism ~6 s"),
    (_T + "test_factorwise_matrix_algebra_transporter.sage", "test_gl_type_transporter_over_a_matrix_algebra"): (6.0, "Commutant centralizers", "levi_transporter, Q[F]|X = M2(Q): index invariant rejects, conjugate found ~6 s"),
    (_T + "test_factorwise_nonmaximal_rank_two.sage", "test_gl_type_transporter_of_rank_two_over_a_nonmaximal_order"): (6.0, "Commutant centralizers", "levi_transporter, multiplier ring Z[sqrt5] index 2 in O_K: ~6 s"),
    (_T + "test_factorwise_order_transporter.sage", "test_gl_type_transporter_over_a_nonmaximal_order"): (6.0, "Commutant centralizers", "levi_transporter, 3 central idempotents, multiplier-ring invariant: ~6 s"),
    (_T + "test_fuchsian_domain.sage", "test_fuchsian_domain_of_the_discriminant_six_maximal_order"): (5.0, "Fuchsian domains", "Voight domain, area 2pi/3, signature (0;2,2,3,3): ~10 sides, enumeration to norm ~50: ~5 s"),
    # Euclidean Picard (Zhao) d=2,7,11 and Gaussian Picard (FFLP); 3x3 product ~1.5 ms.
    (_T + "test_euclidean_picard.sage", "test_source_matrices_are_integral_unitary"): (0.00402, "Explicit unitary groups", "3 fields x 6 membership checks = 0.05 s (unit-corrected)"),
    (_T + "test_euclidean_picard.sage", "test_source_cusp_presentations"): (0.00402, "Explicit unitary groups", "~40 3x3 products = 0.06 s (unit-corrected)"),
    (_T + "test_euclidean_picard.sage", "test_maximal_order_is_used_at_seven"): (0.00335, "Explicit unitary groups", "one membership (unit-corrected)"),
    (_T + "test_euclidean_picard.sage", "test_trace_form_signature"): (0.00402, "Explicit unitary groups", "3 rank-6 trace Grams = 0.05 s (unit-corrected)"),
    (_T + "test_euclidean_picard.sage", "test_words_recover_long_products"): (0.134, "Explicit unitary groups", "18 words x 14 letters, descent <=50 steps x 3x3 = 2 s (unit-corrected)"),
    (_T + "test_euclidean_picard.sage", "test_word_of_a_cusp_element_uses_no_inversion"): (0.0134, "Explicit unitary groups", "3 cusp words = 0.2 s (unit-corrected)"),
    (_T + "test_euclidean_picard.sage", "test_non_element_has_no_word"): (0.00335, "Explicit unitary groups", "membership raises (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_four_source_matrices_lie_in_the_group"): (0.00335, "Explicit unitary groups", "5 memberships (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_heisenberg_square_of_the_source_translation"): (0.00067, "Explicit unitary groups", "one product (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_cusp_element_word"): (0.00268, "Explicit unitary groups", "one cusp word = 0.04 s (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_reduction_recovers_long_products"): (0.0067, "Explicit unitary groups", "9-letter word, descent <=20 steps = 0.1 s (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_non_element_is_rejected_by_membership"): (0.00067, "Explicit unitary groups", "membership (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_trace_coordinates_are_order_linear_isometries"): (0.00402, "Explicit unitary groups", "rank-6 trace lattice, 4 generators = 0.05 s (unit-corrected)"),
    (_T + "test_gaussian_picard.sage", "test_generators_in_a_changed_basis_remain_order_linear_isometries"): (0.0067, "Explicit unitary groups", "4 generators x 36 trace pairings = 0.1 s (unit-corrected)"),
    # Congruence images: |image| bounded by |U_3(O/2O)|.
    (_T + "test_picard_congruence.sage", "test_gaussian_principal_subgroup_of_level_two"): (0.134, "Congruence images", "image order 48 (GAP on Z-coordinates mod 2, inv_images.py): BFS 48x4 products + 145 Schreier memberships ~2 s (unit-corrected)"),
    (_T + "test_picard_congruence.sage", "test_gaussian_preimage_of_a_cyclic_image"): (0.134, "Congruence images", "image 48, <M> image 2, index 24: ~2 s (unit-corrected)"),
    (_T + "test_picard_congruence.sage", "test_zhao_principal_subgroup_of_level_two"): (0.603, "Congruence images", "image order 168 (GAP, inv_images.py): 840 products + 673 Schreier memberships ~9 s (unit-corrected)"),
    (_T + "test_picard_null_transporters.sage", "test_gaussian_images_of_the_cusp_reduce"): (0.0268, "Explicit unitary groups", "8 null descents of height from 12-letter words = 0.4 s (unit-corrected)"),
    (_T + "test_picard_null_transporters.sage", "test_zhao_images_of_the_cusp_reduce"): (0.067, "Explicit unitary groups", "18 null descents = 1 s (unit-corrected)"),
    (_T + "test_nonmaximal_indefinite_hermitian.sage", "test_stabilizer_index_matches_the_finite_orbit"): (0.067, "Congruence images", "GF(2) image of the 6 moves has order 12 (inv_nonmax.py): index <= 12, ~1 s (unit-corrected)"),
    # Bestvina-Savin ocean; module fixtures charged to the first test using them.
    (OCEAN, "test_example_discriminant_and_anisotropy"): (0.05, "Bestvina-Savin ocean", "discriminant and Hilbert symbol"),
    (OCEAN, "test_inputs_outside_the_source_hypotheses_are_rejected"): (0.05, "Bestvina-Savin ocean", "3 hypothesis checks raise"),
    (OCEAN, "test_example_cell_structure_is_the_printed_one"): (5.0, "Bestvina-Savin ocean", "fixture: ocean of Delta=6, one cell orbit, 3 vertex orbits: ~5 s"),
    (OCEAN, "test_example_vertices_are_the_printed_ones"): (0.5, "Bestvina-Savin ocean", "vertex orbits of built fixture"),
    (OCEAN, "test_example_square_has_the_printed_vertex_labels"): (0.1, "Bestvina-Savin ocean", "4 invariants"),
    (OCEAN, "test_vertex_invariants_satisfy_the_edge_relation"): (5.0, "Bestvina-Savin ocean", "fixture: ocean of diag(1,-3), Delta=12, chi=-1/2: ~5 s"),
    (OCEAN, "test_example_euler_characteristic"): (0.05, "Bestvina-Savin ocean", "orbifold sum"),
    (OCEAN, "test_generators_preserve_the_form"): (0.5, "Bestvina-Savin ocean", "side pairings of both oceans"),
    (OCEAN, "test_scalar_kernel_is_not_trivial"): (0.05, "Bestvina-Savin ocean", "one scalar"),
    (OCEAN, "test_diagonal_euler_characteristic_is_the_mass_formula_value"): (0.05, "Bestvina-Savin ocean", "orbifold sum"),
    (OCEAN, "test_diagonal_form_agrees_with_the_quaternion_route"): (0.05, "Bestvina-Savin ocean", "matrix equality"),
    (OCEAN, "test_quaternion_generators_are_words_in_the_ocean_group"): (0.667, "Bestvina-Savin ocean", "quaternion Voronoi generators (~3 s) + ocean reduction ~1 s (unit-corrected)"),
    (OCEAN, "test_reduction_of_an_ocean_generator"): (0.2, "Bestvina-Savin ocean", "one ocean reduction 0.2 s"),
    (OCEAN, "test_images_modulo_split_primes_agree[5-2]"): (0.667, "Bestvina-Savin ocean", "Voronoi generators ~3 s + BFS in GL2(F5) (480 elements) (unit-corrected)"),
    (OCEAN, "test_images_modulo_split_primes_agree[13-5]"): (0.833, "Bestvina-Savin ocean", "Voronoi generators ~3 s + BFS in GL2(F13) (<=26208 elements x ~6 gens, ~0.5 s) (unit-corrected)"),
    (_T + "test_hermitian_parabolic_congruence.sage", "test_level_two_subgroup_of_the_setwise_stabilizer"): (1.0, "Congruence images", "stabilizer_generators (hermitian_lattices.py:592): ~15 lifts x 0.02 s + image mod 2 <= 64 (assumed bound) + Schreier memberships at 0.1 ms"),
    (_T + "test_hermitian_parabolic_congruence.sage", "test_pointwise_stabilizer_fixes_the_marked_vector"): (2.0, "Hermitian parabolics", "two stabilizer generator sets, Levi index 1: ~2 s"),
    (_T + "test_hermitian_parabolic_congruence.sage", "test_gaussian_line_stabilizer_equals_the_source_cusp_group_up_to_scalars"): (4.0, "Hermitian parabolics", "stabilizer gens + FFLP words + Levi closure of mu4 x mu4 (16) ~4 s"),
    # Quaternion Voronoi on Z<1,i,j,k> in (-1,3).
    (_T + "test_quaternion_voronoi.sage", "test_the_order_has_one_perfect_form_class_and_integral_generators"): (3.0, "Quaternion Voronoi", "1 perfect class, |Min| small, facets <=6, ~3 s"),
    (_T + "test_quaternion_voronoi.sage", "test_units_in_a_box_are_words_in_the_generators"): (1.33, "Quaternion Voronoi", "Voronoi 3 s + 11^4 = 14641 nrd evaluations 0.7 s + word ball radius 4 (~10^4) + descents ~4 s (unit-corrected)"),
    (_T + "test_quaternion_voronoi.sage", "test_the_cell_area_equals_the_mass_formula"): (0.5, "Quaternion Voronoi", "Voronoi + cell areas ~3 s (unit-corrected)"),
    (_T + "test_quaternion_voronoi.sage", "test_the_norm_one_extraction_applies_to_the_voronoi_generators"): (0.5, "Quaternion Voronoi", "Voronoi + Schreier index 2 ~3 s (unit-corrected)"),
    (_T + "test_quaternion_voronoi.sage", "test_the_full_orthogonal_group_generators_of_the_anisotropic_lattice"): (1.0, "Quaternion Voronoi", "Voronoi 3 s + normalizer box 2401 (2 s) + pairs (unit-corrected)"),
    (_T + "test_quaternion_voronoi.sage", "test_public_group_generators_of_the_anisotropic_lattice"): (12.0, "Quaternion Voronoi", "public O() dispatch to the quaternion route: recognition + same 6 s work + public wrapping"),
    # Witnesses.
    (_W + "test_bianchi_groups.sage", "test_swan_generators_of_gl2_over_sqrt_minus_five"): (6.0, "Bianchi groups", "Swan, Q(sqrt-5), h=2: one singular point class, ~20 generators, MatrixGroup over F3 (|GL2(F3)|=48) ~6 s"),
    (_W + "test_bianchi_groups.sage", "test_stabilizer_of_the_nonfree_lattice_o_plus_p"): (4.0, "Bianchi groups", "Steinitz class p, conjugated Swan generators ~4 s"),
    (_W + "test_hermitian_levi.sage", "test_split_gaussian_rank_two_levi_index"): (11.0, "Hermitian parabolics", "_levi_cosets (hermitian_lattices.py:693-717): 12 moves, index 3: 106 lifts x 0.1 s (45 value() evaluations of 5x5 at 0.4 ms + residue loop <= 256 x 25 entries in Python ints)"),
    (_W + "test_hermitian_levi.sage", "test_glued_eisenstein_rank_two_levi_index"): (210.0, "Hermitian parabolics", "_levi_cosets: 16 moves, index 8: ~700 lifts x 0.3 s (glue denominator 3 enlarges the residue loop of _admissible)"),
    (_W + "test_hermitian_negative_transporter_and_preimages.sage", "test_two_i_lines_with_isometric_quotients_have_no_transporter"): (0.3, "Hermitian parabolics", "transporter_to (l.738): ~14 lifts x 0.02 s"),
    (_W + "test_hermitian_negative_transporter_and_preimages.sage", "test_gaussian_line_stabilizer_meets_a_joint_residue_and_determinant_preimage[setwise]"): (0.5, "Congruence images", "stabilizer_generators: ~15 lifts x 0.02 s + finite_preimage_generators index 8 x ~10 gens"),
    (_W + "test_hermitian_negative_transporter_and_preimages.sage", "test_gaussian_line_stabilizer_meets_a_joint_residue_and_determinant_preimage[pointwise]"): (0.3, "Congruence images", "pointwise Levi 4 lifts x 0.02 s + preimage index 8"),
    (_W + "test_hermitian_parabolic_stabilizers.sage", "test_flag_stabilizers_based_and_setwise_on_the_split_frame"): (6.0, "Hermitian parabolics", "split Z[i]^5 flag (1,1): ~60 lifts x 0.1 s"),
    (_W + "test_hermitian_parabolic_stabilizers.sage", "test_setwise_line_stabilizer_of_the_glued_lattice_induces_the_hand_levi_image"): (1.2, "Hermitian parabolics", "k=1, n=3: 12 moves, index 2: ~60 lifts x 0.02 s"),
    (_W + "test_hermitian_parabolic_stabilizers.sage", "test_act_eisenstein_cusp_has_the_full_split_levi_image"): (400.0, "Hermitian parabolics", "levi_candidates (l.561): mu6 x U(I_3), U(I_3) from automorphisms of the A2^3 trace form (order 10368) filtered to 1296; levi_image: 7776 lifts x 0.05 s (n=5, k=1, 7 parameters, 36 value() evaluations at 0.4 ms + residue loop)"),
    (_W + "test_split_quaternary_orthogonal_number_field.sage", "test_split_quaternary_orthogonal_group_of_the_norm_form_of_eleven_nineteen_over_sqrt5"): (60.0, "Number-field orthogonal groups", "C0 = B x B, B=(11,19)_K ramified at the 2 primes over 11: Voronoi of a maximal order over Q(sqrt5) (herm dim 6) ~60 s"),
    (_W + "test_ternary_orthogonal_number_field.sage", "test_ternary_orthogonal_group_of_one_minus_eleven_minus_nineteen_over_sqrt5"): (30.0, "Number-field orthogonal groups", "C0 = B: one Voronoi of a maximal order over Q(sqrt5) ~30 s"),
    (_W + "test_li_margulis_box.sage", "test_isotropic_orbits_of_u_plus_u2"): (60.0, "Li-Margulis boxes", "isotropic_pair_grams: 3 Grams (d=1: a=0; d=2: a=0,2); box 1 (81 vectors): ~30 pairs each deduplicated by marked_transporter (pairing_slices.py:68) ~100 calls x 0.5 s (per-call cost assumed from isometry-testing.md) = 50 s"),
    (_W + "test_li_margulis_box.sage", "test_isotropic_orbits_of_u_plus_u"): (15.0, "Li-Margulis boxes", "1 Gram, ~30 pairs, ~30 marked_transporter calls x 0.5 s (assumed) = 15 s"),
    # Li-Margulis boxes (tests).
    (LMB, "test_null_vectors_of_u"): (0.2, "Li-Margulis boxes", "box 2 on rank 2: 25 vectors, 1 orbit"),
    (LMB, "test_two_orbits_of_square_minus_two"): (8.0, "Li-Margulis boxes", "box 2 on U+A1: exactly 10 vectors of square -2 (ab = -1 - c^2), ~15 marked_transporter calls x 0.5 s (assumed)"),
    (LMB, "test_rank_two_tuple"): (8.0, "Li-Margulis boxes", "box 1 on U+A1: 8 primitive null vectors, ~16 hyperbolic pairs, ~16 marked_transporter calls x 0.5 s (assumed)"),
    (LMB, "test_isometry_found_and_absent"): (0.1, "Li-Margulis boxes", "rank 2, boxes 1 and 3: <=49 vectors, tuples <=10^3 = 0.1 s"),
    (LMB, "test_generator_bound_formula"): (0.01, "Li-Margulis boxes", "closed formula"),
    (LMB, "test_isometry_bound_exponents"): (0.01, "Li-Margulis boxes", "closed formula"),
    (LMB, "test_prescribed_null_gram_on_two_planes"): (16.0, "Li-Margulis boxes", "box 1 on U+U: ~32 primitive null vectors, ~32 marked_transporter calls x 0.5 s (assumed)"),
    (LMB, "test_isotropic_pair_grams_of_u_plus_u_three"): (0.1, "Li-Margulis boxes", "divisors of exponent 3: 4 Grams"),
    (LMB, "test_theorem_one_constant_composition"): (0.01, "Li-Margulis constants", "closed formula"),
    (LMB, "test_theorem_two_constant_composition"): (0.01, "Li-Margulis constants", "closed formula"),
    (LMB, "test_embedding_bound_of_a_line"): (0.01, "Li-Margulis boxes", "closed formula"),
    (LMB, "test_definite_form_is_rejected"): (0.01, "Li-Margulis boxes", "signature check raises"),
    (LMB, "test_orthogonal_group_of_u"): (0.1, "Li-Margulis boxes", "box 1 on U: 9 vectors, 4 isometries"),
    (LMB, "test_quaternion_unit_of_issue_40"): (3.0, "Li-Margulis boxes", "box 3 on rank 4: 2401 vectors, shells q=2,-6 of size <=40, tuples <=10^5 pairing checks ~3 s"),
}

for _name in (
    "test_appendix_epsilon", "test_corollary1_and_lemma20", "test_default_constants", "test_haar_ball_and_theorem5",
    "test_harish_chandra_decay_constant", "test_lemma10_radius", "test_lemma11_constant", "test_lemma14_constant",
    "test_lemma15_constant_power", "test_lemma15_radius_and_zeta", "test_lemma24_and_corollaries",
    "test_matrix_coefficient_constant", "test_section6_constant", "test_theorem6_constant",
):
    TABLE[(LMC, _name)] = (0.01, "Li-Margulis constants", "closed rational formulas on integers of <=10^3 digits")

# Commutant groups (commutant_groups.py:1-60 dispatch).
_COMM = {
    "test_unipotent_transporter_to_the_inverse_exists": (1.0, "radical branch A/rad=Q: coset t_r(+-N) integral point ~1 s"),
    "test_unipotent_transporter_to_the_square_is_absent": (1.0, "radical branch, empty coset ~1 s"),
    "test_symplectic_quaternion_centralizer_is_the_norm_one_units": (0.667, "QuaternionCommutant symplectic: O^1 from Voronoi of (-1,3) order (1 class) ~4 s (unit-corrected)"),
    "test_symplectic_quaternion_transporters": (1.0, "symplectic transporter: Voronoi + 3 reduced-norm point searches ~6 s (unit-corrected)"),
    "test_orthogonal_quaternion_centralizer_is_the_real_torus": (2.0, "QuaternionCommutant orthogonal: norm-one torus of Q(c) (one PARI bnfinit deg 2) + anticommuting coset ~2 s"),
    "test_orthogonal_quaternion_transporters": (3.0, "orthogonal transporter: torus + norm equations ~3 s"),
    "test_gaussian_hyperbolic_centralizer_contains_the_scalar_units": (2.0, "CentralFieldCommutant second kind, hyperbolic over Z[i], unitary model U(1,1) ~2 s"),
    "test_real_quadratic_symplectic_centralizer": (3.0, "CentralFieldCommutant alternating n=2 over real quadratic: SL2(O_K) model ~3 s"),
    "test_gaussian_hyperbolic_transporters": (4.0, "hyperbolic model isometry + coset reps (index 1) x 3 calls ~4 s"),
    "test_anisotropic_gaussian_centralizer_is_infinite": (6.0, "ocean branch: Bestvina-Savin ocean of the binary form ~6 s"),
    "test_anisotropic_gaussian_transporters": (8.0, "ocean model isometry + coset search, 3 calls ~8 s"),
    "test_real_quadratic_orthogonal_centralizer": (3.0, "EtaleCommutant: relative norm equations over K^sigma (bnfinit, bnfisnorm) ~3 s"),
    "test_vector_orbit_in_the_hyperbolic_plane": (2.0, "VectorOrbit: marked transporter of a cyclic span ~2 s"),
    "test_real_quadratic_orthogonal_transporters": (4.0, "etale transporters: norm equations, 2 calls ~4 s"),
    "test_product_levi_quotient_transporter": (3.0, "levi_transporter product Levi quotient ~3 s"),
    "test_anisotropic_eisenstein_centralizer_is_infinite": (1.0, "_BinaryHermitianOrder over Z[omega], h=1: quaternion-order Voronoi ~6 s (unit-corrected)"),
    "test_anisotropic_eisenstein_transporters": (8.0, "binary Hermitian isometry + coset search, 3 calls ~8 s"),
    "test_hermitian_quaternion_centralizer_is_infinite": (12.0, "quaternion_hermitian_commutant M2(D): Hermitian over the right order ~12 s"),
    "test_golden_orthogonal_centralizer_rank_four": (12.0, "symmetric rank >= 3 over Q(sqrt5): H+H+L0 split, Eichler transvections + spinor reps ~12 s"),
    "test_split_orthogonal_centralizer_of_minus_one_on_U": (1.0, "SplitSimpleCommutant M2(Q): O(U) finite (order 4) ~1 s"),
    "test_split_orthogonal_centralizer_rank_three_is_infinite": (4.0, "SplitSimpleCommutant M3(Q): O(b) of rank-3 indefinite form ~4 s"),
    "test_split_symplectic_centralizer_is_SL2": (1.0, "SplitSimpleCommutant alternating: Sp2(Z)=SL2(Z) generators ~1 s"),
    "test_two_pell_units_transporters": (4.0, "EtaleCommutant two factors, Pell units via bnfinit: 3 calls ~4 s"),
    "test_swapped_factor_transporter": (3.0, "levi_transporter with swapped adjoint factors ~3 s"),
    "test_split_orthogonal_transporters": (2.0, "split orthogonal: form isometry + coset reps ~2 s"),
    "test_split_symplectic_transporters": (2.0, "split symplectic: form isometry, 2 calls ~2 s"),
    "test_definite_hermitian_quaternion_centralizer_is_infinite": (12.0, "M2(D) definite D, Hermitian plane over right order ~12 s"),
    "test_ramified_dyadic_eichler_generators": (6.0, "eichler_orthogonal with ramified dyadic prime: realizing reflection search ~6 s"),
    "test_skew_hermitian_quaternion_centralizer_is_infinite": (12.0, "skew-Hermitian over quaternion order, M2(D) ~12 s"),
    "test_diagonal_null_vector_beyond_ternary_subforms": (0.4, "qfsolve on diagonal form of rank >= 4 ~0.4 s"),
    "test_dyadic_spinor_norm_groups_of_binary_two_adic_lattices": (0.4, "closed 2-adic spinor norm tables (Earnest-Hsia) ~0.4 s"),
    "test_dyadic_spinor_norm_group_after_maximal_norm_splitting": (0.4, "2-adic splitting + table ~0.4 s"),
    "test_voronoi_perfect_forms_over_real_quadratic_field": (5.0, "SemisimpleOrderVoronoi split branch over Q(sqrt2), herm dim 6: P<=10 classes x F<=30 facets (assumed) x 15 probes x 1 ms Gram+qfminim"),
    "test_voronoi_sl1_of_lipschitz_order_over_sqrt2": (5.0, "same Voronoi (P, F assumed) + norm-one Schreier extraction"),
    "test_anisotropic_skew_hermitian_quaternion_centralizer": (12.0, "anisotropic skew-Hermitian over quaternion order ~12 s"),
    "test_anisotropic_skew_hermitian_quaternion_transporters": (18.0, "transporters: model isometry + coset search, 2 calls ~18 s"),
    "test_maximal_orders_containing_the_lipschitz_order_over_sqrt2": (4.0, "OrderLattices over Q(sqrt2): overorders of index 2^k ~4 s"),
    "test_witt_one_parabolic_generators_over_sqrt2": (12.0, "_witt_one_model over Q(sqrt2), h=1: parabolic generators ~12 s"),
    "test_witt_one_parabolic_generators_mixed_sign_rank_one": (12.0, "_witt_one_model mixed signs, rank-one anisotropic kernel ~12 s"),
    "test_witt_one_pseudo_basis_over_sqrt10": (12.0, "_witt_one_model over Q(sqrt10), h=2: Steinitz pseudo-basis ~12 s"),
    "test_binary_orthogonal_group_of_one_minus_three_over_sqrt5": (6.0, "binary orthogonal: torus K[c]^1 via relative norm equations + reflection coset ~6 s"),
    "test_factorwise_centralizer_and_transporter": (6.0, "factorwise: hyperbolic + Pell factors, centralizer + 2 transporters ~6 s"),
    "test_factorwise_centralizer_with_a_gl_type_pair": (4.0, "factorwise GL pair centralizer ~4 s"),
    "test_factorwise_transporter_with_a_gl_type_pair": (6.0, "factorwise GL pair transporter, 2 calls ~6 s"),
}
for _name, (_s, _b) in _COMM.items():
    TABLE[(COMM, _name)] = (_s, "Commutant centralizers", _b)


def lookup(path: str, function: str) -> tuple[float, str]:
    seconds, section, basis = TABLE[(path, function)]
    return seconds, f"[{section}] {basis}"


def _key(path: str) -> str:
    return path[len("tests/"):] if path.startswith("tests/") else path


def _missing_nodes() -> list[str]:
    path = Path(__file__).resolve().parent / "missing_E.txt"
    return [n for n in path.read_text().splitlines() if n]


def _scope_nodes() -> list[str]:
    lines = (Path(__file__).resolve().parent / "nodes_initial.txt").read_text().splitlines()
    return [n for n in lines if _key(n.split("::")[0]) in MODELS]


def estimates() -> dict[str, tuple[float, str]]:
    result = {}
    for node in _scope_nodes():
        path, func = node.split("::")
        model = MODELS[_key(path)]
        seconds, basis = model(func.split("[")[0])
        result[node] = (seconds, f"[{SECTIONS[model.__name__]}] {basis}")
    for node in _missing_nodes():
        path, func = node.split("::")
        result[node] = lookup(path, func)
    return result


if __name__ == "__main__":
    import statistics

    table = estimates()
    assert set(table) == set(_scope_nodes()) | set(_missing_nodes())
    assert not set(_scope_nodes()) & set(_missing_nodes())
    print(len(table))
    files = sorted({n.split("::")[0] for n in table})
    for file in files:
        values = [s for n, (s, _) in table.items() if n.split("::")[0] == file]
        print(file, len(values), min(values), statistics.median(values), max(values))
