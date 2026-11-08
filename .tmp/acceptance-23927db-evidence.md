# Partial Acceptance reports at 23927db/16485a1

[Acceptance 42412497409](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/42412497409) records port `23927dbb9f3bedc9cb765c2dc69dde7a748420d1` and provider `16485a149e2ac30c660d2eb3083da2b01731fab2`. This is the inspected result set of the completed failed run. It contains push-selected samples; full-file final-pair evidence remains required. A testcase marked passed below has no failure, error or skip in its inspected JUnit report. Remaining shards and final-pair evidence remain required.

| Family | Outcome | Inspected cases |
| --- | --- | ---: |
| boundaries | failure | 2 |
| boundaries | passed | 8 |
| centralizers | passed | 5 |
| double_cosets | passed | 19 |
| isometry | passed | 166 |
| isotropic_orbits | passed | 9 |
| perfect_domains | passed | 92 |
| orthogonal_group | passed | 8 |
| reflective | passed | 5 |
| vector_orbits | passed | 30 |
| vector_stabilizers | passed | 25 |

## Failed assertions

- `boundaries-2`: `test_dawes_buildings_have_the_published_counts[dawes_L2_gamma2]` returns `(2, 2, 3)` instead of `(3, 2, 4)`. The line and flag counts disagree; the plane count agrees. The pinned line-orbit lemma and the independent argument in `references/dawes-L2-boundary.md` establish (2,2,3). The source theorem cited by the fixture contradicts that lemma. Issue #33 comment 6716895523 records the proposed correction; its application is pending a decision against the unchanged-suite requirement. The frozen assertion still fails.
- `boundaries-3`: `test_dawes_building_incidence_graph_is_the_published_figure[dawes_2U2_A2_stable]` fails at `lattice.base_ring()(rank)` in `_isotropic_orbit_representatives`: the rank is a native Sage integer, rejected as a scalar of the owned coefficient ring. The graph assertion is unreached. Rank is an integer index; accepting it must preserve integral-index validation.

## Inspected artifacts without JUnit

- `acceptance-boundaries-11`
- `acceptance-boundaries-12`
- `acceptance-boundaries-13`
- `acceptance-boundaries-14`
- `acceptance-boundaries-15`
- `acceptance-boundaries-16`
- `acceptance-boundaries-8`

- `acceptance-isotropic_orbits-3`
- `acceptance-boundaries-17`
- `acceptance-isotropic_orbits-2`
- `acceptance-orthogonal_group-1`
- `acceptance-isotropic_orbits-12`
- `acceptance-isotropic_orbits-6`

## Recorded case IDs

| Artifact | Case | Outcome |
| --- | --- | --- |
| `acceptance-boundaries-0` | `test_dawes_buildings_have_the_published_counts[dawes_2U_A2_stable]` | passed |
| `acceptance-boundaries-1` | `test_dawes_buildings_have_the_published_counts[dawes_2U_minus2_minus6]` | passed |
| `acceptance-boundaries-10` | `test_dtower_boundary_pictures[lambda_10]` | passed |
| `acceptance-boundaries-2` | `test_dawes_buildings_have_the_published_counts[dawes_L2_gamma2]` | failure |
| `acceptance-boundaries-3` | `test_dawes_building_incidence_graph_is_the_published_figure[dawes_2U2_A2_stable]` | failure |
| `acceptance-boundaries-4` | `test_dawes_last_index_in_the_chain_is_the_discriminant_image_of_o_plus` | passed |
| `acceptance-boundaries-5` | `test_dtower_boundary_component_counts[N_3]` | passed |
| `acceptance-boundaries-6` | `test_dtower_boundary_component_counts[N_4]` | passed |
| `acceptance-boundaries-7` | `test_dtower_boundary_component_counts[N_9]` | passed |
| `acceptance-boundaries-9` | `test_dtower_boundary_pictures[lambda_9]` | passed |
| `acceptance-centralizers-0` | `test_centralizer_images_match_oscar[oscar_lattices_with_isometry_110]` | passed |
| `acceptance-centralizers-0` | `test_centralizer_images_match_oscar[oscar_lattices_with_isometry_118]` | passed |
| `acceptance-centralizers-0` | `test_centralizer_images_match_oscar[oscar_lattices_with_isometry_133]` | passed |
| `acceptance-centralizers-0` | `test_centralizer_images_match_oscar[oscar_enumeration_88]` | passed |
| `acceptance-centralizers-0` | `test_centralizer_images_match_oscar[oscar_enumeration_96]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n55_big5760_sma1440_vf47_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n55_big5760_sma960_vf47_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n56_big2592_sma288_vf58_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n56_big40320_sma5040_vf19_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n56_big768_sma192_vf158_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n56_big768_sma384_vf158_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n56_big768_sma384_vf158_idx1]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n60_big1440_sma288_vf168_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n60_big1440_sma480_vf168_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n60_big34560_sma5760_vf30_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n60_big34560_sma8640_vf30_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n60_big5760_sma960_vf58_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n64_big1536_sma384_vf324_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n70_big80640_sma10080_vf36_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n72_big18432_sma2304_vf150_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n72_big3456_sma576_vf409_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n80_big34560_sma5760_vf156_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_double_cosets_partition_the_rational_group[DoubleCoset_n96_big92160_sma23040_vf414_idx0]` | passed |
| `acceptance-double_cosets-0` | `test_the_recorded_case_families_are_all_present` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism4]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism24]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism44]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism49]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism54]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism55]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism59]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism87]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism89]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism403]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism405]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism427]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism430]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism441]` | passed |
| `acceptance-isometry-0` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism442]` | passed |
| `acceptance-isometry-0` | `test_row_matrix_morphism_roundtrip_uses_coordinate_labels` | passed |
| `acceptance-isometry-0` | `test_eichler_generator_families_preserve_nonlexical_complement` | passed |
| `acceptance-isometry-0` | `test_two_u_full_generating_data_accepts_raw_sage_integer_norm` | passed |
| `acceptance-isometry-0` | `test_large_cell_transporter_invariant_rejection_returns_none` | passed |
| `acceptance-isometry-0` | `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U_A3]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism10]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism36]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism37]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism43]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism72]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism76]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism83]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism409]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism431]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism433]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism439]` | passed |
| `acceptance-isometry-1` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism443]` | passed |
| `acceptance-isometry-1` | `test_column_matrix_morphism_roundtrip_uses_coordinate_labels` | passed |
| `acceptance-isometry-1` | `test_degenerate_isometry_witness_respects_complement_coordinates` | passed |
| `acceptance-isometry-1` | `test_odd_eichler_envelope_uses_declared_coordinate_order` | passed |
| `acceptance-isometry-2` | `test_isometry_verdicts_match_certified_pairs[indefinite_jl_U_2U_2I3]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism1]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism2]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism9]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism12]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism13]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism14]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism15]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism27]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism34]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism35]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism41]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism65]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism67]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism71]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism73]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism92]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism99]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism406]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism416]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism422]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism429]` | passed |
| `acceptance-isometry-2` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism435]` | passed |
| `acceptance-isometry-2` | `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U_A2_A2]` | passed |
| `acceptance-isometry-2` | `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U_2E8]` | passed |
| `acceptance-isometry-3` | `test_isometry_verdicts_match_certified_pairs[hecke_U_vs_U2]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism5]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism23]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism29]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism33]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism45]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism46]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism56]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism60]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism61]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism66]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism77]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism85]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism402]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism408]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism410]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism418]` | passed |
| `acceptance-isometry-3` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism436]` | passed |
| `acceptance-isometry-4` | `test_isometry_verdicts_match_certified_pairs[indefinite_jl_U_I3]` | passed |
| `acceptance-isometry-4` | `test_isometry_verdicts_match_certified_pairs[hecke_U2_A2_by_u]` | passed |
| `acceptance-isometry-4` | `test_isometry_verdicts_match_certified_pairs[hecke_3x3_swapped_basis]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism26]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism32]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism38]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism50]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism51]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism63]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism78]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism80]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism91]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism94]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism95]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism97]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism400]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism404]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism434]` | passed |
| `acceptance-isometry-4` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism437]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism6]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism7]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism17]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism21]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism25]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism30]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism39]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism47]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism57]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism70]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism74]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism75]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism79]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism84]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism93]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism98]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism407]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism411]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism417]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism424]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism425]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism432]` | passed |
| `acceptance-isometry-5` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism438]` | passed |
| `acceptance-isometry-5` | `test_equal_gram_groupoid_identification_preserves_coordinate_labels` | passed |
| `acceptance-isometry-5` | `test_conway_sloane_same_genus_pair_is_not_isometric` | passed |
| `acceptance-isometry-6` | `test_isometry_verdicts_match_certified_pairs[indefinite_jl_U_E8]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism3]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism22]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism58]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism62]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism86]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism412]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism414]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism415]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism420]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism426]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism440]` | passed |
| `acceptance-isometry-6` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism444]` | passed |
| `acceptance-isometry-6` | `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U_A2]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism0]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism8]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism11]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism16]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism18]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism19]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism20]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism28]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism31]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism40]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism42]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism48]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism52]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism53]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism64]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism68]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism69]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism81]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism82]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism88]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism90]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism96]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism401]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism413]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism419]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism421]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism423]` | passed |
| `acceptance-isometry-7` | `test_lorentzian_pairs_are_isometric[LORENTZ_TestIsomorphism428]` | passed |
| `acceptance-isometry-7` | `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U]` | passed |
| `acceptance-isometry-7` | `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_U_E7]` | passed |
| `acceptance-isotropic_orbits-0` | `test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_0]` | passed |
| `acceptance-isotropic_orbits-1` | `test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors[IsotropicCases_15]` | passed |
| `acceptance-isotropic_orbits-11` | `test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank0-primitive_isotropic_vector_orbits]` | passed |
| `acceptance-isotropic_orbits-7` | `test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1]` | passed |
| `acceptance-isotropic_orbits-5` | `test_published_representatives_lie_in_distinct_o_plus_orbits[rank0]` | passed |
| `acceptance-isotropic_orbits-8` | `test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_2]` | passed |
| `acceptance-isotropic_orbits-10` | `test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_2_2]` | passed |
| `acceptance-isotropic_orbits-4` | `test_o_n_maps_onto_o_of_the_discriminant_form` | passed |
| `acceptance-perfect_domains-0` | `test_large_pairing_graph_fallback_returns_exact_transporter` | passed |
| `acceptance-perfect_domains-0` | `test_non_diagonal_positive_shell_uses_the_pinned_cvp_center` | passed |
| `acceptance-perfect_domains-0` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_19]` | passed |
| `acceptance-perfect_domains-0` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_20]` | passed |
| `acceptance-perfect_domains-0` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_28]` | passed |
| `acceptance-perfect_domains-0` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_36]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_4]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_5]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_12]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_18]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_19]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_23]` | passed |
| `acceptance-perfect_domains-0` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_26]` | passed |
| `acceptance-perfect_domains-0` | `test_lazy_neighbor_cones_retain_exact_transformed_rays` | passed |
| `acceptance-perfect_domains-1` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_0]` | passed |
| `acceptance-perfect_domains-1` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_18]` | passed |
| `acceptance-perfect_domains-1` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_22]` | passed |
| `acceptance-perfect_domains-1` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_25]` | passed |
| `acceptance-perfect_domains-1` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_27]` | passed |
| `acceptance-perfect_domains-1` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_0]` | passed |
| `acceptance-perfect_domains-1` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_11]` | passed |
| `acceptance-perfect_domains-1` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_16]` | passed |
| `acceptance-perfect_domains-1` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_22]` | passed |
| `acceptance-perfect_domains-2` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_8]` | passed |
| `acceptance-perfect_domains-2` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_15]` | passed |
| `acceptance-perfect_domains-2` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_31]` | passed |
| `acceptance-perfect_domains-2` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_32]` | passed |
| `acceptance-perfect_domains-2` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_9]` | passed |
| `acceptance-perfect_domains-2` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_10]` | passed |
| `acceptance-perfect_domains-2` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_14]` | passed |
| `acceptance-perfect_domains-2` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_27]` | passed |
| `acceptance-perfect_domains-2` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_32]` | passed |
| `acceptance-perfect_domains-2` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_34]` | passed |
| `acceptance-perfect_domains-2` | `test_facet_orbit_permutation_view_materializes_exact_isometries` | passed |
| `acceptance-perfect_domains-2` | `test_large_cell_cone_uses_exact_normaliz_h_presentation` | passed |
| `acceptance-perfect_domains-3` | `test_positive_shell_survives_a_permuted_gram_basis` | passed |
| `acceptance-perfect_domains-3` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_16]` | passed |
| `acceptance-perfect_domains-3` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_6]` | passed |
| `acceptance-perfect_domains-3` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_15]` | passed |
| `acceptance-perfect_domains-3` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_24]` | passed |
| `acceptance-perfect_domains-4` | `test_positive_shell_scales_rational_direction_pairing_bound` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_1]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_3]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_4]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_5]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_6]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_9]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_10]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_11]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_17]` | passed |
| `acceptance-perfect_domains-4` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_35]` | passed |
| `acceptance-perfect_domains-4` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_7]` | passed |
| `acceptance-perfect_domains-4` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_25]` | passed |
| `acceptance-perfect_domains-4` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_28]` | passed |
| `acceptance-perfect_domains-4` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_29]` | passed |
| `acceptance-perfect_domains-4` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_31]` | passed |
| `acceptance-perfect_domains-4` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_39]` | passed |
| `acceptance-perfect_domains-5` | `test_non_diagonal_positive_shell_respects_explicit_coordinate_labels` | passed |
| `acceptance-perfect_domains-5` | `test_non_diagonal_isotropic_shell_uses_the_exact_sphere_center` | passed |
| `acceptance-perfect_domains-5` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_7]` | passed |
| `acceptance-perfect_domains-5` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_34]` | passed |
| `acceptance-perfect_domains-5` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_1]` | passed |
| `acceptance-perfect_domains-5` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_8]` | passed |
| `acceptance-perfect_domains-5` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_13]` | passed |
| `acceptance-perfect_domains-5` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_37]` | passed |
| `acceptance-perfect_domains-5` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_38]` | passed |
| `acceptance-perfect_domains-6` | `test_non_diagonal_isotropic_shell_respects_explicit_coordinate_labels` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_12]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_13]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_21]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_29]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_30]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_33]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_38]` | passed |
| `acceptance-perfect_domains-6` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_39]` | passed |
| `acceptance-perfect_domains-6` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_3]` | passed |
| `acceptance-perfect_domains-6` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_20]` | passed |
| `acceptance-perfect_domains-6` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_21]` | passed |
| `acceptance-perfect_domains-6` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_33]` | passed |
| `acceptance-perfect_domains-6` | `test_traversal_retains_only_adjacency_facet_orbits` | passed |
| `acceptance-perfect_domains-7` | `test_isotropic_shell_survives_a_permuted_gram_basis` | passed |
| `acceptance-perfect_domains-7` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_2]` | passed |
| `acceptance-perfect_domains-7` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_14]` | passed |
| `acceptance-perfect_domains-7` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_23]` | passed |
| `acceptance-perfect_domains-7` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_24]` | passed |
| `acceptance-perfect_domains-7` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_26]` | passed |
| `acceptance-perfect_domains-7` | `test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_37]` | passed |
| `acceptance-perfect_domains-7` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_2]` | passed |
| `acceptance-perfect_domains-7` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_17]` | passed |
| `acceptance-perfect_domains-7` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_30]` | passed |
| `acceptance-perfect_domains-7` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_35]` | passed |
| `acceptance-perfect_domains-7` | `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_36]` | passed |
| `acceptance-vector_orbits-2` | `test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_30]` | passed |
| `acceptance-vector_orbits-1` | `test_primitive_orbit_counts_in_u_plus_e8_match_dutour_sikiric_hulek[2d_2]` | passed |
| `acceptance-reflective-0` | `test_edgewalk_branch_representative_sample` | passed |
| `acceptance-reflective-0` | `test_edgewalk_period_action_matrix_uses_coordinate_label_order` | passed |
| `acceptance-reflective-0` | `test_binary_cycle_coordinate_order_and_period_certificate` | passed |
| `acceptance-reflective-0` | `test_vinberg_roots_match_oscars_chamber_up_to_isometry` | passed |
| `acceptance-reflective-0` | `test_e10_simple_roots_match_dutour_sikiric_hulek` | passed |
| `acceptance-orthogonal_group-5` | `test_mertens_omega_and_minus_one_generate_o_l` | passed |
| `acceptance-orthogonal_group-5` | `test_zero_form_orthogonal_group_is_the_full_general_linear_group` | passed |
| `acceptance-orthogonal_group-5` | `test_orthogonal_group_of_the_hyperbolic_plane_has_order_four` | passed |
| `acceptance-orthogonal_group-5` | `test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q` | passed |
| `acceptance-orthogonal_group-4` | `test_degenerate_orthogonal_group_embeds_non_diagonal_complement` | passed |
| `acceptance-vector_orbits-0` | `test_anisotropic_binary_vector_classes_use_primitive_embeddings` | passed |
| `acceptance-isotropic_orbits-9` | `test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_2_1]` | passed |
| `acceptance-vector_orbits-3` | `test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[0]` | passed |
| `acceptance-orthogonal_group-2` | `test_degenerate_orthogonal_group_contains_radical_shear_and_gl_action` | passed |
| `acceptance-orthogonal_group-3` | `test_degenerate_orthogonal_group_generates_rank_two_radical_action` | passed |
| `acceptance-orthogonal_group-0` | `test_lorentzian_generators_match_polyhedral_common[LORENTZ_OrthogonalGroup0]` | passed |

| Newly inspected index assertion | Outcome | Seconds |
| --- | --- | ---: |
| I_2_1 | passed | 3002.232 |

The new orthogonal-group reports include LORENTZ_OrthogonalGroup0 and Mertens. A7 and the rank-one-labelled published-representative and I_2_10 plane-class shards supply logs without JUnit; their final assertions remain unproved. All reports above retain the port/provider pair in the heading.


## Additional vector assertions

The following archives retain the exact port/provider commits above. Every inspected JUnit case has no failure, error or skip, and its case ID occurs in the matching shard log. These are push-selected cases; the complete frozen file remains required at the final pair. The successful vector-stabilizer shards establish collection and their listed assertions, while the other shards remain open.

| Artifact | ID | Case | Outcome | JUnit seconds |
| --- | ---: | --- | --- | ---: |
| `acceptance-vector_orbits-4` | 13603890524 | `test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[2]` | passed | 991.473 |
| `acceptance-vector_orbits-5` | 13603485801 | `test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[4]` | passed | 562.383 |
| `acceptance-vector_orbits-6` | 13603975213 | `test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[6]` | passed | 618.115 |
| `acceptance-vector_orbits-7` | 13603293303 | `test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[8]` | passed | 842.340 |
| `acceptance-vector_orbits-8` | 13604031672 | `test_e10_orbit_counts_match_brandhorst_gonzalez_alonso[10]` | passed | 774.395 |
| `acceptance-vector_orbits-9` | 13603598221 | `test_published_representatives_are_one_per_orbit[2]` | passed | 1125.670 |
| `acceptance-vector_orbits-10` | 13604215936 | `test_published_representatives_are_one_per_orbit[4]` | passed | 837.065 |
| `acceptance-vector_orbits-11` | 13603837516 | `test_published_representatives_are_one_per_orbit[6]` | passed | 377.424 |
| `acceptance-vector_orbits-12` | 13604578472 | `test_published_representatives_are_one_per_orbit[8]` | passed | 905.537 |
| `acceptance-vector_orbits-13` | 13604518555 | `test_published_representatives_are_one_per_orbit[10]` | passed | 972.211 |
| `acceptance-vector_orbits-14` | 13604149686 | `test_published_representatives_are_one_per_orbit[12]` | passed | 694.801 |
| `acceptance-vector_orbits-15` | 13604768083 | `test_published_representatives_are_one_per_orbit[14]` | passed | 987.725 |
| `acceptance-vector_orbits-16` | 13605596123 | `test_published_representatives_are_one_per_orbit[16]` | passed | 864.865 |
| `acceptance-vector_orbits-17` | 13605980532 | `test_published_representatives_are_one_per_orbit[18]` | passed | 763.108 |
| `acceptance-vector_orbits-18` | 13605039569 | `test_published_representatives_are_one_per_orbit[20]` | passed | 809.294 |
| `acceptance-vector_orbits-19` | 13605391608 | `test_published_representatives_are_one_per_orbit[22]` | passed | 423.112 |
| `acceptance-vector_orbits-20` | 13605692454 | `test_published_representatives_are_one_per_orbit[24]` | passed | 490.716 |
| `acceptance-vector_orbits-21` | 13606261334 | `test_published_representatives_are_one_per_orbit[26]` | passed | 711.443 |
| `acceptance-vector_orbits-22` | 13606034850 | `test_published_representatives_are_one_per_orbit[28]` | passed | 1576.437 |
| `acceptance-vector_orbits-23` | 13606590529 | `test_published_representatives_are_one_per_orbit[30]` | passed | 1103.091 |
| `acceptance-vector_orbits-24` | 13605519902 | `test_i_2_10_is_transitive_on_norm_minus_one_vectors` | passed | 366.623 |
| `acceptance-vector_orbits-25` | 13606208135 | `test_k3_lattice_has_one_primitive_orbit_per_norm[0]` | passed | 113.482 |
| `acceptance-vector_orbits-26` | 13606283297 | `test_k3_lattice_has_one_primitive_orbit_per_norm[-2]` | passed | 181.965 |
| `acceptance-vector_orbits-27` | 13606581457 | `test_k3_lattice_has_one_primitive_orbit_per_norm[2]` | passed | 190.673 |
| `acceptance-vector_orbits-28` | 13606551625 | `test_k3_lattice_has_one_primitive_orbit_per_norm[4]` | passed | 164.184 |
| `acceptance-vector_orbits-29` | 13607285226 | `test_k3_lattice_has_one_primitive_orbit_per_norm[-4]` | passed | 191.283 |

## Completed vector-stabilizer assertions

Every vector-stabilizer shard now has an inspected JUnit report, matching log and exact port/provider files for the pair above. All reports have no failures, errors or skips. The table gives every selected node and its JUnit call time. The E10 assertions compare the complete multiset of stabilizer-image indices over every vector class at each displayed norm with the published fixture. Reflected-vector tests verify the actual pointed isometry witnesses. The other nodes check the stated integral isotropic actions, containing bounds, markings and negative decisions. The OSCAR-order assertion compares the actual vector stabilizer order with its fixture. These are the run's selected assertions; the full frozen module remains required at the final pair.

| Artifact | Case | JUnit seconds |
| --- | --- | ---: |
| `acceptance-vector_stabilizers-0` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_-1]` | 1917.173 |
| `acceptance-vector_stabilizers-0` | `test_isotropic_pointed_automorphism_generators_preserve_marking_and_form` | 7.579 |
| `acceptance-vector_stabilizers-1` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_1]` | 335.068 |
| `acceptance-vector_stabilizers-1` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_5]` | 6.370 |
| `acceptance-vector_stabilizers-2` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_7]` | 341.216 |
| `acceptance-vector_stabilizers-2` | `test_general_isotropic_pointed_isometry_carries_negative_marking` | 0.553 |
| `acceptance-vector_stabilizers-3` | `test_divisibility_four_isotropic_plane_stabilizer_constructs_exact_containing_bound` | 0.905 |
| `acceptance-vector_stabilizers-3` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_0]` | 524.983 |
| `acceptance-vector_stabilizers-4` | `test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[10]` | 755.403 |
| `acceptance-vector_stabilizers-4` | `test_split_divisibility_one_pointed_isometry_carries_actual_negative_marking` | 0.604 |
| `acceptance-vector_stabilizers-5` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_2]` | 275.880 |
| `acceptance-vector_stabilizers-5` | `test_isotropic_pointed_isometry_rejects_different_divisibility` | 0.018 |
| `acceptance-vector_stabilizers-6` | `test_divisibility_two_isotropic_pointed_isometry_has_exact_witness` | 0.848 |
| `acceptance-vector_stabilizers-6` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_4]` | 151.738 |
| `acceptance-vector_stabilizers-7` | `test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[4]` | 749.733 |
| `acceptance-vector_stabilizers-8` | `test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[8]` | 487.086 |
| `acceptance-vector_stabilizers-9` | `test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[6]` | 502.181 |
| `acceptance-vector_stabilizers-10` | `test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[0]` | 180.472 |
| `acceptance-vector_stabilizers-11` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_8]` | 249.798 |
| `acceptance-vector_stabilizers-12` | `test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[2]` | 912.813 |
| `acceptance-vector_stabilizers-13` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_3]` | 164.087 |
| `acceptance-vector_stabilizers-14` | `test_equivalence_witness_carries_a_representative_to_its_reflection[r_6]` | 120.441 |
| `acceptance-vector_stabilizers-15` | `test_divisibility_four_isotropic_stabilizer_constructs_exact_containing_bound` | 4.547 |
| `acceptance-vector_stabilizers-16` | `test_vector_stabilizer_has_oscars_order` | 87.958 |
| `acceptance-vector_stabilizers-17` | `test_isotropic_plane_and_flag_automorphism_generators_preserve_their_marks` | 8.496 |
