# Run timings

Every local test run and CI run, with its duration and outcome. New local runs are appended
automatically to `.tmp/run-timings.tsv` by `~/.cache/sage-docker-tmp/runtest.sage`.
Machine: 10 cores, 7 GB RAM, Sage in the `sage-research` container.

## 2026-11-11 to 2026-11-12

| Test or run | Code | Duration | Outcome |
|---|---|---|---|
| `test_port_negative_isometry_index_two` (diag(1,-2,-64,128), #68) | port 7882647 | 10 h 15 min, not finished | Killed. Stack stayed in `EichlerEnvelope.right_coset_transversal` (eichler.py:514) the whole time: the orbit of phi(L) on integral structures is too large for det 16384 |
| `test_congruence_isotropic_planes` (#44) | sage-extensions before 8a26ca5 | more than 1 h 45 min, not finished | Killed. Pairwise coset membership in `_orbit_words` |
| `test_congruence_isotropic_planes` | sage-extensions 9ca46e7 | 1 min 37 s | 2 failed (CON-17 regression from the oriented-flag Mor change) |
| `test_congruence_isotropic_planes` | sage-extensions 4ef62d6 | 1 h 07 min | 2 passed |
| `test_congruence_isotropic_planes` | sage-extensions b9bdd1e, publish during the run | 37 min | 2 failed with `tokenize.TokenError` (source changed under the process) |
| `test_congruence_isotropic_planes` | sage-extensions 27be129 | 59 min 39 s | 2 passed |
| `test_isotropic_flag_graded_action` Gamma case | sage-extensions before 8a26ca5 | more than 40 min, not finished | Killed |
| `test_isotropic_flag_graded_action` | sage-extensions 8a26ca5 and later | 2 min 14 s to 3 min 31 s | 3 passed |
| `test_ternary_orthogonal_number_field` (v2) | port 1a65bd6 | more than 1 h, not finished | Killed. p^8 candidate search in `maximal_orders_containing` |
| `test_ternary_orthogonal_number_field` (v2) | port 699e64a | 1 min 46 s | Failed at the `_neighbours` count assertion |
| `test_ternary_orthogonal_number_field` (v2) | port e3c34b7 | more than 35 min, not finished | Killed. Sage ideal membership in the zero-divisor loop |
| `test_mixed_abelian_normalizer` (v2) | port before 78f996f | more than 1 h 50 min, not finished | Killed. In-order walk of {0..r}^n in `rational_intertwiner` |
| `test_mixed_abelian_normalizer` (v2) | port 78f996f | 2 min 07 s | 1 passed |
| `test_factorwise_dedekind_transporter` (v2) | port before 78f996f | more than 18 min, not finished | Killed, same intertwiner walk |
| `test_skew_hermitian_quaternion_transporters` (v2) | port 87adc4d | more than 15 min | Moved to CI witnesses |
| `test_unipotent_transporter_beyond_the_linear_step` (v2) | port 63d27f6 | killed at the 7 GB limit | Moved to CI witnesses |
| `test_definite_hermitian_quaternion_transporters` (v2) | port | 1 h 57 min (7043 s) | passed; moved to CI witnesses |
| `test_arithmetic_normalizer_of_the_full_orthogonal_group` (v2) | port 2a98dbf | 44 min (2634 s) | passed; moved to CI witnesses |
| witness `test_gamma_three_splits_null_lines_..._sixteen_orbits` (v2) | port | more than 5 h, not finished | Killed |
| `test_commutant_groups` full (v2) | port be6354f | 2 h 04 min (7426 s) | 29 passed |
| `test_li_margulis_box` | port | 33 min (1982 s) | 16 passed |
| `test_orbit_counts_basis_change` (#18) | port 6f99942 | 14 min 33 s | 2 passed |
| `test_enriques_line_stabilizers` (#70) | port fda734f | 21 min 12 s | 2 passed |
| `test_dawes_2U_minus2_minus6_stabilizers` (#18) | sage-extensions 64f7a0b | 4 min 25 s | 3 passed (failed in 4 min 13 s before 84d54f0) |
| `test_general_isotropic_markings` | sage-extensions 4ef62d6 and later | 11 min 04 s to 11 min 35 s | 23 passed |
| `test_partial_isotropic_flag_orbits` | sage-extensions | 9 min 27 s to 9 min 59 s | 5 passed |
| `test_oriented_isotropic_flags` | sage-extensions b9bdd1e + 1b617bc | 6 min 49 s | 5 passed, 1 test-literal failure, fixed |
| `test_generator_determinism` (#72) | port c6db72b | 2 min 17 s | 1 passed |
| determinism + pinned Lorentzian + facet contracts | port c6db72b | 2 min 48 s | 75 passed |
| `test_edgewalk_options` (#15) | port fd8b122 | 1 min 21 s | 6 passed |
| `test_unsupported_quotient_group` (#44) | port 59a4a8e | 1 min 04 s | 1 passed |
| `test_flag_invariant_consumer` (#18) | port 79f7e89 | 46 s | 2 passed |
| `test_spacelike_vector_classes_index_one` (#18) | port ba7abdf | 40 s | 1 passed |
| `test_marked_subgroup_transporter` | sage-extensions adc6814 | 31 s | 4 passed (4 failed in 27 s at 9ca46e7 and 51c8167) |
| `test_based_isotropic_orbits` | sage-extensions | 1 min 01 s to 1 min 05 s | 5 passed |
| Acceptance sample, CI run 43114552529 | port b5e5f54, sage-extensions 27be129 | dispatched 2026-11-11 about 19:55; still running 2026-11-12 10:06 | 67 success, 18 failure, 7 cancelled, 14 queued |
| GitHub long-witness jobs, runs at f03ef36 and later | port | 6 h each | cancelled at the GitHub per-job limit |

## Acceptance run 43114552529 failures (port b5e5f54, sage-extensions 27be129)

| Shards | Error | Cause |
|---|---|---|
| boundaries 1, 2, 3, 11; isotropic_orbits 2, 3, 6 (Dawes U ⊕ U+<-2>+<-6> buildings, O_plus_N two lines and two planes) | `the plane's basis, its completion to P^perp and that to L must give a ZZ-basis of L` (determinant -3, -4) | Upstream 7b72780 (#326) completed the plane basis greedily by rank, which does not give a ZZ-basis; replaced by Smith completion of primitive sublattices |
| vector_orbits 15-23; vector_stabilizers 0, 3 (E10 pointed classes, degree 14 and similar) | `PariError: the PARI stack overflows (current size: 1073741824)` in `QuadraticForm.number_of_automorphisms` of a definite lattice, reached from `Isom(L, M).is_empty()` through `lattice_morphisms.py:3357` and `:4246` | The PARI `qfauto` stack limit of the CI job; not a port defect |
