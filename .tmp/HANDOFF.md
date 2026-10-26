# Acceptance repair handoff

## Task

Repair the failed mathematical computations across the frozen acceptance
families, in the port or in research according to operation ownership.
The scope is isometry, boundaries, centralizers, isotropic orbits,
orthogonal groups, reflective lattices, vector orbits and vector stabilizers.
Preserve perfect-domain correctness and every frozen fixture and assertion.

## Completion condition

Every reported failing case must compute its required result. Re-run cases
whose runners were lost and cases lacking reports. Final evidence must identify
the same published port and research revisions across all families. Repairs to
shared research operations require independent mathematical checks at their
owners. A passing helper or selected example leaves the other cases open.

The reference failure set is Acceptance
[40205836899](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40205836899),
port `5c990fd35f20f194afdcd95e0bb3647118298645`, research
`f44e76f18e2e9cec43b09a846e2e74bca5fd50f6`.
Download its `acceptance-summary` and individual shard artifacts for exact case
IDs and tracebacks. Missing JUnit reports require the job logs as well.
The repair scope remains open for every listed failing family.

## Now

Verify the corrected Lorentzian flip at port `5024411` with research
`c55fa8a7d5818ff075a7146e394c70a50368d580`. The frozen U_2U_A3
isometry and Dawes discriminant-image cases both reached zero denominators
because the port reversed GAP's IsSubset(A,B) convention in every flip exit.
The small supporting face x=1 of x^2-y^2 now returns exactly
(1,-1), (1,0), (1,1). The frozen vector-orbit consumer also reached the old
second-scheme loop. Its result remains open.

The current selected frozen run is `.tmp/repair-flip-frozen-consumers.log`,
using immutable port sources `.tmp/port-5024411/src`. It exercises, in order,
the Dawes discriminant-image case, original U_2U_A3 isometry,
`test_published_representatives_are_one_per_orbit[18]`, and
`test_equivalence_witness_carries_a_representative_to_its_reflection[r_0]`.
Read its actual outcomes before selecting further repairs.

Research's exact OSCAR orthogonal-group adapter has independently verified
order 1393459200 and generator actions on E8 + <18>. The definite isometry
witness operation still uses PARI qfisom; inspect an actual failure before
changing that operation.

The isolated provider worktree is `.tmp/research-acceptance-repair`, branch
`fix/acceptance-quotient-cardinality`. Its repairs through c55fa8a7 are also
on research main at `f98317ee80577ce0ca20a529db2e40818396924a`.
Preserve unrelated work in both checkouts. Polarized extension groups must
restrict by O(L^f,h), not the full invariant orthogonal group.

Acceptance still required, in order:

1. Finish the selected frozen consumers above. The original isometry
   completes its finite lattice orbit after the quotient repair; its next
   observed failure was the corrected flip containment. Earlier interrupted
   runs and timeout traces are not passing evidence.
2. Complete the isotropic consumers after exact integer configuration
   pairings. `IsotropicCases_813` has not returned a verified result.
3. Repair all remaining reflective, group, boundary and vector cases using
   exact shard tracebacks, including cases with missing JUnit reports.
4. Obtain full frozen Acceptance evidence at one published port/provider
   pair. Run 40299696962 was cancelled and supplies no further evidence.

The following unchanged frozen cases pass locally with port source
`5024411d3fde70e6e475aa2e1e59b334772d61b5` and provider
`c55fa8a7d5818ff075a7146e394c70a50368d580`:

| Family | Case |
| --- | --- |
| boundaries | test_dawes_last_index_in_the_chain_is_the_discriminant_image_of_o_plus |
| isotropic_orbits | IsotropicCases_813 |
| orthogonal_group | test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q |
| reflective | ListReflect_590 |
| centralizers | oscar_lattices_with_isometry_110; oscar_lattices_with_isometry_118 |

These are selected local results, not hosted full-family acceptance.
Research's independent integer residue-cube, polarized U+U, stable-group
axis-stabilizer and E8 + <18> group examples pass. The full repair condition
remains open.

Other work is advancing the same port branch, including positive-vector
validation and bounded CVP. Inspect current branch tips before changing
those paths; retain their changes. Post-v1 requirements remain in TODO.md
and issues #40/#41.
