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

Finish the unchanged frozen vector consumers with research
`53c48e0ab7c8f74f9c15f2e0970b11f29376dc23` and port source `5024411`.
The repaired definite witness backend passes exact change-of-basis, inverse
and pairing checks for both signs of E8 + <18>, and rejects the independent
non-isometric binary pair. `.tmp/repair-definite-witness-fixed.log` is the
current selected run, with reflected-vector `r_0` followed by degree-18
vector representatives. Both previously exhausted PARI qfisom at 1 GiB;
the new backend uses OSCAR invariant-sublattice decomposition through the
existing adapter. Read actual outcomes before further changes.

The I_1_2 isotropic-stabilizer assertion remains open. Its instrumented
single-case run is `.tmp/repair-isotropic-stabilizer-stack.log`, using the
same source implementation. Repeated samples locate its work in
`LorentzianPerfectLocalBackend.cell_transporter`, matching perfect cells
while constructing the ambient orthogonal group. No new assertion failure
has been observed. The earlier uninstrumented run was interrupted and is
not passing evidence. `.tmp/repair-isotropic-stabilizer-diagnosis.md`
records the diagnostic boundary.

The corrected Lorentzian flip preserves GAP's IsSubset(A,B) convention.
Both source exits return the exact supporting face x=1 of x^2-y^2 in the
independent owner example. Research's OSCAR orthogonal-group adapter has
verified order 1393459200 and exact generator actions on E8 + <18>.

The isolated provider worktree is `.tmp/research-acceptance-repair`, branch
`fix/acceptance-quotient-cardinality`. Its repairs through 2e204df5 are also
on research main at `fdfed1d61c4b2530abb9822c93e457c1d26bc825`.
Preserve unrelated work in both checkouts. Polarized extension groups must
restrict by O(L^f,h), not the full invariant orthogonal group.

Acceptance still required, in order:

1. Finish the selected vector consumers above and the unchanged isotropic
   stabilizer case `I_1_2` in `.tmp/repair-isotropic-stabilizer-stack.log`. The original
   U_2U_A3 isometry now passes its unchanged witness assertion.
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
| isometry | ci_indef_U_2U_A3 (change of basis and witness) |
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
