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

Investigate the remaining reflective failures in Acceptance
[40310936005](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40310936005),
port `6618e63ce1decc6ccdb42d0ddaa9e33836b1bdaf`, research
`53c48e0ab7c8f74f9c15f2e0970b11f29376dc23`. Inspect each completed
shard's JUnit and revision metadata. Every isometry, perfect-domain and
vector-stabilizer shard has a complete passing JUnit report at that pair.
The reflected-vector r_0 and degree-18 orbit consumers pass there.
Remaining families are open. Reflective shards 2, 4 and 7 contain assertion
failures and then exit 137 before writing JUnit; orthogonal-group shard 0
also exits 137. These runs are incomplete evidence.

The first failure in reflective shard 7 maps to ListReflect_549. Its
Vinberg coset normalization used QQ division in place of the pinned integer
QuoInt. The integer-division repair passes that unchanged root-count case.
Use subsequent exact failures to determine which other branches need repair.

The I_1_2 isotropic-stabilizer assertion remains open. A captured perfect-cell
pair has a verified integral witness from the existing research finite-family
transporter. The current port basis-image search exceeds 180 seconds on that
same pair; research returns in 87 seconds. The frozen owner split and pinned
Lorentzian.g:660-696 prescribe the generic configuration transporter.
The diagnostic data and exact certificate are in
`.tmp/probe/configuration-transporter/`; these are optional local evidence,
not required checkout inputs. The full I_1_2 consumer has no verified result.
The unchanged Enriques case_1 boundary consumer is also still unresolved.

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

1. Repair the remaining reflective cases from exact tracebacks. Recover
   missing reports by running selected individual cases.
2. Complete the I_1_2 and Enriques case_1 public consumers; use the existing
   research configuration operation where the pinned algorithm calls it.
3. Repair every remaining group, boundary, orbit and centralizer failure.
4. Obtain full frozen Acceptance evidence at one published port/provider pair.

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
