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
The full repair condition remains open.

## Now

Finish the real consumers with research
`260c6a92c192d879d8a7cf9c988eeef97d4de84b`. The pinned provider now uses
OSCAR normal forms for finite-form orthogonal groups and Bliss certificates
for finite-configuration transport. The split quadratic forms in dimensions
2, 4 and 10 and the exponent-four hyperbolic form pass their independent
order and generator assertions. Finish the selected bilinear U(2) check in
`.tmp/repair-two-discriminant-fixed.log`; a passing quadratic check does not
establish the bilinear value-modulus branch.

Next work, in order:

1. Finish the unchanged `ci_indef_U_U_E7` isometry assertion in
   `.tmp/repair-isometry-UUE7-bliss.log`. It fails with MemoryError in the
   old Sage graph-isomorphism route in Acceptance
   [40320087648](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40320087648)
   (port `7e56222`, provider `53c48e0`). Bliss passes the reindexed A2 family
   and the independent U + E8(-2) full-configuration witness checks.
2. Re-run unchanged I_1_2 isotropic-stabilizer and Enriques case_1 boundary
   consumers. Both reach Sage's exhaustive finite-form orthogonal-group
   search with the old provider. Their prior runs were stopped after that
   call was identified. The repaired dimension-ten operation now passes its
   mathematical assertions; the full consumers are still unverified.
3. Finish the A7 + diag(1,1,-1) discriminant-image assertion. Its initial
   perfect configuration has 49,861 vectors; the old pairing block retained
   51,057,664 exact scalars, and the hosted shard exited 137. Pairings now
   stream by row, and the hash counts the pinned upper-triangular multiset
   directly. The exact x=1 face of x^2-y^2 and a large integral basis change
   verify the invariant. The full A7 consumer is still open.
4. Read all completed artifacts from the newest Acceptance run and repair
   the remaining reflective, boundary, group, orbit and centralizer cases.
   Check JUnit and revision metadata, not job status alone. The Vinberg
   integer-quotient correction passes the unchanged ListReflect_549, 231,
   516, 241 and 880 cases. Other reflective cases remain open. Immediate
   failure reporting now preserves tracebacks before a killed shard exits.
5. Obtain full frozen Acceptance evidence at one published port/provider
   pair, and carry any further provider repairs to the current library owner.

Every isometry, perfect-domain and vector-stabilizer shard has a complete
passing JUnit report at port `6618e63` / research `53c48e0` in
[40310936005](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40310936005).
That does not validate subsequent changes. The reflected-vector r_0 and
degree-18 orbit consumers also pass in that run.

## Provider and scope boundaries

The isolated provider checkout is `.tmp/research-acceptance-repair`, branch
`fix/acceptance-quotient-cardinality`. Repairs through 2e204df5 were merged
into research main at `fdfed1d61c4b2530abb9822c93e457c1d26bc825`.
The newer finite-form and Bliss repairs are published on the provider branch.
Research main is moving library source to
[sage-extensions](https://github.com/dzackgarza/sage-extensions); its checkout
has unrelated merge conflicts. Preserve that work. Verify the replacement
library source before integrating provider changes there.

Polarized extension groups must restrict by O(L^f,h). Generic arithmetic,
finite configurations and finite-form groups belong to the library; the
indefinite searches and specialized group assembly belong to this port.
Preserve unrelated port work, including positive-vector validation and CVP.
Post-v1 requirements remain in TODO.md and issues #40/#41.
