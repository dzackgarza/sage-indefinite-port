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

Repair the definite orthogonal-group computation in research
`LatticeIsometryMor._engine_group`. The frozen consumer
`test_published_representatives_are_one_per_orbit[18]` reaches
`IntegralLattice(gram).orthogonal_group()` on a rank-nine positive definite
form with an orthogonal norm-18 direction and an E8 component. PARI `qfauto`
exhausts its stack while enumerating short vectors. Preserve the complete
orthogonal group and exact live isometries; increasing a resource limit is
not evidence that the required vector computations work.

Use research branch `fix/acceptance-quotient-cardinality`, published revision
`ee96afcf95ccc3e26f7b45cf92d513dc8ae7a061`, for isolated repairs. Its worktree
is `.tmp/research-acceptance-repair`; a cloud checkout can fetch the branch.
It extends the pinned provider with integer quotient cardinality, polarized
full-glue restriction by O(L^f,h), and subgroup stabilizers. Integrate its
new changes into research main after reconciling the active main work;
preserve all unrelated dirty files there. The cardinality change already
exists on main as `1d5d0bb24d09d12c4967a31e97b52db104734902`.
Do not replace the polarized restriction by the unpolarized extension group.

Acceptance still required, in order:

1. Verify the newly pinned provider in hosted Acceptance, including the
   original isometry case, all centralizers and the boundary consumers.
   The isometry case now reaches the finite integral-representation orbit;
   its local run was interrupted, so its result is unresolved.
2. Complete the isotropic consumers after exact integer configuration
   pairings. `IsotropicCases_813` has not returned a verified result.
3. Repair all remaining reflective, group, boundary and vector cases using
   exact shard tracebacks, including cases with missing JUnit reports.

At port `36bf6d5c5ce3731179b12b0f86d732fb5dc0f739` with the provider above,
the unchanged reflective case `ListReflect_590` passes. Research's independent
integer residue-cube, polarized U+U, and stable-group axis-stabilizer examples
pass. These examples leave the full repair condition open.
Fetch branch tips and completed Acceptance evidence before selecting a new
failure. Post-v1 requirements remain in TODO.md and issues #40/#41.
