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

Complete the frozen vector consumers using research revision
`c55fa8a7d5818ff075a7146e394c70a50368d580`. Its definite orthogonal groups
use the existing OSCAR adapter. The independent lattice E8 + <18> has
verified group order 1393459200 and exact generator actions. The unchanged
`test_published_representatives_are_one_per_orbit[18]` and
`test_equivalence_witness_carries_a_representative_to_its_reflection[r_0]`
still need results. The definite isometry witness operation still uses PARI
qfisom; inspect an actual failure before changing that operation.

The isolated provider worktree is `.tmp/research-acceptance-repair`, branch
`fix/acceptance-quotient-cardinality`. Its repairs through c55fa8a7 are also
on research main at `f98317ee80577ce0ca20a529db2e40818396924a`.
Preserve unrelated work in both checkouts. Polarized extension groups must
restrict by O(L^f,h), not the full invariant orthogonal group.

Acceptance still required, in order:

1. Resolve the original isometry case
   `test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U_A3]`.
   The quotient cardinality assertion is repaired, but a completed result
   remains missing. A run with provider ee96afcf produced an incomplete
   timeout traceback in `_ExactCVPEngine._filtered_candidates`.
   `.tmp/repair-isometry-diagnosis.md` records competing explanations;
   `.tmp/repair-isometry-boundaries.log` records a run without that timeout.
   A partial timeout trace does not identify the currently executing operation.
2. Complete the isotropic consumers after exact integer configuration
   pairings. `IsotropicCases_813` has not returned a verified result.
3. Repair all remaining reflective, group, boundary and vector cases using
   exact shard tracebacks, including cases with missing JUnit reports.
4. Obtain full frozen Acceptance evidence at one published port/provider
   pair. Run 40299696962 was cancelled and supplies no further evidence.

With port 36bf6d5/d0a9b3d and provider ee96afcf, unchanged cases
`ListReflect_590` and `oscar_lattices_with_isometry_118` pass. Independent
research examples for integer residue cubes, polarized U+U, and stable-group
axis stabilizers pass. The full repair condition remains open.

Other work is advancing the same port branch, including positive-vector
validation and bounded CVP. Inspect current branch tips before changing
those paths; retain their changes. Post-v1 requirements remain in TODO.md
and issues #40/#41.
