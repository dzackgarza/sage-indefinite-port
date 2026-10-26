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

Resolve research finite presented-module cardinality through the existing
invariant-factor and quotient-ring computations. The isometry consumer
`test_change_of_basis_gives_an_isometric_lattice[ci_indef_U_2U_A3]`
reaches this operation through `IntegralStructureAction.finite_representation`.
The same error occurs in group, boundary, isotropic and stabilizer cases.
Research commit `398887702a` changes quotient cardinality; verify its effect.
Current research also imports `sage.categories.ordered_fields` from
`categories/rings/commutative_algebra.py`; the local Sage image cannot import
that module. Repair this integration defect before evaluating the current
provider and before changing the provider pin.

Next, classify and repair the remaining tracebacks: polarized equivariant
lattice extension subgroups, centralizers, isotropic and reflective branches,
and PARI failures in vector computations. Follow each actual source call;
retain the pinned algorithms and the public mathematical objects.

Fetch current branch tips before editing. The inspected port tip is
`4b066e3f973cc80dad5f774a4b69e14d445cfd29`; its Acceptance run
`40294291516` was pending. Research main has unrelated dirty scheme source;
preserve it. Older MVP statements and helper-only evidence do not discharge
this repair task. Post-v1 requirements remain in TODO.md and issues #40/#41.
