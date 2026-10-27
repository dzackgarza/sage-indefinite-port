# Acceptance repair handoff

## Ambient task

Repair the actual failed mathematical computations in the frozen Acceptance
suite, in sage-indefinite-port or its research provider according to operation
ownership. The required scope is isometry, boundaries, centralizers, isotropic
orbits, orthogonal groups, reflective lattices, vector orbits and vector
stabilizers. Preserve perfect-domain correctness and every frozen assertion.

The repair is incomplete. Work is paused at the user's request. Resume only
when requested. This document replaces the previous handoff.

## Completion condition and verified status

Every reported failing case must compute its required result. Cases with lost
runners or missing reports require execution. Final evidence must identify one
published port/provider pair across all families. Shared provider repairs also
need independent mathematical assertions at their owners. No additional stubs
or weakened expectations are authorized. A helper pass or a selected case does
not discharge its family.

The original failure set is Acceptance
[40205836899](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40205836899),
port `5c990fd35f20f194afdcd95e0bb3647118298645`, provider
`f44e76f18e2e9cec43b09a846e2e74bca5fd50f6`.

The latest inspected full run is
[40339136782](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40339136782),
port `b3f01e016440b405081a61a8517be536475f8ec1`, provider
`07b24a2d3c02921d2a6098c32be0a1a171ce6406`. Its completed JUnit reports establish:

| Required family | Evidence at this pair | Remaining acceptance |
| --- | --- | --- |
| perfect_domains | All 8 shards pass | Repeat at final pair |
| isometry | All 8 shards pass | Repeat at final pair |
| vector_stabilizers | All 6 shards pass | Repeat at final pair |
| centralizers | Entire shard, 5 assertions, passes | Repeat at final pair |
| isotropic_orbits | Shards 8 and 9 pass | Other shards and final pair |
| vector_orbits | Shard 0 passes | Other shards and final pair |
| reflective | Killed shards lack complete reports | All shards |
| boundaries | Full-family evidence open | All shards |
| orthogonal_group | Full-family evidence open | All shards |

Job status alone is insufficient. Download JUnit and both revision files.
Local copies are under `.tmp/repair-40339136782/`; `verified-junit.json`
predates the centralizer download and is not an authoritative complete index.

The current published port is `199168dad31efc17981411d04579a8cf68a1089d`,
branch `implement/lorentzian-completion`. Its justfile pins provider
`53d0db2e7d0c2ddb96617cfb361a43c18235d1d2`. Acceptance
[40354157956](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/40354157956)
was pending at the pause; 40339136782 was still in progress. Refresh those
results before selecting failures. Subsequent handoff-only revisions do not
change this code/provider pair.

## Next work, ranked

### Now: remove the repeated preamble operations in the discriminant transport

The Enriques route is now identified and `_gamma_h` for case_1 computes its
exact result locally in 2609 s. It does not use the generic predicate coset
path, so the finite-group transversal repair above is not what this case
needs. The stages of `.tmp/probe_enriques_stages.py`, provider `53d0db2e7d`:

| stage | seconds |
| --- | --- |
| setup: lattice, transport isometry, polarization | 147 |
| `invariant.O().stabilizer(h).image_under(discriminant_representation)` | 1356 |
| `_coinvariant_glue_inverse_maps` | 332 |
| 19 `_coinvariant_discriminant_from_invariant` transports | 400 |
| `coinvariant_extension_subgroup` (`discriminant_preimage`) | 382 |

The stabilizer image is the port's pointed Aut((M,h)) (`.tmp/probe-enriques-pointed.log`):
O(h^perp), a definite rank-9 lattice with 21 generators, 732 s; the
`_lattice_action` stable bound 92 s; the residual stabilizer 86 s. A warm
definite O(9) call costs 0.1 s (`.tmp/probe_isometry_group_sign.py`), so most
of the 732 s is the first Julia engine call of the process, not the group.

The transports are dominated by per-element preamble work, not by any group
algorithm (`.tmp/probe-enriques-profile.log`: 20 conversions, 2081 s under
cProfile, 10.08M `category.__classcall__`, 20.7M set-element constructions,
9.2M `__contains__`, 2.0M `_own_ring`, 84k `abelian_group_gps.exponents`).
That is roughly a hundred thousand preamble operations per converted
generator for a 12-generator, order-1024 discriminant group.

Next specimen: one conversion of one invariant generator, timed before and
after the change. The unit is filed as
[research #621](https://github.com/dzackgarza/research/issues/621), whose body
carries the same profile and the same specimen. It is not the port's to
repair: `_presented_module_from_morphism`, `Mor.Of`,
`_normalized_form_automorphism` and `_coinvariant_discriminant_from_invariant`
are provider owners. The port now pins provider
`d20472d309dd64d168a6aacd71f074570b81538`; the timings above were taken at
`53d0db2e7d`. The same overhead multiplies the 2395-vector
`IsotropicCases_8716` and the 49,861-vector A7 configuration, so measure the
change on that conversion before rerunning either family. The case also still
needs `_counts(group)`; no case has completed the unchanged assertion.

The same route runs on a rank-4 specimen in seconds
(`tests/lattices/test_polarized_equivariant_lattice.sage`: `U + U`, the
involution exchanging the summands, `h = e + u + f + v`). There
`.tmp/research-acceptance-repair/.tmp/probe_transport_profile.py` profiles two
conversions there: 3.75 s in `_coinvariant_discriminant_from_invariant`, 2.50 s
of it inside `TorsionFormOrthogonalGroup.from_morphism`, and 1.07 s in the
caller's `morphism.inverse()`. `normalization_isometry` dominates the first
conversion only (one cached call); the marginal cost is the presentation work
each conversion repeats: `_normalized_form_automorphism` rebuilds the engine
matrix by lifting every generator image through the framing,
`_presented_module_from_morphism` and `Mor.Of` build a presentation per
morphism, and `_from_engine` transports the automorphism and its inverse into
two more presented morphisms. Measure on this specimen; do not rerun the
43-minute Enriques probe.

The owner checks already passing on the stable provider are nonnormal S3
cosets, the C4 preimage, the anisotropic F3 plane and the split four-dimensional
F2 form. They verify actions, sections, stabilizer orders and Schreier data.
Their sources are `tests/groups/test_finite_group_subgroups_and_quotients.sage`
and `tests/lattices/test_elementary_two_discriminant_groups.sage`.
The older protected `test_schreier_generating_morphisms.sage` uses a retired
spelling; preserve its expectation and resolve any bridge at its owner.

Provider `53d0db2e7d` repairs an unconditional
f-string in `BiproductModules.from_coproduct_cocone`. The owner biproduct
universal-map identity assertion passes. That repair is not full Enriques
acceptance. Earlier local I_1_2 stabilizer-index and D-tower N_6 boundary
assertions pass; their final-pair hosted confirmation remains required.

### Then: complete IsotropicCases_8716 and the killed reflective shards

The current provider uses the owned lattice Gram tensor and family map at the
private graph adapter, lowers their matrix pullback there, shares exact edge
labels and releases temporary family realizations. Bliss supplies index
certificates; their composite is checked against every exact colour, pairing
and relation table. Canonical graph copies are released between the two
certificate computations. These contracts must remain intact.

Unchanged `IsotropicCases_8716` is still unverified. The latest isolated run
at provider `daf20e327f` reached multiple neighbour transporters in the first
perfect cell, then was stopped without a result. Use
`.tmp/repair-isotropic-8716-certificates.log` and
`.tmp/repair-isotropic-8716-current-stack.json` only as diagnostic evidence.
The initial configuration has 2395 vectors. The earlier capped probes are
not hosted-memory certificates: virtual reservations exceeded 6 GiB while
resident memory and swap were much smaller. The latest run's artificial
address-space limit was removed. Do not reinstate that limit as a correctness
gate or treat interrupted argument rendering as a test result.

At the older hosted pair, killed-shard progress logs locate these candidate
active cases: isotropic shards 2/7/10/0 at 8716/8806/8714/8666; reflective
shards 13/4/5 at 7551/7799/7528. They are locators, not completed failures.
Additional killed shards are recorded in `.tmp/repair-current-jobs.json` and
`.tmp/repair-*-hosted.log`; refresh the run before using that snapshot.
The inspected killed logs contain exit 137 without an assertion traceback.

On provider `53d0db2e7d`, unchanged ListReflect_7551 and ListReflect_7799
pass in `.tmp/repair-reflective-hosted-locators.log`. They do not establish
full-shard memory behaviour. The negative A2 pairing obstructions, positive
redundant-family transporter, ListReflect_841 and ListReflect_964 also pass.
Read the final-pair hosted results to determine whether retained graph data
were the remaining full-shard defect.

### Then: finish the A7 discriminant-image computation

The unchanged
`test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q` is still open.
Its initial perfect configuration has 49,861 distinct vectors. Streamed
pairing/hash calculations reach exact cone conversion, but neither the full
PPL consumer nor the separate Normaliz cone probe has a result.

- Full consumer log: `.tmp/repair-orthogonal-a7-flint.log`, provider 07b24a2.
- Exact Normaliz probe: `.tmp/repair_a7_normaliz.py` and
  `.tmp/repair-a7-normaliz.log`.
- Reusable local initial rows: `.tmp/repair-a7-initial-rows.sobj`.

Both computations were stopped. No backend change was made. Sage's installed
Normaliz backend accepts the same rational V-presentation; see its
[official documentation](https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/backend_normaliz.html).
Do not infer that it solves this case from availability alone.

If complete facet construction remains the obstacle, inspect the actual
pinned `lorentzian_perfect.h::DualDescriptionRecordFullDim` call and
`POLY_RecursiveDualDesc.h` before changing the specialized facet-orbit route.
The current `facet_orbit_permutations` materializes all facets and full ray
permutation tuples. Assess its required memory before restarting another
large cone computation. This is an unresolved implementation question, not
a selected replacement algorithm. Preserve exact facets, orbit coverage,
side pairings and the owned cone interfaces.

### Finally: discharge all remaining families at one pair

Read every completed artifact from the refreshed full run. Route each actual
failure to the operation owner, repair it and execute its unchanged consumer.
Keep missing/interrupted shards open. The full repair condition above remains
a requirement even when selected mathematical examples pass.

As of 2026-10-27 10:20 UTC every Acceptance run after 40339136782 is
`cancelled`: successive pushes to `implement/lorentzian-completion` supersede
each other, so no full suite has completed for over fourteen hours. A complete
evidence run requires a quiet branch; read the run of the revision you cite,
not the newest one.

## Checkout and ownership constraints

The stable provider worktree is `.tmp/research-acceptance-repair`, branch
`fix/acceptance-quotient-cardinality`. Its published tip is the pinned
`53d0db2e7d0c2ddb96617cfb361a43c18235d1d2`. Its tracked tree is clean.

The replacement library checkout is `.tmp/sage-extensions-acceptance-repair`.
The repairs are carried to its published main at
`74d34ba59d6581969101556dcef9e8508af10aa2`; its tracked tree is clean.
Do not switch the port provider merely because those changes are carried.
The library still needs verification of presentation-square construction,
lattice-isometry inverse construction and the selected free basis in group
resolutions, tracked in
[sage-extensions #1](https://github.com/dzackgarza/sage-extensions/issues/1).

Research main at `/home/dzack/research` has unrelated migration conflicts;
preserve them. Preserve the pre-existing port paths `.agents`, `.hermes`,
`=2.18`, `cysignals_crash_logs/` and `tests/.report/`.

Generic arithmetic, finite configurations, group sections and gluing belong
to the provider. Indefinite search, perfect-domain traversal, Allcock cycles,
parabolic compatibility-image discovery and specialized group assembly belong
to this port. Polarized extension groups must restrict by O(L^f,h).
Use the pinned upstream routines and their actual callees. Preserve the
post-v1 requirements in TODO.md and issues #40/#41.

## Reproduction constraints

The local Sage wrapper uses container `sage-dev-develop`. A targeted command
from the port root uses:

```sh
direnv exec . env \
  SAGE_JULIA_COMMAND=/home/dzack/.julia/juliaup/julia-1.12.5+0.x64.linux.gnu/bin/julia \
  PYTHONPATH="$PWD/.tmp/research-acceptance-repair/src:$PWD/src:$PWD:/home/dzack/ai-review-ci/tool-artifacts/pytest_plugins" \
  sage -python -m pytest -p qc_sage_session 'tests/acceptance/MODULE.sage::CASE' -q --tb=short
```

Run only selected cases locally. For combined owner and consumer checks,
use `-c pyproject.toml` from the port root. Local probes and saved rows are
optional aids, not cloud-checkout prerequisites. The four local long-running
repair computations were terminated at the pause; none supplies a pass.
Hosted Actions may still run independently. Refresh their state on resume.
