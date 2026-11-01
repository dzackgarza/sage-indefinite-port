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

The local Sage wrapper uses container `sage-research`. A targeted command
from the port root uses:

```sh
direnv exec . env \
  PYTHONPATH="$PWD/.tmp/research-acceptance-repair/src:$PWD/src:$PWD:/home/dzack/ai-review-ci/tool-artifacts/pytest_plugins" \
  sage-docker-python -m pytest -p qc_sage_session 'tests/acceptance/MODULE.sage::CASE' -q --tb=short
```

Run only selected cases locally. For combined owner and consumer checks,
use `-c pyproject.toml` from the port root. Local probes and saved rows are
optional aids, not cloud-checkout prerequisites. The four local long-running
repair computations were terminated at the pause; none supplies a pass.
Hosted Actions may still run independently. Refresh their state on resume.

## Recorded repairs after integration on main

### Frozen case made runnable: `test_zero_form_orthogonal_group_is_the_full_general_linear_group`

Amendment chronology, stated exactly:

- **2026-10-29, upstream record first.** Issue #33 comment
  [6495821128](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6495821128)
  reports this same defect — `rows_of` returns mutable list-of-lists, so set
  membership raises `TypeError` and tuple-pattern checks never match — and says
  the test canonicalizes every generator matrix to tuple-of-tuples in commit
  `d48b2c1`.
- **That commit is not in this repository.** `GET /repos/.../commits/d48b2c1`
  answers 422 "No commit found for SHA"; the SHA is in no fetched ref. The
  published line merged into main (`73ef6713`) still carried the un-normalized
  assertion, and the case failed at the CI pair.
- **2026-10-31, `b8f35104` on main.** The same normalization is applied on
  `main`: `rows = tuple(tuple(tuple(row) for row in rows_of(g)) for g in
  computed)`, matching the sibling `test_zero_form_orthogonal_group_uses_permuted_coordinate_labels`
  and `test_degenerate_orthogonal_group_*` cases.

No assertion, recorded matrix or family changed; the case previously could not
evaluate its comparison. The record above is the #33 route; this note is written
after the change and does not claim otherwise.

Evidence, port `b8f35104`, provider `53c48e0ab7` (the `pyproject` pin — a
different pair from the acceptance baseline below, labelled separately):
`test_orthogonal_group.sage -k "degenerate or zero_form or hyperbolic_plane or
ambiguous_binary"` → 9 passed, 38 deselected (112 s); 8 passed + 1 failed before.

### Provider contract resolved: one pin, `064cd6e8e90c`

`acceptance.yml` runs `just ci-provision-sage`, which fetched a hard-coded
research revision and installed it with `--no-deps`; `pyproject.toml` declared a
different one (`53c48e0a`), so a packaging install and CI ran different provider
APIs (36 provider files differ, the commits are unrelated lines). The CI
revision is the acceptance contract, so `pyproject.toml` now declares
`064cd6e8e90cc9cbad17ce7d9574b497c864b8c0` and the recipe reads the revision
from `pyproject.toml` (`sed -n 's|.*research@\([0-9a-f]\{40\}\).*|\1|p'`) and
fails if it is absent, so the two cannot drift again.

Baseline pair from here on: port `main` + provider `064cd6e8e90c`. Results taken
against `53c48e0ab7` (the former `pyproject` value) are labelled separately and
are not acceptance evidence.

### Resolved by #33 amendment: `test_binary_cycle_period_respects_permuted_coordinate_basis`

Reproduced at the baseline pair (port `f05c8830`, provider `064cd6e8e90c`) in
3.7 s. Two corrections, each recorded on #33 before it was applied:

1. [comment 6548844861](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6548844861)
   — orientation. Frame A `Lattices(ZZ)([[2,0],[0,-1]])`, start `(1,0)`: representatives
   `(1,1)`, `(3,4)`, q = 1, 2, period rows `[[3,4],[2,3]]`, det 1, form preserved.
   Frame B `Lattices(ZZ)([[-1,0],[0,2]], module_generators=("negative","positive"))`,
   start `(0,1)` — same lattice, same start vector, reversed coordinates:
   representatives `(-1,1)`, `(-4,3)`, period rows `[[3,-2],[-4,3]]`, det 1, form
   preserved. The case asserted equality of the swapped Frame B set with the
   Frame A set; it holds only under the reflection `(x, y) |-> (x, -y)`, measured
   as `{(1,-1),(3,-4)}` on both sides. The pinned source is equivariant, not
   invariant: `src_lorentzian/two_dim_lorentzian.h` at `a55fcb7b` defines
   `det_two(r,l) = r(0)l(1) - r(1)l(0)` (line 50) and `Promised` requires
   `det_two(r,l) == 1` in the supplied frame (lines 70-74).
2. [comment 6548850475](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6548850475)
   — owned integer at the comparison. After (1), `period.determinant() == 1`
   failed with `1 == 1`: the provider returns an owned `ZZ` element and the owned
   `ZZ` registers no coercion to or from Sage's `ZZ`, while the preparsed literal
   is `Integer(1)`. Corrected to `int(period.determinant()) == 1`, the convention
   `497fff9a` established in this file.

Evidence after both: `test_reflective.sage -k "binary_cycle or promised or
edgewalk_period"` → 4 passed, 17528 deselected (2.91 s).

What is verified: both frames produce det-1 form isometries as periods, with
primitive positive representatives inside the bound, and the reflection relation
between frames. What is not claimed: any general statement about the cycle under
arbitrary basis changes; only the measured reversal case is established.

### Shared runtime defect: canonical labeling aborts in the pinned image

Owner: the shared runtime, filed as
[ai-review-ci#490](https://github.com/dzackgarza/ai-review-ci/issues/490); the
image build's own tracker (`dzackgarza/sage`) has issues disabled. Not a port
defect and not worked around in the suite.

| item | value |
| --- | --- |
| reference | `ghcr.io/dzackgarza/sage:develop` |
| manifest digest (amd64) | `sha256:c3a9f1637d43b727ab2b24c0ccf1bc9ee513a08a02094d20265ef2d8643dd` |
| `org.opencontainers.image.revision` | `b0cd7d86e21eb450ed90e4c83707920e1d5d429a` (`/sage`) |
| Sage / Python | `11.0.beta5` / 3.14.10 |

Invocation: `/sage/.venv/bin/python repro.py` inside the image, one subprocess per
backend. Measured there, in that image only:

| backend | call | outcome |
| --- | --- | --- |
| bliss | `Graph([[0,1]]).canonical_label(algorithm="bliss")` | SIGSEGV |
| bliss | `graphs.CycleGraph(5).canonical_label(algorithm="bliss")` | SIGSEGV |
| partn_ref | `Graph([[0,1]]).canonical_label(algorithm="sage")` | ok, 1 edge |
| partn_ref + labels | `Graph([[0,1]]).canonical_label(edge_labels=True)` | SIGSEGV |
| partn_ref | `Graph([[0,1],[1,2],[2,3],[3,4],[4,0]]).canonical_label(algorithm="sage")` | ok, 5 edges |

Reproduction scope: the matched dev-image investigation reproduced the bliss
abort **natively**, directly on the image's interpreter, as well as under PRoot,
and isolated build ABI and copied-binary defects in that build. This is therefore
an image-build defect, not an artefact of the downstream runner, and the earlier
"runner blocker" framing is withdrawn.

Consequence: frozen cases that canonicalize root graphs with `edge_labels=True`
abort their shard (`test_mertens_omega_and_minus_one_generate_o_l` surfaced as
`MemoryError` in the edge-labelled path; `IsotropicCases_186` likewise). Those
families cannot be evaluated until the image is fixed; the arithmetic behind them
is untouched by this record. Port work continues independently of the image fix:
the defect is tracked upstream (ai-review-ci#490) and no port module is changed to
accommodate it.


## Library migration to sage_extensions (port-owned)

The port now imports the shared library at its current home. `dzack_research.preamble.*`
became `sage_extensions.*` in 22 files (12 under `src/`, 10 under `tests/`), with
`FramedFreeModules` and `ModulesOverIntegralDomains` taken from their defining modules
because `sage_extensions.all` does not re-export them. `pyproject.toml` declares
`sage-extensions @ ...@7b222e86214524a53802350f3cfda9c513707600` and `ci-provision-sage`
clones that repository, reads the revision from `pyproject.toml`, and installs it
`--no-deps -e`. The pre-move research library is no longer installed or used.

Verified at port `main` + sage-extensions `7b222e86` in the shared image
(`b0cd7d86...`, Sage 11.0.beta5, CPython 3.14.10):

| check | result |
| --- | --- |
| every port module imports | 18/18 |
| `U.Isom(U).an_element()`, `(U + <2>).Isom(...).an_element()` | isometry returned, determinant 1, form preserved |
| `test_vector_orbits::test_anisotropic_binary_vector_classes_use_primitive_embeddings` | pass |
| `test_isometry -k "coordinate or permuted or split"` | 4 pass, 1 library-blocked |
| `test_vector_stabilizers::...plane_and_flag...` | library-blocked |
| `test_centralizers` | 3 pass, 2 library-blocked |

Pre-migration acceptance stays on record as taken at the pre-move pair (port `04eb0934`,
provider `064cd6e8e90c`): 230 passing cases on the active hosted run.

Library blockers, all filed at
[sage-extensions#5](https://github.com/dzackgarza/sage-extensions/issues/5) and not
worked around here: `module_morphism()` raises `_form_module_morphism` on every isometry
including the identity; isometries lack `_element_function`; a degenerate Gram matrix
constructs a `_PresentedModule` without `Isom`. The port references none of those
internals and constructs only through `Lattices`, `Isom`, `Aut` and the homset
constructors.

### Verified canonical pair after the library repair

Port `main` (`7cb71cca`) + sage-extensions `7b222e86214524a53802350f3cfda9c513707600`,
image `b0cd7d86e21eb450ed90e4c83707920e1d5d429a` (Sage 11.0.beta5, CPython 3.14.10), OCI
environment preserved (`JULIA_DEPOT_PATH=/root/.julia:/opt/julia-depot:`). This is the pair the
shared image recipe should provision.

- Every former sage-extensions#5 reproduction passes there: `module_morphism()` on
  `Isom(images)`, `Aut(images)`, `Isom(dict)` and `Isom(L).identity()`; `U.O().cardinality()` = 4;
  `U.O().discriminant_image()`; `U.stable_orthogonal_group().order()` = 4; identity
  product/inverse/determinant; a degenerate Gram constructing a lattice with `Isom`
  (rank 3, signature (2,0)).
- Bounded consumers, all passing: centralizers 5/5 (78 s), plane/flag + divisibility 7/7 (68 s),
  pointed/divisibility 9/9 (68 s), isometry coordinate subset 5/5 (65 s) including the
  degenerate witness, isometry-state regressions 3/3 (2 s), reflective binary-cycle 4/4 (3 s),
  anisotropic binary vector classes 1/1 (1.5 s).
- `test_orthogonal_group` bounded subset still aborts in the image's bliss canonical-label path
  (ai-review-ci#490) — an image defect, distinct from this pair; not rerun further.
- `test_boundaries` remains open under ARITHMETIC/#18; whole-file isometry/stabilizer runs are
  long-tail and also pass through the image defect.

Pre-migration acceptance (port `04eb0934` + research `064cd6e8`) stays recorded separately and
does not certify this pair.

### A_REFLECTIVE branch sample at the canonical pair with the owner Bliss patch

Identities for this run: port `main` `52eeedfa17573823b7687039a6d27f05726bb54c`; provider
sage-extensions `7b222e86214524a53802350f3cfda9c513707600`; base image
`b0cd7d86e21eb450ed90e4c83707920e1d5d429a` (Sage 11.0.beta5, CPython 3.14.10); owner patches
`bliss-repair.patch` `baf4b826e8093c7891de3650789c7f1e13c7d9e82de362a124486956e53f32b3` and
`bliss-repair-v2.patch` `ff86f4a9a3cf02ae47ee81dc221ee048ba7f9fe277cfc3dbeb2f1ae091dbb5e1`;
installed `sage/graphs/bliss.cpython-314-x86_64-linux-gnu.so`
`d941337e3a69b8256771c2fbbff5c2604e724b717dcc41352e12a3f7f9b75a75`, recorded in
`/opt/bliss-repair.json`. This rootfs is the base image plus those owner patches; the published
image is unchanged.

`tests/acceptance/test_reflective.sage::test_edgewalk_branch_representative_sample` → **1 passed
(14.5 s)**. It exercises the starting-corner and cusp branches, `U+<46>` termination, the compact
hexagon under a permuted labelled Gram, the binary branches, and the required cycle/edge routines,
and it performs the exact labelled-graph comparison that the canonicalization repair unblocked.

This is branch coverage for #15, not A_REFLECTIVE acceptance: the file's 17,532 parametrized
cases remain long-tail and unrun, and no family is closed by this record.

### A_PERFECT per-case record (bounded, 180 s per case)

92 frozen cases of `test_perfect_domains.sage`, each run as its own invocation so a
long case cannot swallow the others. Raw table: `~/sage-image/perfect-percase.log`.
Port `52eeedfa` (main) + sage-extensions `7b222e86`, base image `b0cd7d86` with the owner
Bliss patches.

**86 cases passed** (sum of case times 1318 s):

- test_facet_orbit_permutation_view_materializes_exact_isometries — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_0] — 4s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_10] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_11] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_12] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_13] — 5s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_14] — 4s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_15] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_16] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_17] — 4s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_18] — 6s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_19] — 4s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_1] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_20] — 7s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_21] — 11s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_22] — 6s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_23] — 10s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_24] — 16s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_25] — 8s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_26] — 83s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_27] — 5s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_28] — 20s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_29] — 14s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_2] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_30] — 64s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_31] — 32s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_32] — 88s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_34] — 30s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_35] — 32s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_36] — 13s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_38] — 103s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_3] — 4s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_4] — 2s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_5] — 4s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_6] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_7] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_8] — 3s
- test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_9] — 3s
- test_isotropic_shell_survives_a_permuted_gram_basis — 2s
- test_large_cell_cone_uses_exact_normaliz_h_presentation — 12s
- test_large_pairing_graph_fallback_returns_exact_transporter — 3s
- test_lazy_neighbor_cones_retain_exact_transformed_rays — 4s
- test_non_diagonal_isotropic_shell_respects_explicit_coordinate_labels — 1s
- test_non_diagonal_isotropic_shell_uses_the_exact_sphere_center — 2s
- test_non_diagonal_positive_shell_respects_explicit_coordinate_labels — 2s
- test_non_diagonal_positive_shell_uses_the_pinned_cvp_center — 2s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_0] — 4s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_10] — 2s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_11] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_12] — 3s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_13] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_14] — 3s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_15] — 3s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_16] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_17] — 3s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_18] — 6s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_19] — 4s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_1] — 4s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_20] — 8s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_21] — 11s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_22] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_23] — 7s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_24] — 157s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_25] — 8s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_26] — 43s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_27] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_28] — 9s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_29] — 9s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_2] — 3s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_30] — 43s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_31] — 21s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_32] — 88s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_34] — 31s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_35] — 28s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_36] — 15s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_38] — 92s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_3] — 6s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_4] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_5] — 8s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_6] — 5s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_7] — 8s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_8] — 3s
- test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_9] — 3s
- test_positive_shell_scales_rational_direction_pairing_bound — 2s
- test_positive_shell_survives_a_permuted_gram_basis — 1s
- test_traversal_retains_only_adjacency_facet_orbits — 4s

**6 cases exceeded the 180 s cap** — three rank-6 fixtures in both modes:
`Result_Enumeration_33` (total, isotropic), `Result_Enumeration_37` (total, isotropic),
`Result_Enumeration_39` (total, isotropic). Exceeding the cap establishes only that a case
ran longer than 180 s; it is not evidence of nontermination.

The pinned reference (`CI_tests/28B_LorentzianPerfStabEqui/Result_Enumeration`) records
exactly the fixture counts for all three: 33 → total 18 / isotropic 15, 37 → 2 / 2,
39 → 12 / 12, so the fixtures are faithful and only the cost is open. Case 37 is the
anomaly: two perfect cells in total, yet over 180 s, so the cost is not the cell count.

## Case 37: both polyhedral engines are fast; the provider's dual-point hash collapses

Measured on port `52eeedfa` + `sage-extensions 7b222e86`, base image `b0cd7d86` with the owner
Bliss patches — the same identities as the records above. Engine timings here are engine-only: no
owned cone is constructed inside the timed call.

`Result_Enumeration_37` total mode, initial perfect cell: 202 rays, rank 6, max |coordinate| 59.

| call | time | result |
| --- | --- | --- |
| Normaliz `NmzResult(cone, "SupportHyperplanes")` | 0.113 s | 570 facets, 0 equations |
| PPL `Polyhedron(vertices=[0], rays=<202 rows>).inequalities()` | 0.052 s | 570 facets, 0 equations, 112 extreme rays |
| PPL `Polyhedron(ieqs=<570 Normaliz rows>)` | 0.09 s | 112 extreme rays |

The two 570-facet families are the same rows: identical as raw signed primitive integer rows and
identical as sign-canonical sets, 570 distinct on each side.

**The earlier attribution of case 37's cost to PPL's V-to-H conversion is withdrawn.** PPL V-to-H on
those 202 rays is 0.052 s. The 1800 s probe that did not return sat in `cell_cone`'s last step, the
provider's owned-cone construction, which was inside the timed call.

That step is `sage_extensions.categories.polyhedral_cones.RationalPolyhedralCone.__init__` calling
`finite_ordered_set` over ~570 covectors of `lattice.dual_module()`. It is superlinear because
`hash()` of a dual-module point does not depend on the point: 570 distinct covectors share one hash,
so `dict.fromkeys` degenerates to pairwise `==`, each a module-morphism `_richcmp_` costing ~17 ms.

- `finite_ordered_set` over 570 such covectors: not finished in 560 s.
- `_rational_polyhedral_cone` from the 570-facet H-presentation: not finished in 779 s.
- Isolated, on the standard lattice `I6`, no port and no engine: 64 distinct dual points give 1
  distinct hash, and `dict.fromkeys` costs 0.585 / 1.961 / 6.525 / 22.926 s at n = 8 / 16 / 32 / 64.
- Primal points are unaffected: 125 points give 80 distinct hashes and 0.001 s at n = 64.

Bounded faulthandler stack, and the provider-only reproducer: `.tmp/provider-report-dual-point-hash.md`
(untracked, `.tmp/*` is ignored). Reported to `sage-extensions` at the pin; the owner is repairing it.

Consequence for the port: **no routing change is justified.** Both `cell_cone` branches — PPL
`from_rays` and the Normaliz direct path — end in the same `RationalPolyhedralCone.__init__`, so
neither avoids the wall. Case 37 is blocked on the provider, not on port routing.

## All 80 frozen fixture cones agree exactly between Normaliz and PPL

Sweep of every fixture in both modes, engine-only, 11 s total: for all 40 fixtures x 2 modes the
minimal H-presentations coincide exactly — same facet count, same primitive signed rows, same set.
Every cone is pointed (0 equations throughout). Raw table: `~/sage-image/logs/engine-agreement.json`.

Facet counts are what separate the capped cases, and they do not separate them the way the caps
suggested:

- `Result_Enumeration_37` has **570** facets in both modes (202 / 112 rays) — by far the largest, and
  the only fixture near the provider's hash-collapse wall. This is why case 37 is not a cost problem
  but a blocker.
- Next largest: `Result_Enumeration_24` total 104 facets (137 rays), `Result_Enumeration_38`
  isotropic 56, `Result_Enumeration_26` 48, `Result_Enumeration_30/31` 32.
- `Result_Enumeration_33` has **11** facets (16 rays) and `Result_Enumeration_39` **8** (22 rays).
  Their cone construction is milliseconds, so their cost is not the engines and not this defect: it
  is the traversal over cells. Case 33 total is measured at 558 s and passes.

So of the six 180 s caps, one (37, both modes) is a provider blocker and the other four are traversal
cost that a larger cap may simply absorb — 33 total already does so at 558 s. That is the next
measurement worth making, and it needs no code change.

## A_PERFECT: the six 180 s caps reduced to three verdicts and one real blocker

Each capped case rerun unchanged on the same port `52eeedfa` + `sage-extensions 7b222e86` + base image
`b0cd7d86` with the owner Bliss patches, one per invocation, under a 1800 s cap.
Raw table: `~/sage-image/perfect-33-39-longcap.log`.

| case | mode | 180 s cap | 1800 s cap |
| --- | --- | --- | --- |
| `Result_Enumeration_33` | total | cap | **passed, 558 s** |
| `Result_Enumeration_33` | isotropic | cap | exceeded the 1800 s cap |
| `Result_Enumeration_37` | total | cap | blocked on the provider defect below |
| `Result_Enumeration_37` | isotropic | cap | blocked on the provider defect below |
| `Result_Enumeration_39` | total | cap | **passed, 237 s** |
| `Result_Enumeration_39` | isotropic | cap | **passed, 261 s** |

So of the 92 A_PERFECT cases: 86 pass inside 180 s, three more pass inside 1800 s (33 total, 39 total,
39 isotropic), and three remain without a verdict — 33 isotropic over the 1800 s cap, and both 37
modes blocked on the provider `finite_ordered_set` defect. Exceeding a cap is recorded as exactly
that, never as nontermination.

The 39 passes are consistent with the engine sweep: `Result_Enumeration_39` has 8 facets and 22 rays
in total mode, 8 facets and 10 rays isotropic, so its cone construction is milliseconds and the 237 s
/ 261 s is traversal. `Result_Enumeration_33` isotropic (11 facets, 14 rays) is likewise traversal
cost, not cone construction — it is the one remaining cost question that a larger cap may still
absorb, and it needs no code change.

These bounded runs do not substitute for the provider repair: case 37's cost is the dual-point hash
collapse, and it stays blocked until the repaired pin is published and consumed.

## A_GROUPS: the Mertens generating-set case passes on the repaired runtime

`tests/acceptance/test_orthogonal_group.sage::test_mertens_omega_and_minus_one_generate_o_l`
— run unchanged, single node, one invocation — **passes, 7.11 s** on port `52eeedfa` +
`sage-extensions 7b222e86` + base image `b0cd7d86` with the owner Bliss patches
(installed `bliss.cpython-314-x86_64-linux-gnu.so` sha256 `d941337e…`, `/opt/bliss-repair.json`).
It reports PASSED, not SKIPPED, so `require(MARKED_AUTOMORPHISMS)` was satisfied and the body ran.

This case previously failed with `MemoryError` in
`sage/groups/perm_gps/partn_ref/automorphism_group_canonical_label.pyx:513`, so it is frozen-consumer
evidence that the owner Bliss repair addresses the `partn_ref` canonicalization defect — and that the
#15 branch sample did not cover this branch.

The case is the rank-3 form of determinant `-155`: it assembles the full `O(L)` generator family and
checks the finite images of the recorded Mertens Omega generators together with `-I` against the
computed group, via `consumer.py:120–134` `assert_same_finite_images` and `discriminant_image_order`.
A_GROUPS now has one accepted case; the #33 contract for that family is a branch-representative
sample recorded by case id, not the whole file.

Do not read the whole file as cheap: one bounded invocation of
`test_orthogonal_group.sage` with only the known-open `a7_plus_diagonal` node deselected passed six
nodes and then exceeded a 900 s cap. That is a cap on the file as a whole, not a failure of any node,
and the family must stay per-case. The A7 discriminant-image node remains open and is not covered by
this result.

## Provider pin 69907aef unblocks case 37; targeted acceptance rerun

`sage-extensions` moved from `7b222e86` to **`69907aef538a8a4fe08c7431d2228afdd9d62927`**
(dual-point hashing, including primitive-line reflection repair `24885f3`), pinned in
`pyproject.toml`; the justfile provisioning reader extracts the new revision from that file, so the
two cannot drift. Nothing else in the tree carried the old sha (`uv.lock` does not record it).

Out-of-band check of the repair, standard lattice `I6`, no port and no engine: 64 distinct dual points
now give **64 distinct hashes** (previously 1), `dict.fromkeys` 0.142 s at n=64 (previously 22.9 s),
`finite_ordered_set` 0.142 s (previously 49.8 s). Linear, as the contract requires.

**`Result_Enumeration_37` now passes in both modes**, unchanged cases, one invocation each, with the
repaired provider bound at `/repos/sage-extensions` and `inspect.getfile(sage_extensions)` confirming
that path:

| case | mode | before | after |
| --- | --- | --- | --- |
| `Result_Enumeration_37` | total | blocked, `cell_cone` >1800 s | **passed, 36.80 s** |
| `Result_Enumeration_37` | isotropic | blocked, `cell_cone` >1800 s | **passed, 46.61 s** |

That is a 202-ray / 570-facet and a 112-ray / 570-facet cone: the engine work was always 0.05–0.11 s,
so the entire previous cost was the provider container, exactly as the report said. No port edit was
needed or made; routing was never the problem.

A_PERFECT is now **91 of 92 cases with a verdict**: 86 inside 180 s, and 33 total 558 s, 39 total
237 s, 39 isotropic 261 s, 37 total 38 s, 37 isotropic 48 s. The single remaining case is
`Result_Enumeration_33` isotropic, which exceeded a 1800 s cap; its cone has 11 facets, so the provider
repair is not expected to change it and its cost is traversal.

Targeted acceptance for the two provider changes, all passing with the repaired pin:

- `tests/test_perfect_cell_transporter.py` — 46 passed, 25.34 s (owned-cone construction throughout).
- `test_perfect_domains.sage::test_large_cell_cone_uses_exact_normaliz_h_presentation` — passed, 9.00 s.
- `test_reflective.sage` binary-cycle / permuted-basis nodes (reflection repair `24885f3`) —
  `test_binary_cycle_coordinate_order_and_period_certificate` and
  `test_binary_cycle_period_respects_permuted_coordinate_basis` — 2 passed, 2.38 s.

The raw-engine 80/80 agreement stands separately and is unaffected: it is an engine-only statement
about Normaliz and PPL on identical ray sets, and it never depended on the provider container.

## Enriques case_1: instrumented 1800 s run, phase evidence, no hang

Job `jmvxdwf53`, node
`test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]`, cap 1800 s,
provider `e83c4d6` bound at `/repos/sage-extensions`. It reached its cap: `elapsed=1801s`. Raw log with
the stack dumps: `.nucleus/jobs/jmvxdwf53.log`. No job was left running afterwards.

Node identity confirmed against `tests/fixtures/enriques_87_polarizations.json` before running: case 1,
`face_polarization` `[-1,-1,0,0,0,0,0,0,0,0]`, `face_polarization_norm` 2, `group_order` 696729600,
`line_orbits` 5, `plane_orbits` 9, `flag_orbits` 18.

**Instrumentation matters, and the first attempt produced nothing.** `faulthandler` writes to file
descriptor 2, and pytest captures at fd level by default, so a run without `-s` discards every dump and
the invocation ends as an opaque cap — exactly the failure mode this run was meant to replace. The
working form is `python -u -m pytest ... -s` plus a small `PYTHONPATH` plugin whose
`pytest_runtest_setup` arms `faulthandler.dump_traceback_later(60, repeat=True)` and whose teardown
disarms it. 29 dumps were captured.

**The case is progressing, not stuck, and the phase split is not where it first looked.**
`faulthandler` truncates when the stack is deep: only **14 of the 29 dumps carry a complete stack**
(40–74 frames); the other 15 hold 0–7 frames and cannot place the call at all. Attributing only from
the 14 complete dumps:

| window | test frame | meaning |
| --- | --- | --- |
| 60–120 s | `test_boundaries.sage:213`, inside `_gamma_h` | the port's coinvariant / pointed route, including the only Julia-engine minute (`interface.py:358 _request_unlocked` under `lattice_engines.py:1146 definite_isometry_group`) |
| 300–1680 s | `test_boundaries.sage:226` | `_gamma_h` has already returned; the test is computing `target.subgroup_on(tuple(g.discriminant_morphism() for g in group.group_generators()))` |

So `_gamma_h` costs roughly the first two to five minutes, and the remaining ~23 minutes are the
discriminant-form image computation on `U + E8(-1)`. The recorded `group_order` 696729600 is
`|W(E8)|`, and the frames in that window are `fqf_orthogonal.py:112`/`:484`, `abelian_aut.py:203
matrix`, `abelian_group_gap.py exponents`, `torsion_quadratic_module.py:149 _mul_` and
`free_module.py _element_constructor_` — building the subgroup of the finite quadratic form's
orthogonal group and asking for its cardinality, not the port's recursion.

**An earlier reading of this log, recorded a few minutes before this correction, put the whole run
inside `_gamma_h`. That was wrong.** It came from searching all 29 dumps instead of the complete ones;
truncation hides the test frame in a third of them. `_gamma_h` is intermediate, as stated before the
run, and the expensive step is at `:226`.

The re-entered frames, which are where the time actually goes:

- `sage/structure/free_module.py` `_element_constructor_` — at 180, 420, 1260, 1560 s.
- `sage/groups/abelian_gps/abelian_group_gap.py` `exponents` / `:154` — at 240, 540, 600, 660, 1080 s.
- `sage/modules/fqf_morphism`-side and `fqf_orthogonal.py:112 __call__`, `:484 _act_`,
  `abelian_aut.py:203 matrix` — 300, 1620 s.
- Owned category overhead inside the same loop: `set_categories.py` `position_of` / `__call__`,
  `refine.py:362 construction_scope`, `mor_foundation.py:403 operation`, `category.py:706 __contains__`,
  `category_with_axiom.py __classget__`, `category.py __classcall__`, `mor_categories.py Of`,
  `indexed_families.py:424,90`, `torsion_quadratic_module.py:149 _mul_`.

So the cost is a long chain of small constructions — element construction in the formed/torsion-form
modules and GAP abelian-group exponent computations during the discriminant-form action — not a single
long call and not a hang. `_gamma_h` is only the beginning: the test then computes
`image.cardinality()` against the recorded 696729600 and `_counts` against 5/9/18.

This is a cost finding, not an acceptance result: case_1 has **no verdict**. It exceeded a cap, which is
recorded as exactly that.

## Runtime hold and revision coverage

No shared-runtime job was active at the end of this run, and no further PRoot/bwrap/udocker command was
started from this project. The next shared-runtime run waits on the coordinated maintenance boundary.

Revision coverage is mixed and must not be flattened into "91/92 on the current pin":

- At `e83c4d6` (current pin): `Result_Enumeration_37` total **40.12 s**, isotropic **34.51 s**, and the
  Mertens `O(L)` node **6.54 s**, all passing.
- At `69907ae`: 37 both modes (36.80 s / 46.61 s), Mertens 7.11 s, and the targeted suite —
  `test_perfect_cell_transporter.py` 46 passed, the large-cell Normaliz node, and the two reflective
  binary-cycle nodes.
- At `7b222e86` (previous pin): 33 total 558 s, 39 total 237 s, 39 isotropic 261 s, and the 86 cases
  inside the 180 s cap.

`Result_Enumeration_33` isotropic remains the only historically capped A_PERFECT case, and it now needs
a bounded phase diagnosis on the current provider — the same instrumented treatment case_1 just had,
not another opaque cap.

## Runtime resumed: binds verified, targeted acceptance re-taken at the current pin

After the shared-runtime maintenance the baseline `/repos` pins were restored to the published lock,
which is *older* than this project's declared provider. The per-process bind was verified three ways
rather than assumed, because the failure mode here is silent: a run that picks up the restored baseline
provider looks exactly like a run on the declared pin until something slow or wrong happens.

| check | value |
| --- | --- |
| host tree `~/sage-extensions-e83c4d6` `module_morphisms.py` sha256 | `892172a6b10d3acb` |
| same file as seen inside the run at `/repos/sage-extensions` | `892172a6b10d3acb` (match) |
| restored baseline `/repos` copy in the rootfs | `be6b82f51bcd1726` (different — not what runs) |
| `sage_extensions.__file__` | `/repos/sage-extensions/src/sage_extensions/__init__.py` |
| repair behaviour on `I6` | 64 distinct dual points give **64 distinct hashes** |

Targeted acceptance re-taken at `e83c4d6`, so this is no longer mixed with the `69907ae` run:

- `tests/test_perfect_cell_transporter.py` — 46 passed, 25.41 s.
- `test_perfect_domains.sage::test_large_cell_cone_uses_exact_normaliz_h_presentation` — passed, 10.29 s.
- `test_reflective.sage` binary-cycle / permuted-basis nodes — 2 passed, 2.63 s.

Cumulative coverage **at `e83c4d6`** is now: `Result_Enumeration_37` total 40.12 s and isotropic
34.51 s, the Mertens `O(L)` node 6.54 s, and the four targeted items above. Everything else — the 86
cases inside the 180 s cap, 33 total 558 s, 39 total 237 s, 39 isotropic 261 s — still dates from
`7b222e86` and would need re-taking for a single-revision A_PERFECT statement.

The Maxima environment recipe correction found during readback belongs to the image owner and is not
duplicated here.

## Every formerly capped A_PERFECT case now passes at the current pin

All six capped cases, unchanged, one invocation each at `e83c4d6` with the verified bind:

| case | mode | at `7b222e86` | at `e83c4d6` |
| --- | --- | --- | --- |
| `Result_Enumeration_33` | total | passed, 558 s | **passed, 36.17 s** |
| `Result_Enumeration_33` | isotropic | exceeded 1800 s cap | **passed, 56.30 s** |
| `Result_Enumeration_37` | total | blocked, `cell_cone` >1800 s | **passed, 40.12 s** |
| `Result_Enumeration_37` | isotropic | blocked, `cell_cone` >1800 s | **passed, 34.51 s** |
| `Result_Enumeration_39` | total | passed, 237 s | **passed, 34.41 s** |
| `Result_Enumeration_39` | isotropic | passed, 261 s | **passed, 23.52 s** |

The dual-point hash defect was therefore not a case-37 curiosity. `RationalPolyhedralCone.__init__`
does O(facets²) comparisons at ~17 ms each, and a traversal builds an owned cone per cell, so an
11-facet cell costs about 2 s per cone and a traversal touching hundreds of cells reaches the cap.
That is why the whole upper tail of A_PERFECT — not just the 570-facet specimen — was sitting at
558 s, 261 s, 237 s and beyond 1800 s. Fixing the hash took the entire tail to 24–56 s.

`Result_Enumeration_33` isotropic needed no phase diagnosis after all: it finished in 56 s, inside the
first 60 s dump interval, so `stackdump` recorded zero dumps. The instrumented run was set up and never
had to say anything, which is the right outcome.

A per-case sweep of all 40 fixtures in both modes at `e83c4d6` is running (job `jmvxfxg7n`, 300 s cap
per case, one invocation each) so that A_PERFECT can be stated at a single revision instead of the
mixed `7b222e86` / `69907ae` / `e83c4d6` coverage recorded above.

## A_PERFECT at a single revision: the whole file passes locally on `e83c4d6`

`tests/acceptance/test_perfect_domains.sage` at one revision, one provider, verified bind:

- **80 of 80** parametrized cases (`test_perfect_domain_orbit_counts_match_polyhedral_common` and
  `test_isotropic_mode_counts_match_polyhedral_common` for all 40 fixtures), each run as its own
  invocation with a 300 s cap. Raw per-case table: `.nucleus/jobs/jmvxg3dnz.log`. **No failures, no
  errors, no caps, no `NO_VERDICT` lines.**
  - sum of case times 642 s (against 2387 s for the 86 cases that ran inside 180 s on `7b222e86`, plus
    four cases that needed 237–558 s and two that never returned)
  - total mode: 40 cases, 306 s, slowest 46 s (`Result_Enumeration_37`)
  - isotropic mode: 40 cases, 336 s, slowest 65 s (`Result_Enumeration_33`)
- **12 of 12** standalone semantic tests in the same file (`-k "not counts_match_polyhedral_common"`),
  12.96 s, slowest `test_large_cell_cone_uses_exact_normaliz_h_presentation` at 9.32 s.

So 92 of 92 tests in the file pass at `e83c4d6`, and the slowest case in the entire family is 65 s —
a 180 s per-case cap is comfortable for the first time, where before it was exceeded by four cases and
never reached by two.

**This is local evidence, not hosted certification.** TODO.md is explicit that A_PERFECT requires "a
successful hosted run whose head SHA equals the implementation commit being accepted", and the
baseline run 40807218493 does not certify subsequently changed code. No hosted Acceptance run at
`a848fd7c` has been observed from here. The A_PERFECT row stays open until that run exists.

Two further consequences worth carrying:

- The 180 s cap that produced the original "86 of 92" record was measuring the provider defect, not the
  algorithms. The record above should be read with that in mind rather than taken as evidence that the
  port's traversal was ever slow in its own right.
- `Result_Enumeration_33` isotropic was the last case without a verdict and it never needed the phase
  diagnosis I queued for it; it now runs in 65 s.

## Port fix: integral coordinates over a rational lattice (the last known isometry failures)

`tests/acceptance/test_isometry.sage::test_change_of_basis_gives_an_isometric_lattice` — all six
cases failed at HEAD with

```
TypeError: int() argument must be a string, a bytes-like object or a real number,
           not 'OwnedJoinCategory.element_class'
```

raised from `recursive.py:154`, `int(coordinates(label))`.

**This was a latent port defect, not a provider regression.** The same probe run against three
different providers gives identical behaviour — `e83c4d6`, `69907ae`, and the restored baseline
`/repos` pin — so no pin change caused it. The isolating observation:

| lattice over | `to_vector()(label)` type | `int(...)` |
| --- | --- | --- |
| `ZZ` | `_OwnedRingParent_with_category.element_class` | works |
| `QQ` | `OwnedJoinCategory.element_class` | **`TypeError`** |

A lattice over `QQ` states a coordinate as a rational that does not admit `int`. The failing path is
`eichler.py:464` → `_isometry_from_row_matrix(rational_envelope, generator)`: the integer row-matrix
boundary being handed a rational envelope. It became *reachable* only once the Bliss repair stopped
killing these cases earlier with `MemoryError`, which is why the earlier record saw "1 failed, 165
passed" and this round saw six failures.

The repair reuses the helper the port already has rather than adding a near-duplicate:
`edgewalk_scalars._integer`, whose contract is "a rational entry with a denominator is a defect, not a
value to truncate", converting through `numerator()`/`denominator()` and raising `ArithmeticError` on
a genuine fraction. Applied at the four coordinate extractions that cross into a possibly-rational
lattice — `recursive.py:127` (`_column_matrix`), `:155` (`_isometry_from_row_matrix`), `:185`
(`_isometry_from_column_matrix`), and `eichler.py:922,926` (`_pullback_envelope_vector`).

All six cases pass, one invocation each, 900 s cap:

| case | time |
| --- | --- |
| `ci_indef_U_2U` | 3 s |
| `ci_indef_U_2U_A2` | 113 s |
| `ci_indef_U_2U_A3` | 151 s |
| `ci_indef_U_2U_A2_A2` | 122 s |
| `ci_indef_U_U_E7` | 336 s |
| `ci_indef_U_2U_2E8` | 575 s |

`ci_indef_U_2U_2E8` is the case that historically died with `MemoryError`; `ci_indef_U_U_E7` is the
case that failed in hosted run 41242413071. The other 160 cases in the file were not re-run in this
round, so this closes the known failures and does not by itself certify the full file.

Independent confirmation: the library separately reproduced `Γ_h.order()` failing at
`recursive.py:154` on an owned rational element, and published `89a0523d` with regressions and no
library code change.

## Correcting the line-226 attribution for Enriques case_1

The earlier entry attributed the 300–1680 s window to `target.subgroup_on(...)`. That is too narrow.
The traced frames there were `fqf_orthogonal`, `abelian_aut.matrix`, `abelian_group_gap.exponents`,
`torsion_quadratic_module._mul_` and `free_module._element_constructor_`, which is consistent with the
generator request inside that expression — `group.group_generators()` — as well as with the subgroup
construction. The library's readback pins it further: supplied finite-form generators and the subgroup
construction are fast there and do not enumerate the finite group, while Γ_h's retained
character/intersection generating morphism requests the **preimage group's Eichler/orbit-cover family
through this port's `OrbitCoverModel`**. So a meaningful share of that window is port-owned family
construction, and attributing all of it to `subgroup_on` was wrong.

The next measurement is therefore to instrument the retained `OrbitCoverModel` family construction
itself — the real path, not a replacement library subgroup algorithm — with the original Enriques
assertions unchanged.

## Enriques case_1: OrbitCoverModel is not on this path (measured, not assumed)

Job `jmvxi4o1u`, the unchanged node, `OrbitCoverModel`'s 14 callable members wrapped in place
(behaviour-preserving: plain functions and `@staticmethod`s; caching lives in module-level dicts, so no
cache is disabled). Sampling ran to 29 faulthandler dumps (~1740 s) before the cap.

**Result: zero calls to any of the 14 members in that window.** The teardown table printed its header
with no rows.

That is a real negative, not a broken wrapper. The control run confirms the wrapper counts calls,
including calls that raise:

```
instrumented members: 14
call raised as expected: AttributeError
TOTALS: {'source_generating_family': 0.0}
CALLS: {'source_generating_family': 1}
CONTROL OK: the wrapper observes calls
```

So for case 1 the library's readback does not reproduce: Γ_h's retained generating morphism does not
reach this port's `OrbitCoverModel` inside the sampled window. Combined with the earlier stack
evidence — test frame at `test_boundaries.sage:226`, top frames `fqf_orthogonal.py`,
`abelian_aut.matrix`, `abelian_group_gap.exponents`, `torsion_quadratic_module._mul_`,
`free_module._element_constructor_` — the time is in the discriminant-form layer on generators already
in hand (`g.discriminant_morphism()` and the subgroup/cardinality work), not in family construction.
Case 1 remains **unresolved, capped**, and the OrbitCoverModel hypothesis is not confirmed for it.

### Instrumentation defect in my own probe, and the fix

I made the table print from a `SIGTERM` handler that raised `SystemExit`. pytest catches
`BaseException` in `wrap_session`, so the run did not exit: it continued past its 1800 s cap with
faulthandler already cancelled and produced no further evidence until I stopped job `jmvxi4o1u`
myself. Two rules for the next instrumented run:

- Let `timeout` do the killing, with `timeout -k` so a wedged child is still killed, and never raise
  from a signal handler inside pytest.
- Print the partial table from the periodic dump instead of at exit, so a capped run still reports its
  numbers and its stacks together.

### Next exact analysis

Instrument the layers that the stacks actually name, on the unchanged node: `discriminant_morphism` on
each group generator, the FQF `subgroup_on`, `cardinality`, and `group_generators` itself, with the
table emitted on every periodic dump. Original Enriques assertions stay untouched.

## Discriminant-layer instrumentation: design and positive control

The OrbitCoverModel attribution is **refuted on this path** for Enriques case 1: 14 wrapped members,
29 dumps (~1740 s), zero calls, with a control proving the wrapper counts calls. Recorded in
`7af4b0b5`; it stands.

The replacement probe splits the four phases that the stacks actually name, each with its own
wrapper and its own marker:

| phase | wrapped | definition |
| --- | --- | --- |
| generator enumeration | `generator_enumeration` | `OwnedGroups.ParentMethods.group_generators` (`groups.py:5430`) |
| each morphism conversion | `discriminant_morphism` | `LatticeIsometryMethods.discriminant_morphism` (`lattice_morphisms.py:1795`) |
| subgroup admission | `subgroup_admission` | `TorsionFormOrthogonalGroup.subgroup_on` (`torsion_form_modules.py:1504`) |
| cardinality of its result | `result_cardinality` | `cardinality` on the class defining it on the returned group |

Enclosing classes were resolved by AST scan of the pinned provider, not guessed. Markers are
unbuffered; only the first 40 calls of each phase are marked, plus any single call over 35 s.
A daemon thread prints the running table every 60 s so an externally capped run still reports its
numbers, and `faulthandler` dumps stacks on the same cadence. There is no signal handler and no
`SystemExit`: the external `timeout -k 10 1800` is the only stopper, which is what the previous run
lacked.

**The control earned its place.** The first version recorded `subgroup_admission` *after* the call
returned, so a raising call went uncounted — the control caught it as `recorded 0 -> 0`. Fixed with
`try/finally`; the corrected control records `0 -> 1` for all three wrappers, including calls that
raise. Evidence kept at `~/sage-image/logs/discprof-control.log`, and the failed first control is
recorded there too rather than dropped.

The run itself uses the unchanged node and the unchanged assertions. Per instruction, whatever it
identifies gets handed to the owning library or Sage project as a minimal reproduction, and no
hypothesized owner is patched.

## The slow call in Enriques case_1: exact parent, presentation, caller

From the phase table at t=256.7 s of the instrumented run, and then from source — not from the test's
outer line and not from the aggregate:

| phase | calls | total | worst |
| --- | --- | --- | --- |
| `generator_enumeration` | 21 | 177.52 s | **118.67 s** |
| `discriminant_morphism` | 21 | 13.75 s | 11.65 s |
| `subgroup_admission` | 1 | 0.03 s | 0.03 s |
| `result_cardinality` | 957 | 0.01 s | 0.00 s |

**Exact parent.** The stabilizer of the case-1 polarization inside `O(M)`, with `M = U(2) + E8(-2)`:

```
lattice_centralizers.py:852   stabilizer = extension.invariant.O().stabilizer(invariant_polarization)
lattice_centralizers.py:985   invariant_image = invariant_subgroup.image_under(...)   <- requests group_generators
```

**Retained presentation.** A marked-lattice stabilizer subgroup. Its resolution realises through:

```
group_generators -> selected_group_resolution -> realise
  -> lattice_morphisms.vector_stabilizer_generators           :3664
  -> _marking_stabilizer_generators                          :3691
  -> lattices._marked_automorphism_generators                :1874
  -> port marked_lattices.automorphism_group_generators        :109
  -> port recursive.pointed_automorphism_generators            :591
  -> recursive._lattice_action :457 -> orthogonal_group_generators :333
  -> lattices._core_automorphism_generators                  :1849
  -> lattice_morphisms._computed_group_generators :3575 / _engine_group :3456
  -> lattice_engines.definite_isometry_group                 :1146 -> Julia
```

`eichler.OrbitCoverModel` does not appear anywhere on that chain, which is exactly why the separate
OrbitCoverModel probe measured zero calls. The two measurements agree instead of contradicting.

**Initial versus repeated.** `group_generators` is not cached, but `selected_resolution`
(`objects.py:383`) replaces the `_PendingResolution` with the realised resolution in the per-owner
registry, so a second request on the same parent returns the stored object. The 177.52 s therefore
cannot be 21 equal charges: it is a small number of first realisations, one of them 118.67 s, plus
cheap repeats. This is reasoned from source and is deliberately not asserted as measured — the running
probe cannot attribute per-call identity, because my marker index is "prior completions + 1", so
nested calls share a number (the earlier log shows `#20 start` twice).

**Actionable minimal path, published but not run:** [.tmp/probe-gamma-h-stabilizer.py](.tmp/probe-gamma-h-stabilizer.py)
rebuilds exactly that parent the way `_gamma_h` does and times the initial realisation against a repeat
request in one process, so the initial/repeat split becomes a measurement rather than an inference. It
runs after the sole bounded probe finishes; nothing else is launched now.

No owner is patched and the frozen assertions and generator semantics are untouched.

## Correction: the cheap `subgroup_admission` reading was the wrong call

The t=256.7 s table showed `subgroup_admission: 1 call, 0.03 s`, and I read that as "admission is free".
**That was wrong, and the aggregate table is exactly why.** The test frame is still
`test_boundaries.sage:226` at t=540 s, so the `subgroup_on` that line 226 is waiting on has not
returned. The one recorded admission is a *different, nested* admission — `lattice_centralizers.py:987`
calls `subgroup_on` as well — which completed cheaply inside `_gamma_h`'s coinvariant route.

What the dumps and the table actually say together:

| dump | t | test frame | top main-thread frame |
| --- | --- | --- | --- |
| 1–2 | 60–120 s | `:213` `_gamma_h` | `interface.py _request_unlocked`, `libgap_mixin.py is_finite` |
| 4–5 | 240–300 s | `:226` | `abelian_group_gap.py exponents` |
| 7 | 420 s | `:226` | `abelian_aut.py __init__` |
| 8 | 480 s | `:226` | `groups.py representative` |
| 9 | 540 s | `:226` | `owned_category.py _object_of` |

And by t=496.8 s the wrapped totals had **frozen**: `generator_enumeration` 21 calls / 177.52 s and
`discriminant_morphism` 21 calls / 13.75 s unchanged since t=256.7 s, while `result_cardinality` grew
from 957 to **137 242** calls for a total of 0.91 s. So of ~497 s elapsed, only ~192 s is inside the
four wrapped phases; roughly 300 s is unaccounted, and it is inside the still-running line-226
`subgroup_on`, dominated by per-object construction (`_object_of`, `representative`,
`abelian_aut.__init__`) with a huge number of individually trivial cardinality calls.

Two corrections follow, both recorded rather than quietly fixed:

- The slow generator call is real and identified (the polarization stabilizer in `O(M)`,
  `M = U(2) + E8(-2)`), but it is **not the whole case**: it is 177.5 s of the first ~380 s.
- Enriques case_1's remaining cost is inside the outer `subgroup_on` at line 226, which no wrapper of
  mine has yet timed because it has not returned.

`discprof.py` now tracks **in-flight** calls and prints `RUNNING depth=… oldest=…s` for any wrapped
phase still executing, so the next run cannot repeat this misreading — a capping run reports what is
running, not only what has finished. Not yet exercised; it waits for the current probe to end.

## Resolved: the call that never returns is a `description` f-string

The bounded probe reached its cap: `pytest_exit=124 elapsed=1800s`, 29 dumps, process verified gone.
Wrapped totals at t=1757.2 s: `group_generators` 21 / 177.52 s and `discriminant_morphism` 21 /
13.75 s (both frozen since t≈380 s), `cardinality` 1 570 390 calls / 9.95 s, nested `subgroup_on`
1 / 0.03 s. The line-226 `subgroup_on` never returned.

Last complete dump (76 frames, t≈1740 s), innermost first:

```
abelian_group_gap.py:152 exponents
abelian_aut.py:203 matrix
fqf_orthogonal.py:484 _act_ / :288 _preserves_form / :267 _element_constructor_
groups.py:305 _element_from_engine / :525, :472 representative
groups.py:486 _TransversalCoset._repr_
finite_ordered_sets.py:587 _repr_
groups.py:2133 _symmetric_group          <- caller
```

`_symmetric_group` builds `description=(f"Symmetric group S_{size}" if standard else
f"Symmetric group on {points}")`. The f-string is eager, so on the non-standard branch it formats a
finite ordered set of transversal cosets; that `_repr_` iterates every point and each coset's `_repr_`
calls `representative()`, materialising and form-validating a real element through
`_element_from_engine` into Sage's `FqfOrthogonalGroup`. A description string that is usually never
read costs a full materialisation of the coset space.

Handed to sage-extensions with the object path, the exact lines, the measured counters and the
expected contract: `.tmp/handoff-owner-symmetric-group-description.md`. No library code was patched
here.

## A_STABILIZERS coverage: one real pass, and the rest is not "small"

Run on the current declared pin `e83c4d6`, one invocation per node, 600 s cap. Kept strictly separate
from the Enriques blocker.

| node | result |
| --- | --- |
| `test_vector_stabilizer_has_oscars_order` | **passed, 65 s** |
| `test_equivalence_witness_carries_a_representative_to_its_reflection[r_-1]` | capped at 600 s, **no verdict** |

I stopped the remaining fifteen nodes rather than spend ~2.5 h producing more caps, and then measured
*why* instead of guessing. The node's unit is one pointed E10 Isom, repeated over the 61 recorded orbit
representatives. Timing them individually:

| representative norm `q` | one pointed Isom |
| --- | --- |
| 2 | **0.22 s** |
| 6, 8, 10, 12, 18 | 56–74 s each |
| 14, 16 | 36–60 s each |
| 20 | 57–58 s each |

So the cost is strongly norm-dependent and plateaus near ~56 s per representative once `q >= 6`; a
900 s probe cap was reached after ~13 of the 61. Extrapolating, one `r_*` node is roughly
**50–60 minutes**, and the ten of them roughly 8–10 hours — plus the six E10 square nodes, whose
`vector_classes` enumeration was not measured.

Conclusion for #33: the remaining A_STABILIZERS sample is **not a small case set**. The exact ids are
pinned (`r_-1` … `r_8`, squares `[0|2|4|6|8|10]`, and the OSCAR `stabilizer_of_100` node that now
passes), but each root node needs its own budget rather than a bounded sequence, and a per-node ~1 h
budget is a scheduling decision, not a probe detail.

Also worth carrying: the norm dependence is a real property of the pointed route at rank 10, not a
fixture artefact — the same 61 representatives, same reflection, differing only in `q`.

## One instrumented pointed-Isom unit: engine-bound, not representation-bound

Single unit, exact frozen specimen — `r_-1`, orbit representative index 3,
`[-1,-3,0,0,0,0,0,0,0,0]`, `q = 6`, asserted in the probe — uninstrumented 56.12 s, profiled 143.24 s
(so percentages below are attribution, not a runtime claim), with 30 s in-flight dumps.

| self-time | calls | site | judgement |
| --- | --- | --- | --- |
| 58.39 s | 4 | `{readline}` on the Julia pipe | **necessary mathematics**: `definite_isometry_group` (`:1146`) computing the definite isometry group |
| 27.44 s | 1 | `libgap_mixin.py:181 is_finite` | **candidate**, caller `_engine_finiteness` (`groups.py:158`) |
| ~2 s each | 3–14 M | `_from_engine_element`, `__contains__`, `__class__call__`, `isinstance`, `_element_constructor_` | broad owned-object-model churn, spread thin |

**No representation defect in this unit**: total time inside `repr` is 1.38 s over 139 205 calls and
every `_repr_` is ≤ 0.03 s. That is the opposite of the Enriques case-1 finding, and the two must not
be conflated — this plateau is engine mathematics plus object-model overhead, not formatting.

Handed to its owner as a *candidate*, not a conclusion:
`.tmp/handoff-owner-engine-finiteness.md`. The 27.44 s GAP `IsFinite` is measured and attributed to
`_engine_finiteness`, but grep of the pinned tree finds only the definition and the two
`assert _engine_finiteness(engine) is True` sites in `_gap_model_of_engine` (`:191`, `:198`), neither
passing `compute=True`, so which branch fires is the owner's to read. If it is avoidable, it is 27 s
of a 56 s unit.

**No runtime guarantee is extrapolated.** The 61-representative node and the ten `r_*` nodes stay
unresolved; the OSCAR node's 65 s pass remains the only new accepted A_STABILIZERS coverage. The
Enriques case_1 run is not repeated while its provider hand-off is open.

## Full unchanged r_-1 node running (one heavy node, nothing else)

Job **`jmvxlw71n`**. Exact pair, verified before launch: port `6e1be797c8a491d62390a50ba2e4b51d750c82df`
(clean worktree) + provider `e83c4d69db6fba9cadc560d6fd8a518bfe6c5dae` bound at `/repos/sage-extensions`.

- Node: `test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_-1]`,
  **unchanged**, original assertions, all 61 representatives.
- Cap: `timeout -k 30 7200` — 120 minutes with a 30 s kill-after grace; `python -u` and `-s` so output is
  unbuffered.
- Progress: `stackdump` arms `faulthandler` every 300 s, so a 2 h run yields ~24 in-flight stacks instead
  of one opaque wait. The plugin only arms a dumper; it does not touch the test.

This is the acceptance measurement the unit profile justified: the unit costs ~56 s per representative
once `q >= 6`, and the node's assertion covers all 61, so it needs a real budget rather than more unit
samples. No other heavy node is running, and Enriques case_1 is not re-run while its display blocker is
open.

The GAP `IsFinite` candidate stays exactly that — a separate owner question. It is **not** grounds to
skip or weaken a finiteness check, and nothing in this run bypasses one.

## A_REFLECTIVE sample cannot currently be recorded "by case id"

While the `r_-1` node runs, source inspection of the reflective family turned up a contract-level
problem rather than another cost problem.

`test_reflective.sage` parametrizes over `load_allcock_rank3_reflective()` with
`ids=lambda case: case["id"]`. That fixture is not id-unique:

| fixture | entries | distinct ids |
| --- | --- | --- |
| `allcock_rank3_reflective.json` | 8 595 | **176** |
| `reflective_forms_8821.json` | 8 821 | distinct |
| `root_systems_103.json` | 103 | distinct |

In detail: `allcock_rk3_1` occurs 374 times, `allcock_rk3_2` 373, `allcock_rk3_3` 367, `allcock_rk3_4`
363, and only 56 ids occur once. The duplicates are **not** byte-identical — 120 ids have differing
duplicates, i.e. genuinely different lattices sharing a label.

So pytest disambiguates those node ids by appending `0`, `1`, … and the #33 requirement that an
accepted A_REFLECTIVE sample be "recorded by case id" cannot be satisfied for this parametrization by
the id alone: `allcock_rk3_2` names 373 distinct lattices. A sample must be recorded by full node id
(including pytest's index suffix) or by an intrinsic discriminator, and the fixture ids want making
unique.

That is a frozen-suite/contract change and belongs to #33, with its reason recorded before adoption —
not something to patch locally. Recorded here so it is not discovered later as a sampling ambiguity.

Note also that the required A_REFLECTIVE sample is otherwise easy to name at node level:
`test_edgewalk_branch_representative_sample`, `test_u_plus_46_nonreflective_certificate_terminates_after_three_vertices`,
`test_root_systems_match_polyhedral_common[ListSimpleRootSystem_0]`,
`test_e10_simple_roots_match_dutour_sikiric_hulek`, and the four OSCAR Vinberg nodes
(`vinberg_I_1_10`, `vinberg_diag_1_-1_-13`, `vinberg_scharlau_U_A2`, `vinberg_diag_1_-1_-3`).

## Correction: the finiteness caller is real, and my earlier grep was truncated

The owner is right that `e83c4d6` has a `compute=True` caller at `groups.py:7936`, in
`GeneratedSubgroups._finiteness_decision`. My earlier note claimed a grep found only the definition and
two assert sites. **That was my error**: I piped the grep through `head -6`, which hid the other nine
references. The tree has eleven.

Reading each `stats` entry's `callers` map from the existing profile — no re-run — gives the actual
chain for the pointed-Isom unit:

```
libgap_mixin.py:181 is_finite                         1 call   27.44 s
  <- _engine_finiteness (groups.py:158)               1 call   27.44 s   compute=True branch
    <- _finiteness_decision (groups.py:7924)          1 call   27.85 s   @cached_method
      <- is_finite (set_categories.py:2415)         959 calls   27.86 s
        <- invariant_module_construction
           (rational_integral_stabilizers.py:368)            27.85 s
```

So the parent whose construction asks is `rational_integral_stabilizers.invariant_module_construction`
(`:368`). Of 959 `is_finite` calls, 27.85 s is one undecided parent; the placement guards and the cache
are working, and one parent is genuinely unplaced.

**I am not declaring the GAP call avoidable.** The decision is deliberate, cached, and reached only after
`OwnedFiniteGroups()` / `OwnedInfiniteGroups()` placement fails. A finite presentation never proves
finiteness, so any argument for skipping it must show that site already knows the answer. That is the
owner's call, and the corrected hand-off says so: `.tmp/handoff-owner-engine-finiteness.md`.

## Next clean boundary: pin the `_symmetric_group` fix

Provider main `059e666` publishes the eager-description fix (degree-based string, `natural_points` and
the action unchanged, targeted regression passing). The running `r_-1` node stays on `e83c4d6` — no
provider change while it runs. At its next clean boundary: pin `059e666` (or a later owner pin if
supplied), then rerun the original Enriques case_1 with `discprof` phase markers and in-flight
reporting. The old description diagnosis is already established and should not be re-derived.

## Pin decision at the next clean boundary: `9dc89ac`, not `059e666`

Owner supplied a later pin, **`9dc89ac1d9e4ab074da8872386aaf5b9df290f63`** ("test: pin the labelled-ranking
and admission contracts of the symmetric engine"), superseding `059e666`. Verified in the tree rather
than taken on report:

- `_symmetric_group` (`groups.py:2093`) now carries the docstring "The description names the degree and
  never formats the point set".
- `engine = SymmetricGroup(degree)` is unconditional — the canonical Sage `SymmetricGroup` engine for
  labelled finite sets, replacing the `PermutationGroup(())` / transposition-and-cycle construction that
  `e83c4d6` used on the non-standard branch.
- The `else f"Symmetric group on {points}"` f-string that caused the Enriques cap is gone, while
  `natural_points` and the checked ingress remain per the owner's description.

Prepared as a separate tree at `~/sage-extensions-9dc89ac` so the running `r_-1` job keeps its exact
`e83c4d6` pair. **`pyproject.toml` is deliberately untouched while `jmvxlw71n` runs** — the pin flips
only at its result or cap, so the recorded pair stays honest.

At the boundary: pin `9dc89ac`, commit and push, then rerun the original Enriques case_1 with the
existing `discprof` phase instrumentation and in-flight reporting. The description investigation is
settled; the rerun measures the fix, it does not re-derive the diagnosis.

## Provenance of the 92/92 A_PERFECT run (`jmvxg3dnz`) — port SHA

The #33 record needs the port checkout that executed the 92/92 A_PERFECT result. **The SHA was not
captured by the run.** Checked directly rather than assumed:

- `jmvxg3dnz.log` contains no 40-hex string and no revision marker of any kind; the job command never
  echoed one. `grep -nE "[0-9a-f]{40}"` returns nothing, and no line matches `port|commit|rev`.
- Nothing in `.nucleus/jobs/` records a checkout state either.

What *is* established, and how far it goes:

| claim | status | basis |
| --- | --- | --- |
| HEAD at launch was **`a848fd7c`** | **derived, not captured** | job start 06:24:37 UTC falls between commit `a848fd7c` at 06:20:11 and the next commit `40bc5ef2` at 06:35:59; no commit exists in between |
| the executed `src/` and `tests/` equal `40bc5ef2`'s | **provable** | tree hashes are identical at `79de4190`, `a848fd7c`, `40bc5ef2`: `src` = `17144b985368f53b4af41c7d03f93dfe07f02e53`, `tests` = `7d02ecd6ec0493dc076a95491dc760e679a14b07` |
| nothing else changed in the window | **provable** | `git diff --stat 79de4190 40bc5ef2` is `1 file changed, 166 insertions(+)`, all of it `.tmp/HANDOFF.md` |

Residual uncertainty, stated rather than glossed: tree-hash equality is a statement about commits, not
about the live worktree. If an uncommitted edit to `src/` or `tests/` had existed at launch and was
later reverted, the hashes would not show it. The next commit's diff is docs-only, which is evidence
against such an edit but not proof.

**For the #33 record, the honest entry is:** port `a848fd7c` (derived — HEAD at launch by timestamp
bracketing), provider `e83c4d6` (recorded), with the executed source provably identical to
`40bc5ef2` because every commit in the window touched only `.tmp/HANDOFF.md`. Naming `40bc5ef2` alone
is not sufficient, as noted: it is the commit that *recorded* the result, made 11 minutes after the run
finished, and it differs from the launched tree only in documentation.

Future runs echo the port SHA and the provider SHA in the run header, so this derivation is not needed
again.
