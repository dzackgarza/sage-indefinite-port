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
