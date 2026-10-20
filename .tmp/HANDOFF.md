# Handoff: finish the entire sage-indefinite-port

## Immediate continuation

Resume TODO.md "Resume point: CELLS consumer verification, then independent families" using the decision table in AGENTS.md.
This is the continuation point; the historical diagnoses below do not select a different assignment.

The exact pending node is:

```text
tests/acceptance/test_perfect_domains.sage::test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_37]
```

The last broad run is `.tmp/accept-original-full-20261020.log`, with its adjacent source manifest, exit and timing files.
It passed 37 cases before a diagnostic cutoff.
Check for an existing live run before launching anything.
Use a current source snapshot and the exact node, without replaying the prefix or imposing the same unsupported 30-minute limit.
A completed result advances to the affected family; a wrong result or exception requires repairing the precise source mismatch; an interruption remains owed while independent families continue.

The current recursive.py already calls the assigned Witt-extension, pointwise-perpendicular-kernel and rational-lift operations.
Their complete integration and consumers are unverified, not unimplemented.
Historical provider findings below need current reproduction before any provider edit.
Preserve the full completion table; do not replace it with the current specimen.

When the next action changes, replace this continuation point with the exact remaining step, source correspondence and result/job pointer.
Do not append another competing handoff or restart completed investigation on a new turn.

## 1. Ambient task

Finish the entire port of the specific pinned upstream implementation into native Sage/Python, under the original issue tree rooted at #4 and its frozen acceptance contract.
This includes the specialized algorithms, their assigned generic callees, all applicable upstream fixtures, and removal of unfinished markers when justified.
The task is incomplete.
Passing individual cases or repairing a dependency does not replace this obligation.
No numerical completion percentage is established.

Use the leaf issue bodies and `references/source-to-sage-translation-map.md` to recover the full contract; this document does not narrow it to the current repair.
Read `AGENTS.md` and `references/port-design.md`, including concrete implementation assignments, before editing.
Preserve the pinned algorithm's branches, invariants, enumeration scope, termination and certificates.
No replacement algorithms or basis-dependent shortcuts.
Read called helpers as well as the named routine.

### User decisions and execution constraints

- Work in `/home/dzack/gitclones/sage-indefinite-port` and, for the authorized upstream repair interval, only `/home/dzack/research/.tmp/wt-port`. Leave `/home/dzack/research` main's dirty files alone.
  They belong to another writer.

- Stay within these assigned trees.
  Other workers and their trees are not this task.

- Obtain missing pinned reference source yourself.
  The user is not responsible for supplying an alternative algorithm or locating the upstream implementation.

- The pinned references are `Indefinite.jl@374a5ebb` and `polyhedral_common@1592b246`; Allcock edgewalk also uses the vendored paper `references/vendor/arxiv/2110.03784/`. Never build the C++ reference.

- The frozen acceptance suite is the only port suite.
  No added, edited or weakened cases.
  Only capability xfail removal is permitted once its cases pass.

- Run one Sage mathematical test process at a time.
  Long runs must be detached inside `sage-dev-develop`, with output redirected inside the container.

- Substitute the research worktree's `src` in PYTHONPATH. Do not test against `/home/dzack/research/src` or edit `sys.path`.

- Keep the original pending `marked_lattices.py` readiness changes uncommitted until all their cases pass.
  No push was requested.

## 2. Completion conditions

Completion requires source-faithful implementation and the full original frozen consumer obligation at the final source revision.
Older diagnostic evidence is useful for localization but does not establish acceptance of subsequently edited code.
No stubs were permitted by the user.
Planned frozen preamble calls may be authored before their provider exists, but missing runtime behavior remains unresolved work.

| Required family | Acceptance condition | Verified status / remaining obligation |
| --- | --- | --- |
| Boundary quotients, 119 cases | All frozen cases through actual production paths, including Dawes, Enriques, D-tower and degree-2/4 K3 | Full rerun required. Current continuation's boundary run completed no cases before interruption. Recheck coinvariant extension, selected enumeration and the port's earlier Eichler `[-1]` repair. |
| Centralizers, 5 cases | All discriminant-image results and the complete #18 consumer contract | Five-case evidence exists at earlier source; final-source rerun and readiness release remain. See evidence pointer below. |
| Vector stabilizers, 17 cases | Full file, including r_-1 through r_8, OSCAR and E10 squares 0/2/4/6/8/10 | Earlier-source reflection evidence exists; E10 squares 2/4/6/8/10 and full final-source file remain unaccepted. Square 2 was interrupted, not passed. |
| Orthogonal groups | Full frozen file; LORENTZ_OrthogonalGroup0 with all 536 generators checked against recorded images and any explicitly sourced runtime requirement | Verify generators and the complete file. The historical 30-minute budget needs its authoritative source; the 2153-second timing alone does not establish failure. |
| Vector orbits | Full frozen file, including every E10 case and the higher-index/spacelike obligations | E10 timeout problem unresolved. Trace the full pinned dispatcher for nonpositive norms; a restriction in one Lorentzian helper does not establish project-wide absence. |
| Isometry, 158 cases | Full frozen file, including the applicable 145 equivalences and 3 stabilizer datasets | No complete accepted run. Earlier killed run retained no usable full failure diagnosis. |
| Isotropic orbits, 8832 cases | All frozen cases with required line/plane/flag assembly | No complete run; ownership repair below is still required. |
| #15 perfect domains | Cases 24/26/30 and all 40×2 counts with full frozen semantic assertions | Unaccepted; execute full obligation. |
| #15 reflective corpus | All 17,526 cases and all applicable equivalence/stabilizer datasets | Unaccepted; execute full obligation. |
| Frozen operation ownership | Research owns independent generic computations; port owns specialized search, image discovery, recursion and assembly | The assigned extension/kernel/lift calls are present in `recursive.py`; complete source correspondence and consumer verification remain. Presence or passing specimens does not discharge this. |
| Readiness and xfail removal | Every applicable unfinished/xfail marker removed only after its entire obligation passes | Unfinished. Centralizer marker removal is pending locally; original handoff counted 46 markers. Recount against final source. |
| #36 and terminal gates | Satisfy #36's original type-checking obligation and original final verification contracts | Still open. While #36 is open, follow the exact documented deferred-mypy push exception; no broader gate bypass. Existing `[unverified]` obligations remain outstanding. |
| Other leaves rooted at #4 | Every leaf's original obligation and acceptance, with dependency order respected | Reconcile current issue tree and live claim before authoring. Do not interpret this table as permission to omit another leaf. |

For diagnosis only, earlier-source evidence is available in port `.tmp/`: `resume-r0-integrated.log`, `resume-reflection-r_1.log` through `resume-reflection-r_8.log`, `resume-r-minus1-preimage.log`, `resume-centralizers.log`, `resume-stabilizers-oscar.log`, and `resume-stabilizers-e10-0.log`. The corresponding source manifests, where present, bind only those runs; `resume-centralizers-source.sha256` and `resume-r-minus1-preimage-source.sha256` are not manifests for the current research source.
Do not overwrite old evidence or restart the old stabilizer runner as proof of current-source acceptance.

## 3. Remaining work and historical diagnoses

### Provider dependency: verify the historical exact finite-image failure

Select this repair only if a current consumer still needs it, and only within an authorized exclusive research interval.
Reproduce the named failure before editing.
The historical assigned worktree was `port/after-410`; verify its current source.
The following were its pending files at handoff:

- `src/dzack_research/preamble/categories/isotropic_lifts.py`: integral residue classes use the quotient's explicit `.elements().filtered(...)`; the class×period inverse now has the integral lift locus as its domain.

- `src/dzack_research/preamble/categories/sets/set_categories.py`: image membership checks admission to an owned inverse morphism's domain before applying it, preserving undecided admission as a proposition.

These changes are not accepted as a finished provider repair.
The decisive existing regression is in `.tmp/resume-lift-set-regressions.log`: the image of literal `{1,2,3}` under squaring cannot compute its cardinality.
Reproduce both existing nodes:

```
tests/sets/test_injections_surjections_and_subsets_of_finite_sets.sage::test_injection_inverse_on_image_and_characteristic_are_compatible
tests/constructions/test_sets_construct.sage::test_subsets_by_condition_and_by_image
```

The failure is in `_ImageSetEngine` around lines 3817–4040. `_FiniteLiteralSet` retains explicit finite points and iteration, but it is not a selected ordinal `EnumeratedSets()` realization.
`_finite_image_values`, cardinality, injectivity, iteration and membership currently restrict finite computation to that category.
The condition-set engine already recognizes finite literals.
Repair the generic finite-literal route consistently, without assigning a fake ranking, trying to iterate arbitrary merely-finite sets, or treating undecided equality as false.
The inverse guard is not exercised by the failing no-inverse image.

Then verify both integral-lift regimes.
Diagnostic probe `.tmp/probe/resume-reduction-lift-locus.py` checks a nonempty locus, exact integral restriction, the marked vector, and rejection of nonintegral/incorrect lifts.
`.tmp/probe/resume-empty-lift-locus.py` is prepared but unrun and checks an empty integral locus between U and U(2). These probes do not replace frozen acceptance.
The nonempty probe's evidence is `resume-reduction-lift-locus-after-2.log`. Regenerate required preamble references after source settles and use normal gates.
The owning research issue is #401; independent provider work does not reclassify or block the port obligation.
A red gate requires diagnosis before further authoring.

### Algorithm unit: verify and complete the source-derived port assembly

Read `src/sage_indefinite_port/indefinite/recursive.py`, especially `_IsotropicExtension`, `isotropic_stabilizer_action`, `isotropic_subspace_isometry`, `_isotropic_pointed_isometry`, and the split isotropic shortcut.
The current route already calls the assigned research computations; verify that it fully replaces the copied generic calculations.
The specialized recursion, Hom(K,I) shears, compatibility actions, image discovery and assembly stay here.

The source/callee correspondence and bound proof are recorded on the owning issue: https://github.com/dzackgarza/sage-indefinite-port/issues/18#issuecomment-6217015400 Read it before implementing; it is a mathematical obligation, not completion evidence.
The local copy is `.tmp/probe/extension-callee-correspondence.md`.

Exact anchors: pinned `Lorentzian.g:115–249`, `Fundamental.g:471–516,751–794`, and `IndefiniteForms.g:988–1055,1605–1650`. The rational complement is the Euclidean nullspace specified by the helper.
Research's coisotropic Witt extension is available but needs final consumer verification.

Required source correspondence and consumer checks (repair actual omissions only):

1. Verify the existing partial perpendicular isometry from owned morphisms and prescribed source images; call its `witt_extension`. Keep raw engine coordinates inside the tensor owner.
   Do not substitute an arbitrary complement.

2. Verify consumption of `R.pointwise_perpendicular_kernel().group_generators()` for the complete homogeneous direction module.
   For the source's period d, convert `1+T_j` to `1+T_j/d` through owned rational morphisms.
   The proof uses T_i T_j = 0 because each direction kills N and has image in I contained in N. Keep Hom(K,I) shears.

3. In the existing pointed isotropic Isom route, verify primitive markings, recursive reduction Isom, `R.rational_lifts(S,reduced).integral_members()`, and the appropriate integral or rational base point.
   An empty integral fibre for one reduced map does not imply the whole Isom is empty: retain the full target correction assembly.
   Preserve upstream's global-sign normalization against actual marked vectors.
   Repair the general route rather than relying on the even divisibility-one shortcut.

4. Retain the certified invariant-module bound D = 2 d² e² a. Here a is the intrinsic pairing index |coker(L -> I*)|, not the determinant from a coordinate-dependent complement; e is the reduction discriminant magnitude.
   The issue comment proves the denominator bound for line, full plane and flag actions.
   Construct the pairing map using owned dualization and `metric_map`; verify actual contracts before calling them.
   Keep the full orbit-span and containment certificate, with no truncation or heuristic bound.

Also verify `IsotropicInduction.line_automorphism_generators`, `plane_automorphism_generators`, `flag_automorphism_generators`, and the isotropic branch of `EichlerRecursion.pointed_automorphism_generators` against their concrete port-design assignments.
A call to an unresolved parabolic group does not implement them.

### Validation: complete consumers and release justified markers

Run vector stabilizers, boundaries and centralizers against the settled provider and repaired port, one file/case process at a time, retaining exact failures.
Continue through the entire completion table, including orthogonal groups, vector orbits, isometry, isotropic orbits and both #15 corpora.
Diagnose discrepancies at the precise pinned source step; acceptance is not restricted to cases already fast.

The port has two pending tracked files: `src/sage_indefinite_port/indefinite/marked_lattices.py` and `tests/acceptance/test_centralizers.sage`. The former removes three #18 unfinished decorators; the latter removes the centralizer xfail and its unused import.
Commit these only with the verified readiness they claim.
Preserve unrelated untracked `.agents`, `.hermes`, `=2.18`, `cysignals_crash_logs/`, and `tests/.report/`.

For E10 performance, inspect `.tmp/probe/lorog0.log` and profile the current path before assigning a cause.
Relevant generic computations include family sections, fixed-map Smith preimages, isometry equality and finite-action word sections.
The final timed consumer and generator images decide whether performance is fixed.

### Execution reference

The Sage container is `sage-dev-develop`, root, executable `/usr/bin/sage`. Docker is installed via Snap at `/snap/bin/docker` on this host; invoke it through `/usr/bin/sg docker -c '/snap/bin/docker ...'`. The existing `sage-dev-develop` container is running.
Source `.envrc` for the correct Sage, Python, preparser, and Julia bridge launchers.
For port tests use:

```
PYTHONPATH=/home/dzack/ai-review-ci/tool-artifacts/pytest_plugins:/home/dzack/gitclones/sage-indefinite-port/src:/home/dzack/gitclones/sage-indefinite-port:/home/dzack/research/.tmp/wt-port/src:/home/dzack/research/.venv/lib/python3.14/site-packages
```

Launch `docker exec -d -u root -w /home/dzack/gitclones/sage-indefinite-port sage-dev-develop sh -c '...'`, with the export and log redirection inside the shell.
The test invocation inside that shell is:

```
/usr/bin/sage -python -c 'import sage.all,sys,pytest; sys.exit(pytest.main(sys.argv[1:]))' -p qc_sage_session NODE --runxfail -vv --tb=long
```

Use a full file or the full inventory for `NODE`, without fail-fast; for a known failure select its exact node first.
Append an explicit exit code to each log.
A diagnostic cutoff is optional and must not become an acceptance budget or a reason to optimize without evidence.
A timeout is not a mathematical failure.
Research diagnostic regressions used `--no-time-gates`; this does not discharge normal time gates.
Check your own prior job's exit before launching another; leave unrelated processes alone.

Normal commit hooks remain required.
For research, put `TMPDIR=/home/dzack/research/.tmp/wt-port/.tmp` inside the `sg` shell so container preparse can see it.
Read research's contribution contract before committing.
Machine-wide QC serialization can wait tens of minutes; it is not permission to bypass a hook.
No broad research-suite success is established.

### Retired directions

The direct qfisom repair is no longer an untested proposal: inspect the current scaled-dual implementation rather than restarting the old mass/group-based emptiness path.
Similarly, assess the current reduced-dual full-group and cached fixed-map preimage implementations before reproducing obsolete diagnoses.
Their source and guarantees still require final consumer acceptance.

Do not revive the old handoff's claim that one Lorentzian helper establishes a project-wide spacelike source gap.
Trace the pinned full dispatcher, including higher-index Eichler cover paths; obtain sparse missing source rather than inventing a replacement.
Exact pinned C++ headers are available for reading in `.tmp/probe/pinned-combined-algorithms.h`, `pinned-lorentzian-perfect.h`, and `pinned-approximate-models.h`. Reference source is never compiled.

This file is the sole authoritative handoff.
Scratch ledgers and old runner logs are evidence pointers, not alternative work plans or completion claims.
