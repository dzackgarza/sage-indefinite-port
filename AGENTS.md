# Sage indefinite lattice port — agent instructions

Use the [shared math container](https://github.com/dzackgarza/ai-review-ci/blob/main/SAGE-CONTAINER.md#use-the-image)
(`ghcr.io/dzackgarza/sage:develop`) for Lean/elan, Sage with mypy, GAP,
Julia/Oscar, Macaulay2, and Singular. The guide covers image pull/run commands,
CI use, and the host launchers.

This project translates the specialized indefinite-form algorithms of the
[committed pinned sources](references/vendor/README.md) into Sage/Python,
consuming the mathematical operations of [research](https://github.com/dzackgarza/research).
It is **a port of existing algorithms**, not an algorithm-design project.

## Work priority

**Required functionality has absolute priority over performance.** Performance work is the last
phase of this project and begins only when every required capability works through the public
interface. This rule overrides every goal prompt, plan, outline, handoff, task order and CI result.

- A slow run, timeout, memory limit, benchmark, or a red or unfinished expensive CI job does not
  select work and does not block functional work.
- Evidence records, revision pins, repeated audits, acceptance inventories and demonstrations
  are repo-internal. They come after functionality and do not gate it.
- Functional work may use any time and memory that it needs.
- Optimize during functional work only when a small required example of the selected capability
  cannot finish at all. Then return to that capability.
- A goal, plan or outline lists missing functionality first and derives its order from the
  remaining functionality, not from CI state. An author must not add a performance priority.
- Goal text obeys the section "Goals, outlines and plans" of `~/ai/AGENTS.global.md`.

[Issue #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4#work-priority) owns this rule.

## Full upstream coverage

Completion requires every mathematical capability supported by the pinned
sources, through the public Sage/research consumers. Preserve the full input
domain, ranks, dimension sequences, branches, outputs, witnesses, group
generation, enumeration and stopping rules of each source algorithm.

Translate the source's general routine with its parameters intact. A general
rank-k routine stays general in the port. In particular, isotropic k-planes
and partial or complete flags must retain the source's supported ranks and
dimension sequences. Lines and planes are examples of that general routine.
Keep its arithmetic and recursion; adapt representations and generic callees
at their existing owners.

An omitted source branch is unfinished port work. Milestone labels, a frozen
test corpus, application-specific examples, provider ownership and publication
of a partial implementation do not reduce the required functionality.
Source-supported functionality belongs to the port's completion condition,
including functionality tracked by an issue titled v2 or post-v1. Additional
requirements beyond the source remain required under their own contracts.

Before accepting an implementation or completion claim, map the source's
dispatch branches and parameter domain to the production route and its
mathematical evidence. Record every uncovered branch on the issue that owns its
algorithm in the tree rooted at #4. Missing source functionality keeps that algorithm incomplete.
A rank-three isotropic sublattice accepted upstream but rejected by the port
is a missing translation, even when every line and plane example passes.

## Current functional MVP phase

The user's 2026-10-22 MVP directive supersedes older completion and refinement
gates below. Complete every algorithm family and its real public consumers;
verify each implemented mathematical branch on small, independently justified
examples. Publish coherent functional changes without waiting for the complete
acceptance matrix. Preserve the full frozen CI suite as nonblocking long-tail
evidence. Performance, typing (#36), formatting, refactoring, exhaustive
coverage and repeated audits are deferred until functional MVP delivery.
Repair only defects that block the required mathematical computations or their
small correctness examples. Functional publication may be incremental; port
completion requires the full upstream coverage specified above.

The v1 boundary sets delivery priority, not the eventual functional scope.
Cases unsupported by the pinned upstream, including higher-real-index forms
of rational Witt index zero or one, remain explicit **post-v1 requirements**.
Retain their inputs, required outputs and public consumers in the open backlog.
An upstream limitation explains the additional implementation work; it does
not cancel the requirement. These deferred cases do not gate v1 delivery.

## Functional work selection

**Implement the user's entire sage-indefinite-port assignment, not a series of
convenient examples.** There is exactly one committed, well-ordered functional
queue: the numbered obligations below. It refines the order required by #4
and #69; it does not replace their full mathematical contracts. `itree next`
currently selects #69; #44 was closed on 19 November 2026. An obsolete
handoff, issue comment, agent goal, earlier generated plan, test result or
latest local discovery cannot select a different task.

**Selection is deterministic:** work on the *first unchecked number* below.
Complete every mathematical output, input domain, rank, branch, integral
witness, exact group-generation/orbit/negative-decision requirement and public
consumer belonging to that number before checking it or moving to the next.
Its issue's full contract is binding, including requirements not summarized
on the queue line. A successfully implemented special case or a pushed commit
does not check the number off. Closed former issue numbers absorbed into #69
do not imply that their larger retained contracts are complete.

**A real prerequisite is the only exception to executing the first unchecked
number.** It must supply a specifically named mathematical input used by that
number, have a concrete failing public example, and belong to an existing
owner. Implement that prerequisite immediately, then return to the interrupted
number; do not advance the queue or create a parallel discretionary backlog.
Do not reassign priorities based on ease, curiosity, interesting failures,
available tests, CI duration, or which piece of code is already open.
Do not reinterpret a partially implemented algorithm as a new, smaller task
whose completion permits switching to another issue.

**Fixed functional queue (in this order, first unchecked is active):**

1. [ ] **#42, Tier 0 — uniform finite-image arithmetic subgroups.** Complete
   the *actual joint image* and exact owned integral preimages for
   discriminant actions, determinant, normalized real and base-changed field
   spinor, positive-cone component where defined, simultaneous/composite
   congruence levels, nonnormal targets and explicitly generated ambient
   subgroups. Supply full kernel/section words, subgroup generators,
   membership, cosets, marked stabilizers and transporters through the public
   `L.O()`/marked-group consumers. All factors must retain correlations;
   neither finite images alone nor separate character kernels suffice.
2. [ ] **#40/#55, Tier 0 — Witt-index-0/1 constructive orthogonal groups.**
   Implement the complete effective Li–Margulis bound and full finite
   integral generator search, exact positive isometries and certified
   negative integral isometry decisions on higher-real-index lattices of
   rational Witt index below two, without an assumed rational hyperbolic
   plane. Both preserved low-Witt specimens and public `O(L)`/`Isom(L,M)`
   must work. The existing lazy group iterator is not full generation.
3. [ ] **#48, Tier 0 — content, nonprimitivity and odd forms.** Complete all
   nonzero-/zero-norm content strata, actual primitive saturation and
   nonprimitive sublattice/flag equivalence with stabilizers, and form
   doubling with original-lattice discriminant and isometry semantics.
4. [ ] **#62/#53, Tier 0 — general orbit quotients.** Complete exhaustive
   duplicate-free canonical representatives and integral transporters on
   arbitrary finite- and infinite-index `Gamma`-sets, including HNF
   submodules, flags, Hermitian encodings and prescribed pairing slices.
   Finite bounded searches must not be advertised as complete quotients.
5. [ ] **#51/#58 (definite branches), Tier 0 — integral CM/Hermitian unitary
   groups.** Complete all source-supported definite-at-every-place and
   finite-order cyclotomic centralizers on the *original* order/module,
   including mixed-place signs, nonmaximal and nonfree integral structures,
   scalar kernels, actual discriminant images and the public group API.
6. [ ] **#52, Tier 0 — Allcock–Carlson–Toledo `E^(4,1)`.** Complete the
   exact source reflection/translation reduction, full `U` versus index-two
   `Aut+`, projective and mod-theta congruence distinctions, primitive
   Eisenstein-null/root orbit generators, stabilizers and transporters
   through original trace-lattice public consumers.
7. [ ] **#59, Tier 0 — indefinite Hermitian source algorithms.** Complete
   Bestvina–Savin's actual Mendoza-spine ocean, integral/projective lifts,
   Gaussian and all three Euclidean Picard groups, exact word reductions,
   cusp stabilizers and original-order transporters; not just group matrices.
8. [ ] **#60, Tier 0 — integral Hermitian parabolics.** Complete the full
   nonabelian unipotent radical, center and actual integral/glued Levi image,
   rank-general parabolic and partial-flag subgroups, witnesses, finite-index
   intersections and original-module public consumers.
9. [ ] **#43, Tier 1 — one full finite-double-coset algorithm for marked
   actions.** Complete subgroup-correct splitting, setwise/pointwise
   stabilizer generators and positive/negative transporters for *all* named
   arrows, vectors, sublattices, configurations and flags, including
   nonnormal finite images and arbitrary user-supplied `Gamma`.
10. [ ] **#47, Tier 1 — complement and gluing.** Complete indefinite
    integral nondegenerate embedding stabilizers and transporters, all
    complement genera, discriminant-form glue compatibility, infinite
    complement factors and exact restrictions to source-action groups.
11. [ ] **#49, Tier 1 — simultaneous commuting involutions.** Complete
    arbitrary-rank joint eigensublattice integral gluing, entire centralizer
    generators, subgroup restrictions, original-lattice witnesses and exact
    orbit actions, not independent-factor images.
12. [ ] **#41/#61, Tier 1 — degenerate and equivariant marked objects.**
    Complete all original-rank radical GL and shear actions, based/moving
    degenerate Gram domains, nonprimitive mixed frames, arbitrary source
    GL/diagram actions and full stabilizers/negative transport decisions.
13. [ ] **#46, Tier 2 — Gram-constrained configurations.** Complete
    primitive fixed-norm vector families, arbitrary-rank fixed-source
    nondegenerate embeddings, moving/unbased orbit classification,
    integral stabilizers and transporters across all complement classes.
14. [ ] **#56, Tier 2 — bounded primitive Gram embeddings.** Implement the
    Li–Margulis completion bounds and full orbit and isotropic-vector
    reduction with proofs of exhaustion and exact integral witnesses;
    consume #47 gluing rather than narrowing to rank one.
15. [ ] **#40, Tier 2 — all low-Witt vector classes.** Use the completed
    Tier-0 orthogonal group/isometry backend for fixed-norm and isotropic
    vector/line classifications, complete orbit representatives, negative
    equivalence and integral stabilizer generators in Witt indices 0/1.
16. [ ] **#45, Tier 2 — quotient rational Tits building.** Complete flags
    of all supported dimension sequences, all face maps with witnesses,
    simplicial identities, face-stabilizer actions and incidence
    multiplicities under full and specified arithmetic subgroups.
17. [ ] **#50, Tier 2 — exact normalizers.** Complete finite nonabelian
    subgroup and infinite cyclic normalizers, twisted conjugacy transporters,
    full centralizer/inversion cosets and finite-index subgroup restrictions.
18. [ ] **#57, Tier 3 — Clifford arithmetic.** Complete ranks 4, 5 and 6,
    actual quaternion/Clifford orders, reduced-norm-one and full unit
    extraction, rational/integral isometry decisions and full group outputs.
19. [ ] **#51/#58, Tier 3 — remaining noncompact arithmetic unitary and
    commutant-order branches.** Complete all unsupported number fields,
    higher rank/anisotropic Hermitian forms, nonfree modules, centralizer
    groups and exact equivariant transporters with original integral action.
20. [ ] **#54, Tier 4 — general Grunewald–Segal arithmetic fallback.**
    Complete A1/A2 effective group generation, arbitrary simultaneous
    centralizers and tensor/marking stabilizers, exact positive and
    negative integral transporters without semisimplicity or Witt
    restrictions, and finite-index arithmetic subgroup intersections.
21. [ ] **#15 — source-complete Lorentzian/perfect-cell/Allcock algorithms.**
    Complete missing face/edgewalk branches, all prescribed source leaves,
    exact reflective domain and finite-rank public consumers.
22. [ ] **#18 — source-complete indefinite recursion.** Complete all pinned
    source dispatches, general isotropic ranks/partial flags, integral
    stabilizers and branch-wide public orbit/isometry consumers.
23. [ ] **#68 — spinor-genus integral isometry.** Complete genuine
    integral non-isometry decisions and positive witnesses beyond the
    Conway–Sloane example, at the original public `Isom` carrier.
24. [ ] **#72 — deterministic complete automorphism generation.** Retain
    exact generated groups and source semantics without unseeded
    process-dependent algorithmic outputs.
25. [ ] **#70/#73–#77 — required last-phase optimizations.** Only after
    full mathematical functionality: improve the source-exact integral
    stabilizer/group and arithmetic heavy cases without losing kernels.
26. [ ] **#64/#78/#80–#83 — final common-revision acceptance.** Complete
    the *unchanged* published full corpus, stable `O+(N)` orbit counts,
    Enriques 87 modular groups, K3 Type-II and related boundary cases,
    all at one actual released port/provider revision. A test timeout,
    fixture amendment or partial published family is not completion.

**Checkpoint discipline.** At each commit, leave the currently active
number unchecked unless *its entire mathematical contract* has been
implemented through the public interface; completed subcases stay part of
that same number. Only check off a number in a commit containing the actual
finished mathematical implementation and independent public checks (not a
planning-only commit), then proceed to the next number. If source evidence
shows that an obligation was fully implemented already, inspect the actual
public path and the full domain before checking it. If a genuine mathematical
dependency is discovered, implement it in place and record its necessary
input relationship, without silently reordering this list. All 26 boxes,
not just the next one, must be discharged to finish the user's entire task.

A prerequisite must name the dependent capability, the exact operation it
consumes, and a small independently justified example that cannot complete
without the repair. A shared module or a timeout in a larger specimen does not
establish that dependency. Repair the smallest blocking operation, verify that
example, then return to the dependent implementation.

Issues #64 and #70 retain expensive-case and final-release obligations. Schedule
performance work there after functional coverage, except for a demonstrated
blocking prerequisite above. A terminal expensive run contributes its result to
that issue; it does not automatically select another optimization or rerun.
Choose subsequent work by missing functionality. Full release acceptance remains
required after implementation. Publication uses the existing small branch checks.

## Sources of authority

| Source | Responsibility |
| --- | --- |
| [Issue #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4), [#33](https://github.com/dzackgarza/sage-indefinite-port/issues/33), [#15](https://github.com/dzackgarza/sage-indefinite-port/issues/15), [#18](https://github.com/dzackgarza/sage-indefinite-port/issues/18), [#36](https://github.com/dzackgarza/sage-indefinite-port/issues/36) | Original obligation, ownership, issue dependencies and acceptance |
| [Issue tree rooted at #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4), navigated with `uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port` | Remaining algorithm, dependency and acceptance obligations, one work unit per node; a closed node needs source-backed evidence, and the next action lives on the open node |
| [Port design](references/port-design.md), [assignments](references/port-design.md#implementation-assignments), [translation map](references/source-to-sage-translation-map.md) | Frozen API, owner split, exact upstream routine locations |
| [Vendored source and Allcock paper](references/vendor/README.md) | **Algorithm authority**; translation notes never replace source |
| [Acceptance workflow](.github/workflows/acceptance.yml), [oracle manifest](tests/fixtures/oracle_manifest.yaml) | Frozen executable cases and results for exact published revisions |

## Entry and continuity

For a **new cloud checkout**, work from `main`, where
implementation and acceptance now land. Read this document once, then
the next open work unit from `itree next` and the relevant
issue/source sections. Fetch the branch and identify its current commit and
newest *completed* Acceptance run before changing it. If the issue's recorded commit or
selected failure predates newer code/results, record the actual next
unresolved implementation step on the issue; do not enact its stale instructions.

A normal continuation starts with the next open work unit and its precise
source/consumer, rather than repeat source-wide surveys, re-read the entire DAG,
or reselect work. An issue comment is **not** proof that named people or agents remain
active. Separate branch tips can contain unmerged work without an active owner.
Never overwrite or discard another checkout's dirty source. Use a dedicated
worktree/branch for parallel work and integrate only demonstrated relevant diffs.

The issue tree rooted at #4 supplies authoritative contracts and dependencies;
the *single numbered queue above* supplies their fixed execution order.
Do not replace it with a fresh agent plan or select a later issue because its
example seems simpler. #69 is the first open itree work unit; do not
reselect until its active obligations have been completely discharged.
Milestone deadlines are forecasts, not completion evidence. Use
`uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port`
only to confirm issue-tree ownership after completing the active queue unit,
never to reshuffle the numbered list. A unit closes only when its entire accepted obligation is discharged;
releasing a PR, a narrowed example, or a marker alone does not close it.

## Source-faithful implementation

1. Open the **actual pinned routine and its called helpers** before changing
   its specialized port: `references/vendor/Indefinite.jl@374a5ebb`,
   `references/vendor/polyhedral_common@1592b246`, and
   `references/vendor/arxiv/2110.03784` for Allcock's edgewalk.
   Use the translation map for navigation; record exact source symbols/lines
   beside the corresponding implementation.
2. Preserve the source algorithm's cases, mathematics, enumeration, stopping
   rules, certificates and exact witnesses. Translate generic callees into
   existing Sage/research operations. No alternative algorithms, bounded
   heuristic shortcuts, basis-dependent stand-ins, ad-hoc facades or duplicate
   arithmetic primitives in place of the specialized source route.
   Account explicitly for an unavoidable representation/callee substitution
   against its source contract. A mathematical departure needs a pinned-source
   proof or counterexample on the owning issue.
3. **Do not compile or invoke** upstream C++ code. Production is native
   Sage/Python with exact standard backends and the formed-lattice preamble.
   Do not bypass the owned lattice, morphism, homset or tensor interfaces with
   raw matrices. Never introduce `Any` escapes, forged stubs,
   `sys.path` changes or weakened exceptions to silence failures.
4. The frozen public mathematical objects are the research preamble's lattice,
   element, morphism and tensor objects. Pointed markings are fixed in the
   coslice, moving sublattice/configuration markings use the arrow category,
   and flags use `Fun([2], Lattices)`. Their Isom/Aut questions differ.
   Preserve those distinctions and their complete construction data.

**Division of work.** `research` owns generic kernels, fibres, equation and
module arithmetic, partial-isometry extension, rational/integral lift loci,
finite configurations/actions from supplied data, group sections and general
gluing. This port owns indefinite/Lorentzian search, perfect-domain flips and
quotients, Allcock cycles, approximate and Eichler recursion, parabolic
compatibility-image discovery, orbit classification, and specialized group
assembly. The concrete assignment is in
[port design](references/port-design.md#implementation-assignments).
Calling an unresolved generic `.group_generators()` is not an implementation
of a specialized parabolic algorithm.

A defect in a generic research operation is reported at its research issue
with the input contract, failed consumer and exact provider revision.
Continue independent port work against the frozen API rather than changing
research source in the port worktree. Do not claim an upstream research
operation is correct merely because its public name exists.

## Mathematical implementation feedback loop

**The inner loop is one mathematical obligation, not a test module, test
family, CI tier, or acceptance suite.** Identify the pinned source operation,
its exact hypotheses, and the smallest input on which the proposed behavior
is nontrivial. Establish the expected answer independently: prove the orbit
count, subgroup equality, integral obstruction, termination condition or
transporter property as applicable. A passing test checks that example; it
does **not** substitute for a proof of the routine's full advertised domain.

Run one explicitly named pytest node with
`just test-node 'tests/path.py::test_exact_mathematical_claim'`. The recipe
rejects bare files and selections expanding to more than one collected test,
including an unqualified parametrized function; specify its exact parameter ID. Prefer a directly executable
one-case witness of the actual public operation over a large imported
fixture or an expensive integration case. If the smallest valid counterexample
is expensive, keep that example; do not replace the required mathematics
with a weaker claim merely to save time.

For each source correction: establish the failure or missing branch on that
single witness; implement the full source-supported construction; rerun the
same node; then check the *few other nodes for genuinely distinct mathematical
branches* touched by the change. For generator algorithms prove the family
generates the complete group, not merely that each proposed generator
preserves the form. For orbit enumerators prove exhaustiveness and pairwise
inequivalence, not just that returned elements are legal. For transporters
verify existence/nonexistence against an independent invariant, integral
maps, both directions and the requested subgroup when applicable. For
inductive routines prove the rank-general step separately from testing its
smallest instance. Keep every essential hypothesis and do not infer
completeness from lack of further examples in a search.

**Do not run a whole test module, multi-module selection, acceptance family
or CI workflow as a routine inner-loop step.** Broader integration tests run
once when a demonstrated dependency or a completed coherent implementation
needs them; frozen full acceptance is the separate final gate. The existing
suite-complexity/budget records concern that later gate and do not determine
what to implement or what to test in the hot loop. In particular, a long
regression suite is never a prerequisite to making the next mathematical
code change.

## Tests and acceptance

Run-time work follows [the suite triage](references/complexity/triage.md): the probe-backed suite estimate is the problem, and a fix counts only by the hours it removes from it. A saving of seconds counts only when it compounds over enough cases to remove hours. The suite is not run whole while the estimate exceeds 6 hours.

Every collected test has a time budget derived from the complexity theory of
its algorithm, never from a measured runtime. [references/complexity](references/complexity/README.md)
states the rules for a budget and owns the theory; a new or changed test gets
its budget there before it lands.

The [frozen acceptance modules](tests/acceptance/) and
[fixtures](tests/fixtures/) are not edited, deleted or weakened to make a
port case pass. A legitimate fixture/API correction goes through #33 with its
reason recorded before the change. Diagnostic probes belong outside the frozen
suite; their local `.tmp/probe` history is optional and is **not** part of
the cloud checkout's required inputs.

**Accepted criterion (2026-10-22, [#33 amendment](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6343829665)):**
`A_REFLECTIVE`, `A_GROUPS` and `A_ISOTROPIC` use a
**branch-representative sample** of unchanged cases exercising *every*
branch of the pinned algorithm. The other families use the exact criteria
in the acceptance table on [#64](https://github.com/dzackgarza/sage-indefinite-port/issues/64), which identifies where sample evidence applies.
Record selected case IDs, covered and uncovered branches, ranks/signatures,
passing results and the exact commit. A sample must not avoid a genuinely
distinct case branch or alter any assertion. Full-corpus failures outside the
sample remain defects to record, **not** reasons to declare the sample or
the entire port correct. Long-tail exhaustive coverage does not gate the
agreed milestones.

GitHub [Acceptance](.github/workflows/acceptance.yml) runs only by manual
dispatch (`gh workflow run Acceptance -f corpus=tier1`, `-f corpus=tier2` or `-f corpus=full`).
A result does not select the next task. It runs the frozen
suite in sharded Sage jobs using the
`ghcr.io/dzackgarza/sage:develop` research image; the run
installs the published `research` branch and records its commit.
Only the fork's `Research environment image` workflow builds that image,
from commits pushed to `dzackgarza/sage` `develop`. Sage source is never
checked out or built on the host.
Use the exact run/artifact for the revision under discussion. Job status
alone cannot prove a case passed. A full workflow can be red while the
explicitly cited branch-representative samples pass; conversely, an
interrupted/missing shard cannot be counted as passed. The final revision
requires its chosen case IDs to pass under the same revision/provider.

In a capable local Sage environment, run selected small production consumers
with independent expected answers. These establish the exercised MVP branches;
they do not establish unexecuted frozen cases. Do not launch an entire family locally, repeatedly poll queued
CI, or substitute local machine-specific `.tmp` logs for published proof.
If the cloud runtime lacks Sage, use the committed source/fixture data and
GitHub Actions; local Snap Docker, `sage-dev-develop`, the former
`/home/dzack/research` tree and local Zotero files are not prerequisites.

## Commit, publication and diagnostics

- Commit each coherent source change, including WIP progress where required;
  identify the pinned source step or failed consumer. Push so CI tests the exact
  revision. Preserve all unrelated work and claim state.
- Publish working mathematical implementations directly to `main` when the
  repository permits it. A PR/branch containing only a queue, plan, receipt,
  evidence or status change is not a substantive task. Never replace
  completion of the first unchecked unit with a queue-administration session.
- Local `just test-commit` checks only staged Python/JSON/TOML syntax;
  `just test-push` checks only outgoing committed syntax.
  Both delegate to [ai-review-ci](https://github.com/dzackgarza/ai-review-ci)
  when installed. They **do not** run mypy, the full suite, formatters or
  mandatory PR review. An existing invalid source file can escape these
  *changed-files-only* checks; test the precise imported source for syntax
  when an acceptance import reports a parsing error.
- Type debt belongs to [#36](https://github.com/dzackgarza/sage-indefinite-port/issues/36)
  after algorithms/consumers have settled. CI may report types/style/debt,
  but during development such findings are retained, **nonblocking**
  diagnostics; they do not divert the current mathematical implementation.
- At the recorded source revision there were no active `@unfinished`
  wrappers or acceptance xfails. Do not invent marker-removal work or infer
  readiness from their absence. Public capability is judged by actual
  consumer behavior.
- Treat a red **actual local syntax check** as an immediate defect in the
  committed source. Treat other CI failures as evidence to route to their
  algorithm/provider owner; investigate concrete failing cases instead of
  performing broad cleanup, optimization or repetitive diagnostics.
- The final issue #4 closure requires the acceptance evidence on [#33](https://github.com/dzackgarza/sage-indefinite-port/issues/33) and [#64](https://github.com/dzackgarza/sage-indefinite-port/issues/64) and public
  research-to-port consumers, with a recorded port and provider revision.
  `just test-ci` is a nonblocking diagnostic collection, not a substitute
  for the frozen Acceptance workflow.

Do not require a local agent-memory vault, home-directory paths, or a running
coordinator to resume cloud work. Durable mathematical decisions belong in
the owning source/issue; the open work unit records only the next actionable
step, not an activity log.
