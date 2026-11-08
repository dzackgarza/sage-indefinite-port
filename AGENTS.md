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

Implement missing public capabilities before performance refinement and final
corpus certification. Use the functional sequence in issue #4 and the current
handoff. Start with the next concrete missing branch in #44; then implement
ready capabilities in #69 using their actual mathematical prerequisites.

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

The issue tree rooted at #4 is the work plan. Use
`uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port`
only when selecting a new unit after completion or an evidenced dependency
change. A unit closes only when its entire accepted obligation is discharged;
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

## Tests and acceptance

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
dispatch (`gh workflow run Acceptance -f corpus=sample` or `-f corpus=full`).
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
