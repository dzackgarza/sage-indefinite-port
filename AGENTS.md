# Sage indefinite lattice port — agent instructions

Use the [shared math container](https://github.com/dzackgarza/ai-review-ci/blob/main/SAGE-CONTAINER.md#use-the-image)
(`ghcr.io/dzackgarza/sage:develop`) for Lean/elan, Sage with mypy, GAP,
Julia/Oscar, Macaulay2, and Singular. The guide covers image pull/run commands,
CI use, and the host launchers.

This project translates the specialized indefinite-form algorithms of the
[committed pinned sources](references/vendor/README.md) into Sage/Python,
consuming the mathematical operations of [research](https://github.com/dzackgarza/research).
It is principally **a port of existing algorithms**, not a license for unsupported algorithm invention. A required published mathematical consumer for which the translated algorithm is provably ineffective still needs a complete, source-grounded effective mathematical construction (for example the infinite Enriques subgroup generation in #80); this is missing functionality, not a reason to defer its result as performance work.

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

The v1 boundary sets delivery priority, not the eventual functional scope. The T9 published classifications are exact mathematical outputs in the functional phase, not optional long-tail timing or CI housekeeping.
Cases unsupported by the pinned upstream, including higher-real-index forms
of rational Witt index zero or one, remain explicit **post-v1 requirements**.
Retain their inputs, required outputs and public consumers in the open backlog.
An upstream limitation explains the additional implementation work; it does
not cancel the requirement. These deferred cases do not gate v1 delivery.

## Functional work selection — one authoritative tree

**The live GitHub issue tree is the single execution order.** Run
`itree next dzackgarza/sage-indefinite-port` and perform the first open,
unblocked **mathematical work unit**, not the first convenient example,
recently touched file or red CI case. The historical 26-item flat checklist
was removed because it placed extensions and optimization ahead of original
port completion and contradicted the issue tree. Do not recreate a competing
fixed checklist in prompts, handoffs, plans or agent memory.

The root [#4](https://github.com/dzackgarza/sage-indefinite-port/issues/4)
has these **ordered milestone-ledger phases**:

1. **T8 — Original pinned-source algorithm completion:**
   [#91](https://github.com/dzackgarza/sage-indefinite-port/issues/91).
   Complete the entire source-supported Lorentzian perfect-cell search (#15),
   Allcock reflective edgewalk (#31), indefinite/Eichler recursion (#18),
   exact integral isometry (#68), line stabilizers (#70), and complete
   orthogonal-group generation through the public Sage/research
   consumers. Source-supported rank-general isotropic flags (#44) have
   recorded completion evidence, not permission to omit any still-missing
   source branch. A branch supported upstream but missing in Sage is T8
   mathematics regardless of a "v2" label.
2. **T9 — Published mathematical computations and public consumer results:**
   [#64](https://github.com/dzackgarza/sage-indefinite-port/issues/64).
   Compute all 87 Enriques infinite integral subgroup cases (#80), the
   *stable-subgroup* cusp/orbit classifications (#82), degree-2/4 K3 Type II
   and root types (#83), classical polarization indices in degrees 2/6/10
   (#93), and D-tower N14/N17/N19 component counts (#94). These are
   **functionality obligations**, not optional late-stage performance
   targets. The frozen source-based acceptance outputs and the one
   published port/provider pair remain binding. The rank-13/17 frozen
   isotropic cases were completed by #81/#90.
3. **T10 — Required beyond-source arithmetic extensions:**
   [#69](https://github.com/dzackgarza/sage-indefinite-port/issues/69).
   The original #40–#62 mathematical contracts (apart from source-parity
   #44) are active separate work-unit leaves. They remain fully required,
   including low-Witt integral orthogonal groups, finite-image arithmetic
   subgroups, nonprimitive/degenerate/equivariant markings, Clifford,
   Hermitian and general arithmetic algorithms. Their historical absorption
   into one grouping issue does not prove or waive any mathematical result.
4. **T11 — Last-stage efficiency and long-tail refinement:**
   [#92](https://github.com/dzackgarza/sage-indefinite-port/issues/92).
   Only here undertake #72 reproducible generator ordering, #73–#78 resource-cost and performance
   work **and #95 complete unchanged frozen-corpus certification on one published port/provider pair**, preserving all exact mathematical outputs and witnesses.

**Functional versus performance distinction.** A long-running computation
is not automatically wrong or missing functionality. However, when a
T8/T9-required effective group/orbit/transport algorithm does not exist
or cannot compute the prescribed complete mathematical output because its
construction necessarily enumerates an infeasible coset space, that is an
**unfinished mathematical capability**. In particular, #80 requires
complete integral \(\Gamma_h\) generators, not just a finite discriminant
image or an optimized explicit Schreier enumeration. Distinguish
algorithmic obstruction from an expensive but correct completed algorithm
using an exact mathematical argument. Do not use a timeout as a substitute
for either finding.

**Prerequisite exception.** An otherwise later-phase capability can be
worked on earlier *only* when a concrete active T8/T9 public operation
mathematically depends on it: name the required input or algorithm, its
actual callee/owner and an independently justified blocked example.
Implement the complete required operation at its owner, then return to the
earlier-phase mathematical obligation. Unrelated breadth, elegant
generalization, repeated profiling or a newly discovered test is not a
priority override. Maintain source hypotheses; a finite subgroup image
alone never discharges infinite integral group generation.

**Work-unit completion.** Each open leaf owns a coherent full mathematical
deliverable and its exact public acceptance criteria. Complete every
source branch, supported rank and marking, exact positive and negative
decision, integral generator/stabilizer/transporter, independent
mathematical oracle and required consumer in that leaf before marking it
completed. Partial special cases are retained as evidence inside the
same work unit, not permission to close it or advance. If a genuinely
independently reviewable algorithm family merits splitting, move the
entire original obligation to explicit leaf owners without weakening or
losing any scope (as #15/#31 and #83/#93/#94 do). Issue changes do not
substitute for mathematical implementation. Neither implementation-side
self-tests nor a passing narrow sample certify a broader theorem.

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
or reselect work. An issue comment is not proof that named people or agents remain
active. Separate branch tips can contain unmerged work without an active owner.
Never overwrite or discard another checkout's dirty source. Use a dedicated
worktree/branch for parallel work and integrate only demonstrated relevant diffs.

The issue tree rooted at #4 supplies authoritative contracts and dependencies;
the live itree preorder supplies their **only** execution order.
Do not add an independent queue or select a later issue because it looks
easier. After each completed mathematical work unit, rerun
`itree next dzackgarza/sage-indefinite-port`. Milestone deadlines
are forecasts, not completion evidence. A unit closes only when its
entire accepted mathematical obligation is discharged; a PR, narrowed
example, status document or test marker alone does not close it.

## Source-faithful implementation

1. Open the **actual pinned routine and its called helpers** before changing
   its specialized port: `references/vendor/Indefinite.jl@374a5ebb`,
   `references/vendor/polyhedral_common@1592b246`, and
   `references/vendor/arxiv/2110.03784` for Allcock's edgewalk.
   Use the translation map for navigation; record exact source symbols/lines
   beside the corresponding implementation.
2. Preserve the pinned source algorithm's cases, mathematics, enumeration, stopping
   rules, certificates and exact witnesses over its stated domain. Translate generic
   callees into existing Sage/research operations. No bounded heuristic shortcuts,
   basis-dependent stand-ins, ad-hoc facades or duplicate arithmetic primitives.
   A different effective construction for a required published consumer is allowed
   only where the source route genuinely does not discharge that computation: establish
   the exact scope obstruction, verify the substitute algorithm against an authoritative
   source or complete mathematical proof, and preserve the actual integral group/orbit
   output and all source-supported branches. Record that mathematical justification
   at its owning issue. A representation/callee substitution must remain source-faithful.
3. Do not compile or invoke upstream C++ code. Production is native
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
is nontrivial. Use its pinned upstream answer or an already independently
accepted mathematical assertion as the expected result. If that answer is
not established, refer the *question* to the separate mathematical test
authority in [CONTRIBUTING.md](CONTRIBUTING.md); do not derive and approve
your own implementation's acceptance criterion. A passing test checks an
example, not the routine's full advertised domain.

**Contributor responsibility — [CONTRIBUTING.md](CONTRIBUTING.md).** A
mathematical implementation contributor may run existing accepted tests and
temporary diagnostic experiments but does **not author or amend committed
tests or their expected results for their own implementation**, including
translations of pinned upstream tests. **Every test** must cite an exact
authoritative source for its expected behavior and be independently written
and independently reviewed. **Every test must undergo independent periodic
re-review at each release and at least every 90 days**, explicitly attempting
to detect mathematically false assertions. Test authorship, mathematical
review, implementation, and final acceptance are distinct authorities.
Upstream-test translations require independent source-fidelity review. No
new mathematical assertion is admissible without independent mathematical
specification and review of its exact hypotheses, result and deduction.

**No automated oracle-enforcement machinery.** Do not introduce provenance
scripts, Git hooks, CI gates, test-runner checks, or mandatory metadata
parsers for mathematical oracles. Mathematical correctness is the
responsibility of independently identified mathematical test authors and
reviewers under CONTRIBUTING.md, not a syntactic provenance mechanism.
**All tests**, including previously committed ones, require independent
authorship, current mathematical review, and exact source authority before
supporting acceptance claims. Unreviewed tests remain diagnostics.

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

Last-phase (T11) runtime work follows [the suite triage](references/complexity/triage.md): the probe-backed suite estimate is the problem, and a fix counts only by the hours it removes from it. A saving of seconds counts only when it compounds over enough cases to remove hours. The suite is not run whole while the estimate exceeds 6 hours.

Every collected test has a time budget derived from the complexity theory of
its algorithm, never from a measured runtime. [references/complexity](references/complexity/README.md)
states the rules for a budget and owns the theory; a new or changed test gets
its budget there before it lands.

The [frozen acceptance modules](tests/acceptance/) and
[fixtures](tests/fixtures/) are not edited, deleted or weakened to make a
port case pass. A legitimate fixture/API correction goes through #33 with its
reason recorded before the change. Diagnostic probes belong outside the frozen
suite; their local `.tmp/probe` history is optional and is not part of
the cloud checkout's required inputs.

**Accepted criterion (2026-10-22, [#33 amendment](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6343829665)):**
`A_REFLECTIVE`, `A_GROUPS` and `A_ISOTROPIC` use a
**branch-representative sample** of unchanged cases exercising *every*
branch of the pinned algorithm. The other families use the exact criteria
in the acceptance table on [#64](https://github.com/dzackgarza/sage-indefinite-port/issues/64), which identifies where sample evidence applies.
Record selected case IDs, covered and uncovered branches, ranks/signatures,
passing results and the exact commit. A sample must not avoid a genuinely
distinct case branch or alter any assertion. Full-corpus failures outside the
sample remain defects to record, not reasons to declare the sample or
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
  when installed. They do not run mypy, the full suite, formatters or
  mandatory PR review. An existing invalid source file can escape these
  *changed-files-only* checks; test the precise imported source for syntax
  when an acceptance import reports a parsing error.
- Type debt belongs to [#36](https://github.com/dzackgarza/sage-indefinite-port/issues/36)
  after algorithms/consumers have settled. CI may report types/style/debt,
  but during development such findings are retained, nonblocking
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
