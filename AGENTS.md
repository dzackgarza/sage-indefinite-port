# Review Guidelines

These are additional requirements for reviewing agent work.
They do not replace the reviewer’s normal role, repo-specific standards, or technical judgment.
They provide the failure model that should shape the review.

The task is not merely to review a PR. The task is to decide whether a completion claim is true under the original objective.
The standard is full, correct, provable completion against the original requirements and repo guidelines.
Anything less is incomplete work that must not be treated as a win.

## Failure Model

Agents systematically produce impressive non-completion.
Common patterns are: polished summaries that imply finished work, caveats that quietly narrow the goal, reclassification without proof, delegated discovery presented as resolution, process language that substitutes for evidence, merged PRs treated as completion, passing checks treated as semantic proof, and artifacts that look substantial while leaving required work unowned.

Treat the agent’s summary, PR description, closing comment, issue closure, “goal completed” statement, and self-reported validations as untrusted.
They may be diagnostic pointers, but they are not evidence that the work is complete.
The evidence is the original issue or task, the code diff, tests, source/runtime facts, review comments, and produced artifacts.

## Decisive Invariants

Preserve the original success condition.
Read the original issue or task before accepting any restatement of it.
Keep its quantifiers intact: “all,” “complete,” "full subset," “zero remaining,” and similar terms cannot be quietly narrowed to examples, partial coverage, known blockers, or whatever the PR happened to touch.

Nothing required may disappear silently.
A required work family must be implemented, explicitly falsified, or validly reclassified with evidence that satisfies the issue’s own standard.
Partial implementation is not completion.
Future work is not completion.
Count reduction is not completion.
Resolved review threads are not completion.
Passing checks are not completion.
Substantial-looking work is not completion.
“Better than before” is not completion.

Goal substitution is the main thing to detect.
Ask whether the submitted work solves the original problem or merely produces a narrower artifact: cleaner metadata, a partial subset, a better explanation, a new issue, a renamed scope, a local workaround, or proof that someone should investigate later.

Technically correct administrative artifacts can be goal substitution.
A well-written issue, comment, audit note, scope statement, or enumeration of remaining work may be required, but it does not complete implementation, testing, proof, or downstream cleanup.
If the original task requires execution, the artifact is only useful insofar as it drives that execution; it must not become the stopping point.

Treat self-scoped remaining-work lists as a severe completion-laundering pattern.
When an agent is asked to enumerate remaining work, the domain is the original full completion requirement, not the agent’s intended subset, the PR’s current shape, a closeability criterion, or the work left after deferral and reclassification.
A valid enumeration subtracts only artifact-proven completed work from the original contract.
Deferrals, routed follow-ups, owner changes, and truthful incompletion notes remain unresolved work unless the original task explicitly made that administrative routing the whole deliverable.

If an agent repeats a narrowed enumeration after being corrected, treat that as a hard misalignment signal, not as an innocent wording issue.
The reviewer should identify the original full requirement, the scope the agent substituted, and the required work hidden by that substitution.

Silent reclassification is not resolution.
If the PR says remaining work is out-of-scope, research-owned, stub-owned, plugin-owned, downstream-owned, or future-owned, require evidence from the relevant source/runtime behavior, repo boundary, or original acceptance criteria.
A sentence in the PR description is not enough.

Ownership boundaries matter.
The submitting repo must prove its own claimed behavior and do the blocker forensics required by its own issue.
Do not require a receiving or downstream repo to classify another project’s internal uncertainty unless the original issue explicitly made that part of acceptance.
When an external issue is created, it should be written for that receiving repo, not for a reader who already knows the submitting repo’s context.

## Evidence Expectations

Review tests as evidence, not as decoration.
Valid tests exercise the real production path or semantic requirement.
Be skeptical of helper-only tests, tautologies, assertions of the implementation’s own output, bypasses around the runtime/plugin/stub path, example-only coverage where the issue required full coverage, weakened assertions, and missing invalid-nearby cases where the fix could overgeneralize.

For plugin work, the evidence should usually distinguish valid generic behavior from invalid nearby ordinary Python and should not hard-code a downstream consumer.
For stubs work, the evidence should be source-backed: the upstream surface exists, the stub matches public behavior, no fake API is added, no Any/object opacity escape is introduced, and inherited-method inflation is not used unless source exposes that surface.

Watch for code-level laundering: hard-coded consumer names, support for local research abstractions as if they were external API, fake stubs, broad Any/object escapes, line suppressions, diagnostic filtering, deletion of required data, broad type widening, and any move that makes checks pass by weakening the problem instead of solving it.

## When Acting on Review Feedback

A positive disposition requires a commit.

Do not resolve an accepted review comment until the code/proof remediation is committed and the reply cites the commit.

Never reply “accepted,” “aligned,” “fixed,” “addressed,” or “will address” to a review thread unless the remediation is already committed.
A thread cannot be resolved on intent or future work.

Every substantive review item must receive its visible thread- or surface-local disposition and evidence before resolution.
The canonical field contract and state machine live in [[pr-feedback-triage/SKILL|pr-feedback-triage]]. Do not create top-level disposition ledgers or tracked review-log files.
Migrate legacy ledger-only resolutions by posting the canonical disposition and evidence on each affected thread before treating it as closed.

Review comments are not implementation specs.
The worker must translate accepted feedback into first-principles remediation requirements before assigning implementation.

For each comment:
- Identify the concern.
- Identify the proposed fix.
- Decide whether the concern is true under global + repo policy.
- Decide whether the proposed fix preserves those policies.
- If the concern is true but the fix is wrong, apply a policy-compatible remediation.

## Writing the Review

Write nuanced feedback for an intelligent reader.
Do not force a machine-readable template, a mandatory table, or a simplistic pass/fail label when prose communicates the situation better.
Do make the completion judgment clear: whether the original task can be considered complete, what evidence supports that judgment, and which unresolved requirements block completion if any remain.

Do not foreground effort, progress, good intentions, volume of work, or “substantial” partial implementation when required work remains.
Mention completed pieces only when they are necessary to identify the exact remaining blockers or to prevent redoing already-correct work.
Do not compare incomplete work to “no work done” or “completely fake work”; compare it to the expected standard: the task done correctly, completely, and provably.

When required work remains, lead with the incompleteness and the concrete blockers.
Do not make the reader excavate the missing work from beneath praise, context-setting, or a narrative of what did get done.

Nuance belongs in the evidence and blocker analysis, not in softening the completion standard.
The review should make it easy to finish the work, not easy to feel satisfied with less than the original contract required.

<!-- agent-memory:start -->
# Agent memory

This repository uses the central agent memory vault at `/home/dzack/.agent-memory-vault`.

Project memory key: `projects/github.com__dzackgarza__sage-indefinite-port/index`.

Repository `.agents` and `.hermes` paths are symlinks to the same vault-owned project directory.

Before changing architecture, search both project and global memory:

```bash
agent-memory search --scope both "<task or subsystem>"
```

Record durable repo-specific lessons with:

```bash
agent-memory add --scope project --type decision --title <title> --content <content>
agent-memory add --scope project --type trap --title <title> --content <content>
agent-memory add --scope project --type advice --title <title> --content <content>
agent-memory add --scope project --type context --title <title> --content <content>
agent-memory add --scope project --type reference --title <title> --content <content>
```

Plan work is card-backed. Create and update plan cards with `agent-memory plan add` and `agent-memory plan update`, not `agent-memory add --type plan`.

Use `agent-memory retrieve <key>`, `agent-memory update <key>`, and `agent-memory delete <key>` for memory CRUD.

The vault should be committed at all times. Treat staged or unstaged vault changes as an ephemeral error state. Before normal memory work resumes, load the bundled vault-maintenance skill with `agent-memory maintain skill vault-maintenance` and follow its referenced check, repair, and commit workflows.

Move reusable lessons during maintenance with:

```bash
agent-memory maintain move <key> --to global/advice
```
<!-- agent-memory:end -->

# Taking Work

The GitHub issue tree rooted at #4 is the only plan: its work units are the leaf issues, their blocked-by links are the dependency edges, and each unit's body holds its obligation and acceptance. Issue bodies cite "map §N" in [references/source-to-sage-translation-map.md](references/source-to-sage-translation-map.md).
Take the next unit with `uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port`, work it to the acceptance in its issue body, commit, close it with `itree close`, and take the next unit in the same turn.
A turn that ends with a ready unit untouched stops the repository until someone notices.

# Execution: complete the assigned port obligation

The issue's complete source translation and frozen consumers are the assignment.
A continuation, local failure, profile, or passing specimen does not select a new
assignment. TODO.md holds the remaining work; .tmp/HANDOFF.md holds the current
continuation point. Update those existing documents only when the next action
changes. Do not create another plan, receipt, audit, or progress ledger.

## Choose the next action from the actual state

| Observed state | Required next action |
| --- | --- |
| A run is still executing | Reuse its process/log; do not launch a duplicate. Continue independent implementation without editing its source snapshot. |
| An exact consumer failed with a wrong result or exception | Locate the first divergence from the pinned routine and its callees. Repair the complete source step, including its other branches and callers; run that consumer directly. |
| A diagnostic timer interrupted a case | Record it as interrupted. Do not infer a defect or start optimization. Run it detached without that diagnostic cutoff, or retain it as owed while running independent families serially. |
| The exact repaired consumer passes | Run the affected family; do not keep polishing that specimen. |
| A source step exists but has no consumer evidence | Verify its source correspondence and run its frozen consumers. Do not implement it again. |
| A source step or specialized helper is absent | Translate that step and its called helpers from the pinned source. Do not wait for another failing example to reveal each omitted branch. |
| A provider operation is absent | Keep its dependent consumers unresolved, record the exact provider contract at its owner, and continue independent port implementation against that API. |
| The assigned unit's full implementation and consumers are satisfied | Release its justified readiness markers and finish its issue acceptance; take the next existing ready unit. |

Within a selected algorithm unit, use the existing translation map to identify
its complete source routine/helper boundary and frozen consumer files. Account for
all source branches, not only the path taken by a failing specimen. Inspect current
port bodies before listing missing work. Preserve that correspondence as source
comments or a concise update to the existing map when it changes. A helper pass
never replaces the unit's consumer obligation. Do not turn source accounting into
a prerequisite audit of the entire repository.

## Keep validation advancing

A repair's first rerun is the exact failed node. Run earlier passing nodes again
only when changed code or dependencies can affect them, or during required final
acceptance. After the targeted check, run its affected family. Complete inventory
runs omit fail-fast and collect the full outcomes. Separately run independent
families when one remains interrupted; keep every omitted case explicitly owed.
These development runs do not waive the full final-source corpus.

Do not add a wall-clock budget to a mathematical requirement. Use an explicit
owner/issue budget where one exists, naming its source. The runner's 1,800-second
cutoff establishes no bound on perfect-domain enumeration. A profiling session
must answer a concrete source discrepancy or required execution problem; its
output alone cannot assign the next optimization. A slow successful computation
returns to the next missing capability or consumer check. Implement source-faithful
improvements needed for actual failures; do not invent faster algorithms.

## Continue without restarting

At interruption, replace the handoff's next action with the exact outstanding
source step or consumer node, the existing source correspondence, and the result
or running-job pointer. State separately what is implemented, what failed, and
what has not run. On continuation, resume that action. Check only whether its
inputs changed; do not repeat repository discovery, ownership analysis, source
reading, or successful probes whose relevant inputs are unchanged. A changed
input invalidates its affected evidence, not the whole project's accumulated work.

Keep one Sage mathematical process at a time and preserve its exact source and
dependency snapshot. While it runs, work on independent source in an isolated tree
if necessary. Reuse coherent implementation checkpoints; remove readiness markers
only after their complete obligations pass. No source-presence claim, benchmark,
partial result or administrative closure constitutes full acceptance.

# Nonnegotiable: port the existing algorithm; do not design a new one

**This project is a port, not an algorithm-development project.** The specialized algorithms in scope have already been implemented in the pinned open-source reference. The task is to translate **that implementation** into native Sage/Python, replacing its low-level and generic callees with the assigned Sage backends or research preamble operations. Do not search for alternative algorithms, reinvent an algorithm, or describe a known upstream algorithm as an open problem. The no-C++-build rule below prohibits compiling the reference; it does **not** authorize replacing its mathematics.

Before authoring or revising any specialized routine:

1. **Read the pinned source first**, including the routine's called helpers, not just its name or a summary. Use `references/vendor/Indefinite.jl@374a5ebb` and `references/vendor/polyhedral_common@1592b246`, and the vendored Allcock algorithm paper at `references/vendor/arxiv/2110.03784/` for edgewalk, as indexed by `references/source-to-sage-translation-map.md`. Identify the exact source symbol, file and line span for the routine and its relevant dependencies. The translation map is a navigation aid; the source is the algorithm specification. If a referenced implementation is absent from a sparse vendored snapshot, obtain that reference instead of improvising an algorithm.
2. **Translate, don't substitute.** Preserve the source algorithm's branching, mathematical invariants, enumeration scope, termination conditions, certificate checks, and treatment of boundary/exceptional cases. Replace a generic called operation (e.g. lattice algebra, exact CVP, polyhedral faces, finite-action calculations) with the corresponding existing backend/preamble primitive; do not replace the specialized caller with a newly devised search, shortcut, oracle, heuristic, or enumeration bound. In particular, a basis-dependent fast path or a special-case formula that avoids repairing the source algorithm's general path is **not** completion.
3. **Account for each deviation.** When a direct line-by-line translation is impossible because of representations or dependency contracts, identify the specific source step, the assigned callee in Sage/preamble, and why the translation has the *same inputs, outputs and mathematical guarantees*. A proposed mathematical change requires a source-backed proof or a concrete counterexample to the reference, recorded on the owning issue **before** adopting the change. Do not silently change the algorithm to get a fixture passing. If the source itself is incomplete under the issue's full acceptance, that is unresolved implementation, not permission to invent a replacement.
4. **Validate against the original full obligation**, using the frozen consumer paths and every applicable upstream fixture (not an xfail, a sample, a local helper-only check, or counts alone). A discrepancy requires locating the precise divergence from the pinned source and fixing it at that step. Passing on a subset never closes a unit requiring all cases.

**Current source anchors for #15:** `Lorentzian.g:315–429` (`LORENTZ_FindPositiveVectors`, `LORENTZ_SearchInitialVector`), `:705–968` (`GetUpperBound`, `LORENTZ_Kernel_Flipping`), `:971–1232` (`LORENTZ_GetShortPositiveVector`, `LORENTZ_GetOnePerfect`, `LORENTZ_DoFlipping`), and `:1239–1371` (complete perfect-cell quotient traversal), plus `IndefiniteFormsFundamental.g:11–30` (the perturbation helper). Allcock edgewalk and binary-cycle source correspondences belong to the same rule. These anchors must be read before editing the corresponding port bodies. A prose comment citing a source line is **not** evidence of equivalence when the code does something else.

# Working Against the Frozen API

- **The API is frozen, and work is written against it.** The boundary is frozen in #33 and [references/port-design.md](references/port-design.md). The design separates public mathematical owners from computational owners. A preamble method may dispatch to this port. Planned primitive names may be called before implementation, but a missing specialized algorithm remains port work.
- **Preamble deficiencies are never blockers.** A missing, wrong or restricted preamble operation is filed (or found) as a research issue that names the planned API, and port work continues against that name. A call into an unimplemented operation fails at run time. That decides only whether acceptance cases pass, never whether port work proceeds. No port issue is blocked by a research issue; blocked-by edges connect port units only.
- **One active worker owns each working tree.** While `research` has its own active managed workstream, this port's worker does not edit, stage or commit in `/home/dzack/research`, even for a narrow generic preamble fix. Record the exact source-derived defect at its research owner and continue the independent port implementation against the frozen API, retaining the failed acceptance case until its dependency exists. If the current port task truly cannot proceed without that other repository's implementation, the fleet's serial-blocker rule governs a temporary workstream switch; two workers must not compete for the research index or modify one another's uncommitted source. An explicitly assigned serial upstream-repair interval has one worker in that repository and ends before the port worker resumes.
- **Separate operation ownership from implementation routes.** Research accepts a computation only when its implementation needs no knowledge of the specialized upstream algorithms being ported here and it is a reasonably primitive mathematical construction usable independently in other problems. State its given data, hypotheses, output and independent route. Public ownership of a kernel, fibre, section or centralizer does not assign every algorithm answering that operation to research. Generic image/word-section computations consume supplied generators; they do not discover arithmetic parabolic generators or their liftable image. Research owns independent equation solvers, partial-isometry extension and pointwise-perpendicular kernels. The port owns specialized searches, compatibility-action construction, image discovery and recursive assembly. Neither source residency nor indefiniteness alone assigns ownership. See [port design](references/port-design.md#assigning-a-computation).
- **Use the concrete assignments.** [Implementation assignments](references/port-design.md#implementation-assignments) names the actual functions and repair owners on both sides. In particular, implement `IsotropicInduction.line_automorphism_generators`, `plane_automorphism_generators`, `flag_automorphism_generators` and the isotropic branch of `EichlerRecursion.pointed_automorphism_generators` here. Calling research's unresolved parabolic `.group_generators()` does not implement them. Consume research's independent `R.rational_lifts(...).integral_members()` and `R.pointwise_perpendicular_kernel()` computations. Move the copied binary cycle and scaled-shell multiplier search here, and implement Eichler approximate/cover/full-group assembly here; frozen preamble names may dispatch to those implementations. Research's matching instructions are in `/home/dzack/research/AGENTS.md`, “Concrete boundary with sage-indefinite-port”.
- **No new tests.** The acceptance suite (tests/acceptance, built in #33) is the only suite and is frozen. No unit adds, edits or weakens a case. A capability lands by removing the xfail markers that cite its unit. Changing a case, or changing the frozen API, means reopening #33 with the reason.
- **State questions as mathematics.** Orbit, stabilizer and equivalence questions are Isom and Aut of marked lattices: diagrams of lattice morphisms ending at L, with natural isomorphisms as morphisms (references/port-design.md, "Marked lattices"). A pointed lattice ⟨n⟩ → L lives in the coslice ⟨n⟩/Lattices (the marking is fixed pointwise); sublattices and configurations live in the arrow category Ar(Lattices) (the marking may move by an automorphism of its domain); flags live in Fun([2], Lattices). Do not introduce ad-hoc vocabulary ("transporter", "stabilizer data", "full generating data") for constructions that have standard names.

# Commit Gates and Work Claims

- **Red gate protocol**: The first time a commit gate, hook, or QC stage goes red, diagnosing that failure becomes the current task. Stop authoring; root-cause and fix the gate, or report it as a blocker with a reproducer. Never keep writing code behind a red gate, and never accumulate uncommitted work around one. A gate that is red on two consecutive commit attempts is a defect to diagnose, not an environment condition to wait out.

- **Claim freshness — no off-ledger work**: At every claim and every release, reconcile the shared queue/claim state against actual repository state across all branches before selecting work. Never select work from a queue older than your last branch sync. All authoring requires a live claim; batch-committing a body of work authored off-ledger is prohibited.

- **The `[unverified]` commit tag**: `[unverified]` may appear in a commit message only while the repository's verification phase is formally deferred by a named contract, and every `[unverified]` commit must name the contract or terminal phase that discharges it.

- **Deferred type checking (temporary, discharged by #36)**: while #36 is open, the push gate's mypy errors are acknowledged debt and are not paid down. A push may then use `git push --no-verify`, and only under all of these conditions:
  - `just test-push` was run on the commit being pushed and failed in its mypy stage (`_sage-mypy`) and nowhere earlier.
  - The tip commit of the push carries this trailer, verbatim apart from the count: `Type-Debt-Ack: MYPY-DEFER-7Q4K #36 errors=<N>`, where `<N>` is the "Found N errors" count of that run.
  - `<N>` is no greater than the count in the previous `Type-Debt-Ack` trailer in `git log`. A push that would raise the count does not use this exception: fix the new errors first.

  The trailer covers every commit in that push. Commits need no exception, because the commit gate does not run mypy. The exception covers only mypy debt: any other gate failure is fixed, never bypassed. #36's last commit deletes this rule.

# Architecture and Dependencies

- **Sage runtime on this host** (verified 2026-10-17): Docker is installed by
  Snap at `/snap/bin/docker`, not `/usr/bin/docker`. The `sage-dev-develop`
  container is already running. Use
  `/usr/bin/sg docker -c '/snap/bin/docker ...'` to access it; an unsuccessful
  `command -v docker` does not mean Docker is unavailable.
  The wrapper `~/.local/bin/sage` calls `docker` from `PATH`, and the preamble
  Julia bridge requires `/home/dzack/.juliaup/bin/julia`. The tracked `.envrc`
  configures `PATH`, `SAGE_BIN`, `SAGE_PYTHON_BIN`, `SAGE_PREPARSE_BIN`, and
  `SAGE_JULIA_COMMAND` for the container-backed Sage and central QC recipes.
  Verified with `sage --version` (SageMath 11.0.beta1) and a preamble import.
  Keep one Sage mathematical test process at a time.

- **The `research` preamble is the lattice substrate and a co-developed dependency**:
  This repository consumes `dzack_research.preamble` (categories, formed modules, lattices, isometries, discriminant modules) from the checkout at `/home/dzack/research/src`. `.envrc` puts it on `PYTHONPATH` and `MYPYPATH`; `pyproject.toml` declares `dzack-research`.
  - **Use the preamble's own classes and types**: annotate with `Lattice` (and `Lattice.Element`), `LatticeIsometry`, `LatticeIsometryHomset`, and the tensor and presented-module classes. Never add Protocol, stub, or facade types over the preamble, and never edit `sys.path`.
  - **Fix preamble defects upstream**: a missing annotation, a class reachable only from a private module, or a wrong result in the preamble is fixed in `research` under its own QC. Port code is written against the corrected or planned API meanwhile and does not wait (see "Working Against the Frozen API"). No local workaround, no fallback to raw SageMath lattices (their indefinite support is insufficient), no reinvention of lattice, discriminant, or formed-category infrastructure here.
  - **Ownership**: `research` owns mathematical interfaces and independently reusable primitive computations. This repository owns approximate models, coverings, indefinite recursion and traversal, binary edgewalk cycles, scaled first-success searches and arithmetic parabolic assembly. Research may dispatch public methods to these port algorithms. Invariant decomposition and gluing from supplied component groups are research primitives; non-scalar centralizers need an identified computational route and cannot be assigned to research merely by their name. The complete centralizer consumer obligation remains in #18. The operation-by-operation split is in `references/port-design.md`.

- **NEVER BUILD THE C++ CODE**:
  **Do NOT attempt to compile, build, configure, or invoke C++ compilation toolchains** for `polyhedral_common` or any other C++ source in `references/`.
  The upstream C++ codebase in `references/polyhedral_common` is strictly a reference implementation for algorithm extraction, logic translation, and structural understanding.
  All production algorithms in this repository MUST be implemented as native SageMath / Python code delegating low-level operations to standard system backends (FLINT for exact integer linear algebra, Normaliz / cddlib / PPL for polyhedral cones). Research supplies independent primitive computations, including exact module algebra, specified equation fibres, partial-isometry extension, finite configurations and finite actions from supplied effective data. This port supplies specialized indefinite Isom, Aut, orbit, search and assembly algorithms. Reuse generic GAP/Bliss/backend operations through their owners. A claimed primitive gap must state its independent contract and input data; the output of this port cannot be assumed as its own prerequisite.

## Imported Contribution Policies (from `/home/dzack/research/CONTRIBUTING.md`)

### ARC: Mathematical Architecture & Ownership
- `ARC-01`: Own universal properties, categories, morphisms, functors, and adjunctions natively in this repository’s category framework.
- `ARC-02`: Represent subobjects as pairs `(S, iota: S -> M)`. Place predicates, isometries, embeddings, and containment checks on morphism data and hom-sets.
- `ARC-03`: Build structural objects/functors first, then derive numerical invariants.

### ENG: Computational Backend Delegation
- `ENG-01`: Delegate heavy computations to reliable exact backends (SageMath, Singular, OSCAR, Macaulay2, PARI/GP) when available.
- `ENG-02`: Do not hand-roll standard mathematics that mature upstream dependencies already provide.
- `ENG-03`: Keep owned logic minimal; offload heavy numerical work to engine backends through standard bridges.
- `ENG-04`: Prefer native engine implementations (e.g., Julia/OSCAR, Singular) for multi-step heavy computation when bridge overhead is worse.

### BRG: Interoperability and Bridge Boundaries
- `BRG-01`: Use structured, persistent bridge interfaces for external systems; avoid ad-hoc subprocess scripts with temp files for core algorithms.
- `BRG-02`: Validate and wrap data at bridge boundaries to keep backend representations out of public API types.

### ENV: Environment, Execution, and Tooling
- `ENV-01`: Use exact physical paths in shell/tool invocations (e.g., `/home/dzack/...`).
- `ENV-02`: Define project orchestration, gates, and doc generators in the root `justfile`.

### DEV: Development Discipline
- `DEV-01`: Use explicit typing on public APIs; avoid `Any`/`object` in public contracts.
- `DEV-02`: Add concrete falsifiable specimens for every new category, functor, or operation.
- `DEV-03`: Before adding code, run `just preamble-megadoc`, reuse existing constructions, and implement new work in the most general mathematical form before specialization.

## Imported Contribution Policies (from `/home/dzack/gitclones/sage-categories/CONTRIBUTING.md`)

### POL-SCOPE: Implementation Order and Scope
- `POL-SCOPE-001`: Build dependency chains in order and complete required milestones before starting later phases.
- `POL-SCOPE-002`: Execute only active-phase work with accepted prerequisites. Do not start unsupported descendants early.
- `POL-SCOPE-007`: Measure success by categorical ownership, auditable declarations, and explicit mathematical structure.
- `POL-SCOPE-008`: Keep the mathematical declaration as the source of truth. Keep runtime representation private except for authorized SymPy propositions.

### POL-SHADOW: Package-owned Mathematical Surface
- `POL-SHADOW-001`: Build a package-owned categorical replacement for the supported Sage surface.
- `POL-SHADOW-003`: Keep public API closed over package-owned categories, objects, morphisms, and propositions.
- `POL-SHADOW-005`: Do not guarantee compatibility with unsupported arbitrary Sage API behavior.
- `POL-SHADOW-006`: Absorb required Sage constructions into the package-owned architecture.

### POL-ONT: Foundational Ontology
- `POL-ONT-001`: Treat raw Python values as private carriers, not as category members.
- `POL-ONT-002`: Use separate constructor and refinement layers; avoid runtime-type-sniffing constructor overloads.
- `POL-ONT-008`: Do not use catch-all fallback constructors. Route each input mode through explicit constructors and fail hard.
- `POL-ONT-010`: Survey existing generic categorical machinery before adding greenfield implementations.
- `POL-ONT-011`: Do not add standalone procedural helper functions for core predicates.

### POL-MATH: Mathematical Architecture
- `POL-MATH-001`: Identify object, element, morphism, category, functor, and universal-property owners before implementation.
- `POL-MATH-002`: Model named categories as categories, not as utility-like classes.
- `POL-MATH-003`: Model functors with explicit object and morphism maps.
- `POL-MATH-008`: Do not duplicate data that a defining morphism already determines.
- `POL-MATH-014`: Use inspected mathematics to justify construction sites; do not use tests as proof of mathematics.
- `POL-MATH-034`: Every truth question gets one category-owned predicate and a SymPy proposition; `ask()` resolves it.
- `POL-MATH-037`: Construct values in the stated mathematical category; constructors do not certify proofs.

### POL-CAT: Category Ownership and Inheritance
- `POL-CAT-001`: A category owns constructors, local operations, and role types.
- `POL-CAT-002`: Use category-owned `ObjectType`, `ElementType`, and `MorphismType`.
- `POL-CAT-004`: A category level only declares its own structure and operations.
- `POL-CAT-006`: Do not re-export methods owned by another category.
- `POL-CAT-018`: Distinguish property subcategories from “data-carrying” categories.
- `POL-CAT-021`: Keep `Mor(n, C)` and endpoint categories as first-class objects in the categorical layer.

### POL-REP: Semantic Representation
- `POL-REP-001`: Treat Sage vectors/matrices as private representations; own semantic objects in the category layer.
- `POL-REP-002`: Return semantic mathematical objects at public API boundaries.
- `POL-REP-003`: Compare elements and morphisms semantically, not by raw coordinate unpacking.
- `POL-REP-009`: Lower private semantics to computation once and reconstruct semantic outputs before returning.
- Local repository rule: do not use raw matrices as mathematical API objects. Use morphisms for structure and tensor valences for numerics (for example, `(0,2)` Gram tensors and `(1,1)` endomorphism tensors).

### POL-ENGINE: Computation Engine Boundary
- `POL-ENGINE-001`: Public API is owned mathematics; engines provide private realizations and algorithms.
- `POL-ENGINE-002`: Keep engine types private; only authorized SymPy proposition expressions may cross the boundary.
- `POL-ENGINE-007`: Do not introduce selectable backends or replaceable engine abstractions in the public API.
- `POL-ENGINE-015`: Use fixed dependency assignments for private algorithms and reconstruct exact semantic outputs.

### POL-FORM: Forms and Lattices
- `POL-FORM-001`: Model lattices as `R`-modules with a specified form, not as raw free `ZZ` modules.
- `POL-FORM-002`: Encode bilinear forms via Gram tensors (`M ⊗ M -> W`).
- `POL-FORM-004`: Do not assume positivity, freeness, embeddedness, or unimodularity by default.
- `POL-FORM-006`: Define complements, norms, and reflections under correct hypotheses only.
- `POL-FORM-008`: Use exact arithmetic and exact coefficient rings.

# Existing Preamble Lattice API (Obviates Foundational Tasks)

The user's `research` repository (`dzack_research.preamble.categories.lattices`, `lattice_morphisms`, `modules`, and `forms`) owns and provides the formed lattice algebra and categorical infrastructure below. This surface inventory does not certify every listed computational method as implemented or assign its specialized algorithm to research.
**Do NOT reinvent or re-implement any of the following infrastructure; consume it directly from `dzack_research.preamble`. This completely obviates greenfield implementation for Phase 2 and early foundational tasks.**

## Full lattice interface in `preamble-megadoc.md` (`#subsystem-lattices`)

### Construction and category entry points
- `Lattices` as the lattice category over a ring with constructor:
  - `Lattices(cls, *args)`
  - `Lattices(...)(self, data, basis=None, names=None, form=None, module_generators=None)`
- Category constructors:
  - `FiniteRankLattices`, `NondegenerateLattices`, `RationalLattices`, `EvenLattices`, `RootLattices`
- `Lattices(ZZ)` named lattice constructors:
  - `Lattices(ZZ)("U")`, `"U(m)"`, `"An"`, `"Dn"`, `"En"`, `"E8"`, `"II_{1,9}"`, `"Leech"`, and Gram-tensor constructors.
- Direct sum and scaling:
  - `L1 + L2` and `L.twist(scale)`.

### `Lattices` parent object API (full interface)
- Isometry and homset constructors:
  - `Aut`, `O`, `orthogonal_group`
  - `SO`, `special_orthogonal_group`
  - `Emb(codomain)`, `Isom(codomain)`
  - `hom(images, codomain=None)`, `identity_morphism`
  - `similarity_homset(other, scale)`, `similarity(scale, images=None, codomain=None)`
- Pairings and form evaluations:
  - `b(left, right)`, `q(vector)`, `correlation`, `correlation_morphism`
  - `gram_tensor`, `gram_matrix(basis=None)`, `discriminant`, `rank`, `signature_pair`, `level`, `genus`
- Lattice reductions and basis control:
  - `LLL()`, `BKZ(block_size=20)`, `HKZ()`
  - `lll_reduction`, `bkz_reduction`, `hkz_reduction`, `orthogonal_group`, `stable_orthogonal_group`
- Duality and discriminant data:
  - `dual_module`, `dual_lattice`, `metric_dual`, `discriminant_module`, `discriminant_projection`
  - `discriminant_group`, `discriminant_bilinear_form`, `discriminant_quadratic_form`
  - `discriminant_class`, `divided_discriminant_class`, `discriminant_representation`
  - `discriminant_image`, `discriminant_representation_is_surjective`
- Orthogonality and decomposition:
  - `radical`, `radical_quotient`, `orthogonal_complement`
  - `module_generating_set`, `module_generator`, `biproduct_factors`, `decomposition`, `decomposition_names`, `summands`
  - `is_decomposable`, `indecomposable_name`, `primitive_isotropic_subobject`
  - `local_modification(prime, *discriminant_classes)`, `overlattice(*discriminant_classes)`, `even_overlattice_inclusions`
- Arithmetical predicates:
  - `is_finite_rank`, `is_even`, `is_nondegenerate`, `is_unimodular`
  - `is_definite`, `is_positive_definite`, `is_negative_definite`, `is_p_elementary(prime)`
  - `is_locally_isometric(other, prime)`, `is_isometric(other)`, `is_similar(other, scale)`
- Discrete geometry and orbits:
  - `minimum`, `successive_minima`, `shortest_vectors`, `roots`, `roots_of_square(square)`
  - `vectors_of_square(square)`, `vectors_of_square_and_divisibility(square, divisibility)`
  - `kissing_number`, `root_sublattice`, `theta_series(precision=20, variable='q')`
- Isotropic and orbit APIs:
  - `isotropic_flag`, `isotropic_flag_orbit_representatives(rank=2)`
  - `isotropic_line_orbit_representatives`, `isotropic_plane_orbit_representatives`
  - `voronoi_cell(bound=None)`, `voronoi_relevant_vectors`
- Group and positivity APIs:
  - `stable_orthogonal_group`, `spinor_kernel_subgroup`, `positive_cone_subgroup`
  - `twist(scale)`, `reflection(root)`
  - `two_elementary_invariants`, `delta`
- Packing and metric helpers:
  - `babai(target)`, `closest_vector(target)`, `contact_polytope`, `covering_radius`
  - `packing_density`, `packing_radius`, `center_density`, `hadamard_ratio`, `hermite_invariant`

### `Lattices` element methods
- `b(other)`, `norm`, `q()`, `div()`, `divisibility_ideal`
- `divided_discriminant_class()`, `is_root()`
- `to_list`, `to_tuple`, `to_vector`, `monomial_coefficients`

### Subcategory APIs
- `FiniteRankLattices`: `is_finite_rank`
- `RationalLattices`: `fraction_field`, `is_nondegenerate`
- `EvenLattices`: `is_even` inherited by lattice parent
- `RootLattices`:
  - parent methods: `cartan_type`, `coxeter_number`, `highest_root`, `simple_reflections`, `simple_roots`, `fundamental_weights`
  - element methods: `coroot`, `height`, `is_negative_root`, `is_positive_root`

### Morphisms and homsets
- `LatticeMorphism` (`ModuleMorphism`):
  - constructor: `__init__(self, parent, images)`
- `LatticeEmbedding` (`LatticeMorphism`):
  - `discriminant_inclusion`, `is_injective`
- `LatticeIsometry` (`LatticeEmbedding`):
  - `centralizer_discriminant_image`, `cyclic_subgroup`, `determinant`, `discriminant_isometry`
  - `discriminant_morphism`, `formed_coinvariants`, `invariant_lattice`, `inverse`
  - `is_surjective`, `preserves_positive_cone`, `real_spinor_norm_sign`
- `LatticeHomset`:
  - `__init__(domain, codomain)`, `_element_constructor_(images)`
- `LatticeEmbeddingHomset`:
  - `an_element`, `even_overlattice_inclusions`, `is_empty`
- `LatticeIsometryHomset`:
  - `act`, `acting_group`, `compose`, `discriminant_image`, `discriminant_preimage`
  - `group_generators`, `identity`, `is_empty`, `isotropic_equivalence_witness`
  - `isotropic_orbit_representatives`, `isotropic_stabilizer_generators`, `number_of_group_generators`
  - `one`, `order`, `transporter`, `vector_equivalence_witness`, `vector_orbit_representatives`
  - `vector_stabilizer_generators`, `vectors_are_equivalent`

### Lattice objects and named tables
- `CoxeterDiagram`
  - `cardinality`, `connected_components`, `elliptic_subdiagrams`, `graph`, `coxeter_entry`
  - `coxeter_matrix`, `is_connected`, `is_elliptic`, `is_hyperbolic`, `is_parabolic`
  - `is_rooted`, `index_set`, `induced_subdiagram`, `preferred_positions`, `roots`
  - `root_gram_tensor`, `schlafli_tensor`, `signature_pair`, `vertex_names`
- `Genus`
  - `class_number`, `determinant`, `discriminant_form`, `excess`, `exists`
  - `local_symbol`, `level`, `mass`, `representative`, `representatives`, `signature_pair`
- `IsotropicFlag`
  - `basis`, `lattice`, `rank`, `terms`, `top`
- `OrthogonalCharacterQuotient`
  - `image`, `image_keys`, `splitting_isometries`, `stabilizer_image_keys`
  - `subgroup_image_keys`, `witness_meets_subgroup`
- `VectorPrimitiveExtension`
  - `class_of_representative`, `complement_is_definite`, `representative_of`
- registries: `register_indecomposable`, `register_indecomposable_gram`

### Helper functions from the lattice subsystem
- lattice constructors and homset factories:
  - `lattice(...)`, `lattice_embedding_homset`, `lattice_homset`, `lattice_isometry_homset`, `lattice_latex`
- lattice operators:
  - `diagonal_gram`, `orthogonal_sum`, `scale_gram_tensor`, `signature_pair_of_gram`
  - `colimit_lattice`, `discriminant_of_gram`, `signature_pair_of_gram`
- orbits, stabilizers, and equivalence:
  - `definite_complement_extensions`, `isotropic_equivalence_witness`, `isotropic_orbit_representatives`, `isotropic_stabilizer_generators`
  - `subgroup_isotropic_are_equivalent`, `subgroup_isotropic_orbit_representatives`
  - `subgroup_vector_orbit_representatives`, `subgroup_vectors_are_equivalent`
  - `transport_isotropic_object`, `vector_equivalence_witness`, `vectors_are_equivalent` (homset-level)
- reduction and algorithm helpers:
  - `babai`, `bkz_reduction`, `hkz_reduction`, `closest_vector`, `gaussian_heuristic`
  - `generator_pairings`, `gluing_route_discriminant_classes`
  - `packing_radius`, `packing_density`, `center_density`, `kissing_number`
  - `roots`, `roots_of_square`, `vectors_of_square`, `vectors_of_square_and_divisibility`
  - `voronoi_cell`, `voronoi_relevant_vectors`
  - `oscar_centralizer_discriminant_image`, `oscar_even_unimodular_primitive_embedding`, `oscar_rational_spinor_norm_sign`

### Compatibility note
- The above interface is the preamble contract. This repo must use it and should not introduce local fallback lattice mechanics in place of these names.

---

### What Actually Needs to be Built in `sage-indefinite-port`
Consume the existing foundational formed algebra. Development in this repository focuses on:

1. **Lorentzian perfect-domain reduction**: Traversal of the $(1,n)$ reduction complex for $O^\Omega(L)$ and $O(L)$ generators.
2. **$2U$-Eichler approximate models & Higher-Witt-index recursion**: Generating $A(L)$ and lifting covering sets to full $O(L)$ for Witt index $\geq 2$.
3. **Parabolic recursion & Isotropic orbit decomposition**: Compose research's group actions, sections and kernels with recursive reduced-group computations to obtain cusp and flag stabilizers and the frozen rank-at-most-two isotropic orbits.
4. **Centralizer and arithmetic-subgroup consumers**: Supply indefinite orthogonal-group and orbit computations to research's general centralizer, gluing and subgroup operations; verify all frozen consumer cases.
