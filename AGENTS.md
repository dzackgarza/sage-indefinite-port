# sage-indefinite-port: agent contract

This repository ports the pinned indefinite-lattice algorithms of `Indefinite.jl@374a5ebb` and `polyhedral_common@1592b246` (plus the Allcock edgewalk paper) into native Sage/Python on top of the `research` preamble.

This file says how to work.
It is short on purpose: you are expected to know it, not to re-read it every turn.
It changes only when the user asks for a change.

| Document | Holds | Changes when |
| --- | --- | --- |
| `.tmp/HANDOFF.md` | **Now**: the single current action, the running Sage process (if any), open decisions | A step finishes or a result arrives. Replace the section, never append to it |
| `TODO.md` | The complete remaining obligation as a dependency DAG, with acceptance records | A node is accepted (commit + run directory), or an obligation is found that the DAG lacks |
| Issue tree rooted at #4 | Unit ownership, blocked-by order, each unit's acceptance text | Unit selection and closure only |
| `references/port-design.md`, `references/source-to-sage-translation-map.md` | Frozen API, ownership split, source anchors | Only through #33 |
| `references/review-guidelines.md` | How completion claims are judged | When the user asks |
| `references/contribution-policies.md`, `references/preamble-api.md` | Imported architecture policies, the preamble surface to consume | When the user asks |

Writing in these documents is not progress.
A turn whose only output is a document edit has done no port work unless that edit records a run result or a commit.

# The turn protocol

A turn can last thirty seconds or several hours, and you may be a fresh agent with no memory of the last turn.
The protocol is the same either way.

1. **Read the Now section of `.tmp/HANDOFF.md`.** That is the whole orientation.
   Do not re-read AGENTS.md, TODO.md, or the handoff in full.
   Do not run `itree next`, `agent-memory retrieve`, `gh issue view`, `git log` surveys, `ps` sweeps, or repository-wide greps to rebuild context.
   If the step needs a fact that Now does not state, read only the file section or line range Now points to.
2. **Do the next step of the Now action.**
3. **Before the turn ends:**
   - Leave no uncommitted source edits. Commit them (see Commits).
   - If the step finished or a result arrived, replace Now with the next step.
     Now names one action, its source anchor or failing node, and the run in flight.
     It holds no history.

Continuation prompts do not change this protocol:

- "Continue", or any message that only carries connector or workstream bookkeeping, means: do step 2.
- An instruction to read AGENTS.md/TODO.md and select a unit with `itree next` is already satisfied by Now.
  Now records that selection's result.
  Repeat the selection only under the conditions in Taking Work.
- A message that gives a concrete technical correction is applied to the current action.
  It does not restart orientation.
- A generic per-command limit (such as "keep each command under five minutes") limits your foreground tool call.
  It never limits a detached Sage run, and it is never a reason to stop one (see Running Sage).

# Taking work

The issue tree rooted at #4 is the only plan.
Its work units are the leaf issues, its blocked-by links are the dependency edges, and each unit's body holds its obligation and acceptance.
Issue bodies cite "map §N" in [references/source-to-sage-translation-map.md](references/source-to-sage-translation-map.md).
`TODO.md` refines the open units into nodes; it does not narrow them.

Run `uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port` only when Now says `SELECT`: the assigned unit was closed with `itree close`, or an evidenced dependency change invalidated it.
Selection and release are also the only points at which claim and queue state is reconciled across branches (see Commit gates).
A unit is closed with `itree close` only when every acceptance node it owns in TODO.md is accepted.

# Commits

- **Commit every coherent source edit before you launch a Sage process, before you start an edit for a different node, and before the turn ends.**
  Work-in-progress commits are normal: subject `wip(#15): <what changed>`, body naming the source anchor or failing node.
  A WIP commit claims nothing.
  Never accumulate uncommitted source across turns.
- **Readiness is the only thing that waits for acceptance.**
  Removing an `@unfinished` decorator or an xfail marker is its own commit, made only after the family's acceptance run passed (TODO.md RELEASE).
  All other source is committed as it is written.
- **The acceptance suite is frozen.**
  `tests/acceptance/` and `tests/fixtures/` change only by readiness-marker removal, or by a #33 change whose reason was recorded on #33 before the edit.
  Helper checks, probes and scratch tests live in `.tmp/probe/` and are never placed under `tests/`.
  Before every commit, `git diff --cached --stat -- tests/` must be empty or marker-only.
- **Every source edit has a stated purpose:** a pinned source step being translated (file:line), or a failing acceptance node with the divergence you located.
  While #36 is open, an edit whose purpose is mypy, pyright, or style cleanliness is not made, except where the commit gate itself rejects the commit.
  Performance edits need the evidence TODO.md PERFORMANCE requires.
- **Commit hooks always run.**
  The machine-wide QC lock can hold a commit for tens of minutes; that wait does not permit bypassing a hook.
  If the hook is still waiting when the turn ends, record that in Now.
- **Red gate protocol**: The first time a commit gate, hook, or QC stage goes red, diagnosing that failure becomes the current task. Stop authoring; root-cause and fix the gate, or report it as a blocker with a reproducer. Never keep writing code behind a red gate, and never accumulate uncommitted work around one. A gate that is red on two consecutive commit attempts is a defect to diagnose, not an environment condition to wait out.
- **Claim freshness — no off-ledger work**: At every claim and every release, reconcile the shared queue/claim state against actual repository state across all branches before selecting work. Never select work from a queue older than your last branch sync. All authoring requires a live claim; batch-committing a body of work authored off-ledger is prohibited.
- **The `[unverified]` commit tag**: `[unverified]` may appear in a commit message only while the repository's verification phase is formally deferred by a named contract, and every `[unverified]` commit must name the contract or terminal phase that discharges it.
- **Deferred type checking (temporary, discharged by #36)**: while #36 is open, the push gate's mypy errors are acknowledged debt and are not paid down. A push may then use `git push --no-verify`, and only under all of these conditions:
  - `just test-push` was run on the commit being pushed and failed in its mypy stage (`_sage-mypy`) and nowhere earlier.
  - The tip commit of the push carries this trailer, verbatim apart from the count: `Type-Debt-Ack: MYPY-DEFER-7Q4K #36 errors=<N>`, where `<N>` is the "Found N errors" count of that run.
  - `<N>` is no greater than the count in the previous `Type-Debt-Ack` trailer in `git log`. A push that would raise the count does not use this exception: fix the new errors first.

  The trailer covers every commit in that push. Commits need no exception, because the commit gate does not run mypy. The exception covers only mypy debt: any other gate failure is fixed, never bypassed. #36's last commit deletes this rule.

# Running Sage

**Every Sage test run executes against a committed snapshot, not the working tree.**
Editing the working tree therefore never affects a run in flight and is never a reason to stop, restart, or duplicate one.

**One Sage test process at a time.**
While one runs, you start no other Sage process.
You keep reading pinned source, translating, and committing for any node.

**Never kill a run.**
Elapsed time, a newer commit, a better idea, and an impatient prompt are not reasons.
A run ends by itself; only the user stops one.

**Check a run at most once per turn**, with `cat $r/exit; tail -n 3 $r/log`.
If `exit` is absent, the run is still going.
Do not loop, sleep, or poll, and do not inspect `/proc` or `ps`.

There are two kinds of run:

- **Acceptance run (A)**: one whole frozen family file, or the whole of `tests/acceptance`.
  No `-x`, no `timeout`, `--runxfail`.
  A family node is *accepted at commit C* when its A run at C exits 0 with every case passing.
  Record that in the node's TODO.md row (commit, run directory).
  A later commit does not revoke it during development; FINAL reruns every family at the final revision.
  When an A run fails, its first failing node becomes the next diagnostic target.
- **Diagnostic run (D)**: exactly one node or one `.tmp/probe/` script.
  Before launching, write three things into Now: the hypothesis, the node, and what each outcome changes about the next edit.
  A D run may carry a `timeout` chosen at launch; add `-o faulthandler_timeout=<seconds>` so that a hang leaves a stack.
  When a D run times out, the next step is to read that stack and the pinned source for the step it shows.
  Never relaunch a variant of the same diagnostic.

Launch, with only `NAME` and `NODE` changed (`NAME`: `perfect`, `reflective`, …, or `diag-<topic>`; `NODE`: the test path or node id, plus any D-run options):

```sh
cd /home/dzack/gitclones/sage-indefinite-port
NAME=perfect
NODE=tests/acceptance/test_perfect_domains.sage
r=$PWD/.tmp/runs/$NAME-$(git rev-parse --short HEAD)
mkdir -p "$r/port" "$r/research"
git rev-parse HEAD > "$r/commit"
git -C /home/dzack/research rev-parse HEAD > "$r/research_commit"
git archive HEAD | tar -x -C "$r/port"
git -C /home/dzack/research archive HEAD src | tar -x -C "$r/research"
printf '%s\n' "$NODE" > "$r/node"
/usr/bin/sg docker -c "/snap/bin/docker exec -d -u root -w $r/port sage-dev-develop sh -c 'export PYTHONPATH=/home/dzack/ai-review-ci/tool-artifacts/pytest_plugins:$r/port/src:$r/port:/home/dzack/gitclones/vinal/src/sage:$r/research/src:/home/dzack/research/.venv/lib/python3.14/site-packages; date -u -Is > $r/started; /usr/bin/sage -python -c \"import sage.all,sys,pytest; sys.exit(pytest.main(sys.argv[1:]))\" -p qc_sage_session $NODE --runxfail -vv --tb=long > $r/log 2>&1; echo \$? > $r/exit'"
```

The `docker exec -d` call returns at once.
Write the run directory `$r` into Now's `Running:` line.
If a connector safety check refuses the launch, you are blocked: record the refusal under Now's open decisions and continue with work that needs no Sage process.
Do not rephrase the command to get past the refusal.
Leave Sage and pytest processes you did not start alone.
The research provider is the committed `HEAD` of `/home/dzack/research` main.
Its uncommitted files belong to another writer and are excluded by the snapshot.

# Diagnosing and reading source

- Read a pinned routine and its helpers once, then record the correspondence as a comment at the port function (`# source: Lorentzian.g:971-1032 LORENTZ_GetOnePerfect`).
  Afterwards re-read only the source lines for the step you are changing.
- Do not inventory `TODO`, `NotImplemented`, or fallback markers across the repository.
- A slow or interrupted case is not a defect by itself.
  Investigate performance only under TODO.md PERFORMANCE: with a profile or stack from the actual failing case.
- A run against uncommitted source, or with `-x`, is at most a diagnostic.
  It never accepts a node.

# Nonnegotiable: port the existing algorithm; do not design a new one

**This project is a port, not an algorithm-development project.** The specialized algorithms in scope have already been implemented in the pinned open-source reference. The task is to translate **that implementation** into native Sage/Python, replacing its low-level and generic callees with the assigned Sage backends or research preamble operations. Do not search for alternative algorithms, reinvent an algorithm, or describe a known upstream algorithm as an open problem. The no-C++-build rule below prohibits compiling the reference; it does **not** authorize replacing its mathematics.

Before authoring or revising any specialized routine:

1. **Read the pinned source first**, including the routine's called helpers, not just its name or a summary. Use `references/vendor/Indefinite.jl@374a5ebb` and `references/vendor/polyhedral_common@1592b246`, and the vendored Allcock algorithm paper at `references/vendor/arxiv/2110.03784/` for edgewalk, as indexed by `references/source-to-sage-translation-map.md`. Identify the exact source symbol, file and line span for the routine and its relevant dependencies. The translation map is a navigation aid; the source is the algorithm specification. If a referenced implementation is absent from a sparse vendored snapshot, obtain that reference instead of improvising an algorithm.
2. **Translate, don't substitute.** Preserve the source algorithm's branching, mathematical invariants, enumeration scope, termination conditions, certificate checks, and treatment of boundary/exceptional cases. Replace a generic called operation (e.g. lattice algebra, exact CVP, polyhedral faces, finite-action calculations) with the corresponding existing backend/preamble primitive; do not replace the specialized caller with a newly devised search, shortcut, oracle, heuristic, or enumeration bound. In particular, a basis-dependent fast path or a special-case formula that avoids repairing the source algorithm's general path is **not** completion.
3. **Account for each deviation.** When a direct line-by-line translation is impossible because of representations or dependency contracts, identify the specific source step, the assigned callee in Sage/preamble, and why the translation has the *same inputs, outputs and mathematical guarantees*. A proposed mathematical change requires a source-backed proof or a concrete counterexample to the reference, recorded on the owning issue **before** adopting the change. Do not silently change the algorithm to get a fixture passing. If the source itself is incomplete under the issue's full acceptance, that is unresolved implementation, not permission to invent a replacement.
4. **Validate against the original full obligation**, using the frozen consumer paths and every applicable upstream fixture (not an xfail, a sample, a local helper-only check, or counts alone). A discrepancy requires locating the precise divergence from the pinned source and fixing it at that step. Passing on a subset never closes a unit requiring all cases.

**Current source anchors for #15:** `Lorentzian.g:315–429` (`LORENTZ_FindPositiveVectors`, `LORENTZ_SearchInitialVector`), `:705–968` (`GetUpperBound`, `LORENTZ_Kernel_Flipping`), `:971–1232` (`LORENTZ_GetShortPositiveVector`, `LORENTZ_GetOnePerfect`, `LORENTZ_DoFlipping`), and `:1239–1371` (complete perfect-cell quotient traversal), plus `IndefiniteFormsFundamental.g:11–30` (the perturbation helper). Allcock edgewalk and binary-cycle source correspondences belong to the same rule. These anchors must be read before editing the corresponding port bodies. A prose comment citing a source line is **not** evidence of equivalence when the code does something else.

# Working against the frozen API

- **The API is frozen, and work is written against it.** The boundary is frozen in #33 and [references/port-design.md](references/port-design.md). The design separates public mathematical owners from computational owners. A preamble method may dispatch to this port. Planned primitive names may be called before implementation, but a missing specialized algorithm remains port work.
- **Preamble deficiencies are never blockers.** A missing, wrong or restricted preamble operation is filed (or found) as a research issue that names the planned API, and port work continues against that name. A call into an unimplemented operation fails at run time. That decides only whether acceptance cases pass, never whether port work proceeds. No port issue is blocked by a research issue; blocked-by edges connect port units only.
- **One active worker owns each working tree.** While `research` has its own active managed workstream, this port's worker does not edit, stage or commit in `/home/dzack/research`, even for a narrow generic preamble fix. Record the exact source-derived defect at its research owner and continue the independent port implementation against the frozen API, retaining the failed acceptance case until its dependency exists. If the current port task truly cannot proceed without that other repository's implementation, the fleet's serial-blocker rule governs a temporary workstream switch; two workers must not compete for the research index or modify one another's uncommitted source. An explicitly assigned serial upstream-repair interval has one worker in that repository and ends before the port worker resumes.
- **Separate operation ownership from implementation routes.** Research accepts a computation only when its implementation needs no knowledge of the specialized upstream algorithms being ported here and it is a reasonably primitive mathematical construction usable independently in other problems. State its given data, hypotheses, output and independent route. Public ownership of a kernel, fibre, section or centralizer does not assign every algorithm answering that operation to research. Generic image/word-section computations consume supplied generators; they do not discover arithmetic parabolic generators or their liftable image. Research owns independent equation solvers, partial-isometry extension and pointwise-perpendicular kernels. The port owns specialized searches, compatibility-action construction, image discovery and recursive assembly. Neither source residency nor indefiniteness alone assigns ownership. See [port design](references/port-design.md#assigning-a-computation).
- **Use the concrete assignments.** [Implementation assignments](references/port-design.md#implementation-assignments) names the actual functions and repair owners on both sides. In particular, implement `IsotropicInduction.line_automorphism_generators`, `plane_automorphism_generators`, `flag_automorphism_generators` and the isotropic branch of `EichlerRecursion.pointed_automorphism_generators` here. Calling research's unresolved parabolic `.group_generators()` does not implement them. Consume research's independent `R.rational_lifts(...).integral_members()` and `R.pointwise_perpendicular_kernel()` computations. Move the copied binary cycle and scaled-shell multiplier search here, and implement Eichler approximate/cover/full-group assembly here; frozen preamble names may dispatch to those implementations. Research's matching instructions are in `/home/dzack/research/AGENTS.md`, “Concrete boundary with sage-indefinite-port”.
- **No new tests.** The acceptance suite (tests/acceptance, built in #33) is the only suite and is frozen. No unit adds, edits or weakens a case. A capability lands by removing the xfail markers that cite its unit. Changing a case, or changing the frozen API, means reopening #33 with the reason.
- **State questions as mathematics.** Orbit, stabilizer and equivalence questions are Isom and Aut of marked lattices: diagrams of lattice morphisms ending at L, with natural isomorphisms as morphisms (references/port-design.md, "Marked lattices"). A pointed lattice ⟨n⟩ → L lives in the coslice ⟨n⟩/Lattices (the marking is fixed pointwise); sublattices and configurations live in the arrow category Ar(Lattices) (the marking may move by an automorphism of its domain); flags live in Fun([2], Lattices). Do not introduce ad-hoc vocabulary ("transporter", "stabilizer data", "full generating data") for constructions that have standard names.

# Architecture and dependencies

- **Sage runtime on this host**: Docker is installed by Snap at `/snap/bin/docker`, not `/usr/bin/docker`.
  The `sage-dev-develop` container is already running; reach it with `/usr/bin/sg docker -c '/snap/bin/docker ...'`.
  An unsuccessful `command -v docker` does not mean Docker is unavailable.
  The preamble Julia bridge requires `/home/dzack/.juliaup/bin/julia`.
  The tracked `.envrc` configures `PATH`, `SAGE_BIN`, `SAGE_PYTHON_BIN`, `SAGE_PREPARSE_BIN`, `SAGE_JULIA_COMMAND` and `PYTHONPATH` for the container-backed Sage and the central QC recipes.
  Test runs use the launch in Running Sage, not `.envrc`'s working-tree paths.
- **The `research` preamble is the lattice substrate and a co-developed dependency**:
  This repository consumes `dzack_research.preamble` (categories, formed modules, lattices, isometries, discriminant modules) from the checkout at `/home/dzack/research/src`. `.envrc` puts it on `PYTHONPATH` and `MYPYPATH`; `pyproject.toml` declares `dzack-research`.
  - **Use the preamble's own classes and types**: annotate with `Lattice` (and `Lattice.Element`), `LatticeIsometry`, `LatticeIsometryHomset`, and the tensor and presented-module classes. Never add Protocol, stub, or facade types over the preamble, and never edit `sys.path`.
  - **Fix preamble defects upstream**: a missing annotation, a class reachable only from a private module, or a wrong result in the preamble is fixed in `research` under its own QC **by its exclusive repository owner**; the Sage port worker does not edit, stage, or commit in research while its separate workstream is active (see "Working against the frozen API"). Port code is written against the corrected or planned API meanwhile and does not wait. No local workaround, no fallback to raw SageMath lattices (their indefinite support is insufficient), no reinvention of lattice, discriminant, or formed-category infrastructure here.
  - **Ownership**: `research` owns mathematical interfaces and independently reusable primitive computations. This repository owns approximate models, coverings, indefinite recursion and traversal, binary edgewalk cycles, scaled first-success searches and arithmetic parabolic assembly. Research may dispatch public methods to these port algorithms. Invariant decomposition and gluing from supplied component groups are research primitives; non-scalar centralizers need an identified computational route and cannot be assigned to research merely by their name. The complete centralizer consumer obligation remains in #18. The operation-by-operation split is in `references/port-design.md`.
  - The preamble surface to consume, not reinvent, is listed in [references/preamble-api.md](references/preamble-api.md).
    The architecture policies every edit obeys are in [references/contribution-policies.md](references/contribution-policies.md).
- **NEVER BUILD THE C++ CODE**:
  **Do NOT attempt to compile, build, configure, or invoke C++ compilation toolchains** for `polyhedral_common` or any other C++ source in `references/`.
  The upstream C++ codebase in `references/polyhedral_common` is strictly a reference implementation for algorithm extraction, logic translation, and structural understanding.
  All production algorithms in this repository MUST be implemented as native SageMath / Python code delegating low-level operations to standard system backends (FLINT for exact integer linear algebra, Normaliz / cddlib / PPL for polyhedral cones). Research supplies independent primitive computations, including exact module algebra, specified equation fibres, partial-isometry extension, finite configurations and finite actions from supplied effective data. This port supplies specialized indefinite Isom, Aut, orbit, search and assembly algorithms. Reuse generic GAP/Bliss/backend operations through their owners. A claimed primitive gap must state its independent contract and input data; the output of this port cannot be assumed as its own prerequisite.

# What is built here

Consume the existing foundational formed algebra. Development in this repository focuses on:

1. **Lorentzian perfect-domain reduction**: Traversal of the $(1,n)$ reduction complex for $O^\Omega(L)$ and $O(L)$ generators.
2. **$2U$-Eichler approximate models & Higher-Witt-index recursion**: Generating $A(L)$ and lifting covering sets to full $O(L)$ for Witt index $\geq 2$.
3. **Parabolic recursion & Isotropic orbit decomposition**: Compose research's group actions, sections and kernels with recursive reduced-group computations to obtain cusp and flag stabilizers and the frozen rank-at-most-two isotropic orbits.
4. **Centralizer and arithmetic-subgroup consumers**: Supply indefinite orthogonal-group and orbit computations to research's general centralizer, gluing and subgroup operations; verify all frozen consumer cases.

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

## Technical continuity

An automatic continuation does not reset the current mathematical task. Reuse
already established issue, source, and test context rather than reconstructing
it with repeated status queries, source reads, and memory calls. Retrieve a
specific project-memory entry only when a required fact is genuinely missing
or its recorded source assumptions have changed.

Record a durable checkpoint when a new mathematical result, source
correspondence, or failure materially changes the next implementation action;
do not create notebook entries as receipts for ordinary tool calls or probes.
Investigations must resolve a question needed for the next source correction
or complete consumer validation, not become a recurring activity of their own.
While a long test runs, continue independent pinned implementation work; its
result is valid only for the exact source revision and evidence recorded by the
runner required in TODO.md.

Plan work is card-backed. Create and update plan cards with `agent-memory plan add` and `agent-memory plan update`, not `agent-memory add --type plan`.

Use `agent-memory retrieve <key>`, `agent-memory update <key>`, and `agent-memory delete <key>` for memory CRUD.

The vault should be committed at all times. Treat staged or unstaged vault changes as an ephemeral error state. Before normal memory work resumes, load the bundled vault-maintenance skill with `agent-memory maintain skill vault-maintenance` and follow its referenced check, repair, and commit workflows.

Move reusable lessons during maintenance with:

```bash
agent-memory maintain move <key> --to global/advice
```
<!-- agent-memory:end -->
