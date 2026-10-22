# Sage indefinite lattice port — agent instructions

This project translates the specialized indefinite-form algorithms of the
[committed pinned sources](references/vendor/README.md) into Sage/Python,
consuming the mathematical operations of [research](https://github.com/dzackgarza/research).
It is **a port of existing algorithms**, not an algorithm-design project.

## Sources of authority

| Source | Responsibility |
| --- | --- |
| [Issue #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4), [#33](https://github.com/dzackgarza/sage-indefinite-port/issues/33), [#15](https://github.com/dzackgarza/sage-indefinite-port/issues/15), [#18](https://github.com/dzackgarza/sage-indefinite-port/issues/18), [#36](https://github.com/dzackgarza/sage-indefinite-port/issues/36) | Original obligation, ownership, issue dependencies and acceptance |
| [TODO.md](TODO.md) | Remaining algorithm, dependency and acceptance obligations; a checked node needs source-backed evidence |
| [`.tmp/HANDOFF.md`](.tmp/HANDOFF.md) | One concrete **Now** action, its source/consumer and latest known remote revision; replace obsolete Now rather than append a history |
| [Port design](references/port-design.md), [assignments](references/port-design.md#implementation-assignments), [translation map](references/source-to-sage-translation-map.md) | Frozen API, owner split, exact upstream routine locations |
| [Vendored source and Allcock paper](references/vendor/README.md) | **Algorithm authority**; translation notes never replace source |
| [Acceptance workflow](.github/workflows/acceptance.yml), [oracle manifest](tests/fixtures/oracle_manifest.yaml) | Frozen executable cases and results for exact published revisions |

## Entry and continuity

For a **new cloud checkout**, work from
`implement/lorentzian-completion` (not default `main`).
Read this document once, then `.tmp/HANDOFF.md` **Now** and the relevant
TODO/issue/source sections. Fetch the branch and identify its current commit and
newest *completed* Acceptance run before changing it. If the handoff commit or
selected failure predates newer code/results, replace Now with the actual next
unresolved implementation step; do not enact its stale instructions.

A normal continuation starts with the current Now action and its precise
source/consumer, rather than repeat source-wide surveys, re-read the entire DAG,
or reselect work. The handoff is **not** proof that named people or agents remain
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
in the TODO acceptance table, which identifies where sample evidence applies.
Record selected case IDs, covered and uncovered branches, ranks/signatures,
passing results and the exact commit. A sample must not avoid a genuinely
distinct case branch or alter any assertion. Full-corpus failures outside the
sample remain defects to record, **not** reasons to declare the sample or
the entire port correct. Long-tail exhaustive coverage does not gate the
agreed milestones.

GitHub [Acceptance](.github/workflows/acceptance.yml) runs the full frozen
suite in sharded Sage jobs on every push using the
`ghcr.io/dzackgarza/sage-indefinite-port-sage:py314` image; the run
installs the published `research` branch and records its commit.
Use the exact run/artifact for the revision under discussion. Job status
alone cannot prove a case passed. A full workflow can be red while the
explicitly cited branch-representative samples pass; conversely, an
interrupted/missing shard cannot be counted as passed. The final revision
requires its chosen case IDs to pass under the same revision/provider.

In a capable local Sage environment, run one focused node or diagnostic
to locate the next failing source step. That is diagnostic evidence, not
acceptance. Do not launch an entire family locally, repeatedly poll queued
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
  mandatory PR review. No type-debt exception trailer or `--no-verify`
  workaround is required for these old gates.
- Type debt belongs to [#36](https://github.com/dzackgarza/sage-indefinite-port/issues/36)
  after algorithms/consumers have settled. CI may report types/style/debt,
  but during development such findings are retained, **nonblocking**
  diagnostics; they do not divert the current mathematical implementation.
- There are currently no active `@unfinished` wrappers or acceptance
  xfails. Do not invent marker-removal work or infer readiness from their
  absence. Public capability is judged by actual consumer behavior.
- Treat a red **actual local syntax check** as an immediate defect in the
  committed source. Treat other CI failures as evidence to route to their
  algorithm/provider owner; investigate concrete failing cases instead of
  performing broad cleanup, optimization or repetitive diagnostics.
- The final issue #4 closure requires the TODO acceptance evidence and public
  research-to-port consumers, with a recorded port and provider revision.
  `just test-ci` is a nonblocking diagnostic collection, not a substitute
  for the frozen Acceptance workflow.

Do not require a local agent-memory vault, home-directory paths, or a running
coordinator to resume cloud work. Durable mathematical decisions belong in
the owning source/issue; the tracked handoff records only the next actionable
step, not an activity log.
