# Cloud handoff — sage-indefinite-port

This is the tracked entry point for the
[`implement/lorentzian-completion`](https://github.com/dzackgarza/sage-indefinite-port/tree/implement/lorentzian-completion)
branch. It is a **single-action locator**, not a coordinator reservation,
activity journal, substitute for [TODO.md](../TODO.md), or claim that older
local worktrees or agents are still active. Full ownership and source-fidelity
rules are in [AGENTS.md](../AGENTS.md).

## Now

**First, repair a confirmed source importability defect on the published branch.**
At the last source-bearing commit `69a456e`, parsing
`src/sage_indefinite_port/indefinite/lorentzian_cells.py` fails at line
2052: `except ArithmeticError, OverflowError, RuntimeError, TypeError, ValueError:`
is invalid Python syntax. The exception types require parentheses. This was
confirmed against the published code at `2777a40` using `ast.parse` without
running Sage or touching the frozen tests. Apply the local syntax repair,
check that the file parses, and commit/push it before attempting an acceptance
case. Do **not** reinterpret this one-line syntax error as a missing algorithm.

**Selected obligation: #15 EDGEWALK / A_REFLECTIVE.** Use the *current*
branch and its latest **completed** Acceptance report to isolate a failing
branch-representative reflective case in
`tests/acceptance/test_reflective.sage`, starting at the
`test_simple_root_counts_match_polyhedral_common` consumer when still failing.
Trace that case through
`src/sage_indefinite_port/indefinite/edgewalk.py` and
`edgewalk_extensions.py` to the
[translation map](../references/source-to-sage-translation-map.md) and the
pinned `Lorentzian.g` / Allcock binary-cycle source.
Correct the precise algorithm/callee divergence; keep the frozen test unchanged.

Before changing EDGEWALK, compare the remote `port/edgewalk` branch with the
integration branch **for the affected routines only**. It has an unmerged
source correction at `386f9b8`; do not discard or blindly merge it. The
port algorithm needs starting corners, cusp/root candidates, next-wall selection,
oriented binary-cycle stopping and reflectivity/termination, with exact
source-correspondence. Record a representative set covering the pinned
algorithm's ranks, signatures, and branches, not simply a passing subset.

**Diagnostic starting point, not current verdict:** the fully reported
[Acceptance run 39667397208](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/39667397208)
on port `30799df`, research `0a70597c04` exposed three reflective error
families: domain-construction `ValueError`, string-formatting `TypeError`,
and nonsquare-matrix `ValueError`. Source and provider have changed since;
check whether each failure remains on the latest committed revision before
editing. Fix one demonstrated route rather than infer thousands of distinct
algorithm bugs from the number of repeated cases.

## Published starting state

- Last source-bearing integration revision before the documentation handoff:
  `69a456e419a2be42da682039addc5fbe79b861eb` (2026-10-22).
  Subsequent documentation commits may advance the branch without changing
  that source. Fetch the **current remote tip** before editing; if a newer
  source revision has repaired the import defect or changed the selected
  failure, use the new evidence instead of repeating the old step.
- [Latest branch Acceptance runs](https://github.com/dzackgarza/sage-indefinite-port/actions/workflows/acceptance.yml):
  the recorded latest run for `69a456e` is
  [39832125180](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/39832125180).
  It was pending when the handoff was prepared; **do not count pending or
  missing shards as passed**. The newest *completed* run known at the time was
  [39702377287](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/39702377287)
  on earlier source and failed. The detailed older source-backed failure
  classifications above are from `39667397208`.
- `A_PERFECT` is [accepted at `30799df`](https://github.com/dzackgarza/sage-indefinite-port/actions/runs/39667397208)
  (all eight shards, 80 cases). Earlier reports show 158 isometry cases with
  no case-level failure, **but one isometry shard reported failure**, so
  `A_ISOMETRY` is not recorded as accepted. Other acceptance rows in TODO.md
  remain open until their evidence is recorded; do not promote diagnostic
  counts to a closure.
- The pinned [C++/Julia/GAP sources and papers](../references/vendor/README.md),
  complete [oracle manifest](../tests/fixtures/oracle_manifest.yaml),
  [public API](../references/preamble-api.md) and
  [issue/contract details](https://github.com/dzackgarza/sage-indefinite-port/issues/33)
  are tracked and available remotely. None of the ignored local
  `.tmp/probe` files, old shell logs, Zotero copies, local agent-memory vault,
  or host-specific Docker paths is needed to locate the next source step.

## Remaining work and acceptance

The complete outstanding dependency graph is [TODO.md](../TODO.md) and the
[issue #4 tree](https://github.com/dzackgarza/sage-indefinite-port/issues/4).
After EDGEWALK, unresolved port obligations include the #18 parabolic,
recursive, orbit, boundary and centralizer consumers and research-owned
`FINITE`/`LIFTS` provider contracts. **Do not assume an outstanding TODO
node means its entire implementation is absent.** Identify the precise
remaining failing branch before authoring.

The [October 22 amendment to #33](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6343829665)
replaces exhaustive-case completion for large families with evidence from a
**branch-representative sample**. Preserve every frozen fixture and assertion;
include IDs of cited passing cases, branches exercised, cases/branches not
covered, their commit, and provider revision. The CI workflow still runs
full long-tail coverage and may fail outside the cited sample. Successful
sample coverage does not excuse a known uncited defect, and full-suite failure
alone is not the completion criterion.

Several *remote* branches contain commits not ancestors of the integration
tip: `port/cells`, `port/edgewalk`, `port/eichler`, `port/parabolic`,
`finish/handoff-remainder`, and `implement/pinned-indefinite-reduction`.
These are **potential source deltas, not active-worker ownership claims**.
Compare only the code required by the selected obligation, and integrate
applicable proven work instead of restarting it. The `port/recursive` tip
was an ancestor of `69a456e` at preparation time.

## Remote execution boundary

Use GitHub Actions for complete acceptance and its retained run artifacts;
source is published, and CI provisions the image and current `research` main.
The local `just test-commit` and `just test-push` commands are **syntax-only**
development checks. These check the staged or outgoing changes respectively,
not necessarily preexisting syntax defects elsewhere in the source tree.

After a real source correction and a pushed revision, update **Now** to the
next failing source/consumer with its exact commit and Acceptance run; replace
the existing action rather than append logs, dated status summaries or a
second plan. Never end at a probe or a report when the assigned algorithm
still needs implementation.
