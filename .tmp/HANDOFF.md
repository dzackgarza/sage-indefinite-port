# Cloud handoff — sage-indefinite-port

This is the tracked entry point for the
[`implement/lorentzian-completion`](https://github.com/dzackgarza/sage-indefinite-port/tree/implement/lorentzian-completion)
branch. It is a **single-action locator**, not a coordinator reservation,
activity journal, substitute for [TODO.md](../TODO.md), or claim that older
local worktrees or agents are still active. Full ownership and source-fidelity
rules are in [AGENTS.md](../AGENTS.md).

## Now

**The functional MVP is implemented and published against the pinned
upstream algorithms.** The latest code revisions are port
`666d5a50cd3501aca5ee07eb04f7e2bcabfcc672` and research
`ec87ab857183f91d5fa03ed1a6e2678667763411`. The continuation repaired actual
public-consumer omissions: Vinberg starting corners, radical Aut/Isom,
marked and arithmetic stabilizer generators, character-intersection cosets,
coslice subgroup transporters, and an isotropic transporter that lost its
isometry parent. Small native Sage examples pass, including the compact
six-wall hexagon, a nonreflective certificate, radical group images
12/72/48, and the `SO(2U) cap O_plus(2U)` boundary counts 1/2/2 with
stabilizer image 72. Exact evidence and limits are in
[TODO.md](../TODO.md#functional-mvp-evidence-from-the-current-implementation-work).

Do not restart the completed source-dispatch investigation or require the
full frozen acceptance matrix before recognizing this MVP delivery.
The retained source-shared extensions are open post-v1 requirements:
[low rational Witt index, #40](https://github.com/dzackgarza/sage-indefinite-port/issues/40)
and [degenerate ambient vector markings, #41](https://github.com/dzackgarza/sage-indefinite-port/issues/41).
Their explicit inputs and outputs remain required. Full frozen regression
acceptance and #36 refinement also remain open. Fetch current tips before
any subsequent implementation work; the broader project is not closed by
the MVP milestone.

## Historical published starting state

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
The functional evidence table records the implemented consumers and provider
contracts. The low-Witt functionality remains an outstanding requirement in
the explicit post-v1 backlog, including groups, integral witnesses, vector
and line classes, stabilizers, transporters and their downstream consumers.
These are functional implementation requirements, separate from the deferred
regression, performance, typing and refactoring work.
**Do not assume an outstanding TODO node means its entire implementation is
absent.** Identify the precise remaining branch before authoring.

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
