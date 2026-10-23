# Cloud handoff — sage-indefinite-port

This is the tracked entry point for the
[`implement/lorentzian-completion`](https://github.com/dzackgarza/sage-indefinite-port/tree/implement/lorentzian-completion)
branch. It is a **single-action locator**, not a coordinator reservation,
activity journal, substitute for [TODO.md](../TODO.md), or claim that older
local worktrees or agents are still active. Full ownership and source-fidelity
rules are in [AGENTS.md](../AGENTS.md).

## Now

**The functional MVP is delivered for the pinned source algorithms.** Its
small public-consumer examples are recorded in
[TODO.md](../TODO.md#functional-mvp-evidence-from-the-current-implementation-work).
The former instruction to invent a higher-index route for rational Witt
index less than two was incorrect and is withdrawn.

The complete pinned Julia/GAP/executable chain has the same two-rational-plane
requirement as the port. The reference paper's Theorems 4.4–4.6 establish the
general Aut/Isom recursion for `max(p,q) >= 5`, which ensures that every
higher-index recursive stage has the two required planes. The rank-four
forms `diag(2,2,-6,-6)` and `U+<2,-4>` lie outside that guarantee. The reference
contains no separate branch for them; its `QmatExt` remains in the same
rational quadratic space. The port's reduction, envelope, covering families,
recursive stabilizers and integral-witness corrections are present. See the
[source-backed applicability correction](../TODO.md#upstream-applicability-correction-to-the-former-mvp-blocker)
for the public dispatch and the separate vector-orbit hypotheses.

The native centralizer gluing, general CM image, primitive embedding/spacelike,
finite-image, coisotropic lift, cusp equality, zero-rank root and native
Witt-index repairs are published in research `2024cb764e4e28ceb818f82f2b1075291ce83a68`
and its ancestors. Port `59227592439cdccd859e9282bc02ced2a977a3c4` contains the
spacelike integration, empty plane/flag decision and earlier syntax,
signed-view, parabolic and anisotropic-zero repairs. The source comparison
was made at port `bccf3dcae47bfb16a70dce73e55b7247865546f7`; subsequent
documentation commits do not change those algorithms. Fetch current tips
before editing. Full frozen regression acceptance and #36 refinement remain
separate, deferred work; neither is claimed complete by the MVP delivery.

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
contracts. The former low-Witt entry is not an outstanding port branch.
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
