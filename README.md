# Sage indefinite lattice port

Native Sage/Python implementations of the indefinite-lattice isometry, orthogonal-group, orbit, Lorentzian perfect-domain and Allcock edgewalk algorithms from pinned upstream sources. The public mathematical lattice objects and operations are supplied by [research](https://github.com/dzackgarza/research).

## Cloud-agent entry point

The active integration branch is [`implement/lorentzian-completion`](https://github.com/dzackgarza/sage-indefinite-port/tree/implement/lorentzian-completion), **not** default `main`.

```bash
git clone https://github.com/dzackgarza/sage-indefinite-port.git
cd sage-indefinite-port
git switch --track origin/implement/lorentzian-completion
```

Read [AGENTS.md](AGENTS.md) for ownership and source-fidelity constraints; then read the single action in [the tracked handoff](.tmp/HANDOFF.md). Follow its relevant [TODO.md](TODO.md) node and the [issue tree rooted at #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4). The handoff is a starting locator, not evidence that its recorded commit is still the remote tip: fetch and compare GitHub before editing. Do not infer another agent owns a file from an old handoff or a branch name.

## Exact source and acceptance material

- [Port design and implementation assignments](references/port-design.md)
- [Pinned routine-to-Sage translation map](references/source-to-sage-translation-map.md)
- [Committed upstream source trees and papers](references/vendor/README.md)
- [Oracle manifest](tests/fixtures/oracle_manifest.yaml), [fixture provenance](tests/fixtures/DATA_SOURCES.md), and [frozen acceptance tests](tests/acceptance/)
- [Academic citations](CITATIONS.md) and [review criteria](references/review-guidelines.md)
- [Acceptance workflow](.github/workflows/acceptance.yml) and [Acceptance runs](https://github.com/dzackgarza/sage-indefinite-port/actions/workflows/acceptance.yml)

The CI workflow provisions Sage and installs the published `research` dependency; it runs the full frozen fixture suite on pushes. Each run is tied to its own port commit and recorded research revision. The accepted development criterion for large families is **branch-representative sample evidence**, not completion of all long-tail cases; see [issue #33's October 22 amendment](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6343829665) and TODO.md. No fixture or assertion is deleted to form a sample.

## Local development checks

The shared [ai-review-ci](https://github.com/dzackgarza/ai-review-ci) checkout is used only when the following delegated `just` recipes run locally:

```bash
just test-commit  # staged Python/JSON/TOML syntax
just test-push    # outgoing committed Python/JSON/TOML syntax
```

Neither command performs a full test suite, typing pass, format pass or PR review. Complete acceptance is in GitHub Actions; additional quality reports are nonblocking development diagnostics. A cloud checkout need not reproduce the original Linux machine's local Docker, Zotero, agent-memory vault or `.tmp/probe` logs. Use the committed pinned sources and GitHub CI instead.
