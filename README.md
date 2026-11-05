# Sage indefinite lattice port

Native Sage/Python implementations of the indefinite-lattice isometry, orthogonal-group, orbit, Lorentzian perfect-domain and Allcock edgewalk algorithms from pinned upstream sources. The public mathematical lattice objects and operations are supplied by [sage-extensions](https://github.com/dzackgarza/sage-extensions).

## Run it

Two even lattices of signature $(2,1)$ with the same discriminant form are isometric, and the port finds an isometry:

```python
>>> from sage_extensions.all import ZZ, Lattices
>>> A = Lattices(ZZ)([[0, 1, 0], [1, 0, 0], [0, 0, 2]])
>>> B = Lattices(ZZ)([[2, 1, 0], [1, 0, 0], [0, 0, 2]])
>>> A
Integral lattice of rank 3 and signature (2, 1)
>>> A.is_isometric(B)
True
>>> A.Isom(B).an_element()
Morphism from Integral lattice of rank 3 and signature (2, 1) to Integral lattice of rank 3 and signature (2, 1)

```

The block runs as a doctest, locally with `just readme-example` and in CI on every push as the
`readme` job of the [acceptance workflow](.github/workflows/acceptance.yml).

## Completion scope

The port must preserve every mathematical capability of the pinned sources,
including their general input domains, ranks, flags, branches, outputs and
certificates. General source routines are translated with their parameters
intact. Full upstream isotropic k-plane and flag functionality is required;
the current line/plane implementation leaves that translation incomplete.

Milestones and test samples order delivery and verification. They do not
authorize a smaller port. Every missing source-supported capability remains
required until implemented through its public consumers. See
[full upstream coverage](AGENTS.md#full-upstream-coverage) and
[the remaining parity obligation](https://github.com/dzackgarza/sage-indefinite-port/issues/44).

## Cloud-agent entry point

The active integration branch is `main`; the former `implement/lorentzian-completion` line is merged into it.

```bash
git clone https://github.com/dzackgarza/sage-indefinite-port.git
cd sage-indefinite-port
git switch --track origin/main
```

Read [AGENTS.md](AGENTS.md) for ownership and source-fidelity constraints; then take the next open work unit from the [issue tree rooted at #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4) with `uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port`. An issue is a starting locator, not evidence that its recorded commit is still the remote tip: fetch and compare GitHub before editing. Do not infer another agent owns a file from an old issue comment or a branch name.

## Exact source and acceptance material

- [Port design and implementation assignments](references/port-design.md)
- [Pinned routine-to-Sage translation map](references/source-to-sage-translation-map.md)
- [Committed upstream source trees and papers](references/vendor/README.md)
- [Oracle manifest](tests/fixtures/oracle_manifest.yaml), [fixture provenance](tests/fixtures/DATA_SOURCES.md), and [frozen acceptance tests](tests/acceptance/)
- [Academic citations](CITATIONS.md) and [review criteria](references/review-guidelines.md)
- [Acceptance workflow](.github/workflows/acceptance.yml) and [Acceptance runs](https://github.com/dzackgarza/sage-indefinite-port/actions/workflows/acceptance.yml)

The CI workflow provisions Sage and installs the published `research` dependency; it runs the full frozen fixture suite on pushes. Each run is tied to its own port commit and recorded research revision. The accepted development criterion for large families is **branch-representative sample evidence**, not completion of all long-tail cases; see [issue #33's October 22 amendment](https://github.com/dzackgarza/sage-indefinite-port/issues/33#issuecomment-6343829665) and the acceptance table on [#64](https://github.com/dzackgarza/sage-indefinite-port/issues/64). No fixture or assertion is deleted to form a sample.

## Local development checks

The shared [ai-review-ci](https://github.com/dzackgarza/ai-review-ci) checkout is used when the following delegated `just` recipes run locally:

```bash
just test-commit  # staged Python/JSON/TOML syntax
just test-push    # outgoing committed Python/JSON/TOML syntax
```

Neither command performs a full test suite, typing pass, format pass or PR review. Complete acceptance is in GitHub Actions; additional quality reports are nonblocking development diagnostics. A cloud checkout uses the committed pinned sources and CI without any private local machine state.
