# sage-indefinite-port

Port of indefinite lattices and orthogonal group algorithms to SageMath.

## Documentation and References

- [CITATIONS.md](file:///home/dzack/gitclones/sage-indefinite-port/CITATIONS.md): Formal academic citations, literature references, and BibTeX entries.
- [oracle_manifest.yaml](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/oracle_manifest.yaml): Acceptance criteria and oracle manifest covering 13 test suites.
- [references/README.md](file:///home/dzack/gitclones/sage-indefinite-port/references/README.md): Overview of upstream software subtrees (`polyhedral_common`, `Indefinite.jl`).

## Quality Control

Run quality control checks with `just`:

```bash
just test-commit  # Fast commit-tier QC (syntax, types, formatting)
just test-push    # Full test suite and lint checks
just test-ci      # CI acceptance checks
```

