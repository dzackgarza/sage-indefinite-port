# sage-indefinite-port

Port of indefinite lattices and orthogonal group algorithms to SageMath.

## Quality Control

Run quality control checks with `just`:

```bash
just test-commit  # Fast commit-tier QC (syntax, types, formatting)
just test-push    # Full test suite and lint checks
just test-ci      # CI acceptance checks
```
