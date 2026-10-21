# Handoff

The rules for using this file are in AGENTS.md "The turn protocol".
Read **Now** and act on it.
Replace Now when a step finishes or a result arrives; never append history.
The full remaining obligation is in TODO.md, and Now never narrows it.

## Now

**Coordinator:** the Claude session that rewrote this protocol is running the port with parallel agents, one per file: edgewalk (`edgewalk*.py`), recursion (`recursive.py`), parabolic (`isotropic_induction.py`), Eichler (`eichler.py`), and perfect cells plus orthogonal groups (`lorentzian_cells.py`, `indefinite_reduction.py`, `groupoid_cache.py`, `marked_lattices.py`). Do not edit those files while it runs.

**Research fixes in progress** (uncommitted, in `/home/dzack/research/.tmp/wt-port-fix`): small-family `engine_permutation`; general-codimension `witt_extension`; `invariant_module_construction` without `is_finite()` on infinite groups; setwise sublattice stabilizer on `O_plus()`. Still open: `right_coset_section` re-admission (`g_sets.py:682`), non-spanning families in `_coordinate_basis_data`, E10 orbit-span speed.

**Acceptance:** CI Acceptance on each push, in `ghcr.io/dzackgarza/sage-indefinite-port-sage:py314`.

## State

- **Committed:** all inherited source, in `df7d45b` (WIP, claims nothing).
  The helper checks formerly in the frozen perfect-domain file are in `.tmp/probe/perfect_domain_helper_checks.sage`; the four scratch tests from `tests/` are in `.tmp/probe/`.

- **Best evidence for A_PERFECT:** `.tmp/accept-perfect-final-20261021g.log` (source equal to `df7d45b` apart from removed type-ignore comments, `-x`): total-mode 0–39 and isotropic 0 passed, isotropic 1 failed as in step 1.

- **CI:** `.github/workflows/acceptance.yml` runs every family on every push; `review-pr.yml` still fails (`SAGE_BIN must be set`) and is not the acceptance signal.

- **Provider:** CI installs research `main` from GitHub; local checks use `/home/dzack/research/src`.

## Open decisions (user)

- The local push gate (`just test-push`, `_sage-pytest` in `~/ai-review-ci/justfiles/sage.just`) still runs the whole acceptance suite before every push.
  Removing it there was deferred because that file holds another writer's uncommitted edits.

- Whether the moved helper checks in `.tmp/probe/` should become frozen cases through #33, or stay probes.
