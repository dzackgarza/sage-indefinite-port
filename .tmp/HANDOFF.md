# Handoff

The rules for using this file are in AGENTS.md "The turn protocol".
Read **Now** and act on it.
Replace Now when a step finishes or a result arrives; never append history.
The full remaining obligation is in TODO.md, and Now never narrows it.

## Now

**Unit:** #15. **Node:** CELLS → A_PERFECT (TODO.md row A_PERFECT). **Selection:** already made; do not run `itree next` until #15 is closed.

**Acceptance:** the first CI Acceptance run is for the push of `implement/lorentzian-completion` made on 2026-10-21 (AGENTS.md "Running Sage" says how to find and read it). Until it finishes, work from the known failure below.

**Steps, in order:**

1. Repair the known A_PERFECT failure.
   `test_isotropic_mode_counts_match_polyhedral_common[Result_Enumeration_1]` raises `ArithmeticError: the exact-distance isotropic shell returned a nonisotropic vector` in `_ambient_vectors` (`lorentzian_cells.py`, called from `vectors_at_multiplier` / `first_shell` / `_search_initial_vectors` / `initial_cell`); traceback in `.tmp/accept-perfect-final-20261021g.log`.
   Hypothesis: the rewritten exact-distance shell (`_fixed_distance_vectors` and its LLL transport) diverges from the pinned isotropic initial-vector step, since the Oct 20 source passed isotropic cases 0–38 (`.tmp/accept-original-full-20261020i.log`).
   Read the pinned step `_search_initial_vectors` translates, run this one node locally (AGENTS.md "Running Sage"), locate the divergence, fix it there, commit, push.

2. When the Acceptance run finishes, download `acceptance-summary` and write its per-family counts here.
   Every A_PERFECT failure becomes the next step of this node.
   A failure in another family is that family's diagnostic pointer in its TODO.md row; it does not displace this node.

3. While no local check is needed, take up EXTEND (TODO.md row EXTEND, #18).
   Compare `_IsotropicExtension` in `src/sage_indefinite_port/indefinite/recursive.py` with items 1–2 of `.tmp/probe/extension-callee-correspondence.md`. Commit each correction and record the correspondence as source comments.

## State

- **Committed:** all inherited source, in `df7d45b` (WIP, claims nothing). The helper checks formerly in the frozen perfect-domain file are in `.tmp/probe/perfect_domain_helper_checks.sage`; the four scratch tests from `tests/` are in `.tmp/probe/`.

- **Best evidence for A_PERFECT:** `.tmp/accept-perfect-final-20261021g.log` (source equal to `df7d45b` apart from removed type-ignore comments, `-x`): total-mode 0–39 and isotropic 0 passed, isotropic 1 failed as in step 1.

- **CI:** `.github/workflows/acceptance.yml` runs every family on every push; `review-pr.yml` still fails (`SAGE_BIN must be set`) and is not the acceptance signal.

- **Provider:** CI installs research `main` from GitHub; local checks use `/home/dzack/research/src`.

## Open decisions (user)

- The local push gate (`just test-push`, `_sage-pytest` in `~/ai-review-ci/justfiles/sage.just`) still runs the whole acceptance suite before every push. Removing it there was deferred because that file holds another writer's uncommitted edits.

- Whether the moved helper checks in `.tmp/probe/` should become frozen cases through #33, or stay probes.
