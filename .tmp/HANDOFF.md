# Handoff

The rules for using this file are in AGENTS.md "The turn protocol".
Read **Now** and act on it.
Replace Now when a step finishes or a result arrives; never append history.
The full remaining obligation is in TODO.md, and Now never narrows it.

## Now

**Unit:** #15. **Node:** CELLS → A_PERFECT (TODO.md row A_PERFECT). **Selection:** already made; do not run `itree next` until #15 is closed.

**Running:** container PID 582706. Started 2026-10-21 06:07 UTC by the old runner `.tmp/accept-perfect-final-20261021g.runner.sh` against uncommitted source, with `-x`. Its log has had no output since 06:40, after total-mode `Result_Enumeration_37`. It predates the snapshot protocol and can accept nothing.
Whether to stop it is the user's decision (Open decisions).
Until its `.tmp/accept-perfect-final-20261021g.exit` exists, start no Sage process.

**Steps, in order:**

1. Restore the frozen suite.
   `tests/acceptance/test_perfect_domains.sage` holds 79 uncommitted lines of helper checks (`_integral_positive_vector_allmeth`, `_fixed_distance_vectors`). Copy those lines verbatim into `.tmp/probe/perfect_domain_helper_checks.sage`, then make the acceptance file match its committed content again with a forward edit.
   Move the untracked `tests/test_binary_split_vectors.sage`, `tests/test_edgewalk_exact_coordinates.sage`, `tests/test_edgewalk_positive_ray.py` and `tests/test_edgewalk_rank2_promised.sage` into `.tmp/probe/` unchanged.

2. Checkpoint-commit the inherited source in one WIP commit, so that a snapshot can run it:

   - every modified file under `src/` and `.envrc`;

   - the untracked modules `indefinite/indefinite_reduction.py`, `indefinite/groupoid_cache.py` and `indefinite/edgewalk_scalars.py`, which the modified modules import;

   - the pinned source `references/vendor/polyhedral_common@1592b246/src_{indefinite,isotropy,lorentzian}/` and `references/vendor/README.md`.

   Subject: `wip(#15,#18): checkpoint inherited perfect-domain, edgewalk and recursion source`. `git diff --cached --stat -- tests/` must be empty.
   Leave `.agents`, `.hermes`, `=2.18`, `cysignals_crash_logs/` and `tests/.report/` untracked.

3. Once PID 582706 has exited, and not before, launch an A run of `tests/acceptance/test_perfect_domains.sage` at that commit (AGENTS.md "Running Sage"). Write its run directory into **Running**.

4. While it runs, take up EXTEND (TODO.md row EXTEND, #18), which needs no Sage process.
   Compare `_IsotropicExtension` in `src/sage_indefinite_port/indefinite/recursive.py` with items 1–2 of `.tmp/probe/extension-callee-correspondence.md`. Commit each correction and record the correspondence as source comments.

5. When the A run's `exit` appears:

   - Exit 0 with all 80 cases passing: write the commit and run directory into TODO.md A_PERFECT `Accepted at`, then make the next CELLS-dependent family the next step.

   - Otherwise: the first failing node becomes a D-run target.
     Write its hypothesis here before launching it.

## State

- **Implemented, uncommitted until step 2:** the perfect-domain traversal in `lorentzian_cells.py` (+1,475 lines since `48670ca`, Oct 17); edgewalk, Eichler and recursion edits; the extension, kernel and lift calls in `recursive.py`. The recursion calls are present but unverified (TODO.md EXTEND, POINTED, BOUND).

- **Best evidence for A_PERFECT:** `.tmp/accept-original-full-20261020i.log`, on earlier uncommitted source: total-mode 0–39 and isotropic 0–38 passed.

- **Not run on any committed source:** every acceptance family.

- **Provider:** `/home/dzack/research` main `HEAD`. The old `research/.tmp/wt-port` worktree is gone; `.envrc` now points at `/home/dzack/research/src`.

## Open decisions (user)

- Stop PID 582706 or let it finish.
  It cannot accept anything, and it blocks every Sage run until it exits.

- Whether the moved helper checks in `.tmp/probe/` should become frozen cases through #33, or stay probes.
