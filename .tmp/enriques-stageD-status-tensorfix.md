# Stage-D after the tensor-hash repair: terminal and instrumentation records (run state: atom-B live)

## Provenance of what ran

- Provider `sage-extensions` `50de3c1` (contains the tensor repair `c2c81a3`; its recorded
  `src/sage_extensions/tensors/tensor.py` sha256 `8da116e10f939e86eb3f3ef849d8f541b3c1f763f254042efd703543eaa1d029`,
  runtime import `/repos/sage-extensions/src/sage_extensions/tensors/tensor.py` under the proot binding;
  the attribution probes also print `lattice_morphisms.py sha256[:16] = 00f0c682880cc022` per run).
- Port `3f24ab905d56bad2e009bacb38b9f98dc0b318e6` (`refs/heads/main` of the port remote, verified
  2026-11-03; the earlier working note that it was local-only was stale). The stage-D witness-scope
  correction `d5507b07` (det +1, proper subgroup, index unproved) stays evidence.
- Launch/provenance records: provider `.tmp/launch_stageD_tensorfix.sh` (tracked in `5c6ad9b`),
  `.tmp/launch-stageD-cost-atoms.sh`, `.tmp/launch-stageD-atomB.sh` here; every run under
  `PEAKRSS_LIMIT_MIB=6144`, one heavy slot, probes bound read-only as `/probe.py`.

## What is measured where

- `.tmp/enriques_stageD_hashfix_tensorfix.log` — terminal log of the original call
  `action.stabilizer().group_generators()` on N = U+U(2)+E8(-2), line I_1_1 (job jmw0qw891):
  isotropic_stabilizer_action 17.0 s, then the stabilizer call hit the inner 2400 s cap (exit 124).
  PEAKRSS peak tree 2332 MiB of the 6144 MiB guard: **time-bound, no memory wall**.
- `.tmp/stageD-cost-atoms.log` (probe `.tmp/probe-stageD-cost-atoms.py`) — first separation:
  **ATOM A finished: `finite_representation()` 154.8 s with an honest orbit of exactly 1024 points**
  (poller: 1024/1024 positions stable from ~t=175 s). The later stabilizer stage did not end by its
  960 s cap. Before the repair the orbit collapsed to 1 point and hid all of this.
- `.tmp/stageD-atomB.log` (probe `.tmp/probe-stageD-atomB.py`, run jmw0t49hp still live with its cap)
  — disaggregation of the post-BFS stages. Measured so far: `permutation_representation` 0.0 s;
  `finite_stabilizer` 1.8 s; then `_orbit_of`'s second orbit traversal (with word retention) in-band
  ~340 s before the words comprehension; words tuple dimensions **1024 points x 25 letters = 25600
  rows** at ~30-37 ms/row by the first 100 rows. The whole cost is **not** assigned to one internal
  operation until this raw log is complete.

## Boundary

Semantic content is repaired (port `3f24ab90`: on U+U(2) every supplied generator's actual action and
its converted integral output compare correctly, the stabilizer family carries the -e element, small
model). **Small-model success is not full acceptance**: on N the BFS finishes (1024 points) and the
later stabilizer stage is time-bound; 5/10/19 orbit acceptance and S1 stay open by name. The
stabilizer stage's cost is an owning-library item (`rational_integral_stabilizers.py` /
`g_sets.py`), which the atom-B log is measuring before any change.
