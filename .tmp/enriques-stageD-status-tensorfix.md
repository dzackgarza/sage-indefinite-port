# Stage-D after the tensor-hash repair: terminal and instrumentation records (run state: atom-B ended, capped)

## Provenance of what ran

**Launch-time evidence, i.e. what the runs themselves recorded.** Both launchers bind the library checkout with proot's
plain `-b` (`PROV=/workspace/nucleus/work/projects/p-8c69fc…/repo` → `/repos/sage-extensions`, PYTHONPATH
`/repos/sage-extensions/src`). A plain `-b` is **not** a read-only mount (proot mounts read-write by default);
the recorded content is the launch content because nothing wrote through that path, not because writing was
impossible. The same holds for the `/probe.py` binds of the probe sources themselves, and each run header prints the runtime module path plus
`lattice_morphisms.py sha256[:16] = 00f0c682880cc022`. **No commit OID was recorded at launch**, and the logs do not
print a `tensor.py` digest. A run is therefore pinned by those recorded digests (and the bound path), not by a commit
name.

**Repair identity is a repository object, not a launch reading.** `c2c81a3` is the tensor-repair commit, and the
`tensor.py` digest quoted here — sha256 `8da116e10f939e86eb3f3ef849d8f541b3c1f763f254042efd703543eaa1d029` — is the
object *at that commit*, verified by hashing `c2c81a3:src/sage_extensions/tensors/tensor.py`. It identifies the code
that must be present, not the instant a job started.

**Later record commits are not launch pins.** `50de3c1` **postdates the original job launch**: it is the later
equivalent record commit whose tree contains the same tensor repair. An earlier version of this note named it as the
provider "that ran", which is wrong; where a commit name is wanted for the *repair*, use `c2c81a3`, and where the pin
of a *specific run* is wanted, quote that run's own header digests. Likewise the port commit
`3f24ab905d56bad2e009bacb38b9f98dc0b318e6` (`refs/heads/main`, remote-verified) is the artifact commit, not a
launch reading; the stage-D witness-scope correction `d5507b07` (det +1, proper subgroup, index unproved) stays
evidence on its own terms.

## What is measured where

- `.tmp/enriques_stageD_hashfix_tensorfix.log` — terminal log of the original call
  `action.stabilizer().group_generators()` on N = U+U(2)+E8(-2), line I_1_1 (job jmw0qw891):
  isotropic_stabilizer_action 17.0 s, then the stabilizer call hit the inner 2400 s cap (exit 124).
  PEAKRSS peak tree 2332 MiB of the 6144 MiB guard: **time-bound, no memory wall**.
- `.tmp/stageD-cost-atoms.log` (probe `.tmp/probe-stageD-cost-atoms.py`) — first separation:
  **ATOM A finished: `finite_representation()` 154.8 s with an honest orbit of exactly 1024 points**
  (poller: 1024/1024 positions stable from ~t=175 s). The later stabilizer stage did not end by its
  960 s cap. Before the repair the orbit collapsed to 1 point and hid all of this.
- `.tmp/stageD-atomB.log` (probe `.tmp/probe-stageD-atomB.py`, job **jmw0t49hp, ended**) — the completed
  disaggregation of the post-BFS stages, with its **actual** figure in place of the earlier sample estimate:

  | atom | measured | note |
  | --- | --- | --- |
  | A `finite_representation()` | 156.6 s | orbit of 1024 points |
  | B1 `permutation_representation` | 0.0 s | codomain `Symmetric group on 1024 points` |
  | B2 `finite_stabilizer` | 1.8 s | — |
  | B3 `dims` | 1024 points x 25 letters = **25600 words** | letter images cached |
  | B3 `orbit_schreier_morphism` | **1520.2 s**, finished at **t = 2042.3 s** | supersedes the ~662 s first-200-row estimate (25.8 ms/row at row 200 was an early sample, not the average) |
  | B4 `preimage_subgroup` construct | 0.0 s | at t = 2042.3 s |
  | B4 generator realization | **did not finish by the 2160 s cap** | terminal dump in Sage parent/classcall/`coerce_dict` construction |

  The run ended exit **124** after **36m1s**; the run manager stopped its **1 descendant**. Guard block from the
  same raw log: `child waited=124 shell_code=124 samples=1080`, **peak tree RSS 1660 MiB of the 6144 MiB guard**,
  `low-water MemAvail 10943 MiB` — time-bound throughout, not a memory wall, and the low-water reading is the one
  that matters for that claim.

  The in-band row-split accounting it also printed (`act` 18.0 ms, `section(act_p)` 0.0 ms, `section(p)` 0.0 ms,
  word ops 0.0 ms per row) shows the cost is in the action evaluation itself, not in the word bookkeeping around it.

**Who owns what now.** The stage that remains is `g_sets.py` (`orbit_schreier_morphism` and the realization that
follows it), and the central repair there is being implemented in `sage-extensions` *after* this job ended. It is not
duplicated in this port, and no further heavy probe is started here. The port's consumer integration resumes once
that repair is measured and available.

## Boundary

Semantic content is repaired (port `3f24ab90`: on U+U(2) every supplied generator's actual action and
its converted integral output compare correctly, the stabilizer family carries the -e element, small
model). **Small-model success is not full acceptance**: on N the BFS finishes (1024 points), B3 is measured at
1520.2 s, and the run is capped inside B4 generator realization. 5/10/19 orbit acceptance and S1 stay open **by
name**; nothing here substitutes for them, and no partial stage is reported as an acceptance.

The remaining stage is an owning-library item in `g_sets.py`, whose central repair is now being implemented in
`sage-extensions` after this job ended. This port does not duplicate that repair and starts no further heavy probe;
its consumer integration resumes once the repair is measured and available.

## Terminal of the retained-transition acceptance rerun (`jmw0urc9n`, cap unchanged)

- `action.stabilizer()` **completed: 1898.5 s, 25600 generators** (the full $1024 \times 25$ Schreier family): B3 words through the retained transition table, then B4 realization of every preimage word into a group element. No second sweep, no per-row re-acts.
- The run then entered the probe's own conversion loop `tuple(N.integral_isometry(g) for g in ...)` and hit its **unchanged** 2400 s inner cap there (exit 124), PEAKRSS peak tree 3637 MiB of the 6144 MiB guard: again **time-bound, no memory wall**. The first boundary the conversions own (`integral_restriction -> factor_through_or_none -> _preimage_or_none`) is thus measured real, not from a sample: 25600 per-element conversions, one preimage computation each; the conversion stage, not the stabilizer realization, is now the acceptance frontier.
- The raw log `.tmp/enriques_stageD_schreierfix_6675f1f.log` is published with this note. The stabilized Schreier family size (25600) is itself the cost driver of that stage; deducing semantics from it duplicates nothing.
