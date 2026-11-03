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
- The raw log `.tmp/enriques_stageD_schreierfix_6675f1f.log` is published with this note, together with its launcher
  `.tmp/launch-stageD-schreierfix.sh` (`PEAKRSS_LIMIT_MIB=6144`, inner cap `timeout -s KILL 2400`, outer 2700). The
  stabilized Schreier family size (25600) is itself the cost driver of that stage; deducing semantics from it
  duplicates nothing.
- **Cleanup, from the run manager rather than the log:** one descendant was cleaned up after the KILL. `peakrss`
  records the process tree's peak and wait status, not the cleanup, so this fact is attributed to the manager's
  report and is not read out of the log.

### Provenance scope of this run, stated precisely

- The launcher binds the library working tree (`PROV=…/p-8c69fc…/repo` → `/repos/sage-extensions`), and the log's
  header prints the runtime module path plus `lattice_morphisms.py sha256[:16] = 00f0c682880cc022`. **That digest
  does not identify the repaired file**: this repair is in `g_sets.py`, and `lattice_morphisms.py` is unchanged by
  it. The log therefore does not by itself prove which `g_sets.py` ran; the run's repair identity is the
  operator-recorded tree `6675f1f`, which is why it is encoded in the log's filename. A future run should print the
  `g_sets.py` digest for the same reason the earlier ones print `lattice_morphisms.py`.
- **`6c1500ff` is not this job's pin.** It is the subsequently published provider repair (exact retained
  augmentation, one-shot-safe seeds, distinct Schreier-image presentation). It is what the changed-source guarded
  probe `jmw0whgrg`, already running in the library session under the same 2400 s / 6144 MiB bounds, is being
  measured against. Nothing in this note claims `jmw0urc9n` imported it.

### The exact consumer continuation (prepared, not yet run)

When `jmw0whgrg` reports, the port-side continuation is one chain and one check, with no new provider code here:

```
IsotropicInduction._stabilizer_generators            (port: isotropic_induction.py)
  -> recursion.isotropic_stabilizer_action(isotropic, radical_images)      # port, unchanged
  -> .stabilizer().group_generators()                                     # provider (repaired)
  -> port consumer: lattice.integral_isometry(g) per generator            # the measured frontier
```

Expected data from the provider: a finite ordered set of isometries fixing the lattice, of size 1024 x 25 = 25600 on
this input. The port's acceptance checks on it, exactly the ones already used:

- every returned generator is an integral unimodular isometry of `N` (`M G M^T = G`, integer rows);
- each fixes the line **setwise**, and the family contains at least one element acting as `-1` on it -- the
  sign-reversing witness `diag(-1,-1,1,...,1)` is an element of `Aut((N, Z e1))` and must be *in the generated
  group*, which a family acting `+1` throughout is not (index unproved, as recorded in
  `.tmp/handoff-port-enriques-stageD-hashfix.md`).

Two port-side obligations follow from the measurement, and neither duplicates the provider repair: the conversion
stage is the frontier (25600 per-element `integral_isometry` conversions, one preimage computation each), so the port
consumer must decide what it actually needs from that family instead of realizing all 25600 unconditionally; and the
original acceptance -- 5/10/19 orbits and S1 -- stays open on its own terms until the consumer verdict exists.

## Terminal of the deduped-presentation rerun (`jmw0whgrg`, cap unchanged): the in-line dedupe is cost-negative at N scale

- The rerun with provider `6c1500ff` (stabilizer presenting each distinct Schreier image once; the log filename
  carries the 7-character abbreviation `6c1500f`) hit its unchanged 2400 s inner cap **inside**
  `action.stabilizer().group_generators()` (exit 124) after **40m1s**, peak tree 3675 MiB of 6144: the stabilizer
  itself did not complete. That is checked rather than inferred: the raw log is scanned for the completion line
  `action.stabilizer() <seconds>: <n> generators`, and **no such line exists** there — the last stage line is the
  call's own announcement. The run manager reports **one descendant cleaned up**; `peakrss` records the tree's peak
  and wait status, not cleanup, so that fact is attributed to the manager.
- What that measures: the in-line dedupe is **cost-negative at N scale** — it did not finish where the verbatim
  family did (1898.5 s). Its termination frames land in lattice-morphism tensor construction
  (`_tensor_view -> tensor.from_morphism -> from_matrix -> matrix_entry`). That is the attributed stack at
  termination and **not** an isolated hashing cost: no per-component split against image evaluation was measured, so
  no such split is claimed.
- Corrected at the owner (the corrections this section makes are carried in provider `77da8a3`, which withdrew the
  unsupported timing and causal estimates): provider `71f5c13` makes `lattice_stabilizer` return the verbatim Schreier
  family again (whose measured completion on N is 1898.5 s, from `jmw0urc9n`), and keeps the dedupe
  as `_deduplicated_generating_morphism`, with its contract stated - image family and subgroup
  preserved, not the free domain or word enumeration - and its unhashable boundary stated. The raw
  log of this rerun is `.tmp/enriques_stageD_dedupefix_6c1500f.log`, published with this note.
- What `jmw0urc9n` establishes about the conversion stage, and no more: the probe finished
  `action.stabilizer()` at t=1919.0 s and hit the unchanged 2400 s cap ~481 s later inside its
  per-element `N.integral_isometry(g)` loop (frames `integral_restriction -> factor_through_or_none
  -> _preimage_or_none`); no processed conversion count and no per-conversion rate was logged, so the
  previously quoted "~512 s implied" total is withdrawn. Likewise the parenthetical attributing the
  `jmw0whgrg` cap to coordinate hashing is beyond the evidence: its terminal frames land at
  `_tensor_view -> tensor.from_morphism -> from_matrix -> matrix_entry`, the attributed stack at
  termination, not a measured cause split against image evaluation; that attribution is unmeasured.
  The per-image integral-conversion stage, not the stabilizer construction, remains the acceptance
  frontier; the same capped run is not repeated before owner work on it.
