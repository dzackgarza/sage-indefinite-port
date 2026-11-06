# #64 family acceptance at one pair: what runs, what does not, and why

> **Superseded in part, 2026-11-06.** The tables below were measured on the pre-2026-11-05 runtime and
> are history. The current per-family figures, measured on the adopted runtime at port `45c22a7c` +
> provider `cdafca84`, are in `.tmp/acceptance-coverage.md` under "Full rerun on the adopted runtime,
> 2026-11-05/06", together with the per-case `orthogonal_group` table and the `boundaries` guard stop.
> Read that section for the state of #64; the text below is kept for its method and its attributions.

Recorded pair for every row: port at the commit named, provider **`cdafca84`** (the `pyproject.toml`
pin), corrected runtime, guard `PEAKRSS_LIMIT_MIB=6144` process-tree. Each accepted row is a run whose
per-case `junit.xml` and log sit in `.tmp/ci/results/<family>-0/`.

| family | cases | result | record |
| --- | --- | --- | --- |
| `centralizers` | 5 | **5 passed**, 78.71 s | `centralizers-0` |
| `perfect_domains` | 92 | **92 passed**, 357.79 s | `perfect_domains-0` |
| `isometry` | 166 | **166 passed**, 1279.75 s | `isometry-0` |
| `double_cosets` | 19 | **19 passed**, 1.33 s | `double_cosets-0` |
| `vector_stabilizers` | 31 | **guard hit at 15/31**: peak tree RSS 6147 MiB > 6144 | `vector_stabilizers-0/log` ends at 15 progress marks, no summary |

## The `vector_stabilizers` result, stated precisely

The whole file cannot be run as one shard at this pair. `peakrss` recorded the tree at **6147 MiB
(hwm 6159)** against the 6144 MiB guard at t = 752 s and the run stopped there; the pytest log ends
after fifteen progress marks with no summary line, so **no case failed** — the run was stopped by the
guard, and fifteen cases are all that had completed.

Two things follow, and neither is a verdict on the family:

* the CI matrix runs this family as **six shards** (`acceptance.yml`, `"vector_stabilizers": 6`), which
  is why whole-file shard 0/1 is not the CI arrangement. My local runner omitted `--shard-id/--num-shards`
  because `pytest-shard` is CI-provisioned and absent here; that deviation was harmless for the four
  families above but is **not** harmless here, where the file exceeds the guard.
* the honest way to accept this family from this repository is bounded subsets — individual cases or
  small groups — each finishing under the guard, exactly as the standing "one heavy node at a time" rule
  requires. A guard stop is not a case verdict and not a cost measurement.

## Bounded subsets of `vector_stabilizers` at this pair

Following the sharding conclusion, the cases after the guard stop were run one process each:

* cases 1–15 of 31 completed in the whole-file run before the guard stopped it (no per-case verdict
  lines, because `-rA` prints its summary at the end and the run never reached it);
* the six E10 square cases (`test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[0,2,4,6,8,10]`)
  were then run individually, each with its own 380 s cap. The first one already took the tree to
  **6186 MiB (hwm 6758)** at t = 214 s, over the 6144 MiB guard, and the run stopped there; `pgrep`
  shows no surviving process.

**Attributed, by holding everything but the provider fixed.** The earlier sweep's figure came from a
different harness, so it proved nothing on its own. The same case was therefore re-run at the *same port
commit*, with the *same command*, under the *same runner and guard*, changing only the provider tree:

| provider | result |
| --- | --- |
| `cdafca84` (this repo's pin) | **guard stop**: tree 6186 MiB (hwm 6758) at t = 214 s |
| `42d7361c` | **1 passed in 134.99 s**, peak tree RSS **3702 MiB**, low-water 6841 MiB |

Provider lineage, from `git merge-base --is-ancestor` on the provider repository:

```
cdafca84 (2026-11-01, our pin) -> 42d7361c (2026-11-02) -> ba4efaa (2026-11-04)
```

So the pin this repository is on is **an ancestor** of the revision at which the case passes, and the
difference is provider-attributable: ~2.5 GiB of process-tree memory on this one case. That is not a
port regression (the port commit is identical in both runs) and it is not a harness artefact (the
harness is identical in both runs).

What it does mean for `#64`: `vector_stabilizers` at this pin cannot be accepted from this repository
under the 6144 MiB guard — not because the file needs sharding (CI's six shards partition work but do not
lower a single case's peak), but because **one case needs ~2.5 GiB more at this pin than at a later
provider revision**. The two provider-side items now measured are therefore, together:

1. `ba4efaa` (Nov 4) fixes the coset transversal index overflow that fails `case_87`/`case_50` of the
   Enriques family — absent at the pin (0 occurrences of `_gap_list_entry`/`ELM_LIST`), present there (14);
2. a revision between the pin and `42d7361c` (Nov 2) brings the E10 stabilizer case under the guard —
   6186 MiB at the pin against 3702 MiB there, same port, same harness.

Pin adoption is the owner's decision; both measurements are recorded so it can be made on evidence.

## `vector_stabilizers` bounded run at `42d7361c`: 21 of 31, zero failures

With the attribution above established, the whole file was run once at provider `42d7361c` (port
unchanged, `-v` so each verdict prints as it lands, inner cap 3200 s):

```
passed = 21   failed = 0
```

The 21 are the fifteen base cases plus **all six E10 squares**
(`test_e10_stabilizer_indices_match_brandhorst_gonzalez_alonso[0,2,4,6,8,10]`) — the very cases that
stopped the pin run on the guard. The ten `test_equivalence_witness_carries_a_representative_to_its_reflection`
cases had not completed when the bounded run's cap was reached; `r_-1` alone has taken ~1271 s in an
earlier measurement, so this is a horizon question, not a failure.

Final figures of that run, from its own files: `junit.xml` records **21 tests, 0 failures, 0 errors**;
the log's last verdict is `...alonso[10] PASSED [ 67%]`, so the cap was reached entering the `r_*` block;
`peakrss` records **peak tree RSS 5139 MiB** of the 6144 MiB guard with low-water MemAvail 3265 MiB.
The 5139 MiB is the accumulated high-water of 21 cases in one process, well above the 3702 MiB of the
single E10 case measured alone — which is what a shard boundary is for, and why the CI arrangement runs
this family as six shards rather than one process.

Per-case lines and the in-flight `junit.xml` are in `.tmp/ci/results/vector_stabilizers-42d7361c/`.
This is **bounded branch coverage at a non-pin provider**, and it is labelled that way on purpose: the
pin remains `cdafca84`, where the family cannot be accepted at all because of the memory difference. It
is not the family's acceptance and does not become one by being quoted.

## What is *not* claimed

Nothing here says a family is complete because its file exists, because its markers are absent, or
because a sample passed. `#64`'s obligation is the remaining families at one pair; three of the nine now
have that and one is measured as needing sharding. The Enriques-dependent parts of `boundaries` are
blocked on provider items recorded in `.tmp/handoff-owner-enriques-coset-overflow.md`, which is a
provider matter and not a port defect.

## Interrupted at the runtime checkpoint (NOT acceptance)

Two runs of the `boundaries` family were stopped for runtime maintenance. Neither is evidence of
passing, and neither is reported as such.

**1. `dawes`/`dtower`/`k3` subset, 25 selected cases — stopped by me after 30 min on its first case.**
The first selected case, `test_dawes_building_incidence_graph_is_the_published_figure[dawes_2U2_A2_stable]`,
had not finished in ~30 minutes; it builds the same Dawes building whose counts case has capped
historically. Verdict: **interrupted, zero cases completed**. A case that has not finished may not be
called expensive; it is simply unfinished.

**2. `dtower` + `k3` subset, 23 selected cases — interrupted by the runtime checkpoint at 2 passed.**

```
test_dtower_boundary_component_counts[N_3] PASSED [  4%]
test_dtower_boundary_component_counts[N_4] PASSED [  8%]
```

Two cases completed, zero failed, and the run was killed with the rest unrun. That is **2 of 23 at the
pin with a per-case record**, not 23, and not the family. No `junit.xml` was written (pytest writes it at
the end), so the record here is the per-case `-v` lines and the log itself, which is retained.

Both runs are preserved under `.tmp/ci/results/boundaries-subset-0/` and
`.tmp/ci/results/boundaries-fast-0/` with their `commit` and `research_commit` files.

## Runtime checkpoint: quiescence evidence

At the stop instruction the rootfs was checked before and after termination:

```
my processes in the rootfs (test_boundaries / peakrss.py / run.sh, this project) : 3 before, NONE after
any rootfs user, all projects (proot-x86_64 / sage-image/rootfs)               : NONE
```

Three PIDs (`1294461`, `1294469`, `1294470`) were sent `TERM`, then re-checked; `NONE` remained. No new
rootfs job was launched, and none will be until a verified resume. Pins, guard (`6144 MiB`), and
ownership are unchanged: `pyproject.toml` still pins `cdafca84`, the provisioned provider checkout is
untouched, and no dependency repository was edited.

## State of #64 at 2026-11-06, port `45c22a7c` + provider `cdafca84`

Measured on the adopted runtime (`/sage` at `9f1680d2`), 6144 MiB tree guard, one family per process.
Full detail and the per-case tables are in `.tmp/acceptance-coverage.md`.

| family | cases run / in criterion | result | what #64 still needs |
| --- | --- | --- | --- |
| double_cosets | 19/19 full file | **19 passed** | nothing; now also scheduled by CI (plan row added in `45c22a7c`) |
| centralizers | 5/5 full file | **5 passed** | nothing |
| perfect_domains | 92/92 full file | **92 passed** | nothing |
| isometry | 166/166 full file | **166 passed** | nothing |
| reflective | 5/5 sample | **5 passed** | nothing for the amended node |
| orthogonal_group | 9 sample cases, one process each | **7 passed; 2 with no verdict at the 1800 s cap** | `LORENTZ_OrthogonalGroup0` and `a7_plus_diagonal` need a horizon decision (both under the guard) |
| isotropic_orbits | 2 of 13 at the 2400 s cap | 2 passed, 0 failed | 11 sample cases unrun |
| vector_orbits | 7 of 30 at the 2400 s cap | 7 passed, 0 failed | 23 sample cases unrun |
| vector_stabilizers | 14 of 25 at the 3400 s cap | 14 passed, 0 failed | 11 sample cases unrun (`alonso[10]`, `r_-1..r_8`) |
| boundaries | 1 passed then guard stop at 6146 MiB, no JUnit | **no verdict** | the Enriques provider item, and the Dawes/classical-index subgroup primitive |

Nothing in this table is an acceptance of a family that did not finish, and no cap or guard stop is
counted as a pass or as a failure. The hosted requirement is unchanged: `ci-acceptance` in the
Acceptance workflow is the only source of hosted evidence, and no hosted run at this head has been
observed from this workspace (the connector exposes no Actions runs).
