# #64 family acceptance at one pair: what runs, what does not, and why

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
