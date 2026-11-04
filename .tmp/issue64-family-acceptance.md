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

**Not attributed, and deliberately so.** An earlier sweep of these same node ids at provider `42d7361c`
recorded peaks of 3735–3846 MiB. This run, at port `7b6102e7` + provider `cdafca84` under the gate
runner, reached 6186 MiB before the same case finished. The two runs differ in provider revision and in
harness, so the difference is **measured but not attributed** — it is neither a provider regression nor
a port regression on this evidence, and calling it either would be exactly the inference the standing
rules forbid.

What it does mean for `#64`: `vector_stabilizers` at this pair cannot be accepted from this repository
under the 6144 MiB guard as the file stands, and the single E10 case that was sampled does not fit
either. Establishing the family needs either the CI six-shard arrangement (which partitions work but
does not by itself lower a single case's peak) or a case-level memory investigation that starts by
reproducing the `42d7361c` figure under this exact runner.

## What is *not* claimed

Nothing here says a family is complete because its file exists, because its markers are absent, or
because a sample passed. `#64`'s obligation is the remaining families at one pair; three of the nine now
have that and one is measured as needing sharding. The Enriques-dependent parts of `boundaries` are
blocked on provider items recorded in `.tmp/handoff-owner-enriques-coset-overflow.md`, which is a
provider matter and not a port defect.
