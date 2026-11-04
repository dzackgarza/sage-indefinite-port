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

## What is *not* claimed

Nothing here says a family is complete because its file exists, because its markers are absent, or
because a sample passed. `#64`'s obligation is the remaining families at one pair; three of the nine now
have that and one is measured as needing sharding. The Enriques-dependent parts of `boundaries` are
blocked on provider items recorded in `.tmp/handoff-owner-enriques-coset-overflow.md`, which is a
provider matter and not a port defect.
