# Parked: schedule the `double_cosets` acceptance family in CI

**Why parked.** `.github/workflows/acceptance.yml` carries an explicit `shards` dict of exactly the
nine families it runs. The new `tests/acceptance/test_double_cosets.sage` (19 cases, all passing
through `just ci-acceptance`'s recipe body) is therefore not scheduled by CI until one line is added
to that dict. Editing a workflow file requires the GitHub App `workflows` permission, which this
connection does not have; `#66` records the same boundary for the image pin. Nothing here bypasses it.

**The patch.** One line, after `"centralizers": 1,`:

```yaml
              "double_cosets": 1,
```

**Why no other change is needed.** `just ci-acceptance <family> <shard> <shards> [selection]` already
handles an unlisted family through its `*` branch (`tests=(tests/acceptance/test_$family.sage)` with
`--shard-id/--num-shards`), so the recipe, the summary job and the artifact names all work unchanged.
Only the matrix that *schedules* the family needs the line.

**Verified before parking.** `run-ci-acceptance.sh double_cosets 0 1 ''` on the corrected runtime,
pinned pair `c746b913` + `cdafca84`: **19 passed in 1.33 s**, `junit.xml` = 19 tests, 0 failures,
0 errors, peak tree RSS 397 MiB of 6144.
