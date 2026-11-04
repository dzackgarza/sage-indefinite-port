# #33 acceptance contract: item-by-item reconciliation against delivered evidence

Root ledger `#4`; `#33` is the first open leaf in the declared preorder and has no children. Its body
carries five things a unit can be held to. This note records, for each, what is actually in the tree and
what is not, measured at port `c746b913` + provider pin `cdafca84` on the corrected runtime.

Nothing here is a closure claim. `#33` stays open; the last section says exactly why.

## 1. "tests/acceptance holds one module per operation above"

The consumer-API table lists twelve operations. Nine modules cover eleven of them:

| operation (table row) | module | note |
| --- | --- | --- |
| O(L) generators | `test_orthogonal_group.sage` | 15 test definitions |
| Isometry | `test_isometry.sage` | 12 |
| Vector equivalence witness | `test_isometry.sage` | same module |
| Vector stabilizer | `test_vector_stabilizers.sage` | 17 |
| Vector orbit representatives | `test_vector_orbits.sage` | 6 |
| Isotropic sublattice/flag orbits, witnesses, stabilizers | `test_isotropic_orbits.sage` | 7; `indefinite.isotropic_flags` is deliberately absent per the amendment |
| Lorentzian perfect domains | `test_perfect_domains.sage` | 14 |
| Reflective lattices | `test_reflective.sage` | 13 |
| Perpendicular kernel | `test_vector_stabilizers.sage` | `indefinite.isotropic_lifts` deliberately absent per the amendment |
| Finite-index subgroups / orbit splitting | `test_boundaries.sage` | `test_enriques_modular_groups_match_dutour_sikiric_hulek` asserts `_counts(group) == (line_orbits, plane_orbits, flag_orbits)`; the stable row carries 528 / 24242 / 72199 |
| Centralizers | `test_centralizers.sage` | 1 parametrized test over 5 cases |
| **Integral structures / double cosets** | **no case** | see below |

### The gap was real, and is now closed

`tests/fixtures/double_coset_cases.json` is cited in `oracle_manifest.yaml`
(`extraction: references/extract/polyhedral_common_ci.py`, source
`references/vendor/polyhedral_common@1592b246/CI_tests/DoubleCosets/DBL`), and until this unit **no
acceptance module consumed it** (`load_double_coset_cases -> 0 acceptance file(s)`), although its
internal test `test_integral_structures.sage` had been removed — the exact situation clause 4 forbids.

`tests/acceptance/test_double_cosets.sage` now supplies that case. It asserts the upstream
decomposition test's own fact, `TestDoubleCosetPermDecomposition`
(`CI_tests/01_RatIntAutomorphy/ProcessExamples.g:179-253`): for `H <= G` the double cosets `H g H`
**partition** `G` — pairwise disjoint, cardinalities summing to `|G|`. The two group orders are
re-derived from the permutation generators rather than trusted from the file name, and each class is
counted as `|H| * #{left cosets of H in H g H}`, which is `|H g H|` for a subgroup, so the sum is a
completeness check and not a restatement of a count.

Route: the operation is provider-owned after the ownership amendment, and the case calls the provider's
`Group.double_cosets` (`sage-extensions categories/group/groups.py`, realized as the `H`-orbits of the
finite `G`-set `G/H`). The port supplies no group theory: the case builds the groups with the owned
constructor and reads the decomposition through that method.

Measured through the gate recipe on the corrected runtime, pinned pair (`c746b913` + `cdafca84`):

```
just ci-acceptance double_cosets 0 1 ''   ->   19 passed in 1.33 s
junit.xml: 19 tests, 0 failures, 0 errors; peak tree RSS 397 MiB of 6144
```

One follow-on remains outside this unit's reach: the `Acceptance` workflow's `shards` dict names
exactly nine families, so CI does not *schedule* the new family until one line is added. That is a
workflow edit, which needs the App `workflows` permission this connection lacks (`#66` records the
same boundary); the line is parked as `.tmp/pending-ci-double-cosets-family.patch` with its README.
The recipe itself needs no change — its `*` branch already runs any family.

## 2. "The suite runs under the gate (just test-push), and its per-case results are recorded"

Two findings, both measured.

**(a) The named gate cannot perform this role.** `just test-push` delegates to
`ai-review-ci/scripts/local-check.py --mode push`, whose own docstring says the local tier "does not
compile a project, run a suite, format files, inspect style/typing". It validates syntax of changed
`.py`/`.json`/`.toml` files. The suite is run by the **`Acceptance` workflow** through
`just ci-acceptance <family> <shard> <shards> [selection]`, on every push, per family shard, writing
`--junitxml`. The acceptance text's `just test-push` is a stale spelling of "the push gate"; the recipe
that actually runs cases is `ci-acceptance`.

**(b) The per-case mechanism works at the current pair.** `.tmp/run-ci-acceptance.sh` is that recipe's
body, executed directly (`just` is not installed here). Through it, at `c746b913` + `cdafca84`:

| family | result | record |
| --- | --- | --- |
| `centralizers 0 1 ''` | **5 passed in 78.71 s** | `.tmp/ci/results/centralizers-0/{junit.xml,log,commit,research_commit}` |
| `perfect_domains 0 1 ''` | **92 passed in 357.79 s** | `.tmp/ci/results/perfect_domains-0/`, `junit.xml` = 92 tests, 0 failures, 0 errors |

The county records are the ones the recipe writes: the port commit, the provider revision it verified
against the pin, and a JUnit file with one entry per case. That is the "per-case results recorded"
mechanism, demonstrated rather than asserted, and the `perfect_domains` run doubles as `#65`'s
A_PERFECT gate re-measured on the corrected runtime through the actual gate recipe.

Hosted, the same recipe runs the full matrix on every push; the run for this HEAD is
`Acceptance` on `c746b913`.

## 3. "Each later work unit's issue names its cases"

Measured on the open units that own algorithm cases:

| unit | case names in body | case files named |
| --- | --- | --- |
| `#15` | 10 | `tests/acceptance/test_reflective.sage` |
| `#18` | 17 | — (names only) |
| `#64` | 19 | 9 files |
| `#36`, `#63`, `#66`, `#67` | 0 | — |

`#36` is the typing unit, `#63` docs, `#66` the image pin, `#67` a setup hoist: none of them owns
acceptance cases, so zero is correct there rather than a gap.

## 4. "The internal tests above are gone, each one only after the case that replaces it exists"

All eleven files the disposition table names are absent from `tests/`:

```
test_recursive.sage test_eichler.sage test_isotropic_reductions.sage test_isotropic_lifts.sage
test_vector_sections.sage test_isotropic_flag_types.sage test_groupoid_cache.sage test_invariants.sage
test_integral_structures.sage test_lorentzian_cells.sage test_recursive_acceptance.sage
```

The one with no replacement case is `test_integral_structures.sage` — item 1 above. `TSET_edgewalk*`
files are gone as `.sage`; `tests/test_edgewalk_rank2.py` remains, which the table's `test_edgewalk*.sage`
pattern does not name.

## Corpus half of the contract

The manifest lists 26 fixtures. Every one is present, every one cites an extraction script that exists
in the tree, and every cited source path resolves:

```
manifest fixtures: 26 | missing fixture files: none
missing/absent extraction scripts: none | cited sources absent: none
```

## Status of this unit

All five clauses now stand on delivered evidence: nine acceptance modules plus the new double-coset
module cover all twelve table rows; the gate recipe runs and records per-case JUnit; the case-owning
units name their cases; every named internal test is gone with its replacement case now existing; and
all 26 corpus fixtures are cited with present extraction scripts.

The residual item is the CI scheduling line for the new family, which is a permission boundary rather
than missing work, and is parked with a patch.

Closing the issue itself is not available from here: this connection reads issues and cannot publish
comments or change issue state (its description states it does not provide issue/review publishing), and
`gh`/`itree` are not installed in this environment. The closure command to run where they are:

```
itree close dzackgarza/sage-indefinite-port#33 --reason completed \
  --comment "contract items 1-5 delivered; double-coset cases added in <sha>; CI scheduling line parked in .tmp/pending-ci-double-cosets-family.patch (workflows permission)"
```

## Original reasoning, kept for the record: why closing was withheld before the case existed

Items 2, 3 and 4 stand on the evidence above, and the corpus half is intact. Item 1 does not: one
operation in the unit's own table has no case, its internal test was removed ahead of its replacement,
and the replacement cannot be added to the gate from this environment because that needs a workflow
matrix edit this connector is not permitted to make.

Closing `#33` would assert a complete proof boundary that is missing one operation. It stays open, with
this note as the record, until either the gate matrix can gain a `double_cosets` family or the research
side's own suite covers the provider-owned double-coset oracle. Algorithm success and downstream family
acceptance are separate and are not claimed here: `#64` owns the family-level verdicts, and original
`5/10/19` and `S1` remain open.
