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

### The one real gap, and why it is not closable here

`tests/fixtures/double_coset_cases.json` is cited in `oracle_manifest.yaml`
(`extraction: references/extract/polyhedral_common_ci.py`, source
`references/vendor/polyhedral_common@1592b246/CI_tests/DoubleCosets/DBL`), and **no acceptance module
consumes it**:

```
load_double_coset_cases -> 0 acceptance file(s)
```

The disposition table says `test_integral_structures.sage → double-coset cases`; that internal file is
gone and its replacement cases were never written, so one capability currently has no suite coverage —
the exact situation `#33`'s item 4 forbids ("removed once that case exists").

Three facts decide the disposition:

1. Under the ownership amendment the operation is **provider-owned**: this port has no
   `integral_structures` module (correctly — the amendment says that port entry point is not created),
   and the oracle is the provider's double-coset decomposition.
2. A tenth acceptance module **cannot be run by the gate as it stands**: `acceptance.yml` carries an
   explicit `shards` dict of exactly the nine families, and adding one is a workflow edit, which this
   connector cannot make (missing App `workflows` permission; `#66` records the same boundary).
3. The fixture is pure permutation-group data (degree 55, group orders 5760/1440), so it cannot be
   folded into an existing lattice-case family without misdescribing it.

So the gap is precise and blocked, not overlooked: provider-owned operation, oracle cited and unused,
gate matrix not extensible from here.

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

## Why #33 is not being closed

Items 2, 3 and 4 stand on the evidence above, and the corpus half is intact. Item 1 does not: one
operation in the unit's own table has no case, its internal test was removed ahead of its replacement,
and the replacement cannot be added to the gate from this environment because that needs a workflow
matrix edit this connector is not permitted to make.

Closing `#33` would assert a complete proof boundary that is missing one operation. It stays open, with
this note as the record, until either the gate matrix can gain a `double_cosets` family or the research
side's own suite covers the provider-owned double-coset oracle. Algorithm success and downstream family
acceptance are separate and are not claimed here: `#64` owns the family-level verdicts, and original
`5/10/19` and `S1` remain open.
