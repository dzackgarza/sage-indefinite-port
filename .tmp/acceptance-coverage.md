# #33 acceptance coverage, tied to the pair each result actually ran on

Every row names the port checkout and provider revision that produced it. `captured` means the run
header echoed both SHAs; `derived` means they were reconstructed from commit timestamps and tree
hashes, with the reasoning recorded in `.tmp/HANDOFF.md`. No row is attributed to a revision it did not
run on, and no cap is counted as a pass or a failure.

## A_PERFECT

| what | port | provider | result | provenance |
| --- | --- | --- | --- | --- |
| full file, one invocation per test: all 92 | `87fbf863` (src tree) | **`cdafca84`** | **92/92 pass**, zero failures; max peak 506 MiB; slowest `33` isotropic 77 s | captured in the run header + bind-verified |
| earlier full file: 80/80 parametrized + 12/12 standalone | `a848fd7c` | `e83c4d6` | 92/92 pass; sum 642 s, slowest 65 s | derived (log had no SHA) |
| formerly capped six, rerun at 1800 s | `a848fd7c` | `e83c4d6` | 33 total 36.17 s, 33 iso 56.30 s, 37 total 40.12 s, 37 iso 34.51 s, 39 total 34.41 s, 39 iso 23.52 s | derived |

This is **local** evidence. A_PERFECT requires a hosted run at an accepted head SHA; none has been
observed, so the row stays open.

## A_STABILIZERS

| node | port | provider | result |
| --- | --- | --- | --- |
| OSCAR `stabilizer_of_100` | `7a2692f` | `e83c4d6` | passed 65 s (derived) |
| `r_-1` | `6e1be797` | `e83c4d6` | passed 1850.35 s (captured) |
| `r_0` | `fc6a022f` | `9dc89ac` | passed 696.33 s (captured, bind-verified) |
| `r_1` | `92991467` | `9dc89ac` | passed 402.82 s (captured) |
| `r_2` | `43a49697` | `9dc89ac` | passed 200.44 s (captured) |
| `r_3` | `64cf8c6f` | `9dc89ac` | passed 129.70 s (captured) |
| `r_4` | `ac1f1a8f` | `9dc89ac` | passed 109.04 s (captured) |
| `r_5` | `1a2ead82` | `9dc89ac` | passed 6.07 s (captured) |
| `r_6` | `685b26fd` | `9dc89ac` | passed 80.75 s (captured) |
| `r_7` | `4feccb6f` | `9dc89ac` | passed 488.51 s (captured) |
| `r_8` | `16bfcdc6` | `9dc89ac` | passed 203.62 s (captured) |
| square `[0]` | src `87fbf863` / `e18e04ee`-era | `9dc89ac` | **killed, no verdict** — exit 137 at 373 s against a 7200 s cap |
| squares `[2] [4] [6] [8] [10]` | — | — | unresolved, not failed |
| **square `[0]`** | src `87fbf863` | **`565fc6b9`** | **passed 160 s, peak 3 893 MiB** — unchanged assertions, 6 144 MiB guard, cap 7 200 s |
| **square `[2]`** | src `87fbf863` | **`565fc6b9`** | **passed 176 s, peak 3 773 MiB** |
| **square `[4]`** | src `87fbf863` | **`565fc6b9`** | **passed 222 s, peak 3 795 MiB** |
| **square `[6]`** | src `87fbf863` | **`565fc6b9`** | **passed 385 s, peak 3 283 MiB** |
| **square `[8]`** | src `87fbf863` | **`565fc6b9`** | **passed 256 s, peak 3 720 MiB** |
| **square `[10]`** | src `87fbf863` | **`565fc6b9`** | **passed 238 s, peak 3 834 MiB** |

Eleven accepted cases at the older pins, none a cap, spanning 6.07 s to 1850.35 s. Four (`r_-1`, `r_0`,
`r_1`, `r_7`) would have been lost to the 600 s cap used earlier.

**The six square nodes are now all accepted** on the maintained provider `565fc6b9`, which carries the
composite-disjunct-first reorder in `IntegralStructureAction.stabilizer()` (plus the two earlier
`containing is None`-first reorders). `[0]` previously never finished: 413 s kill at 8 127 MiB tree, and a
1 200 s run with nine stack dumps in one dead `is_finite` query. Now the slowest square is 385 s and the
largest peak is 3 893 MiB, so the guard is never approached. Bind-verified in-run:
`rational_integral_stabilizers.py` sha256 `4b0bff018f332d3e`, all three repairs present.

## A_ISOMETRY

| what | port | provider | result |
| --- | --- | --- | --- |
| `test_change_of_basis_gives_an_isometric_lattice`, all six | `ab140b70` | `e83c4d6` | passed: 3 s / 113 s; 151 s / 122 s; 336 s / 575 s |
| **full file, one case per invocation: all 166** | src tree `87fbf863` (= `HEAD:src`) | **`c546077d`** | **166/166 pass**, zero non-zero exits; max peak 3 845 MiB; slowest `ci_indef_U_2U_2E8` 577 s / 3 758 MiB (the historical `MemoryError`); total wall 4 501 s |
| remaining 160 cases in the file | — | — | **done — see the 166/166 row above; this row is superseded** |

The 166-case row is read off the per-node lines of the sweep log (job `jmvxurihq`), every one of which ends
`verdict=1 passed`; the header records `port 75dc17f6-equivalent src`, and `git rev-parse 75dc17f6:src` is
`87fbf863c413…`, identical to `HEAD:src`, with zero commits touching `src/` between them — so that run is the
current source tree, not an equivalent of it.

The six were failing at HEAD before `ab140b70` (rational-coordinate `int()` crossing), which was fixed
in the port, not the provider. The earlier hosted 158-case baseline does not certify current code.

## A_REFLECTIVE

| what | port | provider | result |
| --- | --- | --- | --- |
| `test_edgewalk_branch_representative_sample` | pre-session | `7b222e86`-era | passed 14.5 s (branch coverage only) |
| `test_binary_cycle_coordinate_order_and_period_certificate`, `…_respects_permuted_coordinate_basis` | `a848fd7c`-era | `e83c4d6` | passed 2.38 s together |

**Contract blocker:** `allcock_rank3_reflective.json` has 8 595 entries and only **176 distinct ids**
(120 ids have non-identical duplicates), so #33's "recorded by case id" cannot be satisfied for that
parametrization. That belongs to #33, not to a local fixture edit.

## A_GROUPS

| node | port | provider | result |
| --- | --- | --- | --- |
| `test_mertens_omega_and_minus_one_generate_o_l` | `a848fd7c`-era | `e83c4d6` | passed 7.11 s / 6.54 s (two pins) |
| nodes 1–6 (mertens + five Hecke/ambiguous binary forms) | src `87fbf863` | `cdafca84` | 6/6 passed, 4–46 s |
| **nodes 10–45 and 47 (37 nodes)** | src `87fbf863` | **`565fc6b9`** | **37/37 passed**, 0 failures; sum 1 780 s; slowest 96 s; max peak 3 753 MiB |
| node 46 `a7_plus_diagonal` | — | — | **excluded, documented open** |
| nodes 7–9 `LORENTZ_OrthogonalGroup{0,1,2}` | src `87fbf863` | `cdafca84` / `565fc6b9` | **open — see below** |

The 37 passing nodes are the degenerate/radical family (10–17) and the OSCAR discriminant-image family
(18–45, 47).

### The three Lorentzian nodes are open, with a measured cause

Generation of the group is **12.9–20.8 s** (2 quotient cells, 3 adjacencies). The rest of a 1 200 s run is the
acceptance assertion's own mod-p image-order comparison: `reduction_image_order` is
`MatrixGroup([matrix(GF(p), rows)…]).order()` (`consumer.py:111`), and `comparison_primes` returns
`{2,3,5,7}` for all three cases (det −1024, −256, −1024). Measured on `LORENTZ_OrthogonalGroup0`:

| p | image order | cost |
| --- | --- | --- |
| 2 | 46 998 591 897 600 = `|O^+_10(2)|` | 0.06 s |
| 3 | 5 158 051 199 765 220 556 800 = `|O^+_10(3)|` | 4.94 s |
| 5 | not obtained | **> 1 000 s**, flat 535 MiB, GAP `Size` on a matrix group over `F_5` |
| 7 | not reached | action space `7^10 = 282 475 249` points |

Both sides agree at p=2 and p=3. Mutual containment was tested as a cheaper criterion and **also fails at
p=5** (900 s cap, 3 862 MiB), so it is not a fix. Restricted-prime runs are diagnostics only and are **not**
counted as all-prime acceptance. Detail: `.tmp/investigation-lorentzian-image-order.md`.

## A_VECTORS (new: first run on current pins)

One invocation per node, src `87fbf863`, provider `cdafca84`, guard 6 144 MiB, cap 1 200 s.

| what | result |
| --- | --- |
| 60 of 64 nodes executed | **58 passed**, max peak 4 341 MiB, sum 5 619 s, slowest 569 s |
| `test_published_representatives_are_one_per_orbit[14]` | **FAILED** (exit 1, 652 s) — open, message not yet captured |
| `test_k3_lattice_has_one_primitive_orbit` | **exceeded its 1 200 s cap** (killed at 1 325 s) — open |
| the four `test_k3_lattice_has_one_primitive_orbit_per_norm` nodes | **not run** — rank-22 K3 family, deferred under host memory pressure (3.3 GiB free vs ~4 GiB nodes) |

## A_CENTRALIZERS, A_VECTORS, A_ISOTROPIC, A_BOUNDARIES

A_CENTRALIZERS: 5/5 accepted at src `87fbf863` + `c546077d` (2.98 s/417 MiB, 2.56 s/274 MiB,
80.18 s/3 739 MiB, 39.46 s/2 786 MiB, 38.99 s/2 801 MiB). A_VECTORS: see its own section above.

### A_STABILIZERS, the pointed-isometry nodes (14)

All 14 nodes of the non-E10 pointed family passed at src `87fbf863` + **`565fc6b9`**: 14/14, sum 146 s,
slowest 68 s, max peak 3 581 MiB (nodes 1–14 of `test_vector_stabilizers.sage`).

- **A_BOUNDARIES**: in flight on `565fc6b9` at the time of this record, one node landed with a **failure** —
  `test_dawes_buildings_have_the_published_counts[dawes_2U_A2_stable]`, exit 1 after 172 s. Open; the
  assertion message has not yet been captured, and the remaining nodes have no verdict. Enriques `case_1`
  remains the documented blocker (1800 s on `7b222e86`, 3600 s on `9dc89ac`, 3600 s on `c546077d`); its
  display blocker was fixed by `9dc89ac`.
- **A_ISOTROPIC**: **not run**. The file collects **8 832** nodes (8 821 of them
  `test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors`, ranks 3–18,
  6 476 with `has_isotropic=False`). The earlier "186 cases" figure was wrong. A systematic stride sample,
  not an unbounded chain, is the intended coverage.

## Not port acceptance

The 80/80 Normaliz-vs-PPL engine agreement is an **engine-only** statement about identical ray sets. It
is not owned-cone cost and not acceptance, and it never depended on the provider.
