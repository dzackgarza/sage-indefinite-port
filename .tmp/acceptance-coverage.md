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
| square `[0]` | `e18e04ee` | `9dc89ac` | **killed, no verdict** — exit 137 at 373 s against a 7200 s cap |
| squares `[2] [4] [6] [8] [10]` | — | — | unresolved, not failed |

Eleven accepted cases, none a cap, spanning 6.07 s to 1850.35 s. Four (`r_-1`, `r_0`, `r_1`, `r_7`)
would have been lost to the 600 s cap used earlier.

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
| rest of the file | — | — | not run; one bounded file invocation passed six nodes then exceeded 900 s |

## A_CENTRALIZERS, A_VECTORS, A_ISOTROPIC, A_BOUNDARIES

Historical baseline only (run 40807218493) except as noted below; nothing on current pins.

- `A_BOUNDARIES`: Enriques `case_1` **unresolved** — capped at 1800 s on `7b222e86`, 3600 s on
  `9dc89ac`, rerunning on `c546077d`. Its display blocker was fixed by `9dc89ac`; the surviving cost is
  the orbit-hash path, handed to the provider.
- `test_isotropic_orbits`: 186 cases passed at `7b222e86`-era; not re-run.

## Not port acceptance

The 80/80 Normaliz-vs-PPL engine agreement is an **engine-only** statement about identical ray sets. It
is not owned-cone cost and not acceptance, and it never depended on the provider.
