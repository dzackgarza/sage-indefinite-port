# #64 shard decisions from measured per-case cost

Pair: port `893cb233` + provider `cdafca84`, adopted runtime (`/sage` at `281da26ffd`). Guard 6144 MiB on
the **process tree**, inner cap per shard 345 min in CI.

## The rule these decisions follow

One matrix row is one pytest process, and the guard bounds that process's whole tree — Sage plus the
Julia bridge, whose resident baseline alone is ~3.0-3.3 GiB once a case touches Oscar. So a sample cannot
be partitioned by case count; it is partitioned by whichever quantity binds it, and the two are different
in different families:

- **memory-bound**: the *sum* of the cases' retained peaks reaches the guard. Partition so that no shard
  holds two heavy cases. Splitting such a family by count alone leaves a shard with two 3.7 GiB cases and
  it dies at 7.4 GiB.
- **time-bound**: the peak stays well under the guard and the *sum of durations* exceeds the cap. There,
  round-robin over an array ordered by cost balances the shards, and memory is not a constraint.

A guard stop is not a case verdict and not a cap: it says the arrangement, not the case, exceeded a bound.

## Per-family evidence and decision

### `orthogonal_group` — memory-bound; 9 sample cases → 6 shards, implemented

Per case, each run alone in its own process (peak tree RSS, and JUnit duration):

| case | peak | duration |
| --- | --- | --- |
| `[LORENTZ_OrthogonalGroup0]` | **5490 MiB** | exit 124 at the 1800 s cap |
| `test_degenerate_orthogonal_group_embeds_non_diagonal_complement` | 3089 MiB | 81.0 s |
| `test_degenerate_orthogonal_group_generates_rank_two_radical_action` | 2989 MiB | 77.7 s |
| `test_degenerate_orthogonal_group_contains_radical_shear_and_gl_action` | 2959 MiB | 87.1 s |
| `test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q` | 1858 MiB | exit 124 at the 1800 s cap |
| `test_mertens_omega_and_minus_one_generate_o_l` | 381 MiB | 9.9 s |
| `test_signed_discriminant_image_of_u_plus_minus_three_is_all_of_o_q` | 368 MiB | 4.2 s |
| `test_zero_form_orthogonal_group_is_the_full_general_linear_group` | 326 MiB | 1.6 s |
| `test_orthogonal_group_of_the_hyperbolic_plane_has_order_four` | 291 MiB | 0.9 s |

The undivided sample reached **6236 MiB at t = 618 s** with five cases resident (2026-11-06, this pair; the
same arrangement reached 6154 MiB on the previous runtime). Summing two 3 GiB degenerate cases already
exceeds the guard, so the partition is by case, not by count:

| shard | cases | expected peak |
| --- | --- | --- |
| `orthogonal_group_0` | `[LORENTZ_OrthogonalGroup0]` | 5490 MiB |
| `orthogonal_group_1` | `a7_plus_diagonal` | 1858 MiB |
| `orthogonal_group_2` | `degenerate…contains_radical_shear_and_gl_action` | 2959 MiB |
| `orthogonal_group_3` | `degenerate…generates_rank_two_radical_action` | 2989 MiB |
| `orthogonal_group_4` | `degenerate…embeds_non_diagonal_complement` | 3089 MiB |
| `orthogonal_group_5` | the four light cases | ≤ 1500 MiB |

Implemented in `justfile` as explicit selection names, and each collects exactly its cases (verified:
1+1+1+1+1+4 = 9, a partition of the sample, no duplicates).

### `vector_stabilizers` — memory-bound; 25 sample cases → 18 shards, implemented

Seventeen of the 25 peak above 3.6 GiB **alone** (roots `r_0..r_8` 3673-4708 MiB, e10 squares `[0,2,4,6,8,10]`
3735-3846 MiB, `divisibility_four…exact_containing_bound` 3651 MiB, `vector_stabilizer_has_oscars_order`
3608 MiB) and eight peak under 440 MiB. Measured per case in jobs `jmvz58wdp`, `jmvzyjrg9`, `jmvz3a0di`
(the peaks are runtime-stable across the two runtimes we have: the 14-case shard peaked at 4405 MiB here and
5139 MiB for 21 cases earlier, consistent with retention, not with release).

The family's own run is **time**-stoppable as well — 14 of 25 cases in 3400 s — but a partition chosen for
time would put two heavy cases in one shard and lose it to the guard. So: the selection array is ordered by
measured peak, descending, and the matrix runs **18 shards**. Round-robin over that order gives each shard
one heavy case plus at most one light one; the worst shard is `r_-1` (4708) + `divisibility_two` (297) =
**5005 MiB**, under the guard. The ordering comment in the `justfile` says so, because sorting that array
alphabetically would silently reintroduce the failure.

### `vector_orbits` — time-bound; 30 sample cases → 3 shards, implemented

The seven cases that completed in the 2400 s cap took 510, 471, 466, 463, 458, 0.3 and 0.2 s, and the
family's peak over those seven was **1281 MiB** — the guard is not close. The binding quantity is time, so
three round-robin shards are safe and give ~800 s each on the measured cases. Twenty-three sample cases are
still unmeasured at this pair; the shard count is a floor for them, not a claim about them.

### `isotropic_orbits` — undecided, one shard, with the measurement named

Two cases complete in 2.2 s and 0.0 s, and the third (`test_o_plus_n_has_two_line_and_two_plane_orbits`)
runs past 2400 s with the tree at 4844 MiB for the pair. I have no per-case peak for the long ones, so the
partition is **not** decided here: the next bounded measurement is the peak of that single case alone, which
says whether it is a memory partition or a horizon question. Left at one shard deliberately rather than
guessed.

### `boundaries` — undecided, one shard, with the measurement named

The sample passed one case and was killed by the guard at **6146 MiB**, so it is memory-bound, but the
per-case peaks of its 18 cases are unmeasured and two of its rows are provider-blocked anyway (Enriques
needs a provider revision past the pin; Dawes and the classical indices share the subgroup primitive). Not
partitioned by guess.

## What the matrix carries now

`push` grows from 24 to **48** rows; `workflow_dispatch` stays at 74 (full files, unchanged shard counts).
`reflective` keeps one shard (5 cases, 2654 MiB). Every sampled selection now accepts `--shard-id/--num-shards`,
which is what lets the matrix partition a sample at all; with `shards = 1` it is the identity, so the
already-passing arrangements are unchanged.

The `justfile` half is committed and pushed. The matrix half is a **workflow file** (`.github/workflows/acceptance.yml`)
and is committed locally as a workflow commit for the orchestrator to land.
