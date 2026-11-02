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
| **`r_-1`** | src `65ea68a68af0` | **`42d7361c`** | **passed 1 271 s, peak 4 708 MiB** |
| **`r_0`** | src `65ea68a68af0` | **`42d7361c`** | **passed 498 s, peak 3 988 MiB** |
| **`r_1`** | src `65ea68a68af0` | **`42d7361c`** | **passed 352 s, peak 4 121 MiB** |
| **`r_2`** | src `65ea68a68af0` | **`42d7361c`** | **passed 186 s, peak 3 839 MiB** |
| **`r_3`** | src `65ea68a68af0` | **`42d7361c`** | **passed 136 s, peak 3 714 MiB** |
| **`r_4`** | src `65ea68a68af0` | **`42d7361c`** | **passed 126 s, peak 3 830 MiB** |
| **`r_5`** | src `65ea68a68af0` | **`42d7361c`** | **passed 10 s, peak 318 MiB** |
| **`r_6`** | src `65ea68a68af0` | **`42d7361c`** | **passed 106 s, peak 3 673 MiB** |
| **`r_7`** | src `65ea68a68af0` | **`42d7361c`** | **passed 454 s, peak 4 060 MiB** |
| **`r_8`** | src `65ea68a68af0` | **`42d7361c`** | **passed 254 s, peak 3 780 MiB** |
| **square `[0]`** | src `49e3965f` tree | **`42d7361c`** | **passed 152 s, peak 3 782 MiB** |
| **square `[2]`** | src `49e3965f` tree | **`42d7361c`** | **passed 116 s, peak 3 735 MiB** |
| **square `[4]`** | src `49e3965f` tree | **`42d7361c`** | **passed 208 s, peak 3 807 MiB** |
| **square `[6]`** | src `49e3965f` tree | **`42d7361c`** | **passed 182 s, peak 3 786 MiB** |
| **square `[8]`** | src `49e3965f` tree | **`42d7361c`** | **passed 198 s, peak 3 804 MiB** |
| **square `[10]`** | src `49e3965f` tree | **`42d7361c`** | **passed 242 s, peak 3 846 MiB** |
| **OSCAR `stabilizer_of_100` (node 15)** | src `49e3965f` tree | **`42d7361c`** | **passed 70 s, peak 3 608 MiB** |

Eleven accepted cases at the older pins, none a cap, spanning 6.07 s to 1850.35 s. Four (`r_-1`, `r_0`,
`r_1`, `r_7`) would have been lost to the 600 s cap used earlier.

The roots family was re-baselined on the current provider `42d7361c` (src `65ea68a68af0`) because that
pin changes what `_group_cardinality` asks. **All ten root nodes now pass at the current pins** (sweep
`jmvyzjrg9`, src `65ea68a68af0`, provider `42d7361c`, cap 3 600 s, 6 144 MiB guard): sum 3 293 s, slowest
`r_-1` at 1 271 s, max peak 4 708 MiB. Same assertions; no node reached either limit.

The comparison to the older rows is a **re-baseline, not a controlled one** — every earlier figure also
differs in the port src, which is not held fixed anywhere. Read strictly: two nodes dropped a lot
(`r_-1` 1 850.35 s at port `6e1be797`/provider `e83c4d6` → 1 271 s; `r_0` 696.33 s at `fc6a022f`/`9dc89ac`
→ 498 s), and the seven comparable-to-`9dc89ac` nodes are **mixed within noise**: `r_2` 200 → 186, `r_7`
489 → 454, but `r_3` 130 → 136, `r_4` 109 → 126, `r_5` 6 → 10, `r_6` 81 → 106, `r_8` 204 → 254. So the
provider pin is not a systematic speedup; the two large drops are the two heaviest nodes and the rest move
both ways. What the sweep establishes is acceptance at the current pins, and that the family is bounded by
`r_-1`'s memory, not by a cap.

**The E10 squares and the OSCAR node are now re-baselined at the same pin.** All seven pass, one pytest
invocation per node, cap 1 800 s, guard 6 144 MiB, unchanged assertions (job `jmvz3a0di`): sum **1 168 s**
against **1 502 s** for the same seven at `565fc6b9`. Node by node, old → new: `[0]` 160 → 152, `[2]` 176 →
116, `[4]` 222 → 208, `[6]` 385 → **182**, `[8]` 256 → 198, `[10]` 238 → **242**, OSCAR 65 → **70**. The large
drop is `[6]`; two nodes ticked up slightly. So the pin is not a uniform speedup here either, but nothing
regressed materially and nothing came near a limit.

That leaves the 14 pointed-isometry nodes (1–14) as the only part of A_STABILIZERS still carrying a
`565fc6b9` row; they cost 146 s together, so re-baselining them is a short run, not a heavy one.

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
| node 46 `a7_plus_diagonal` | — | — | **open — blocked on cone conversion** of a 49 861-vector perfect configuration; see the dependency section below |
| nodes 7–9 `LORENTZ_OrthogonalGroup{0,1,2}` | src `87fbf863` | `cdafca84` / `565fc6b9` | **open — see below** |

The 37 passing nodes are the degenerate/radical family (10–17) and the OSCAR discriminant-image family
(18–45, 47).

### The three Lorentzian nodes are open, with a measured cause

**They are capped, not wrong.** Every prime the harness can finish agrees on both sides, for all three cases
(`tests/test_lorentzian_reduction_image_orders.py`, 3 tests, 87 s):

| case | p=2 | p=3 |
| --- | --- | --- |
| `LORENTZ_OrthogonalGroup0` | 46 998 591 897 600 both sides | 5 158 051 199 765 220 556 800 both sides |
| `LORENTZ_OrthogonalGroup1` | 45 660 871 065 600 both sides | 5 158 051 199 765 220 556 800 both sides |
| `LORENTZ_OrthogonalGroup2` | 46 998 591 897 600 both sides | 5 158 051 199 765 220 556 800 both sides |

At p=5 the run does not finish; p=7 is never reached. `LORENTZ_OrthogonalGroup1` has determinant −256, so its
p=2 reduction is degenerate and its image is *not* a full orthogonal group there (`2^30·3^5·5^2·7`); where
`p ∤ det` the image is exactly `|O^+_10(p)|`.

Generation is 12–21 s. The measured cost of every route to the p=5/p=7 *orders*:

| route | p=5 | p=7 |
| --- | --- | --- |
| `MatrixGroup(...).order()` (the harness) | > 1 000 s, flat 535 MiB | not reached |
| mutual containment (`g in H`) | > 880 s, 3 862 MiB | — |
| `libgap.IsSubgroup(GO, G)` | 5 660 MiB, killed at 1 200 s | — |
| `libgap.Size(libgap.GO(+1, 10, q))` — the *ambient* group only | **0.00 s** | **0.00 s** |
| GAP recognition (`recog`, `genss`) | **packages not installed** in GAP 4.14.0 | — |

The provider route does **not** stop at the ring map: `PrimeFields()(p).quotient_map()` is the canonical
retained `Z -> F_p` map, and with it `orthogonal_group_base_change` → `base_change(...).Aut()` →
`subgroup(images)` completes in ≤ 0.2 s at p=2 *and* p=5. What is missing is one decision: the subgroup's
`is_finite()` is already `True`, but `cardinality()` falls through to `_group_cardinality`
(`groups.py:383-394`), which counts only groups *placed* in `OwnedFiniteGroups()`, and no
`_selected_cardinality_decision` is installed. Exact contract and the two things the owner must supply:
`.tmp/handoff-owner-lorentz-image-order.md`.

Nothing here concludes anything about the image from the ambient order or from p=2/p=3 agreement.

### Provider pin moved to `42d7361c` (admission only)

The maintained library published `42d7361cfda3d8505f365b26d9aa687447d6df9d` — "a decided finiteness is
placed, and cardinality asks it" — adopted through the maintained checkout and verified against the actual
remote (it is that repository's `origin/main`; staged by `git archive` of the commit). It closes the
admission gap recorded here. Verified without re-running the engine BSGS:

```
[p=2] is_finite True 0.10s | in OwnedFiniteGroups now True (was False)
[p=2] cardinality 46998591897600 in 0.02s == reduction_image_order at p=2
[p=5] is_finite True 0.11s | in OwnedFiniteGroups now True (was False)
[p=5] cardinality NOT called
```

So the provider's own semantics now reproduce the harness number at p=2. At p=5 the fix routes
`cardinality()` to `_finite_order` → `as_permutation_group()`, i.e. the same BSGS; **the p=5/p=7 obligation is
unchanged and those three nodes stay open.** The six accepted square nodes were accepted at `565fc6b9` and
have not been re-run at `42d7361c`; their rows keep the pin they ran on.

## A_VECTORS (new: first run on current pins)

One invocation per node, src `87fbf863`, provider `cdafca84`, guard 6 144 MiB, cap 1 200 s.

### The K3 per-norm family is now **accepted, 4/4**, at the original limits

src = the `49e3965f` tree (memo + immutable shared Gram), provider `42d7361c`, one pytest invocation per node,
**cap 1 800 s and guard 6 144 MiB unchanged**, original assertions, each node with its own transcript and
60 s stack dumps:

| node | result | peak |
| --- | --- | --- |
| `test_k3_lattice_has_one_primitive_orbit_per_norm[-2]` | **passed 462 s** | 3 948 MiB |
| `[2]` | **passed 466 s** | 3 960 MiB |
| `[4]` | **passed 448 s** | 3 938 MiB |
| `[-4]` | **passed 430 s** | 3 957 MiB |

Sum 1 806 s, no node within 1 200 s of its cap. Before the memo two of these capped at 1 800 s (3 579 /
3 564 MiB) and the other two had never been started, so this is a family that moved from "not established" to
accepted on a **port fix**, with no cap raised and no assertion weakened. The stack dumps put the bulk of each
node inside `right_coset_transversal()` — the finite-representation coset space, not the vector dedup.

**What this does and does not say.** It accepts one parametrised test — the four K3 per-norm nodes — and
nothing more. It is not A_VECTORS acceptance: the other five test functions in
`test_vector_orbits.sage` (`primitive_orbit_counts_in_u_plus_e8`, `published_representatives_are_one_per_orbit`,
`e10_orbit_counts`, `i_2_10_is_transitive`, `anisotropic_binary_vector_classes`) are separate obligations with
their own rows, and the `[14]` failure record there is untouched by this. It is not port acceptance either:
A_BOUNDARIES, the A_ISOTROPIC tail, A_GROUPS LORENTZ and `a7_plus_diagonal` all remain open below.

### Unit suite at the current pins, over all of `tests/`

src = the `38d1c9f1` tree, provider `42d7361c`, guard 6 144 MiB:

```
198 passed, 2 warnings in 140.26 s          peak tree RSS 3 995 MiB
```

This run covers every file in `tests/`, including `test_acceptance_ci.py`, which the earlier 182-test figure
did not: running it exposed one **stale test** (`test_batch_runner_continues_after_malformed_batch_junit`
still passed `--sage`, removed from `ci/pytest_batches.py` in `0e544edd`), repaired in `38d1c9f1` by changing
the invocation and no assertion. So the immutability change is validated at suite scale, and the suite is now
green over its whole directory rather than over ten named files.

| what | result |
| --- | --- |
| 60 of 64 nodes executed | **58 passed**, max peak 4 341 MiB, sum 5 619 s, slowest 569 s |
| `test_published_representatives_are_one_per_orbit[14]` | **FAILED** at that pin pair (exit 1, 652 s) — does **not** reproduce at current pins (see below) |
| `test_k3_lattice_has_one_primitive_orbit` (no `_per_norm` suffix) | **unreliable earlier record, not an obligation.** No function of that name exists in the file at any commit, and the log line carrying it also carries hand-written prose, so it is a reconstruction rather than a collected node id. Nothing is claimed either way about it. |
| the four `test_k3_lattice_has_one_primitive_orbit_per_norm` nodes | **not run in that sweep** — deferred under host memory pressure; **all four now accepted**, see the K3 section above |

### `[14]` does not reproduce at current pins

Run again on provider `565fc6b9` with the port at src `965041fd` (`fd63f677`), unchanged assertions:

| run | result |
| --- | --- |
| first | **passed 153.0 s**, peak 4 012 MiB |
| `PYTHONHASHSEED=0` | **passed 193.5 s** |
| `PYTHONHASHSEED=1` | **passed 162.6 s** |
| `PYTHONHASHSEED=2` | **passed 197.1 s** |

4/4 passes across four hash seeds, so it is also stable there — not a flake at the current pins. The
old-pin failure stands as an observation at provider `cdafca84` with src `87fbf863` (exit 1 after 652 s);
its cause is **not established**, and that pin pair is no longer current. The recorded A_VECTORS result
above stays at its own pins and is not rewritten.

## A_CENTRALIZERS, A_VECTORS, A_ISOTROPIC, A_BOUNDARIES

A_CENTRALIZERS: 5/5 accepted at src `87fbf863` + `c546077d` (2.98 s/417 MiB, 2.56 s/274 MiB,
80.18 s/3 739 MiB, 39.46 s/2 786 MiB, 38.99 s/2 801 MiB). A_VECTORS: see its own section above.

### A_STABILIZERS, the pointed-isometry nodes (14)

All 14 nodes of the non-E10 pointed family passed at src `87fbf863` + **`565fc6b9`**: 14/14, sum 146 s,
slowest 68 s, max peak 3 581 MiB (nodes 1–14 of `test_vector_stabilizers.sage`).

- **A_BOUNDARIES**: partially run on `565fc6b9`, then paused by decision (not by failure) after four nodes,
  to give the blocked research operation and the observed failures priority over extending the sweep.
  Every Dawes building node that ran **failed**, all with one cause, now repaired (see below).

### A_BOUNDARIES: the Dawes failures, their cause, and the repair

Four nodes ran before the pause:

| node | before repair | after repair |
| --- | --- | --- |
| `test_dawes_buildings_have_the_published_counts[dawes_2U_A2_stable]` | FAILED exit 1, 172 s | **passed 274 s, peak 3 941 MiB** |
| `[dawes_2U_minus2_minus6]` | FAILED exit 1, 214 s | ran 1 672 s and **exceeded the 6 144 MiB guard** (peak 6 147 MiB, exit 143) — open |
| `[dawes_L2_gamma2]` | FAILED exit 1, 262 s | runs; **exceeded the 1 800 s cap** (2 279 s, peak 6 094 MiB) — open |
| `[dawes_2U2_A2_stable]` | **CAPPED** at 1 800 s (1 858 s, peak 3 903 MiB) | not re-run |

Not a wrong count, and not a Dawes-specific problem: `IsotropicInduction._integral_image`
(`isotropic_induction.py:151`) converted a coordinate with `int()`, and that coordinate belongs to the
*rational* space, so it is a base-ring element:

```
TypeError: int() argument must be a string, a bytes-like object or a real number,
not 'OwnedJoinCategory.element_class'
```

This is the same defect class as the old-pin T_Co traceback, and the port already had the right helper for
it — `edgewalk_scalars._integer`, which was introduced by `ab140b70` for `recursive.py` and never applied
here. **Repair:** that one conversion now reads through `_integer`, which keeps the integrality invariant
(a non-integral entry raises `ArithmeticError` instead of being truncated). No behaviour changes on the
integral path.

Direct evidence that the repair targets this defect, produced on the same object:

```
coordinate container : sage_extensions.owned_category.OwnedJoinCategory.element_class
pre-fix  int(entry)  -> TypeError: int() argument must be ... not 'OwnedJoinCategory.element_class'
fixed    _integer(entry) -> 1
```

**Regression:** `tests/test_isotropic_integral_image.py` — 2 tests, 1.20 s: an identity of the rational space
must carry an integral vector to itself through `_integral_image` (this raised `TypeError` before), and a
genuinely non-integral image must raise `ArithmeticError` rather than be truncated.

**Provenance note.** The sweep driver prints `git rev-parse HEAD:src` when it starts. The Dawes and T_Co runs
above executed the *working tree*, whose content is commit `06bc9562` with src tree
`965041fd64e6cad6d0a7f4cebff69ef28f06b0f5` — i.e. `87fbf863` plus this one repair. A header line reading
`87fbf863` for those runs is the pre-commit `HEAD:src`, not the tree under test.

## The research T_Co handoff: does not reproduce on current pins

Specimen reproduced verbatim from `dzackgarza/research@a145cea4`
`computations/reports/coble_tco_surjectivity_specimen.py`, including its construction
(`Lattices(ZZ)([list(row) for row in GRAM.rows()])`, not a Sage matrix) and both same-session controls.
Provider `565fc6b9`, port `44dc4d6c` (src `87fbf863c413`), Sage 11.0.beta5, guard 6 144 MiB.

```
construction: rank 11, determinant -2048, signature (2, 9), even True   <- matches the artifact exactly
call: L.discriminant_representation_is_surjective()
  returned True
controls: S_Co -> True,  M = S_En + <delta> -> True
RESULT: returned True
```

29 min, peak tree 4 239 MiB (python 1 161, julia 3 306), host low-water 9 071 MiB.

**The old-pin `TypeError` does not persist.** It was raised at the engine's `recursive.py:154`
`int(coordinates(label))`, and that site has since been reading through `_integer` (`ab140b70`), so no repair
is required for this operation and none was made. Old-pin record: provider `c546077d`, Sage 11.0.beta5.
- **A_ISOTROPIC**: covered by a **documented systematic sample**, not the 8 832-node chain. The file collects
  8 832 nodes: 8 821 `test_isotropic_line_orbits_exist_exactly_when_polyhedral_common_finds_isotropic_vectors`
  (ranks 3–18; 6 476 with `has_isotropic=False`), plus 11 others. The earlier "186 cases" figure was wrong.

  | coverage | result |
  | --- | --- |
  | every 100th node (`stride 100`, 89 cases across the whole index range) | **89/89 passed**, 0 failures, max peak 749 MiB, slowest 548 s, sum 968 s |
  | the 11 non-parametrised nodes | 1 `CAPPED` (`test_o_plus_n_has_two_line_and_two_plane_orbits`, 1 800 s cap, peak 4 604 MiB); the other 10 not run |

  So the isotropic-line family is healthy everywhere it was sampled, and the heavy tail is the
  Enriques-adjacent group-orbit set, which is where the known blocker lives.

## What still depends on what (read from source, not assumed)

An inventory of the remaining obligations and what actually blocks each:

| obligation | blocker | evidence |
| --- | --- | --- |
| A_VECTORS `test_k3_lattice_has_one_primitive_orbit_per_norm[-2,2,4,-4]` | **accepted 4/4** — the cost was a port-side recomputation defect, not a provider gap | 462 / 466 / 448 / 430 s at peaks 3 948–3 960 MiB, cap 1 800 s, guard 6 144 MiB |
| A_VECTORS `test_k3_lattice_has_one_primitive_orbit` (no `_per_norm`) | **not an obligation** — no such node exists in the file at any commit | see the note in the A_VECTORS table |
| A_STABILIZERS roots (10 nodes) | **resolved** — all ten accepted at `42d7361c`, 3 293 s serial, no cap, max peak 4 708 MiB | sweep `jmvyzjrg9` |
| A_STABILIZERS E10 squares (6) + OSCAR node 15 | **resolved** — all seven accepted at `42d7361c`, 1 168 s serial | job `jmvz3a0di` |
| A_STABILIZERS pointed-isometry nodes 1–14 | **cheap re-baseline pending** — accepted at `565fc6b9` only, 146 s together | ~300 s serial with startup |
| A_BOUNDARIES (87 Enriques + ~28 non-Enriques) | **heavy slot**; Enriques `case_1` is a known block at 3 600 s | only 4 Dawes nodes ever ran |
| A_ISOTROPIC tail (10 non-parametrised + 1 capped) | **heavy slot**; the capped one is the Enriques-adjacent set | 89/89 stride sample passed |
| A_GROUPS LORENTZ 0/1/2 | **external owner** — GAP lacks `recog`/`genss`, and no recognition route exists for the generated matrix group | `.tmp/handoff-owner-lorentz-image-order.md` |
| A_GROUPS node 46 `a7_plus_diagonal` | **cone conversion**, not a cap — see below | `.tmp/HANDOFF.md` |
| BOUND/#18 finite-representation transporter | **unverified since the provider moved**, not known-failing | see below |
| FINITE, LIFTS | **research #401** | TODO |

### BOUND/#18: the note predates `42d7361c` and is now unverified

The TODO says the provider's finite-representation transporter "still fails when it requests cardinality of a
represented finite torsion module". That was written before the provider's cardinality-admission fix. A probe
instrumenting `IntegralStructureAction.transporter` on the divisibility-four isotropic `U(4)+U(4)` specimen
finds a witness **without the transporter branch being entered at all** (the rational transport is already
integral), so that specimen cannot settle it either way. The claim therefore stays *unverified* rather than
either resolved or refuted, and the instrumentation to settle it exists
(`.tmp/probe-transporter-path.py`): it needs a pointed pair whose rational transport is non-integral, which
the A_STABILIZERS pointed family may supply when that family is re-run at the current pin.

### K3: measured attribution, and a withdrawn inference

**Call graph** (`_pointed_classes`, `marked_lattices.py:150-170`): K3 = `U^3 + E8 + E8` has positive index 3,
no rule above `case 1` applies, so the case goes to `RECURSION.primitive_vector_classes`
(`recursive.py:1057-1065`), which enumerates `OrbitCoverModel.primitive_covering_vectors(square)` and
deduplicates by pairwise `pointed_isometry`.

That is the graph. The **measured attribution**, one stage at a time, each stage printed before it starts:

| stage | measured |
| --- | --- |
| imports + Gram + `lattice(gram)` | 1.3 s |
| `lattice.signature_pair()` → (3, 19) | **0.0 s** |
| `AttackProfile.from_lattice` → positive_index 3, sign 1 | **0.0 s** |
| `EichlerRecursion._model(signed_view)` | 9.4 s |
| **`primitive_covering_vectors(-2)`, first 200 vectors** | **> 410 s, fewer than 200 yielded** |

So the profile is *not* the cost; the cover enumeration's **yield rate** is. Memory stayed at the Julia
baseline (3 580 MiB tree), so this is time, not memory.

**Withdrawn inference.** I had written that the other three norms "are no smaller" shells, and therefore add
no information. That was an inference from norm magnitude, not a measurement, and it does not follow: the
enumeration could behave differently per square. `[-2]` and `[2]` were run and both capped (3 579 and
3 564 MiB); `[4]` and `[-4]` were **not** started, and no cap was raised.

**Resolved: `[-2]` passes, and the "yield rate" reading was itself incomplete.** The stage table above was
taken *before* the Gram memo, and the missing 410 s was not the cover's yield rate. `_isometry_from_row_matrix`
asks for the source and the target Gram on every generator, and `stabilizer_of_lattice` builds **136**
generators, so 272 calls to a 25 s pure function of an immutable rank-22 rational lattice — about 6 800 s —
sat in front of the first cover vector. Memoizing `_gram_in_label_order` (commit `829793bf`, hardened to
return an immutable shared matrix in `49e3965f`) removes that, and the node now **passes**:

| node | limit | result | where the time went |
| --- | --- | --- | --- |
| `test_k3_lattice_has_one_primitive_orbit_per_norm[-2]` | 1 800 s cap, 6 144 MiB guard | **passed 462 s, peak 3 948 MiB** | assembly done by 60 s; seven consecutive 60 s stack dumps inside `envelope.right_coset_transversal()` (`eichler.py:489`), reached from `primitive_covering_vectors` (`:1000`) |

### The recomputation defect class: searched, one instance, now exhausted

The K3 fix was one instance of a general fault shape — a **pure function of an immutable lattice recomputed
per call**. I checked every other hot function of that shape in the port rather than fixing them one cap at a
time:

| function | memoized? | where |
| --- | --- | --- |
| `groupoid_cache._gram_in_label_order` | **was not** — fixed | the K3 defect |
| `invariants._gram` | yes, `_GRAMS` by identity holding the lattice | invariants |
| `indefinite_reduced_presentation` | yes, `@cache` | 9 call sites, all share one reduced presentation per lattice |
| `perfect_complex` | yes, "memoized as IndefiniteForms.g:1173-1299 memoizes its groups" | `lorentzian_cells.py:415` |
| `OrbitCoverModel.source_generator_matrices` | yes, `_GENERATOR_MATRICES` by model identity | `eichler.py:661` |
| `Envelope.stabilizer_of_lattice`, `right_coset_transversal` | yes, `@cached_method` | `eichler.py:461,479` |

So there is no second instance of the class among the paths these families take, and the remaining caps are
not this fault. In particular the A_GROUPS node 46 cone blocker is *not* a recomputed perfect complex: the
complex is already memoized and its cost is the configuration's, not a repeat.

So the K3 finite-representation locus reached from the vector side *completes* for this configuration. That is a
positive instance for the BOUND/#18 row below, which until now had only capped runs and a stale defect note:
the same `right_coset_transversal` computation does terminate here, inside 1 800 s, and it is the bulk of the
time. It does not by itself accept BOUND — different specimen, different statement — but it moves that row
from "no successful instance" to "one successful instance on the K3 envelope".

### K3 and BOUND share one mathematical locus

The K3 covering does not fall back to a naive model enumeration. Measured on the signed view:

```
model built: envelope is None? False        # the Eichler envelope exists
similarity.cokernel().exponent() = 1        # in 0.11 s -- so c runs over divisors of 1 only
```

`primitive_covering_vectors` (`eichler.py:980-1000`) therefore takes its **envelope branch**, whose first step is
`envelope.right_coset_transversal()`. That call goes `eichler.py:488` →
`stabilizer_of_lattice()` → `right_cosets()` → `IntegralStructureAction.lattice_stabilizer()`
(`rational_integral_stabilizers.py:228`) → `permutation_representation().preimage_subgroup(...)`, whose
construction is `_FiniteIntegralRepresentation` and its `M/dM` cardinality assertion (`:112-125`).

So the A_VECTORS K3 nodes and BOUND/#18's "represented finite torsion module" obligation are the **same
mathematical locus**, not two unrelated problems: the cover cannot be produced until the lattice stabilizer's
finite representation exists for that lattice. The successful `U + <-2>` transporter specimen exercises a
different action and does not close this one.

**Correct but costly, not missing.** Instrumenting the same call on the K3 specimen, with each step printed
before it runs:

| step | measured |
| --- | --- |
| model + envelope | 10.1 s |
| `OrbitCoverModel.source_generator_matrices(envelope.model)` | **79 s → 136 matrices** |
| assembling those 136 as owned rank-22 rational isometries | **fewer than 20 in the next ~140 s** (≥ 7 s each ⇒ ≥ 950 s for all 136) |
| `rational_envelope.O().subgroup(...)`, then `right_cosets()` (the residue orbit) | not reached |

Nothing raised: no missing operation, no incorrect semantic — the construction is correct and its cost is
`136 approximate generators x owned-isometry assembly at rank 22`, paid *before* the residue orbit that the
finite representation actually needs. So the K3 acceptance nodes are infeasible by construction cost at this
rank, not blocked by a defect, and no port-side mathematical repair is indicated by this evidence.

What that leaves as the real question is a different one, and it is the provider's: whether an action's
lattice stabilizer can be reached without materializing every approximate generator as an owned rational
isometry. That is a semantic-reuse question about the retained contracts, not a port algorithm — the port
already follows the pinned source's generator list.

**Harness lesson, recorded because it cost 900 s.** The first cover probe reported *zero output* and I read that
as "the pre-cover stage exceeded 900 s". It was not: the command piped through `head`, so the probe's output was
truncated by SIGPIPE — the same pitfall already recorded once in this project. The staged probe above, without
the pipe, shows the real attribution.

### `a7_plus_diagonal`: blocked on cone conversion, with the artifact

The node needs `O(L).discriminant_image()` for `A7 + diag(1,1,-1)`. `.tmp/HANDOFF.md:163-172` records the
concrete block: its initial perfect configuration has **49 861 distinct vectors**; the streamed pairing/hash
passes reach exact cone conversion, but neither the full PPL consumer (`.tmp/repair-orthogonal-a7-flint.log`,
provider `07b24a2`) nor the separate Normaliz probe (`.tmp/repair_a7_normaliz.py`,
`.tmp/repair-a7-normaliz.log`) produced a result. So "excluded" is not a placeholder: the blocker is the cone
conversion of a ~50k-vector configuration in two independent engines.

### The transporter: specimen found, branch reachable, and it works at `42d7361c`

A backlog note recorded that the provider's finite-representation transporter "still fails when it requests
cardinality of a represented finite torsion module", written before the provider's cardinality-admission fix.
The first specimen I tried (divisibility-four isotropic `U(4)+U(4)`) never enters the branch — it finds its
witness on the integral path — so it could not settle the note.

A scan of small lattices for a pair whose extended rational transport is **non-integral** found one:

```
L = U + <-2>,   source (-2,-2,-1)  ->  target (2,2,-1),   q = 6 both, content 1 both
```

Instrumenting `IntegralStructureAction.transporter` shows the correction branch is entered and returns a
witness: events `[('ok', True)]`, no exception, and the witness carries the source onto the target. So the
transporter is **reachable and working** at `42d7361c`, and the blanket "still fails" is not true for this
kind of action.

What this does *not* settle: whether the *torsion-module* variant — an action whose finite module `M/dM` is
genuinely torsion rather than free — still hits a cardinality problem. That needs a specimen whose
`finite_module()` is torsion, which this one is not.

Pinned as `tests/test_finite_representation_transporter.py`: it asserts the specimen is genuinely
rational-nonintegral (the transporter must actually be called, else the test fails as an invalid specimen),
that the witness exists and acts correctly, and that every transporter call returned rather than raised.

### A green unit suite is not a proof of no source gaps

182 unit tests passing says the non-acceptance code paths that have tests still behave. It says nothing about
untested branches — for instance `primitive_vector_classes` at positive index >= 3, which is exactly the K3
case above and has no unit coverage.

### Shared GAP capability: a handoff to coordinate, not a stop

The missing `recog`/`genss` packages are recorded as a concrete request in
`.tmp/handoff-owner-lorentz-image-order.md` (exact contract, measured routes, what a recognition build would
unblock). That is an item to coordinate with the provider owner, not a reason to consider LORENTZ closed or
impossible.

### Keyboard-level: nothing else is inventable

Everything above is either waiting on the running job, on the heavy slot that job occupies, or on an external
owner. The non-acceptance unit suite is green (182 passed), so there is no latent unit-level defect to work
while the slot is busy.

## Unit-suite failures: two real test defects, repaired

Running the non-acceptance unit suite on provider `42d7361c` first gave `2 failed, 141 passed`. Both
failures were **defects in the tests themselves**, not in the port and not environmental:

```
tests/test_edgewalk_rank2.py::test_edgewalk_period_action_matrix_uses_coordinate_label_order
  TypeError: base_ring (=Integer Ring) must be a ring or a semiring
tests/test_pinned_lorentzian_cells.py::test_positive_vector_shell_respects_explicit_coordinate_labels
  TypeError: unable to coerce <class 'dict'> to an integer
```

Both tests do `from sage_extensions.all import ZZ`, which shadows Sage's `ZZ`, and then hand that provider ring
object to a **Sage** constructor — `matrix(ZZ, ...)` and `vector(ZZ, ...)`. The provider ring is a different
object, so Sage refuses it. Repairing the second site exposed the reciprocal half of the same mixup: once the
matrix is a Sage matrix its entries are Sage integers, and the provider's `linear_combination` rejects those
("construct it from a ring of this session") while accepting a Python `int` — which is what every other
coefficient in the port and tests uses, and what `consumer.element` does with `int(entry)`.

I had first attributed these to the environment after reproducing them at `3c2425a1` in a detached worktree.
That attribution was wrong: an older commit failing identically shows the defects are old, not that they are
environmental. After the repair the whole non-acceptance suite is green:

```
182 passed in 223.6 s   (peak tree RSS 3 963 MiB)
```

covering `test_edgewalk_rank2`, `test_lorentzian_flipping`, `test_perfect_cell_transporter`,
`test_pinned_lorentzian_cells`, `test_configuration_pairings`, `test_cell_faces_facet_orbits`,
`test_allcock_fixture_identity`, plus the three targeted files of this unit.

## Not port acceptance

The 80/80 Normaliz-vs-PPL engine agreement is an **engine-only** statement about identical ray sets. It
is not owned-cone cost and not acceptance, and it never depended on the provider.
