# Measurements taken for direction d6 (exact certificate layer)

Environment: this VM, `python3` 3.12.14 with the SageMath install at `/opt/math`
(`import sage.all` works in plain `python3`), 4 cores, 15 GB. No A7 specimen was constructed and
the 49 861-row `SupportHyperplanes` call was **not** re-run.

Reproduce: `cd .ultramath/smvzd6hvn/m7/code && python3 bench.py > ../results/bench.json`.
Raw output: `results/bench.json`; primitives: `code/primitives.py`; harness: `code/bench.py`.

## 1. Correctness of (P3) normal-from-incidence and (P1) facet certificate

Ground truth: Sage's own exact facet enumeration of the same cone (and, for the cube cones, the
closed-form facet family `x_0 ± x_j ≥ 0`).

| cone | rays | dim | facets | normal recovered exactly from incidence alone | (P1) certificate accepted |
| --- | --- | --- | --- | --- | --- |
| cone over 9-cube | 512 | 10 | 18 | 18 | 18 |
| cone over 3-cube | 8 | 4 | 6 | 6 | 6 |
| random cone over 20 points | 20 | 4 | 18 | 18 | 18 |
| random cone over 30 points | 30 | 5 | 92 | 92 | 92 |
| random cone over 40 points | 40 | 6 | 365 | 365 | 365 |

499 facets, all recovered **from the incidence set alone** (`O(n^4)` cofactor construction, no
dependence on `|I|`), all certified, and all matching the engine's normals. The dim-10 case is the
specimen's dimension and has 256-ray incidence sets, i.e. non-simplicial facets of the size that
makes a rational nullspace over `|I|` rows expensive.

## 2. Correctness of the rank certificate

`rank_lower_certificate` vs `matrix(QQ, rows).rank()` on 120 random integer matrices (dim 2–10,
entries ≤ 15, dependent rows included): **0 mismatches**. Each returned rank is backed by an exact
Bareiss determinant of a square minor, not by the modular reduction.

## 3. Correctness of (P2) adjacency

Cone over the 9-cube: 18 facets, 153 pairs, analytic answer 144 adjacent (the 9 non-adjacent pairs
are the opposite facets `x_0 - x_j`, `x_0 + x_j`). (P2) returns 144; it read at most 65 rows of any
joint incidence set. Cone over the 3-cube: 15 pairs, analytic 12, (P2) returns 12.

**A dead end found by this harness and kept on the record:** the "no ray is negative on both
normals" test — which an earlier draft used as the decision procedure — passes all 153 pairs on
the 9-cube cone, i.e. it wrongly calls all 9 opposite-facet pairs adjacent. The condition is
*necessary* (if a ray were negative on both, `y + y'` would leave the dual cone) but **vacuous**,
because facet normals are nonnegative on every ray. (P2) must be the joint-rank test.

## 4. Cost at the specimen's ray count (49 861 rays, dimension 10)

| operation | measured | note |
| --- | --- | --- |
| one pairing pass of a covector over all 49 861 rays (`numpy int64`, certified bound) | **29 ms** | `matrix(ZZ, …) * vector` in Sage: 119 ms |
| one incidence / zero-support bitset per facet orbit | **15 ms** | 6.1 KB per orbit |
| bitset disjointness test (sound pre-filter, Lemma 5 — prunes nothing) | 0.44 µs | kept only as a cheap sanity query |
| exact (P2) adjacency decision, per facet-orbit pair | **181–208 µs** | includes the bitset AND and the certified minor |
| exact (P3) normal from an incidence set | 0.05–0.3 ms | independent of `|I|` |
| certified rank lower bound, 256 rows | 0.35 ms | Sage `QQ` rank: 0.45 ms |
| certified rank lower bound, 1 024 rows | 0.88 ms | Sage `QQ` rank: 1.87 ms |
| certified rank lower bound, 4 096 rows | 2.99 ms | Sage `QQ` rank: 7.32 ms |
| certified rank lower bound, 16 384 rows | 11.5 ms | Sage `QQ` rank: 31.2 ms |

The certified rank is *faster* than `matrix(QQ, …).rank()` from 256 rows upwards because it exits
as soon as `n` pivots are found: its cost is `O(n^3)`-ish in the found pivot count rather than
`O(|I| n^2)` rational arithmetic.

## 5. Derived closure cost (the cost model of §5 of `solution.md`)

From the per-pair 181–208 µs and 15 ms per orbit support:

| facet orbits `o` | support builds | all `o(o-1)/2` exact adjacency decisions | state size |
| --- | --- | --- | --- |
| 50 | 0.8 s | 0.22 s (1 225 pairs) | ≤ 0.3 MB |
| 200 | 3.0 s | 3.6 s (19 900 pairs) | ≤ 1.2 MB |
| 600 | 10.4 s | 37.5 s (179 700 pairs) | ≤ 3.7 MB |

Compare the blocked Normaliz route: peak 5 405 MiB and **CAPPED at 1 800 s** on one
`SupportHyperplanes` call. The certificate layer is two to three orders of magnitude under the cap
even at 600 facet orbits.

## 6. What these numbers do and do not support

* They support: the three primitives are exact and correct against an independent engine at the
  specimen's dimension; a facet-orbit-sized state is a few MB, never a facet list; the per-orbit
  and per-pair costs are small enough that the layer cannot be the bottleneck.
* They do **not** support: any claim about the A7 cell's facet-orbit count `o`, about the
  enumeration's total time, or about the ridge subproblems (facet enumeration in dimension 9 over
  `|I_F|` generators). The costs in §4–5 use synthetic covectors with the specimen's ray count and
  dimension, not A7 covectors.
* No A7 run was made, as required. No acceptance assertion, cap or memory guard was touched; the
  provider library was not modified (the code here is a new module in this member's directory).
