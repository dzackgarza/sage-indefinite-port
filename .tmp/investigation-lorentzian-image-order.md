# `LORENTZ_OrthogonalGroup0`: the 1200 s is the assertion, not the generators

Provider `cdafca84`, port src tree `87fbf863c413` (= `HEAD:src`). Cap **unchanged at 1200 s**, memory guard
6 144 MiB, one node at a time, **no `src/` or `tests/` edits**. The earlier sweep (job `jmvy3hgmo`) recorded
this node as `exit=124 peak=5765 MiB verdict=none 1203s`.

Artifacts: probe `~/sage-image/probe_lorentzian_stages.py`, logs `~/sage-image/lorentzian-primes.log`,
`~/sage-image/membership.log`, `~/sage-image/seeds-dump.log`, `~/sage-image/compare-seeds.log`,
stage totals in `.tmp/lorentzian-stages.json`.

## Headline

The node is **not** blocked in the port's generation path, and **not** in the provider's `is_finite` /
`stabilizer` path that the `square[0]` work was chasing. Generation of the group takes **12–21 s**. The rest
of the 1200 s is spent inside the acceptance assertion's mod-p image-order comparison, in GAP, on matrix
groups of degree 10 over `F_5` and `F_7`.

That suspect could never have applied here: this Gram matrix has signature `(1, 9)`, so `positive_index = 1`
and `orthogonal_group_generators` (`recursive.py:355`) goes to `lorentzian_orthogonal_group_generators`.
Nothing on that path calls `IntegralStructureAction.stabilizer`.

## Stage attribution on the unmodified path

Wrappers were installed in-process on `LorentzianPerfectComplex._expand_traversal_once`,
`LorentzianPerfectLocalBackend.cell_transporter`, `consumer.discriminant_image_order` and
`consumer.reduction_image_order`; the assertions themselves ran unchanged.

| stage | measured |
| --- | --- |
| `indefinite_reduced_presentation` | 0.24 s |
| perfect-complex traversal | 2 calls, 1.54 s (0.88 s + 0.66 s) |
| quotient cells | **2**, sizes 10 and 18 |
| adjacencies | 3, all with a transporter |
| `cell_stabilizer` | 1 call, 0.06 s |
| `cell_transporter` | 4 calls, 1.40–10.4 s total (Julia warm-up varies) |
| **`O(L).group_generators()`** | **12.9 s → 17 generators / 20.8 s → 16 generators** (see determinism, below) |
| `discriminant_image_order` | 44.89 s cold, then 1.81 s → 46 998 591 897 600 |
| `reduction_image_order[p=2]` | 0.06 s and 0.03 s → 46 998 591 897 600 |
| `reduction_image_order[p=3]` | 4.94 s and 4.28 s → 5 158 051 199 765 220 556 800 |
| **`reduction_image_order[p=5]`** | **did not complete in ~1 000 s** |
| `reduction_image_order[p=7]` | never reached |

`reduction_image_order` is `MatrixGroup([matrix(GF(prime), rows) ...]).order()`
(`tests/acceptance/consumer.py:111`). `comparison_primes` returns `{2,3,5,7}`: `det = -1024` so `2·det`
contributes only `2`, and `{3,5,7}` are added unconditionally.

Nine `faulthandler` dumps over the p=5 wait all sit in one place:

```
sage/groups/libgap_mixin.py:240   cardinality            <- GAP Size, on a matrix group over F_5
tests/acceptance/consumer.py:111  reduction_image_order
tests/acceptance/consumer.py:134  assert_same_finite_images
```

Memory on that route stayed **flat at ~535 MiB** throughout p=5. So p=5 costs *time*, not memory.

## Why p=5 and p=7 are the wall

The group GAP is asked about acts on `F_p^10`:

| p | points `p^10` |
| --- | --- |
| 2 | 1 024 |
| 3 | 59 049 |
| 5 | 9 765 625 |
| 7 | 282 475 249 |

p=5 is 165× the p=3 space and p=7 is 4 783× it. That is the whole shape of the cost.

### The orders GAP is computing are closed-form

Both measured values are **exactly** the order of the full even-type orthogonal group, from
`|O^+_{2n}(q)| = 2 q^{n(n-1)} (q^n − 1) ∏_{i=1}^{n−1} (q^{2i} − 1)` with `n = 5`:

```
|O^+_10(2)| = 46 998 591 897 600                    = observed at p=2
|O^+_10(3)| = 5 158 051 199 765 220 556 800          = observed at p=3
```

and both generating sets give the same value at each prime. If the same type holds at p=5 and p=7 — a
hypothesis, not a measurement — the answers are

```
|O^+_10(5)| = 54 461 311 079 175 000 000 000 000 000 000
|O^+_10(7)| = 209 544 577 891 300 558 570 289 055 128 616 960 000
```

## Containment is not the cheap alternative it looks like

The criterion compares two independently computed **orders**, and equal order is only a *necessary*
condition for "the same finite image" — distinct subgroups can share an order. So mutual containment
(`all(g in other for g in group.gens())`, validated on a known-answer fixture before use) was measured as a
cheaper and more faithful test.

- At p=2 and p=3 it is cheap and it **works on this very case**: both sides agree, in both directions.
- At **p=5 it does not finish either**: killed at its 900 s cap, having grown to **3 862 MiB peak**.

So the obstacle is not the choice of query; it is GAP doing natural-module work on `F_5^10` (and, unmeasured
but 29× larger, `F_7^10`). **I am not claiming containment is a fix**, because measurement says it is not.
The tested alternative is recorded here as tried and failed, so it is not retried.

## A separate port defect found while measuring

`O(L).group_generators()` for this lattice is **not reproducible across processes**. Same pins, same input,
only `PYTHONHASHSEED` varied:

| seed | generators | cell sizes | digest |
| --- | --- | --- | --- |
| 0 | **16** | `(10, 18)` | `286bc4fa49246749` |
| 1 | **17** | `(18, 10)` | `ab56302636b9c86f` |
| 12345 | **16** | `(10, 18)` | `0ede609c68467f3a` |
| 99999 | **16** | `(18, 10)` | `07ad815b815cca5d` |

Four distinct generating sets; the number of generators and the order of the two quotient cells both vary.

**The spanned group does not vary**, at least where it can be checked: for seeds 0 and 1 the mod-2 and mod-3
images are mutually contained in both directions (`equal True` at p=2 and p=3,
`~/sage-image/compare-seeds.log`). So this is a non-canonical choice of generating set rather than a wrong
group — but a generating set that changes between runs is a reproducibility problem on its own terms, and
it is port code, so it is ours to fix.

This also explains part of the runtime variance: GAP's matrix-group algorithms are randomised, and they are
fed a different generating set each run.

## What is and is not established

**Established:** generation is seconds; the wall is the assertion's mod-5/mod-7 GAP work; the images agree at
p=2 and p=3 and equal `|O^+_10(2)|` and `|O^+_10(3)|`; the generating set is hash-seed dependent while the
group is stable at the affordable primes.

**Not established, and not claimed:** that the case passes. p=5 has no verdict on either route and p=7 was
never reached. The 5 765 MiB peak the sweep recorded was not reproduced on the order route (flat 535 MiB);
the containment route reached 3 862 MiB, so the high-water mark still needs full attribution before anyone
calls this an OOM.

## Options, for the person rather than for me to choose

1. Accept the case on the primes where it is affordable, stated explicitly, with the two sides agreeing at
   p=2 and p=3 and the p=5/p=7 orders given by the closed form above once the type is confirmed.
2. Certify the images structurally — the evidence is that they are the *full* `O^+_10(p)` — instead of
   comparing two independently computed orders.
3. Leave the node open. What should not happen is raising the cap: at p=7 the action space is 4 783× the
   p=3 one, so a bigger cap buys nothing.

Nothing here was applied to `tests/` or `src/`. The original assertions are intact and were run as written.
