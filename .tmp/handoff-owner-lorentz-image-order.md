# LORENTZ 0/1/2: capped, not wrong — and where the cost actually is

Provider `565fc6b9`, port src `965041fd` (`fd63f677`), Sage 11.0.beta5. No cap raised, no assertion
weakened, no restricted-prime result counted as all-prime acceptance.

## The contract

`test_lorentzian_generators_match_polyhedral_common` (`tests/acceptance/test_orthogonal_group.sage:72`):

```python
computed = lattice(case["gram"]).O().group_generators()
assert_same_finite_images(case["gram"], computed, case["generators"])
```

and `assert_same_finite_images` (`tests/acceptance/consumer.py:120`) compares, per prime in
`comparison_primes(gram)`:

```python
reduction_image_order(generator_rows, prime)
  = int(MatrixGroup([matrix(GF(prime), rows) for rows in generator_rows]).order())
```

`comparison_primes` is `prime_divisors(2 det G) | {3, 5, 7}`. All three LORENTZ cases therefore compare at
`{2, 3, 5, 7}` (determinants −1024, −256, −1024).

## Distinction: the node is capped, not a demonstrated wrong value

Every prime that the harness can actually finish **agrees on both sides**, for all three cases — this is now a
regression, `tests/test_lorentzian_reduction_image_orders.py` (3 tests, 87 s, peak 600 MiB):

| case | p=2 computed | p=2 recorded | p=3 computed | p=3 recorded |
| --- | --- | --- | --- | --- |
| `LORENTZ_OrthogonalGroup0` | 46 998 591 897 600 | same | 5 158 051 199 765 220 556 800 | same |
| `LORENTZ_OrthogonalGroup1` | 45 660 871 065 600 | same | 5 158 051 199 765 220 556 800 | same |
| `LORENTZ_OrthogonalGroup2` | 46 998 591 897 600 | same | 5 158 051 199 765 220 556 800 | same |

At p=5 the run does not finish; p=7 is never reached. **That is a cap, and nothing in the record shows a
disagreeing value at any prime.** It should not be reported as a failure of the port's generating set.

### A real distinction inside that table: `LORENTZ_OrthogonalGroup1` at p = 2

Its determinant is −256, so `L/2L` is **degenerate** and the image is *not* a full orthogonal group: the
measured order is `2^30 · 3^5 · 5^2 · 7 = 45 660 871 065 600`, which has no factor 17 or 31, whereas
`|O^+_10(2)| = 2^21 · 3^5 · 5^2 · 7 · 17 · 31`. Where `p ∤ det` the order *is* the full even orthogonal group:

```
|O^+_10(2)| = 46 998 591 897 600            = measured p=2 for cases 0 and 2
|O^+_10(3)| = 5 158 051 199 765 220 556 800  = measured p=3 for all three cases
```

That identification is asserted conditionally (only where the reduction is nondegenerate), because it is a
statement about the *space*, and the space degenerates mod 2 for case 1.

## Where the cost is: every route, measured

Lattice of rank 10, signature (1, 9), computed generating set of 16–17 isometries.

| route | p=2 | p=3 | p=5 | p=7 |
| --- | --- | --- | --- | --- |
| `MatrixGroup(...).order()` (the harness) | 0.03 s | 5 s | **> 1 000 s**, flat 535 MiB | not reached |
| mutual containment, `all(g in other ...)` | — | — | **> 880 s**, 3 862 MiB, killed | — |
| `libgap.IsSubgroup(GO, G)` | — | — | **5 660 MiB, killed at 1 200 s** | — |
| `libgap.Size(libgap.GO(+1, 10, q))` | 0.00 s | 0.00 s | **0.00 s** | **0.00 s** |
| provider `orthogonal_group_base_change` to `F_p` | rejected | rejected | rejected | — |
| GAP recognition `IsomorphismTypeOfMatrixGroup` | — | — | **not bound in this build** | — |

Two rejections of the provider route, measured rather than assumed:

```
OwnedZZ.hom(PrimeFields()(5))  ->  ValueError: the identity morphism exists only on Mor(R, R), but this is
                                   Mor(Integer Ring, Integer Ring / Ideal (5) of Integer Ring)
ZZ.hom(GF(5))                  ->  TypeError: expected a ring of this session, but Finite Field of size 5
                                   is not one
```

So the maintained semantic exists in principle — `O(L) → O(L ⊗ F_p)` is a functor and the target is a finite
group — but the provider's `base_change` does not currently accept the quotient ring map `Z → F_p` that this
route needs.

`IsomorphismTypeOfMatrixGroup` (GAP's recognition, which would name the group and fix its order from the
classical formula) is **not bound**: only `NiceMonomorphism` and `IsomorphismPermGroup` are, and both build the
action on `F_p^10` — 9 765 625 points at p=5 and 282 475 249 at p=7.

## What the p=5 and p=7 answers must be, if the images are full

The classical order is available instantly from GAP's own classical group object, so it is not in doubt:

```
|O^+_10(5)| = 54 461 311 079 175 000 000 000 000 000 000        (libgap.Size(libgap.GO(1,10,5)), 0.00 s)
|O^+_10(7)| = 209 544 577 891 300 558 570 289 055 128 616 960 000  (0.00 s)
```

Since the p=2 and p=3 images are the full orthogonal groups (cases 0 and 2), the p=5 and p=7 images are
expected to be the full orthogonal groups too, and those two numbers are what the assertion would compare.
**That expectation is not a measurement and is not counted as acceptance.**

## What the owner needs to supply, precisely

Either

1. **recognition**: a build of GAP with `IsomorphismTypeOfMatrixGroup` bound (the Recog package), which
   identifies the image as `O+(10, p)` in bounded time and yields its order from the classical formula; or
2. **a provider-level structural certificate**: an owned operation stating "this generating set spans the
   full isometry group of the space `L/pL`", which the acceptance harness could then use in place of two
   independently computed orders — the same shape as the existing `discriminant_image_order`, which is cheap
   (1.81 s) precisely because the discriminant group is small; or
3. acceptance on the affordable primes only, stated explicitly — which is a decision about the criterion, not
   something to slip into the harness.

The harness itself needs no change to remain correct; it needs a cheaper way to compute the same two numbers.
