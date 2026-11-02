# LORENTZ 0/1/2 mod-5/7 image order: the provider route works, the last decision does not exist

Provider `565fc6b9`, port src `965041fd`, Sage 11.0.beta5, GAP 4.14.0. No cap raised, no assertion weakened,
no restricted-prime result counted as acceptance, and nothing here infers image equality from the ambient
group's order or from p=2/p=3 agreement.

## Correction to my earlier note

I previously reported that the provider's base-change route was "rejected" and cited

```
OwnedZZ.hom(PrimeFields()(5)) -> ValueError: the identity morphism exists only on Mor(R, R) ...
```

as evidence. **That was the wrong invocation, not a missing route.** The canonical retained map is the
quotient's own `quotient_map()`, whose domain is the lattice's base ring *by identity*:

```python
ring_map = PrimeFields()(5).quotient_map()      # Integer Ring -> Integer Ring / Ideal (5)
ring_map.domain() is L.base_ring()              # True
```

With that map the whole chain runs, at p=5 as well as p=2:

| step | p=2 | p=5 |
| --- | --- | --- |
| `L.orthogonal_group_base_change(ring_map)` | 0.04 s | 0.03 s |
| `L.base_change(ring_map).Aut()` | 0.00 s | 0.00 s |
| images of the 16–17 generators | 0.11 s | 0.09 s |
| `O(V).subgroup(images)` | 0.18 s | (same path) |

## The gap, stated exactly

`O(V) = O(L (x) F_p)` is finite — `V.base_ring().is_finite()` is `True` — and so is the generated subgroup:

```
subgroup.is_finite()                      -> True      (engine-decided, groups.py:7933)
subgroup in OwnedFiniteGroups()           -> False     (never placed)
subgroup.has_selected_group_resolution()  -> True
subgroup.cardinality()                    -> AssertionError at groups.py:393
```

The refusal is `_group_cardinality` (`categories/group/groups.py:383-394`), which counts a group only when it
*is placed* in `OwnedFiniteGroups()`; failing that, only a chosen free basis or a known-infinite placement
avoids the final

```
assert False, "cannot compute the cardinality of the group ...: it is not known to be finite, ..."
```

So the provider decides finiteness with its engine and then declines to count the group because nothing
*placed* it as finite and no `_selected_cardinality_decision` was installed (`groups.py:7927-7930`). That
channel exists — `_cardinality_decision` first consults `self._selected_cardinality_decision` — it is simply
never selected for a generated matrix subgroup.

## Required semantic operation (minimal contract)

```
input   L          a Lattices(ZZ) object; here rank 10, signature (1, 9), det -1024, even
        p          a prime with p not dividing det L        (here 5 and 7)
        G          a finite tuple of isometries of L        (here 16-17 generators of O(L))
        V = L (x) F_p, nondegenerate                        (guaranteed by p not dividing det L)
        subgroup   O(V).subgroup(tuple(image of g along L.orthogonal_group_base_change for g in G))

output  |subgroup|, the order of the subgroup of the finite classical group O(V) that G generates,
        as an exact integer, computed without a BSGS on the natural module F_p^10 (9 765 625 points at
        p = 5, 282 475 249 at p = 7).
```

Two separate things are missing, and the first alone is not enough:

1. **placement / decision**: a generated matrix subgroup whose `_finiteness_decision()` is `True` should
   either be placed in `OwnedFiniteGroups()` or have a `_selected_cardinality_decision` installed, so that
   `cardinality()` does not fall through to `_group_cardinality`. This is a small, local change.
2. **an order algorithm that is not the engine's BSGS**: if (1) is fixed by routing to `_finite_order(group)`,
   the engine still builds that BSGS. Measured here: `MatrixGroup(...).order()` > 1 000 s at p=5 at a flat
   535 MiB, and `libgap.IsSubgroup` against `GO(+1,10,p)` 5 660 MiB before being killed at 1 200 s.

What *is* available cheaply is the **ambient** order, from GAP's classical group object:

```
libgap.Size(libgap.GO(1, 10, 5)) -> 54 461 311 079 175 000 000 000 000 000 000    (0.00 s)
libgap.Size(libgap.GO(1, 10, 7)) -> 209 544 577 891 300 558 570 289 055 128 616 960 000   (0.00 s)
```

Those are the ambient groups, **not** the image, and this note does not use them to conclude anything about
the image. They are listed because they show what a recognition-based decision could fall back on.

## GAP capability: unavailable, not unloaded

```
GAP 4.14.0
IsBound(IsomorphismTypeOfMatrixGroup)        false
LoadPackage("recog")                         fail      (absent from GAPInfo.PackagesInfo)
LoadPackage("genss")                         fail      (absent)
LoadPackage("forms")                         fail      (absent)
installed: autpgrp aclib alnuth atlasrep congruence cryst crystcat ctbllib fga gapdoc hap polycyclic
           primgrp radiroot smallgrp tomlib transgrp design factint float grape guava io laguna nq polenta
           sonata toric utils
```

The two packages that provide matrix-group recognition — `recog` and `genss` — are **not installed**, so this
is a capability absence, not a symbol that happens to be unbound at start-up. Only `NiceMonomorphism` and
`IsomorphismPermGroup` are bound, and both act on the natural module.

## What is measured, and what is not claimed

**Measured:** the provider chain to the subgroup (≤ 0.2 s at p=2 and p=5); its refusal to count; the exact
locus; that finiteness is already decided `True`; the cost of the engine order and of membership; the absence
of the two GAP recognition packages; the ambient classical orders.

**Not claimed:** that the image is or is not the full orthogonal group at p=5/p=7; that no route exists. The
routes measured here are finite in number and the provider's own ontology has a channel
(`_selected_cardinality_decision`) that nothing currently fills for this object.

## Why this keeps the all-prime scope

The acceptance node still compares computed and recorded image orders at `{2, 3, 5, 7}` and is still open. The
affordable primes agree on both sides for all three cases and are pinned by
`tests/test_lorentzian_reduction_image_orders.py`, which also identifies the order with
`|O^+_10(p)|` wherever `p` does not divide `det L` — a regression, not a substitute for the p=5/p=7 comparison.
