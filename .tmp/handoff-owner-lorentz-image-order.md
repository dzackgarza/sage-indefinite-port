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

---

## Pin update: `42d7361c` closes the admission gap (not the engine cost)

Adopted through the maintained checkout at the next safe pin boundary, verified against the actual remote:
`42d7361c` is `origin/main` of `dzackgarza/sage-extensions`, and the staging is `git archive` of that commit.

The commit is `fix: a decided finiteness is placed, and cardinality asks it`, and it closes exactly the gap
reported here: `_group_cardinality` now asks `group.is_finite()` after the placement and free-basis cases
decline, and a decided finiteness is placed. Verified without re-running the engine BSGS:

```
[p=2] is_finite() True in 0.10s | in OwnedFiniteGroups() now True      (was False)
[p=2] cardinality 46998591897600 in 0.02s | harness 46998591897600 | agree True
[p=5] is_finite() True in 0.11s | in OwnedFiniteGroups() now True      (was False)
[p=5] cardinality NOT called
```

Two things worth recording. First, the provider route now *reproduces the harness number at p=2 exactly*
(46 998 591 897 600) through its own semantics — an independent cross-check of `reduction_image_order` at the
prime where both are affordable. Second, at p=5 only admission was checked on purpose: the fix places the
group and routes `cardinality()` to `_finite_order`, which hands the finite matrix group to GAP through
`as_permutation_group()` — the same BSGS cost as before, measured earlier at > 1 000 s. **The p=5/p=7
obligation is unchanged and the nodes stay open.**

## The GAP capability: absence in the measured image, and the exact dependency

The statement "GAP 4.14.0 has `recog` and `genss` not installed" is a fact about the pinned image, not about GAP.
What the Lorentz nodes actually need, and where it comes from, measured inside that image:

* **Capability.** The order of the subgroup generated by the stabilizer's matrix images over `F_p`, i.e. a
  Schreier–Sims (or group-recognition) order computation for a matrix group, used to an exact value at the
  original primes. The contract being reproduced is the recorded image order at `p = 2, 5, 7`, whose prime
  factorisation is `2^21 · 3^5 · 5^2 · 7 · 17 · 31` for the p=2 case — the primes and exponents the acceptance
  assertion names, unchanged.
* **Runtime in the image.** GAP is a *system* installation, not Sage's: `/usr/bin/gap`, root `/usr/share/gap`,
  packages in `/usr/share/gap/pkg`, currently AClib, Alnuth, AtlasRep, AutPGrp, Congruence, Cryst, CrystCat,
  CtblLib, FGA, GAPDoc, HAP, Polycyclic, PrimGrp, RadiRoot, SmallGrp and the rest of the Debian set — no `recog`,
  no `genss`.
* **Supported reproducible route for `genss`.** Sage's optional package `gap_packages` (spkg version 4.16.1,
  `/sage/build/pkgs/gap_packages/spkg-install.in`) installs `genss` into `$SAGE_LOCAL/lib/gap/pkg` as part of its
  documented package list. That is a maintained install (`sage -i gap_packages`), owned by the Sage build, and it
  is the narrowest supported dependency: `genss` supplies randomized Schreier–Sims for matrix groups, which is
  precisely the order computation these nodes lack.
* **`recog` is not in that bundle.** `grep -c '^    recog' spkg-install.in` returns 0: the Sage bundle ships
  `genss` but not `recog`. `recog` (group recognition over finite fields, `RecognizeGroup`) comes from the GAP
  package distribution as a standalone package, installed into the GAP root's `pkg` directory by whoever owns the
  image's GAP (`/usr/share/gap/pkg`). That is a different owner from the Sage build.

So there are two names to route, in this order: the Sage optional `gap_packages` for `genss` (a reproducible
install that also gives a checkable stabilizer chain), and the image's GAP owner for `recog` if recognition
rather than a Schreier–Sims order is required. Neither requires a provider patch, and neither changes the
assertion: the nodes still have to reproduce the recorded order at the original primes.

## What the maintained contracts actually permit for an exact order

Read from the packages' own manuals, because the distinction decides whether a run is evidence.

**genss, chapter 3 (`Stabiliser Chains`), verbatim:** the `DeterministicVerification` component

> "Set this component to true to switch on a deterministic verification routine after the randomised
> Schreier-Sims procedure. **This is not yet implemented.**"

and on the size it returns:

> "the resulting size need not be a divisor of the actual group order. In particular, equality of the computed
> size and the actual group order does not by itself prove correctness unless the group order was known to the
> algorithm during the computation. If the group object G for some reason already knows the group order, then
> this procedure keeps adding random elements until the known size is reached, and returns a correct and proven
> stabiliser chain for G."

So, concretely:

* **Monte Carlo route (not evidence).** `StabilizerChain(G, rec(ErrorBound := e))` and `SizeMC` give a bound on
  the error probability, not a proof. Feeding that number into the acceptance assertion would be probability
  relabelled as proof, and it is excluded.
* **Known-size route (circular, not evidence).** `rec(Size := N)` — or a group object that already stores its
  order — yields a genuinely *proven* chain, but only because `N` was assumed. If `N` is the value the acceptance
  is meant to establish, the argument is circular. It is usable only in the opposite direction: given an
  independently computed `N`, the proven chain *verifies* it.
* **The gap.** `DeterministicVerification` being unimplemented means the package offers no way, today, to turn a
  Monte Carlo chain into a proven one without being handed the answer.

**So the exact order needs a certificate, and the certificate belongs to the maintained library.**

Two candidate routes have been named. Both are **proposals that require actual source and API verification**, and
neither is established evidence today:

1. **Deterministic Schreier-Sims verification of a candidate chain** -- the classical certificate and what
   `DeterministicVerification` would do: check orbit closure at each level and that every Schreier generator sifts
   to the identity, after which `|G| = product of orbit lengths` is exact, with no probability and no assumed
   order. **This is general group-order and certificate machinery, so it is not to be implemented in this port.**
   It is a capability request to the maintained library, routed after that owner's current Enriques priority.
2. **An "independent exact computation" through `reduction_image_order`** -- **retracted as a route.** Checked
   against the source: `reduction_image_order` is not a provider API at all. It is
   `tests/acceptance/consumer.py:109`,

   ```python
   def reduction_image_order(generator_rows, prime):
       return int(MatrixGroup([matrix(GF(prime), rows) for rows in generator_rows]).order())
   ```

   a *harness* helper that calls Sage's plain matrix-group BSGS -- the same algorithmic route as GAP's plain
   `Size`, which is the one measured at over 1000 s on the degree-10 case. It is the reference used to record the
   expected numbers, so it cannot also serve as independent evidence for them, and its p=2 agreement with the
   provider `cardinality()` is an implementation cross-check between two BSGS implementations, not a
   certificate. No number derived from it is admissible as proof, and the expected order is never an input to
   whatever computes the order.

What the capability request therefore asks for, precisely: given generators as matrices over `F_p` (degree 10 in
the measured case), an operation that returns the **exact** order together with a checkable certificate, or an
explicit `fail` -- never a Monte Carlo figure, never a figure that required the expected order as input. The
acceptance keeps its original primes -- `2, 5, 7`, with `2^21 * 3^5 * 5^2 * 7 * 17 * 31` at p = 2 -- and stays
open until such an operation exists and reproduces them.

`recog`'s `RecognizeGroup` returns recognition nodes carrying `Size` fields assembled from directly recognised
leaves, with `immediate verification` checks in the tree. That is a third, weaker candidate: it is constructive
recognition rather than a proof of order, and it can return `fail`. It is recorded here as a fallback to
cross-check against, not as the certificate.
