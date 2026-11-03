# Stage-D after the provider hash fix: the wall is gone, and the returned group is index two short

Continues `.tmp/handoff-port-enriques-stageD-attribution.md`. Measured on provider
`e609d03221be54240a510f29c7fd6950f7125965` (contains `46b2c392`, native tensor hash), port `3e2bbd51` + this
commit, same published input: direct-gram `N` from `unpolarized_enriques.json` (rank 12, signature (2,10), det 1024),
line `I_1_1` `Z e1`. Guard 6144 MiB process-tree. No 40-minute loop repeated, no cap extended, no competing suite.
The project's job list showed `running: 0` before the run, so the heavy slot was free.

## Provenance of the tree that actually ran

```
provider package:                 /repos/sage-extensions/src/sage_extensions/__init__.py
lattice_morphisms.py sha256[:16]  00f0c682880cc022
```

The bound tree is the fetched `e609d032` checkout, and the digest is of the file the run imported — not of the
port's `pyproject` pin, which is unchanged. Stages run once each, in order, so the reduction and `O(K)` are **warmed**
before the action; `stabilizer()` is the same call cold or warm.

## The measurement

| stage | before (`0dad201c`) | after (`e609d032`) |
| --- | --- | --- |
| `isotropic_reduction` | 1.1 s | 1.1 s |
| `orthogonal_group_generators(K)`, K rank 10 (1,9) | 2.1 s, 22 gens | 2.1 s, 22 gens |
| `isotropic_stabilizer_action` (whole) | 16.7 s | 24.6 s |
| `action.stabilizer().group_generators()` | **never returned (2400 s cap)** | **31.7 s, 25 generators** |
| peak tree RSS | 1430 MiB | 496 MiB |

So the attributed hot path was the wall: with `LatticeIsometryMor.__hash__` off the repr-built cache key, the orbit
enumeration that never finished now finishes in half a minute. Confirmed as the reason: the previous run's dumps were
all in `from_finite_orbits -> _orbit_words -> _act_on_position -> quotient_automorphism -> LatticeIsometryMor.__hash__`;
those frames no longer appear.

## What the returned generators are — and what they are not

All 25 are genuine integral isometries of `N` (`M G M^T = G` verified through the port's `is_isometry`) and all fix
the line **setwise**. But their action on the line is

```
line action per generator: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
```

— every one is in the **pointwise** stabilizer. `Aut((N, Z e1))` is strictly larger: the element `M = diag(-1,-1,1,
...,1)` on the fixture basis is integral, `M G M^T = G`, `det M = -1` and `M e1 = -e1` (checked directly against the
fixture gram, which is block diagonal in the `U, U(2), E8(-2)` splitting). A generating set whose elements all act as
`+1` generates a subgroup of the pointwise stabilizer, so the returned 25 generators cannot generate the setwise line
stabilizer: the delivered group is index two short.

That is the remaining owning-library call, and it is a **semantic** one, not a cost one:
`IntegralStructureAction.stabilizer()` (`rational_integral_stabilizers.py:459`) returned through
`lattice_stabilizer` (:242), i.e. `G ∩ O(L)` for the constructed rational group `G`. Either the `+̶1` image family that
`isotropic_stabilizer_action` passes as `radical_images=((-e1,),)` never reaches `G`, or `lattice_stabilizer` drops it.
The port's own contract for `line_automorphism_generators` ("`Gamma = {+-1}` on the line") is therefore not delivered
on this input, and `case_3`'s consumer needs it.

This is not something the port may paper over: an ambient `Isom(L, L)` or a decision oracle cannot stand in for this
subgroup action / orbit / stabilizer-generator contract, and no such substitution was made.

## Residual cost, named but not blocking

The dumps during the (now bounded) 31.7 s land in `orbit_schreier_morphism` (`g_sets.py:520`) →
`_orbit_schreier_morphism` (:75) → `finite_ordered_sets.position_of` → `contains` → `LatticeIsometryMor.__eq__`
(`lattice_morphisms.py:1482`), which compares morphisms by linear combination. That is a linear scan per orbit point,
so the remaining cost is quadratic in the orbit; it is a cost item for the library, not a wall, and it no longer
prevents the call from returning.

## Kept open, unchanged

- `case_3` acceptance: line/plane/flag counts **5/10/19** and generators that are genuine integral isometries — the
  count obligation is untouched by this measurement, and the generator obligation is now **partly** met (integral,
  form-preserving) and **not** met (setwise line action).
- The Enriques family's other stages and the Z/2-shortfall above are port-consumer obligations, not discharged here.
- No provider patch was made in this repository; the semantic finding is routed to `sage-extensions`, which owns
  `rational_integral_stabilizers.py`.
