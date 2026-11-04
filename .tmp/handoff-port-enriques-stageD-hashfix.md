# Stage-D after the provider hash fix: the wall is gone, and the returned group is proper (index unproved)

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

— every one is in the **pointwise** stabilizer. `Aut((N, Z e1))` is strictly larger, and the sign-reversing
semantic witness is explicit. On the fixture basis (block diagonal in the `U, U(2), E8(-2)` splitting), take

```
M = diag(-1, -1, 1, ..., 1)          # 12 x 12, two negative diagonal entries, the rest +1
```

with, all checked exactly against the fixture gram in this repository and no run required:

| property | value |
| --- | --- |
| integral | yes |
| `M G M^T = G` (row convention) | **true** |
| `det M` | **+1** (`(-1)(-1) = +1`; orientation-preserving) |
| `M e1` | `-e1` |

**Correction to the previous version of this note:** it asserted `det M = -1`. That was wrong, and it was written from
a hardcoded print rather than a computation: with two negative diagonal entries the determinant is `+1`. The
conclusion does not depend on the determinant — what matters is that `M` is a genuine integral isometry of `N`
negating the line vector, and that is unaffected. It does refine the statement: the sign-reversing element is in
`SO(N)`, so no orientation-preserving restriction can exclude it.

A generating set whose elements all act as `+1` on the line generates a subgroup of the **pointwise** stabilizer, while
`M` lies in the setwise stabilizer and not in that subgroup. Therefore the delivered 25 generators do not generate
`Aut((N, Z e1))`: the returned subgroup is **proper**.

**What that does *not* prove.** It does not establish the index. Exact index two would additionally require that the
returned subgroup is the *full* pointwise stabilizer, and no such completeness proof is available here: the only
route to the pointwise stabilizer in this port is
`pointwise_stabilizer_generators` → `pointed_automorphism_generators` → the same
`isotropic_stabilizer_action(...).stabilizer()` construction, so comparing the two would not be an independent check.
The previous version of this note said "index two short"; that is retracted, and the claim published here is
**proper, index unproved**. An index would need an independent generator set or order for `Stab_{O(N)}(e1)`, which
would have to come from the owning layer.

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

## Retained raw measurement and guard output

The full captured output of the run is tracked verbatim at `.tmp/enriques_stageD_hashfix_e609d032.log` (4018 bytes):
the provenance lines above, every stage line, the 25-generator result, the per-generator sign list, the
form-preservation and setwise checks, and the `peakrss` block. Its guard figures are the ones quoted here — peak tree
RSS **496 MiB** against the 6144 MiB limit, `child waited=0`, `shell_code=0` — together with the `low-water MemAvail`
line the supervisor wrote. The three stack dumps it contains are the residual-cost frames named below; no dump reaches
the old hash path.

## Kept open, unchanged

- `case_3` acceptance: line/plane/flag counts **5/10/19** and generators that are genuine integral isometries — the
  count obligation is untouched by this measurement, and the generator obligation is now **partly** met (integral,
  form-preserving) and **not** met (the returned subgroup is proper: it misses the explicit sign-reversing witness).
  The size of the shortfall is unproved, and the setwise group contract stands as stated.
- The Enriques family's other stages and the Z/2-shortfall above are port-consumer obligations, not discharged here.
- No provider patch was made in this repository; the semantic finding is routed to `sage-extensions`, which owns
  `rational_integral_stabilizers.py`.
