# Port-owner attribution of the case_3 stage-D block: it is `integral_stabilizer.stabilizer`, not the induction

Measured on the exact published input, at provider `0dad201c` (fetched tree `.scratch/prov/sage-extensions-main`);
port `ee4c2901` plus the working tree of this commit; direct-gram `N` from
`tests/fixtures/unpolarized_enriques.json` (rank 12, signature (2,10), det 1024) and the published line
representative `I_1_1` (`Z e1`, basis `[1,0,...,0]`); guard 6144 MiB process-tree; **no 40-minute run repeated and
no cap extended**.

## The sub-stage table (focused, distinguishing measurement)

Reproduce with `.scratch/probe-stage-attrib.py`; `STACK_EVERY=20` arms `faulthandler` so a cap still names its frame.

| stage | measured | verdict |
| --- | --- | --- |
| `N` built, `I_1_1` element | 1.5 s | — |
| `primitive_isotropic_subobject(v).isotropic_reduction()` | **1.1 s** | not the block |
| quotient `K = v^perp/<v>` | rank 10, signature (1,9), det 1024 | — |
| `AttackProfile.from_lattice(K)` | positive index 1 | Lorentzian leaf selected |
| `orthogonal_group_generators(K)` | **2.1 s, 22 generators** | not the block |
| `reduction.pointwise_perpendicular_kernel()` | **0.7 s, 0 generators** | not the block |
| `isotropic_stabilizer_action((v,), ((-v,),))` (whole construction) | **16.7 s** | not the block |
| `action.stabilizer().group_generators()` | **does not return** (300 s cap here; 2400 s in the retained run) | **the block** |

So the non-returning operation is *inside* `.stabilizer()`. The induction, the reduction, the Lorentzian group and
the kernel are all fast on this input: the earlier framing "the innermost
`line_automorphism_generators(N, e1)` never returns" is true of the *call*, but the intractability is not in the
induction and not in `O(K)`.

## The frame, verbatim (repeated dumps, one location)

```
rational_integral_stabilizers.py:459  stabilizer
rational_integral_stabilizers.py:437  finite_representation
rational_integral_stabilizers.py:132  __init__           (from_finite_orbits)
group/g_sets.py:172                   from_finite_orbits
group/g_sets.py:156                   _orbit_words        (breadth-first over the orbit)
rational_integral_stabilizers.py:179  _act_on_position
rational_integral_stabilizers.py:419  quotient_automorphism
lattice_morphisms.py:1505             LatticeIsometryMor.__hash__
tensors/tensor.py:1949                _repr_
tensors/tensor.py:2177                _coordinate_component_repr
rings/ring_foundation.py:1218         _repr_
rings/ring_foundation.py:6446         one
rings/ring_foundation.py:5498         __init__           (owned engine element construction)
```

and, at the same `_act_on_position`, a second dump through

```
rational_integral_stabilizers.py:192  _act_on_position
modules/.../finitely_presented_modules.py:736  subobject_on
modules/.../finitely_presented_modules.py:3399 _presented_module_from_morphism
modules/pure/modules.py:2498          _framing_lift
```

## Attribution, stated as the contract the caller depends on

`IntegralStructureAction.stabilizer()` takes its **finite-representation branch** (line 437) because the port supplies
`with_stable_containing_bound(...)` *and* the rational group has a selected resolution and finite generating set.
That branch builds `_FiniteIntegralRepresentation`, whose `from_finite_orbits` runs the BFS of `g_sets._orbit_words`
over the orbit, one word per point. Each step calls `_act_on_position`, which (a) constructs a presented submodule
(`subobject_on` / `_presented_module_from_morphism`) and (b) hashes the resulting `LatticeIsometryMor`, whose
`__hash__` (`lattice_morphisms.py:1505`) reaches a **tensor `_repr_`** and constructs fresh owned ring elements
inside it. The per-element cost is therefore not arithmetic but representation and module-bookkeeping, and the
number of elements is the orbit size.

Two things follow, and they are the port-side contract ask to the owning layer:

1. **`LatticeIsometryMor.__hash__` must not build a display string.** `_coordinate_tensor()._cache_key()` currently
   reaches `tensor._repr_` → `_coordinate_component_repr` → ring `_repr_` → owned element construction. This is the
   same defect class already recorded for `TorsionFormAutomorphism.__hash__`; here it is on the hot path of every
   orbit step. Hashing the tensor's coordinates (or an interned key) without formatting fixes it.
2. **The finite representation needs a lazy or decision-based orbit, or a caller-supplied exact orbit relation.**
   The port decides `Isom(L, L)` exactly (`EichlerRecursion.isometry`) and would happily supply
   `orbit_relation`, but `_FiniteIntegralRepresentation` only accepts an enumerated orbit from `_orbit_words`. Since
   the *stabilizer* is what the caller wants, and the orbit is materialised only to take a stabilizer in it, a
   relation-driven route returns the same object without the enumeration.

Until either lands, `action.stabilizer().group_generators()` is the operation that does not return on this input, and
that is a property of the finite-representation enumeration, **not** evidence that `Aut((N, Z e1))` is intractable.

## Corrections to the published handoff

- Its second bullet says "flat at 1430 MiB RSS ... host MemAvail ~11 GiB throughout ... no memory component". The
  retained log it cites carries `low-water MemAvail : 526 MiB` (final samples 11054 MiB). The low-water reading is
  the one that matters: a host that came within 0.5 GiB of exhaustion during a 2400 s run means memory pressure
  cannot be excluded as a contributing factor, and "~11 GiB throughout" is an inference from the last samples only.
  Recorded here as a correction, not as a claim of an OOM.
- Its "attribute inside `isotropic_stabilizer_action` vs `stabilizer()` vs `integral_isometry`" ask is now answered:
  `isotropic_stabilizer_action` 16.7 s, `stabilizer()` non-returning, `integral_isometry` never reached.
- Its isolation scoping and pins stand as recorded; the launch `c2288184` provenance is a live-bind reading, not an
  immutable snapshot, and is used here only as provenance.

## What is not claimed, and stays open

- No provider patch was made here: the two items above are routed to `sage-extensions`, which owns that layer.
- No smaller stabilizer, no finite-index substitute and no reduced contract were adopted. The original case_3
  acceptance — line/plane/flag counts 5/10/19 and generators that are genuine integral isometries of `N` — remains
  open on its own terms.
- The 40-minute run was not repeated and the 2400 s cap was not extended; the 300 s probe above is the distinguishing
  measurement, and it names the frame rather than timing the wall again.
