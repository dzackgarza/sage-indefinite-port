# sage-extensions: one 27 s GAP finiteness decision inside a pointed Isom

**Corrected version.** An earlier version of this note said a grep of the pinned tree "finds only the
definition and the two assert sites, neither passing `compute=True`". That was **false**, and the cause
was my own command: I piped the grep through `head -6`. The tree has eleven references, including the
`compute=True` caller at `groups.py:7936`. The owner's correction is right, and the identification
below is what the profile actually supports.

Pinned `sage-extensions` **`e83c4d69db6fba9cadc560d6fd8a518bfe6c5dae`**; line numbers its own. No
library code was changed here.

## Specimen

`test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_-1]`,
orbit representative index 3, vector `[-1,-3,0,0,0,0,0,0,0,0]`, `q = 6`. Uninstrumented 56.12 s;
under `cProfile` 143.24 s, so the shares are attribution, not a runtime claim.

## The decision and its actual caller chain

Measured, innermost first, with call counts and cumulative time from `cProfile`:

```
libgap_mixin.py:181 is_finite                          1 call   27.44 s   (GAP IsFinite)
  <- _engine_finiteness (groups.py:158)                1 call   27.44 s   compute=True branch
    <- _finiteness_decision (groups.py:7924)           1 call   27.85 s   @cached_method
      <- is_finite (set_categories.py:2415)          959 calls   27.86 s
        <- invariant_module_construction
           (rational_integral_stabilizers.py:368)             27.85 s
```

So the answer to "which parent's construction asks" is
**`rational_integral_stabilizers.invariant_module_construction` (`:368`)**. It asks `is_finite` on a
generated rational-matrix subgroup whose placement is not already decided, so
`GeneratedSubgroups._finiteness_decision` falls through its guards to
`_engine_finiteness(..., compute=True)`.

The 959-versus-1 split matters: `is_finite` is asked 959 times for 27.86 s total, of which **27.85 s is
the single undecided parent**. The placement guards and `@cached_method` are doing their job; one parent
is genuinely unplaced.

## What is not being claimed

I am **not** declaring the GAP call avoidable. `_finiteness_decision`'s own contract says finiteness of
a subgroup is independent of its ambient group's, and that engine computation happens only when asked,
never while constructing or placing the subgroup. Both hold here: the call is made once, cached, and
only after the retained finite/infinite placements failed to answer.

In particular, **a finite presentation never proves finiteness** — the subgroup is supplied by
generators, and being finitely generated says nothing about its order. Any argument that the decision
could be skipped has to show that `invariant_module_construction` already knows the answer for this
parent, from a retained placement or from the mathematics of the construction itself. That is the
owner's determination, and the site is `rational_integral_stabilizers.py:368`.

## Suggested owner reading

1. Does `invariant_module_construction` need finiteness of that subgroup for the module it builds, or is
   it asking defensively? If it needs it, the 27 s is real work and the question becomes whether a
   cheaper sufficient criterion exists for this parent.
2. If a placement is available — because the construction already knows the subgroup sits in a finite or
   infinite owner — retaining it as `OwnedFiniteGroups()` / `OwnedInfiniteGroups()` placement would let
   the existing guard answer without the engine.
3. If neither, then no action: the cost is the honest price of a genuinely undecided subgroup.

## Reproduction

`~/sage-image/probe_pointed_unit_profile.py` (specimen pinned and asserted), profile at
`~/sage-image/logs/pointed-q6.prof`, `pstats` transcript at `~/sage-image/logs/pointed-q6.out`. Caller
queries read each `stats` entry's `callers` map, which is how the chain above was obtained rather than
inferred from source.
