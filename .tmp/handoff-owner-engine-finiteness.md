# sage-extensions: a 27 s GAP finiteness decision inside one pointed Isom (candidate, not confirmed)

Found by profiling one acceptance unit exactly as the frozen node builds it. Pinned
`sage-extensions` **`e83c4d69db6fba9cadc560d6fd8a518bfe6c5dae`**; line numbers its own. No library
code was changed.

## Specimen (unchanged, exact)

`tests/acceptance/test_vector_stabilizers.sage::test_equivalence_witness_carries_a_representative_to_its_reflection[r_-1]`,
orbit representative at index 3, vector `[-1, -3, 0, 0, 0, 0, 0, 0, 0, 0]`, `q = 6`. Uninstrumented
duration for this exact specimen: **56.12 s**. Under `cProfile` the same unit takes 143.24 s, so the
percentages below are relative attribution, not a runtime claim.

## What dominates — and what does not

Self-time (`tottime`) over the profiled unit, 143.2 s total:

| self | calls | site | judgement |
| --- | --- | --- | --- |
| **58.39 s** | 4 | `{readline}` on the Julia bridge pipe | **necessary mathematics** — `lattice_engines.definite_isometry_group` (`:1146`) computing the definite isometry group in the engine, one call, cum 58.39 s |
| **27.44 s** | 1 | `sage/groups/libgap_mixin.py:181 is_finite` | **candidate** — a GAP finiteness decision, caller `_engine_finiteness` (`groups.py:158`) |
| 2.57 s | 4 154 252 | `ring_foundation.py:5867 _from_engine_element` | volume churn in the owned object model |
| 2.30 s | 3 247 511 | `objects.py:158 __contains__` | " |
| 2.19 s | 4 615 498 | `sage/categories/category.py:436 __classcall__` | " |
| 2.03 s | 14 376 823 | `{isinstance}` | " |
| 1.74 s | 5 052 819 | `ring_foundation.py:5896 _element_constructor_` | " |

**There is no eager-representation defect in this unit.** Total time inside `repr` is **1.38 s** over
139 205 calls, and every `_repr_` method is ≤ 0.03 s. This is the opposite of the Enriques case-1
finding, and it is worth stating plainly so the two are not conflated: this plateau is engine-bound
mathematics plus broad object-model overhead, not string formatting.

## The candidate

`_engine_finiteness` (`groups.py:158`) is:

```python
def _engine_finiteness(engine, *, compute=False):
    if engine in SageFiniteGroups():
        return True
    if engine.category().is_subcategory(SageGroups().Infinite()):
        return False
    if compute and isinstance(engine, MatrixGroup_generic) and engine.base_ring() in (ZZ, QQ):
        return bool(engine.is_finite())
    return None
```

The 27.44 s call reaches `libgap_mixin.is_finite`, i.e. GAP's `IsFinite` on a GAP-backed group. The
in-provider entry points that reach a finiteness assertion are the two `assert _engine_finiteness(engine) is True`
sites in `_gap_model_of_engine` — `groups.py:191` for `AbelianGroup_class` and `groups.py:198` for
matrix groups, both immediately before `libgap(... .permutation_group())`.

**What I can state from measurement:** one GAP finiteness decision costing 27.44 s occurs inside a
single pointed Isom, and the profiler attributes it to `_engine_finiteness`.

**What I cannot state without the owner's reading:** whether the expensive branch is the explicit
`compute=True` call or a category-membership refinement inside `engine in SageFiniteGroups()`. A
grep of the pinned tree finds only the definition and the two assert sites, neither of which passes
`compute=True`, so the trigger is not obvious from source alone and I am not going to assert one.

**Why it matters if avoidable:** the group handed to GAP there is being converted to a permutation
representation, which only makes sense for a finite group — so finiteness is already implied by the
caller's intent. A construction-time flag or a cheaper argument would remove 27 s of a 56 s unit.

## Reproduction

`~/sage-image/probe_pointed_unit_profile.py` (specimen pinned and asserted in the script) with the
profile kept at `~/sage-image/logs/pointed-q6.prof` and the full `pstats` output at
`~/sage-image/logs/pointed-q6.out`. Run with the port repository bound at `/work` and the provider at
`/repos/sage-extensions`, then `pstats.Stats(...).sort_stats("tottime").print_stats(16)`.
