# sage-extensions: hashing a torsion-form automorphism costs a GAP round-trip

Pinned `sage-extensions` **`9dc89ac1d9e4ab074da8872386aaf5b9df290f63`**; line numbers its own. No
library code changed here.

## Context: the eager-description fix is confirmed working

`_symmetric_group`'s degree-based description removed the earlier blocker. In a 3600 s instrumented run
of `test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]` with 29 two-minute stack dumps,
**zero dumps land in any `_repr_`**, against a pre-fix run whose stacks sat in
`finite_ordered_sets.py:587 _repr_` for ~1500 s. That defect is closed.

## What replaces it

The case still does not finish: `pytest_exit=124 elapsed=3600s`, test frame still
`test_boundaries.sage:226`. The phase table at t=3543 s:

| phase | calls | total | note |
| --- | --- | --- | --- |
| `group_generators` | 69 | 162.85 s | plus **one call RUNNING 3374 s** |
| `discriminant_morphism` | 210 | 37.03 s | |
| `cardinality` | 806 691 | 5.01 s | |
| `subgroup_on` (nested) | 1 | 0.02 s | |

So one `group_generators` realisation has been running for 56 minutes. **14 of the 16 complete dumps**
have their top frame inside the orbit-word computation:

```
pointed_orbit            g_sets.py:190
  _orbit_words           g_sets.py:156
    <lambda>             predicate_subgroups.py:82
      act                g_objects.py:1452
        action_of        g_sets.py:319
          __call__       groups.py:2982
            _permutation_from_point_map   g_sets.py:1222
              <lambda>   g_sets.py:1283
                act      groups.py:528
                  __hash__   torsion_form_modules.py:1164
                    matrix   abelian_aut.py:203
                      __call__ fqf_orthogonal.py:114 / abelian_aut.py:166  -> GAP
```

## The concrete cost

`TorsionFormAutomorphism.__hash__` (`torsion_form_modules.py:1164`):

```python
def __hash__(self) -> int:
    return hash((
        id(self.parent()),
        tuple(int(entry) for entry in self._engine().matrix().list()),
    ))
```

`self._engine()` is a Sage `FqfOrthogonalGroup` element, and `.matrix()`
(`abelian_aut.py:203`) converts the GAP automorphism to a matrix, driving
`abelian_group_gap` (`exponents`, `_element_constructor_`). So **every hash of a torsion-form
automorphism performs a GAP conversion.**

That matters because `_permutation_from_point_map` builds a mapping over the action's points, so the
cost is paid per point per orbit step — which is exactly where 14 of 16 dumps sit.

The hash is *correct*: it is derived from the same engine data as `__eq__`, so equality and hashing
agree. This is a cost question, not a correctness one.

## Question for the owner

Can the hash be computed from data the element already holds, instead of forcing a GAP matrix
conversion — for example from the element's engine coordinates captured at construction, or from a
memoised matrix? If the element's `_engine()` cannot give a stable cheap digest, then the cost is
inherent and the orbit computation needs a different point key; either way the decision is at
`torsion_form_modules.py:1164` and `g_sets.py:1222`.

I am **not** asserting this is avoidable. It is the frame the run sits in, consistent across 14 of 16
samples, with a per-element GAP conversion visible in the source.

## Reproduction

Run `test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]` with the
`discprof` plugin (phase table every 60 s, stacks every 120 s); the log is
`~/sage-image/logs/`-external, reproduced from `.nucleus/jobs/jmvxn1hj8.log`. A cheaper direct probe:
time `hash(automorphism)` for a few thousand elements of a `TorsionFormOrthogonalGroup`.
