# sage-extensions: two `is_finite()` calls whose answer cannot change the branch

Pinned `sage-extensions` **`c546077dce79a76eb32e92d080aa9c44e4c67552`**; line numbers its own. No library
code changed here. Found from a port profile of one acceptance unit, and it is a **general-owner** matter:
the port supplies the containing bound that makes both queries dead; the ordering is the provider's.

## The two sites

`categories/rational_integral_stabilizers.py:371` — `IntegralStructureAction.invariant_module_construction`:

```python
group = self.rational_group()
containing = self.stable_containing_bound()
if group.is_finite() is not True and containing is None:      # <- line 371
    assert False, (
        f"computing the orbit-generated integral module for {group} requires a known finite group "
        "or a justified finite stable containing module; supply the latter with "
        "with_stable_containing_bound(...)"
    )
return GObjects(group, Modules(...)).orbit_generated_submodule(..., containing=containing)
```

`categories/modules/orbit_spans.py:44`:

```python
if group.is_finite() is True and containing is None:          # <- line 44
    span = ambient.subobject_on(tuple(act(vector) for act in acting_maps(group) for vector in generators))
else:
    ...  # bounded closure, requires containing and a finite generating set
```

Both put `group.is_finite()` on the **left** of a conjunction whose other operand is
`containing is None`. Python evaluates left to right, so the finiteness query runs first, and when a
stable containing bound **is** supplied the right operand is `False` — the conjunction is `False` at both
sites whatever the group's order is, and the `else` branch of the second is taken regardless.

**The answer cannot affect the outcome, yet it is computed.**

## Why this is the port's normal path, not an edge case

`IndefiniteRecursion._lattice_action` (`src/sage_indefinite_port/indefinite/recursive.py:445-468`) always
supplies the bound:

```python
containing = Modules(lattice.base_ring()).Mono(split, space)(lambda label: space.wrap(...))
return IntegralStructureAction(rational.O().subgroup(extensions), inclusion).with_stable_containing_bound(containing)
```

So in this port `stable_containing_bound()` is never `None` on the recursion path, and every
`invariant_module_construction` / orbit-span call made beneath it pays for a query it cannot use.

## Measured

From a `cProfile` of one pointed-isometry unit (root `r_-1`, representative `q = 6`, 56.12 s
uninstrumented):

```
libgap_mixin.py:181 is_finite                     1 call   27.44 s   (GAP IsFinite)
  <- _engine_finiteness (groups.py:158)           1 call   27.44 s   compute=True branch
    <- _finiteness_decision (groups.py:7924)      1 call   27.85 s   @cached_method
      <- is_finite (set_categories.py:2415)     959 calls   27.86 s
        <- invariant_module_construction
           (rational_integral_stabilizers.py:368)          27.85 s
```

27.44 s of a 56 s unit, for an answer the branch above discards. Of 959 `is_finite` calls in that unit,
this one undecided group accounts for 27.85 s; the placement guards and caching are doing their job.

`_finiteness_decision` (`groups.py:7924`) is correct as written — it returns the retained
`OwnedFiniteGroups()` / `OwnedInfiniteGroups()` placement when one exists and only falls through to
`_engine_finiteness(_engine_group(self), compute=True)` when none does. The defect is not that decision;
it is that **two consumers ask a question they cannot use**.

## Avoidable lost premise, not a required computation

This is the distinction that matters for the owner:

- **Avoidable:** at both sites, when `containing` is supplied, `group.is_finite()` is dead. Reordering
  the conjunction to put `containing is None` first preserves the truth value exactly — a pure
  evaluation-order change — and makes the query run only where it can change the result.
- **Genuinely required:** the same query *is* required when no bound is supplied, because then it
  decides between unbounded orbit closure and refusing the computation. That path is untouched by the
  reorder.

Proposed, minimal:

```python
# rational_integral_stabilizers.py:371
if containing is None and group.is_finite() is not True:
```
```python
# orbit_spans.py:44
if containing is None and group.is_finite() is True:
```

This is not a shortcut around `IsFinite` and it introduces no finiteness assumption: the predicate is
unchanged, only deferred until it can matter. The port keeps its own responsibility — supplying the
justified stable containing bound — which it already does.

## Smallest reproduction

Provider-only, no port code, and it is a call counter rather than a timing claim:

1. Wrap `OwnedGroups.ParentMethods.is_finite` with a counter.
2. Build `IntegralStructureAction(G, i)` for a small rational lattice with a subgroup `G` of a couple of
   supplied isometries, and attach a stable containing bound with
   `.with_stable_containing_bound(bound)`.
3. Call `invariant_module_construction()` twice — once with the bound attached, once without
   (`containing is None`).
4. Expect: **counter increments on the bounded call even though the branch outcome is identical** with
   `is_finite` forced `True` or forced `False`; and the unbounded call genuinely needs the answer.

A port-level instance of the same query is
`IndefiniteRecursion.isotropic_stabilizer_action((v,), ())` on `U` with `v = e0 + e1`, whose returned
action has `stable_containing_bound()` non-`None` by construction.

I have not executed the counter reproduction — one heavy acceptance run is in flight — so this is a
source-derived claim with a measured call stack behind it, not a measured counter.

## Provenance (added: the earlier note alone did not carry this)

- **Artifact**: `.tmp/handoff-owner-finiteness-shortcircuit.md`, commit
  **`75dc17f6be36559742ffb7dbfa7beeb18646977a`** ("docs: finiteness handoff -- two consumers ask an
  unusable question; first centralizer case passes"), an ancestor of the current `main`.
- **Measurement**: `~/sage-image/logs/pointed-q6.prof` with the transcript
  `~/sage-image/logs/pointed-q6.out`, produced by `~/sage-image/probe_pointed_unit_profile.py`. The unit
  is root `r_-1`, orbit representative index 3, vector `[-1,-3,0,0,0,0,0,0,0,0]`, `q = 6`:
  **56.12 s uninstrumented**, 143.24 s under `cProfile`, of which **27.44 s** is GAP `IsFinite`.

## The minimal owned group and its construction

The group the discarded query is asked about is built in the port at
`src/sage_indefinite_port/indefinite/recursive.py:456-468`, `IndefiniteRecursion._lattice_action`:

```python
perpendicular = vector.orthogonal_complement().inclusion()
complement    = perpendicular.domain()
extensions    = tuple(generator.extension_across(vector, vector)
                      for generator in self.orthogonal_group_generators(complement))
rational      = lattice.base_change(lattice.base_ring().fraction_field_map())
inclusion     = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()
split         = lattice.subobject_on((vector, *(perpendicular(point) for point in complement.module_generators())))
denominator   = split.inclusion().cokernel().exponent()
scale         = rational.base_ring().one() / rational.base_ring()(denominator)
containing    = Modules(lattice.base_ring()).Mono(split, space)(lambda label: space.wrap(...))
return IntegralStructureAction(rational.O().subgroup(extensions), inclusion).with_stable_containing_bound(containing)
```

So: an owned **generated subgroup** `rational.O().subgroup(extensions)` of the rational orthogonal
group, acting on the generic fibre of `L`, with a **stable containing bound supplied** — which is exactly
the premise that makes the two `is_finite()` calls dead.

That is the smallest instance to reproduce with: `U` with `v = e0 + e1` (`q = 2`) reaches the same
construction through `isotropic_stabilizer_action`/`_lattice_action`, and its
`stable_containing_bound()` is non-`None` by construction.
