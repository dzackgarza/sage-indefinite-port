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

---

## Follow-up: a third site, and the repair does not cover it

Observed on the repaired pin `cdafca84`, running the `square[0]` stabilizer case as a full consumer check:

- It **hit a 6 144 MiB memory guard at 191 s** (peak tree RSS 6 193 MiB, wrapper exit 143). The case is
  unchanged from the one that previously died at ~9.8 GB, so **the bounded-finiteness repair does not
  change this case's memory profile.**
- The 120 s stack was still inside the GAP query, but through a **different site**:

```
libgap_mixin.py:195 is_finite
  groups.py:170 _engine_finiteness            (compute=True branch)
  groups.py:7943 _finiteness_decision
  set_categories.py:2423 is_finite
  rational_integral_stabilizers.py:444 stabilizer      <-- not covered by 00bcc86
  work/src/sage_indefinite_port/indefinite/recursive.py:593 pointed_automorphism_generators
  marked_lattices.py:109 automorphism_group_generators
  ...
```

`IntegralStructureAction.stabilizer()` (`:436-446`) has the same shape:

```python
if (
    rational_group.is_finite() is True
    or (
        rational_group.has_selected_group_resolution()
        and rational_group.has_selected_finite_group_generating_set()
        and self.stable_containing_bound() is not None
    )
):
    return self.finite_representation().lattice_stabilizer()
```

If the second disjunct holds, the condition is `True` whatever `is_finite()` answers — so the query is
again dead work, and the fix is the same pure reorder, moving the composite disjunct first:

```python
if (
    (
        rational_group.has_selected_group_resolution()
        and rational_group.has_selected_finite_group_generating_set()
        and self.stable_containing_bound() is not None
    )
    or rational_group.is_finite() is True
):
```

**What the owner still needs to check, which I have not established by measurement:** whether the second
disjunct is actually `True` for the group this path builds. It requires a selected group resolution, a
selected finite generating set, and a non-`None` containing bound. The port supplies the bound
(`_lattice_action`), so the first and third are plausible — but "plausible" is not measured, and the
reorder only helps if all three hold.

Direct check: instrument `IntegralStructureAction.stabilizer` to print the three conjuncts for the group
reached from `recursive.pointed_automorphism_generators` on the `square[0]` specimen.

## Why this is routed, not fixed here

Both the `or` and the two conjunctions are provider code; the port's contribution is the bound, which it
already supplies. And the failure is a **memory** failure, not a wrong answer — the case never reaches its
assertions, so there is no port-side assertion to repair.

The two acceptance records are untouched by this: A_PERFECT 92/92 at `87fbf863` + `cdafca84` and
A_ISOMETRY 166/166 at `87fbf863` + `c546077d` stand on their own pins and were not re-run.

---

## A rank-2 reproduction of the memory growth

While checking whether the `stabilizer()` disjunct holds, the probe never reached its first assertion:
`EichlerRecursion._lattice_action` on the **hyperbolic plane** `U = [[0,1],[1,0]]` with `v = e0 + e1`
(`q = 2`) exceeded a 3 GiB tree-RSS guard after ~48 s:

```
timelike v q = 2
[peakrss] STOPPED_AT_MEMORY_LIMIT: tree RSS 3233 MiB exceeded 3072 MiB
[peakrss] peak tree RSS : 3233 MiB   low-water MemAvail : 9023 MiB
```

That is rank 2 — the smallest Lorentzian lattice there is — and the blow-up is inside the port's
`_lattice_action` before any E10 or discriminant-form machinery is touched. So the memory behaviour that
sinks `square[0]` is **not** a rank-10 phenomenon, and a reproduction does not need the 9.8 GB specimen or
a 6 GiB budget.

`_lattice_action` (`recursive.py:445-468`) builds, in order: `orthogonal_group_generators(v^perp)` for the
rank-1 complement (cheap: `O(<-2>) = {±1}`), `extension_across` for each generator, the base change to
`QQ`, `generic_fibre_map`, `lattice.subobject_on(...)`, `split.inclusion().cokernel().exponent()`, and a
`Modules(...).Mono` bound built by a lambda over `split`'s generators. One of those is allocating hundreds
of megabytes on a two-dimensional lattice; the suspects worth instrumenting first are the cokernel exponent
and the `Mono` over the generic fibre, because both build module objects over the *rational* space.

I have not isolated which line it is — the probe was stopped by my own guard, which is the correct
behaviour on a shared host, and I am not going to guess at the owner's site. The value here is that the
reproduction is now cheap: a two-line lattice, one vector, one call.

Next step for whoever takes it: run `_lattice_action` on `U` with `v = e0 + e1` under a line-level memory
profile (`VmRSS` deltas around each of the six construction steps above) — the whole thing is under a
minute.

---

## RETRACTION: the rank-2 "memory blow-up" was the Julia child process

My previous note claimed the memory growth "reproduces on rank 2, inside `_lattice_action`". **That was
wrong**, and the error was mine: the sampler summed RSS over the process *tree*, and the Julia bridge is a
child whose baseline footprint is ~3 GiB.

Instrumenting each of the nine construction steps of `_lattice_action` on `U` with `v = e0 + e1`, printing
the process's own `VmRSS` after every step:

```
1 perpendicular inclusion              0.11s  273 ->   277 MiB  (+4)
2 complement domain                    0.00s  277 ->   277 MiB  (+0)
3 orthogonal_group_generators [ext]   62.87s  277 ->   383 MiB  (+106)   <- the Julia engine call
4 base_change(QQ)                      0.00s  383 ->   383 MiB  (+0)
5 generic_fibre_map                    0.02s  383 ->   384 MiB  (+1)
6 subobject_on(v, perp generators)     0.01s  384 ->   385 MiB  (+1)
7 split cokernel exponent              0.03s  385 ->   385 MiB  (+0)
8 stable containing bound (Mono)       0.00s  385 ->   385 MiB  (+0)
9 IntegralStructureAction(...)         0.06s  385 ->   387 MiB  (+2)
FINAL own RSS 387 MiB after 63 s
```

Running the same probe under the tree sampler: **FINAL own RSS 387 MiB, peak tree RSS 3 576 MiB.** The
difference is the Julia bridge, not the port. Step 3 is the only expensive step and it is expensive in
*time* (62.87 s), not in memory (+106 MiB).

So the port allocates ~114 MiB in total on this path at rank 2. **There is no rank-2 memory reproduction,
and I should not have said there was.**

### What this means for `square[0]`

It also means my `square[0]` figures need re-reading. The sampler reports `peak single VmHWM` by summing
each process's `VmHWM` across the tree — it is **not** a single-process maximum, despite the label. So
"peak single VmHWM 9 801 MiB" for `square[0]` does not establish that Python reached 9.8 GB; part of it,
possibly most, may be the Julia engine.

**Unresolved and now clearly mis-specified:** whether `square[0]`'s growth is Python's or Julia's. Nobody
should route that to the provider, or to the port, on the strength of a summed figure. The sampler must
report per-process peaks *by executable* before the question is asked again.

That is the next instrumentation task, and it is small: record `VmHWM` per pid together with
`/proc/<pid>/comm`, and report the largest Python and the largest Julia separately.
