# Measured: Enriques case_87 fails in the provider coset/Schreier path, not in the port

**Owner:** `sage-extensions` (`categories/group/predicate_subgroups.py` and the coset machinery).
**Port side:** none. The port's `test_boundaries.sage` only calls the public consumer route.

## What happens

`test_enriques_modular_groups_match_dutour_sikiric_hulek[case_87]` at port `1bedc833` + provider
`cdafca84` fails in ~140 s (not a cap) at

```python
image = target.subgroup_on(tuple(g.discriminant_morphism() for g in group.group_generators()))
```

where `group = _gamma_h(polarization)` ends in `.intersection(anti_invariant.O_plus())`, i.e. the
subgroup `{g in O(N) : g lies in every selected subgroup}`. Asking **that** object for
`group_generators()` traces:

```
groups.py:5438            group_generators
groups.py:5416            selected_group_resolution -> selected_resolution(OwnedGroups())
objects.py:383/116        _PendingResolution.realize -> resolution factory
predicate_subgroups.py:47   generating_morphism(subgroup)
predicate_subgroups.py:1100 _restricted_preimage_generating_morphism
predicate_subgroups.py:1138   restricted.schreier_generating_morphism()
predicate_subgroups.py:821      coset_space(generators=None)
predicate_subgroups.py:811        _preimage_coset_space
predicate_subgroups.py:83           _free_cover(ambient, generators=None)
predicate_subgroups.py:62             ambient.selected_group_resolution().augmentation()
```

and terminates in

```
E   KeyError: ((Subgroup of Orthogonal group of Finitely presented module on 12 module generators
                over Integer Ring with invariant factors [0 ↦ 2, ..., 9 ↦ 2],), ())
E   OverflowError: value too large to convert to int
E     (raised in repr() of a _TransversalCoset; sage/libs/gap/element.pyx:3166)
```

The `KeyError` comes from Sage's `cached_method` and is then superseded; the terminal exception is the
`OverflowError` while a `_TransversalCoset` is being rendered for the message. The substantive failure
is the preimage/coset construction for an **intersection** subgroup of an infinite orthogonal group
whose discriminant module has nine invariant factors `2`.

## Why this is not the port's, and not a cap

Two runs, same input, differ only in whether the port's `#67` caching change is present:

| build | result |
| --- | --- |
| hoisted `_enriques_transport()` (port `8c1c9925`) | **1 failed in 141.62 s** |
| genuine pre-hoist copy (`git show 8c1c9925^:…`, 0 occurrences of `_enriques_transport`) | **1 failed in 136.76 s** |

Identical `KeyError` → `OverflowError` in both. The caching is exonerated by this comparison, and the
defect predates it. Raw logs: `.tmp/enriques-case87-coset-overflow.log` (hoisted, `--tb=long`) and
`.scratch/prehoist-case87-real.log` (pre-hoist).

It is also **not** a timeout: 140 s against a 1860 s cap. Whatever the earlier "Enriques is capped"
narrative recorded, this case at this pair returns an error, and the two must not be confused.

## What the port needs from the owner

The contract the port's consumer uses: `G.group_generators()` for `G` an intersection of a
discriminant-preimage subgroup with `O^+(N)` must return generators (or raise a typed refusal naming the
missing capability). The current path instead dies inside `_free_cover`'s recursive
`selected_group_resolution().augmentation()` with an untyped `OverflowError`, which no caller can act on
and which `#33`'s "typed errors" invariant forbids.

A minimal reproducer is the shipped case: `_gamma_h` at `case_87`'s polarization, then
`group.group_generators()`.

## The family is heterogeneous, measured

Four cases sampled at the same pair, each with its own 500 s inner cap:

| case | outcome |
| --- | --- |
| `case_87` | **failed** at ~140 s (the provider `OverflowError` above, full traceback captured) |
| `case_50` | **failed** at 141.31 s with the **same** traceback, captured afterwards with `--tb=long` (see below) |
| `case_1` | exceeded 500 s (killed by the sampling cap) |
| `case_2` | exceeded 500 s (killed by the sampling cap) |

**`case_50` shares case_87's cause — now measured, not inferred.** The first sample used `--tb=line` and
captured only `1 failed in 151.04s`, so the exception was left unattributed at that point. A second run
with `--tb=long` (`.tmp/enriques-case50-coset-overflow.log`) shows the identical signature:
`OverflowError: value too large to convert to int` raised in the `repr` of a `_TransversalCoset`, at
`sage/libs/gap/element.pyx:3166`, `1 failed in 141.31 s`. Same failure, 141 s rather than 151 s.

**Correction, recorded rather than quietly fixed.** The first version of this table said `case_50`
exceeded 500 s. That was inferred from two `Killed` lines in a log whose output was buffered behind a
pipe, before the run finished; the completed log shows `1 failed in 151.04s`. Two cases error quickly
and two outrun the sampling cap — inferred outcomes from partial logs are not measurements, and this
note now only carries what the finished logs show.

So the family contains at least two distinct failure modes: a fast, reproducible provider error, and
cases that need more than 500 s. Those must not be merged into one story: an error is a defect to fix, a
slow case is a cost question, and only a bounded run that finishes can call any case "expensive".

## The published fix, located and compared (source-grounded)

The provider line already contains the repair for this exact failure:

```
ba4efaa  fix: coset positions are read at GAP's index width, not a machine int's
```

* It is **155 commits downstream of this repo's pin** `cdafca84` (`git merge-base --is-ancestor`
  confirms `cdafca84` is an ancestor of `ba4efaa`).
* The fix is **absent at our pin and present at `ba4efaa`**, by source inspection:
  `grep -c '_gap_list_entry\|ELM_LIST' src/sage_extensions/categories/group/groups.py` gives **0** in
  `cdafca84` and **14** in `ba4efaa`.
* Its own diagnosis matches ours exactly: `GapElement_List.__getitem__` binds its index through
  `cdef int`, and the right transversal of `O(A_N)` has index `|O(A_N)| = 46998591897600 > 2**31`, so a
  valid position raises `OverflowError: value too large to convert to int` at
  `sage/libs/gap/element.pyx:3166` — the line in both of our tracebacks. The repair reads transversals
  through GAP's `ELM_LIST`, which takes the index as a GAP integer.

What this did **not** do: it did not change this repo's pin, and it did not edit the provider. The
measurement of whether `ba4efaa` makes the port's failing cases pass is run against a per-process bind
of a tree extracted from that commit, which is how this repo consumes a provider revision in a run.

## Measured against `ba4efaa`: the overflow clears, a second provider error appears

Both cases were re-run at the same port with the provider tree extracted from `ba4efaa` bound in place
of the pin (`--tb=line`, each with its own 650 s inner cap):

| case | at pin `cdafca84` | at `ba4efaa` |
| --- | --- | --- |
| `case_87` | `OverflowError` at ~140 s | **no error**: ran past the 650 s cap |
| `case_50` | `OverflowError` at 141.31 s | **different failure at 130.35 s** |

So the transversal fix does what it says: the `_TransversalCoset` index overflow is gone, and `case_87`
stops failing fast and becomes a long computation instead. That is progress, and it is not acceptance —
`case_87` has no verdict yet, only "longer than 650 s".

`case_50` at `ba4efaa` exposes the **next** provider defect (`.tmp/enriques-case50-ba4efaa-next-error.log`):

```
E   KeyError: (({g in Orthogonal group O(rank 10, (1,9)) : g fixes -11*e_0 - 30*e_1 + ...},), ())
During handling of the above exception, another exception occurred:
E   ValueError: Pcgs([f1, ..., f10]) -> [f1, f1*f2*f5*f7*f9*f10, ...] is not in the group
    Group of isometries of Finite quadratic module over Z with invariants (2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
```

i.e. a pcgs/permutation-representation image that is not recognized as a member of the target isometry
group of the discriminant module. That is a provider-layer membership/representation matter, in the same
family of machinery but a distinct defect from the transversal index.

**Reported to supervision as the exact missing provider acceptance:** at `ba4efaa`, `case_50` needs the
pcgs image membership to resolve; `case_87` needs a bounded run long enough to finish (650 s was not).
Neither is something this repository may patch, and neither is a reason to weaken or skip the case.

## Status kept honest

`case_87` is one of the 87; this note establishes one failing case and its owner, not the family's state.
The family's acceptance remains `#64`'s, and the 5/10/19 counts and `S1` are untouched by anything here.
