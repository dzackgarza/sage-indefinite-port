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

## Status kept honest

`case_87` is one of the 87; this note establishes one failing case and its owner, not the family's state.
The family's acceptance remains `#64`'s, and the 5/10/19 counts and `S1` are untouched by anything here.
