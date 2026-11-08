# #67: the case-independent Enriques setup is computed once per run (measured)

`tests/acceptance/test_boundaries.sage` calls `_gamma_h(face_polarization)` once per Enriques case, 87
times per family run, and again from `test_classical_polarization_indices`. Before the change each call
rebuilt three case-independent objects: `Involutions.I_En`, its primitive extension's invariant lattice,
and `lattice(E10).twist(2).Isom(invariant).an_element()`. The handoff's instrumentation timed the last
two at 4.3 s and 19.5 s at provider `42d7361c`.

After the change `_enriques_transport()` owns them under `@cache` and `_gamma_h` reads it.

## Measured, not inferred

A scratch copy of the module (`.tmp/probe67-enriques-cache.sage`; the suite file itself is unchanged
apart from the hoist) ran one probe test under the same plugin the acceptance suite uses:

```
python -m pytest .tmp/probe67-enriques-cache.sage -p qc_sage_session \
  -k probe_enriques_setup_is_cached
-> 1 passed, 119 deselected in 14.92 s
```

The probe asserts object identity, not just equality, across two calls:

```
first[0] is second[0]     # the involution
first[1] is second[1]     # the invariant lattice
first[2] is second[2]     # the transport isometry
```

All three held. So the cold setup ran once in that process (~15 s end to end, provider `cdafca84`,
which is faster than the 23.8 s recorded at `42d7361c`), and every later `_gamma_h` call in the same
run reuses the same objects.

## What this does and does not establish

It establishes the unit's acceptance sentence for the caching itself: the involution, its primitive
extension and the transported isometry are computed once per acceptance run instead of once per case,
and the objects handed to the polarization stage are the same ones every time, so no assertion can see
a different setup than before the change.

It does not establish the family: whether the 87 cases pass is `#64`'s obligation, and the Enriques
family remains the capped case it was. Collection is unchanged at 119 tests, and the assertion text in
both callers is untouched.
