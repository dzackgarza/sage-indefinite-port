# Owner handoff: the Enriques `case` bottleneck — exactly what is measured, and what is not

Provider of record: `dzackgarza/sage-extensions` at `42d7361cfda3d8505f365b26d9aa687447d6df9d`
(`~/sage-extensions-42d7361c`). Port of record: `sage-indefinite-port` at src
`2cdac947`-tree / `a74fb2ca`-tree (identical `src/`), delivered as
`~/sage-image/node_run.py` and `~/sage-image/probe_enriques_stages.py`.

**Nothing in `sage-extensions` was patched, and no ambient order or group was substituted.** The acceptance
assertions are unchanged and were not run against any replacement object.

## 1. The obligation, unchanged

`tests/acceptance/test_boundaries.sage:219-230`
(`test_enriques_modular_groups_match_dutour_sikiric_hulek`), one node per case, 87 cases from
`tests/fixtures/enriques_87_polarizations.json` (source: Dutour Sikirić–Hulek, arXiv:2302.01679,
`table_subgroups1`).

Input, per case: the **face polarization** `h` in the coordinates of `U + E8(-1)` (the fixture's
`face_polarization`, 10 integers), transported to the invariant lattice

```
M = U(2) + E8(-2)        (the Enriques-involution invariant lattice of the K3 lattice)
```

Output, per case, two claims and no others:

1. `int(image.cardinality()) == case["group_order"]`, where `image = target.subgroup_on(tuple(g.discriminant_morphism() for g in group.group_generators()))`
   and `target = N.discriminant_group().orthogonal_group()`, `N` the **anti-invariant** lattice of the
   involution (rank 12, signature (2, 10)) and `group = Γ_h = O^+(N) ∩ π_N^{-1}(π_M(O(M, h)))`.
2. `_counts(group) == (line_orbits, plane_orbits, flag_orbits)` — `Γ_h`-orbits of isotropic lines, isotropic
   planes and line-in-plane flags of `N` (`_counts`, `test_boundaries.sage:35-41`).

Two specimen rows, both from the fixture, both with the smallest published cost proxy
(`(L + 10P + 100F) · log2 |Γ_h|`) among the 87:

| case | `group_order` (published) | line / plane / flag orbits | signature of `N` |
| --- | --- | --- | --- |
| `case_3` | `92 897 280` | 5 / 10 / 19 | (2, 10) |
| `case_1` | `696 729 600` | 5 / 9 / 18 | (2, 10) |

`case_1` is the case that previously "blocked at 3 600 s"; at the current pins it is **still capped**,
measured this session: `test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]`, cap 1 800 s,
guard 6 144 MiB → **CAPPED**, killed at 1 800 s, peak 4 182 MiB, 30 stack dumps landing in **26 distinct
innermost frames** (diffuse; not one hot call).

## 2. Measured stage timing — established

`probe_enriques_stages.py` reproduces `_gamma_h` (`test_boundaries.sage:200-215`) and stage 1 of the
assertion, printing before each stage. One bounded run, provider `42d7361c`, cap 2 400 s, guard 6 144 MiB,
job `jmvz5p00r`:

| stage | time | scope |
| --- | --- | --- |
| fixtures + `sage_extensions.catalogue.Involutions.I_En.primitive_extension().invariant` | **4.3 s** | case-independent |
| `lattice(E10).twist(2).Isom(invariant).an_element()` | **19.5 s** | case-independent; note `2r` in the `.sage` file preparses to the integer `2` |
| `involution.polarized(h).coinvariant_extension_subgroup()` | **249 s** | per case |
| `∩ N.O_plus()` | 0.1 s | per case |
| **bundled stage**: `N.discriminant_group()`, `.orthogonal_group()`, `discriminant_morphism()` per `Γ_h` generator | **> 2 127 s and unfinished** — the probe ended at its own 2 400 s cap, still inside this stage | per case |

### 2b. The bundle is **not** the `O(q_N)` group — measured by a separate acceptance node

`tests/acceptance/test_isotropic_orbits.sage::test_o_n_maps_onto_o_of_the_discriminant_form` asserts
`n.discriminant_group().orthogonal_group().cardinality() == 46 998 591 897 600` (and the same order for
`n.O().discriminant_image()`), for the unpolarized Enriques lattice of the same rank 12, signature (2, 10)
shape. At the same pins it **passes in 468 s, peak 3 868 MiB** (job `jmvz75d9x`), and **6 of its 7 stack
dumps are inside `n.O().discriminant_image()`**, not inside `orthogonal_group()`.

So for this form: `discriminant_group()`, `orthogonal_group()` and its `.cardinality()` are **cheap**; the
~400 s of that node is the discriminant image of the infinite `O(N)`. Consequently the > 2 127 s unfinished in
the Enriques bundle is **not** `O(q_N)` materialization, and the earlier "materialize the whole group" framing
is dead twice over (first the laziness, now this).

### 2c. The bottleneck is **generating `Γ_h` itself** — measured, bounded, 900 s

A second bounded probe (`probe_enriques_morphism_cost.py`, cap 900 s, guard 6 144 MiB, job `jmvz7ma89`) builds
`case_3` and then reads the generators, timing every step:

| step | time |
| --- | --- |
| involution, invariant, transported isometry | 13.6 s |
| `involution.polarized(h).coinvariant_extension_subgroup()` | → 218.9 s |
| `N` = rank 12, signature (2, 10); `group = … ∩ N.O_plus()` | → 218.9 s (the intersection itself is 0.1 s) |
| **`group.group_generators()`** | **> 681 s, unfinished at the 900 s cap** — the generator count was never reached |
| `discriminant_morphism()` per generator, `cardinality()`, orbit counts | **never reached** |

So the > 2 127 s of the longer probe is accounted for at its *first* sub-step: computing a **generating set of
the intersection `Γ_h = coinvariant_extension_subgroup(h) ∩ O^+(N)`**. That is the thing to make cheaper or
lazy, and it is a different object from `O(q_N)`, from the count assertion, and from the orbit counts. The
per-generator `discriminant_morphism()` cost remains **unmeasured** because the generator list never arrived.

So the per-case cost, in measured order: **~219 s** to build the coinvariant subgroup, then **> 681 s and
unfinished** to produce a generating set of `Γ_h`, then an unmeasured number of `discriminant_morphism()` calls,
then the cardinality, then the three orbit counts. The longer probe's bundle exceeded **2 127 s** unfinished
(peak 4 943 MiB tree: julia 3 177, python 1 754). Everything before the orbit counts already exceeds the
1 800 s cap, which is why raising the cap would not help.

## 3. What is *established* about the bundle, and what is *not*

The bundle is three calls and the probe prints only around all three, so the split inside it is
**unmeasured**. In particular the following is established from source:

- `orthogonal_group()` on the descriptor is a delegating accessor:
  `form_modules.py:3759` → `self.automorphism_group()`, a `@cached_method` at `form_modules.py:3755` →
  `_torsion_form_automorphism_group(form, quadratic=True)` (`torsion_form_modules.py:1606-1612`).
- That constructor **does not materialize the group**: `TorsionFormOrthogonalGroup.__init__`
  (`torsion_form_modules.py:1208-1235`) ends with
  `self.select_group_generators(self._computed_group_generators)`, and
  `select_group_generators` (`groups.py:5531-5548`) is documented and implemented as **lazy**: "selecting
  computes nothing, and the function is called when `group_generators()` or the selected resolution is first
  read."
- What *would* materialize it: `_computed_group_generators` (`torsion_form_modules.py:1237-1242`) is
  `finite_ordered_set(tuple(self._from_engine(g) for g in self._engine_group().gens()))`, and
  `_engine_group` (`:1288-1291`) → `_engine_orthogonal_group` (`:1752-…`), which builds `O(q_p)` per primary
  part — GAP `GO(ε, n, p)` for elementary odd parts, OSCAR's normal-form algorithm otherwise — and crosses
  every generator to a live automorphism.
- Subgroup construction from ambient elements requires ambient-parented elements:
  `subgroup_on` (`torsion_form_modules.py:1523-1536`) rejects any generator whose `parent() is not self`,
  then calls `self.subgroup(supplied)`.

**Correction to an earlier note of mine in this repo.** I wrote in `.tmp/acceptance-coverage.md` that
`TorsionFormOrthogonalGroup.__init__` eagerly materializes `O(q_N)` and that this was the > 700 s term. That
was wrong: `select_group_generators` is lazy, so the > 700 s is somewhere in the bundle {`discriminant_group()`,
`orthogonal_group()` (lazy effect only), `discriminant_morphism()` × generators}, and the whole-group
generator family is only forced when something reads `group_generators()`/the resolution — most plausibly the
`subgroup_on` call, which the probe had not reached. Recorded here rather than left standing.

What is **not** established, and should not be assumed by the owner:

- which of the three calls holds the ≥ 1 090 s;
- whether `discriminant_morphism()` per `Γ_h` generator is cheap, expensive, or itself forces the ambient
  engine group;
- whether a route exists in the maintained contract that decides `|image|` without the ambient generating
  family. I have **not** shown such an API is missing. `_engine_subgroup_from_generators`
  (`torsion_form_modules.py:1309-1313`) builds a subgroup from generators through `self._engine_group()`, but
  it is a method *on* the ambient group object and still needs that engine group; whether the engine group
  parent (as opposed to its `.gens()`) is the expensive part is unmeasured.

## 4. What the owner is asked to consider (not a claim that it is possible)

Decide `|π(Γ_h)| ≤ O(q_N)` for a subgroup **given by its generators in `O^+(N)`** without materializing the
ambient `O(q_N)` generating family, if the maintained contract permits it. The discriminating
measurement is three separate timings, which needs provider-side instrumentation:

```
T1 = N.discriminant_group()                                  # the torsion quadratic module on the Smith engine
T2 = T1.orthogonal_group()                                   # expect ~0 by the laziness above
T3 = [g.discriminant_morphism() for g in Γ_h.group_generators()]
T4 = T2.subgroup_on(T3);  int(T4.cardinality())
```

Whichever of `T1`…`T4` holds the time decides whether this is (a) a construction-cost matter inside this
specimen's `q_N`, (b) a lazy-vs-eager question in the ambient group, or (c) a genuinely missing subgroup-order
route.

Note, explicitly **not** proposed as a substitute: `discriminant_group().subgroup_on(...)`
(`discriminant_modules.py:109-111` → `_discriminant_subgroup`, `:624-644`) builds a **submodule of `A_N`** on
the unformed module's Smith engine. That is a different object from a subgroup of `O(q_N)` and is not
interchangeable with the assertion; it is listed only because its engine route is cheap and the owner may know
whether an analogous `O(q_N)`-side route is available.

## 5. Why this is an owner item and not a port patch

Nothing in the port is on the measured path here: the port's contribution is `Γ_h`'s generators
(`automorphism_group_generators`) plus the three orbit counts; the ≥ 1 090 s is inside the provider's
discriminant-form machinery. Project rule: no local patches to `sage-extensions`. The port-side obligations
that remain are unaffected and are tracked separately in `.tmp/acceptance-coverage.md`.

Secondary, and cheap for whoever touches the acceptance harness: **23.8 s of every case is case-independent**
(the involution, its primitive extension, and the transported isometry — rows 1 and 2 above), recomputed 87
times = ~34 min across the family. That is 3 % of one case's cost, so it is not the reason any case caps.
