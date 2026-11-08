# #70 blocked on one provider operation: the Schreier factorization of the stabilizer family

**Status: no port patch is landed here.** The unit as specified cannot be implemented port-side without
either reaching into provider-private state or changing the group. Both are excluded by AGENTS.md and by
the user's instruction, and both are excluded here for measured reasons, not stylistic ones. The port's
three conversion sites are unchanged.

Pair for every measurement below: port `893cb233` (`src/` and `tests/` identical to `45c22a7c`), provider
pin `cdafca84065405cb810bba047f8a965299b03762`, adopted runtime (`/sage` at `281da26ffd`, Sage
11.0.beta5, GAP 4.16.1).

## 1. What the port does today, exactly

Three call sites convert every element of the provider's Schreier family:

- `src/sage_indefinite_port/indefinite/isotropic_induction.py:298-304` — `_stabilizer_generators`,
  used by `line_automorphism_generators` (`:306`), `plane_automorphism_generators` and
  `flag_automorphism_generators`;
- `src/sage_indefinite_port/indefinite/recursive.py:1053-1057` — `line_stabilizer_generators`;
- `src/sage_indefinite_port/indefinite/recursive.py:1559-1571` — `pointed_automorphism_generators`'s
  isotropic branch.

All three are the same expression:

```python
return tuple(lattice.integral_isometry(generator) for generator in action.stabilizer().group_generators())
```

## 2. Why the family is that size, read from the provider at the pin

`IntegralStructureAction.stabilizer()` (`src/sage_extensions/categories/rational_integral_stabilizers.py:439`)
takes the finite branch and returns `_FiniteIntegralRepresentation.lattice_stabilizer()`
(`:228-243`), which is

```python
self.permutation_representation().preimage_subgroup(
    self.finite_stabilizer(),
    generating_morphism=self.g_set().orbit_schreier_morphism(self._selected_point()),
)
```

so the stabilizer is a `PreimageSubgroups` object whose selected generating morphism is Schreier's lemma on
the lattice orbit. That construction has exactly one implementation,
`_orbit_schreier_morphism` (`src/sage_extensions/categories/group/g_sets.py:60-85`), and it writes one word
per (orbit point, free generator):

```python
words = tuple(
    (~word_section(orbit.act(augmentation(generator), point))) * generator * word_section(point)
    for point in orbit
    for generator in free.group_generators()
)
```

i.e. `(~t_{s·p}) · s · t_p`, in the provider's own docstring wording. The family size is therefore
`|orbit| × |letters|`, an identity of the construction and not a property of the mathematics: on the
apex Enriques line of #18 it is `1024 × 25 = 25600`.

## 3. Why a port-side subset cannot fix it

### 3a. Selecting by finite image is wrong, by measurement

The finite representation sees the *image* of an element; the lift the provider returns is one specific
element of that image's fibre, and the fibre can contain the sign-reversed element. The recorded
measurement on the #18 apex line (sage-extensions `TRAPS.md` at `c4c0469e50`, entry "Issue #18 measured on
the adopted runtime") converts the first 2000 family elements and finds line actions `+1 → 1920`,
**`-1 → 80`**, first `-1` at prefix position 14, with the form preserved by 2000 of 2000. Two elements with
equal finite image can therefore differ by a sign on the marked line. A subset chosen by image (or by
"does it change the image subgroup") would drop sign-reversing generators from the returned group and
silently change `Aut((L, Zu))` — the acceptance cases assert exactly that marking content.

### 3b. The alternatives the provider already exposes do not reduce the family

Probed on the adopted runtime (log: `.tmp/issue70-evidence/probe70-routes.log`, reproduced verbatim in §6a).
The stabilizer object carries `kernel_image_generating_morphism`
(`src/sage_extensions/categories/group/predicate_subgroups.py:848`) and `schreier_generating_morphism`
(`:815`, `g_sets.py:661`). On `2U` with the line `Ze1` the family has 7 elements and both alternatives
return 7 words, the same mod-3 image order (72) as the family. On `U(4) + <-8>` with the line `Ze1` the
family has 6 elements (mod-3 image order 12) and both alternatives raise, inside the provider:

```
TypeError: membership of Morphism from Lattice of rank 3 and signature (1, 2) to Lattice of rank 3 and
signature (1, 2) in Subgroup of Orthogonal group O(Lattice of rank 3 and signature (1, 2)) is undecided;
sup...
```

So neither alternative is a route to fewer generators: where it works it reproduces the family, and where
the family is not the trivial one it refuses. Neither returns a generating *set* distinct from the
enumerated family.

### 3c. The factorization is in provider-private state

The words need `word_section` (`p |-> t_p`) and `augmentation` (`F ->> G`). On the stabilizer object those
live behind `self._coset_word_section` — the name is asserted by the provider itself at
`g_sets.py:663` — and `_selected_schreier_augmentation(words)`. `probe70-names.log` lists the public
surface; the only public accessors are `coset_space()`, `right_cosets()`, `free_cover(...)`,
`selected_group_resolution()` (a category object, not the instance's morphism) and
`group_generators()`. Reading `_coset_word_section` from the port is the private-state access the project
forbids; rebuilding `(~t_{s·p})·s·t_p` at the port from those privates is a duplicate of the provider's
Schreier construction, which `references/port-design.md` assigns to research.

## 4. The measurement that prices the fix

From the same `TRAPS.md` entry, on the #18 apex line at this kind of pair: `isotropic_stabilizer_action`
20.3 s; `stabilizer().group_generators()` 1643.0 s for 25600 words; conversion 0.143 s per element
(286.7 s for the first 2000), i.e. ≈ 3668 s for all of them, ≈ 5311 s total against the 2400 s inner cap.
Converting 1049 factors (25 letters + 1024 transversal values) at the same rate is ≈ 150 s, which fits.
The waste is a factor of ≈ 24 in the conversion stage, and it is exactly the factor `|orbit|`.

## 5. The provider operation this unit needs

I need one operation on the integral-structure action that hands the port the Schreier factorization
without exposing private state: call it `IntegralStructureAction.stabilizer_schreier_factors()`; it takes
no argument beyond the action itself (the action already carries the rational group, its selected
generating set, the lattice orbit and its section), and returns two tuples of owned automorphisms of the
same rational space — `letters`, the images of the free cover's generators in `G`, and `transversals`,
one element `t_p` per point `p` of the lattice orbit in the orbit's own point order, with `t_p` carrying
the base point to `p` — so that the port forms each Schreier word through provider multiplication only;
and it returns as its certificate the guarantee that the composites `(~t_{s·p})·s·t_p` taken in that
order equal `stabilizer().group_generators()` element for element, which the port can verify directly on
any case small enough to materialize both and which is what makes the substitution sound rather than
plausible. With that operation the port converts `|letters| + |orbit|` elements (1049 on #18, ≈ 150 s)
instead of `|letters| × |orbit|` (25600, ≈ 3668 s), composes the family with provider multiplication, and
returns the identical group and the identical marking content; without it, the only port-side change
available is to relabel identity elements, which is not a fix.

The same need is what dzackgarza/sage-extensions#273 solved on the torsion-form side
(`stabilizer_generators(element)`, returning 6 generators instead of enumerating `O(A_T_Co)`, certified by
`|<gens>| = |O|/|orbit|`); this unit needs the integral-structure counterpart, where the certificate cannot
be an order because the stabilizer is infinite, and is instead the elementwise equality with the family the
provider already builds.

## 6. Probe logs, verbatim

### 6a. `.tmp/issue70-evidence/probe70-routes.log`

```
[2U line] action built in 2.7s
[2U line] family size 7 built in 0.0s
[2U line] converted all in 0.2s (7 conversions)
[2U line] mod-3 image order of the full family: 72 (0.0s)
[2U line] kernel_image_generating_morphism: domain generators 7 in 1.6s
[2U line] kernel_image_generating_morphism: evaluated 7 elements in 0.6s
[2U line] kernel_image_generating_morphism: converted 7 in 0.2s
[2U line] kernel_image_generating_morphism: mod-3 image order 72
[2U line] schreier_generating_morphism: domain generators 7 in 0.5s
[2U line] schreier_generating_morphism: evaluated 7 elements in 0.6s
[2U line] schreier_generating_morphism: converted 7 in 0.2s
[2U line] schreier_generating_morphism: mod-3 image order 72
[U(4)+<-8> line] action built in 61.3s
[U(4)+<-8> line] family size 6 built in 0.0s
[U(4)+<-8> line] converted all in 0.2s (6 conversions)
[U(4)+<-8> line] mod-3 image order of the full family: 12 (0.0s)
[U(4)+<-8> line] kernel_image_generating_morphism: TypeError: membership of Morphism from Lattice of rank 3 and signature (1, 2) to Lattice of rank 3 and signature (1, 2) in Subgroup of Orthogonal group O(Lattice of rank 3 and signature (1, 2)) is undecided; sup
[U(4)+<-8> line] schreier_generating_morphism: TypeError: membership of Morphism from Lattice of rank 3 and signature (1, 2) to Lattice of rank 3 and signature (1, 2) in Subgroup of Orthogonal group O(Lattice of rank 3 and signature (1, 2)) is undecided; sup
[peakrss] run: child waited=0 shell_code=0 samples=37
[peakrss] peak tree RSS      : 3630 MiB
[peakrss] low-water MemAvail : 5232 MiB
[peakrss] peak single VmHWM  : 3653 MiB
[peakrss] peak VmHWM by executable (this is the per-process attribution):
[peakrss]   julia                      3192 MiB peak
[peakrss]   python                      452 MiB peak
[peakrss]   timeout                       7 MiB peak
[peakrss]   proot-x86_64                  3 MiB peak
```

### 6b. `.tmp/issue70-evidence/probe70-names.log`

```
relevant names: ['_engine_subgroup_from_generators', '_finite_generation_decision', '_finite_magma_generation_decision', '_left_coset_construction', '_predicate', '_preimage_morphism', '_selected_schreier_augmentation', 'contains_character_kernel', 'coset_space', 'defining_predicate', 'double_cosets', 'finite_image_lifts', 'finite_image_section', 'finite_index_schreier_generating_morphism', 'finite_index_schreier_generators', 'finite_words', 'free_cover', 'generated_submonoid', 'group_generating_set', 'group_generator', 'group_generator_morphism', 'group_generators', 'group_generators_are_computable', 'has_computed_group_generators', 'has_finite_coset_space', 'has_selected_finite_group_generating_set', 'has_selected_group_resolution', 'has_selected_left_cosets', 'has_selected_resolution', 'image_set', 'image_under', 'intersection', 'is_finitely_generated', 'is_finitely_generated_as_magma', 'isotropic_orbit_representatives', 'kernel_image_generating_morphism', 'left_cosets', 'number_of_group_generators', 'predicate', 'predicate_subgroup', 'predicate_submonoid', 'preimage_morphism', 'right_cosets', 'schreier_generating_morphism', 'select_group_generators', 'select_group_resolution', 'select_left_cosets', 'selected_group_resolution', 'selected_resolution', 'selected_resolution_generating_set', 'selected_resolution_generator', 'vector_orbit_representatives']
  _selected_data_registry: CachedMethodCallerNoArgs -> {Category of groups: <sage_extensions.categories.abstract_categories.objects._PendingResolution obje
  _selected_schreier_augmentation: raised TypeError: PredicateSubgroups.ParentMethods._selected_schreier_augmentation() missing 1 required positional argument: 'words'
  has_selected_finite_group_generating_set: method -> True
  has_selected_group_resolution: method -> True
  has_selected_left_cosets: raised TypeError: OwnedGroups.ParentMethods.has_selected_left_cosets() missing 1 required positional argument: 'subgroup'
  has_selected_resolution: raised TypeError: Objects.ParentMethods.has_selected_resolution() missing 1 required positional argument: 'owner'
  selected_group_resolution: method -> <sage_extensions.categories.abstract_categories.resolutions.Resolutions.parent_class object at 0x7d3
  selected_resolution: raised TypeError: Objects.ParentMethods.selected_resolution() missing 1 required positional argument: 'owner'
  selected_resolution_generating_set: raised TypeError: Objects.ParentMethods.selected_resolution_generating_set() missing 1 required positional argument: 'owner'
  selected_resolution_generator: raised TypeError: Objects.ParentMethods.selected_resolution_generator() missing 2 required positional arguments: 'owner' and 'label'
[peakrss] run: child waited=0 shell_code=0 samples=4
[peakrss] peak tree RSS      : 419 MiB
[peakrss] low-water MemAvail : 3641 MiB
[peakrss] peak single VmHWM  : 419 MiB
[peakrss] peak VmHWM by executable (this is the per-process attribution):
[peakrss]   python                      408 MiB peak
[peakrss]   timeout                       7 MiB peak
[peakrss]   proot-x86_64                  3 MiB peak
```

## 7. What stays true meanwhile

The three call sites are correct; they are simply priced by the family's size. No frozen case, fixture or
assertion was touched, and the families that use these routes keep their recorded results (see the
acceptance rerun of 2026-11-06 in `.tmp/acceptance-coverage.md`). #70 stays open until the operation above
exists and the port's conversion count is bounded by the factors rather than by the family.


## Current original O_plus_N failure

At immutable port 2e36cff/provider 7960e1ae, the unchanged original line/plane assertion is terminal under its original 2400-second bound. The final stack is `_integral_stabilizer_family` yielding `inverse_transversals[image_position] * letter * realized_transversals[position]`, through owned isometry composition and formed-module-arrow construction. Its process group is empty. The complete requested assertion remains open.

The supplied finite action already retains every Schreier factor and its exact generator permutation. The next repair lowers those retained factors into the same rational framing once, forms every exact Schreier matrix product there and applies the existing integral Nielsen reduction. A product is admitted over ZZ only after the rational factors compose. Each surviving full integral action is raised in the actual lattice isometry homset. The complete family, including integral kernels invisible in finite quotients, must survive as the same generated group.

Verification must compare every factored full forward/inverse action with the owned source family on small noncommuting and nonintegral-factor examples, then establish the unchanged original O_plus_N assertion. A reduced family size or isolated matrix timing cannot establish that assertion.
