# A7 `a7_plus_diagonal`: the facet requirement, the reference route, and what this image's engines offer

Status: **open**. Nothing in `sage-extensions` was patched. The acceptance assertion is unchanged and no limit
was raised. This document replaces the earlier "no port-side lever remains" note, which was inferred from a
missing keyword and is wrong.

## 1. The exact requirement

`LorentzianPerfectLocalBackend.cell_cone(cell)` (`src/sage_indefinite_port/indefinite/lorentzian_cells.py:1011`)
returns the owned polyhedral cone of a perfect cell, from which the traversal reads

- `halfspace_covectors()` — the cell's facet (support) covectors (`:1121`, `:1129`, `:1142`, `:1152`, `:2830`),
- `primitive_rays()` — the cell's extreme rays, used for isotropic marks (`:1166`).

It is called for **every cell** in the quotient traversal (`:2830`), and the result is retained per cell
(`cell_cones[id(cell)]`). Above 10 000 rows the cone is built from Normaliz's raw H-data
(`_normaliz_cell_h_presentation`, `:1032-1067`); below, from `RationalPolyhedralCones(...).from_rays(...)`
(`:1013`).

Downstream, the port already knows it only needs **one facet per stabilizer orbit** — its own comment at
`:2834-2840` says so ("The upstream Lorentzian perfect algorithm computes one facet per stabilizer orbit
(lorentzian_perfect.h::f_adj); the Python traversal has already retained those exact orbit representatives and
words"). It nevertheless materializes the *complete* H-description and then keeps only representatives.

## 2. The measured specimen and its measured cost

From `load_oscar_lattice_oracles()["discriminant_images"]`, id `A7_plus_diag_1_1_-1`:

| quantity | measured |
| --- | --- |
| Gram | `A7 ⊕ <1> ⊕ <1> ⊕ <-1>` — rank 10, signature (9, 1) |
| `initial_cell(...)` | **49 861 rows**, row length 10, built in **4.5 s** |
| row rank | **10 = ambient dimension** — the cone is full-dimensional. Pointedness does **not** follow from this (see the correction in §5) and is certified separately |
| max coordinate | **15** (small entries; not an arithmetic-size problem) |

Node run at src `fad2c5bb`-tree, provider `42d7361c`, cap 1 800 s, guard 6 144 MiB (job `jmvzblu48`):
**CAPPED**, peak 5 405 MiB (below the guard ⇒ time-bound), and **all 29 stack dumps inside one call**:

```python
hyperplanes = NmzResult(normaliz_cone, "SupportHyperplanes")     # lorentzian_cells.py:1051
```

So the cost is Normaliz's facet enumeration of a full-dimensional cone with 49 861 generators in dimension 10.

## 3. What the reference does for this same requirement

The reference never builds an H-description of the 49 861-vector configuration, and never calls Normaliz
(the only "normaliz" strings in `polyhedral_common@1592b246` are `PositivityNormalizeGcdDot` &c.).

- **Initial cell** — `LORENTZ_GetOnePerfect` (`src_lorentzian/lorentzian_perfect.h:969-1100`): start from a
  critical set (`LORENTZ_SearchInitialVector`), then a **rank-growing loop**: while `rank(CritSet) < n`, take
  the nullspace of the expanded critical set, pick an outside ray (`LORENTZ_GetOneOutsideRay`), and flip
  (`LORENTZ_Kernel_Flipping`). It returns `{CritSet, eNSPbas, CentralVectT}` — critical rays plus **one facet
  normal**, never a facet list. At most `n` steps.
- **Facets and adjacency** — `f_adj` (`:1433-1472`): compute the cell's stabilizer
  (`LORENTZ_ComputeStabilizer`), then

  ```cpp
  vectface TheOutput = DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, data.rddo);
  ```

  i.e. a dual description **taken modulo the stabilizer group** (`GRPperm`), giving one facet per orbit; each
  orbit representative is then flipped (`LORENTZ_DoFlipping`) to obtain the neighbouring cell.
- **Where the primary definition lives** — `DualDescriptionRecordFullDim` and `RecordDualDescOperation` are
  declared in the vendored header but **defined outside this tree** (`rddo` is constructed at `:1542`). The
  upstream repository's own layout names the home: `src_dualdesc` ("for computing dual description on serial
  computers"), and the routine is the **adjacency decomposition** method, documented by the author (BIRS 21w5011,
  "permutalib/polyhedral: Advanced methods for lattice and group computations"):

  > Input: the vertex-set of a polytope P and a group G acting on P. Output: O, the orbits of facets of P.
  > Compute some initial facet F (by linear programming) and insert the corresponding orbit into O as undone.
  > For every undone orbit O of facet: take a representative F of O; find the ridges contained in F, i.e. the
  > facets of the facet F (this is a dual description computation); for every ridge R find the corresponding
  > adjacent facet F0 such that R = F ∩ F0; for every adjacent facet found test if the corresponding orbit is
  > already present in O; if no, insert it as undone. Mark the orbit O as done. Terminate when all orbits are
  > done.

  A maintained Python binding of the same routine exists upstream: `py_polyhedral.dual_description(EXT, GRP)`.
  GAP counterparts are vendored here: `references/vendor/Indefinite.jl@374a5ebb/indef/lib/DualDescriptionSimple.g`
  and `DualDescriptionAndADM.g`.

  So the algorithm to reproduce is documented, and its two moving parts are exactly the ones the port already
  owns: the initial facet normal (the port's `initial_cell` computes `base_normal` by the source's bisection)
  and the per-ridge flip (the port's `_kernel_flipping` / `flip_across`). What the port lacks is the
  orbit bookkeeping *between* them — the adjacency decomposition itself — not the geometry.

## 4. What this image's engines document, measured

PPL route: `RationalPolyhedralCones(lattice).from_rays(rays)` — the port's own comment at `:1015` records it as
the previous A7 bottleneck ("PPL's generic V-to-H conversion is the current A7 bottleneck: its initial perfect
configuration has 49,861 rays"), which is why the Normaliz branch exists. Not re-measured this session.

Normaliz route: the recognised input types in **this** build were enumerated by offering each keyword to
`NmzCone` (`probe_normaliz_options.py`):

| keyword | this build |
| --- | --- |
| `vertices`, `cone`, `subspace`, `congruences` | **accepted** |
| `inequalities`, `equations`, `excluded_faces`, `open_facets`, `grading`, `offset` | **recognised** (rejected only for my deliberately mismatched dimensions / semantics) |
| `number_field` | recognised (needs a 3-entry list) |
| `interior`, `internal_point`, `incidence`, `weights` | **`Unknown type` — not implemented** |

## 5. Bounded alternatives investigated, and their outcome

| alternative | outcome |
| --- | --- |
| `interior` / `internal_point` point | **not implemented in this build** (measured, above). This is what the earlier note wrongly generalised into "no port-side lever". |
| `subspace=` (declare lineality) | **inapplicable, but not for the reason I first gave.** I wrote "full-dimensional, so no lineality"; that inference is invalid — `R = {(1,0),(-1,0),(0,1)}` in `Q^2` has rank `2 = n` and lineality `{(x,0)}`. What is needed is pointedness, and for these cells it is certified a-posteriori by `p := Σ_F a_F` over the facet orbits (the normals span `V*`, so no ray lies in every facet, so `<p,r> > 0` on every ray and on `C \ {0}`) — measured in `tests/test_facet_orbit_contracts.py`. Note that the port's `_positive_direction` is *not* such a certificate: it is the source's `INDEFINITE_GetShortPositiveVector`, the cell's central vector, an element of positive norm, and an interior point of a cone says nothing about its lineality. `subspace=()` is correct **given** that certificate, which is now asserted in `tests/test_facet_orbit_contracts.py`. |
| `congruences` / `excluded_faces` / `open_facets` | recognised, but none of them reduces a full-dimensional double description; not pursued as a cost lever. |
| PPL `from_rays` | the reason the Normaliz branch exists; documented in the port as the earlier bottleneck; not re-measured at the current pins. |
| **source route: symmetry-reduced dual description, one facet per stabilizer orbit** | **not implemented in the port.** The ingredients exist (`cell_stabilizer`, `_kernel_flipping` returning the flipped normal, `flip_across`), and the traversal already retains facet-orbit representatives and Schreier words — but `cell_cone` does not use them; it computes the full H-description and discards the non-representatives. |

### 5b. Measured on cells both engines can handle — and the extrapolation I must withdraw

Two specimens, built the way the traversal builds them (`initial_cell`):

| lattice | rays | facets | ray orbits under the stabilizer | stabilizer generators | `cell_cone` | `from_rays` |
| --- | --- | --- | --- | --- | --- | --- |
| `A2 + <1> + <-1>`, rank 4, signature (3,1) | 27 | 6 | **4** | 4 | 0.11 s | 0.04 s (6, agree) |
| `A3 + <1> + <-1>`, rank 5, signature (4,1) | 89 | 24 | **6** | 5 | 0.27 s | 0.26 s (24, agree) |

**Withdrawn.** I wrote that the facet/ray ratio of 0.75–1.0 here implies the A7 cell's H-description has "an
order of magnitude of 50 000 covectors". Two points at ranks 4 and 5 do not establish anything at rank 10, and
neither the ratio nor the A7 facet count has been measured — the A7 count is precisely what the blocked call
would produce. The measured facts are only those in the table above.

What the table *does* support is the size of the symmetry reduction: 89 rays collapse to **6** orbits, 27 to
**4**. That is the quantity the reference exploits (`f_adj` keeps one facet per orbit), and it is the reason the
orbit route is worth implementing. It is not evidence about the facet count.

### 5c. Why the reduction cannot be done by shrinking the input

It is tempting to feed the engines one ray per orbit instead of all rays. That is wrong in general: with
`R = {e1, e2}` in `R^2` and the swap acting, the ray orbits are `{e1}, {e2}` and `cone(e1) != cone(e1, e2)`. A
cone generated by orbit representatives is a subcone, not the cone. The reduction has to happen **inside** the
dual description, where the processed ray set is kept stabilizer-invariant (rays inserted orbit by orbit), so
that pruning the facet list to orbit representatives is sound at every step. That is the structure the reference
has and the port does not.

### 5d. The contracts, now pinned

`tests/test_facet_orbit_contracts.py` (4 passed, 3.9 s, peak 429 MiB) asserts, on both specimens, the four
properties an orbit-reduced enumeration depends on: a facet covector annihilates exactly a rank-`n-1` incidence
and is non-negative on every other ray; each stabilizer generator permutes the ray set; each sends facets to
facets; every facet lies in exactly one orbit. Any implementation of the orbit route has to keep passing these.

### 5e. Corrections to the certificate layer, from review

Three of these are mine, and two of them I had asserted as fact:

1. **A nonzero minor modulo a prime *does* prove a rational rank lower bound.** `det != 0` in `Z/pZ` implies
   `det != 0` in `Z`, hence in `Q`. I wrote the opposite ("a nonzero reduction never proves a lower bound"),
   which is wrong. What modular arithmetic cannot do is certify an *upper* bound, since `det == 0 mod p` says
   nothing about `det` over `Z`. The port's certificate layer now uses Sage's exact `rank()` for the decision
   and a determinant for the certificate, so the question no longer arises.
2. **`max_entry^n * sqrt(n!)` is not a Hadamard bound.** Hadamard's inequality gives
   `|det| <= prod_i ||row_i||_2 <= (sqrt(n) * max_entry)^n`, and `sqrt(n!)` is far smaller than `n^{n/2}`. The
   bound was only ever used to size a prime; with a too-small prime the modular rank can be artificially
   lowered, which would *miss* a genuine facet or ridge. That whole path is gone: the module composes Sage's
   `rank`/`pivots`/`right_kernel`/`det` and keeps no private elimination, no modular reduction and no primality
   test.
3. **A sample of random examples cannot certify enumeration completeness.** Agreement on 120 random matrices is
   a check of a primitive, not evidence that a facet enumeration is complete — a missed rank silently drops a
   genuine facet or ridge. Completeness comes only from the closure argument over the facet graph.
4. **`orient_facet_normal` accepted all-zero pairings** while its docstring said it refused them: every row
   giving `<normal, r> == 0` satisfies "non-negative on all rays". It now raises, and a test pins both that and
   the impossible-normal case.

### 5f. A7 pointedness, proved on the specimen instead of inferred

From the constructor hypothesis, evaluated on the real configuration in 0.26 s: `rank[1|A] = 10` for the
49 861 × 10 row set, its kernel is spanned by `w = (w_0, w')`, and here
`w' = (1, 0, 0, 0, 0, 0, 0, 0, 0, 1)` takes the value **1 on every one of the 49 861 rows** (min = max = 1, one
distinct value). So the rows lie in the single affine hyperplane `x_1 + x_10 = -w_0`, which does not contain the
origin, the cone is the cone over a 9-polytope, and it is pointed and full-dimensional. That is a proof for this
cell and, by the same argument, for every cell the port builds — not an extrapolation from the two bounded
specimens, and not the false inference "row rank = ambient dimension, so no lineality".

## 6. Established vs unmeasured

Established: the specimen and its dimensions; that the cost is one `SupportHyperplanes` call; that the cone is
full-dimensional with small entries; which Normaliz input types exist in this build; that the reference's route
is orbit-reduced and structurally different; that the port's own comment states one facet per orbit suffices.

Unmeasured: how long the source's orbit-reduced route would take on this specimen; whether an orbit-wise
implementation in the port would fit 1 800 s; whether PPL at the current pins is faster or slower than the
documented stall; what `DualDescriptionRecordFullDim` actually does internally.

## 7. The port-side step this points to

Replace the per-cell *complete* H-description with the source's **orbit-wise** enumeration:

1. compute the cell's stabilizer orbit decomposition of facets (the traversal's facet-orbit data already
   exists at `:2815-2860`);
2. for each orbit representative, obtain its normal by the source's route (`SearchInitialVector` normal, then
   `_kernel_flipping` per representative) rather than by `SupportHyperplanes` on all 49 861 rows;
3. validate against the current cone on cells small enough for both — `tests/test_large_cell_cone_route.py`
   now provides exactly that comparison (route equality and exact ray recovery), so the change is checkable
   without the A7 node;
4. only then re-run A7, once, as a verification run rather than as an exploration.

No provider patch is required for any of this, and it does not change the assertion, the guard or the cap.
