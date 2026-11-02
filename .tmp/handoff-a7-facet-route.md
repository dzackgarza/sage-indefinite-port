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
| row rank | **10 = ambient dimension** — the cone is full-dimensional, so there is no lineality to declare |
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
- **What is not available** — `DualDescriptionRecordFullDim` and `RecordDualDescOperation` are declared in the
  vendored header but **defined outside this tree** (`rddo` is constructed at `:1542`). Their algorithm cannot
  be compared line by line from what is vendored here. What *is* documented by the call signature is that the
  routine consumes the stabilizer permutation group, so its contract is orbit-reduced.

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
| `subspace=` (declare lineality) | **inapplicable**: the cell's row rank is 10 = ambient dimension, so the cone is full-dimensional and there is no lineality to declare (measured). |
| `congruences` / `excluded_faces` / `open_facets` | recognised, but none of them reduces a full-dimensional double description; not pursued as a cost lever. |
| PPL `from_rays` | the reason the Normaliz branch exists; documented in the port as the earlier bottleneck; not re-measured at the current pins. |
| **source route: symmetry-reduced dual description, one facet per stabilizer orbit** | **not implemented in the port.** The ingredients exist (`cell_stabilizer`, `_kernel_flipping` returning the flipped normal, `flip_across`), and the traversal already retains facet-orbit representatives and Schreier words — but `cell_cone` does not use them; it computes the full H-description and discards the non-representatives. |

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
