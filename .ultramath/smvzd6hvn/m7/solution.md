# Orbit-representative facet enumeration for perfect-cell cones — exact certificate layer

Swarm `smvzd6hvn`, member m7, direction d6. Files: `code/primitives.py`, `code/bench.py`,
`results/bench.json`, `measurements.md`.

## 0. Scope, and how this fits the reference route

The enumeration itself is the reference's orbit-wise adjacency decomposition
(`references/vendor/polyhedral_common@1592b246/src_lorentzian/lorentzian_perfect.h:1433-1472`,
`f_adj` → `DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, rddo)`): one facet per
stabilizer orbit, the ridges of a representative facet, the adjacent facet across each ridge,
orbit saturation under exact transport, termination when every orbit is done.

This document does not re-derive that control flow (directions d1 and d4 own it). It delivers
what that control flow *consumes* and cannot work without: three exact primitives, each with a
machine-checkable certificate, plus the two lemmas that make them cheap, plus a cost model that
is **measured** at the specimen's shape rather than asserted:

| primitive | question | answer |
| --- | --- | --- |
| (P1) facet certificate | is `y` the facet normal of this incidence set? | `is_facet_certificate` |
| (P2) adjacency | do two facet normals share a ridge? | `adjacent_from_zero_sets` |
| (P3) normal from incidence | which normal belongs to this incidence set? | `primitive_normal_from_incidence` |

Everything is integer arithmetic. No LP solver, no perturbation, no tolerance, no probabilistic
rank test. Measured: 499 facets over five cones (including a 10-dimensional cone of the
specimen's dimension with 256-ray incidences) are recovered exactly and certified against
Sage's own facet enumeration; the rank certificate agrees with exact `QQ` rank in 120/120 trials;
one (P2) decision at 49 861 rays costs 181–208 µs. See `measurements.md`.

## 1. The construction, stated exactly

### 1.1 Data and conventions

* `L` a lattice of rank `n` and signature `(n-1,1)`, Gram matrix `G`; `cell.rows = (a_1,…,a_m)`
  the cell's coordinate rows in `lattice.module_generating_set()` order, each `a_i ∈ Z^n`.
  This tuple **is** the `EXT` input of the reference route: rows are lattice vectors, one ray of
  the cone per row. On the target specimen `m = 49 861`, `n = 10`, `max |a_i(j)| = 15`.
* `C = cone(a_1,…,a_m) ⊆ R^n`.
* `G_cell ≤ O(L)` the finite cell stabilizer (`cell_stabilizer`, `lorentzian_cells.py:1228`),
  given by generators; each generator acts on rows by a position permutation `π`,
  `g(a_i) = a_{π(i)}` (`_facet_mask_data` builds these packed permutations, `:1330`).
* A facet record is `(y, I)`: `y ∈ Z^n` primitive, `I ⊆ {1..m}` with `y·a_i = 0 ⇔ i ∈ I` and
  `y·a_i > 0` for `i ∉ I`. `I` is the incidence set; the port encodes it compactly as
  `FacetMaskKey` (`:105`) and decodes it as `FacetIncidence` (`:123`). Note the port's key
  stores the *non-incident* side by default (`(True, positions)` = non-incident, `(False,
  positions)` = incident, bare `int` = non-incident bitset); (P1)/(P2) below work on the
  *incident* bitset, so the adapter is one complement against `(1 << m) - 1`.

### 1.2 The cone is a cone over a polytope (and therefore pointed and full-dimensional)

`LorentzianPerfectCell.__post_init__` (`:602`) certifies `rank[1 | A] = n`. Hence the row space
`{1} × R^n` meets the row set in a unique (up to scale) functional `w = (w_0, w')` with
`w_0 + w'·a_i = 0` for every `i`, and `w_0 ≠ 0` because `rank A = n`. So **every cell's rows lie
in one affine hyperplane** `{x : w'·x = -w_0}`, and

> **Lemma 1.** `C` is the cone over the polytope `P = conv(a_1,…,a_m)`, which is
> `(n-1)`-dimensional; `C` is pointed and full-dimensional; and `F ↦ F ∩ P` is a bijection
> between the facets of `C` and the facets of `P`.

*Proof.* Every `x ∈ C` is `Σ λ_i a_i`, and `w'·x = -w_0 Σ λ_i`; scaling with `κ := Σ λ_i` gives
`C ∩ {w'·x = -w_0} = P`, so faces correspond bijectively (every face of `C` is `C ∩ (face of P)`
because both are cut by linear functionals; a `k`-face of `C` meets that hyperplane in a
`(k-1)`-face of `P`). `C ⊆ {w'·x ≥ 0}` or `≤ 0` by sign of `-w_0`, so `C` contains no line, and
`dim C = 1 + dim P = n`. ∎

This is why no hypothesis about full-dimensionality or lineality has to be *checked* for a
perfect cell: it is a consequence of the cell invariants. It also explains the specimen's numbers
— 49 861 rays in rank 10 means a 9-polytope with 49 861 vertices.

### 1.3 The lemmas the primitives rest on

> **Lemma 2 (a facet is its incidence set).** Let `y, y'` be facet normals with `I_y = I_{y'}`.
> Then `y' = λ y` with `λ > 0`. Consequently the map *facet ↦ incidence set* is injective, and
> orbit bookkeeping may run on incidence keys alone.

*Proof.* `F_y = C ∩ y^⊥ = cone{a_i : i ∈ I_y}` depends only on `I_y`, and `F_{y'} = F_y`; a
facet is the intersection of `C` with a supporting hyperplane, and the hyperplane supporting a
facet with nonempty relative interior is unique up to positive scaling. ∎

> **Lemma 3 (facet test, with a free upper bound).** Let `y ∈ Z^n`, `y ≠ 0`, primitive, with
> `y·a_i ≥ 0` for all `i` and `I = {i : y·a_i = 0}`. Then `F_y` is a facet of `C` **iff**
> `rank{a_i : i ∈ I} = n-1`. In that case `rank ≤ n-1` holds for free, since `y` annihilates
> `I`; only the lower bound needs a certificate.

> **Lemma 4 (adjacency is a joint rank, with a free upper bound).** Let `y ≠ y'` be facet
> normals. The smallest face of the dual cone `C*` containing both rays is `ann(I_y ∩ I_{y'}) ∩
> C*`, of dimension `n - rank{a_i : i ∈ I_y ∩ I_{y'}}`. Hence `y, y'` span a 2-face — they are
> *adjacent* — **iff** `rank{a_i : i ∈ I_y ∩ I_{y'}} = n-2`. Since both normals annihilate the
> joint set and are independent, `rank ≤ n-2` is free; only the lower bound needs a certificate.

*Proof.* Faces of `C*` are `ann(J) ∩ C*` for `J ⊆ {1..m}`, and the polarity
`dim(ann(J) ∩ C*) = n - rank{a_i : i ∈ J}` for a pointed full-dimensional `C`; the face spanned
by a set of rays is the intersection of the facets (of `C*`) containing them, i.e. `ann(J) ∩ C*`
with `J` the intersection of their incidence sets. ∎

> **Lemma 5 (a dead end — recorded so it is not repeated).** "Some ray is strictly negative on
> both normals" is a *necessary* condition for non-adjacency, but it is **vacuous**: facet
> normals are nonnegative on every ray, so the condition is never satisfiable and filters
> nothing. Measured: on the cone over the 9-cube the filter prunes 0 of 153 facet pairs while 9
> of them are genuinely non-adjacent. An earlier draft of this layer used that condition as the
> decision procedure; `code/bench.py` caught it (144 of 153 pairs decided correctly against the
> analytic truth 144). (P2) must be the exact joint-rank test.

> **Lemma 6 (transport).** Let `g` be a cell stabilizer generator with position permutation `π`,
> i.e. `g(a_i) = a_{π(i)}`. Then `g` maps the facet `(y, I)` to the facet `(y_g, π(I))` where
> `y_g` is the unique primitive normal of the hyperplane spanned by `{a_j : j ∈ π(I)}`. So a facet
> orbit is transported either by the integer action of `g` on `y` or by re-deriving the normal
> from the permuted incidence via (P3); the two must agree exactly, and that agreement *is* the
> transport certificate (cost `O(n^4)`), not an assumption.

## 2. The primitives (implemented, in `code/primitives.py`)

### (P3) `primitive_normal_from_incidence(incidence_rows, dimension) -> tuple[int, ...]`

1. Extract `n-1` rows of the incidence set that are independent over `Q`, with an exact
   certificate (`rank_lower_certificate`, below). If there are not `n-1`, the set is not a facet
   incidence and the primitive refuses (`ArithmeticError`) — a facet must own a codimension-one
   span.
2. Cofactor construction on those `n-1` rows: `y_i = (-1)^i det(M_i)`, `M_i` = the matrix with
   column `i` deleted. Cost `O(n^4)` integer operations, **independent of `|I|`**.
3. Divide by the entry gcd (primitivity) and fix the sign by the first nonzero entry.
4. Certify: `y·a_i = 0` for **every** row of the incidence set, not only the chosen ones.

This replaces a rational nullspace of an `|I| × n` matrix (the port's
`_pinned_nullspace_matrix(selected.transpose())`, `:1727`, cost `O(|I| n^2)` rational) by an
`O(n^4)` integer computation plus one certified rank check. Same answer, one canonical form:
because the result is primitive and sign-normalized, two routes that agree agree *exactly* — the
"up to positive scaling" comparisons the acceptance tests perform become equalities.

### (P1) `is_facet_certificate(normal, incidence, ray_rows) -> bool`

Exact conjunction of: `y ≠ 0` primitive; pairings vanish exactly on `incidence`; all other
pairings strictly positive; `rank` of the incidence rows `= n-1` via a certified independent
subset (Lemma 3: the upper bound is free). Every accepted facet is therefore a facet, verifiably,
with no appeal to the enumeration's correctness.

### (P2) `adjacent_from_zero_sets(first_zeros, second_zeros, ray_rows, dimension)`

`J = I ∩ I'` is a big-integer AND of two incidence bitsets (0.4 µs at 49 861 rays). Then rows of
`J` are reduced modulo a prime above the Hadamard bound until `n-2` pivots appear; the resulting
minor is re-certified exactly by a Bareiss determinant (Lemma 4: the upper bound is free). The
loop exits as soon as `n-2` independent rows are found, which for adjacent facets is within the
first `n-2` rows; non-adjacent pairs consume all of `J`, which is the honest worst case. Returns
the adjacency and the number of rows examined.

### Certified rank: `rank_lower_certificate(rows, wanted, *, full_rank=False)`

* Modular reduction with a prime `p > Hadamard bound` of the `n × n` minors only *finds candidate
  pivots* — it is never the certificate.
* The certificate is the exact Bareiss determinant of the square minor on the chosen rows' own
  pivot columns (upper triangular mod `p`, so it cannot vanish mod `p` unless it vanishes in `Z`;
  the Bareiss check settles it).
* `full_rank=True` additionally returns the rank of the whole set modulo the prime; used only by
  the test-suite cross-check against `matrix(QQ, …).rank()`.

### Bounded exact pairings: `pairings(normal, ray_rows)`

`int64` is used only after proving `|y|_1 · max|a_i(j)| < 2^63`; otherwise the exact
Python-integer path is taken. The check is a proof obligation, not a guard against silent
wrap-around. (For the specimen: entries ≤ 15 and cofactors below the Hadamard bound
`(√10·15)^{10} ≈ 6 · 10^{17}`; one covector's `|y|_1` is bounded by `n` times that, so the
fast path is chosen whenever `|y|_1 · 15 < 2^63`, and the fallback is exact when it is not.)

## 3. Completeness of the enumeration, with the hypothesis ledger

The following is the argument the layer certifies; each hypothesis is named and, where the port
supplies it, located.

**H1** `rank A = n` and `C` pointed — supplied by Lemma 1 from the cell invariants.
**H2** `G_cell` is finite and permutes the rows (so it acts on incidence keys) — supplied by
`cell_stabilizer` and `_facet_mask_data`'s packed permutations, which raise if an image is not a
row (`:1330`).
**H3** facets are exactly the certified records (P1) and facets ↔ incidence sets is injective
(Lemma 2).
**H4** adjacency is exactly the joint-rank condition (Lemma 4).
**H5** a codimension-2 face is contained in exactly two facets. *Proof:* let `R` be a `(n-2)`-face
of the pointed full-dimensional cone `C`. `R`'s span `V` has dimension `n-2`; `C` has a
supporting hyperplane through `V` on each side, and each supporting hyperplane through `V`
contains `R`; two independent supporting hyperplanes through `V` exist (otherwise `C ⊆ V`, false),
and a third would leave `R` in a smaller face. Hence exactly two. ∎
**H6** the facet graph of `C` (facets as vertices, ridges as edges) is connected. *Proof:* slice
`C` by an affine hyperplane `{ℓ = 1}` for a linear `ℓ > 0` on `C ∖ {0}`; the slice is a polytope
whose vertices are exactly the facets of `C` and whose edges are exactly the ridges (H5), so the
facet graph is a polytope's vertex-edge graph. That graph is connected: for a generic functional
`f`, every vertex that does not maximize `f` has an edge to a vertex with a strictly larger
`f`-value (if `f` were in the interior of the vertex's normal cone, the vertex would maximize
`f`), so following maximizing steps reaches the `f`-maximal face, which is a polytope of lower
dimension where the argument repeats; hence every vertex reaches the maximizing set, so the graph
is connected. ∎

> **Theorem (completeness of a closed orbit set).** Let `O` be a finite, nonempty set of facet
> orbits of `C` such that for every facet `F` in a member of `O` and every ridge `R ⊂ F`, the
> facet adjacent to `F` across `R` lies in a member of `O`. Then `O` is *all* facet orbits.

*Proof.* Let `S` be the union of `O` and suppose some facet `F* ∉ S`. Pick `F_0 ∈ S` and, by H6,
a facet path `F_0 = F_{0}, …, F_{k} = F*`. Let `j` be the first index with `F_{j} ∉ S`; then
`F_{j-1} ∈ S` and the ridge `R = F_{j-1} ∩ F_j` is contained in `F_{j-1}`. By closure the
facet adjacent to `F_{j-1}` across `R` is in `S`, and by H5 that facet is `F_j` — contradiction. ∎

This is precisely the reference's termination condition (`insert any orbit not already present;
mark done`) with the middle step made explicit, and it makes the algorithm's correct thing to
check a *local* one: for each enumerated orbit representative, all of its ridge neighbours are
covered. The queue never needs to know more.

**Termination.** The facet set is finite (H3), `O` only grows, and each iteration completes one
orbit; so the loop runs at most `|O|` times. The ridge subproblem for a representative reduces the
dimension by one (a facet cone inside an `(n-1)`-space), so the recursion is well founded with
depth at most `n-1` and bottoming out at `n = 2`, where a 2-dimensional pointed cone has exactly
two facets, both of which are its extreme rays.

## 4. Exactness: why no step can be a heuristic or a bound

* **Normals** are produced by cofactors from exact incidence rows (P3), verified against *all*
  rows of the incidence set, and normalized to a canonical primitive integer vector. There is no
  LP, so there is no optimality tolerance; there is no perturbation, so there is no degenerate
  tie to break.
* **Facets** are accepted only with (P1), whose rank lower bound is an exact nonzero minor and
  whose rank upper bound is a proof (the normal annihilates the set), not an estimate.
* **Adjacency** is decided by (P2) with the same two-sided structure; the only prime in the code
  is used to *propose* pivot rows and every proposal is re-verified exactly.
* **Transport** is the integer action of a stabilizer generator (`g`'s row matrix) composed along
  a Schreier word, and is verified by re-deriving the image facet's normal from the image
  incidence (Lemma 6). A wrong word or a wrong key raises `ArithmeticError` instead of silently
  producing a wrong orbit.
* **The two "bounds" in the code are proofs, not approximations**: the Hadamard bound sizes the
  certificate data, and the `int64` bound is *checked* before the fast path is taken, with an
  exact Python-integer fallback behind it.
* **Canonicalization removes scale comparisons.** Because normals are primitive with a fixed sign
  convention, matching an enumerated facet against a full-route covector is an equality test, not
  a search over positive rational multiples.

## 5. Termination and cost model (measured at the specimen's shape)

Notation: `o` = number of facet orbits of the cell, `r` = number of ridges of a facet (bounded by
the facet count of the `(n-1)`-dimensional facet cone), `m = 49 861`, `n = 10`.

* **Number of enumeration iterations** = `o` (Theorem, §3). Each iteration performs one ridge
  subproblem for its representative, so the number of dual-description calls is `o` and the total
  number of ridge reconstructions is at most `Σ_orbit r`.
* **Per-orbit cost of this layer** (measured, `results/bench.json`):
  one pairing pass over all rays 29 ms; one zero-support bitset 15 ms (a fraction of that if the
  incidence is sparse); one (P3) normal 0.05–0.3 ms (independent of `|I|`); one certified rank
  lower bound 0.35 ms at 256 rows, 2.9 ms at 4 096 rows, 11.5 ms at 16 384 rows.
* **Closure / orbit bookkeeping** (measured): exact (P2) decision 181–208 µs per facet-orbit pair.
  For `o = 200`: 3.0 s support builds + 3.6 s for all 19 900 pairs. For `o = 600`: 10.4 s + 37.5 s
  for all 179 700 pairs. Both are two to three orders of magnitude below the 1 800 s cap.
* **Why the full 49 861-facet materialization never happens.** The layer's state is one primitive
  normal (`n` integers), one incidence bitset (`m` bits = 6.1 KB) and one Schreier word per facet
  *orbit*; there is no per-facet ray-indexed object and no facet list. At `o = 600` that is ≤ 4 MB,
  against the 5 405 MiB peak of the blocked Normaliz call. The only per-ray work is the pairing
  pass of one covector per orbit, which produces a bitset and no facet.
* **Honest limits.** `o` and `r` on the A7 specimen are *unknown* — they are exactly what the
  blocked `SupportHyperplanes` call would have produced — so the enumeration's own cost is bounded
  by `o · r`, not measured. The measurements above are of the certificate layer at the specimen's
  `m` and `n`, with synthetic covectors; they establish the layer's cost, not the route's total.

## 6. Implementation plan: exact signatures and the places that change

New module `src/sage_indefinite_port/indefinite/facet_certificates.py`:

```python
def rank_lower_certificate(rows: tuple[tuple[int, ...], ...], wanted: int,
                           *, full_rank: bool = False) -> tuple[tuple[int, ...], int | None]: ...
def exact_minor_is_nonzero(rows: Sequence[Sequence[int]]) -> bool: ...
def independent_rows(rows, wanted: int) -> list[tuple[int, ...]] | None: ...
def primitive_normal_from_incidence(incidence_rows, dimension: int) -> tuple[int, ...]: ...
def orient_normal(normal: tuple[int, ...], outside_ray) -> tuple[int, ...]: ...
def zero_support(normal: tuple[int, ...], ray_rows) -> int: ...      # incident bitset
def adjacent_from_zero_sets(first_zeros: int, second_zeros: int, ray_rows,
                            dimension: int) -> tuple[bool, int]: ...
def is_facet_certificate(normal, incidence: Iterable[int], ray_rows) -> bool: ...
def pairings_bounded(normal, ray_rows) -> "np.ndarray": ...
```

Wiring, with the current call sites:

1. **`flip_across`, `lorentzian_cells.py:1708`** — replace the rational nullspace at `:1727`
   (`_pinned_nullspace_matrix(selected.transpose())`, cost `O(|I_F| n²)` with rationals) by
   `orient_normal(primitive_normal_from_incidence(selected_rows, n), vectors[nonincident])`. The
   `facet_nullspace.nrows() != 1` check becomes the `ArithmeticError` of `independent_rows`. The
   `base_normal` computation at `:1737` stays exactly as it is: its nullspace is over the *whole*
   row set with a unique solution and is not a facet computation.
2. **`_facet_mask_data`, `:1330`** — the inner `facet_keys()` generator currently consumes
   `_iter_cell_halfspace_pairing_rows(...)`, i.e. the engine's complete H-list, and converts each
   covector into a key with an `O(m · #facets)` blockwise pairing. On the orbit route the keys
   arrive from the enumeration instead: one key per orbit plus its Schreier orbit, produced by
   `zero_support` on the orbit representative and `_permuted_facet_key` for the rest. Everything
   downstream (`_iter_facet_mask_orbits` at `:1518`, `_facet_incidence`, `facet_orbit_permutations`
   at `:1279`) is unchanged — it already speaks keys and permutations, which is exactly the
   interface this layer produces.
3. **`cell_cone`, `:1011`** — keep the current branch unchanged for cells below the existing
   size threshold, so the current behaviour is preserved where it is already correct and both
   existing tests keep exercising it. For large cells return an orbit-indexed cone object exposing
   the same consumer surface (`halfspace_covectors()`, `primitive_rays()`, pairing queries) but
   backed by records:

   ```python
   @dataclass(frozen=True)
   class FacetOrbitRecord:
       normal: tuple[int, ...]        # primitive, oriented, = (P3) of `incidence`
       incidence: int                 # bitset over cell.rows, = `zero_support(normal, rows)`
       word: SchreierWord | None      # carries the representative to this facet
       stabilizer_generators: tuple[int, ...]   # indices of generators fixing the orbit

   def cell_facet_orbits(self, cell: LorentzianPerfectCell) -> tuple[FacetOrbitRecord, ...]: ...
   ```
4. **`perfect_domain_traversal`, `:2790`** — the per-cell block beginning `cell_cone = backend.
   cell_cone(cell)` consumes `cell_cone.halfspace_covectors()` and pairs each covector against the
   whole ray matrix. On the orbit route it iterates the records and evaluates each `normal` once
   against the ray matrix (`pairings_bounded`), which is the same incidence data at `o` passes
   instead of `#facets` passes. `_expand_traversal_once` (`:2564`) needs no change: it already
   takes `orbit[0][0]` as the representative key per orbit.
5. **Caps and guards**: nothing changes. No assertion is relaxed, the 1 800 s cap and the
   6 144 MiB guard stay, and the A7 node is still a single verification run — this layer is
   measured on the two bounded specimens and on synthetic cones of the same shape.

## 7. Validation

**Run here (this VM's Sage, `code/bench.py`, `results/bench.json`).** Five cones, 499 facets:

| cone | facets | normal recovered from incidence alone | certificate accepted | matches Sage's facet list |
| --- | --- | --- | --- | --- |
| cone over the 9-cube (dim 10, 512 rays, 256-ray incidences) | 18 | 18 | 18 | 18 (closed form) |
| cone over the 3-cube (dim 4, 27-ray scale) | 6 | 6 | 6 | 6 (closed form) |
| random dim 4, 20 rays | 18 | 18 | 18 | 18 |
| random dim 5, 30 rays | 92 | 92 | 92 | 92 |
| random dim 6, 40 rays | 365 | 365 | 365 | 365 |

`rank_lower_certificate` against `matrix(QQ, …).rank()`: 120/120 trials agree. Adjacency by
(P2) on the 9-cube cone reproduces the analytic answer (144 of 153 pairs adjacent, the 9
non-adjacent pairs being the opposite facets `x_0 - x_j`, `x_0 + x_j`), and (P2) had to read at
most 65 rows of any joint incidence set.

**Tests that must keep passing** (unchanged, and quoted from their own docstrings):
`tests/test_large_cell_cone_route.py` — the two cone routes agree and recover exactly the
configuration rays; `tests/test_facet_orbit_contracts.py` — a facet covector annihilates exactly a
rank-`(n-1)` incidence and is non-negative (strictly positive) elsewhere, stabilizer generators
permute rays and send facets to facets, and every facet lies in exactly one orbit. Both are
exactly what (P1) and Lemma 2–4 assert locally, so the new route can satisfy them per facet
instead of only globally.

**New differential test (bounded, no A7 run).** Cones over sub-configurations of the A7 cell's
rows: take `k ≤ 2 000` rows of the A7 cell, keep the rows' rank at `n = 10`, and compare the full
route (`RationalPolyhedralCones.from_rays` / `_normaliz_cell_h_presentation`) against the orbit
route on that cone. A sub-configuration is not itself a perfect cell, so the test constructs the
cone directly from rows; `LorentzianPerfectCell` is not involved. Assert (i) every enumerated
orbit passes (P1); (ii) every full-route facet matches exactly one orbit normal, as an equality
after the canonical form; (iii) every orbit member's normal is reproduced by the Schreier
transport. This exercises rank 10, entries ≤ 15 and heavy degeneracy without touching the 49 861
row call.

## 8. What this layer does not establish

* The facet-orbit count `o` of the A7 cell, and therefore the enumeration's total cost: `o` is
  what the blocked engine call would have produced. No A7 run is made here.
* The ridge subproblem's cost: that is a facet enumeration in dimension `n-1 = 9` over `|I_F|`
  generators per orbit, and its size depends on `o`/`r` as above. Directions d1/d4 own it.
* Port integration beyond the signatures above: direction d5 owns the seams.
* One caveat on `zero_support` as the sole adjacency state: for a *full* closure certificate one
  also needs, per orbit, that its ridge neighbours were computed, which requires the ridge
  subproblem; (P2) decides adjacency but does not produce the ridge.
