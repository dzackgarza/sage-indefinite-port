# Orbit-representative facet enumeration for perfect-cell cones: construction, proof, exactness, cost, implementation and validation plan

Scope. This document specifies how `LorentzianPerfectLocalBackend.cell_cone`'s expensive call
`NmzResult(normaliz_cone, "SupportHyperplanes")` (`src/sage_indefinite_port/indefinite/lorentzian_cells.py:1051`)
is replaced, on the traversal path, by an enumeration of the cell cone's facets **modulo the cell stabilizer**:
one representative per orbit, with exact transport words, no heuristic pruning, and no materialization of the
complete facet family. It also states the exact contracts, the completeness/termination proof, the exactness
argument for every quantity (normal, incidence, transport, adjacency witness), the cost and memory ledger, the
function-level implementation plan, and the validation plan on bounded specimens.

Everything below is stated for the port's data types. Nothing here proposes a provider (`sage-extensions`) patch;
every routine named is in the project's own tree.

---

## 1. Data, conventions and contracts

### 1.1 The cell and its cone

A cell is `LorentzianPerfectCell(lattice, rows, allowed_squares)` (`:587`): `rows` is a tuple of `m` integer
row tuples of length `n`, `n = lattice.module_rank()`, and the cell cone is

```
C = cone(rows) = { sum_i lambda_i rows[i] : lambda_i >= 0 }  subset  R^n .
```

Hypotheses used in the proof, all checked on the target specimen:

* **(H1) Full-dimensionality**: `rank(rows) = n`. Measured: rank 10 = ambient dimension on the A7 specimen, so
  there is no lineality and `C` is a pointed, full-dimensional cone.
* **(H2) Finite stabilizer**: `G = <cell_stabilizer(cell)> <= GL_n(Q) ∩ GL_n(Z)` is a finite group of integral
  isometries that permutes the rows. `cell_stabilizer` (`:1228`) returns generators only (not the whole group).
* **(H3) Rational data**: all rows are integral; all arithmetic below is in `ZZ`/`QQ`, never floating point.
* **(H4) Non-degeneracy of the affine slice**: the rows lie in one affine hyperplane not containing the origin
  (`__post_init__`, `:602`: the augmented family `[[1, *row]]` has rank `n`). This is what makes the port's
  affine-normal computations (`flip_across`, `:1708`) well defined and gives the convention dictionary of §1.4.

### 1.2 Row-matrix convention for `EXT`, and what the reference hands the dual description

Upstream `f_adj` (`references/vendor/polyhedral_common@1592b246/src_lorentzian/lorentzian_perfect.h:1433-1472`)
computes the stabilizer of the *row matrix* `x.EXT` and then calls

```cpp
vectface TheOutput = DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, data.rddo);
```

with `EXT_T = UniversalMatrixConversion<T,Tint>(x.EXT)` (`:1443`) and `x.EXT = MatrixFromVectorFamily(ListTotal)`
(`LORENTZ_GetEXT`, `:964`). So the contract is exactly:

* **input**: the `m × n` matrix whose *rows* are the configuration vectors — **not** augmented with a leading
  `1`, and **not** transposed — together with a permutation group on the `m` row positions;
* **output**: `vectface`, i.e. a family of `Face` values, each a 0/1 vector of length `m`; entry `i` is `1`
  exactly when row `i` lies on the facet.

That is the same object as the port's `FacetIncidence = frozenset[int]` (`:123`) and its packed key
`FacetMaskKey` (`:105`, `_facet_complement_key`), where a key encodes either the *incident* or the
*non-incident* position set, whichever is smaller (`:97-121`), and `_facet_incidence(key, ray_count)` (`:123`)
recovers the exact incident set. **The output is incidence-level; no covector is part of the reference's contract.**
The port's traversal already only consumes incidences: `_expand_traversal_once` (`:2564`) calls
`flip_across(source, facet)` with a `FacetIncidence`, and `flip_across` recomputes the facet's normal itself from
the incident rows. This document keeps that shape: normals are computed exactly but are *certificates* used
inside the enumeration, not part of the traversal-facing contract.

Because the rows span an affine hyperplane not through the origin (H4), the two possible readings of `EXT` — "these
are the vertices of a polytope" (so the output lists facets of `conv(rows)`) and "these are generators of a cone"
(so the output lists facets of `cone(rows)`) — give **the same incidence family**: every facet hyperplane of the
polytope, homogenised, is a facet hyperplane of the cone, and conversely every facet of the cone meets the affine
slice in a facet of the polytope. Lemma A (§1.4) maps the corresponding normals onto each other exactly, so no
information is lost or invented when the reference's affine data is read by the port's homogeneous machinery.

### 1.3 Stabilizer permutation convention

`cell_stabilizer(cell)` returns integral isometries (generators). For a generator `g` let
`A_g = matrix(ZZ, [coordinate_row(g(lattice.module_generator(label))) for label in labels])`, so that for a row
vector `v` the image is the row vector `v * A_g`. The port's permutation of positions is

```
pi_g[i] = j   <=>   rows[i] * A_g == rows[j]          (point_permutations_from_rows, :1246; _facet_mask_data, :1362-1390)
```

i.e. `pi_g` is the *image position array* (`pi_g[i]` is where position `i` goes), consistently with
`_permuted_bitset(mask, pi)` (`:181`), which maps source bit `p` to target bit `pi[p]`. `_facet_mask_data`
builds these arrays blockwise (`:1377-1389`) and stores them packed as `array` with typecode `H`/`I`/`Q`.

Schreier words are the port's nested-tuple `SchreierWord = (parent_word, generator_index) | None` (`:172`).
`_schreier_generator_indices` (`:172`) yields generator indices in *application order*: the composition loop

```python
permutation = array(typecode, range(ray_count))
for generator_index in _schreier_generator_indices(word):
    permutation = array(typecode, (generator[position] for position in permutation))
```

(`facet_orbit_permutations`, `:1305-1316`) applies `pi_{g_{i1}}, pi_{g_{i2}}, ...` from left to right, so the word
`(i1, i2, ..., ik)` denotes the composed image map `pi_{g_{ik}} o ... o pi_{g_{i1}}`. Any implementation and any
test must use exactly this order; §6 pins it with an explicit composition check.

### 1.4 The exact dictionary between the reference's affine normal and the port's homogeneous covector

`LORENTZ_DoFlipping` (`:1075-1100`) works with **affine** functionals on the augmented rows `(1, v)`: it takes the
unique `eNSPbas` vanishing on all of them (the "height" functional) and the facet's normal direction
`eNSPdir = (0, TheDir)`. The port's facets are **homogeneous** covectors `u` of the cone. The two conventions are
related exactly:

**Lemma A (affine/homogeneous dictionary).** Let `[alpha, h]` be the unique (up to scale) augmented functional with
`alpha + h·v = 0` for every row `v` (this is what `_pinned_nullspace_matrix(matrix([[1,*row] for row in rows]).transpose())`
returns in `flip_across`, `:1737-1743`), set `c0 := -alpha = h·v` (independent of `v`, and `c0 != 0` by (H4)).
Let `[beta, w]` be the augmented functional vanishing on a set `S` of rows (`beta + w·v = 0` for `v in S`).
Then

```
u = w - (beta/alpha) * h        (equivalently u = w + (beta/c0) * h)
```

satisfies `u·v = beta + w·v` for **every** row `v`, so `u` vanishes on exactly `S` and has exactly the signs of the
affine functional on all rows. Conversely, for a homogeneous covector `u` and a chosen row `v0`,
`w = u - (u·v0/c0) * h` and `beta = u·v0` recover an augmented functional agreeing with `u` on all rows.

*Proof.* `u·v = w·v - (beta/alpha)(h·v) = w·v - (beta/alpha)(-alpha) = w·v + beta`. Signs and zero sets are
therefore identical, and the inverse identity is `beta + w·v = u·v0 + u·v - (u·v0)(h·v)/c0 = u·v` because
`h·v = c0`. ∎

Consequence: any upstream incidence, produced with affine normals, is *the same set* the port needs, and the
primitive integral homogeneous normal is obtained from `(beta, w)` by clearing denominators and dividing by the
content of `u`. This is the only place where the two conventions can be confused; the port's `flip_across` already
performs the equivalent computation in one direction (`:1725-1745`), and §5 reuses it rather than duplicating it.

### 1.5 Output contract of the enumeration

`iter_cell_facet_orbits(cell)` yields, one orbit at a time, a list

```
orbit = [ (key_0, word_0), (key_1, word_1), ... ]        word_0 is the empty word
```

with:

1. **Completeness per orbit.** `key_0` is the orbit's canonical representative (the lexicographically least key in
   the orbit, see §2 step 5); every orbit member occurs exactly once; the orbit is closed under all generator
   permutations.
2. **Exact transport.** For each `j`, applying `_schreier_generator_indices(word_j)` to `key_0`'s incidence yields
   exactly `_facet_incidence(key_j, m)`.
3. **Orbit list.** The union of all yielded orbits is exactly the set of facets of `C`; distinct orbits are
   disjoint.
4. **Witness set** (optional, `retain_witnesses=True`): per representative the primitive integral normal `u_0`
   with `u_0 · rows >= 0`, `#{i : u_0·rows[i] = 0} = |incidence| >= 1`, `rank(rows[incidence]) = n-1`; per ridge the
   pair (ridge incidence, adjacent-facet key) and the ratio-test row `y*` of §2 step 4.

The compact-key type is the existing `FacetMaskKey`; no new incidence representation is introduced.

---

## 2. The construction

Notation: `Facets(C)` = the set of facets (each a subset of row positions), `Ridges(F)` = the facets of the cone
`F = cone(rows[i] : i in F)` inside `span(F)`, which are exactly the codimension-2 faces of `C` contained in `F`
(Lemma 2).

**Step 0 — stabilizer and ray permutations.** Compute `generators = cell_stabilizer(cell)` and their position
permutations exactly as `_facet_mask_data` does today (`:1330-1394`). This step is unchanged; it is factored out
into `_facet_ray_permutations(cell)` (§5) so the new route and the old route share one implementation. Complexity
`O(|generators| * m * n)` integer operations.

**Step 1 — one seed facet.** Produce a primitive integral `u` with `u·rows >= 0` and `rank(rows[{i : u·rows[i]=0}]) = n-1`.
Two admissible exact procedures; the second is the one specified for the port because it needs no lattice search:

* *(1a) Reference route.* `initial_cell`'s first flip: start from the augmented height functional
  `[alpha,h]` and any direction normal from the affine nullspace of the augmented rows, then run the port's pinned
  bisection `_upper_bound`/`_kernel_flipping` (`:5526`, `:5795`) once. This is the mechanism of
  `LORENTZ_SearchInitialVector` + `LORENTZ_Kernel_Flipping` and of `flip_across`.
* *(1b) Ridge-interval seed (specified).* Take any `n-2` rows spanning an `(n-2)`-dimensional space `R0`
  (grow it greedily from the row list by exact rank: `O(m n^2)`), compute a basis `n0, n1` of the 2-dimensional
  integer nullspace of the `(n-2) x n` matrix of those rows, and compute the cone
  `K = { (a,b) in Q^2 : a (n0·y) + b (n1·y) >= 0 for every row y }` by exact half-plane intersection in the plane
  (each row contributes one closed half-plane through the origin; `O(m)` rational operations, no search, no
  bisection). If `dim K = 2`, its two extreme rays give the candidate normals `a n0 + b n1`; normalise each
  primitively, orient by `u·rows >= 0`, and **certify it by Lemma 1** (nonnegativity plus tight rank `n-1`). If
  `dim K <= 1` or a candidate fails certification, retry with another `(n-2)`-subset. For a cone with at least one
  facet and `n >= 3` a certifying subset exists (take `n-2` independent rows inside any facet), so the retry loop
  terminates; the greedy rank growth finds it after `O(m)` tries in the worst case, and in practice on the first
  (see §6 for a bounded-specimen assertion that it does).

**Step 2 — worklist over orbits.** Maintain
`seen : list[FacetMaskKey]` (canonical key per *orbit only*, not per facet) and a FIFO `pending` of
`FacetMaskKey` values. Insert the seed's orbit (§2 step 5) and loop:

1. Pop `F` (a facet incidence key) and compute its primitive normal `u := primitive_facet_normal(cell, F)` as the
   exact integer nullspace vector of the incident rows, oriented by `u·rows >= 0`. Recomputing it from the incidence
   (one nullspace of an `(n-1)`-rank row set) keeps the loop free of transported covectors: the permutation action is
   applied to incidences only, where it is exact by §1.3.
2. **Ridges of `F`**: compute the facets of the cone `cone(rows[i] : i in F)` inside its own span, with an exact
   dual description on that sub-configuration (small input: §3 Theorem 2 and §6.5). Each returned covector `v`
   gives the ridge incidence `R_v = {i in F : v·rows[i] = 0}`; keep those with `rank(rows[R_v]) = n-2` and
   de-duplicate by `_facet_complement_key`. Only the zero set of `v` on `rows[F]` is used, so the choice of
   representative modulo the sub-cone's equation space (which the owned cone returns as
   `equation_covectors`) cannot affect the result.
3. **Adjacent facet per ridge (exact ratio test).** For each ridge incidence `R`:
   * compute the integer nullspace of the `|R| x n` matrix `rows[R]`; it is 2-dimensional and contains `u`;
     choose any second generator `w` (exact, `_pinned_nullspace_matrix` on the transposed int matrix);
   * compute `q_y = w·rows[y]` for all rows (one C-level matrix-vector product), and
     `s* = min_{y: q_y > 0} (u·rows[y]) / q_y` over `y not in R`, a rational, or `+inf` if the set is empty;
   * if `s* = +inf`, the ridge is a boundary ridge: record "no adjacent facet" (this is the port's existing
     boundary case, `flip_across` returning `None`, `:1721-1729`) and continue;
   * otherwise set `u' := primitive_integral(u + s* w)`; certify `u'·rows >= 0`,
     `rank(rows[{i : u'·rows[i] = 0}]) = n-1`, and `R` subset of the tight set. Then `F' = {i : u'·rows[i]=0}` is
     the *unique* facet of `C` other than `F` containing `R` (Lemma 3), and `(F', u')` is a new candidate.
4. **Novelty test (§2 step 5).** Compute the orbit of `F'` under the generator permutations by BFS with words.
   If any member's key equals a stored representative key, the orbit is already enumerated — discard.
   Otherwise append the canonical member (the minimum key) to `pending`, and store that key in `seen`.
5. **Orbit construction and canonical form.** The orbit of a key is the BFS closure under
   `_permuted_facet_key(key, pi_g)` for every generator `g`; words are extended as `(word, generator_index)`
   exactly as `_iter_facet_mask_orbits` does (`:1556-1585`), so the word convention of §1.3 is inherited
   unchanged. The canonical member is the minimum key under a total order on keys (compare `(0, int)` before
   `(1, tuple)` forms after converting to the incident-bitmask integers; the port's `_facet_complement_key` gives a
   canonical form for a *given* mask, and the minimum is taken over the orbit's masks in that canonical form).
6. **Emission.** Yield the completed orbit list. If the caller asked for witnesses, attach the normal `u_0` and the
   ridge/adjacency witnesses of the representative.

**Step 3 — termination.** The loop pops one orbit per iteration, marks it seen, and only pushes an orbit whose
orbit-closure is not represented in `seen`. Since `C` has finitely many facets and hence finitely many orbits, the
worklist empties; Theorem 1(iv) makes this exact.

**Step 4 — what is never computed.** No call ever passes `rows` (49861 vectors) to a V→H engine: the only dual
description call of the route receives the incidence rows of a single facet (`|F|` vectors, rank `n-1`), and its
output is the ridge family of that facet. The complete facet family of `C` is never materialized, not even
transiently, and no facet key outside `seen` (one per orbit) is retained.

---

## 3. Proofs

Throughout, `C = cone(rows)` is a pointed full-dimensional polyhedral cone of dimension `n` (H1), `G` a finite group
of integral isometries permuting the rows (H2).

**Lemma 1 (exact facet certificates).** For `u in Z^n \ {0}` primitive:
`u` is a facet normal (i.e. `C ∩ {u·x = 0}` is a facet) **iff** `u·rows[i] >= 0` for all `i` and
`rank(rows[{i : u·rows[i] = 0}]) = n-1`.

*Proof.* If the two conditions hold, `T = cone(rows[{i : u·rows[i]=0}])` is a face of `C` (it is `C ∩ u^⊥`
because any nonnegative combination with `u·x = 0` has all its `u`-positive terms zero), contained in the hyperplane
`u^⊥`, and `dim T = rank(rows[T-indices]) = n-1`; so `T` is a facet. Conversely a facet is `C ∩ u^⊥` for its
supporting `u`, the signs hold by support, and `dim = n-1` gives the rank condition. ∎

Because every facet has a supporting `u`, and scaling by a positive rational does not change the sign pattern or the
zero set, each facet has exactly one *primitive integral* such `u` (up to the sign fixed by the sign condition),
so "primitive integral normal" is a well-defined certificate. The rank and sign tests are exact integer
computations; there is no tolerance, no floating point, no rounding.

**Lemma 2 (ridges of `C` in `F` are the facets of the facet cone).** Let `F` be a facet of `C` with normal `u`, and
let `R subset F` be a face of `C` with `dim R = n-2`. Then `R` is a facet of the cone `F` (computed in
`span(F) = u^⊥`), and conversely every facet of the cone `F` is an `(n-2)`-face of `C`. Moreover each such `R` lies
in exactly two facets of `C`.

*Proof.* `F` is a pointed full-dimensional polyhedral cone in the `(n-1)`-dimensional space `u^⊥`, and its faces are
exactly the faces of `C` contained in `F` (face-lattice property of a face). Hence the facets of `F` (`dim n-2`) are
the `(n-2)`-faces of `C` inside `F`.

For the last claim let `R` be an `(n-2)`-face of `C` and put `S = { mu in span(R)^⊥ : mu·rows >= 0 }` — the normal
cone of `C` along `R`, restricted to `span(R)^⊥`, which is 2-dimensional. Then `S` is the set of supporting
functionals vanishing on `R`, and the facets of `C` containing `R` correspond bijectively to the extreme rays of
`S`. `S` is a cone (`mu·rows >= 0` is homogeneous), closed, and **pointed**: if `mu, -mu in S` then `mu` annihilates
every row, and `rows` spans, so `mu = 0`. It is nonzero because `R` is a face, and `dim S >= 1`. If `dim S = 1`
then `S = R_+ u` for the normal `u` of a single facet, so `R = C ∩ u^⊥` would be a facet, contradicting
`dim R = n-2 < n-1`; hence `dim S = 2`. A pointed 2-dimensional polyhedral cone has exactly two extreme rays, and
each of them is a facet normal: for an extreme ray `R_+ mu` the face `C ∩ mu^⊥` contains `R` and is strictly larger
than `R` (else `dim S = 1`), hence equals `n-1` and is a facet. Distinct facets have distinct normals, so `R` lies
in exactly two facets. ∎

**Lemma 3 (exact adjacent-facet oracle).** Let `F` be a facet with primitive normal `u`, `R` an `(n-2)`-face
contained in `F`, `w in Z^n` any integer vector in `span(R)^⊥` with `w notin span(u)` (such a `w` exists and
`span(R)^⊥ = span(u, w)` exactly), and

```
s* = min { (u·y)/(w·y) : y a row with y notin R and w·y > 0 }   (rational, possibly +inf)
u' = primitive_integral( u + s* w ).
```

Then: (i) `u'` is a facet normal; (ii) `R subset {i : u'·rows[i] = 0}` and `F' := {i : u'·rows[i]=0}` is a facet
`F' != F` containing `R`; (iii) if `s* = +inf` then no facet other than `F` contains `R`; (iv) consequently, by the
uniqueness in Lemma 2, `F'` is *the* facet adjacent to `F` across `R`, independently of the choice of `w`.

*Proof.* For `y in R`: `u·y = 0` (as `R subset F`) and `w·y = 0` (as `w in span(R)^⊥`), so `u'·y = 0`; hence
`R` is contained in `F'`. For `y notin R` compute `u'·y = u·y + s*(w·y)`:

* if `w·y <= 0`: `u·y >= 0`, `s* >= 0` and `w·y <= 0` give `u'·y >= 0`;
* if `w·y > 0`: `s* <= (u·y)/(w·y)` by minimality, so `u'·y >= u·y + (u·y/(w·y))(w·y) = 2 u·y >= 0`.

(Note `s* >= 0` because every candidate ratio `(u·y)/(w·y)` has a nonnegative numerator and a positive denominator.)
So `u'·rows >= 0`, and dividing a nonnegative vector by the positive content of its coordinates preserves that; the
sign certificate of Lemma 1 therefore holds after primitivisation. The tight set of `u'` contains `R` (rank `n-2`)
and the minimising row `y*`; and `y* notin span(R)`, because `w` annihilates `span(R)` while `w·y* > 0`. Hence the
tight set has rank `n-1`, so by Lemma 1 `F'` is a facet, and `R subset F'`. Finally `u'` is not proportional to `u`:
it equals `u + s* w` with `s* > 0` (the minimising set is nonempty) and `w notin span(u)`. So `F' != F`, and `F'` is
a facet other than `F` containing `R`. If `s* = +inf` there is no row outside `R` with `w·y > 0`, while `w` vanishes
on `R`; hence `w·rows <= 0`, so for any `mu = a u + b w` with `a, b in Q` and `b != 0` one of `mu`, `-mu` is negative
on some row, and the only supporting functionals in `span(R)^⊥` are the positive multiples of `u`: only `F` contains
`R`. Since by Lemma 2 there are at most two facets through a ridge, the construction in the finite case returns *the*
other facet through `R`, independent of which second generator `w` of the nullspace was chosen. ∎

*Remark.* Lemma 3 is a ratio test — the homogeneous simplex pivot — and is the exact replacement for the reference's
"find the adjacent facet `F0` with `R = F ∩ F0`" linear program: it is one matrix-vector product, one rational
minimum, one gcd, and one sign/rank certificate per ridge. It cannot be a heuristic: its output is certified by the
two conditions of Lemma 1, and its uniqueness claim makes it canonical.

**Lemma 4 (orbit novelty test is exact, and representative-only state suffices).** Let `G` act on keys by
`key |-> _permuted_facet_key(key, pi_g)`. Orbits are equal or disjoint. If `S` is a set of *canonical orbit
representatives* and `O(F)` the BFS closure of `F`'s orbit, then `O(F)` meets `S` iff `F`'s orbit is already
enumerated. Hence testing candidates against one key per orbit — instead of a `visited` set containing every facet —
is sound and complete.

*Proof.* Group orbits partition the set; closure under the generators is the orbit (standard). Therefore
`O(F) = G·F`. If `G·F = G·F0` for a stored representative `F0`, the BFS of `F` necessarily meets `F0`; conversely a
meeting point shows the orbits coincide. ∎

**Theorem 1 (completeness and uniqueness).** Assume (H1), (H2), Steps 1 and 2 are executed exactly, and the
ridge computation of Step 2.2 returns *all* facets of the facet cone `F` (induction hypothesis: the dual description
on the sub-configuration is complete). Then the enumeration terminates after exactly `N` iterations, where `N` is
the number of `G`-orbits of facets of `C`, and the emitted orbits are pairwise disjoint and cover `Facets(C)`.

*Proof.* Define the graph `Adj` on `Facets(C)` with an edge `{F, F'}` whenever `F ∩ F'` is an `(n-2)`-face
(Lemma 2 says this is the usual ridge-adjacency).

*(i) `Adj` is connected.* Since `C` is pointed and full-dimensional, `C°` has nonempty interior; pick
`lam in int(C°)`, which is strictly positive on `C \ {0}`. Then `P = {x in C : lam·x = 1}` is a polytope of dimension
`n-1` whose facets are exactly `F ∩ {lam·x = 1}` for `F in Facets(C)`, and two facets of `P` share a ridge iff the
corresponding facets of `C` share an `(n-2)`-face (Lemma 2). The facet graph (dual graph) of a full-dimensional
polytope is connected: for facets `F`, `F'` take `p in relint(F)`, `p' in relint(F')`; the segment `[p,p']` lies in
`P` by convexity, and after an arbitrarily small generic perturbation of `p'` inside `relint(F')` it meets no
`(d-2)`-face except at isolated crossing points, each of which lies in the relative interior of `F ∩ F'` for the
facet being left and the facet being entered; hence consecutive facets in the crossing order are adjacent. This is
the standard connectedness of the dual graph (Ziegler, *Lectures on Polytopes*, the dual graph of a polytope is
connected), so `Adj` is connected.

*(ii) Every facet is reached.* Let `F` be any facet, `F0` the seed. By (i) there is a chain
`F0 = F_0, F_1, ..., F_k = F` with `F_{j} ∩ F_{j+1}` a ridge. Inductively, when the orbit of `F_j` is processed,
Step 2.2 enumerates all its ridges (hypothesis), and Step 2.3-2.4 either inserts the orbit of `F_{j+1}` or finds it
already represented. Since orbits are closed under `G` and `G` maps facets to facets (H2), the whole orbit of
`F_{j+1}` is thereby either already enumerated or queued. So every facet orbit is enumerated.

*(iii) Exactly once.* Step 2.4 discards a candidate whose orbit meets `seen` (Lemma 4: that is exactly the case of an
already enumerated orbit), so distinct emitted orbits are disjoint and no orbit is emitted twice. Each insertion
marks the orbit seen, so no orbit is ever re-queued.

*(iv) Termination.* Each iteration removes one element from the finite worklist and adds only orbits not in `seen`;
`Facets(C)` is finite because `C` is a polyhedral cone with `m` generators (every facet is `C ∩ u^⊥` for one of
finitely many primitive normals — the finitely many that arise from the finitely many `(n-1)`-subsets). Hence at most
`N` insertions, and the worklist empties. ∎

**Corollary (the reference's route is recovered).** Steps 1 and 2 with (1a) as the seed and the pinned bisection
`_kernel_flipping` in place of Lemma 3 is precisely the documented upstream adjacency decomposition; the port-side
version specified here differs only in using the exact ratio test (Lemma 3) and in the representative-only state
(Lemma 4). Completeness and termination are Theorems 1 and 2 above either way, because both substitutions leave the
*set* of enumerated facets and their adjacency relations unchanged (Lemma 3 gives the same unique neighbour as the
reference's LP, and `_kernel_flipping` is exact by construction).

**Theorem 2 (cost and memory ledger).** Let `N` be the number of facet orbits, `r` the maximum number of ridges of a
facet, `m` the ray count, `n` the rank. Then the route performs

* exactly `N` dual-description calls, each on a facet incidence (`|F|` vectors) at rank `n-1` — never on all `m`
  vectors at rank `n`;
* exactly `N * r` ratio tests (Lemma 3), each `O(m n)` integer operations plus one `O(|R| n^2)` nullspace;
* exactly `N * r` novelty tests, each a BFS of one orbit: `O(|generators| * |orbit|)` key permutations, each
  permutation touching the set bits of one key;
* at most `N` retained facet keys at any time (Lemma 4) plus one orbit buffer at a time;

so the memory is `O(N * keysize + |orbit| * keysize)` instead of the `O(#Facets * keysize)` of the H-presentation
route, and the number of *expensive* (polyhedral) calls is `N`, independent of `#Facets` and of `m`. On the A7
specimen `m = 49861`, `n = 10`; `N` is the number of stabilizer orbits of facets (the A2 and A3 specimens have 6 and
24 facets and 4 and 6 *ray* orbits, which is the measured scale of the symmetry reduction; the facet-orbit count is
not measured for A7 and must not be guessed). The `Normaliz SupportHyperplanes` call on 49861 generators is not on
this path at all.

*Proof.* Count the operations of §2. The claims on memory follow from Step 2.5-2.6 (orbit buffers are transferred to
the consumer and released; see `_iter_facet_mask_orbits`'s ownership discipline at `:1560-1585`, which the plan
keeps) and Lemma 4 (one key per orbit in `seen`). The `N` calls follow from one ridge computation per popped orbit;
the `N r` ratio tests from one per ridge. ∎

---

## 4. Exactness: why nothing is a heuristic or a bound

* **Normal**: primitive integral, obtained as `u = w + (beta/c0) h` (Lemma A) or as `primitive_integral(u + s* w)`
  (Lemma 3); certified by Lemma 1's two conditions. There is no inequality-tolerance test anywhere.
* **Incidence**: `{i : u·rows[i] = 0}` computed in `ZZ`; packaged with `_facet_complement_key` and unpacked with
  `_facet_incidence`, exactly the port's existing exact representation.
* **Rank witnesses**: `rank` over `QQ` of the incident rows; the certificate is the rank itself, `n-1` for facets and
  `n-2` for ridges.
* **Adjacency witness**: the pair (ridge incidence `R`, minimising row `y*`) with the exact rational value
  `s* = (u·y*)/(w·y*)`, plus the two Lemma-1 certificates for `F` and `F'`. Verifying a witness is a finite exact
  computation; falsifying it is decidable. No LP, no floating point, no step size, no bound on a norm.
* **Transport witness**: the Schreier word is a word in the *given* stabilizer generators; its effect is the exact
  permutation obtained by composing the stored `array` permutations (§1.3). "The found facet maps onto the
  representative" is therefore an exact statement about a specific composition of specific integral matrices, and is
  checked by comparing keys, not by comparing covectors.
* **Isometry witness**: when a caller needs a categorical isometry rather than a permutation (the port's
  `facet_orbits`, `:1669`, materializes them), the composed product `A_{g_{ik}} ... A_{g_{i1}}` is integral and is
  certified by `A G A^T = G` plus row-permutation equality — the same certificate `cell_transporter` uses.
* **No hidden materialization**: the only polyhedral call receives one facet's incidence. `cell_cone` keeps its
  current behaviour for callers that genuinely need the complete H-presentation (`halfspace_covectors`,
  `cell_primitive_rays`, `cell_halfspace_covectors`), and those callers stay on the current route.

---

## 5. Implementation plan (exact signatures and call sites)

All new code goes in `src/sage_indefinite_port/indefinite/lorentzian_cells.py`, in
`LorentzianPerfectLocalBackend`. Reused owners (no parallel machinery): `cell_stabilizer` (`:1228`),
`_facet_complement_key` (`:105`), `_facet_incidence` (`:123`), `_permuted_facet_key` (`:140`),
`_permuted_bitset` (`:181`), `_schreier_generator_indices` (`:172`), `_pinned_nullspace_matrix` (`:5586`),
`_source_remove_fraction` (`:4196`), `_facet_mask_data`'s permutation block (`:1330-1394`), `flip_across` (`:1708`)
for the seed's cross-check, `_iter_facet_mask_orbits` (`:1518`) for the small-cell route, and
`_facet_orbit_transversal_masks`/`_facet_orbit_transversal_words` (`:1621`, `:1648`), `facet_orbits` (`:1669`),
`_expand_traversal_once` (`:2564`) as integration points.

```python
# 1. Factored-out permutation construction (behaviour identical to the block inside _facet_mask_data).
def _facet_ray_permutations(self, cell: LorentzianPerfectCell) -> tuple[array, ...]: ...

# 2. Exact certificates.
def facet_normal_certificate(self, cell, normal_row: tuple[int, ...]) -> tuple[FacetMaskKey, int]:
    """Return (packed incidence key, tight rank); raise unless the row is nonnegative on all rays."""

def primitive_facet_normal(self, cell, incidence: FacetIncidence) -> tuple[int, ...]:
    """Primitive integral homogeneous normal of the facet (Lemma A + normalisation); cached per call only."""

# 3. The two exact steps of the route.
def _cell_seed_facet(self, cell) -> tuple[FacetMaskKey, tuple[int, ...]]: ...
    # Step 1b (ridge-interval seed); falls back to Step 1a via _upper_bound/_kernel_flipping.

def _facet_ridges(self, cell, incidence: FacetIncidence) -> tuple[FacetIncidence, ...]: ...
    # One exact dual description on the incidence sub-configuration (rank n-1); de-duplicated ridge incidences.

def _adjacent_facet_across_ridge(self, cell, facet_normal_row, ridge: FacetIncidence
                                 ) -> tuple[FacetMaskKey, tuple[int, ...]] | None: ...
    # Lemma 3 ratio test; None exactly on boundary ridges (s* = +inf).

# 4. Orbit-level construction and the route selector.
def iter_cell_facet_orbits(self, cell, *, retain_words: bool = True, retain_witnesses: bool = False
                           ) -> Iterator[list[tuple[FacetMaskKey, SchreierWord]]]: ...

def _facet_orbit_source(self, cell, *, retain_words: bool
                        ) -> tuple[Iterator[list[tuple[FacetMaskKey, SchreierWord]]], tuple[array, ...]]:
    """Route selector: orbit enumeration for cells at/above the existing large-cell threshold
    (len(cell.rows) >= 10_000, the same predicate cell_cone uses at :1012), H-presentation otherwise."""
```

Call-site changes, minimal and additive:

1. **`_facet_mask_data` (`:1330`)**: extract its permutation block into `_facet_ray_permutations`; keep the
   H-presentation key generator unchanged for the small route and for callers that ask for `use_owned_cone=True`.
2. **`_facet_orbit_source` (new)**: returns `(orbit_iterator, permutations)` in the same *shape*
   `_facet_orbit_transversal_masks` (`:1621`) already returns, but with the orbits produced by
   `iter_cell_facet_orbits` on the large route and by `_iter_facet_mask_orbits(_facet_mask_data(...))` on the small
   route. This keeps one downstream code path.
3. **`_expand_traversal_once` (`:2564-2600`)**: replace the four lines that build `facet_masks`,
   `generator_permutations`, `canonical_facet_keys` and then iterate `_iter_facet_mask_orbits` with
   ```python
   orbits, generator_permutations = backend._facet_orbit_source(source, retain_words=self._retain_facet_orbits)
   for orbit in orbits:
       facet = _facet_incidence(orbit[0][0], len(source.rows))
       ...
   ```
   Everything after `facet = ...` (neighbour construction, `cell_transporter`, bucket lookup, adjacency record) is
   unchanged, including the `self._traversal_facet_orbits[source_id][facet] = orbit` bookkeeping at `:2630`.
4. **`_iter_facet_mask_orbits` (`:1518`)**: unchanged; it remains the exact orbit grouper for the small-cell
   H-presentation route. Its `canonical_keys` dictionary is only used there; the large route never builds it.
5. **`facet_orbits` (`:1669`), `_facet_orbit_transversal_words` (`:1648`)**: switch their source to
   `_facet_orbit_source`, keeping their public return types (`(FacetIncidence, LatticeIsometryMethods)`), and
   materialize the isometry from the word as today.
6. **`cell_cone` (`:1011`) and `_normaliz_cell_h_presentation` (`:1034`)**: unchanged. They are still correct for
   consumers that need the full H-family (`cell_halfspace_covectors` `:1144`, `cell_primitive_rays` `:1155`), and the
   acceptance assertions of the existing tests keep exercising them. The change is that the *traversal* no longer
   calls the H-presentation for large cells.
7. **Optional new public API for orbit-level consumers**:
   ```python
   def cell_facet_orbit_covectors(self, cell) -> tuple[tuple[tuple[int, ...], FacetIncidence, SchreierWord], ...]:
       """One primitive integral normal per facet orbit, with its incidence and the word carrying it."""
   ```
   This is what a caller that used `halfspace_covectors` for *orbits* (rather than all facets) should use; it does
   not materialize the full family.

Behaviour preservation: small cells (`len(rows) < 10_000`) take exactly today's path, so every existing test that
compares `cell_cone` with `RationalPolyhedralCones.from_rays` (the 27-ray and 89-ray specimens) is unaffected; the
large-cell behaviour changes only in *how the traversal obtains facet incidences*, not in what `cell_cone` returns.

---

## 6. Validation plan

**6.1 Both-route agreement on bounded specimens** (`tests/test_facet_orbit_contracts.py` fixtures:
`A2_plus_1_plus_minus1`, Gram `[[2,-1,0,0],[-1,2,0,0],[0,0,1,0],[0,0,0,-1]]`, 27 rays, 6 facets, 4 ray orbits;
`A3_plus_1_plus_minus1`, 89 rays, 24 facets, 6 ray orbits). For each specimen, with `H = cell_cone(cell).halfspace_covectors()`:

1. **Facet-orbit count** equals the number of orbits of `H` under the generator permutations (computed by the
   current `_iter_facet_mask_orbits`), as a set of orbit *sizes*, not merely as a count.
2. **Each enumerated representative is a facet of the full cone**: its incidence equals the zero set of exactly one
   `h in H` (exact cross-multiplication: `h` and the enumerated normal are positive multiples), and the match is a
   bijection between enumerated orbits and `H`-orbits.
3. **Completeness**: `set(union of all yielded orbit members) == set(key of every h in H)`; no extras.
4. **Witnesses**: for each member `(key_j, word_j)`, `_permuted_facet_key(key_0, composed(word_j)) == key_j`, where
   `composed` is the exact loop of §1.3 — this pins the Schreier-order convention.
5. **Adjacency witnesses**: for each ridge `R` of a representative, `F ∩ F' == R` (set equality) and `F'` occurs in
   the enumerated output; for boundary ridges `_adjacent_facet_across_ridge` returns `None` and the corresponding
   `flip_across(cell, R-containing facet)` also returns `None`.
6. **Exactness of the ratio test**: recompute `u'` from the witness `(R, y*, s*)` by exact arithmetic and check
   `u'` is the normal of the matched full-cone facet, up to positive scaling.

**6.2 Existing tests that must keep passing** (unchanged assertions):
`tests/test_large_cell_cone_route.py` (the two cone routes agree and recover exactly the configuration rays, with
exact positive-scaling matching) and `tests/test_facet_orbit_contracts.py` (covectors annihilate exactly a
rank-`(n-1)` incidence and are nonnegative elsewhere; generators permute rays and send facets to facets; every facet
lies in exactly one orbit). Both are small-specimen tests and stay on the existing route, so they are the regression
net for §5 items 1, 4 and 6.

**6.3 Bounded rank-5 sweep.** `A3 + <1> + <1> + <-1>`-style rank-5 lattices (89 rays, 24 facets, 6 ray orbits). Same
assertions as 6.1; this is where the ridge step first has non-simplicial facets, so it is the first real test of
Lemma 3's boundary case (`s* = +inf`).

**6.4 State and memory assertions.** On the A2 specimen, assert `len(seen) == number of facet orbits` (not
`number of facets`) by inspecting the internal counter, and assert with `tracemalloc` that the enumeration's peak
allocation is below a fixed multiple of `#orbits * key_size`. This is the test that pins the memory claim of
Theorem 2 and guards §1.5's contract against a regression to a full `visited` set.

**6.5 The A7 specimen.** *Not* run speculatively. It is a single verification run after 6.1-6.4 pass: the acceptance
assertions, the 1800 s cap and the 6144 MiB guard are unchanged; the run should report (a) the number of facet orbits,
(b) the number of ridge dual-description calls, (c) peak RSS, and it must not contain any
`NmzResult(..., "SupportHyperplanes")` on the traversal path. No claim about its runtime is made here: the cost
ledger of Theorem 2 bounds the operation count, and the only measured quantities available are the handoff's
(A7: 49861 rows, rank 10, 4.5 s to build the cell; A2: 6 facets/4 ray orbits; A3: 24 facets/6 ray orbits).

---

## 7. What is established here, what is inherited, and what remains unverified

* **Proved here**: Lemma A (affine/homogeneous dictionary, §1.4), Lemma 1 (exact facet certificate), Lemma 2
  (ridges of a facet are the facet cone's facets; each ridge in exactly two facets), Lemma 3 (exact ratio-test
  adjacency oracle with a canonicality/uniqueness statement), Lemma 4 (representative-only novelty test),
  Theorem 1 (completeness, uniqueness, termination, with connectivity of the facet graph proved for cones),
  Theorem 2 (call-count and memory ledger), and the exactness statements of §4.
* **Inherited from the vendored reference and the port** (not re-proved): that `cell_stabilizer` returns generators
  of the full stabilizer; that `flip_across`/`_kernel_flipping` implement the source's flip exactly; the existence
  and correctness of `DualDescriptionRecordFullDim`/`py_polyhedral.dual_description`; and the measured specimen
  sizes (27/6/4 and 89/24/6, A7 49861 rays rank 10).
* **Hypothesis that needs a small check before production**: that the ridge dual description on a facet incidence
  (rank `n-1`, `|F|` generators) is cheap on A7-size cells — it is a strictly smaller input than the blocked call in
  every dimension, but its cost is unmeasured. If it is not, the fallback is the engine-free recursion
  "facets of the facet cone by the same routine", which costs more calls but keeps every input below the full
  configuration; the completeness proof (Theorem 1) covers either implementation, because it only uses that the
  ridge set is complete.
* **Not verified in this VM**: no code was executed against the port. The installed Sage is 10.8 and
  `import sage_extensions` fails in `sage_extensions.rings.real` (`ValueError: base must be a subring of this ring`)
  while the acceptance workflow pins `ghcr.io/dzackgarza/sage:develop`, so the port and its tests could not be run
  here. The arithmetic evidence quoted in the swarm notes (facet covectors of shell configurations with Gram-inverse
  squares 2/3 and 5/3) comes from plain-Sage/PPL computations on shells I reconstructed, not from a port run.

---

## 8. Appendix: pseudocode of the route

```python
def iter_cell_facet_orbits(self, cell, *, retain_witnesses=False):
    m = len(cell.rows)
    perms = self._facet_ray_permutations(cell)        # §1.3 convention, exact
    seed_key = self._cell_seed_facet(cell)            # Step 1b (or 1a); incidence only
    seen_keys = []                                    # one canonical key per enumerated orbit
    pending = [seed_key]

    while pending:
        key = pending.pop(0)
        orbit = orbit_bfs(key, perms)                 # list[(FacetMaskKey, SchreierWord)], word_0 = ()
        rep_key = min((member_key for member_key, _ in orbit), key=canonical_order)
        seen_keys.append(rep_key)
        yield orbit                                   # ownership moves to the consumer; buffer released after use
        F = _facet_incidence(rep_key, m)
        u = self.primitive_facet_normal(cell, F)      # exact nullspace of rows[F], oriented by u·rows >= 0
        for R in self._facet_ridges(cell, F):         # one exact dual description, rank n-1, input = rows[F]
            found = self._adjacent_facet_across_ridge(cell, u, R)   # Lemma 3 ratio test + Lemma 1 certificates
            if found is None:                         # boundary ridge (s* = +inf); the port's flip_across also
                continue                              # returns None for exactly these faces
            next_key = found                          # incidence of the adjacent facet
            if orbit_meets(orbit_bfs(next_key, perms), seen_keys):   # Lemma 4 membership test
                continue
            pending.append(next_key)
```

`orbit_bfs(key, perms)` is the generator BFS with word extension of `_iter_facet_mask_orbits` (`:1556-1585`);
`canonical_order` is the total order of §2 step 5; `orbit_meets` tests the BFS closure against the stored
representatives only. `primitive_facet_normal` is `_pinned_nullspace_matrix(rows[F].transpose())` normalised to a
primitive integer row and oriented by `u·rows >= 0` — it is the same computation `flip_across` performs
(`:1725-1746`), reused rather than duplicated.

Note that the loop stores *incidences* and recomputes the representative's normal from its incidence, so no
covector ever has to be transported through a permutation: the permutation action is used only on keys, where it is
exact by construction (§1.3).

