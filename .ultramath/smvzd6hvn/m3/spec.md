# Orbit-representative facet enumeration for the perfect cell cones

Complete construction, correctness proof, exactness analysis, cost model and
implementation plan for replacing the per-cell complete H-description in
`src/sage_indefinite_port/indefinite/lorentzian_cells.py` by the reference's
orbit-wise dual description (`DualDescriptionRecordFullDim(EXT_T, GRPperm, rddo)`,
`lorentzian_perfect.h:1451-1452`).

Contents

1. Conventions and input/output contracts
2. Hypotheses, and the exact certificates that discharge them
3. The construction
4. Completeness proof
5. Exactness: normals, incidences, transport, Schreier words
6. Termination and cost model (with measurements)
7. Implementation plan: signatures, change sites, what stays as it is
8. Validation plan and the measurements taken here
9. What is not verified here

Companion files (all under `.ultramath/smvzd6hvn/m3/`):

* `prototype/orbit_addm.py` — the algorithm, exact, engine-free (the artifact
  the proof below is about);
* `prototype/validate_specimens.py` — agreement with exact complete dual
  descriptions on 17 specimens, dimension 3–10;
* `prototype/scaling.py`, `prototype/primitive_costs.py` — the cost
  measurements;
* `logs/validation.log`, `logs/scaling.log`, `logs/lp_seed.log`,
  `logs/pencil.log` — the raw measurement output.

---

## 1. Conventions and input/output contracts

**Ambient data.** `V = Q^n` with the lattice's coordinate-label basis fixed.
`cell.rows = (r_1, …, r_m)` are the cell's rays, `r_i ∈ Z^n`, pairwise distinct.
`C = cone{r_1, …, r_m} ⊆ V`.

**Dual pairing.** For `u ∈ V` (written as a coordinate row, i.e. a covector in
the *same* coordinates) and a ray `r`, `⟨u, r⟩ = Σ_k u_k r_k ∈ Z`. The lattice
Gram is never used by the enumeration: only integer pairings, kernels and
ranks. (The Gram enters only when the enumeration's output is *interpreted* —
see §5 — and in the port's stabilizer computation, which is upstream.)

**Facet representation.** A facet of `C` is stored as its **incidence**: the
`frozenset` of ray positions `F = {i : ⟨u, r_i⟩ = 0}` for a supporting covector
`u` (that is, `⟨u, r_i⟩ ≥ 0` for every `i`), together with the compact key the
port already uses (`FacetMaskKey`: sparse position tuple or packed bitset, plus
the complement flag — `_facet_complement_key`, `lorentzian_cells.py:105`). A
facet's **normal** is the unique primitive integral `u ≥ 0` on the rays with
`{i : ⟨u, r_i⟩ = 0} = F` (§5, Lemma 7).

**Group.** `G = ⟨σ_1, …, σ_k⟩` with `σ_j ∈ Sym(m)` the **position permutation**
of the rays: `g_j r_i = r_{σ_j(i)}` for the integral isometry `g_j` that the
port has already certified. Only the permutations are used. `G` acts on facets
by `σ(F) = {σ(i) : i ∈ F}` (`_permuted_facet_key`, `:140`).

**Reference contract, matched bit for bit.** In `f_adj`
(`lorentzian_perfect.h:1433-1472`) the reference calls

```
vectface TheOutput = DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, data.rddo);
```

with `EXT_T = UniversalMatrixConversion<T,Tint>(x.EXT)` and
`x.EXT = LORENTZ_GetEXT(eRec) = MatrixFromVectorFamily(eRec.ListTotal)`, i.e.
**the rows of `EXT_T` are the cell's rays**; `GRPperm` is the stabilizer's
**permutation action on those rows** (`LORENTZ_ComputeStabilizer`); the result
is a `vectface`, a list of **incidences** (index sets of the ray rows), one per
orbit. So the contract is exactly

```
(rows = ray matrix, generators = row permutations)
        ↦ facet incidences, one per stabilizer orbit.
```

`DualDescriptionStandard` at `edgewalk.h:1726` is the transpose convention
(it takes the covector matrix `ListIneq`), which is why the two routines share
one engine but are called with opposite data; the port needs the ray-matrix
version. The vendored GAP implementations
(`references/vendor/Indefinite.jl@374a5ebb/indef/lib/DualDescriptionAndADM.g`,
`__ListFacetByAdjacencyDecompositionMethod`, and
`Fundamental.g`, `__ProjectionLiftingFramework_Rational`, `FuncLift`) are the
readable source of the adjacency decomposition used below; `FuncLift` is
exactly the two-facet pencil of §3 (P2), including the corollary that the
adjacent facet is the unique pencil element that does not vanish on the current
facet.

The reference works with the *homogenized* configuration (`eNSPbas` is the
affine supporting functional; `LORENTZ_GetOnePerfect` returns
`{CritSet, eNSPbas, CentralVectT}`, `lorentzian_perfect.h:1026`), so its facet
cones carry a "facet at infinity" (`TestConicness` in `FuncLift`). The port's
`cell_cone` instead owns the *cone* `C = cone(rows)` in `V`; the translation is
exact and is needed only if the port ever wants the polytope picture: an
affine facet functional `(b, v)` (with `b + ⟨v, x⟩ ≥ 0` on the cell) becomes the
cone covector `v − b·h`, where `h` is the exact affine-hull functional
`⟨h, r_i⟩ = 1` that `flip_across` already computes. The construction below is
stated in the cone picture and needs no `h` except as the pointedness/duality
certificate and the LP normalisation of §3 (P3).

---

## 2. Hypotheses, and the exact certificates that discharge them

* **(H1) Full dimension, no lineality.** `rank(r_1, …, r_m) = n`.
  *Certificate:* one exact rank computation (`matrix(QQ, rows).rank()`), already
  performed by the port (`facet_orbits`, `:1673`, raises
  "does not span its ambient lattice"). Measured on the target specimen: rank
  10 = ambient dimension, so there is no lineality and no equation covector.

* **(H2) Pointed.** `C ∩ (−C) = {0}`.
  *Certificate:* the exact covector `h` with `⟨h, r_i⟩ = 1` for all `i`. It is
  computed by the port inside `flip_across` (`:1737-1743`) as the left kernel of
  the lifted rows `(1, r_i)`, and `flip_across` already raises unless that
  kernel is one-dimensional. Then every nonzero `x ∈ C` has
  `⟨h, x⟩ = Σ λ_i > 0`, so `C ∩ {⟨h, ·⟩ ≤ 0} = {0}`: `C` is pointed and `h` is
  strictly positive on `C \ {0}`. (Note `h` is a *covector*, so `⟨h, r_i⟩ = 1`
  is a statement about pairings; it is not the affine-hull functional of the
  polytope `conv(rows)`, and no affine geometry is used below.)

* **(H3) The group permutes the rays.** For each generator, `σ_j` is a
  permutation of `[m]` and `g_j r_i = r_{σ_j(i)}` with `g_j` an integral
  isometry.
  *Certificate:* `_facet_mask_data` (`:1330-1420`) builds each `σ_j` from the
  exact row action and raises `ArithmeticError("a cell-stabilizer generator
  does not preserve its rays")` otherwise; the isometry identity is checked in
  `perfect_domain_traversal` (`:3001-3007`).

* **(H4) Only orbits are wanted.** The traversal consumes, per cell, one facet
  per stabilizer orbit plus exact Schreier words
  (`_expand_traversal_once`, `:2564-2635`; the port's own comment at
  `:2834-2840`). Nothing downstream needs the complete facet list *as an input*
  — see §6 on what is still retained.

Two further hypotheses are **internal** to the enumeration and are proved, not
assumed: dual-description correctness on the facet sub-cones (it is the same
routine at one dimension lower, §4) and ridge adjacency correctness (§4, L1–L2).

---

## 3. The construction

Notation: for a face cone `Q = cone(S)` with `S ⊆ {r_1,…,r_m}` a ray subfamily,
`Facets(Q)` is its set of facets, each stored as the set of positions of `S`
whose rays it annihilates. `dim Q = rank(S)` because `Q` is pointed (§4, L0).

### P1. Facet normal from an incidence

```
facet_normal(S, A) -> u
```
for `A ⊆ S` with `rank(A) = dim Q − 1`: the kernel `{w ∈ V : ⟨w, r⟩ = 0 ∀ r ∈ A}`
is one dimensional; take its primitive integral generator and the sign with
`⟨u, r⟩ ≥ 0` for all `r ∈ S`; then check `{r : ⟨u, r⟩ = 0} = A` exactly.
Complexity `O(|A| n^3)` for the kernel and `O(|S| n)` for the check.
Implemented: `orbit_addm.py::facet_normal`.

### P2. The two-facet pencil of a ridge (the adjacent-facet operator)

```
ridge_pencil(S, span(S), R) -> (u, v)
```
for `R ⊆ S` with `rank(R) = dim Q − 2`. Let

```
W = { w ∈ span(S) : ⟨w, r⟩ = 0  for every r ∈ R }            (dim W = 2)
N_R = { w ∈ W : ⟨w, r⟩ ≥ 0 for every r ∈ S }                  (a pointed 2-cone)
```

**Lemma 1 (pencil lemma).** The extreme rays of `N_R` are exactly the normals
of the facets of `Q` that contain `R`; there are exactly two of them.
*Proof.* `N_R` is the normal cone of the face `R`; it is a face of the dual cone
`Q*`, of dimension `dim Q − dim R = 2`, and it is pointed because
`w, −w ∈ N_R` forces `⟨w, x⟩ = 0` on `C ⊇ Q`, i.e. `w ∈ Q^⊥ = {0}` by (H1).
A pointed 2-dimensional rational cone has exactly two extreme rays. If `w` is
extreme in `N_R` and `w = w_1 + w_2` with `w_i ∈ Q*`, then `⟨w_i, r⟩ ≥ 0` for
`r ∈ R` and `⟨w, r⟩ = 0` force `⟨w_i, r⟩ = 0`, so `w_i ∈ N_R` and `w_i ∥ w`:
`w` is extreme in `Q*`, i.e. a facet normal of `Q`. Conversely a facet normal
`u` with `R ⊆ Z(u)` lies in `N_R` and is extreme there (a 2-cone's extreme rays
are its boundary rays, and the boundary of `N_R` is where one more ray of `Q`
becomes tight; `u` is tight on all of `R` with `rank R = dim Q − 2`, so it is
a boundary point). ∎

*Computation, exact and `O(|S| n)`.* Take an integral basis `p, q` of `W` (the
left kernel of the `dim Q × |R|` matrix `span(S)·R^T`; `dim = 2` by
`rank R = dim Q − 2`). Put `a_i = (⟨p, r_i⟩, ⟨q, r_i⟩) ∈ Z²` for every ray of
`S`, and drop zeros. Then `y ∈ Q² ↦ y_1 p + y_2 q` identifies `N_R` with
`K = {y : ⟨y, a_i⟩ ≥ 0 ∀ i}`. Since `K` is a pointed 2-cone, the `a_i` span an
angular sector of width below 180°, on which `cross(b, a) < 0` is a total order;
the angular minimum and maximum are one pass each, and the two extreme rays of
`K` are the perpendiculars of those two extreme constraints:
`y₁ = (a_max₂, −a_max₁)`, `y₂ = (−a_min₂, a_min₁)`. Both are then **verified**
(`⟨y, a_i⟩ ≥ 0` for every `i`, exactly), mapped back to `V`, normalised to
primitive integral and sign-fixed by P1's test. No LP, no tolerance, no
iteration bound. Implemented: `orbit_addm.py::facet_normals_through_ridge`.

### P3. One facet of a cone (the seed)

Two exact implementations; both are *certified* before use, so any of them can
be replaced without affecting correctness (§4, L6).

* **P3a (inherited).** If `R` is known to be a ridge of `Q` and `F` a facet of
  `Q` containing `R`, then `R` is a facet of the facet cone `C_F`, i.e. a seed
  for the inner run. Free: `O(|R|)` set operations.
* **P3b (dual-slice LP).** For the cone with no known ridge:
  `S = {u ∈ V : ⟨u, r_i⟩ ≥ 0 ∀ i, ⟨h, u⟩ = 1, u ∈ span(S)}` is a *bounded*
  polytope (with the annihilator equations present only when `dim Q < n`;
  `h` is (H2)'s covector, strictly positive on `Q`). Its **lexicographic
  minimum** (minimise `u_1`, then `u_2`, …) is a unique point, hence a vertex of
  `S`; and the vertices of `S` are exactly the extreme rays of `Q*` (each
  extreme ray of `Q*` meets the slice `⟨h, ·⟩ = 1` once, and a 0-face of the
  slice is contained in a minimal face of `Q*`, which is 1-dimensional because
  `h` is interior). Hence the optimum *is* a facet normal of `Q`, and its
  incidence is returned after the `rank = dim Q − 1` **and** nonnegativity
  checks. No genericity assumption is needed (a total order is used, not a
  random objective); it is `dim Q` exact LPs over one feasible set, each
  pinning the previous coordinate with an equality.
  Implemented: `orbit_addm.py::seed_facet` (Sage `MixedIntegerLinearProgram`
  with the exact `PPL` backend, `QQ` values).

### P4. Setwise stabilizer, from the orbit transversal

```
setwise_stabilizer(generators, F) -> generators of Stab_G(F) on F's positions
```
Schreier's lemma applied to the action on the *set* `F`: with transversal
elements `t_h` (words) carrying `F` to the members `hF` of its orbit, the
elements `t_{s h}^{−1} · s · t_h` fix `F` and generate its setwise stabilizer.
Restricted to `F` these are permutations of `F`'s own positions.
Implemented: `orbit_addm.py::setwise_stabilizer_generators`. In the port this is
*the same machinery it already runs*: `_iter_facet_mask_orbits` (`:1518`)
already computes every transversal word of every facet orbit; the Schreier
generators are `|orbit|·|generators|` word compositions, exactly the bound the
port already pays when it materializes an orbit (`facet_orbits`, `:1669-1705`).

### The main routine

```
facet_orbits(S, H, seed) -> [ (representative incidence F, exact normal u, {member ↦ word}) ]
```
by induction on `dim Q`, where `H ≤ Sym(S)` is a group of permutations of `S`.

* `dim Q = 1`: return the single orbit `{∅}` — the origin is the only facet of
  a ray. (No normal: the covector would not be intrinsic.)
* `dim Q ≥ 2`:
  1. Seed: `F₀` from P3a (`seed`/`seed_ridge` argument) or P3b. Verify it.
  2. Bookkeeping: `where : member incidence ↦ orbit index`, a list of orbit
     records, and a queue of undone orbits (`insert`, guarded by `where`; the
     insertion also computes the orbit closure with words by BFS under `H`, so
     every inserted orbit is saturated).
  3. `insert(F₀)`; then while some orbit is undone:
     * take the orbit's representative `F`, and its discovery ridge `R_F`
       (the ridge through which it was inserted, or `None`);
     * `H_F = Stab_H(F)` restricted to `F` (P4), and the local index map;
     * **inner run**: `facet_orbits(F's rows, H_F, seed = R_F localized)` — that
       is, the facets of the facet cone `C_F`, i.e. the ridge orbits of `C`
       inside `F`, one per `H_F`-orbit;
     * for each inner representative `R`: `(u, v) = P2(R)`; the two facets of
       `Q` containing `R` are `Z(u)`, `Z(v)`; keep the one that is not `F` and
       `insert` it with discovery ridge `R`;
     * mark the orbit done.

`Q`'s facets are `Facets(Q)`; the top call has `Q = C`, `S = rows`,
`H = G =` the cell stabilizer.

Complexity-controlling devices, both exact: (i) the *bank* (`memo` in the
prototype) keyed by the sub-configuration `(rows, generators)` — two nodes with
the same ray set and the same restricted group are the same subproblem, which
is exact equality of the inputs; (ii) the inherited seed (P3a), which is what
makes the number of LPs `O(n)` rather than `O(#nodes)` (§6, Lemma 8).

---

## 4. Completeness proof

Throughout: `C` full dimensional and pointed (H1–H2), `G ≤ Sym(m)` permuting the
rays (H3), `Q` a face cone of `C` with ray list `S`, all faces pointed
(they are faces of the pointed cone `C`), so `dim F = rank F` for every face.

**L0.** `dim cone(S) = rank(S)`, and for a facet `F` of `Q`,
`rank(F) = dim Q − 1`. *(Faces of a pointed cone contain no line, so their
dimension equals the dimension of their linear span.)*

**L1 (ridge pencil).** Every ridge `R` of `Q` (a face with
`rank R = dim Q − 2`) is contained in exactly two facets of `Q`, and they are
the two outputs of P2. *Proof:* Lemma 1.

**L2 (ridges of a facet are its facets; bijection with neighbours).** Let `F`
be a facet of `Q`. Then

```
Φ: Facets(Q) \ {F} → Facets(C_F),   F' ↦ F ∩ F'
```
is a bijection, and `F ∩ F'` is a ridge of `Q`.
*Proof.* A facet of `C_F` is a codimension-one face of `C_F`, i.e. a face `R` of
`Q` with `R ⊆ F` and `rank R = dim F − 1 = dim Q − 2`; so it is a ridge of `Q`
contained in `F`. By L1 it lies in exactly one facet of `Q` other than `F`, call
it `F'`, and then `F ∩ F' ⊇ R`. Since `dim(F ∩ F') ≤ dim Q − 2 = rank R`, we
get `F ∩ F' = R`, so `Φ(F') = R`, i.e. `Φ` is surjective. Injective: if
`F ∩ F' = F ∩ F'' = R` then both `F'` and `F''` contain the ridge `R`; by L1 at
most two facets contain `R`, and `F` is one of them, so `F' = F''`. ∎
(L2 is the statement that the ridge enumeration of a facet cone *is* the
adjacency structure of the cone; it is why the recursion is not a detour.)

**L3 (`G` acts on facets, and inserted orbits are saturated.)** For
`σ ∈ G` induced by the isometry `g`, `F ↦ σ(F)` maps facets of `Q` to facets of
`Q`: if `u` is `F`'s normal then `u ∘ g^{−1}` is supporting and has zero set
`σ(F)` — an exact integer identity, not an approximation. The BFS closure
(`orbit_closure`) therefore returns the whole orbit, with a word per member.

**L4 (neighbour orbits from ridge orbits).** Let `O` be an orbit of `H` on
`Facets(Q)` with representative `F`. Then

```
{ H-orbits of facets other than F's orbit adjacent to some member of O }
   = { orbit of the second facet at R : R a representative of the H_F-orbits of Facets(C_F) }.
```
*Proof.* If `F'` is adjacent to `hF` along the ridge `hR`, then applying `h^{−1}`
gives a facet adjacent to `F` along `R`, whose `H`-orbit is the same as `F'`'s.
Conversely a ridge `R` of `F` gives the facet `F'` from L1. So it suffices to
run over `Facets(C_F)` modulo `Stab_H(F)`, which is exactly the inner run. ∎

**L5 (the facet graph is connected).** The graph with vertex set `Facets(C)` and
one edge per ridge is connected, `n ≥ 2`.
*Proof.* `P := C ∩ {⟨h, ·⟩ = 1}` is a polytope of dimension `n − 1` (H1, H2:
bounded because `h` is interior to the dual cone, and full dimensional because
`C` is). Facets of `C` correspond bijectively to facets of `P`, and ridges of
`C` to ridges of `P` (i.e. edges of the dual polytope `P*`). The 1-skeleton of a
polytope of positive dimension is connected; `P*` has dimension `n − 1 ≥ 1`.
For `n = 2`: `C` is a 2-cone, its two facets are its two extreme rays and the
only ridge is `{0}`, so the graph is a single edge. ∎
(The hypothesis used is exactly (H2): pointedness of `C` is what makes the
slice a polytope. Without it the facet graph of a cone with lineality can be
disconnected — e.g. a line has no facets at all.)

**L6 (the BFS enumerates every orbit).** Consider the quotient graph on the
`G`-orbits of `Facets(C)`, with an edge from `O` to the orbits of L4. It is the
image of a connected graph (L5), hence connected. The main routine performs
exactly a breadth-first search on this quotient graph: every dequeued orbit is
expanded once (it is marked done and cannot be re-inserted because `insert` is
guarded by `where`), and by L4 its expansion finds *all* orbits adjacent to it.
Hence from the seed's orbit every orbit is reached.
*Independence of the seed:* the routine computes the smallest set of orbits that
contains the seed's orbit and is closed under the adjacency relation of L4;
by L5 that set is everything. Therefore the choice of seed (P3a or P3b, or any
other exact facet that passes the same codimension and nonnegativity tests) can
affect the *route* but never the *result*.

**L7 (termination).** `Facets(C)` is finite, so the quotient graph has finitely
many nodes and the while loop halts after at most `#orbits` expansions (each
expansion marks exactly one unfinished orbit done and inserts finitely many
new ones). The recursion halts because `dim C_F = dim Q − 1` (L0) and the base
case is `dim = 1`, so the depth is at most `n`. The bank does not change
termination (it only replaces an identical subproblem by its stored answer).

**Theorem (completeness and partition).** Let `C` be a full-dimensional pointed
polyhedral cone with ray list `rows`, `G ≤ Sym(m)` permuting the rays, and let
`seed` be any incidence that passes the verification of P1/P3. Then
`facet_orbits(rows, G, seed)` terminates and returns a finite list of orbits
such that

1. every returned representative is a facet of `C` (P1/P2 verification plus
   L1); 
2. the returned orbits are pairwise disjoint (guarded by `where`);
3. the union of the returned orbits is `Facets(C)` (L5, L6);
4. every returned member carries an exact word in the generators taking the
   representative onto it (BFS construction), and every member's normal is the
   exact transport of the representative's normal (§5).

*No heuristic enters: the only choices are the seed (verified, and L6 shows the
result does not depend on it), the order of the queue, and the pivot order of
the linear algebra — none of which can change the output.*

---

## 5. Exactness

**(E1) Normals are canonical and exact.** For a facet incidence `F`, the
supporting covectors form a 1-dimensional cone: if `u, u'` are both `≥ 0` on
`C` and vanish exactly on `F`, then `u' = λu` with `λ > 0` (both lie in the
relative interior of the normal cone of `F`). The primitive integral generator
is therefore unique, and P1's sign rule (nonnegative on all rays) selects it.
It is obtained by an exact rational kernel and integer content/primitive
division; every step is integer/rational arithmetic with no bound.

**(E2) Incidences are exact.** `⟨u, r_i⟩` is a `Σ u_k r_{ik}` sum of Python/Sage
integers; `F = {i : ⟨u, r_i⟩ = 0}` is a zero test, not a threshold. The port's
key encoding (`_facet_complement_key`) is a lossless re-encoding.

**(E3) Transport of a normal along a word.** Let `w = (j_1, …, j_l)` be a word,
applied in that order, and `σ_w = σ_{j_l}∘…∘σ_{j_1}`. Let `M_w` be the row
matrix of the isometry, recovered *exactly* from the images of `n` independent
ray positions `i_1, …, i_n` (the pivots of `rows^T`):

```
M_w = B^{−1} · [ r_{σ_w(i_1)} ; … ; r_{σ_w(i_n)} ],      B = [ r_{i_1} ; … ; r_{i_n} ],
```

i.e. `x·M_w = g_w(x)` for row vectors (this is the port's existing pattern at
`:1687-1697` and `:2950-2960`). Because the rays span `V` (H1), the action on
them determines `M_w`; because the images of a lattice basis are lattice
vectors, `M_w` is integral — the port *checks* this
(`ArithmeticError("a Schreier facet permutation does not induce an integral
lattice action")`) together with `|det M_w| = 1` and
`M_w^T·Gram·M_w = Gram`. Then the member's normal is

```
u_w = u_0 · (M_w^{−1})^T          (equivalently the column M_w^{−1}·u_0),
```

which satisfies `⟨u_w, r_j⟩ = ⟨u_0, r_{σ_w^{−1}(j)}⟩`, so
`Z(u_w) = σ_w(Z(u_0)) = F` — the identity the validation checks for every
member. Nothing is scaled or rounded; `u_w` is integral because `M_w` is
unimodular.

**(E4) Words are exact.** A word is a nested tuple of generator indices (the
port's `SchreierWord`); composition is integer permutation evaluation. The
member-to-representative direction is the inverse permutation (an exact array
inverse, or the reversed word with inverse generators). The transport witness
for an adjacency is therefore `(orbit index, word, incidence)` with no search:
the isometry is `_owned_isometry_from_row_matrix(lattice, M_w) · transporter`.

**(E5) Why no step can be a heuristic.** The algorithm never discards a facet
(no pruning), never merges two facets without an explicit permutation and word,
and never evaluates a facet criterion in floating point. Its only "search-like"
element is the seed, and L6 shows the result is invariant under it; a wrong seed
can only fail the explicit `rank = dim − 1` / `⟨u, r⟩ ≥ 0` test, never produce a
silently different answer.

---

## 6. Termination and cost model

**Lemma 8 (how many seeds).** In the run of §3, the number of P3b (LP) seeds is
at most the number of codimensions, `≤ n`; every other node inherits its seed by
P3a.
*Proof.* Every node `Y` of the recursion (a face of `C` at some codimension) is
either (i) the *seed* of the inner run it belongs to, or (ii) discovered by the
pencil at a ridge `Z` of that inner run, in which case `Z ⊂ Y` and `Z` is a
facet of `C_Y` (L2 applied inside the parent cone, whose dimension is one
higher), so `Y`'s own inner run inherits its seed by P3a; or (iii) the root,
whose seed comes from outside. So an LP seed is needed only by the seeds of the
seed chain: the root, then the seed of its inner run, then the seed of *that*
node's inner run, and so on — one per codimension, at most `n`. In the
measurements, `dim = 10` runs used exactly 9 LP seeds and 31–33 inherited ones
(`logs/validation.log`). ∎
(The chain is unavoidable in general: finding one facet of a cone given by rays
is dual to finding one extreme ray of a cone given by inequalities, i.e. an LP.
This is precisely the reference's "some initial facet F (by linear
programming)".)

**Cost per unit, measured in this VM** (Sage 10.8, exact PPL backend, Python
prototype; `logs/lp_seed.log`, `logs/pencil.log`):

| operation | size | measured |
| --- | --- | --- |
| one ridge pencil: `2×n · n×m` exact integer product | `m = 49 861`, `n = 10` | **2.0 ms** |
| one ridge pencil: `O(m)` angular-extreme scan (Python) | `m = 49 861` | **11.2 ms** |
| one `O(m)` exact feasibility pass per candidate | `m = 49 861` | `~10 ms` |
| P3b seed, `dim` sequential exact LPs | `m = 144` | 168 ms |
| | `m = 2 034` | 2.4 s |
| | `m = 7 392` | 10.3 s |
| | `m = 22 608` | 64.8 s |
| orbit closure with words, 17 generators, frozenset form | 49 861 points, 10 000 members | 2.5 s |
| same, packed-bitset form (the port's `_permuted_facet_key`) | — | strictly less (bit shifts instead of tuple rebuilds) |

**Measured end-to-end scaling** (`logs/scaling.log`; the specimens are
perfect-cell-shaped: the vectors of a fixed norm in the affine hyperplane
`x_d = 1` of the Lorentzian lattice of Gram `I_{d-1} ⊕ ⟨−1⟩`, symmetry = signed
coordinate permutations; `full` is PPL's exact complete H-enumeration of the
same ray matrix):

| rays | dim | facets | facet orbits | nodes | LP seeds | pencils | this route | complete H-desc |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 3 | 4 | 1 | 3 | 2 | 2 | 0.01 s | 0.01 s |
| 12 | 4 | 14 | 2 | 5 | 3 | 5 | 0.02 s | 0.00 s |
| 24 | 5 | 24 | 1 | 7 | 4 | 9 | 0.03 s | 0.00 s |
| 90 | 6 | 232 | 3 | 14 | 5 | 28 | 0.10 s | 0.02 s |
| 252 | 7 | 508 | 3 | 21 | 6 | 56 | 0.27 s | 0.13 s |
| 574 | 8 | 1 108 | 3 | 28 | 7 | 86 | 0.76 s | 0.80 s |
| 1 136 | 9 | 2 416 | 3 | 35 | 8 | 116 | 2.05 s | 5.39 s |
| 2 034 | 10 | 5 264 | 3 | 42 | 9 | 146 | **5.40 s** | **42.31 s** |

The two curves separate exactly as the algorithm predicts: the complete
description costs one pivot per (facet × incident ray) *globally*, the orbit
route costs one pencil per *ridge orbit* plus one closure pass per facet.

**Projection to the target specimen (49 861 rays, `n = 10`).** With `ρ` the
number of pencils the run performs and `λ_k` the ray count of the `k`-th seed
chain face,

```
time  ≈  ρ · 15 ms                     (pencils, measured above at m = 49 861)
       + Σ_{k=0}^{9} LP(λ_k)            (≤ 10 exact LPs, Lemma 8)
       + #facets · #generators bit ops  (closure: 2.5 s per 10 000 members measured)
       + k₀ · O(m n)                    (one exact normal per orbit)
```

All three constants are measured at the target's own scale; only `ρ`, `k₀` and
`λ_k` depend on the specimen, and they are exactly the quantities the reference
also pays (`GetInitialRays` per recursion node, one pencil per ridge orbit).
For orientation: the dim-10 measurement needed `ρ = 146` pencils, and even
`ρ = 20 000` would be `300 s`. The blocked alternative needs one Normaliz
call that does not finish in 1 800 s on the same input.

**Why the complete H-description is never triggered.** The route asks an engine
nothing at all. Its inputs are the ray rows and the ray permutations. The only
global object it builds is the set of facet *incidences* (compact keys) with
their words — which `_expand_traversal_once`/`_traversal_facet_orbits` already
retain today, and which the port needs anyway to attach side pairings. It never
builds: a facet normal for a non-representative (transported on demand, §5 E3),
a facet × ray pairing matrix (the pencil's `2 × m` product is the only bulk
pairing), a polyhedron, or an `NmzCone`. If even the key set must be avoided,
the only change is the membership test "is this facet in a known orbit?", which
can be replaced by an exact canonical-image test under `G` (GAP is already
imported by the module: `libgap`); the rest of the algorithm is unchanged. This
is an optimisation, not a correctness requirement.

---

## 7. Implementation plan

### 7.1 New module-level engine (pure, in `lorentzian_cells.py`)

```python
def _facet_primitive_covector(covector: tuple[int, ...], rows, rank: int) -> tuple[int, ...]
    # content division + the P1 sign rule; raises unless ⟨u, r_i⟩ ≥ 0 and the
    # zero set has rank rank-1

def _facet_normal_of_incidence(rows, rank: int, incidence: FacetIncidence) -> tuple[int, ...]
    # P1; raises unless the incidence really is a facet of cone(rows)

def _ridge_pencil_normals(
    rows, span_basis, ridge: FacetIncidence, rank: int, *, ray_transpose=None
) -> tuple[tuple[int, ...], tuple[int, ...]]
    # P2 (Lemma 1); raises unless exactly two verified facet normals are found

def _lexicographic_slice_vertex(rows, height: tuple[int, ...], rank: int) -> tuple[int, ...]
def _one_facet_incidence(rows, height: tuple[int, ...], rank: int) -> FacetIncidence
    # P3b

def _setwise_stabilizer_generators(
    generator_permutations: tuple[array, ...], incidence: FacetIncidence, ray_count: int
) -> tuple[array, ...]
    # P4 (Schreier), restricted permutations of the incidence's own positions

@dataclass(frozen=True)
class _FacetOrbitRecord:
    key: FacetMaskKey                  # compact representative incidence
    normal: tuple[int, ...] | None     # exact primitive supporting covector
    members: dict[FacetMaskKey, SchreierWord]

def _orbit_dual_description(
    rows: tuple[tuple[int, ...], ...],
    generator_permutations: tuple[array, ...],
    rank: int,
    *,
    height: tuple[int, ...],
    seed: FacetIncidence | None = None,
    seed_ridge: FacetIncidence | None = None,
    bank: dict | None = None,
) -> tuple[_FacetOrbitRecord, ...]
```

`rows`, `generator_permutations` and the returned keys are the port's existing
types (`cell.rows`, `tuple[array, ...]` from `_facet_mask_data`, `FacetMaskKey`),
so no representation is introduced twice.

### 7.2 Backend method

```python
@cached_method
def facet_orbit_description(self, cell: LorentzianPerfectCell) -> tuple[_FacetOrbitRecord, ...]
```

* reads `cell.rows`, `self.cell_stabilizer(cell)` (already cached),
  `self._affine_height_covector(cell)` — a small new helper that lifts the
  left-kernel computation out of `flip_across` (`:1737-1743`) so the same exact
  `h` is used in three places instead of being recomputed;
* builds the packed ray permutations exactly as `_facet_mask_data`
  (`:1370-1415`) does — the same loop, moved into a helper
  `_cell_ray_permutations(cell)` that both call sites share;
* calls `_orbit_dual_description(...)` with `seed` from
  `self._traversal_seed_incidence` (see 7.3) or `None`;
* validates and returns. Retained per cell in an explicit cache keyed like
  `_normaliz_cell_cones` (so it can be dropped as soon as the traversal has
  consumed it).

### 7.3 Change sites

**`_iter_facet_mask_orbits` (`:1518`)** — one new keyword,

```python
def _iter_facet_mask_orbits(
    facet_keys, permutations, *, retain_words=True, canonical_keys=None,
    seeds_are_representatives: bool = False,
) -> Iterator[list[tuple[FacetMaskKey, SchreierWord] | None]]
```

* `seeds_are_representatives=True` (the new orbit route): the incoming keys are
  *orbit representatives*, so (i) the per-seed "already visited ⇒ skip" branch
  is removed — the first seed always opens a new orbit; (ii) the BFS still uses
  one `visited` set, but a BFS that *reaches* a facet already claimed by a
  previous orbit now raises `ArithmeticError("two facet-orbit representatives
  are stabilizer-equivalent")` instead of silently merging (this is the exact
  partition check); (iii) the final `visited_count != facet_count` assertion
  becomes "every claimed facet was produced by exactly one orbit".
* `seeds_are_representatives=False` keeps today's semantics bit for bit.

**`_facet_mask_data` (`:1330`)** — one new keyword,

```python
facet_source: Literal["cone", "orbit_description"] = "cone"
```

In the `"orbit_description"` branch it does not touch
`_iter_cell_halfspace_pairing_rows` at all: it yields the representative keys
from `facet_orbit_description(cell)` (same `FacetMaskKey` form), returns
`canonical_keys=None`, and reports `seeds_are_representatives=True` to its
caller. Everything else (the packed permutations, the duplicate check on rows,
the cache-dropping of `_configuration_invariant_cache`) is unchanged.

**`_expand_traversal_once` (`:2564`)** — the only structural change:

```python
large = len(source.rows) >= _ORBIT_DESCRIPTION_RAY_THRESHOLD   # new, default 1_000
facet_masks, generator_permutations, canonical = backend._facet_mask_data(
    source, use_owned_cone=False,
    facet_source="orbit_description" if large else "cone",
)
for orbit in backend._iter_facet_mask_orbits(
    facet_masks, generator_permutations,
    retain_words=self._retain_facet_orbits,
    canonical_keys=canonical,
    seeds_are_representatives=large,
):
```

and, when a flip produces a new representative, record its **seed incidence**
(`self._traversal_seed_incidence[id(neighbour)] = shared_incidence`). To get it
exactly, add to `flip_across` a variant that also returns where the flipped
facet's rays sit in the neighbour:

```python
def flip_across_with_incidence(
    self, cell, facet
) -> tuple[LorentzianPerfectCell | None, FacetIncidence | None]
```

returning `facet`'s rows re-indexed in `neighbour.rows`. This is exact and
`O(|facet|)` by a dict lookup; `_kernel_flipping` guarantees the flipped cell
*contains* the critical set (`lorentzian_cells.py:5845` and `:5871`, the
reference's EXIT 1 / EXIT 2 conditions), so the lookup cannot fail; if it ever did, the
helper raises rather than falling back. It is this seed that makes Lemma 8's
`≤ n` LP bound apply to the *whole traversal*: only the first cell of the
complex ever needs a P3b seed.

**`cell_cone` (`:1011`)** — keep both branches, replace the large branch's
*source*:

```python
@cached_method
def cell_cone(self, cell):
    if len(cell.rows) < 10_000:
        return RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
    description = self.facet_orbit_description(cell)          # no engine call
    return _rational_polyhedral_cone(
        cell.lattice,
        self._iter_orbit_covectors(cell, description),         # every orbit member
        equation_covectors=(),                                 # full dimension certified
    )
```

where `_iter_orbit_covectors` transports each representative normal over its
orbit (§5 E3) — the complete H-family is then *reconstructed* from `k₀` exact
normals instead of computed from 49 861 rays. `cell_halfspace_covectors` and
`perfect_domain_traversal`'s `cell_cone(...).halfspace_covectors()` then work
unchanged; `cell_primitive_rays` keeps its Normaliz `ExtremeRays` branch
(`:1155-1190`), which is not on the blocked acceptance path.

### 7.4 What stays exactly as it is

* `cell_cone` for `len(rows) < 10_000` — the small-cell route is untouched
  (`tests/test_large_cell_cone_route.py` pins it), and `_facet_mask_data`'s
  default `facet_source="cone"` keeps the small-cell traversal bit-identical.
* `_iter_facet_mask_orbits` with the default `seeds_are_representatives=False`.
* `flip_across`, `_kernel_flipping`, `initial_cell`, `cell_stabilizer`,
  `_iter_cell_halfspace_pairing_rows`, all acceptance assertions, the 1 800 s
  cap and the 6 144 MiB guard — nothing here is modified.
* No provider (`sage-extensions`) change: everything added is in
  `sage_indefinite_port` and uses only Sage/PPL and the provider's existing
  public interfaces.

---

## 8. Validation plan, and what was measured here

### 8.1 Tests to add (port side)

`tests/test_orbit_facet_enumeration.py`, parametrized over the specimens of
`tests/test_facet_orbit_contracts.py` plus the task's `A2+⟨1⟩+⟨−1⟩`
(rank 4, 27 rays → 6 facets in 4 stabilizer orbits) and the rank-5
`A3+⟨1⟩+⟨1⟩+⟨−1⟩`-style case (89 rays → 24 facets in 6 orbits):

1. **route equality**: `facet_orbit_description(cell)` and the full-cone route
   (`cell_cone(cell).halfspace_covectors()` + `_facet_mask_data(cell)`) give the
   same *set of orbits*: same orbit count, same multiset of orbit sizes, and
   for every full-cone facet, its incidence lies in exactly one orbit route
   orbit;
2. **normals**: every representative normal is primitive, `⟨u, r_i⟩ ≥ 0` for all
   `i`, `rank(Z(u)) = n − 1`, and equals the cone route's covector for that
   incidence up to a **positive rational** scale;
3. **transport**: for every member, its word reproduces the member, the induced
   `M_w` is integral with `|det| = 1` and `M_w^T Gram M_w = Gram`, and
   `u_0·(M_w^{−1})^T` is the member's exact normal;
4. **trivial group**: with `G = 1` the orbit route returns exactly one orbit per
   facet, so the facet *count* agrees with the cone route — the sharpest
   completeness check;
5. **no engine**: monkey-patch `NmzCone`/`NmzResult` to raise and assert the
   description still succeeds (proves the large route is engine-free).

Keep passing, unchanged: `tests/test_large_cell_cone_route.py` (both cone routes
equal; the Normaliz cone's rays are exactly the configuration rows under exact
positive scaling) and `tests/test_facet_orbit_contracts.py` (facet covectors
annihilate exactly a rank-`(n−1)` incidence and are nonnegative elsewhere;
stabilizer generators permute rays and send facets to facets; every facet is in
exactly one orbit).

### 8.2 What was actually run here

The port **cannot be imported in this VM**: `sage_extensions` fails at import
time under the local Sage 10.8 (`sage_extensions/rings/real.py` →
`Modules(ring)(presentation)` → `ValueError: base must be a subring of this
ring`), so the pytest suite above was **not** run here. What was run is the
algorithm itself, on cones of exactly the same mathematical form (a ray matrix
plus integral isometry permutations, each specimen's hypotheses H1–H3 checked
explicitly), against exact complete dual descriptions:

* `logs/validation.log`: **17 specimens, dimension 3–10**, including isotropic
  ray sets (all rays of Gram-norm 0, the port's `ZERO_SQUARES` mode),
  non-signed-permutation symmetry, and the symmetric group omitted entirely.
  For each one the script asserts: every complete-description facet lies in
  exactly one enumerated orbit, the union of orbits *equals* the complete facet
  set (no extras, none missing), every representative normal is exact and
  canonical, and every member's word transports the representative's normal to
  the member's exact normal with a verified integral isometry. All pass.
* `logs/scaling.log`: the scaling table of §6, with the same equality assertion
  up to 2 034 rays / 5 264 facets / dimension 10 — where the orbit route takes
  5.4 s against 42.3 s for the complete description.
* `logs/lp_seed.log`, `logs/pencil.log`: the per-primitive costs of §6,
  including at 49 861 rays in dimension 10.
* The **A7 specimen was not run** (forbidden), and the 44 640-ray synthetic
  configuration used for the `m = 49 861` timings is not it: it is the
  norm-locus configuration `{p ∈ Z⁹ : |p|² = 10} × {1}` of the same rank,
  signature and entry size.

### 8.3 Residual risks

1. **`ρ` (pencils) and `k₀` (facet orbits) for the A7 cell are unmeasured.**
   They are precisely the numbers the blocked Normaliz call would have produced.
   The dim-10 measurement (`ρ = 146`, `k₀ = 3`) and the shape of the cost model
   are the evidence; the guarantee is the one-shot verification run.
2. **The `cell_cone` reconstruction in 7.3** (complete H-family from `k₀`
   transported normals) is exercised by `perfect_domain_traversal`, not by the
   blocked acceptance node, which only needs the traversal's facet orbits. If
   transport over large orbits proves expensive there, that path can keep
   Normaliz; the acceptance path does not touch it.
3. **The seed for the first cell** is the one genuinely LP-shaped step. If the
   first cell's LP is found too slow at 49 861 rays (measured proxies: 10.3 s at
   7 392 constraints, 64.8 s at 22 608 for *ten* sequential LPs), the exact
   alternative is the cutting-plane variant: solve P3b on `n + 2` rays, verify
   the result against all rays in `O(m n)`, and add one violated ray per round;
   termination is guaranteed because each round adds a distinct ray, and the
   final round (all rays present) is P3b itself. This changes the *route* only.
