# Hypothesis, exactness and witness contracts for the orbit-representative facet enumeration

Scope: this is the proof/contract backbone (deliverable 1's data contracts, deliverable 2's
hypothesis discharge, deliverable 3's exactness and witness contract). It does not re-derive the
adjacency-decomposition search itself; it states exactly what that search may assume, proves each
assumption, and states the verification predicate every produced record must satisfy.

Everything here is exact: integers or `fractions.Fraction`/`sage.rings.QQ`. No step is allowed a
floating-point comparison, and §8 proves that none of the primitives needs one.

---

## 1. Objects and data contracts

### 1.1 Ambient space, rows, pairing

* `EXT` (the row matrix) is `cell.rows`: a tuple of `N` integer tuples of length `n`,
  `n = cell.lattice.module_rank()`. The coordinates are those of `lattice.module_generating_set()`
  (`_coordinate_row` is the inverse map). `LorentzianPerfectCell.__post_init__`
  (lorentzian_cells.py:587-625) already enforces: rows are integral, every row's normalized square
  lies in `allowed_squares`, augmented rank `rank((1,v)) = n`, signature `(1,n)` or `(n,1)`.
* The pairing used for *everything* in the facet route is the **coordinate dot product**
  `⟨a, r⟩ = Σ a_i r_i`, not the Gram pairing: `_iter_cell_halfspace_pairing_rows`
  (:1105-1160) evaluates a covector on the coordinate rows and `tests/test_facet_orbit_contracts.py`
  reproduces the values the same way. The Gram matrix enters the cell's *geometry*, never the
  incidence computation. Any new route must keep this pairing; mixing in the Gram pairing silently
  changes the cone.
* `C = cone(R) = {Σ λ_r r : λ_r ≥ 0}` in `V = Q^n`.

### 1.2 The stabilizer action

`cell_stabilizer(cell)` (:1228) returns isometries `g`; each `g` is applied to a *row* by
`row_matrix * action` where `action` is the integer matrix whose rows are `g(e_i)`, i.e.
`g(x) = x * action` for a row vector `x` (`_facet_mask_data`, :1330-1420).
`_facet_mask_data` converts each generator into a **position-image permutation**
`σ_g = (σ_g[0], …, σ_g[N-1])` with `σ_g[i] = ` the position of `g(row_i)`; it raises
`ArithmeticError` if `g` does not permute the rows. Packed `array` typecodes `H/I/Q` by size.
Contract: the orbit machinery may assume `⟨a, g^{-1} r⟩ = ⟨g^{T} a, r⟩` only through the position
permutation: the action on facet incidences is
`σ_g · I = {σ_g[i] : i ∈ I}` and on covectors it is the contragredient action, which the route never
needs to form explicitly (§5).

### 1.3 Facet records

A facet record is a triple

```
(a, I, w)  with  a ∈ Z^n primitive,  I = {r ∈ R : ⟨a,r⟩ = 0},  w a transport word
```

* `a` is the **primitive integer** facet normal (gcd of coordinates 1) with the sign fixed by
  `⟨a, r_j⟩ > 0` for the least index `j ∉ I` (well defined by T1).
* `I` is the **maximal** incidence, as a `frozenset` of positions, or the port's compact
  `FacetMaskKey` (`_facet_complement_key`, :103), which is an exact injective encoding of the same
  set (`_facet_incidence`, :123, is its exact inverse). **Maximality is a contract, not a
  convenience** — see §6.2 and the measured failure in §8.
* `w` is the transport word (Schreier data) carrying the record's facet onto its orbit
  representative; §5 states its contract and how it is verified.

---

## 2. Standing hypotheses

* **(H1) Finite duplicate-free configuration.** `R` is finite with pairwise distinct rows;
  `_facet_mask_data` raises if the rows repeat (`positions` dictionary size check).
* **(H2) Full-dimensionality.** `rank R = n`.
* **(H3) Pointedness.** There is `h ∈ V*` with `⟨h, r⟩ > 0` for every `r ∈ R`.
* **(H4) Group action.** `G ≤ Sym(R)` is induced by integral isometries with `g(C) = C`; the
  generators are exactly what `cell_stabilizer` returns, and `_facet_mask_data` already validates
  that they permute the rows.
* **(H5) Exact arithmetic.** Every value is an integer or a rational; `Q` arithmetic only.
* **(H6) Farkas/exposedness.** Every face of a finitely generated cone is exposed: for a face `F`
  there is `a ∈ V*` with `⟨a,x⟩ ≥ 0` on `C` and `F = {x ∈ C : ⟨a,x⟩ = 0}`. (Standard; the one
  external lemma used below.)
* **(H7) Perfect-cell structure** is used only through (H1)-(H4); nothing below needs the Lorentzian
  form, so the package applies verbatim to any rational pointed cone with a permutation group.

---

## 3. Facet normal uniqueness and canonicity (T1, T2)

**T1.** Assume (H1), (H2), (H6). Let `F` be a facet of `C`. Then

1. `dim F = n-1` and `F = cone(I_F)` with `I_F = {r ∈ R : r ∈ F}`, `rank I_F = n-1`;
2. the set `N(F) = {a ∈ V* : ⟨a,x⟩ ≥ 0 ∀x ∈ C, ⟨a,x⟩ = 0 ∀x ∈ F}` is a **ray** (1-dimensional
   cone), and the valid normal of `F` is unique up to a positive scalar;
3. `span(I_F) = span F` is a rational subspace, so `N(F)` contains a unique-up-to-sign primitive
   integer vector; exactly one of the two signs is valid, and it is the one positive on any ray
   outside `F`. Hence the primitive integer facet normal is **unique and canonical**;
4. `F` is recovered from `a` as `cone({r : ⟨a,r⟩ = 0})`.

*Proof.* (1) `F = cone(I_F)`: `F` is a face, so `F = C ∩ {a = 0}` for some `a` (H6), hence
`F = cone({r : ⟨a,r⟩ = 0}) = cone(I_F)`. Since `dim F = n-1` and `dim cone(S) = rank S` for any
finite `S` (the linear span of `cone(S)` is `span S`), `rank I_F = n-1`, so `span F = span I_F` is
the rational span of integer rows.

(2) Take `a ∈ N(F)` nonzero. Its zero face `Z(a) = {x ∈ C : ⟨a,x⟩ = 0}` is a face containing `F`,
and `Z(a) ≠ C` because `a` vanishes on `C` only if `a` vanishes on `span R = V` (H2), i.e. `a = 0`.
Since `F` is a *maximal proper* face, `Z(a) = F`. So `a` vanishes on `span F`, a hyperplane: `a` lies
in the 1-dimensional `span(F)^⊥`. Hence `N(F)` is contained in a line; it is a cone, so it is `{0}`,
a ray, or a line. It is not `{0}` (H6 gives a generator) and it is not a line: the two elements of a
line are `±a`, and validity on any `r ∉ I_F` (which exists since `rank I_F = n-1 < n`) forces
`⟨a,r⟩ > 0` for exactly one sign. So `N(F)` is a ray.

(3) `span(F)^⊥` is the annihilator of a rational subspace, so it is a rational line and contains a
primitive integer vector `a_0`, unique up to sign; the valid sign is the one positive on a ray
outside `F`. (4) is (1) read backwards. ∎

**T2 (facet ⟺ certified zero set).** Assume (H1), (H2). Let `a ∈ Z^n` be primitive with
`⟨a,r⟩ ≥ 0` for all `r ∈ R`, and let `Z(a) = cone({r : ⟨a,r⟩ = 0})`. Then `Z(a)` is a facet of `C`
**iff** `rank{r : ⟨a,r⟩ = 0} = n-1`; and in that case `a` is `Z(a)`'s canonical normal (T1).

*Proof.* `Z(a)` is a face (exposed by `a`), and `dim Z(a) = rank` of its ray set. A face of a
full-dimensional cone is a facet iff its dimension is `n-1`. For the last claim, both `a` and the
canonical normal of `Z(a)` are valid primitive vectors in the 1-dimensional `span Z(a)^⊥` with the
same sign, hence equal. ∎

**Corollary (facet identity).** Two certified facet records `(a,I)`, `(a',I')` describe the same
facet iff `I = I'` iff `a = a'`. Therefore the incidence set (or its compact key) is a sound facet
key, and `a` is a sound canonical form. This is what makes the "already found?" test of the search
exact and hashable.

---

## 4. The recursive object: incidence cones and ridges (T3, T4)

Let `F` be a facet, `a_F` its canonical normal, `I_F` its incidence, `C_F = cone(I_F) ⊆ span F`.

**T3 (the facet cone).** Assume (H1)-(H3).

1. `C_F` is full-dimensional in `span F` (`dim C_F = rank I_F = n-1 = dim span F`) and pointed as a
   subcone of the pointed cone `C`. Hence `C_F` obeys (H1)-(H3) *in dimension `n-1`*: the recursion's
   hypothesis is inherited.
2. The facets of `C_F` in `span F` are exactly the `(n-2)`-faces of `C` that are contained in `F`.
   Each is a **ridge** of `C`.
3. (Lift lemma) If `b` is a functional on `span F` that is valid on `I_F`, put
   `b' = b + M a_F` with
   `M ≥ M_0 := max{ -⟨b,r⟩ / ⟨a_F,r⟩ : r ∈ R \ I_F, ⟨b,r⟩ < 0 }` (a finite rational, `0` if the set
   is empty). Then `b'` is valid on all of `R` and `⟨b',r⟩ = ⟨b,r⟩` for `r ∈ I_F`; in particular
   `cone({r ∈ R : ⟨b',r⟩ = 0}) = cone({r ∈ I_F : ⟨b,r⟩ = 0})` **exactly**. Every face of `C_F` is
   therefore a face of `C`.

*Proof.* (1) is immediate from §3(1) plus `C_F ⊆ C`, and subcones of pointed cones are pointed
(a lineality vector of `C_F` would be one of `C`). (2) A subset `G ⊆ F` is a face of `C` iff it is a
face of `C_F`: one direction is `G = C ∩ {a = 0}` intersected with `span F`; the other is (3). Faces
of `C_F` of dimension `n-2 = dim C_F - 1` are exactly its facets. (3) For `r ∉ I_F`:
`⟨a_F, r⟩ > 0` by maximality of `I_F` (a facet normal is strictly positive off its zero set, T1), so
`⟨b', r⟩ ≥ ⟨b,r⟩ + M_0⟨a_F,r⟩ ≥ 0`; strict unless `⟨b,r⟩ = -M_0⟨a_F,r⟩` and the max is attained —
enlarge `M` by `1` in the definition of `M_0` to make it strict on `R \ I_F`. On `I_F`, `a_F`
vanishes. ∎

**T4 (every ridge lies in exactly two facets; its normal cone is 2-dimensional).** Assume (H1)-(H3),
(H6). Let `G` be a ridge of `C` (a face of dimension `n-2`) with incidence
`I_G = {r ∈ R : r ∈ G}` and

    N_G = { x ∈ V* : ⟨x, r⟩ = 0 ∀ r ∈ I_G,  ⟨x, r⟩ ≥ 0 ∀ r ∈ R }.

Then

1. `rank I_G = n-2` and `N_G` is a **2-dimensional pointed cone** (i.e. exactly two extreme rays);
2. `G` is contained in exactly two facets `F`, `F'`, and the extreme rays of `N_G` are `a_F` and
   `a_{F'}`; in particular `F ∩ F' = G` (the ridge is *doubly covered*);
3. given `I_G`, both extreme rays are computed exactly: take a basis `(u_1,u_2)` of `(span I_G)^⊥`,
   form the images `(⟨u_1,r⟩, ⟨u_2,r⟩)` for `r ∈ R`, and intersect the resulting halfplanes; the two
   boundary directions are the two facet normals (as covectors, up to the positive scale used);
4. `a_F`'s image is the boundary ray that `F` occupies; the other boundary ray is the *adjacent
   facet*'s normal. There is no third choice and no sign ambiguity.

*Proof.* (1) `rank I_G = n-2` by §3(1) applied to `G` (or directly: `dim cone(I_G) = rank I_G`).
`(span I_G)^⊥` is `2`-dimensional. `N_G` is that 2-plane intersected with the halfspaces
`⟨·,r⟩ ≥ 0` (`r ∈ R`). By T3(3), the lift `b + M a_F` of any facet normal `b` of `C_F` through `G` is
strictly positive on `R \ I_G`, so `N_G` has nonempty interior in the 2-plane; it is pointed because
it is contained in the 2-plane and `N_G ∩ (-N_G) = {0}` (a lineality vector would be valid and also
negatively valid on `C`, hence vanish on `span R = V`). A pointed 2-dimensional cone has exactly two
extreme rays. (2) `N_G` is the normal cone of the face `G`; its extreme rays expose faces
`Z(e) ⊇ G`. If `x ∈ N_G` is in the *relative interior* then `Z(x) = G` (any strictly interior point of
a normal cone exposes the face), so an extreme ray `e` exposes a face `Z(e) ⊋ G`; being a face it has
dimension `≥ n-1`, and `< n` (it is proper, as `e ≠ 0`), so `dim Z(e) = n-1`: `Z(e)` is a facet. Two
extreme rays, and `Z(a_F) = F`, `Z(a_{F'}) = F'` are two distinct facets containing `G`. Conversely
any facet `H ⊇ G` has `a_H ∈ N_G` with `Z(a_H) = H ≠ G`; a normal-cone element exposing more than `G`
cannot be interior, so it lies on the boundary, and in a 2-dimensional cone every boundary point lies
in one of the two extreme rays. Hence exactly two facets, they are `F` and `F'`, and their normal rays
are the two extreme rays. Finally `F ∩ F' ⊇ G`, and `F ∩ F'` is a face of dimension `≤ n-2` (a face
containing `G` properly would have dimension `≥ n-1` and would hence be a facet — but facets containing
`G` are only `F` and `F'`, and `F ≠ F'`), so `F ∩ F' = G`. (3)-(4) restate (2) in coordinates: a
boundary ray of the intersection of halfplanes satisfies one constraint with equality, i.e. it is
`±rot90(c)` for one constraint vector `c`, and (1) says there are exactly two. ∎

**T4 consequence (the algorithm may not skip ridges).** Doubly-coveredness is exactly what makes
"the search is closed" a completeness certificate (§5). It is also exactly what fails without
pointedness: for `R = {(1,0), (-1,0), (0,1)}` in `Q^2` the cone is the closed half-plane
`{y ≥ 0}`, it has exactly one facet, and the zero face lies in *one* facet, not two — measured in
§8. This is why (H3) may not be dropped, and why (H3) may not be inferred from `rank R = n` (§6.1).

---

## 5. Completeness, transport, and the algebra of the search (T5, T6, T7)

**T5 (closure certificate).** Assume (H1)-(H3), (H6). Let `S` be a nonempty set of facets of `C`
such that

> for every `F ∈ S` and every ridge `G` of `F`, the adjacent facet `F'` (T4) lies in `S`.

Then `S` is the set of all facets of `C`.

*Proof.* Suppose `H` is a facet not in `S`. Pick `x` in the relative interior of `H` and `y` in the
relative interior of any `F_0 ∈ S`. Perturb `x` and `y` within those relative interiors so the
segment `[x,y]` meets no face of `C` of dimension `≤ n-3` (finitely many closed sets of measure
zero, so such a choice exists). Let `F_1, …, F_k` be the facets met by the segment in order
(`F_1 = H`, `F_k = F_0`). Consecutive facets `F_i, F_{i+1}` meet in a face that the transversal
segment crosses; the local cone of `C` at a transversal crossing point has exactly the two local
facets `F_i, F_{i+1}`, so `dim(F_i ∩ F_{i+1}) = n-2`: they share a ridge. Now `F_k ∈ S` and
`S` is closed under adjacency *along ridges*, so `F_{k-1} ∈ S`, …, `F_1 = H ∈ S`. Contradiction. ∎

**T5' (group version, the form the search actually uses).** Assume (H1)-(H4), (H6). Let
`S` be a nonempty `G`-invariant set of facets (`G·S = S`). Suppose that for every `G`-orbit
`O ⊆ S`, for every `F ∈ O` and every ridge `G_r` of `F`, the `G`-orbit of the adjacent facet lies
in `S`. Then `S` is the set of all facets.

*Proof.* `G` preserves `C`, so it maps faces to faces and facets to facets of equal dimension, and it
maps ridges to ridges. Hence `S` is closed under adjacency in the sense of T5: for `F ∈ S` and a ridge
`G_r`, the adjacent facet `F'` lies in the orbit of the adjacent facet computed from any
`G`-translate of `F`, which is in `S` by hypothesis, and `S` is `G`-invariant. Apply T5. ∎

*Consequence for the search.* One representative per orbit, *with all its ridges*, suffices:
adjacency for the other members follows by symmetry, and closure then certifies completeness. The
"all its ridges" clause is not removable: enumerating only some ridges leaves a strictly weaker
statement (T5's hypothesis is about every ridge).

**T6 (facet stabilizer = incidence set stabilizer).** Assume (H1)-(H4). For `g ∈ G` with position
permutation `σ_g` and a facet `F` with maximal incidence `I_F`:

    g(F) = F   ⟺   σ_g(I_F) = I_F .

Hence `Stab_G(F)` = the set stabilizer of the position set `I_F` = the subgroup computed by
`PermutationGroup(...).stabilizer(...)` on the ray action, and `|G·F| = |G| / |Stab_G(F)|`
(orbit-stabilizer).

*Proof.* `g` is a linear automorphism with `g(C) = C`, so it permutes the faces of `C`; the faces of
`C` are exactly `cone(I)` for incidence sets `I`, and `g(cone(I)) = cone(σ_g(I))`. Faces are equal iff
their maximal incidences are equal (`F = cone(I_F)`, T1(4)+(1)), so `g(F) = F ⟺ σ_g(I_F) = I_F`. ∎

**T7 (transport/witness contract).** With `G = ⟨g_1,…,g_m⟩` as position permutations and
`σ_w` the composite of a word `w = (i_1,…,i_k)` (`σ_w = σ_{i_1} ∘ … ∘ σ_{i_k}`, applied to
positions left to right; this is the order `_schreier_generator_indices`
(lorentzian_cells.py:165-200) yields — it collects the nested word from the innermost entry outward
and reverses, i.e. yields generator indices in application order).

A record `(a, I, w)` for the orbit representative `F_rep = (a_rep, I_rep)` of `F` is a **valid
transport witness** iff all of the following hold, each by exact arithmetic:

* **P1 (canonicity of the normal).** `a` is primitive, `⟨a,r⟩ ≥ 0` for all `r ∈ R`,
  `rank{r : ⟨a,r⟩ = 0} = n-1`, and `a`'s sign obeys T1(3).
* **P2 (incidence).** `I = {i : ⟨a, row_i⟩ = 0}` exactly, as sets of positions.
* **P3 (ridge completeness of the record's facet).** the ridge list of `I` produced by T4(3) is
  complete: every facet of `C_F` computed by the recursive call at dimension `n-1` is present.
* **P4 (transport).** recomputing `σ_w` over the packed arrays gives `σ_w(I) = I_rep` as position
  sets; and, as a consequence, the transported normal `σ_w·a` equals `λ a_rep` with `λ > 0` — a
  relation that is *checkable* exactly (`σ_w·a` vs `a_rep` proportional with positive factor), and
  that must be checked rather than inferred, because it also certifies that the word was composed in
  the documented order.
* **P5 (orbit membership).** membership of a found facet in an already-known orbit is decided by
  `∃ g ∈ G : σ_g(I) = I_rep`; a decision procedure may use the Schreier transversal of `⟨G, I_rep⟩`
  (i.e. `σ_g` for `g` in a transversal of `Stab_G(F_rep)`), never a numerical comparison of normals
  up to tolerance.

**Lemma (witness existence).** If `F' = g(F)` for some `g ∈ G` then `g(I_F) = I_{F'}` and `g·a_F` is
a positive multiple of the canonical normal `a_{F'}`; hence a word `w` with `σ_w(I_F) = I_{F'}` exists
by Schreier's lemma, and P3-P4 hold for it. ∎

---

## 6. Two corrections to the recorded reasoning

### 6.1 `rank R = n` does **not** imply "no lineality"

The A7 handoff infers pointedness from "row rank 10 = ambient dimension (so the cone is
full-dimensional, no lineality)". Full-dimensionality follows; pointedness does not. Counterexample:
`R = {(1,0), (-1,0), (0,1))} ⊂ Q^2` has `rank = 2 = n` and `C = {y ≥ 0}` with lineality `{(x,0)}`.
The augmented-rank axiom checked in `__post_init__` does not give it either (`(1,1,0)`, `(1,-1,0)`,
`(1,0,1)` has rank `2 = n`).

What *is* needed is (H3), and it is cheap to certify in two ways:

1. **a-priori:** exhibit one `p ∈ Q^n` with `⟨p, r⟩ > 0` for every row (the reference's central
   vector plays this role; the port has `_positive_direction` and `_search_initial_vectors`);
2. **a-posteriori, and self-validating:** once facet normals are known,
   `p := Σ_{facets F} a_F` satisfies `⟨p,r⟩ > 0` for every ray (a ray lies in some facet, and it
   cannot lie in all of them). This certifies pointedness *and* the covering property from the
   produced data, and it costs one `O(f_orb · |G|)` accumulation over the orbit data.

Without (H3) the doubly-covered property fails (§4 consequence), and the closure certificate T5 has
no content.

### 6.2 Facet identity must be keyed by the *maximal* incidence

For a facet `F` and a proper subset `I ⊊ I_F` with `span I = span I_F`, `cone(I) = F` but
`Stab_G(I) ⊆ Stab_G(F)` can be strict, and the same facet then appears under several keys with
different stabilizers, so orbit bookkeeping is inconsistent. Measured on the cube cone in `Q^4`
(8 rays, 6 facets, `|I_F| = 4 = n`): `|Stab_G(F)| = 2` while a 3-of-4 spanning sub-incidence has set
stabilizer `1` (§8, C-D). The port's incidences are maximal by construction
(`_iter_cell_halfspace_pairing_rows` tests every label), and the contract must survive any
optimization that stores "a spanning subset": it may not.

---

## 7. Exactness of the primitives (T9, T10)

**T9 (no rounding).** Under (H5), every primitive below is an exact sign test or an exact solve, so
no output can be affected by rounding, and every output has an integer certificate checkable by
P1-P5.

| primitive | exact algorithm | certificate |
| --- | --- | --- |
| `rank S` | Gaussian elimination over `Q` (or fraction-free on `Z`) | echelon form |
| `(span I)^⊥` basis | exact kernel of integer rows | `rows · basis = 0`, `rank` count |
| primitive integer normal | clear denominators by `lcm`, divide by `gcd`; sign per T1(3) | normals of `a` = 1, T1 sign |
| validity `⟨a,r⟩ ≥ 0`, incidence | integer dot products | the list of values |
| adjacency across a ridge | T4(3): exact 2-dimensional halfplane intersection | the two normal rays satisfy all constraints |
| orbit membership / transport | permutation composition on packed arrays | P4 recomputation |

**T10 (the 2-dimensional primitive).** Let `K = {t ∈ Q² : ⟨c_i, t⟩ ≥ 0, i = 1..m}` with nonempty
interior (guaranteed by T3(3): the lift of a facet normal of `C_F` is interior to `N_G`).

1. *Correctness.* The extreme rays of `K` are boundary rays; a boundary ray satisfies some `c_i`
   with equality, and `⟨c_i, t⟩ = 0` in `Q²` forces `t ∝ ±rot90(c_i)`. Hence
   `ext(K) ⊆ {±rot90(c_i) : i} ∩ K`, and conversely every element of that set that satisfies all
   constraints with equality on some `i` is extreme. Since `K` is pointed and 2-dimensional,
   `|ext(K)| = 2` exactly (T4(1)).
2. *No floats.* All tests are exact sign evaluations of `c_i · t` and `2×2` determinants, so no
   `atan2`, no `float`, no tolerance appears. Sorting by angle is never needed.
3. *Cost.* A naive candidate sweep is `O(m²)`. Linear time suffices and should be specified: an
   interior point `t_0` is known (the lifted normal), and a 2-dimensional halfplane intersection
   from an interior point maintains the two current boundary rays, updating them in `O(1)` per
   constraint, so `{ext(K)}` is computed in `O(m) = O(N)` exact operations. Per ridge the cost is
   therefore `O(N)` exact integer operations, independent of the facet count.

---

## 8. What was verified exactly here, and what could not be run

Environment: `python3` + Sage (`/opt/math`) is present, but the port could not be imported in this VM:
`sage_extensions` (vendored, read-only) fails at import with
`ValueError: base must be a subring of this ring` (`sage_extensions/rings/real.py:1078`, reached from
`ExactRealField()` at module import). No provider patch is permitted, so no specimen could be run in
this VM. The package is instead verified on self-contained exact cones that exercise exactly the
structural claims (`verify_toy.py`, `verify_ridges.py`, `fractions.Fraction` only):

| cone | rays | n | facets | `|I_F|` | ridges | ridge multiplicity | T4 extreme rays |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hexagonal cone (`z = 1`) | 6 | 3 | 6 | 2 = n-1 | 6 | every ridge in exactly 2 facets | 2 = the two facets |
| cube cone (`(±1,±1,±1,1)`) | 8 | 4 | 6 | 4 = n | 12 | every ridge in exactly 2 facets | 2 = the two facets |

* `min over rows of (sum of facet normals) = 6 > 0` in both cases: the a-posteriori pointedness
  certificate of §6.1 works and is exact.
* Ridge structure: for every facet, the ridges were recomputed by changing coordinates into
  `span F` (RREF, coordinates read off pivot columns), enumerating the facets of `C_F` **inside
  `span F`** (this is where the coordinate change is essential: in the ambient annihilator a facet
  normal of a 2-dimensional `C_F` is only determined modulo `(span F)^⊥`), then lifting with
  `b + M a_F` as in T3(3). The lifted normals' zero sets have rank `n-2` and every ridge is doubly
  covered.
* T4(3): for each ridge, the ambient normal cone was intersected in a 2-dimensional basis of
  `(span I_G)^⊥`; it had exactly two extreme rays, and their zero sets were exactly the two facets
  through the ridge. No float was used anywhere in these computations.
* C-C (lineality, `R = {(1,0),(-1,0),(0,1)}`): `1` facet, `min(sum of normals) = 0` — no strictly
  positive functional exists — and the zero face lies in one facet, not two: the failure mode of §4.
* C-D (cube cone, `G = ⟨swap x,y; flip z⟩`): `|Stab_G(F)| = 2` for the facet `x = 1`, while a 3-of-4
  spanning sub-incidence of it has set stabilizer `1` (§6.2).

Files: `.ultramath/smvzd6hvn/m6/verify_toy.py`, `verify_ridges.py` (runnable with plain `python3`),
`verify_contracts.py` (the port-based version of the same checks; needs a working `sage_extensions`).

---

## 9. Obligations this imposes on the implementation

Any implementation of the orbit-representative route must satisfy, and should assert in tests:

1. rows duplicate-free, integer; `rank rows = n`; and one explicit `p` with `⟨p,r⟩ > 0` (§6.1);
2. every produced record satisfies P1-P5; `a` primitive and canonically signed;
3. incidences maximal; facet keys are the compact keys of §1.3; facet identity by §3 Corollary;
4. per representative facet, the ridge list is complete in the sense of T3(2) (computed on
   `C_F` *inside* `span F`, with the coordinate change), not a partial list — otherwise T5 does not
   apply;
5. `Stab_G(F)` computed as the set stabilizer of `I_F` under the ray permutation action (T6), never
   as a covector stabilizer with tolerance;
6. the closure test `S` is `G`-invariant and closed under ridge adjacency (T5'), and that closure is
   reported as the completeness certificate, together with the a-posteriori positivity certificate
   of §6.1;
7. the only per-orbit size objects retained are `(a, I, w)` records: `f_orb` records, never the `f`
   facet normals that the cancelled Normaliz route materializes.

---

## 10. The construction these contracts pin (statement, not the search's design)

The contract layer is agnostic about *how* the facet graph is walked, but it fixes the data flow,
the order of operations and the certificates, so that any walk satisfying §9 is sound:

**Phase 0 — admission (cheap, `O(N n)`).** `EXT = cell.rows`; assert rows integer and duplicate-free;
`rank EXT = n`; exhibit/verify `p` with `⟨p,r⟩ > 0` for every row (§6.1); build the position
permutations `σ_{g_1},…,σ_{g_m}` of `cell_stabilizer(cell)` by `_facet_mask_data`'s convention
(§1.2), and the permutation group `G = ⟨σ_i⟩` on positions `0..N-1`.

**Phase 1 — seed.** One facet record `(a_0, I_0)` satisfying P1. The source's own route produces a
supporting covector on the way (`initial_cell` maintains `base_normal` and updates it from
`_kernel_flipping`'s `flipped_normal`, in the augmented form `[-scalar, gram·central]`), but that is a
*supporting hyperplane of the cell data*, not automatically a homogeneous facet normal of
`C = cone(cell.rows)`; the seed must be extracted from it (or from any other source) and then
**re-verified** by P1: primitive, nonnegative on every row, zero set of rank `n-1`, canonical sign.
`I_0` is by construction the maximal zero set. Insert the `G`-orbit of `I_0` into the worklist as
*undone*.

**Phase 2 — worklist over facet orbits.** While an undone orbit exists, take its representative
`F = (a_F, I_F)`:

1. `H := Stab_G(F)` = set stabilizer of `I_F` under the position action (T6).
2. Recursion: `ridges(F) := facet_orbit_reps(I_F, H, n-1)` — the same routine one dimension down,
   applied to the coordinate image of `I_F` in `span F` (§4(iii)) so that `cone(I_F)` is
   full-dimensional there (T3(1)).
3. For each ridge orbit representative `G_r` (incidence `I_{G_r} ⊆ I_F`, `rank = n-2`): compute the
   two extreme rays of `N_{G_r}` exactly (T4(3), `O(N)` exact operations); one is `a_F`'s ray, the
   other is the adjacent facet's normal `a'` with incidence `I'`. If `I'`'s orbit is not yet in the
   worklist, insert it as undone and record the transporter `g` with `σ_g(I') = I_F`'s
   representative-normalized form, or `σ_g(I_F) = I'`-form, per the chosen transversal convention,
   together with its Schreier word and the P4 verification.
4. Mark the orbit done.

Base cases: `n = 2` — a pointed 2-dimensional cone's facets are its two extreme rays, computed by the
2-dimensional primitive of T10; `n ≤ 1` — no facets. (Each step lowers the rank of the working
configuration from `n` to `n-1`, so the recursion depth is `≤ n-2`.)

**Phase 3 — certificate and output.** Emit

* one record `(a, I, w)` per facet orbit, P1-P5 verified;
* the closure certificate: the search terminated because the undone set emptied, i.e. `S` (the union
  of the found orbits) is nonempty, `G`-invariant and closed under ridge adjacency — by T5' this is
  exactly `S = all facets`;
* the pointedness/covering certificate `p := Σ_{F ∈ S} a_F` (sum over all orbit elements, i.e. weighted
  by `|G|/|Stab|`) with `⟨p, r⟩ > 0` on every row, recomputed exactly.

**Record count.** `f_orb` records plus, per orbit, its ridge list. Nothing in Phases 0-3 stores or
touches the `f = Σ_O |G|/|Stab_G(F_O)|` individual facet normals; consumers who need the whole orbit
apply the stored permutations on demand. That is the sense in which the route never materializes the
complete H-description that the cancelled Normalizer call at `lorentzian_cells.py:1051` produces.

**Interface obligations towards the existing port.** `cell_cone` keeps its current behaviour for
cells where it is already correct (`len(cell.rows) < 10_000`, the `from_rays`/Normaliz branch), and
the orbit route is selected for the large cells; both routes must agree on bounded specimens
(`tests/test_large_cell_cone_route.py`), and the four contract properties of
`tests/test_facet_orbit_contracts.py` must keep holding. Those acceptance assertions, the 1800 s cap
and the 6144 MiB guard are untouched by anything in this document.
