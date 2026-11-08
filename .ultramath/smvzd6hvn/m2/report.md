# Orbit-representative facet enumeration for perfect-cell cones

Reconstruction of the reference route, its proof, its exactness layer and its port seam.

Status of every claim in this document is marked as one of:

* **[source]** — read directly in the vendored tree (file:line quoted).
* **[proof]** — proved here from the stated hypotheses.
* **[measured]** — executed in this VM on synthetic specimens (numbers given); the
  artefact is `orbit_facets_proto.py`, results in `results_full.txt`,
  `results_balinski_only.txt`, `results_cutoff.txt`.
* **[unverified]** — depends on something I could not run here (see §0).

---

## 0. What could not be verified, and why

**The port cannot be imported or executed in this VM.** Sage here is 10.8 / py3.12
(`/opt/math/bin/sage`). `sage_extensions` calls `base_ring().free_module(1)`
(`categories/modules/native_modules.py:178`); Sage 10.8 requires the explicit base
and raises `ValueError: base must be a subring of this ring`:

```
$ sage -c "print(ZZ.free_module(1))"        # ValueError
$ sage -c "print(ZZ.free_module(ZZ,1))"     # works
```

`lorentzian_cells.py:19 from sage_extensions.all import ...` therefore fails before
any test runs. Every checkout present in `$HOME` fails identically
(`sage-extensions-{42d7361c,565fc6b9,69907aef,9dc89ac,c546077d,cdafca84,e83c4d6}`),
with `PYTHONPATH=src:<ext>/src` and the repo `.venv` (which does carry
`tree_sitter`/`sageparse`). `pyproject.toml` wants `python>=3.14`; the VM has 3.12.

Consequences, stated plainly:

* No number below is a measurement of the A2 / A3+<1>+<-1> / A7 specimens. The A7
  figures used anywhere are the task statement's (49861 rows, rank 10, max entry 15,
  >1800 s), not mine.
* `tests/test_large_cell_cone_route.py` and `tests/test_facet_orbit_contracts.py`
  were **not** run. Whether they still pass after the change in §6 is
  **[unverified]**; §7 says exactly what to run.
* What *is* executable is Sage itself, so the algorithm below is validated on
  synthetic symmetric cones where the complete facet list comes from an independent
  engine (Sage's polyhedron backend) — §7.

---

## 1. What the reference actually does [source]

`lorentzian_perfect.h:1433-1472` (`f_adj`) is the contract:

```
ResultStabilizer           res_stab = LORENTZ_ComputeStabilizer(data.LorMat, x.EXT, os);
vectface                   TheOutput = DualDescriptionRecordFullDim(EXT_T, res_stab.GRPperm, data.rddo);
x.GRP = res_stab.GRPperm;
for (auto &eInc : TheOutput)  ListAdj.push_back({eInc, LORENTZ_GetEXT(LORENTZ_DoFlipping(..., eInc, ...))});
```

So: input = the row matrix of the cell's rays (`EXT_T`, `T` = rationals) plus a
permutation group acting on the **rows**; output = a `vectface`, i.e. a list of
incidences (sets of row indices), used one-for-one to build the neighbouring cells.
Nothing else is consumed. One incidence per facet orbit is therefore *sufficient*,
and `lorentzian_perfect.h:1263,1452` show the `RecordDualDescOperation rddo` is a
memo object threaded through the calls, not an input description.

`DualDescriptionRecordFullDim`, `DualDescriptionStandard` (`edgewalk.h:1726`) and
`RecordDualDescOperation` are **not defined** in the vendored tree. The definition
that *is* vendored is the GAP ADM,
`Indefinite.jl@374a5ebb/indef/lib/DualDescriptionAndADM.g:77-188`, and it is the
same algorithm; the following are its load-bearing details, all read from source.

### 1.1 The object being enumerated is an incidence set, not a normal

A facet is represented by its **incidence**: the set of row indices whose rays it
annihilates (`DualDescriptionAndADM.g:144` `eInc := RPL.FuncRecord(SelectedOrbit)`,
`:155` `Ladj := __ListFacetByAdjacencyDecompositionMethod(EXT{eInc}, RedStab, ...)`).
`EXT{eInc}` is the row submatrix on the incidence. Nothing in the control flow ever
stores a normal.

### 1.2 The recursion is on the incidence's own row set, with the set stabilizer

`DualDescriptionAndADM.g:148`: `RedStab := Data.GroupFormalism.Stabilizer(EXT, WorkingSymGroup, eInc)`
— the **setwise** stabilizer of the incidence inside the working symmetry group;
`:155` recurses on the subconfiguration with that stabilizer. The returned
sub-incidences are ridges of the representative, and `:159-162` lifts each of them:

```
for fInc in Ladj do RPL.FuncInsert(RPLift.FuncLift(fInc)); od;
RPL.FuncPutOrbitAsDone(SelectedOrbit);
```

So one lift per **orbit of ridges**, and the lifted facets are inserted into the
parent's orbit registry; the parent does its own orbit bookkeeping.

### 1.3 The lift is a 2-dimensional exact computation, not an LP and not a bound

`Fundamental.g:605-690` `__ProjectionLiftingFramework_Rational(EXT, OneInc)`. With
rows projected to a coordinate frame (`:518 __ProjectionFrame`, `ColumnReduction` to
`RankMat` columns), the lift of a ridge `TheInc ⊆ OneInc` is:

```
LINCR  := EXTproj2{TheInc};                     # ridge rows
VMA    := NullspaceMat(TransposedMat(LINCR));   # annihilator of the ridge
if Length(VMA) <> 2 then Error("incoherence") fi;
for eExt in EXTproj do VO := [VMA[1]*eExt, VMA[2]*eExt];  # 2-dim projection of every row
...  angular sweep over the VO to get the two extreme directions LV ...
S[iV] := eV[1]*VMA[1] + eV[2]*VMA[2];          # back to the ambient space
VCE   := S[Position(test, 0)];                  # the candidate that does NOT kill all of OneInc
return Filtered([1..Length(EXTproj)], x -> EXTproj[x]*VCE = 0);
```

Three things matter for the port. (i) The ridge's annihilator is computed by exact
nullspace and must be **2-dimensional** — that is the only correctness check the
source performs, and it is a check. (ii) The two candidate facets containing the
ridge are the two extreme rays of the cone `{y in annihilator : y . row >= 0}`; the
parent facet is one of them (`test[iV]=1`, it annihilates all of `OneInc`), the
target is the other. (iii) The result is the **full incidence by an exact
`Filtered` over all rows**, so no sign, scaling or feasibility assumption enters.
No linear programme appears anywhere in this path; `Data.GetInitialRays` does use an
LP (`LinearProgramming.g:20-31`) but only to *seed*.

### 1.4 The registry counts elements, and the stopping rule is Balinski's bound

`SetFunctionality.g:204-330` `OrbitGroupFormalism`: `FuncInsert` deduplicates by a
group invariant plus an isomorphy test (for the standard formalism the invariant is
`Minimum(Orbit(group, inc, OnSets))`, `:117-121`); `OrbSize := GroupOrder/Order(Stabilizer)`
(`:295-303`) and **`nbUndone` accumulates `OrbSize`, i.e. it counts facets, not
orbits** (`:283`, `:307`).

The loop's exit test, `DualDescriptionAndADM.g:131-142`, is

```
eSetUndone := RPL.ComputeIntersectionUndone();
nbUndone   := RPL.FuncNumberUndone();
if RPL.FuncNumberOrbitDone() > 0 then
   if nbUndone <= TheDim-1  or  Length(eSetUndone) > 0 then break; fi;
fi;
```

with `TheDim := RankMat(EXT)-1` (`:114`). For a full-dimensional cone of ambient
dimension `n = rank`, `TheDim = n-1`, so the test is **`#unexpanded facets <= n-2`**.
That is exactly Balinski's connectivity of the facet graph, and it is sound — §3.2.
The name of the vestigial local `RbalinskiRank` at `SetFunctionality.g:205` is the
same argument.

The second disjunct (`all unexpanded facets share a full group orbit of rows`,
`OrbitIntersection` at `SetFunctionality.g:2`) fires often in my runs (§7) and never
lost a facet there, but **I found no proof of it**, so §6 does not port it.

### 1.5 Three accelerators in the source, only one of which the port should take

* `BankRecording` (`DualDescriptionAndADM.g:12-74`, used at `:79-82` and `:177-185`):
  a memo of completed sub-descriptions keyed by a group invariant plus an exact
  isomorphy test. **Sound, and worth porting** — it turns the recursion tree into a
  DAG. Measured effect in §5.
* `TestNeedMoreSymmetry` (`:84-100`): when `Length(EXT) > RankMat(EXT)+4` the working
  group is enlarged by the subconfiguration's own automorphism group. Sound **only
  because** `:173` refines the answer back to the caller's group with
  `LiftingOrbits` = `GlobalLiftingOrbitsOnSets` (`SetFunctionality.g:95-108`) via
  `DoubleCosets(BigGroup, Stab, SmallGroup)`. Optional for the port; harmless if the
  refinement is kept, wrong if it is dropped.
* `IsRespawn` → `DualDescriptionFunction` (`Lorentzian.g:1268-1282,1299-1308`;
  `DualDescriptionSimple.g:1-8` = `POLY_dual_description_group(method="lrs_ring", ...)`):
  when the group is large (`>=50000` and depth `<=2`, or `>=100` with more than
  `rank+20` rows and depth not 2) the ADM is abandoned for a **complete dual
  description**. This is the hidden full materialization. **Do not port it**: it is
  exactly the `SupportHyperplanes`-class cost the port must avoid on the A7 cell, and
  the port has no fallback budget.

---

## 2. The construction, with contracts

Everything below is stated for one cell; the port's recursion applies the same
procedure to subconfigurations.

**Input.**

* `rows = cell.rows`: `m` integer tuples of length `n`, duplicate-free, `rank = n`
  (full dimension), and the cone `C = cone(rows)` **pointed** (no lineality).
  For the A7 cell this is the task's 49861 x 10 matrix.
* `permutations`: one packed `array` per generator of a finite group `G` acting on
  the **positions** of `rows` — exactly what `_facet_mask_data` already produces
  (`lorentzian_cells.py:1330-1408`); each generator comes from an exact lattice
  isometry preserving the ray set (`cell_stabilizer`, `:1228`). `G` stabilizes `C`.

**Output.** A generator of **facet orbits**: for each orbit, a list
`[(mask, word), ...]` covering the orbit exactly once, where `mask` is the port's
facet key (incidence or its complement, `lorentzian_cells.py:123-163`) and `word` is
the exact Schreier word carrying the orbit representative to that member. One
representative per orbit plus its word is all the traversal consumes
(`_expand_traversal_once`, and the side-pairing transport at `:2963-2999`).

**Procedure.**

```
EnumerateFacetOrbits(rows, G):
  n := rank(rows)                                     # = len(rows[0]) after projection
  if n <= 2: return the 2 extreme rays (n=2) as the facets        # base case
  P    := projection of rows to pivot columns of the row space     # exact, injective on span
  seed := OneFacet(P)                                             # exact, any facet
  registry := { orbit-canonical(seed): {incidence: seed, done: false} }
  while exists an unexpanded orbit O:
      if  (#unexpanded facets <= n-2) and (at least one orbit expanded):      break   # Balinski
      F := representative of the orbit of minimum incidence size      # deterministic choice
      mark O expanded
      H  := setwise stabilizer generators of F inside G               # Schreier, exact words
      for R in EnumerateFacetOrbits(P[F], H):                        # RECURSION: ridges of F
          F' := Lift(P, F, R)                                        # 2-dim annihilator, exact
          insert orbit-canonical(F') if new
  return registry
```

`Lift(P, F, R)`:

1. `A := nullspace(transpose(rows[R]))`, exact rational; assert `dim A = 2`.
2. For every row `x`, `p(x) := (A[1].x, A[2].x)` in `Q^2`; drop `p(x) = 0`.
3. `K := {y in Q^2 : y.p(x) >= 0 for all x}`; assert its extreme rays are exactly two
   (a genuine ridge supports exactly two facets).
4. Map each extreme ray back: `a := y[1]*A[1] + y[2]*A[2]`. Of the two, exactly one
   annihilates **all** of `F` (that is `F`'s own normal); take the other.
5. `F' := {i : a . rows[i] = 0}`, then verify `a . rows[i] >= 0` for all `i` and
   `rank(rows[F']) = n-1`. Store `a` normalised to a primitive integer vector.

Three details that the port must get right and that the sources pin down:

* **Projection.** At depth `k` the ambient space is `span(rows)`, not `R^n`, so
  `nullspace` is taken after projecting to pivot columns (`__ProjectionFrame`,
  `Fundamental.g:518,607`). Projecting is injective on the span, so incidences
  computed in projected coordinates are incidences in the original rows.
* **Orientation.** Incidence sets are unordered; `A` is only determined up to scale,
  so both extreme rays of `K` must be tested for feasibility, not signed.
* **Facet-stabilizer.** The recursion's group is the **setwise** stabilizer of `F`,
  not the pointwise one. Its generators are obtained from the orbit words by
  Schreier's lemma — no group enumeration is required, and the port already builds
  exactly such words in `_iter_facet_mask_orbits` (`:1518-1618`).

---

## 3. Completeness and termination

### 3.1 Hypotheses

* **H1** `rank(rows) = n` and `C = cone(rows)` is pointed, so `C` is a full-dimensional
  pointed polyhedral cone.
* **H2** `G` is finite, permutes the rows, and maps `C` onto itself; therefore it maps
  facets to facets and commutes with incidences (`I(gF) = g·I(F)`).
* **H3** `OneFacet` returns a genuine facet (an incidence of rank `n-1` with a
  supporting form nonnegative on all rows). In the port that is `_kernel_flipping`,
  which is exact.
* **H4** `Lift` implements: for a ridge `R ⊂ F` with `rank(rows[R]) = n-2`, the two
  facets containing `R` are the two extreme rays of `{y : y . rows[R] = 0, y . rows >= 0}`.
* **H5** The recursion applied to `rows[F]` with the set stabilizer returns exactly
  one representative per orbit of facets of `cone(rows[F])` (induction hypothesis).
* **H6** Balinski: the facet-adjacency graph of an `n`-dimensional pointed cone is the
  vertex graph of an `(n-1)`-polytope, hence `(n-1)`-connected.

### 3.2 The Balinski stopping rule is sound [proof]

Let `E` = expanded facets, `U` = recorded-but-unexpanded, `X` = not yet recorded.
Expanding a facet records all its neighbours, so **there is no edge between `E` and
`X`** — this is an invariant of the loop, not an assumption.

*Lemma (Balinski).* The facet graph of a full-dimensional pointed cone in `R^n` is
`(n-1)`-connected. *Proof.* Choose `l` strictly positive on `C \ {0}` (H1); then
`P := {y in C* : l(y) = 1}` is a polytope of dimension `n-1`, and its vertices are the
extreme rays of `C*`, i.e. the facets of `C`, with adjacency matching facet adjacency
of `C`. Balinski's theorem gives that the vertex graph of a `d`-polytope is
`d`-connected. ∎

*Claim.* If `|U| <= n-2` then `X = ∅`. Suppose `X != ∅`. Since `|U| <= n-2 <= (n-1)-1`,
`facetgraph(C) - U` is still connected; it contains the nonempty sets `E` and `X`, so
some edge of the connected graph joins the `E` side to the `X` side, i.e. there is an
edge `E`-`X` (edges inside `E` or inside `X`, and edges `E`-`U`, were deleted with
`U`). Contradiction with the invariant. ∎

So stopping at `|U| <= n-2` loses nothing. The source's `nbUndone <= TheDim-1` with
`TheDim = rank-1` is exactly this, and it counts elements — which is required: the
bound is about facets, not orbits.

### 3.3 The recursion returns exactly the ridges of the representative [proof]

Given H1, `cone(rows[F])` is a cone of dimension `n-1`, full-dimensional in `span(F)`
and pointed. Its facets are, by definition, its `(n-2)`-faces, and a `(n-2)`-face of
`F` is an `(n-2)`-face of `C` contained in `F`, i.e. a ridge of `C` in `F`. The map is
a bijection between facets of `cone(rows[F])` and ridges of `C` inside `F`. H5 then
gives one representative per **orbit** of ridges under the set stabilizer of `F`.

### 3.4 Every lift produces a genuine new facet [proof]

Let `R` be a ridge in `F`, `rank(rows[R]) = n-2`. Put
`A = {y : y . rows[R] = 0}`, `dim A = 2`. `K := {y in A : y . rows >= 0}` is pointed
(if `y, -y in K` then `y` annihilates every row, so `y = 0` because `rank = n`) and
2-dimensional (it contains `F`'s normal and, since `R` is a ridge, the other facet's
normal). Hence `K` has exactly two extreme rays `a, a'`, and they are the normals of
`F` and of the unique other facet through `R`. `a'` is therefore a nonzero functional
nonnegative on all rows. Its zero set `Z` contains `R` and cannot equal `R` up to
span: if `a'` vanished on all of `F` then `a'` would lie in `(span F)^perp = span(a)`,
contradicting that `a, a'` are two distinct extreme rays. So some row `r in F` has
`a'.r > 0`, and `F' := {x in C : a'.x = 0}` is a proper face containing the ridge `R`;
being cut out by a single nonzero form it has dimension exactly `n-1`, i.e. it is a
facet. ∎

`R` lies in exactly the two facets `F`, `F'` — this is what "ridge" means, and it is
why one lift per ridge is the complete list of neighbours of `F`.

### 3.5 Completeness [proof]

Let `S` be the returned orbit set. Every element is a facet by H3 and §3.4 and the
invariant `is_facet` check. Conversely, suppose `F` is a facet not in `S`. The facet
graph is connected (§3.2 lemma), so take a path from an expanded facet to `F`; let
`F_0` be the last vertex of `S` on it. `F_0` was expanded (the loop only stops when
every recorded orbit is expanded or §3.2 applies), and expansion records all
neighbours of `F_0` (§3.3 + §3.4), so the next vertex of the path is in `S` —
contradiction, unless the loop stopped by §3.2, which is exactly the case excluded by
§3.2's claim. Hence `S` contains every facet. Each orbit is inserted once because
insertion keys on an orbit-invariant (§2, registry), and `G` maps facets to facets
(H2), so distinct orbits stay distinct. ∎

### 3.6 Termination [proof]

The facet set of `C` is finite. Each iteration either expands an orbit (strictly
decreasing the finite set of unexpanded orbits, never inserting a duplicate) or
stops. The recursion decreases `rank` by exactly 1 each level (`n -> n-1`), so it
terminates at the explicit `n <= 2` base case after at most `n-2` levels. The
registry and the memo are finite. ∎

---

## 4. Exactness

* Every object is an integer or rational vector; no floating point is used anywhere:
  `nullspace`, the 2-dimensional extreme-ray test and the incidence filter are
  `Q`-exact comparisons against 0 and against each other.
* A facet normal is emitted as the **primitive integer** vector of the rational
  nullspace (the port already has `_source_remove_fraction` and
  `_pinned_nullspace_matrix`), never as a floating approximation, and never with a
  tolerance on `a . x >= 0`.
* The lift's correctness check is the source's own: the ridge annihilator must be
  exactly 2-dimensional, and exactly one of the two extreme rays must fail to
  annihilate `F`. Both are equalities, not thresholds; a failure is an exception, not
  a silently dropped facet.
* **The witness is exact.** For each orbit member the enumerator keeps the Schreier
  word as a nested pair (generator index, shorter word) — the representation already
  used by `_iter_facet_mask_orbits` (`:1518`) and consumed at `:2963-2999` — so the
  isometry carrying the representative to any member is a finite composition of
  exact lattice isometries. Nothing is "up to a bound": the word is a certificate, and
  applying it to the representative's normal is verifiable by exact transport
  (`M^{-1} a` for an integral unimodular `M`, the mechanism at `:2985-2999`).
* **The adjacency witness is exact**: the ridge orbit together with the lift's
  identification of which extreme ray is `F`'s own normal. There is no LP, no
  numerical rank decision, and no heuristic choice.

---

## 5. Cost model, and why the 49861-facet materialization never happens

Let `F` = number of facet orbits of `C`, `R_O` = number of ridge orbits of the
representative of orbit `O`, `I(F)` = the incidence of `F`.

**Per level.** One subcone enumeration per expanded orbit, and one lift per ridge
orbit of its representative. A lift costs one rational nullspace of an
`(n-2) x n` integer matrix plus `m` exact dot products and a `Filtered` — i.e.
`O(n^3 + m n)` per lift. The subcone enumeration on `I(F)` costs a rank-`(n-1)`
instance with `|I(F)|` rows.

**Total (level 0)** = `#expanded orbits` subcone calls + `sum_O R_O` lifts, each lift
`O(m n)`. Crucially, `m` is the **cell's** ray count, but the *subcone* work at level
`>=1` is on `|I(F)|` rows only.

**Nothing materializes the facet list.** Every facet enters the state as one integer
(its incidence key) inside one orbit buffer; the loop never holds the H-presentation,
never asks an engine for support hyperplanes, and the only object proportional to the
A7 ray count is the ray matrix itself. `_expand_traversal_once` already consumes one
orbit at a time and drops it.

**Measurements [measured]** — `orbit_facets_proto.py`, Sage 10.8, exact `Fraction`
arithmetic, complete facet list from `Polyhedron.facet_normals` as oracle. Ten
pointed symmetric cones, `dim` 4..6, 5..18 rays, `|G|` up to 720; all ten returned
**exactly** the oracle's orbit partition, and identically under three different seeds:

| specimen | rays | dim | \|G\| | facets | facet orbits | max \|I(F)\| | level-0 exp/lifts | total subcone calls | total lifts | recipe |
|---|---|---|---|---|---|---|---|---|---|---|
| simplex4 | 5 | 5 | 120 | 5 | 1 | 4 | 1/1 | 4 | 3 | full |
| simplex5 | 6 | 6 | 720 | 6 | 1 | 5 | 1/1 | 2 | 2 | full |
| S4 (1,1,0,0)+(1,1,1,0) | 10 | 4 | 24 | 8 | 2 | 4 | 2/2 | 3 | 3 | full |
| S4 (1,2,0,0)+(1,1,1,1) | 13 | 4 | 24 | 8 | 2 | 6 | 2/3 | 4 | 5 | full |
| S4 (2,1,1,0)+(0,1,1,1) | 16 | 4 | 24 | 14 | 3 | 4 | 3/4 | 8 | 8 | full |
| S4 (1,0,1,0)+(1,2,2,0) | 18 | 4 | 24 | 8 | 2 | 6 | 2/2 | 3 | 3 | full |
| S3xS2 | 9 | 5 | 12 | 6 | 3 | 6 | 3/8 | 14 | 32 | full |
| cyclic6 | 12 | 6 | 6 | 30 | 5 | 7 | 5/29 | 127 | 557 | full |
| cyclic6b | 12 | 6 | 6 | 20 | 4 | 7 | 4/23 | 105 | 501 | full |
| cyclic6c | 12 | 6 | 6 | 42 | 7 | 6 | 7/39 | 119 | 687 | full |
| lorentzian_S3 | 10 | 4 | 6 | 8 | 4 | 6 | 4/6 | 9 | 13 | full |

Three facts the numbers establish, in the order they matter for the port:

1. **Level-0 cost is exactly what the model says**: `#facet-orbit` expansions and
   `sum_O R_O` lifts. For `cyclic6`, 5 expansions and 29 lifts for 5 orbits (= 5.8
   ridge orbits per facet orbit).
2. **The recursion, not the top level, is the cost centre, and it multiplies.**
   `cyclic6` totals 127 subcone calls and 557 lifts for 5 top-level orbits — a
   recursion tree of depth 4. Memoizing subconfigurations (the source's
   `BankRecording`) took it from 515/1229 to 127/557; the Balinski stop took it
   further to **52/201**. Both are needed; neither is optional.
3. **`|I(F)|` stays small while the facet count does not**: `cyclic6c` has 42 facets
   but `max |I(F)| = 6`, so the level-1 instances are `6 x (n-1)` problems. This is
   the mechanism that makes the A7 route viable where one `SupportHyperplanes` call is
   not: cost tracks `sum_O R_O` and the per-level row counts, not the facet count.

**Extrapolation to the A7 cell [unverified].** With the task's `m = 49861`, `n = 10`,
level 0 costs `#facet orbits` subcone calls plus `sum_O R_O` lifts at
`O(m n) = 5e5` exact dot products each; level `k` instances have row counts bounded
by the level-`(k-1)` incidences. The number this hinges on — the A7 facet-orbit count
and its ridge-orbit count — is exactly what I could not measure, and it must be
instrumented on the first real run (§7) rather than assumed.

---

*(continued in §6: port seam, and §7: validation)*

---

## 6. Implementation plan for the port

All line numbers are `src/sage_indefinite_port/indefinite/lorentzian_cells.py` as it
stands now.

### 6.1 The seam, in one sentence

`_facet_mask_data` currently does two jobs — build the ray permutations from the
stabilizer, and *obtain* the facet keys from Normaliz. Split the first out, and let
the second be supplied either by Normaliz (unchanged, small cells) or by the new
orbit-wise enumerator (large cells). `_iter_facet_mask_orbits` stays the only owner of
orbit BFS and Schreier words, so the public views need no change.

### 6.2 New functions (additions only)

```python
_ORBIT_FACET_ROUTE_THRESHOLD = 10_000     # the switch cell_cone:1012 already uses

def _facet_ray_permutations(self, cell: LorentzianPerfectCell) -> tuple[array, ...]:
    """Packed ray-position permutations of the cell-stabilizer generators.

    Pure extraction of lorentzian_cells.py:1344-1408 (the loop that turns
    ``cell_stabilizer`` into packed position arrays).  Convention: a generator's
    permutation p satisfies ``p[i]`` = the position of the image of ray i.
    """

def _use_orbit_facet_route(self, cell: LorentzianPerfectCell) -> bool:
    return len(cell.rows) >= _ORBIT_FACET_ROUTE_THRESHOLD

def _facet_orbit_mask_stream(
    self,
    cell: LorentzianPerfectCell,
    permutations: tuple[array, ...],
    *,
    retain_words: bool = True,
) -> Iterator[list[tuple[FacetMaskKey, SchreierWord | None] | None]]:
    """The single dispatch point: Normaliz-derived keys, or orbit-wise enumeration."""

def _facet_orbit_representative_masks(
    self,
    cell: LorentzianPerfectCell,
    permutations: tuple[array, ...],
    *,
    normals: dict[FacetMaskKey, tuple[int, ...]] | None = None,
) -> Iterator[FacetMaskKey]:
    """One facet key per stabilizer orbit, enumerated without any engine call.

    ``normals``, when given, is filled with the primitive integral normal of each
    yielded representative (the hook that later replaces ``cell_cone`` in
    ``perfect_domain_traversal``; see 6.5)."""

def _subconfiguration_facet_orbits(
    self,
    rows: tuple[tuple[int, ...], ...],
    generators: tuple[array, ...],
    *,
    depth: int,
    cache: dict[tuple, tuple[tuple[int, ...], ...]],
) -> tuple[tuple[int, ...], ...]:
    """One facet incidence per orbit of cone(rows) under <generators>, as positions."""

def _mask_orbit_stabilizer_generators(
    self,
    orbit_masks: tuple[FacetMaskKey, ...],
    orbit_permutations: tuple[array, ...],
    permutations: tuple[array, ...],
) -> tuple[array, ...]:
    """Schreier generators sigma_i * g * sigma_l^{-1} of the setwise stabilizer."""

def _ridge_adjacent_positions(
    self,
    rows: tuple[tuple[int, ...], ...],
    facet: tuple[int, ...],
    ridge: tuple[int, ...],
) -> tuple[int, ...]:
    """The lift of Fundamental.g:614-688 (see 2.3): the unique facet != facet on `ridge`."""

def _primitive_normal(self, rows, incidence) -> tuple[int, ...]:
    """Primitive integer supporting form of `incidence`; raise unless it is one."""
```

### 6.3 The four edited places

1. **`_facet_mask_data` (:1330)** — pure refactor: its permutation-building block
   moves into `_facet_ray_permutations`, which it then calls. No behaviour change; the
   Normaliz key stream and `canonical_keys` stay exactly as they are.
2. **`_iter_facet_mask_orbits` (:1518)** — add `verify_presentation: bool = True`
   around the final `visited_count != facet_count` consistency check (that check
   assumes the input stream is the *complete* H-presentation; on the orbit route the
   stream is one seed per orbit and the check must compare orbit count instead). No
   other change: it already yields one list per orbit with exact words, and with
   `canonical_keys=None` it BFSes a generator orbit directly, which is exactly what an
   orbit-of-representatives stream needs.
3. **`_facet_orbit_transversal_masks` (:1621)** — replace its body with
   ```python
   permutations = self._facet_ray_permutations(cell)
   return self._facet_orbit_mask_stream(cell, permutations), permutations
   ```
   `facet_orbits` (:1669), `_facet_orbit_transversal_words` (:1648) and
   `facet_orbit_permutations` (:1279) then inherit the route unchanged.
4. **`_expand_traversal_once` (:2564)** — replace
   ```python
   facet_masks, generator_permutations, canonical_facet_keys = backend._facet_mask_data(source, use_owned_cone=False)
   for orbit in backend._iter_facet_mask_orbits(facet_masks, generator_permutations, retain_words=..., canonical_keys=canonical_facet_keys):
   ```
   with
   ```python
   generator_permutations = backend._facet_ray_permutations(source)
   for orbit in backend._facet_orbit_mask_stream(source, generator_permutations, retain_words=self._retain_facet_orbits):
   ```
   Everything after the loop (flip, representative bucketing, one adjacency per orbit,
   `_traversal_facet_orbits`/`_traversal_facet_generators` bookkeeping) is untouched.

`cell_cone` (:1011) is **not** touched. That is deliberate: it keeps
`tests/test_large_cell_cone_route.py` a genuine test of the Normaliz route, and it
keeps the small-cell behaviour bit-identical.

### 6.4 The enumerator, in port terms

* Projection: `rank`-many pivot columns from the same elimination
  `_pinned_independent_column_indices` (:5603) / `_pinned_nullspace_matrix` (:5586)
  already provide. At depth `k` the rows are the parent's rows restricted to the
  incidence, re-projected.
* Seed: `_kernel_flipping` (:5795) on the configuration gives one critical set; its
  incidence is the seed facet. (The source uses an LP here,
  `LinearProgramming.g:20-31`; the port's flip search is exact and needs no engine,
  and by §3.5 the choice of seed does not affect the result — I verified
  seed-independence on all ten specimens.)
* Registry: `dict[FacetMaskKey, tuple[FacetMaskKey, SchreierWord | None]]` plus a
  `done` set; insert via the orbit-BFS that `_iter_facet_mask_orbits` already
  implements, one seed at a time, so canonicalisation and word bookkeeping are shared
  code rather than a second implementation.
* Stabilizer: `_mask_orbit_stabilizer_generators` as above. Cost is
  `O(|orbit| * #generators)` array compositions — the same order as the BFS that
  built the orbit, so it does not change the asymptotic cost.
* Memo: key on the **exact** subconfiguration `(rows, generators)` first. That is
  sound (the answer is a function of that pair alone) and cheap. The reference's
  stronger bank keys on an *isomorphism class*
  (`BankRecording`, `DualDescriptionAndADM.g:12-74`); the port can reach that later
  with `_integral_configuration_invariants` (:1753) + `cell_transporter` (:1856), and
  it must keep the exact isomorphy test, not just the invariant, or the memo is wrong
  in exactly the way `FuncRetrieveObject` guards against.
* Balinski stop: keep it (proved, §3.2), counting **elements**:
  `sum(len(G) // |Stab(F)| for F in unexpanded) <= n - 2`. Do **not** port the
  `ComputeIntersectionUndone` disjunct (§1.4) and do **not** port `IsRespawn` (§1.5).

### 6.5 The remaining full-materialization, and how to remove it

`perfect_domain_traversal` (:2790) calls `backend.cell_cone(cell)` for **every**
quotient cell (:2831) and consumes `cell_cone.halfspace_covectors()`. On the A7 cell
that is the 1800 s call. The orbit route does not remove it: it only removes the
*duplicate* facet enumeration inside the traversal.

The removal is mechanical once representatives carry normals
(`_facet_orbit_representative_masks(..., normals=...)`): a representative's covector is
the enumerator's exact normal, and every other member's covector is `M^{-1} a` for the
unimodular `M` of its Schreier word — the transform already implemented and argued at
:2985-2999. Two things must be re-expressed at the same time:

* `facet_keys_by_cell` / `by_incidence` (built from `cell_cone`) become the orbit keys
  from `_traversal_facet_orbits`; the completeness assertion at :3178
  (`paired_facets[id(cell)] != expected_paired_facets`) is then a comparison of orbit
  keys, which is *equivalent* because the orbit set is complete (§3.5).
* `cell_cones[...].halfspace_covectors()` is only needed for representatives; the
  side-pairing witnesses and the containment checks
  (`_iter_cell_halfspace_pairing_rows`, `common_face_cache`) accept covectors, not
  cones.

I did not write that change here, and it is [unverified] in the strongest sense: it is
the piece that decides whether the A7 cell finishes, and it needs a run.

### 6.6 Small cells keep the current behaviour

`_use_orbit_facet_route` is `len(cell.rows) >= 10_000`, the same threshold and the same
sense as `cell_cone:1012`, so:

* every cell the current suite exercises (`_isotropic_cell` with 6 rows, A2 with 27,
  A3+<1>+<-1> with 89) takes the existing Normaliz-key path, byte for byte;
* `cell_cone` is unchanged, so `tests/test_large_cell_cone_route.py` continues to
  compare PPL against Normaliz exactly as today;
* only the 49861-row configuration takes the new path, and the two paths are
  comparable by construction (test T2 below), which is the whole point.

---

## 7. Validation

### 7.1 What I ran here [measured]

`orbit_facets_proto.py` (in this directory) is a standalone implementation of §2 in
plain Python + Sage, no port imports. For each specimen it

1. computes the **complete** facet list from the independent engine
   (`Polyhedron(rays=...)`, `Hrepresentation()` normals) and its orbit partition,
2. runs the algorithm of §2, and asserts the two orbit partitions are **equal as sets**,
3. asserts every enumerated incidence *is* a facet (exact normal, nonnegative on all
   rows, `rank = n-1`, zero set exactly the incidence),
4. repeats with two further random seeds and asserts the same result (seed
   independence),
5. reports the cost counters of §5.

Ten specimens, all passing, results in `results_full.txt` (full saturation),
`results_balinski_only.txt` (Balinski stop only), `results_cutoff.txt` (both upstream
stopping rules, the second emulated for comparison). What the specimens cover:
non-simplicial facets (`max |I(F)|` up to 7 of 12 rays), facet stabilizers of order up
to 120, group orders up to 720, recursion depth up to 4, and a Lorentzian-flavoured
cone in signature (3,1). What they do **not** cover: 49861 rays, and the port's own
data types — that is the gap §0 describes.

### 7.2 What must be run before this is accepted

1. The two existing files, unchanged, on a Sage/py version where the provider imports:
   `PYTHONPATH=src:$HOME/sage-extensions-42d7361c/src .venv/bin/python -m pytest
   tests/test_facet_orbit_contracts.py tests/test_large_cell_cone_route.py -q`
   (currently blocked, §0).
2. New test `tests/test_orbit_facet_route.py`, on the two bounded specimens the task
   names — `A2 + <1> + <-1>` (rank 4, 27 rays, 6 facets, 4 ray orbits) and
   `A3 + <1> + <-1>` (rank 5, 89 rays, 24 facets, 6 ray orbits) — asserting:
   * **T1** the orbit partition of `cell_cone(cell).halfspace_covectors()` incidences
     equals the partition produced by `facet_orbits(cell)`; and for each orbit, the
     representative's enumerator normal is a positive rational multiple of the full
     route's covector for the same incidence, compared exactly;
   * **T2** route agreement: forcing `_use_orbit_facet_route` to `False` and to `True`
     gives the same set of canonically-minimal orbit keys;
   * **T3** the rank-2 base case is exercised (any cell with a rank-3 facet cone, and
     `rank 2` directly), because that is the only branch with no induction;
   * **T4** cost identities, not time thresholds: `#subcone calls == #expanded orbits`
     and `#lifts == sum over expanded orbits of #ridge orbits`, and the total number of
     facets covered equals the full route's facet count.
3. The A7 cell, **once**, instrumented (not on this swarm's budget): assert completion,
   that every representative normal is exact and nonnegative on all 49861 rows, the
   orbit-stabilizer identity `sum(orbit sizes) == #facets enumerated`, and record
   `#facet orbits`, `sum_O R_O`, `max |I(F)|`, recursion depth and peak memory against
   the 6144 MiB guard. Those five numbers decide whether the plan in §5 holds at full
   size; they are the only things I could not establish.

---

## 8. Open items and risks, stated plainly

1. **`|I(F)|` is not bounded a priori.** Level 0 cost is fine, but a facet incident to
   a large fraction of the 49861 rays makes the level-1 instance large. Nothing in the
   source bounds it either; the source's answer is `IsRespawn` (abandon the ADM), which
   §1.5 rejects. Instrument it first (§7.2.3) before promising a budget.
2. **The recursion tree multiplies.** Measured 5x-25x the top-level work here; the memo
   and the Balinski stop are what make it affordable. Both must be in from the start.
3. **The upstream `|intersection of undone facets| > 0` cutoff is unproved** and fires
   often in my runs. It is excluded from §6 deliberately; if someone wants it back it
   needs a proof, not a test.
4. **`perfect_domain_traversal` still calls `cell_cone`** (§6.5). Until that call site
   changes, the A7 specimen still hits the 1800 s wall regardless of the facet
   enumeration. This is the single largest remaining piece of work and it is not
   in this document as a finished design.
5. **The port's mask convention stores the smaller of incidence and complement**
   (:1420-1445). The enumerator must carry *incidences* internally (the lift needs
   them) and may only apply that compression at the boundary where keys are yielded.
6. **`flip_across` recomputes the facet normal** by nullspace + `_kernel_flipping`
   (:1708-1751). On the orbit route the traversal already holds the exact normal of the
   orbit representative, so an optional `facet_normal=` argument avoids recomputing it
   per orbit; that is a pure speed change and should be done after T1-T4 pass.
