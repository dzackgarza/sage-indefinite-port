# Orbit-representative facet enumeration for a Lorentzian perfect cell cone

Swarm `smvzd6hvn`, member m4, direction **d8**.
Target: replace `LorentzianPerfectLocalBackend.cell_cone`'s complete H-description of the
49 861-ray A7 cell by an exact enumeration of **one facet per cell-stabilizer orbit**.

Files of this directory:

| file | content |
| --- | --- |
| `proof.md` | this document: construction, hypotheses, proofs, exactness, cost model, implementation plan, validation plan |
| `notes.md` | what was measured, what could not be measured, and why |
| `proto/cone_orbit_engine.py` | exact standalone implementation of the algorithm (Sage 10.8 + pplpy) |
| `proto/specimens.py`, `proto/driver.py` | five bounded specimens and the brute-force comparison |
| `proto/lp_timing.py` | cost of the exact LP primitive at m = 49 861 rows, n = 10 |

---

## 1. The construction

### 1.1 Data and conventions

* $V=\mathbb{Q}^n$, $n=10$; the cell cone is $C=\operatorname{cone}(r_1,\dots,r_m)\subseteq V$ with
  $r_i\in\mathbb{Z}^n$ the cell's coordinate rows, $m=49\,861$ for the A7 specimen,
  $\max_i\lVert r_i\rVert_\infty=15$.
* Pairing: $\langle c,x\rangle=\sum_i c_ix_i$ (plain integer dot product on coordinate tuples;
  the Lorentzian metric plays no role in the H-description, only in the geometry that produces the rows).
  The dual cone is $C^{*}=\{c\in V:\langle c,r_i\rangle\ge 0\ \forall i\}$.
* The stabilizer $G\le \mathrm{GL}_n(\mathbb{Z})$ is given by generators $A_1,\dots,A_k$ acting on
  rows by $x\mapsto xA$ (the port's row convention: `g(e_i)` is row $i$ of `action`). Each generator
  permutes the ray multiset: there is a permutation $\pi_A$ of $[m]$ with $r_{\pi_A(i)}=r_iA$.
  `_facet_mask_data` already packs exactly these $\pi_A$ (`permutations`, unsigned words).
* The contragredient action on facet normals is $c\mapsto c\,(A^{-1})^{T}$; since $A$ is unimodular it is
  again integral. On incidences it is the same $\pi_A$: with $Z_c=\{i:\langle c,r_i\rangle=0\}$,
  $$Z_{c(A^{-1})^T}=\{\,i:\langle c,r_iA\rangle=0\,\}=\pi_A^{-1}(Z_c).$$
  The port's mask action (`_permuted_facet_key`) is this permutation; nothing new is needed on that side.
* Everything below is stated for a *cell*; the same code runs for every quotient cell.

### 1.2 Hypotheses (all used, each stated where it enters)

* **(H1) $C$ is full-dimensional**: $\operatorname{rank}\{r_i\}=n$. Measured for the A7 specimen
  (row rank 10) and for every `initial_cell` output (the rows span the ambient lattice — the port asserts
  this in `facet_orbits`: "a perfect-cell configuration does not span its ambient lattice").
  Used for: facet $\leftrightarrow$ extreme ray of $C^{*}$; $C^{*}$ pointed; $P_C$ bounded.
* **(H2) $u=\sum_i r_i$ is an interior point.** This follows from (H1) and is not an extra assumption;
  the remark is recorded because $P_C$ is defined by it (§1.3). **No pointedness hypothesis is used
  anywhere**: the covering lemma (Lemma 3) and the forced-progress lemma (Lemma 4) need only that
  $P_C$ is a non-empty polytope, i.e. bounded and containing no line, which holds because $P_C$ lies in
  the affine hyperplane $\langle u,\cdot\rangle=1$.
* **(H3) $G$ permutes the row multiset and stabilizes $C$.** Checked exactly at load time by
  `_facet_mask_data` (it raises "a cell-stabilizer generator does not preserve its rays").
* **(H4) The configuration is duplicate-free.** Asserted by the port ("a perfect-cell configuration
  contains repeated rays").
* **(H5) Exact LP is available** over $\mathbb{Q}$ with integer input; no tolerance appears in any
  decision (see §3).

### 1.3 The polar polytope

Let $u=\sum_{i=1}^{m}r_i$. Under (H1), $\langle c,u\rangle>0$ for every $0\ne c\in C^{*}$
(if $\langle c,u\rangle=0$ then all $\langle c,r_i\rangle=0$, so $c\perp\operatorname{span}\{r_i\}=V$),
hence $u\in\operatorname{int}C$ and

$$P_C:=\{y\in V:\ \langle r_i,y\rangle\ge 0\ \forall i,\ \langle u,y\rangle=1\}\subseteq C^{*}$$

is a non-empty polytope: it is a section of the pointed cone $C^{*}$ by an interior functional.

> **Lemma 1 (facets are vertices).** Assume (H1). The map $f\mapsto F_f:=C\cap\{f=0\}$ is a bijection
> from the vertex set of $P_C$ onto the set of facets of $C$.
>
> *Proof.* $C^{*}=\operatorname{cone}$ of its own extreme rays; a vertex of $P_C$ is the normalization
> $\langle u,\cdot\rangle^{-1}$ of an extreme ray of $C^{*}$, and for $0\ne f\in C^{*}$ the zero set
> $F_f$ is a facet iff $f$ is an extreme ray of $C^{*}$ (standard; $C$ full-dimensional closed convex).
> Distinct extreme rays give distinct facets, distinct facets give distinct normalized extreme rays.
> $\square$

> **Lemma 2 (a facet is its incidence).** Assume (H1),(H4). For a facet $F$, $F=\operatorname{cone}\{r_i:i\in Z\}$
> with $Z=\{i:\langle f,r_i\rangle=0\}$ for its normal $f$, the rows indexed by $Z$ span an $(n-1)$-dimensional
> subspace, and the primitive integer generator $f_Z$ of $\operatorname{span}\{r_i\}_{i\in Z}^{\perp}$,
> sign-normalized by $\langle f_Z,r_i\rangle\ge0$ for all $i$, is the unique primitive integer normal of $F$.
> Conversely, if $Z\subseteq[m]$ with $\operatorname{rank}\{r_i\}_{i\in Z}=n-1$ admits an annihilator that is
> non-negative on every row, then $Z$ is the incidence of exactly one facet.
>
> *Proof.* A face of $\operatorname{cone}\{r_i\}$ is spanned by the rays it contains, so $F=\operatorname{cone}\{r_i\}_{i\in Z}$;
> $\dim F=n-1$. The annihilator of $Z$ is one-dimensional over $\mathbb{Q}$ and contains an integer vector
> (integer data), whose primitive generator is unique up to sign; the sign is fixed by non-negativity.
> Uniqueness of the incidence: if two facets had the same $Z$ they would be the same spanned face.
> $\square$

So **the incidence is a complete facet key**, and the *canonical normal* is computed from the incidence by
integer linear algebra — never read off the LP's floating point output. This is the exactness mechanism of §3.

### 1.4 The algorithm

Two primitives.

*`vertex(S, ψ)`* — for a tuple $S$ of ray positions (a chain) and an integer objective $\psi$:

$$\Psi_S:=\{y\in P_C:\langle s,y\rangle=0\ \text{for every } s\in S\},\qquad
  \text{return a vertex of }\Psi_S\text{ maximizing }\psi .$$

$\Psi_S$ is a face of $P_C$, hence a polytope; a vertex is obtained by one exact LP followed, if the
returned optimal point is not a vertex, by the lexicographic cascade
$\psi,e_1,\dots,e_n$ (pin the achieved value each time; after $n+1$ stages the feasible set is a single
point, which is a vertex). The exact vertex test is: the rows vanishing at the point have rank $n-1$.

*`normal(Z)`* — the canonical primitive integer normal of Lemma 2, or `None` if $Z$ is not a facet.

The driver, `enumerate(S)` for a chain $S$ with $d=\operatorname{rank}(S)$:

1. **Node search (forced progress).** Keep $B=\emptyset$ (normals found at this node). While
   $\dim\bigl(\operatorname{span}S+\sum_{f\in B}\operatorname{span}(G f)\bigr)<n$: take
   $0\ne\psi$ in its orthogonal complement, try the two signs, and run `vertex(S, ±ψ)`; if the resulting
   incidence is new, register it and add its normal to $B$, else stop.
2. **Base case $d=n-1$.** $\Psi_S$ is a single vertex; the facet is $\operatorname{cone}\{r_i: r_i\in\operatorname{span}S\}$,
   pure exact linear algebra, no LP.
3. **Base case $d=n-2$.** Every facet containing $S$ contains exactly one more ray: for a transversal of the
   $H$-orbits ($H$ = setwise stabilizer of $\operatorname{span}S$) of the rows outside $\operatorname{span}S$,
   register the incidence $\{i: r_i\in\operatorname{span}(S\cup\{t\})\}$; `normal` silently rejects the
   candidates whose hyperplane is not supporting.
4. **Recursion, $d\le n-3$.** For a transversal $\{t\}$ of the $H$-orbits ($H$ = setwise stabilizer of
   $\operatorname{span}S$) of the rows outside $\operatorname{span}S$, recurse into $S\cup\{t\}$.

"Register" means: if the incidence is already in the visited set, do nothing; otherwise add its whole
$G$-orbit to the visited set and output it as the representative of a new orbit. The visited set is closed
under the generator permutations, so at termination the output holds **exactly one incidence per orbit**.

Cost-relevant detail: `vertex(S,ψ)` solves an LP whose constraint matrix is the full row list. The
specification's implementation strategy (§5.3) is *prune-and-verify*: a float LP on the same rows chooses an
active subset, the exact LP is solved on that subset, and the candidate is verified against all rows exactly;
a mismatch triggers a fallback exact LP on the full list. No decision ever uses a tolerance.

---

## 2. Completeness and termination

> **Lemma 3 (covering).** Assume (H1). Let $S$ be a chain of ray positions with
> $d=\operatorname{rank}(S)\le n-2$ contained in at least one facet, and let $H\le G$ preserve
> $\operatorname{span}S$ setwise. Then every facet containing all rays of $S$ lies in the $G$-orbit of a
> facet containing $\operatorname{span}S$ and some ray $r_t$ with $t\notin\operatorname{span}S$, and it
> suffices to let $t$ range over a transversal of the $H$-orbits of the rows outside $\operatorname{span}S$.
>
> *Proof.* Let $F$ be a facet containing $S$. $\dim F=n-1>d=\dim\operatorname{span}S$, and
> $F=\operatorname{cone}$ of the rows it contains, so $F$ contains a row $r_t\notin\operatorname{span}S$.
> Pick $h\in H$ with $h(t)$ the chosen representative; $h$ preserves $\operatorname{span}S$, so
> $h(F)$ contains $\operatorname{span}S$ and $h(t)$, and $h(F)$ is a facet because $h$ stabilizes $C$. $\square$

> **Lemma 4 (forced progress).** Assume (H1),(H2). Let $B$ be a finite set of facet normals and
> $S$ a chain of ray positions. If $0\ne\psi$ annihilates $\operatorname{span}S$ and
> $\bigcup_{f\in B}Gf$, then for one of the two signs $\sigma=\pm1$, every vertex of $\Psi_S$ maximizing
> $\sigma\psi$ lies outside $\operatorname{cone}(B)$, and hence lies in a facet orbit that has not been
> found yet.
>
> *Proof.* $\max_{\operatorname{cone}(B)\cap\Psi_S}\psi\le 0$, and $\langle\psi,y\rangle=0$ for $y\in\Psi_S$
> is a proper condition unless $\psi$ annihilates all of $\Psi_S$, i.e. all of $C^{*}\cap\operatorname{span}S^{\perp}$
> (impossible for $\psi\ne0$: the cone $C^{*}\cap\bigcap_{s\in S}\ker s$ is full-dimensional in
> $\bigcap_{s\in S}\ker s$, since $S$ is a face chain of the pointed cone $C$). So suppose
> $\psi\le 0$ on that cone; then $-\psi$ lies in $C$ by the bipolar theorem, and as $-\psi$ also
> annihilates $\operatorname{span}S$ we get both $\psi$ and $-\psi$ in $C$ unless $\psi=0$, contradicting
> (H2). Hence one sign has a strictly positive maximum $v>0$ on $\Psi_S$; every vertex attaining it
> satisfies $\sigma\psi>v\cdot$… precisely: $\langle\sigma\psi,y\rangle=v>\max_{\operatorname{cone}(B)}\sigma\psi$,
> so $y\notin\operatorname{cone}(B)$; coinciding orbits would give $y\in\operatorname{span}(Gf)$ for some
> $f\in B$, i.e. $\langle\psi,y\rangle=0$, contradiction. $\square$

> **Theorem 5 (completeness and termination).** Assume (H1),(H2),(H3),(H4). The recursion of §1.4,
> run from $S=\emptyset$, registers every facet of $C$ in exactly one orbit and halts. The output has
> exactly one representative per $G$-orbit of facets.
>
> *Proof.* Termination: the recursion depth is at most $n-2$ and every node has finitely many children;
> within a node the loop of step 1 runs at most $n$ times because each round strictly enlarges
> $\operatorname{span}S+\sum_{f\in B}\operatorname{span}(Gf)$ while step 2/3 are one-shot. (Lemma 4
> guarantees that a round which does not stop the loop has produced a genuinely new orbit.)
> Completeness: prove by downward induction on $d=n-1-\operatorname{rank}S$ that for every chain $S$
> contained in a facet, every facet containing $S$ lies in a registered orbit. For $d\le1$ this is
> step 2/3 by construction (Lemma 2 says every facet containing $S$ is $\operatorname{cone}$ of the rows
> in some hyperplane $\operatorname{span}(S\cup\{t\})$, and the transversal covers all $t$ up to $H$,
> which by Lemma 3 is enough). For $d\ge2$: let $F\supseteq S$; by Lemma 3 there is $t\notin\operatorname{span}S$
> with $F$ in the $G$-orbit of a facet $F'$ containing $\operatorname{span}S$ and $r_t$; the child
> $S\cup\{t\}$ is visited (transversal), and $F'$ contains $S\cup\{t\}$, so the induction hypothesis at
> depth $d-1$ applies. At $S=\emptyset$ every facet contains $S$, so all facets are registered.
> Exactly one per orbit: the visited set is closed under the generators and a seed is only registered when
> unseen, so representatives lie in pairwise distinct orbits, and Theorem 5's first part puts every facet in
> some orbit; orbits partition the facet set. $\square$

**What is *not* claimed.** Completeness does **not** rest on a cheap "cone equality" certificate. The
obvious sufficient test $C^{*}=\operatorname{cone}(f_1,\dots,f_k)$ (equivalently $R_B\subseteq C$ for
$R_B=\{x:\langle f_i,x\rangle\ge0\}$) is sound, exact and cheap, and the prototype implements it — but it is
**too strong**: measured on all five bounded specimens it fails even when the enumeration is complete
(§`notes.md`), because one representative per orbit does not generate the whole dual cone. It is therefore
used as a *validation instrument*, never as a stopping rule. This is the trap to avoid: a stopping rule
based on it either never fires or must be weakened into something unsound.

---

## 3. Exactness: normals, incidences, adjacencies, witnesses

1. **Normals.** Every normal is `normal(Z)`: the primitive generator of $\operatorname{span}\{r_i\}_{i\in Z}^{\perp}$
   computed by an integer kernel over $\mathbb{Z}$ (`_pinned_nullspace_matrix` in the port, or `right_kernel`
   of an integer matrix), scaled by $\gcd$ and sign-fixed by $\langle f,r_i\rangle\ge0$. It is an integer
   vector and it is *derived from the incidence*, not from the LP's floating-point vertex. Two facets are
   equal iff their normals are equal up to positive scaling iff their canonical normals are equal.
2. **Incidences.** $Z=\{\,i:\langle f,r_i\rangle=0\,\}$ by exact integer arithmetic ($m\cdot n$ integer
   multiplications per facet). Two independent checks make a wrong incidence impossible to accept silently:
   (a) `normal(Z)` must reproduce a non-negative primitive annihilator, and (b) $\operatorname{rank}$ of the
   rows in $Z$ must be $n-1$. A candidate that fails either is rejected (and, in the port, counted; a
   non-zero reject count on a specimen whose full cone is known is a hard test failure).
3. **The LP result is a hint, not an oracle.** The candidate vertex from the LP fixes an incidence; the
   *output* is the exact integer object derived from that incidence. If a float LP mis-guesses the active
   set, the exact verification step rejects the candidate (or the exact LP on the active subset disagrees)
   and the fallback full exact LP is run. Hence no produced object depends on a tolerance, a bound, or a
   random seed. (The ψ's in step 1 come from an exact kernel computation, so the search itself is
   deterministic given $u$ and the row order; nothing is sampled.)
4. **Adjacency / transport witnesses.** For each registered orbit the port already has the machinery:
   `_iter_facet_mask_orbits` builds a **Schreier word** $w=(g_{i_1},\dots,g_{i_r})$ with
   $\pi_{w}(\text{rep})\in$ orbit, per facet; the transporting isometry is the product
   $A_{i_1}\cdots A_{i_r}$ of the corresponding exact integral matrices (computed from the word, not from a
   polyhedral computation), and it is *verified* by the port's existing checks
   (`symmetry_action * gram * symmetry_action.transpose() == gram`, $\lvert\det\rvert=1$, integrality of the
   action on a rational basis). The adjacency witness carried into the traversal is the pair
   (orbit representative, word), exactly as `_traversal_facet_orbits` already stores it.
5. **Flips.** `flip_across(cell, facet)` consumes an incidence; it recomputes the facet's normal from the
   incidence's span by exact nullspace and asserts a one-dimensional annihilator. It is therefore
   independent of §1.4's LP and cross-checks the incidence: a wrong incidence would break the degenerate
   `NullspaceMat` assertion. Keep that call path unchanged.

---

## 4. Termination and cost model

Write $n$ for the ambient rank, $m$ for the row count, $N_o$ for the number of $G$-orbits of facets,
$N_f$ for the number of facets, $N_{ro}$ for the number of $G$-orbits of rays.

* **LP count.** Each node of the recursion performs at most $n$ exact LPs in its search loop; the base
  levels ($d=n-1$, $d=n-2$) perform none. Hence
  $$\#\text{LP} \le n\cdot\#\{\text{nodes}\},$$
  and nodes are $G_{H}$-classes of chains $r_{t_1},\dots,r_{t_j}$ with $j\le n-2$, i.e. bounded by the number
  of orbits of chains of faces of $C$. Measured node/LP counts on the five bounded specimens:
  3/5, 2/3, 27/52, 72/143, 145/303. The count is *not* proportional to $N_f$: the 24-cell specimen has
  $N_f=24$, $N_o=1$ and needs 143 LPs, while the wedge specimen has $N_f=12=N_o$ and needs 303.
* **Cost of one LP.** Measured on synthetic data with the A7 shape ($m=49\,861$, $n=10$, entries in
  $[-15,15]$): prune-and-verify totals **0.9 s** per LP (0.22 s float HiGHS + 0.68 s exact re-solve on the
  734-row active set and exact verification against all rows); building the full PPL system for the same
  LP costs 3.0 s at $m=20\,000$, extrapolating to ≈7.6 s at $m=49\,861$. So a run of a few hundred LPs is
  60–300 s, inside the 1800 s cap; thousands of LPs are marginal.
* **The unmeasured quantity** is the node count on A7, which depends on the number of orbits of chains of
  faces — the same kind of quantity the reference's adjacency decomposition pays per facet orbit
  (`f_adj`: one dual description of the representative facet's ridge set per orbit). It is bounded by the
  same $N_o$ in the reference and by the chain-orbit count here; both are unknown for A7 and neither is
  bounded by the facet count.
* **No full-facet materialization.** At no point is an engine given the 49 861 rays and asked for an
  H-description; no structure indexed by *all* facets is built by the enumeration itself: the only
  facet-indexed objects are one incidence per orbit and the per-node normal lists. The A7 ray list is
  touched (i) by one integer dot-product pass per discovered facet, and (ii) by the LP's verification pass.
  **One caveat that must be stated, because it is a real cost:** the port's orbit bookkeeping
  (`_iter_facet_mask_orbits`) saturates each registered orbit by a Schreier BFS, which visits every facet
  of every discovered orbit — at termination, all $N_f$ facets. That traversal is what produces the exact
  Schreier words the traversal's side-pairing API needs, it exists in the port today *after* Normaliz, and
  this change removes the Normaliz call but not that traversal. If the consumer only needs one facet per
  orbit (which is all the reference's `f_adj` needs for `LORENTZ_DoFlipping`), the saturation can be made
  lazy: replace the transversal construction by a Schreier–Sims membership test on demand, and store one
  word per *orbit* instead of one per facet. That is a smaller, separable change with its own proof
  obligation and is listed in §5.4.

---

## 5. Implementation plan in the port

All paths below are under `src/sage_indefinite_port/indefinite/lorentzian_cells.py` unless stated.

### 5.1 New code (one new module, no changes to the vendored provider)

`src/sage_indefinite_port/indefinite/orbit_facets.py`:

```python
@dataclass(frozen=True)
class FacetOrbitEnumeration:
    normals: tuple[tuple[int, ...], ...]        # one canonical primitive normal per orbit
    incidences: tuple[FacetIncidence, ...]      # the matching facet keys, same order
    lp_calls: int
    nodes: int
    rejected_candidates: int                    # must be 0 on validated specimens
    certificate: bool                           # cone(normals) == C* : sufficient, may be False

def polar_vertex(rows, interior, objective, equalities) -> tuple[Fraction, ...] | None
def incidence_at(rows, point) -> FacetIncidence
def canonical_facet_normal(rows, incidence) -> tuple[int, ...] | None
def cone_equality_certificate(normals, rows) -> tuple[bool, tuple[int, ...] | None]
def enumerate_facet_orbit_normals(
    rows, generator_actions, *, row_permutations=None, limit=None
) -> FacetOrbitEnumeration
```

`rows` are the cell's integer coordinate rows, `interior` is $u=\sum r_i$ (primitive-normalized),
`generator_actions` are the stabilizer generators as integer matrices (row convention), and
`row_permutations` are the $\pi_A$ already produced by `_facet_mask_data` (recomputed if omitted).
`enumerate_facet_orbit_normals` is §1.4 verbatim: it needs no polyhedral engine, only exact LP, exact
integer kernels, and the generator matrices; the setwise stabilizer of a chain is obtained from a GAP
permutation group on the ray positions
(`libgap.Group([...])`, `libgap.Stabilizer(group, positions, OnSets)` — the API and its pitfalls are recorded
in `.tmp/facet-orbit-substitution-draft.md`, including that a returned stabilizer must be regrouped with
`subgroup(...)` before use).

### 5.2 Backend entry point

```python
# class LorentzianPerfectLocalBackend
@cached_method
def cell_facet_orbit_enumeration(self, cell) -> FacetOrbitEnumeration:
    """One exact facet representative per cell-stabilizer orbit, for large cells."""
```
It resolves `cell_stabilizer(cell)` (existing), converts each generator to its integer action matrix once
(the same loop `_facet_mask_data` runs), and calls the new module. Threshold: `len(cell.rows) >= 10_000`,
the same constant `cell_cone` uses today.

### 5.3 Changes at the four named sites

1. **`_facet_mask_data` (`:1330`)** — signature unchanged. For large cells the key stream `facet_keys()` is
   fed from `cell_facet_orbit_enumeration(cell).incidences` instead of
   `self._iter_cell_halfspace_pairing_rows(...)` (the Normaliz `SupportHyperplanes` call at `:1051`):
   the incidences are converted to the existing `FacetMaskKey` encoding with the existing
   `_facet_complement_key` block logic. `use_owned_cone=True` (small cells) keeps its present path
   verbatim, so `cell_cone`/`RationalPolyhedralCones.from_rays` are untouched for the 27- and 89-ray
   specimens. Because only one key per orbit is produced, the `canonical_keys` bookkeeping is dropped on the
   large path; `_iter_facet_mask_orbits` already deduplicates overlapping seeds through its `visited` set,
   so a *superset* of orbit representatives (or a repeated one) is harmless — the enumeration may therefore
   also emit keys redundantly without breaking anything.
2. **`_iter_facet_mask_orbits` (`:1518`)** — unchanged. Its final assertion
   (`visited_count != facet_count` → "cell-stabilizer permutations generate incidences outside the
   H-presentation") must be *skipped* on the orbit path, because `facet_count` was the H-presentation's
   facet count, which no longer exists; the exactness it used to certify is now carried by
   `rejected_candidates == 0`, by Theorem 5, and by the tests of §6. Keep it for the small-cell path.
3. **`_expand_traversal_once` (`:2564`)** — unchanged in structure: it still iterates orbits, takes
   `orbit[0][0]` as the representative incidence, calls `flip_across`, and stores the orbit under
   `_traversal_facet_orbits`. The only adjustment is that the seed stream comes from (1); `flip_across`
   recomputes the facet normal from the incidence exactly, so nothing else changes.
4. **`cell_cone` (`:1011`)** — *unchanged*, still the small-cell route and still the reference H-description
   for callers that genuinely want it. What changes is that `perfect_domain_traversal` (`:2790`) must stop
   calling `cell_cone` for large cells: its three uses are
   (a) `representative_keys` / `adjacency_facet_keys` — now comes from `_traversal_facet_orbits` directly
   (already present);
   (b) `boundary_facet_keys` — now "orbits whose `flip_across` returned `None`", recorded by
   `_expand_traversal_once` instead of being recomputed by subtracting key sets;
   (c) `facet_covectors_by_incidence` — now `_integral_covector(cell.lattice, normal)` for each enumerated
   representative, which is exactly the owned covector the code builds today from Normaliz rows.
   `neighbor_factory`'s `target_halfspace_cache` (it transports the *complete* H-family of the target) is the
   one consumer that cannot be served from orbit data; for large cells it should be replaced by the exact
   single-covector certificate "the transporter's image of the source facet covector is a supporting covector
   of the target, and the two facets have the same incidence image", which is $O(mn)$ instead of
   $O(N_f n)$. That change is separable; until it is made, the large-cell path keeps needing the target's
   H-description and the win is confined to the source side. `cell_primitive_rays` (`:1155`) is a different
   query (`NmzResult(normaliz_cone, "ExtremeRays")`) and is untouched; note it is a *second* full-cone
   engine call and is on the isotropic-mark path only.

### 5.4 Optional follow-up (do not bundle)

Lazy orbit saturation (§4): store one Schreier word per orbit and answer "which group element carries facet
$Z$ to the representative" by a Schreier–Sims membership test, so that the traversal no longer visits all
$N_f$ facets. This is the only remaining term that scales with the facet count.

---

## 6. Validation plan

Bounded specimens where **both** routes work (recorded in `.tmp/handoff-a7-facet-route.md`):

| lattice | rank | signature | rays | facets | ray orbits | stabilizer generators |
| --- | --- | --- | --- | --- | --- | --- |
| `A2 + <1> + <-1>` | 4 | (3,1) | 27 | 6 | 4 | 4 |
| `A3 + <1> + <1> + <-1>`-style | 5 | (4,1) | 89 | 24 | 6 | 5 |

New test, `tests/test_facet_orbit_enumeration.py`, asserting on both specimens and both modes:

1. **Orbit completeness.** The `halfspace_covectors()` of `cell_cone` partition into orbits under the
   stabilizer's ray permutations (by exact closure, the way `test_facet_orbit_contracts.py` does it); the
   enumerated representatives number exactly that partition's size, and each enumerated incidence lies in
   exactly one of those orbits. *Every* full-cone facet is in exactly one enumerated orbit.
2. **Exact normals.** Each enumerated normal, as an owned covector, is a positive rational multiple of
   exactly one `halfspace_covectors()` entry, with the exact ratio computed over `QQ` (the pattern of
   `_positive_ratio` in `tests/test_large_cell_cone_route.py`, including its sign-sensitivity assertions).
3. **Exact witnesses.** For each orbit member, the Schreier word's matrix product carries the
   representative's normal to that member's normal (exact integer equality), is an isometry
   ($A G A^T=G$, $\lvert\det A\rvert=1$), and permutes the rays.
4. **No rejected candidates**, and `certificate in {True, False}` with the measured value pinned (currently
   false on both specimens: the document says why).
5. **Kept green, unchanged:** `tests/test_large_cell_cone_route.py` (it calls
   `_normaliz_cell_h_presentation` and `cell_cone` directly on a 6-ray cell, so it exercises the small-cell
   route) and `tests/test_facet_orbit_contracts.py` (covectors annihilate exactly a rank-$(n-1)$ incidence
   and are non-negative elsewhere; generators permute rays and send facets to facets; every facet is in
   exactly one orbit).
6. **Standalone cross-check already run** (`proto/`, this directory): five specimens, exact brute-force
   facet enumeration as ground truth, engine orbit counts equal to the true orbit counts (2, 1, 6, 1, 12 of
   2, 1, 6, 1, 12, counted under the stabilizer subgroup used in each test case).
7. **Then, and only then**, one A7 verification run at the existing cap, no raised limits, no re-runs.

---

## 7. Rejected alternatives (with the reason, so nobody repeats them)

* **Feed one ray per orbit to an engine.** Wrong: cone$\{e_1,e_2\}$ with the swap acting has ray orbits
  $\{e_1\},\{e_2\}$ and $\operatorname{cone}(e_1)\ne\operatorname{cone}(e_1,e_2)$. Any reduction has to be
  inside the enumeration.
* **Stop when a certificate of the form $\operatorname{cone}(\text{representatives})=C^{*}$ passes.** Sound
  but too strong; measured to fail on every bounded specimen even when complete (§2, `notes.md`).
* **Orbit-saturated certificate $\operatorname{cone}(G\cdot\text{representatives})=C^{*}$ as a runtime
  short-circuit.** Complete and sound, but checking it needs the orbit elements, i.e. every facet of every
  found orbit — that is the full-facet materialization the task forbids.
* **Generic engine on the full ray set with a cheaper input format.** No proof of cost on 49 861 rays; the
  measured Normaliz cap is exactly that failure.
* **Randomised facet hunting with a certificate as the only guard.** The certificate is not a generator:
  a failed check yields a separating point $y\in R_B\setminus C$, and turning $y$ into a *new* facet normal
  requires an objective that separates it from $\operatorname{cone}(B)$; that is the forced-progress lemma,
  which is what §1.4 uses — but on its own it stops after $n$ orbits, because $\operatorname{span}(Gf)$
  saturates. Hence the covering recursion is load-bearing.
