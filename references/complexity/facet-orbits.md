# Facet orbits of a cone under a finite group

This document gives the cost of the facet-orbit computation of a pointed,
full-dimensional rational cone $C = \operatorname{cone}(v_1,\dots,v_r)
\subset \mathbb{Q}^d$ under a finite group $G$ that permutes the rays. It also
covers the stabilizer and Schreier-generator steps and the LLL-ordered
solution step that feed it. Unit costs are in [Unit costs](#unit-costs).

## Source and translation

| Step | Pinned source | Port |
| --- | --- | --- |
| Facets of a cell modulo its stabilizer | `polyhedral_common@1592b246/src_lorentzian/lorentzian_perfect.h:1433-1472` (`LORENTZ_ComputeStabilizer`, `DualDescriptionRecordFullDim(EXT_T, GRPperm, rddo)` at :1452) | `src/sage_indefinite_port/indefinite/facet_orbits.py:406` `_facet_orbit_records` |
| Split rule | `POLY_Heuristics.h:320-328` `StandardHeuristicSplitting` (cited in the port; not in the vendored subset) | `facet_orbits.py:461-471` |
| Direct dual description | `POLY_RecursiveDualDesc.h:1465-1474` `DirectFacetOrbitComputation` | `facet_orbits.py:472-492` (Normaliz `SupportHyperplanes`) |
| Adjacency decomposition | `POLY_RecursiveDualDesc.h:837-855`; `SetFunctionality.g:304-323` | `facet_orbits.py:493-548` |
| Balinski stopping test | `Balinski_basic.h:55-60`; `DualDescriptionAndADM.g:133-145` | `facet_orbits.py:503-515` |
| Orbit equivalence without words | `SetFunctionality.g:174-201,234-260` | `facet_orbits.py:436-448` (`RepresentativeAction`, `OnSets`) |
| Ridge closure (complete family) | none (port certificate route) | `facet_orbits.py:261` `cone_facets` |
| Facet certificates | none | `facet_certificates.py` |
| Cell stabilizer | `Lorentzian.g:1318-1322` | `lorentzian_cells.py:1279` `cell_stabilizer` via `cell_transporter` |
| Cell faces and flags | `Lorentzian.g` cell faces | `lorentzian_cells.py:3670` `_CellFaces` |
| LLL-ordered CVP solutions | Lorentzian.g dual-LLL closest-vector step | `lorentzian_cells.py` `_dual_lll_ordered_solutions` |
| Integral Schreier family | `IndefiniteForms.g` stabilizer of a lattice in an envelope | `recursive.py` `_integral_stabilizer_family`, `_integral_stabilizer_action_columns` |

## Invariants that set the cost

- $r$: number of rays. $d$: dimension. $F$: number of facets.
- $|G|$: order of the permutation group on the rays. $N_{\mathrm{orb}}$: number
  of facet orbits. For a representative $f$, $|G_f| = |G| / |Gf|$.
- $\ell$: recursion level. The cell cone of a rank-$n$ lattice enters at
  $\ell = 1$, and each facet step lowers the dimension by one.

Bounds on $F$:

- Upper bound theorem (McMullen 1970; Ziegler, *Lectures on Polytopes*,
  Thm. 8.23): a $(d-1)$-polytope with $r$ vertices has at most
  $f_{d-2}(C(r,d-1))$ facets, which is $O(r^{\lfloor (d-1)/2 \rfloor})$.
- For a simplicial cone $F = d$. For the $d$-cube $F = 2d$. For the
  hypersimplex $\Delta(k,2)$ with $k \ge 5$, $F = 2k$.

Orbit counts: $N_{\mathrm{orb}} \ge F / |G|$ (orbit–stabilizer). Equality
holds when every facet stabilizer is trivial.

## Exact invariants of the test specimens

These values were computed with the production routes in the Sage container:
`LorentzianPerfectLocalBackend.initial_cell(lattice(gram), NONNEGATIVE_SQUARES)`,
`cell_stabilizer_permutations`, `cone_facets`,
`_cell_facet_orbits_from_permutations`, and
`perfect_complex(..., retain_facet_orbits=True)` for the quotient.

| Lattice | $n$ | $r$ | $F$ | rays per facet | $\lvert G\rvert$ | ray-action generators | facet orbits | quotient cells | adjacencies | $O(L)$ generators | branch |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $A_2 \oplus \langle 1\rangle \oplus \langle -1\rangle$ | 4 | 27 | 6 | 9 | 48 | 4 | 1 (size 6) | 1 | 1 | 5 | direct: rule 3, because $r < 40$ |
| $A_3 \oplus \langle 1\rangle \oplus \langle -1\rangle$ | 5 | 89 | 24 | 11 | 384 | 5 | 1 (size 24) | 1 | 1 | 6 | adjacency decomposition: rule 3; the facet cone (11 rays, $d=4$) is direct at level 2 |

For $A_3 \oplus \langle 1\rangle \oplus \langle -1\rangle$, the facet
stabilizer has order $384/24 = 16$.

`Result_Enumeration_0` from `tests/fixtures/lorentzian_perfect_domains.json`
has rank 3, 4 rays, 4 facets, $|G| = 8$, one facet orbit of size 4, facet
stabilizer order 2 and one maximal flag orbit. These values are asserted in
`tests/test_cell_faces_facet_orbits.py`.

The cube and hypersimplex specimens have closed-form invariants:

- $d$-cube: $r = 2^d$, $F = 2d$, $|B_d| = 2^d d!$, $|S_d| = d!$.
- $\Delta(8,2)$: $r = 28$, $F = 16$. A cyclic shift of the moving axes
  together with a transposition generates $S_7$ (order 5040) or $S_8$
  (order 40320).

Counts outside this document:

- Edgewalk vertex orbits and wall counts for $U\oplus\langle 46\rangle$,
  $U\oplus\langle 4\rangle$, $U\oplus\langle 2\rangle$ and
  $U\oplus\langle 2\rangle\oplus\langle 2\rangle$: `vinberg-and-edgewalk.md`,
  edgewalk branches.
- Covering vectors and Eichler-step counts for
  $U(2)\oplus U\oplus\langle 2\rangle$: `orthogonal-group-generators.md`,
  Eichler generation step.
- Quotient traversal and cell transporters:
  `perfect-domains-and-lorentzian-cells.md`.

## Branch selection

The split rule at `facet_orbits.py:461-471`, in rule order:

1. $r - d < 20$: **direct**.
2. $|G| > 5000$ and $\ell < 3$: **adjacency decomposition**.
3. Otherwise: adjacency decomposition iff $r \ge 40$, $\ell \le 1$ and
   $|G| \ge 100$; else **direct**.

$d = 2$ has a separate branch (`facet_orbits.py:451-459`): the two boundary
rays of the pencil, with no dual description.

The rule selects how the computation runs. It does not bound the search.
Both branches return the complete facet set modulo $G$.

## Cost of each branch

### Direct branch

One exact dual description, then one certification per facet, then orbit
grouping.

$$T_{\mathrm{direct}} = T_{\mathrm{Nmz}}(r,d) + F\,(c_{\mathrm{rank}} + d\,r\,c_{\mathrm{py}}) + T_{\mathrm{group}}$$

- With `retain_words=True`, $T_{\mathrm{group}} = F\,|G\text{-gens}|\,|f|\,c_{\mathrm{py}}$
  for the orbit closures (`_orbit`).
- With `retain_words=False`, $T_{\mathrm{group}} \le F\,N_{\mathrm{orb}}\,c_{\mathrm{gap}}$:
  one `RepresentativeAction` per pair of a candidate and a known
  representative of equal incidence size, plus one `Stabilizer` per new
  orbit.

Normaliz is double description with Fourier–Motzkin elimination in the
worst case (Bruns–Ichim, *J. Algebra* 324 (2010) §3). For $r \le 64$ and
$d \le 8$, intermediate facet lists stay below a few thousand rows.

### Adjacency decomposition

Start from one seed facet. Take the pending representative of smallest
incidence. Compute the facet orbits of that facet's cone under the facet
stabilizer $G_f$ (recursion at level $\ell+1$, dimension $d-1$). Rotate each
ridge into its adjacent facet (`_normal_cone_through_ridge`). Register the
result only if its orbit is new.

$$T_{\mathrm{adm}}(C,G) = \sum_{\text{orbit reps } f} \Big[T(C_f, G_f) + R_f\,(c_{\mathrm{rank}} + c_{\mathrm{gap}})\Big] + N_{\mathrm{orb}}\,T_{\mathrm{stab}}$$

Here $R_f$ is the number of ridge orbits of $f$ under $G_f$, and
$T_{\mathrm{stab}}$ is one GAP set stabilizer. Each representative is
expanded at most once. The Balinski test stops the loop early: Balinski's
theorem (Ziegler, Thm. 3.14) says the facet graph of a $(d-1)$-polytope is
$(d-1)$-connected. So when the pending orbits have total size at most
$d-2$, or when they share a full point orbit, they cannot disconnect the
known facets from an unknown facet.

The cost is linear in $N_{\mathrm{orb}}$, not in $F$. This is the reason for
the branch: a large symmetric cone such as a cell with $|G| \gg F$ is reduced
to a few recursive calls on smaller cones.

### Ridge closure (`cone_facets`)

This is the complete facet family, for certificates and tests. It has no
group. Each facet is expanded once, and each of its ridges gives one pencil
rank computation:

$$T_{\mathrm{closure}} = \sum_f R(f)\,(c_{\mathrm{rank}} + r\,d\,c_{\mathrm{py}}),\qquad \sum_f R(f) = 2\,f_{d-3}(C)$$

because each ridge lies on exactly two facets.

## Group steps

- **Cell stabilizer** (`cell_transporter(..., all_images=True,
  generators_only=True)`). This is a backtrack over images of a rational
  basis chosen from the rays. Level $i$ has at most $r$ candidates. These
  are pruned by integral pairing invariants and by the pointwise stabilizer
  of the images already assigned. The search returns $O(\log_2 |G|)$
  generators: each new generator at least doubles the subgroup (Lagrange).
  Cost: $O(n\,r^2)$ pairing checks for the invariants, plus
  $O(\log_2|G|)$ certified isometries of $n^2$ pairings each.
- **Ray action** (`ray_permutations_from_stabilizer`). There are
  $|\text{gens}|\cdot r$ owned isometry applications.
- **Schreier generators.** For an action with $m$ right cosets and $k$
  group generators, Schreier's lemma (Seress, *Permutation Group
  Algorithms*, Lemma 4.2.1) gives at most $m\,k$ generators of the
  stabilizer. Each is a word of length $O(m)$. Each is realized on the
  rank-$n$ basis with $2n$ owned images (forward and inverse). Examples:
  $U\oplus U$ with the trivial bound has $m = 1$; $\langle 1\rangle^3$ with
  the $S_3 \times \{\pm1\}$ action on an index-6 sublattice has $m = 6$,
  $k = 3$, and hence 18 generators.
- **Flags of a cell** (`_CellFaces.flag_simplices`). One maximal flag is
  produced per orbit of $G$ on maximal flags. The descent computes one set
  stabilizer per level, $d$ levels.

## LLL-ordered solutions

`_dual_lll_ordered_solutions` sorts $s$ given solutions by their recursive
coordinates under a fixed $n\times n$ transform. The cost is
$O(s \log s \cdot n^2)$ exact rational operations. No enumeration happens
here. The enumeration volume belongs to the fixed-distance search that
produces the solutions; see `perfect-domains-and-lorentzian-cells.md`.

## Worst cases and termination

- Each branch terminates: the facet set is finite and each orbit is
  registered once.
- The direct branch on a cone with many rays and little symmetry has the
  upper-bound-theorem worst case $O(r^{\lfloor (d-1)/2\rfloor})$ facets. The
  split rule sends such a cone to the direct branch only when the group is
  too small to help.
- Adjacency decomposition on a group with large facet stabilizers recurses
  to level 3 before rule 2 stops applying. The recursion depth is at most
  $d-2$.
- A non-pointed input (nontrivial lineality space) has no seed slice. The
  walk raises (`facet_orbits.py`, seed and pencil certification) and does
  not return a partial family.

## Unit costs

| Primitive | Assumed cost on a 4-vCPU runner |
| --- | --- |
| $c_{\mathrm{rank}}$: FLINT exact rank/kernel, $\le 64\times 8$ | 0.3 ms |
| $T_{\mathrm{Nmz}}$: Normaliz dual description, $r \le 64$, $d \le 8$ | 50 ms |
| owned `from_rays` cone (PPL), $r \le 30$, $d \le 5$ | 100 ms |
| $c_{\mathrm{gap}}$: libGAP call, degree $\le 64$ | 2 ms |
| owned lattice element operation (`sage_extensions` dispatch) | 2 ms |
| pure-Python set or tuple step $c_{\mathrm{py}}$ | 1 µs |
| one `O(L).transporter` of a vertex configuration, rank 3–4 | 0.5 s |
| one `O(L).stabilizer` of a vertex configuration, rank 3–4 | 0.5 s |
| one edgewalk `edge_step`, rank 3–4 | 0.3 s |
| owned lattice construction with its first `O()` parent | 0.1 s |
| owned isometry from images, rank $\le 7$ | 50 ms |
| fresh child interpreter importing Sage and the port (only where a test spawns one) | 8 s |

A budget is the sum of these primitive costs over the work of the selected
branch. It has no multiplicative factor and no floor. A budget that a run
exceeds therefore shows extra work or a slower primitive.
