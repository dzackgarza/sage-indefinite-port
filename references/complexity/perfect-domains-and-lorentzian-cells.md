# Perfect domains and Lorentzian cells: complexity

This document gives the cost of the Lorentzian perfect-domain algorithm and of
its local operations. The costs follow from the mathematics of the input. They
are not measurements. Facet-orbit enumeration and Schreier generators are in
`facet-orbits.md`. Isotropic-vector search is in `isotropic-vector-search.md`.
Generic isometry testing is in `isometry-testing.md`.

## Algorithm and pinned source

Let `L` be an even or odd integral lattice of signature `(1, n)` and rank
`r = n + 1`. Fix a locus `S` of allowed squares (`ZERO_SQUARES` or
`NONNEGATIVE_SQUARES`). A *Lorentzian perfect cell* is the finite family `V`
of primitive vectors `v` with `v^2 ∈ S` that attain the minimum of the linear
pairing with an interior timelike vector. The cells tessellate the positive
cone, and `O(L)` acts on the tessellation with finitely many cell orbits.
The algorithm is Voronoi's algorithm transported to the Lorentzian cone:
find one perfect cell, enumerate its facets up to its stabilizer, flip across
each facet, and test each new cell for equivalence with the known
representatives. The stabilizers and the flip transporters generate `O(L)`.

| Step | Pinned source | Port |
| --- | --- | --- |
| Positive vectors of a pairing shell | `polyhedral_common@1592b246/src_lorentzian/lorentzian_perfect.h:107`, `:397` | `_PositiveVectorEnumerator`, `_fixed_distance_dual_lll_solutions` (`lorentzian_cells.py:4235`), `_bounded_distance_dual_lll_solutions` (`:4287`) |
| Outside ray and initial cell | `lorentzian_perfect.h:840`, `:969` | `perfect_complex` (`lorentzian_cells.py:420`) |
| Flip across a facet | `lorentzian_perfect.h:666` | `LorentzianPerfectLocalBackend.flip_across` (`lorentzian_cells.py:1809`) |
| Cell stabilizer | `lorentzian_perfect.h` `LORENTZ_ComputeStabilizer`; `Indefinite.jl@374a5ebb/indef/lib/Lorentzian.g` | `cell_stabilizer` (`lorentzian_cells.py:1279`) |
| Cell transporter | `lorentzian_perfect.h:1176` (`LinPolytopeIntegral_Isomorphism_GramMat`); `Lorentzian.g:660-696` | `cell_transporter` (`lorentzian_cells.py:1975`) |
| Group generators | `lorentzian_perfect.h:1531`, `:1574`, `:1591` | `lorentzian_orthogonal_group_generators` (`lorentzian_cells.py:445`) |
| Cell cone | none (port representation) | `cell_cone` (`lorentzian_cells.py:1061`) |

## Notation for the invariants

| Symbol | Meaning |
| --- | --- |
| `r` | rank of `L` |
| `N` | number of vectors in the cell configuration, `|V|` |
| `F` | number of facets of the cell cone |
| `G_C` | integral stabilizer of the cell in `O(L)` |
| `h` | number of `O(L)`-orbits of perfect cells |
| `w` | Witt index of `L ⊗ Q` (`w = 1` for an isotropic Lorentzian lattice, `w = 0` otherwise) |

## Positive-vector shell

The shell is the set of `v ∈ L` with `v^2 ∈ S` and `⟨v, t⟩ = m` for a fixed
timelike `t` and the least attained `m`. On the hyperplane `⟨v, t⟩ = m` the
form restricted to `t^⊥` is negative definite. So the shell is the set of
lattice points of an affine lattice of rank `r - 1` inside an ellipsoid. This
is a closest-vector problem.

Branches (selected by `S` and by `t^⊥`):

- **Fixed distance (`S = ZERO_SQUARES`, isotropic route).** The points lie on
  one sphere of the definite form. Cost: the Fincke-Pohst enumeration tree.
  The inner coordinate is solved directly, so the leaf level costs one
  square-root test per interior node.
- **Bounded distance (`S = NONNEGATIVE_SQUARES`).** The points lie in a ball.
  Cost: the number of lattice points in the ball plus the interior nodes of
  the Fincke-Pohst tree.

Cost bound. For a definite lattice `M` of rank `k = r - 1` and radius `R`, the
number of points is at most `vol(B_k(R)) / covol(M)` plus a boundary term
(Gauss counting). Fincke and Pohst (Math. Comp. 44 (1985), 463-471) bound the
tree by the sum over levels of such counts in projected lattices. A dual-LLL
reduced basis (Lenstra, Lenstra and Lovász, Math. Ann. 261 (1982), 515-534)
makes the projected covolumes within a factor `2^{O(k^2)}` of each other.
For `k ≤ 3` and small `m`, the tree has at most a few hundred nodes.

The initial bound is advanced without a rescan
(`test_initial_vector_search_advances_integral_bound_without_rescanning`).
Non-termination: none. The minimum `m` exists because `L` is discrete and `t`
is timelike.

## Initial cell

The initial cell needs one interior timelike vector and one positive-vector
shell. Cost: one exact LLL of rank `r`, one symmetric diagonalisation, and one
shell enumeration. These are polynomial in `r` and in the bit size of the
Gram matrix.

## Cell cone

`cell_cone` (`lorentzian_cells.py:1061`) has two branches:

- **Ray route, `N < 10 000`.** The owned cone is built from rays. The
  backend computes the H-description by double description. The cost is
  output-sensitive: `F` facets, and `F` is bounded by the upper bound
  theorem, `F = O(N^{⌊(r-1)/2⌋})` (McMullen, Mathematika 17 (1970),
  179-184).
- **Normaliz route, `N ≥ 10 000`.** Normaliz returns support hyperplanes
  and equations. The port then makes one linear pass per output covector.
  Duplicated input rays are removed by Normaliz. The cost is the Normaliz
  convex-hull cost plus `O((N + F) r)` exact integer operations.

Both routes give the same cone
(`tests/test_large_cell_cone_route.py`). The threshold is a cost choice, not
a mathematical branch.

## Facets and flips

The facet incidences are the zero sets of the `F` facet covectors on `V`. Cost
`O(F N r)` integer operations. For `N > 4096` the port uses packed bitsets and
the smaller incidence side (`_facet_mask_data`, `lorentzian_cells.py:1417`).
Facet orbits under `G_C` are in `facet-orbits.md`.

`flip_across` computes the kernel of the facet configuration, a direction off
the facet, and the next positive-vector shell along that direction
(`lorentzian_perfect.h:666`). Cost: one exact nullspace of an `N_f × (r+1)`
matrix (`N_f` = facet size) and one positive-vector shell. A flip can leave
the cone at an isotropic ray. Then the facet lies on the boundary at a cusp,
and there is no neighbour. This happens only for `w = 1`.

## Cell stabilizer and transporter

`cell_transporter(source, target)` (`lorentzian_cells.py:1975`) decides if an
integral isometry maps `V_source` onto `V_target`. Branches, in order:

| Branch | Condition | Work |
| --- | --- | --- |
| Rejection | `N`, rank, `det`, signature, Smith form, moment determinant, row-sum content or norm differ (`:1999-2109`) | `O(N r^2)` exact operations plus one Smith form of rank `r` |
| Identity | same lattice and same rows (`:2020`) | `O(1)` |
| Small graph | `N ≤ 256` (`_PAIRING_GRAPH_CANONICAL_LABEL_LIMIT`, `:910`), witness only, no cached stabilizer | coloured pairing-graph isomorphism (Lorentzian.g:660-696, AutomorphismPolytope.g:288-339) |
| Basis search | all other cases, including `all_images=True` and stabilizers | colour refinement (when `N ≤ 4096`, `:2142`), then backtracking over images of a source basis |

**Small graph.** The graph has `N` vertices and `N^2` integer edge weights.
Canonical labelling (McKay and Piperno, J. Symbolic Comput. 60 (2014), 94-112)
is exponential in the worst case. For vertex-transitive cell graphs it is
polynomial in practice, with work `O(N^2 log N)` per refinement round.

**Basis search.** Choose `r` source vectors that span `V ⊗ Q`. Each leaf
assigns `r` target vectors and is certified by one congruence test and one
row-set test, `O(N r^2)`. Colour refinement leaves only colour-compatible
images. When all images are requested, every element of `G_C` (or of the
transporter coset, of size `|G_C|`) is one accepted leaf. So the leaf count is
at least `|G_C|`. With an equitable partition that separates the basis
vectors, the leaf count is `O(|G_C|)`. This is the Plesken-Souvignier bound
for the automorphism search with a fingerprint (Plesken and Souvignier,
J. Symbolic Comput. 24 (1997), 327-334). For `N` vectors in one colour class,
the worst case is `N^r` leaves.

Stabilizer orders for the cells in the tests:

| Cell | `N` | `|G_C|` | Reason |
| --- | --- | --- | --- |
| `U`, rows `e, f` | 2 | 2 | swap of the two isotropic lines |
| `A2 ⊕ A1(-2)`, roots of `A2` lifted | 6 | 12 | `W(A2) × {±1}` |
| `A3 ⊕ A1(-2)`, roots of `A3` lifted | 12 | 48 | `W(A3) × {±1}`; 8 triangular and 6 square facets |
| `3A1 ⊕ A1(-1)`, `t ± e_i` | 6 | 48 | signed permutations of 3 axes, `2^3 · 3!` |
| `diag(5,-1,-2)`, shell `(1,x,y)` | 11 | 4 | sign changes of `diag(-1,-2)` |
| rank-10 simplicial cell | 10 | 1 | a simplicial cell with distinct colours |

## Traversal

The traversal visits `h` cell orbits. For each orbit it does one stabilizer,
`f_C` facet-orbit flips, and up to `h` transporter tests per flip (bucketed by
`_perfect_form_hash_key`, `lorentzian_cells.py:2860`). Total work:

`Σ_C [stab(C) + f_C · (flip + b · transporter)]`,

where `b` is the bucket size. Finiteness of `h` follows from reduction theory
of the self-dual positive cone (Koecher, Math. Ann. 141 (1960), 384-432; Ash,
Mumford, Rapoport and Tai, *Smooth Compactifications of Locally Symmetric
Varieties*, 2nd ed., Cambridge 2010, Chapter II). No a-priori bound on `h` in
terms of `r` and `det L` is cited here. For rank 2 and 3 with small
determinant, `h` is in the single digits.

The traversal terminates because `h` is finite. Non-termination is possible
only from a defect: a missed transporter makes a known orbit look new, and
the traversal then repeats the orbit without end.

## Large-configuration route

When `N > 4096` (`_PAIRING_GRAPH_REFINEMENT_LIMIT`, `lorentzian_cells.py:903`),
the hash key and the transporter use only linear invariants: pairing diagonal,
gcd colours, row-sum data. Cost `O(N r)` per invariant. The quadratic pairing
graph is not built. Known instance: the `A7` perfect configuration with
49 861 rays.

## Worked tessellations

### `U ⊕ A1(-2)`, `S = NONNEGATIVE_SQUARES` or `ZERO_SQUARES`

Write `v = a e + b f + c g` with `e·f = 1`, `e^2 = f^2 = 0`, `g^2 = -2`. Then
`v^2 = 2ab - 2c^2`. The vectors with `v^2 ≥ 0` in the closure of the positive
cone are those with `a, b ≥ 0` and `ab ≥ c^2`.

Take `t = 2e + 2f + g`, with `t^2 = 6 > 0`. Then `⟨v, t⟩ = 2(a + b - c)`. From
`c ≤ √(ab) ≤ (a + b)/2` it follows that `a + b - c ≥ (a + b)/2 > 0`. So
`a + b - c ≥ 1`. Equality needs `a + b ≤ 2`. If `a + b = 1`, then `c = 0`,
which gives `e` or `f`. If `a + b = 2`, then `c = 1` and `ab ≥ 1`, which gives
`e + f + g`. So the cell is `C = {e, f, e + f + g}`. All three vectors are
isotropic, so both loci give the same cell. A brute-force check over
`|a|, |b|, |c| ≤ 6` (plain Python, integers only) gives the same minimal set.

The three vectors have determinant 1 in the basis `(e, f, g)`, so they form a
basis of `L`. In that basis the Gram matrix has zero diagonal and all
off-diagonal entries equal to 1. So every permutation of the three vectors is
an integral isometry. Negation maps the positive cone to the negative cone.
Thus `G_C ≅ S_3`, `|G_C| = 6`, and `G_C` acts transitively on the three
facets. The flip across the facet `{e, f}` gives `{e, f, e + f - g}`. This is
the image of `C` under `g ↦ -g`. So there is `h = 1` cell orbit and one facet
orbit. This is the Farey tessellation of the hyperbolic plane by ideal
triangles: each vertex is a cusp.

Traversal work: one initial shell, one stabilizer search (`N = 3`,
`|G_C| = 6`), one flip, and one transporter test.

### `U`, `S = ZERO_SQUARES`

The only primitive isotropic vectors in the positive cone are `e` and `f`. The
cell is `{e, f}`, with `|G_C| = 2` (the swap). Both facets `{e}` and `{f}` are
isotropic rays, so each flip exits at a cusp and has no neighbour. Thus
`h = 1` and the traversal has no flips.

## Helpers without search

Some tests check a helper that does no search: block iterators, Schreier word
streams, bitset permutations, pairing diagonals, packed fibres, colour
invariants, sign flips, basis reflections, and the invariant cache. Each such
helper is one linear pass over its input. Cost `O(N r)` exact operations,
plus `LAT` per owned lattice that the test builds.

## Unit costs of backend primitives

These are the costs assumed for small inputs (`r ≤ 10`, entries of a few
digits) on a 4-vCPU runner with Sage 10.x:

| Primitive | Cost |
| --- | --- |
| owned lattice construction (sage-extensions parent, Gram check) | 0.1 s |
| one libgap call (group order, `NullspaceMat`, LLL), degree or dim `≤ 16` | 0.1 s |
| one PyNormaliz cone, rank `≤ 4`, few rays, with conversion | 0.1 s |
| owned `Isom.transporter`, graph on `≤ 16` vertices | 0.2 s |
| configuration invariants (Smith form, adjugate), `r ≤ 10` | 0.05 s |
| one positive-vector shell, rank `≤ 3`, a few hundred tree nodes | 0.2 s |
| one flip (facet kernel and one shell) | 0.3 s |
| one basis-search leaf certification | 5 ms |
| one exact small-matrix routine (LLL, kernel, diagonalisation), dim `≤ 5` | 10 ms |
| one row in a Python linear pass | 30 µs |
| one Sage `QQ` operation in a Python generator loop | 3 µs |

Library start-up (libgap, Normaliz) happens once per session and is not part
of the work of one test. A budget is the derived work of its branch, with no
multiplicative factor.
