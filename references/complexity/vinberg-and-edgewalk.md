# Vinberg's algorithm and Allcock's edgewalk

This file states the cost of computing the simple roots of a hyperbolic lattice
`L` of signature `(n-1, 1)` and of deciding reflectivity, by Allcock's edgewalk.
The stabilizers and transporters of vertex configurations in `O(L)` are inputs
here; `isometry-testing.md` owns their cost. The discriminant data `A_L`,
`O(q_L)` are owned by `discriminant-groups-and-congruence-orbits.md`.

No measured runtimes are used. Each cost is a count from the mathematics times
a stated unit cost of a backend primitive.

## Sources

| Step | Pinned source | Port |
| --- | --- | --- |
| Driver, root norms | `polyhedral_common@1592b246/src_lorentzian/edgewalk.h:55-67, 2092-2238` | `src/sage_indefinite_port/indefinite/edgewalk.py:58-102, 249-275` |
| Initial vertex, isotropic start (`isotropic_vinberg`) | `edgewalk.h:1880-2055, 2119` | `edgewalk.py:343-464` (`initial_isotropic_vertex`) |
| Initial vertex, anisotropic start (Vinberg increasing prefixes) | `edgewalk.h:2092-2135`, `vinberg.h:139-141` | `edgewalk.py:465-604, 710-908` |
| Edge directions | `edgewalk.h:525-539` (`GetAdjacencyDirection`) | `edgewalk.py:1358-1380` |
| Next vertex along an edge | `edgewalk.h:966-1061` | `edgewalk.py:1186-1356`, `edgewalk_extensions.py` |
| Main walk, finiteness test, `f_try_terminate` | `edgewalk.h:1677-1792` | `edgewalk.py:1509-1618` |
| Rank 2 (oriented binary cycle) | `two_dim_lorentzian.h` | `edgewalk_rank2.py`, `edgewalk.py:1484-1507` |
| Algorithm | Allcock, arXiv:2110.03784 (`references/vendor/arxiv/2110.03784`) | |

## Upstream check

`polyhedral_common@1592b246/.github/workflows/ci_20_reflective.yml` runs one
job with `timeout-minutes: 120`. The job installs GAP, builds the C++ and runs
`CI_tests/20_Reflective/TestReflectivity.g`, whose `FullTest` computes the
simple roots and the isotropy of all 8821 `ListReflect` lattices (8602 of rank
3). The C++ edgewalk therefore needs less than `120 * 60 / 8821 = 0.82 s` per
lattice on average.

This bound is a check on the derived budgets, not a multiplier: the family
average must be below it. Each step is priced at the compiled cost of its
operation (Cython over FLINT, PARI and GAP; vectorized FLINT matrices), which
is the work the C++ does. Interpreted per-element work in a hot loop of the
port is not priced in; a run above the budget exposes it.

## Unit costs of backend primitives

The units are the cost of each operation at the fixtures' sizes (rank `<= 18`,
entries `< 2^25`) in one shared Sage session.

| Symbol | Operation | Unit cost | Basis |
| --- | --- | --- | --- |
| `OP` | One FLINT or PARI operation on a rank `<= 18` integer vector or matrix | 3e-6 s | At most 324 word-size entries; `O(n^2)` to `O(n^3)` word operations |
| `STAB(n)` | Stabilizer or transporter of a coloured family of `<= n+1` vectors | `2e-4 * max(1, n^2/9) * (1 + log2|O(q_L)| / 64)` s | Plesken–Souvignier backtrack (PARI `qfauto`/`qfisom`) over `O(10)` short vectors of a rank-3 definite form, in C |
| `LP` | One halfspace-ray incidence in a double description | 1e-6 s | cdd in C arithmetic, Fukuda–Prodon 1996 |
| `BASE` | Smith form, root norms, isotropy decision, one isotropic vector | 2e-3 s | PARI `matsnf`, `qfsolve` (Simon 2005) |

**Compiled work per step** (from `edgewalk.py`):

- per vertex orbit: the coloured configuration of `|S| + 1` vectors and its key
  (edgewalk.py:1389-1402), `(n + 2) OP`;
- per transporter comparison rejected by invariant (edgewalk.py:1584-1588):
  one Gram-key comparison, `OP`;
- per edge and root norm (edgewalk.py:1340-1356): the binary search of length
  `s = log2(2 exp(A_L) |A_L|)`, one rank-`n` matrix-vector product per step,
  plus setup: `(30 + 2 s n) OP`.

Each budget is the derived work of the branch with these units.

## Root norms

A root `r` of `L` is primitive with `r^2 > 0` and `2 b(r, x)/r^2 in Z` for all
`x in L`. Then `r/ (r^2/2)` lies in `L^*`, so `r^2` divides `2 exp(A_L)`
(Vinberg 1972, §1; edgewalk.h:55-67). The port takes
`lattice.possible_root_lengths()` (edgewalk.py:262). The count of root norms
`N = tau(2 exp(A_L))` multiplies every per-edge search.

## Edgewalk branches

The branch is selected by the arithmetic of `L`, not by its size.

| Branch | Condition | Port, source |
| --- | --- | --- |
| Binary | `n = 2` | `binary_walk`, edgewalk.py:1484 |
| Cusp start | `L` is isotropic: always for `n >= 5` (Meyer); for `n = 3, 4` iff the Hasse–Minkowski local conditions hold at every `p | 2 det L` | edgewalk.py:716, edgewalk.h:1880-2055 |
| Cocompact start | `L` is anisotropic (`n <= 4` only) | edgewalk.py:717-908 |
| Reflective stop | the root set closed under the found symmetries bounds a polyhedron (`root_chamber_certificate`) | edgewalk.py:1556-1567, edgewalk.h:1700 `f_try_terminate` |
| Nonreflective stop | the symmetries found fix no common nonzero vector (an infinite group) | edgewalk.py:1528-1544, edgewalk.h:1697-1792 |

### Binary branch

`L` of rank 2 and signature `(1, 1)`. The walk is the reduction cycle of the
binary form (Gauss reduction). Its length is at most the period of the
continued fraction of the roots of the form, `O(sqrt(D) log D)` with
`D = -det L` (Lagarias 1980; Buchmann–Vollmer, *Binary Quadratic Forms*,
Thm. 6.10.2 for the cycle length bound). Each step is `O(1)` binary operations.
Isotropic binary forms (`D` a square) have no roots (Vinberg 1972); the walk
returns empty.

### Cusp branch

The walk starts at a primitive null vector `k`. The definite reduction
`k^perp/k` of rank `n-2` gives the affine root system at the cusp. One shell
per root norm is enumerated in `k^perp/k`; the number of vectors per shell is
the theta coefficient, bounded by `O(norm^((n-2)/2) / sqrt(det))`
(Conway–Sloane, SPLAG, Ch. 2 §1.4). For the fixtures the reduction has small
minima and the shells are of size `O((n-2)^2)`.

Ideal vertices recur at each cusp of the polygon or polyhedron. The number of
cusp orbits of a reflective `L` equals the number of `O(L)`-orbits of
primitive isotropic lines, which `isotropic-vector-search.md` bounds.

### Cocompact branch

The start uses Vinberg's increasing prefixes from a timelike `x_0`: roots are
enumerated by the height `b(r, x_0)^2 / r^2` (Vinberg 1972, §3). The prefix
search ends at the first corner. The number of shells searched grows with the
distance from `x_0` to the nearest corner, bounded by the diameter of the
fundamental polygon; per norm this gives a factor `O(sqrt(|A_L|))` in the
enumeration volume (the covolume of `L` is `sqrt(|det L|)`).

### Main walk

Let `R` be the number of simple roots (walls) and `P` the Coxeter polyhedron.

**Vertices.** `P` is acute-angled. Its finite vertices are the elliptic
subdiagrams of rank `n-1` of its Coxeter diagram, and its ideal vertices are
the parabolic subdiagrams of rank `n-2` (Vinberg, *Hyperbolic reflection
groups*, Russian Math. Surveys 40 (1985), Thm. 3.1 and Prop. 4.2). When the
simple roots are recorded (Allcock's rank-3 table, `root_systems_103`, the
OSCAR chamber, E10), the vertex set is enumerated from the root Gram matrix:
an `(n-1)`-subset with positive-definite Gram is a finite vertex; an
`(n-1)`-subset with positive-semidefinite Gram of rank `n-2` lies in an ideal
vertex, identified by the null ray of its kernel combination. For `n = 3` this
gives `V = R` (a polygon); for E10 it gives the 10 vertices of the simplex,
one of them ideal.

**Vertex orbits.** The walk visits one vertex per orbit of
`Sym_L(P) = O^+(L) / W`, the chamber automorphisms that are integral on `L`
(Vinberg 1972, §1: `O^+(L) = W ⋊ Sym_L(P)` for reflective `L`). `Sym_L(P)` is
the group of automorphisms of the labelled Coxeter diagram (root norms and
pairings) whose induced linear map is integral on `L`. It is computed once per
recorded root system. For `root_systems_103` its order is 1, 2, 4 or 12. The
vertex orbit count `V` is the number of orbits of this group on the vertex set.

When no roots are recorded (the 112 `ListReflect` forms of rank `>= 6` outside
`root_systems_103`), the estimate takes `V <= R` vertex orbits. This is an
assumption, not a theorem; it fails for a chamber with a small symmetry group
and more than `R` vertices.

**Edges.** A vertex with wall set `S` has `|S|` outgoing edges in the walk
(`edge_directions`, edgewalk.py:1358-1380: one per facet of the cone of `S`).
`E = sum over vertex orbits of |S_v|`.

**Transporter tests.** edgewalk.py:1584-1588 compares the neighbour of each
edge with the known vertex configurations in order and stops at the first
match. A known configuration with a different invariant (number of walls, the
multiset of rows of the root Gram) is rejected by comparing Gram keys
(`OP`). A full transporter (`STAB`) happens only against configurations of the
same invariant class. So the walk makes at most
`sum_v |S_v| * |class(v)|` full tests and `E V` cheap rejections.

Per vertex orbit: one `STAB(n)` for the stabilizer of the coloured
configuration (edgewalk.py:1574).

Per edge: one `finite_extensions` and one `extension_roots` call per root norm
(edgewalk.py:1340-1356). The source takes the first integral candidate in a
binary search (edgewalk.h:966-1061) of length `s`. The cost per norm is the
compiled work listed under the unit costs.

Per vertex orbit: at most one dual description with `R` halfspaces and `V`
rays, `50 OP + LP R V` (edgewalk.py:1556-1567). The certificate can close
before all orbits are visited; the estimate does not take that credit.

Total:

```
T = BASE + T_start
  + V (STAB(n) + (n + 2) OP) + (sum_v |S_v| |class(v)|) STAB(n) + E V OP
  + E N (30 + 2 s n) OP
  + V (50 OP + LP R V)
```

with `T_start = (1 + cusps) N 40 OP (n-2)^2` on the cusp branch and
`T_start = N (40 OP n^2 + 1e-5 sqrt|A_L|)` on the cocompact branch.

The dominant terms are the full transporter tests and the edge searches.
`|A_L|` enters through `N`, the binary search length and `log2 |O(q_L)|` in
`STAB`. A rank-3 lattice with large `|A_L|` and many primes in `exp(A_L)` costs
more than a unimodular lattice of higher rank with the same `R`.

The unrecorded rank-18 unimodular form (`ListReflect_8820`, `R = 22`) has the
largest bound: `STAB(18) = 7.2e-3 s` and up to `2 (n-1) R` full tests. With its
Coxeter diagram (not recorded in the fixtures) the vertex orbits and classes
would be exact.

## Worst cases and nontermination

- Vinberg's algorithm does not terminate on a nonreflective lattice: the
  root heights are unbounded (Vinberg 1972, §4). The edgewalk terminates on
  every hyperbolic lattice: it stops when the found symmetries generate an
  infinite group (Allcock 2021, Thm. 1.1). The nonreflective stop needs only
  the vertex orbits visited until that point; `U + <46>` stops after three
  vertex orbits.
- With `early_termination=False` (`edgewalk_isometry`, edgewalk.py:105-140)
  a nonreflective walk visits every vertex orbit of a chamber whose symmetry
  group is infinite; the count is finite because `O(L)` has finite covolume,
  but it is bounded only by the covolume over the minimal vertex-orbit volume.
- Reflective hyperbolic lattices exist only for `n <= 22` (Esselmann 1996),
  and none have `n = 21`. The largest fixtures have `n = 18`.

## References

- Allcock, *Root systems for Lorentzian Kac–Moody algebras in rank 3*, arXiv:2110.03784.
- Vinberg, *On the groups of units of certain quadratic forms*, Math. USSR Sb. 16 (1972) 17–35.
- Conway–Sloane, *Sphere Packings, Lattices and Groups*, 3rd ed., Ch. 2.
- Esselmann, *Über die maximale Dimension von Lorentz-Gittern mit coendlicher Spiegelungsgruppe*, J. Number Theory 61 (1996).
- Plesken–Souvignier, *Computing isometries of lattices*, J. Symb. Comput. 24 (1997).
- Fukuda–Prodon, *Double description method revisited* (1996).
- Buchmann–Vollmer, *Binary Quadratic Forms*, Springer 2007.
