# Discriminant groups, congruence subgroups and orbit splitting

This file states the cost of the operations that act through a finite quotient
of `O(L)`: the discriminant group `A_L = L^*/L`, its orthogonal group
`O(q_L)`, discriminant images of stabilizers, congruence representations
`O(L) -> Aut(L/mL)`, their preimages, and the splitting of one `O(L)`-orbit of
markings into orbits of a finite-index subgroup. The generators of `O(L)` and
of parabolic stabilizers are inputs here. Their cost is owned by
`orthogonal-group-generators.md` (§Branches, §Eichler route, §Lorentzian route,
with cell counts from `perfect-domains-and-lorentzian-cells.md`) and by
`boundary-components.md`. Double cosets in
general are owned by `double-cosets.md`; this file covers only their use in
finite images.

No measured runtimes are used. Each cost is a count from the mathematics times
a stated unit cost of a backend primitive.

## Unit costs of backend primitives

| Symbol | Primitive | Upper unit cost on a 4-vCPU runner | Basis |
| --- | --- | --- | --- |
| `c_mor` | One evaluation, composition or construction of a lattice morphism or marked object in the sage-extensions categorical layer (Python objects over a rank `n <= 12` free module) | 0.03 s | `O(n^2)` integer operations plus Python object construction of a few hundred attribute lookups |
| `c_lat` | Construction of a formed lattice from a Gram matrix with its cached invariants (determinant, signature, Smith form) | 0.5 s | Smith form of an `n x n` integer matrix, `n <= 12`, in PARI/FLINT |
| `c_ss(d)` | GAP randomized Schreier–Sims order of a permutation group of degree `d` with `k <= 50` generators | `1 + 1e-7 * k * d * log2(d)^2` s | Seress, *Permutation Group Algorithms*, Sec. 4.5 and Thm. 4.5.5: nearly linear time `O(k d log^c |G|)` for groups with small base; `log |G| = O(log^2 d)` for the classical groups here |
| `c_mat(p, n)` | Order of a matrix group in `GL_n(F_p)` | `c_ss(p^n)` for `p^n <= 10^3`; otherwise `t_reduction` of `orthogonal-group-generators.md` §Finite images | For small `p^n` GAP acts on the vectors. For large `p^n` the vector action has `(p^n - 1)/(p - 1)` points, and constructive recognition is polynomial in `n` and `log p` |
| `c_tqf(A)` | `O(q_L)` for a finite quadratic form, Sage `TorsionQuadraticModule.orthogonal_group` | `2 + 0.02 * |A_L|` s per `p`-primary part | Sage `fqf_orthogonal.py` lifts generators of `O(L_p)` (Brandhorst–Veniani, *Hermitian lattices and K3 automorphisms*; Miranda–Morrison, *Embeddings of integral quadratic forms*, Ch. IV); the permutation action on `A_L` has degree `|A_L|` |

A safety factor of 3 is applied to every total in
`.tmp/probe/estimates/F.py`.

## Discriminant group and its orthogonal group

`A_L` is computed from the Smith form of the Gram matrix: one `c_lat`.
`|A_L| = |det L|`. `O(q_L)` is a product over primes `p | det L` of the groups
`O(q_{L,p})` (Nikulin 1980, Prop. 1.2.1, the orthogonal splitting of `A_L`
into `p`-parts). Its cost is `c_tqf` of each part.

Concrete orders used by the tests:

| Lattice | `A_L` | `O(q_L)` | Reason |
| --- | --- | --- | --- |
| `2U + <-2> + <-6>` | `Z/2 + Z/6` | order 2 | The three elements of order 2 have distinct `q`; the 3-part has `O = {+-1}` |
| `U + <6>` | `Z/6` | order 2 | Cyclic; `O(q) = {+-1}` |
| `U + U + U(3)` | `(Z/3)^2`, hyperbolic | order 4 | `O^+_2(3)` is dihedral of order 4 |
| `2U(2) + A2` | `(Z/2)^4 + Z/3` | `72 * 2 = 144` | 2-part is `O^+_4(2)` of order `2 * 2^2 (2^2-1)(2^2-1) = 72` (Conway et al., *ATLAS*, Ch. 2, orders of classical groups) |
| `N = U + U(2) + E8(-2)` | `(Z/2)^10`, even type `+` | `|O^+_10(2)| = 46998591897600` | `2 q^{n(n-1)} (q^n - 1) prod_{i<n} (q^{2i} - 1)` at `n = 5, q = 2` |

Worst case: `|A_L|` large with a `p`-part of high `p`-rank. The
permutation degree of `O(q_L)` is `|A_L|`, so `c_tqf` grows linearly in
`|A_L|` for a bounded number of generators. A rank-three lattice with
`|A_L| = 2048` costs more here than the rank-twelve lattice `N`, whose
`|A_L| = 1024`.

## Discriminant image of a set of isometries

`discriminant_image_order(L, gens)` maps each generator to `O(q_L)` by its
action on the `rank(A_L)` Smith generators (one `c_mor` each) and computes the
order of the generated subgroup: `c_ss(|A_L|)`.

For `N`, the image of a line stabilizer has index 1 or 527 in `O^+_10(2)`.
527 is the number of nonzero isotropic vectors of the even-type form of
dimension 10 over `F_2`: `(2^5 - 1)(2^4 + 1) = 527`. The order computation is
one Schreier–Sims of degree at most 1023.

`stabilizer_discriminant_image()` of a primitive isotropic sublattice `I`
(sage-extensions `isotropic_parabolics.py:748`) computes the image of
`Stab(I)` from the Levi and unipotent generators of the parabolic. Its cost is
the parabolic generator cost (owned by `boundary-components.md`) plus one
`c_ss(|A_L|)`.

## Congruence representation and preimages

`congruence_representation(m)` (sage-extensions `lattice_morphisms.py:4735`)
reduces each generator of `O(L)` to `GL_n(Z/m)` and forms the image group. The
image is finite. Its order is the cost driver:

| Lattice, modulus | Image | Order |
| --- | --- | --- |
| `2U`, `m = 3` | `O^+_4(3)` | 1152 |
| `2U`, `m = 2` | `O^+_4(2)` | 72 |
| `3U`, `m = 2` | `O^+_6(2)` | 40320 |
| `U + <0>`, `m = 3` | `O(U/3) x {+-1} x F_3^2` | 72 |
| `U + <6>`, `m = 6` | joint image in `Aut(L/6L)` | `|im mod 2| * |im mod 3|` |

Branch: the image is computed in `GL_n(Z/m)` by GAP. Cost
`k * n^2 * c_mor + c_mat(m, n)`, with `k` the number of `O(L)` generators.

`congruence_preimage(m, S)` is `rho^{-1}(S)` (`lattice_morphisms.py:4756`). It
stores the predicate and costs nothing until it is queried.

A stabilizer of `x` inside `Gamma = rho^{-1}(S)` is obtained by Schreier
generators of `Stab_G(x) cap Gamma` (Seress, Lemma 4.2.1). The number of
Schreier generators is `[Stab_G(x) : Stab_Gamma(x)] * k_x`, and the index
equals `|rho(Stab_G(x))| / |rho(Stab_G(x)) cap S|`. For `Gamma(3)` and
`x = Ze` in `2U` this index is 72 (Witness in `witnesses/test_orbits.sage`).
Each Schreier generator costs a membership test in `Gamma`, which is a
reduction modulo `m` (one `c_mor`).

## Orbit splitting across markings

`split_orbit(marked, H)` (`lattice_morphisms.py:4638`) has three branches.

- **`H = O(L)`.** One representative; cost zero.
- **Character-kernel branch** (`:4669-4707`). `H` contains the kernel of the
  finite character map `chi = (rho, det, sn)` to a quotient `Q` of order at
  most `2 * 2 * |O(q_L)|`. For a primitive isotropic line in the arrow
  category, the Eichler certificate (`isotropic_line_stabilizer_image`) gives
  the image of the line stabilizer without parabolic generators. Otherwise
  the stabilizer generators of the marking are computed and the double cosets
  `chi(H) \ Q / chi(Stab)` are formed. Cost: stabilizer generators plus
  `c_ss(|Q|)` plus one `c_mor` per representative.
- **Finite preimage branch** (`:4708-4727`). `H = rho^{-1}(S)` for
  `rho: O(L) ->> F`. The `H`-orbits in `G/P` are the double cosets
  `S \ F / rho(P)`. Their number is at most `|F| / |rho(P)|` and equals it
  for `S = 1`. Cost: stabilizer generators of the marking, `c_mat` for
  `rho(P)`, a double coset enumeration of `|F| / |rho(P)|` cosets (GAP
  `DoubleCosets`, `O(|F|)` group operations), and one section lift and one
  apex action per double coset: `2 * c_mor` each.

The counts for `2U` modulo 3 follow from the order 1152 of `O^+_4(3)` and the
stabilizer orders 36, 72, 144, 36 and 16 of a vector, line, plane, flag and
hyperbolic summand. They give 32, 16, 8, 32 and 72 orbits. A based pair with
prescribed Gram `U` has stabilizer image of order 4: 288 orbits. Complete
rank-three flags in `3U` modulo 2 have Borel image of order 64 in
`O^+_6(2)`: 630 orbits. These counts are proved in the docstring of
`witnesses/test_congruence_orbit_splitting_across_markings.sage`.

Worst case: the finite preimage branch is linear in `|F| / |rho(P)|`, and
`|F|` grows like `m^{dim O(n)}`. For `n = 6` and `m = 3`, `|O^+_6(3)| > 10^{8}`.
Splitting is not attempted through an enumeration of `F`; only the double
coset count enters the per-representative term.

## Marking classes up to `O(L)`

`vector_classes`, `isotropic_classes` and `flag_classes` for a lattice that
contains `U + U` use the Eichler criterion (Gritsenko–Hulek–Sankaran,
*Abelian surfaces with (1,2)-polarization*, Lemma 3.5): the orbit of a
primitive vector is determined by its square and `v / div(v)` in `A_L`. The
number of classes is at most the number of elements of `A_L` with the right
`q` value, modulo `O(q_L)`. Cost: `|A_L|` evaluations of `q` plus one
`c_tqf`. Planes and flags add one parabolic stabilizer each, owned by
`boundary-components.md`.

## Non-termination and refusal

No operation in this file enumerates an infinite set. A refusal arises in the
callee that computes `O(L)` generators: positive index at least 2, rational
Witt index below 2 and no recognized rank-four model raise
`UnfinishedCapability` (`src/sage_indefinite_port/indefinite/recursive.py:1175-1182`,
source `Indefinite.jl@374a5ebb` `IndefiniteForms.g:830`). That refusal costs
one Witt-index computation: Hasse invariants at the primes dividing
`2 det L`, a bounded number of `c_lat`.
