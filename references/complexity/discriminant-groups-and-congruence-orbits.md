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

Every unit is a measured primitive in `ci/test_budget_models/unit_costs.json`,
named by its key, or an operation count times such primitives.

| Symbol | Primitive | Unit cost | Basis |
| --- | --- | --- | --- |
| `c_add` | Addition of two lattice elements | `lattice_element_add` | measured |
| `c_pair` | One pairing `b(x, y)` of lattice elements | `lattice_element_pairing` | measured |
| `c_mor` | One morphism evaluation on rank `n <= 12` | `12 * lattice_element_add` | `n` element operations |
| `c_lat` | Construction of a formed lattice of rank `n <= 12` | `16 * lattice_construction_rank_3` | the work scales with the `n^2` Gram entries: `(12/3)^2` rank-3 constructions |
| `c_ss(d)` | GAP Schreier–Sims order of a permutation group of degree `d` with `k` generators | `libgap_call + k * d * log2(d)^2 * gap_permutation_image_527` | Seress, *Permutation Group Algorithms*, Sec. 4.5 and Thm. 4.5.5: nearly linear `O(k d log^c |G|)` permutation images |
| `c_mat(p, n)` | Order of a matrix group in `GL_n(F_p)` | `c_ss(p^n)` for `p^n <= 10^3`; otherwise `t_reduction` of `orthogonal-group-generators.md` §Finite images | for large `p^n` the vector action has `(p^n - 1)/(p - 1)` points |
| `c_tqf(A)` | `O(q_L)` of a finite quadratic form | `t_oq(|A_L|, log2 |A_L|)` of `orthogonal-group-generators.md` | permutation degree `|A_L|` |
| `C_GAPZ` | Order or finiteness of a small matrix group in GAP | `gap_order_matrix_group_so_plus_10_2` | measured |
| `C_MM` | One product of dense matrices of size at most 8 | `zz_matrix_product_8` | measured |
| `C_OP` | One equivariant or normalizer step after the centralizer | `100 * zz_matrix_product_8 + 10 * qq_rank_8` | at most 100 products and 10 eigenspace ranks |
| `C_WITT` | Rational Witt index | `4 * qq_rank_12` | Hasse invariants at the at most 4 primes involved, one rational diagonalization each |
| `PURE` | Pure-Python check | `100 * libgap_call` | at most 100 interpreter-level calls |
| fixture parse | `tests/fixtures/allcock_rank3_reflective.json` | `8595 * 20 * zz_dot_12` | 8595 cases of about 20 values |

A budget is the derived maximum work of its branch, with these unit costs
and no multiplicative factor (`ci/test_budget_models/F.py`).

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

## Markings on degenerate lattices

`L = C + R` with `R` the radical of rank `t` and `C` nondegenerate. `O(L)`
consists of the block matrices `[A B; 0 D]` with `A in O(C)`,
`B in Hom(C, R)` and `D in GL_t(Z)`
(`src/sage_indefinite_port/indefinite/marked_lattices.py`, `_degenerate_pointed_classes`
and `_radical_residues`; source `IndefiniteForms.g:1496-1511`, `ExtendIsometryGroup`).
The generators are those of `O(C)` (branch of `orthogonal-group-generators.md`
for `C`), `t * rank(C)` elementary shears and the generators of `GL_t(Z)`:
two for `t >= 2`, one sign for `t = 1`.

A primitive vector `v = (c, s)` with `c` of content `d` has orbit invariant
(class of `c` under `O(C)`, class of `s` in `(Z/d)^t` under the image of
`GL_t(Z)`). For `t = 1` that image is `{+-1}`, so the residue classes are
`(Z/d)^* / {+-1}`, of size `phi(d)/2`. For `t >= 2`, `SL_t(Z) -> SL_t(Z/d)` is
onto and transitive on unimodular residues: one class. A transporter is one
congruence lift of a matrix in `GL_t(Z/d)` (strong approximation for `SL_t`,
Newman, *Integral Matrices*, Thm. VII.21), computed by a Smith form of size
`t`: `O(t^3 log d)` integer operations.

Cost: one generation of `O(C)` plus `O(t * rank(C) + t^2)` morphism
operations per stabilizer or transporter. A lazy enumeration of classes of
square `n` with content bound `h` visits the divisors `d <= h` and, for each,
`phi(d)/2` residues: `O(h^2)` candidates, each a canonicalization of constant
cost. The quotient for `n = 0` is infinite (one class per content), so only a
prefix is ever computed.

## Markings by isotropic induction

`primitive_marking_classes`, `basis_marking_classes`, `flag_invariant` and the
transporters between isotropic markings use the induction of
`IsotropicInduction` (`isotropic_induction.py`; source
`IndefiniteForms.g:1700-1760`, `INDEF_FORM_GetOrbit_IsotropicKplane`). A
rank-`k` isotropic marking is built from a rank-`(k-1)` marking `I` by the
classes of isotropic lines in `I^perp / I` under the image of `Stab(I)`.
The number of classes at each step is at most the number of isotropic
classes of `A_{I^perp/I}` modulo the Levi image, bounded by `|A_L|`
(Nikulin 1980, Prop. 1.15.1, for the gluing of `I^perp / I`). For unimodular
`L` every step has one class, so a rank-`k` flag in `mU` has one `O(L)`-class
and the induction costs `k` parabolic stabilizers.

Cost: `sum_{j<k} (stabilizer of a rank-j marking + classes of lines in a
rank-(n-2j) quotient)`. The stabilizer of a rank-`j` marking is the Levi
`GL_j(Z) x O(I^perp/I)` and the unipotent radical of rank
`j(n-2j) + j(j-1)/2`. When `k` exceeds the Witt index the induction stops at
the first step with no isotropic line, after one Witt-index computation.

The splitting of a marking class under `SO(L)` is the character-kernel branch
of `split_orbit` with `|Q| = 2`.

## Pointwise stabilizer of a null vector and its shell decision

For `e` primitive isotropic in `L = U + U`, the pointwise stabilizer
`Gamma = Stab(e)` is `O(e^perp/e) x` the translations `e^perp/e`
(Eichler transvections). `O(e^perp / e) = O(U)` has order 4. A vector `v` in
the isotropic cone has `Gamma`-invariant `(content, b(e, v), class)`; the
decision of `Gamma`-equivalence is one invariant evaluation and one
composite witness `u(rho, gam) . S` with `S` among the four symmetries.

Cost per decision: `O(n^2)` element operations for the invariant plus at
most four morphism compositions for the witness. A test over a set of `N`
vectors that decides every pair costs `N(N-1)/2` decisions. The box
`[-2, 2]^4` contains 129 isotropic vectors of `U + U` (direct count of
`am + cd = 0`), so 8256 pairs.

## Normalizers of subgroups of `O(L)`

`normalizer_generators` (`finite_normalizers.py:1318-1343`) dispatches on the
generated group `H`:

| Branch | Condition | Cost driver |
| --- | --- | --- |
| Finite | every generator of finite order and `H` finite | `C(H)` from the commutant (owned by `centralizers.md`), plus one twisted transporter per automorphism of `H` that preserves the characters: at most `|Aut(H)|` |
| Infinite cyclic | one generator of infinite order | Unit group of the commutant order `Z[f]`: one regulator computation of a real quadratic order, plus one inversion transporter |
| Abelian | commuting generators | Torus algebra (rank of the Lie algebra of the Zariski closure), decomposition normalizer of that algebra |
| Arithmetic | the congruence kernel modulo 3 is not abelian | Image of the supplied group in `GL_n(F_3)`, closure of size at most `|O_n(3)|`, and Schreier generators of the preimage |
| Congruence kernel | kernel abelian, not semisimple | As above, with the unipotent kernel supplied |
| Virtually abelian | otherwise | Torus normalizer and the finite quotient |

The decomposition normalizer of an orthogonal splitting
`L = L_1 + ... + L_r` is the product of the `O(L_i)` and the permutations of
isometric summands: order `prod |O(L_i)| * prod (multiplicities)!` for
definite summands.

## Equivariant markings

`EquivariantMarking(L, X, action, F)` computes the stabilizer and transporters
of a matrix `X` under the centralizer of a finite family `F` of isometries.
The centralizer is computed by `centralizers.md`. The stabilizer is the
subgroup of the centralizer that preserves `X` and its glue in
`A_K` for `K` the orthogonal complement of the image of `X`. Its cost is one
centralizer, one glue stabilizer of degree `|A_K|`, and for an infinite
centralizer one unit-group computation of the commutant order.

## Definite lattices: route of the generators and of the isometry test

The rule follows from the code of sage-extensions:

- `O(L)` generators of a definite lattice over `ZZ`:
  `_computed_group_generators` (`categories/lattice_morphisms.py:4447-4448`) calls
  `_engine_group` (`:4278-4330`). Rank one uses the roots of `a u^2 = a`
  (`:4294-4311`). Every higher rank uses PARI `qfauto` through
  `QuadraticForm.automorphism_group` (`:4313-4330`). No Julia or Oscar call is made.
  `RECURSION.orthogonal_group_generators` (`recursive.py:1135-1148`) reaches the
  same engine through `lattice.O().group_generators()`.
- Isometry of two definite lattices over `ZZ`: `_witness_matrix`
  (`lattice_morphisms.py:5271-5330`) uses PARI `qflllgram` and `qfisom`
  (`QuadraticForm.is_globally_equivalent_to`). The Oscar bridge
  (`_oscar_lattices.definite_isometry`, `:5327-5331`) is reached only when a base
  ring is not `ZZ`.

Every definite lattice in this file's scope is over `ZZ`, so no budget here
includes the Julia start-up cost. The definite cost is B's
`t_definite(n) = 0.5 + 0.02 n^2` (`orthogonal-group-generators.md`, Definite kernel),
without the bridge.

## Derived counts used by the budgets

- **Shell decision.** An invariant `(content, b(e, v), class)` reads 4
  coordinates, forms the 16 products of `b(e, v)` and takes 4 gcd steps: 24
  element operations. A decision evaluates two invariants. A witness is at most
  four compositions of a morphism plus a check on the 4 basis vectors, which is
  64 element operations. The box `[-2, 2]^4` has 129 isotropic vectors. At most
  1150 pairs share `(content, b(e, v))`, and so at most 1150 pairs share a
  class.
- **Canonicalization.** The canonical form of a vector of height `h` on
  `<2> + 0` scans the members of the same square with coordinates of absolute
  value at most `h`: `x = +-h`, `y` in `[-h, h]`, `2(2h + 1)` candidates, each
  one transporter.
- **Radical congruence image.** `radical_congruence_image(L, m)` is a subgroup of
  `GL_t(Z/m)` acting on the `m^t` residues: one Schreier–Sims of degree `m^t`.
- **Normalizers.** The finite branch costs the centralizer
  (`centralizers.md`, Cost) plus one twisted transporter per automorphism of
  `H` that the test exercises: 6 for a Klein four-group (`|Aut| = |S_3|`), 2 for
  `Z/6`. Each transporter is an isometry of a rank-2 commutant piece (B's
  binary model). The congruence branch closes the image modulo 3 inside
  `O^+_4(3)`, of order at most 1152. On `3U` the closure of the supplied group
  (one swap, two signs, two commuting transvections of order 3) has at most
  `2^3 * 9 * 8 = 576` elements.

## Parabolic generators

The generators of the stabilizer of a primitive isotropic sublattice `I` of
rank `j` in a lattice of rank `n` are the Levi generators of `GL_j(Z)` (two) and
the generators of the unipotent radical (at most `j(n - 2j) + j(j - 1)/2 <= n^2/4`).
Each is an isometry given by `n` images. Cost: one construction of the
reduction lattice `I^perp / I` (`c_lat`) plus `(n + n^2/4) * c_mor`. The
generators of `O(I^perp / I)` are a separate term from
`orthogonal-group-generators.md`.
