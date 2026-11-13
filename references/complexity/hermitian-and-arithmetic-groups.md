# Complexity of Hermitian lattice groups, commutant centralizers and arithmetic-group generators

These routines serve issues #51, #54, #56 and #61. The pinned trees under
`references/vendor/` contain no routine for them. The authority for each
route is the source cited in its module docstring. Costs are counts derived
from the mathematics. Unit costs of backend primitives are in the last
section. The Clifford, quaternion norm lattice, number-field genus and
semisimple Voronoi routes are in `clifford-and-hermitian-lattices.md`.

## Explicit unitary groups

Sources: `act_eisenstein.py`, `euclidean_picard.py`, `gaussian_picard.py`,
`binary_hyperbolic_unitary.py`, `eisenstein_hyperbolic_unitary.py`.

The generators are the printed matrices of the source:

- Allcock-Carlson-Toledo (J. Alg. Geom. 11 (2002), Sec. 7) for
  `E^{4,1}` over `Z[omega]`.
- Zhao (Geom. Dedicata 2012) for the Picard groups `PU(2,1; O_d)` with
  `d = 2, 7, 11`.
- Falbel-Francsics-Lax-Parker (arXiv:0911.1104) for `SU(2,1; Z[i])`.

**Branch condition.** Word reduction (`word`, `reduce_null_vector`,
`null_transporter`) is a descent on the height of the null vector `g e`.
Each step applies a cusp translation chosen by nearest-integer division in
`O_d`, then the inversion. The descent terminates when `O_d` is
norm-Euclidean. The norm-Euclidean imaginary quadratic orders are `d = 1, 2,
3, 7, 11` (Lemmermeyer, Expo. Math. 13 (1995), Sec. 4.1). For any other `d`
the step raises.

**Count.** The covering radius `rho` of `O_d` in `C` is less than 1, so each
step multiplies the height by at most `rho^2 < 1`. The number of steps is
`O(log height)`. For a product of `k` generators the height is
`exp(O(k))`, so the word length is `O(k)`.

**Spinor character.** For `act_eisenstein.spinor_character` on the trace
lattice of rank 10, Cartan-Dieudonne factorization uses at most 10
reflections.

## Congruence images

Sources: `picard_congruence.py`, `hermitian_parabolic_congruence`, and
`OrderHermitianModule.stabilizer_of_module`.

**Method.** The index of a congruence subgroup equals the order of the finite
image (Schreier's lemma; Seress, *Permutation Group Algorithms*, Lemma 4.2.1).
The routine enumerates the image by breadth-first search. With `s`
generators it returns `[G : H](s - 1) + 1` Schreier generators.

**Count.** The image orders are invariants of the generator sets. They were
computed in GAP from the `Z`-coordinate matrices modulo 2
(`.tmp/probe/estimates/inv_images.py`, `inv_nonmax.py`):

| Group | Image modulo 2 |
| --- | --- |
| FFLP generators of `SU(2,1; Z[i])` | order 48 (`<M>`: order 2) |
| Zhao generators, `d = 7` | order 168 |
| `U(1,1; Z[i])` moves on the `Z[2i]^2` frame | order 12 |

The setwise stabilizer of the Gaussian line is
`P = l(a, u) n(Z, T)` with `a, u` in `mu_4`, `Z` in `(1 + i) Z[i]` and
`T` in `iZ` (hand derivation in the witness docstring). Its image modulo 2
has order at most 64. This bound is not computed.

**Cost.** (image order) x (generators) x (one product of `3 x 3` matrices
over the field), plus one membership test per Schreier generator.

## Hermitian parabolics

Sources: `HermitianWittFrame` (sage-extensions) and the port's parabolic
witnesses.

**Structure.** The stabilizer of a Witt flag is `L N`. `N` is the
Heisenberg unipotent group `n(Z, T)`. `L` is the Levi group: pairs
`(A, u)` with `A` in `GL_k(O)` and `u` in `U(L_0)` on the anisotropic
kernel.

**Branch condition.** A Levi pair lifts exactly when it preserves the gluing
of the frame. The gluing is the discriminant of the frame modulo the
conductor. The lift test solves one linear system over `O`.

**Count.** The Levi image has index (number of gluing classes) x (unit
condition):
- split Gaussian, rank 2: `|P^1(F_2)| = 3`;
- glued Eisenstein, rank 2: `|P^1(F_3)| * 2 = 8`;
- glued Eisenstein line: `36 / 18 = 2`.

**Trace** (`sage_extensions/categories/hermitian_lattices.py`):

- `levi_lift` (l.727) calls `_lift_to` (l.731), which calls `_admissible`
  (l.~395). `_admissible` solves for the parameter lattice with `_solutions`
  (Hermite forms over `Z`). It evaluates the quadratic entry map at
  `1 + 2s + s(s-1)/2` sample points. Each point is one matrix of size `n x n`
  over the field. Here `s = 2mk + k^2` is the parameter count. It then loops
  over `modulus^{|section|}` residues. Per-lift work is therefore `O(s^2 n^3)`
  field operations plus the residue loop. The modulus is
  `2 lcm(denominators)`, which grows with the glue (`theta` in the
  Eisenstein cases).
- `_levi_cosets` (l.693-717): the moves are the generators of `GL_k(O)` (or
  of its Borel for a flag) and `U(I^perp/I)`. For each representative and each
  move, the routine runs lift tests against known representatives until one
  succeeds. The expected count is `index * |moves| * (index + 1) / 2`, plus
  one lift per Schreier generator in the tests.
- `quotient_isometries_to` (l.524) uses the automorphisms of the definite
  trace form of the quotient through `QuadraticForm.automorphisms`.
- `stabilizer_generators` (l.592) returns the kernel generators from
  `integral_unipotent_kernel` and one lift per Levi generator.

**Candidates.** When the frame is unimodular the candidate set is the whole
finite group. For `E^{4,1}` this is `mu_6 x (mu_6 wr S_3)`, of order 7776.
Each candidate gets one lift test.

## Commutant centralizers

Source: `commutant_groups.py:1-60`, which gives the dispatch. The branch is
selected by the Wedderburn type of the rational commutant `A` and by its
adjoint involution:

| Branch | Selecting condition | Dominant work |
| --- | --- | --- |
| radical | `A / rad = Q` | integral point in the coset `t_r (+-N(Q))`: one Hermite form |
| quaternion, symplectic | `A` indefinite division quaternion algebra over `Q`, canonical involution | Voronoi unit generators of `O` (section Quaternion Voronoi), then norm-one Schreier extraction |
| quaternion, orthogonal | involution `x -> c^{-1} xbar c` | norm-one torus of `Q(c)`: one `bnfinit` of degree 2 and its fundamental unit |
| Levi | `A / rad = prod M_{n_k}(Q)` | Morita reduction to form isometries or GL pairs on the Levi modules. Multiplier-ring and index invariants `[O_K X : X]` decide the absent cases. |
| central field | `A` central simple over `K != Q` | the unitary model: hyperbolic, ocean, `_BinaryHermitianOrder`, symmetric rank 2 (norm equations), symmetric rank `>= 3` (Eichler), Witt index one, or alternating |
| split simple | `A = M_n(Q)` | `O(b)` or `Sp(b)` of the recovered integral form |
| etale | `A = prod K_j` | relative norm equations (`bnfisnorm`), the regulator and units of each `K_j` |

The transporter cost is (model isometry) + (number of coset representatives
of `C(Z)` in `U(Lambda, h)`) x (one integrality test). The coset count is the
index of the order `O` in its maximal order, measured on unitary units.

**Symmetric rank `>= 3` over a totally real `K` with `h(K) = 1`.** The route
is `eichler_orthogonal`. The lattice splits as `H + H + L_0`. `O(Lambda)` is
generated by Eichler transvections and spinor-norm representatives
(Kneser, Arch. Math. 7 (1956); Hsia, J. Number Theory 1975). Spinor norms at
a ramified dyadic prime come from the Earnest-Hsia tables and need one
realizing reflection each.

## Quaternion Voronoi

Sources: `quaternion_voronoi.py` and `quaternion_voronoi_number_field.py`.
The method is BCNS (2015) Sec. 2-4 and Opgenorth (2001).

**Indefinite algebra over `Q`.** For `(-1, 3)_Q` and `Z<1, i, j, k>` the cone
is `Sym_2(R)` of dimension 3. There is one perfect class (the test asserts
it). The number of facets is at most `|Min|`.

**Covolume.** Shimizu's formula gives the covolume of the maximal order of
discriminant 6 as `(pi/3) (2-1)(3-1) = 2 pi / 3` (Voight, *Quaternion
Algebras*, Thm. 39.4.13). The index-2 order has local unit index 3 at 2, so
the covolume of `Lambda^1` is `2 pi`.

**Totally definite branch.** `norm_generator` enumerates
`nrd(x) <= nrd(I)` with one `qfminim` in rank 4 (Kirschmer-Voight, SIAM J.
Comput. 39 (2010), Sec. 6).

**Over `Q(sqrt 2)`.** Two real places give a cone of dimension 6. The number
of perfect classes is not derived here. The budget assumes at most 10 classes
and at most 30 facets, on the scale of the `Q_{2,3}` example in Nebe's 2014
slides, which has 3 classes. Computing the count is the computation under
test.
Sources searched without a count for an order over a real quadratic field:
- BCNS (2015) Sec. 6;
- Nebe's slides (2014);
- the RWTH tables "Unit groups of orders", which list quaternion algebras over
  imaginary quadratic fields;
- Braun-Coulangeon, Math. Comp. 2015, on imaginary quadratic fields;
- Scheckelhoff-Thalagoda-Yasaki (arXiv:2105.00593).

The same assumption applies to `(11, 19)_K` over `Q(sqrt 5)`.

## Bestvina-Savin ocean

Source: `hermitian_ocean.py`; Bestvina-Savin (arXiv:1104.1474), Thm. 5.1
and Thm. 6.1.

**Input condition.** The input is an anisotropic indefinite binary Hermitian
form over `Z[i]`. An isotropic, definite or non-integral form raises.

**Count.** The number of cell orbits is finite. It is determined by the
orbifold Euler characteristic:
- `Delta = 6`: one square cell; vertex stabilizers of orders 4, 3 and 4;
  `chi = -1/6` (Example 7.1).
- `diag(1, -3)`: `chi = -1/2`, by the mass formula through the quaternion
  route.

The ocean is built by enumerating vectors of `Z[i]^2` with
`|f(v)| <= Delta` up to the cell orbits.

## Fuchsian domains

Source: `fuchsian_domain.py`; Voight, *Computing fundamental domains for
Fuchsian groups*, J. Theor. Nombres Bordeaux 21 (2009).

**Count.** For discriminant 6 the area is `2 pi / 3`. The signature is
`(0; 2, 2, 3, 3)`. The exterior domain has `2(2g + e) = 8` paired sides up to
elliptic fixed points. The enumeration stops when the hyperbolic area of the
domain equals the Shimizu area.

## Bianchi groups

Source: `bianchi_groups.py`; Swan, Adv. Math. 6 (1971).

**Count.** For `Q(sqrt -5)`, `h = 2`. The singular points modulo `O` are one
class, `(1 + w)/2`, for the nontrivial ideal class `p = (2, 1 + w)`. Swan's
generators come from the hemispheres over the fundamental rectangle of `O`.
The number of hemispheres is bounded by the elements of norm at most
`|d|/4 + 1` up to translation.

## Number-field orthogonal groups

Sources: `ternary_orthogonal_number_field.py` and
`quaternary_orthogonal_number_field.py`.

**Route.** `SO(M)` of a ternary form is `O^1 / +-1` for the maximal order of
`C_0(M) = B`. In the quaternary case with square discriminant,
`C_0(M) = B x B`. For `B = (11, 19)_K`, `K = Q(sqrt 5)`, `B` ramifies at the
two primes over 11 and at no real place. Its unit group is computed by
Voronoi over two real places, with cone dimension 6, so the cost of the
Voronoi section applies. The number of perfect classes is not derived here.

## Li-Margulis boxes

Sources: `li_margulis.py`; Li and Margulis, *Effective estimates on integral
quadratic forms*.

**Box enumeration.** A box of size `b` on rank `r` holds `(2b + 1)^r`
vectors. The routine filters them by the prescribed norms and pairings. Each
new candidate is deduplicated with a marked transporter, which is the
indefinite route of `vector-orbits.md` (owner B).

**Trace.** `embedding_orbits_in_box` (`li_margulis.py:112`) and
`isotropic_vector_orbits_in_box` (l.~640) deduplicate every candidate with
`pairing_slices.marked_transporter` (l.68). That function returns at once
when the Gram matrices differ. For a mixed marking it uses
`_mixed_marking_transporter`. Otherwise it calls `isometry` in the marking
category, which is the route priced in `isometry-testing.md` (owner B). The
budgets take 0.5 seconds per call on rank at most 4. This is an assumption
until that file prices these inputs.

Hand counts of the candidates:
- `U + A_1`, box 2: 10 vectors of square -2.
- `U + A_1`, box 1: 8 primitive null vectors.
- `U + U(2)`: `isotropic_pair_grams` gives 3 Gram matrices: `d = 1, a = 0`
  and `d = 2, a = 0, 2`.

**Constants.** These are closed rational formulas. The largest integer is
`(4 * 27 * 216)^{18}` (about 100 digits).

## Unit costs assumed

| Primitive | Cost |
| --- | --- |
| Quadratic number-field scalar operation | 2 microseconds |
| `3 x 3` matrix product over a quadratic field | 0.1 milliseconds |
| `5 x 5` matrix product over a quadratic field | 0.4 milliseconds |
| One `levi_lift` on rank 3 | 0.02 seconds |
| One `levi_lift` on rank 5, split or unimodular | 0.05 to 0.1 seconds |
| One `levi_lift` on rank 5, glued with denominator 3 | 0.3 seconds |
| PARI `bnfinit` or `bnfisnorm`, degree 2 | 10 milliseconds |
| PARI `qfminim`, rank at most 9 | 10 milliseconds |
| Quaternion Voronoi over `Q`, one perfect class | 0.5 seconds |
| Voronoi over a real quadratic field, cone of dimension 6 | 5 seconds (assumed class count) |
| Marked transporter, rank at most 4, indefinite | 0.3 to 0.8 seconds |

## Scope G pricing rule

The budgets in `.tmp/probe/estimates/G.py` are the derived maximum work of the branch.
A pure computation on one model instance is charged once. Example: one
Voronoi tessellation per order. The code sometimes repeats such a
computation. `QuaternionVoronoi.perfect_forms`, `unit_generators` and
`reduce_unit` each call `_enumerate` again. `QuaternionOrderVoronoi.transporter`
recomputes the minimal vectors of each known form for each neighbour. A budget
does not include this repeated work, so a test that repeats it can exceed its
budget.

## Voronoi over Q: enumeration count

Source: `quaternion_voronoi.py`, order `Z<1,i,j,k>` of `(-1, 3)`.

**Start form.** For `x = t + x i + y j + z k`, the matrix is
`[[t + y s, -x + z s], [x + z s, t - y s]]` with `s = sqrt 3`. So
`phi_I(x) = Tr(x^t x) = 2t^2 + 2x^2 + 6y^2 + 6z^2`. Its minimal vectors are
`+-1, +-i`, and both have the ray `(1, 0, 1)`. The ray rank is 1. Each flip
adds at least one to the rank, so at most 2 flips reach a perfect form.

**Unit costs.** One `2 x 2` product over `Q(sqrt 3)` is about 0.35 ms. One
`_gram` is 16 pairings, about 11 ms. One `_short` is a Gram plus one PARI
`qfminim` of rank 4, about 15 ms. One `_minimal_vectors` is about 50 ms. One
`_flip` is 1 to 3 bound doublings, about 0.3 s.

**Count.** The perfecting flips cost about 0.6 s. Each facet costs one flip,
one minimum and one `_equivalences` (2 m candidates, m minimal vectors), about
0.5 s. The test asserts one perfect class. The number of facets is at most the
number of distinct rays; the budget assumes at most 6. Then
`V_Q = 0.6 + 6 x 0.5 = 3.6 s`. The `(-3, 2)`-type order of
`_BinaryHermitianOrder` (Eisenstein plane) is assumed to have the same size.

**Box test.** In the box `[-5, 5]^4` there are 308 elements with
`|t^2 + x^2 - 3y^2 - 3z^2| = 1` (stdlib count). Each `reduce_unit` descent
takes about 4 steps of 12 moves at about 0.6 ms, so about 30 ms. The word ball
of radius 4 on at most 12 moves has at most 1.8e4 elements. The test computes
each product twice, so about 3.8e4 quaternion products.

**Public route.** `L.O().group_generators()` for `diag(2, 2, -6, -6)` reaches
`recursive._witt_zero_generators`. It calls `li_margulis.is_isometric_in_box`
with radius 3 (`7^4 = 2401` vectors) and then `orthogonal_group_generators`.

## Ocean tests: cell counts

Source: `hermitian_ocean.py`. A cell is a pair of lax vectors, so
`_transporters` tests at most `2 x 2 x 4 = 16` candidates. Each candidate
costs about 0.5 ms (`is_unitary` on a `2 x 2` matrix over `Q(i)`). The Delta =
6 example has 1 cell orbit (Example 7.1). For `diag(1, -3)`,
`chi = -1/2` and the square cells give at most 4 cell orbits. Then
`_quotient` costs at most `4 x 4 x 4 x 8 ms = 0.5 s`. The image tests do a
breadth-first search in `GL_2(F_p)`: 480 elements for `p = 5` and at most
`168 x 156 = 26208` for `p = 13`.

## Fuchsian domain: enumeration radius

Source: `fuchsian_domain.py`, maximal order of discriminant 6. The pass of
radius `r` enumerates the vectors with `N(x) <= 2 cosh r`. The count of
`gamma` in `O^1` with `d(0, gamma 0) <= r` is about
`2 pi (cosh r - 1) / (2 pi / 3) = 3 (cosh r - 1)`. The count of all lattice
vectors in the rank-4 ellipsoid is about `(pi^2 / 2) B^2 / sqrt(det)`. The
budget takes `det = 36` (the discriminant of the trace form), which gives
about `0.4 B^2` vectors. Each pass enumerates again from zero. Assumption:
the Dirichlet domain has radius at most 2.5, so the passes stop at `r <= 5`.
Then there are about 1.1e4 vectors at 0.3 ms and about 200 norm-one reductions
at 5 ms. Work is about 5 s.

## Swan covering for Q(sqrt -5)

Source: `bianchi_groups.py`. The singular point `s = (1 + w)/2` lies on
hemispheres with `N(mu) = 4` (`mu = 2`) and `N(mu) = 6` (`mu = 1 + w`). So the
punctured squares need `b >= 6`. For `b < 4` no hemisphere passes through `s`,
and `_punctured_radius` returns at once. The `lambda` loop of
`_hemisphere_representatives` makes `sum N(mu)^2` ideal tests over the `mu`
with `N(mu) <= b` (stdlib count):

| b | elements `mu` | ideal tests |
| --- | --- | --- |
| 6 | 10 | 228 |
| 12 | 16 | 714 |
| 20 | 24 | 2810 |
| 50 | 72 | 66376 |

Assumption: the floor-height bound gives `1/h^2 <= 12`. The hemispheres of
`N(mu) = 6` have radius `1/sqrt 6`. The number of generators is then about
60. The box covering is assumed to be about 5 s. For `O + p`, `alpha = 2` and
`O/2O = F_2[e]/(e^2)`, so the orbit index is at most 6.

## Voronoi over Q(sqrt 5) for (11, 19): perfect-class estimate

Source: `ternary_orthogonal_number_field.py`,
`quaternary_orthogonal_number_field.py`, `quaternion_voronoi_number_field.py`.
`B = (11, 19)_K` is unramified at both real places, so
`quaternion_unit_engine` returns `QuaternionOrderVoronoi` with
`dim Herm = 6`.

**Class count.** Gunnells-Yasaki (arXiv:0908.3130, Table 1) list 2
`GL_2(O)`-classes of perfect binary forms over `Q(sqrt 5)` (Ong 1986). The
orbifold Euler characteristic of `O^1 / +-1` for a maximal order is
`Phi(D) = (11 - 1)(11 - 1) = 100` times that of `PSL_2(O_K)` (Shimizu). The
Ong forms have stabilizers of order about 24 to 60. Then the mass of `B` is
about `100 x (1/24 + 1/60) = 6`. With stabilizers of order 2 to 4 this gives
`P = 20` classes (assumption). The budget takes `f = 15` facets per class.

**Unit costs.** One `gram` has 64 entries, each two traces of products of
three quaternions over `K`, about 50 ms. One `_step` is about 8 Grams, so
0.4 s. With cached minima, one neighbour costs about 0.6 s.

**Count.** Ternary: `P f = 300` neighbours x 0.6 s = 180 s. The normalizer
breadth-first search (orbit length assumed at most 20, about 170 generators,
5 ms per `act`) costs about 17 s. The Schreier conjugations cost about 17 s.
Work is about 220 s.

Quaternary: two engines (`R_1`, `R_2`). The `_conjugator` engine is the
engine of `R_1`. Then there are the pair Schreier generators (index at most
4), the lattice-orbit breadth-first search on about 1400 actions, and about
2.8e4 stabilizer checks. Work is about 600 s.

## Commutant tests: branch classification

The dispatch is `commutant_groups.centralizer_generators` and
`equivariant_transporter` (l.2111-2180). The `IsometryCommutant` of rank `r`
solves an integral kernel with `r^2` unknowns. It costs about 0.1 s for
`r <= 4`, 0.5 s for `r = 8` and 2 to 3 s for `r = 16` (dimension 32 gives
1024 products of `16 x 16` matrices).

Branch corrections against the module-docstring classification:

- The real quadratic orthogonal plane is `CentralFieldCommutant` (symmetric,
  rank 2, torus of a degree-4 field). It is not `EtaleCommutant`.
- The swapped-factor test has a commutative commutant of dimension 4, so it
  takes `EtaleCommutant`. It does not take `levi_transporter`.
- The anisotropic skew-Hermitian plane takes `clifford_commutant`. Its
  Clifford algebra `C` is over `Z = Q(z)`, `z^2 = 3/2`, which is `Q(sqrt 6)`.
  `C` ramifies at the primes over 2 and 3 and at both real places
  (`.tmp/probe/estimates/inv_clifford_skew.py`). So `C` is totally definite,
  `Herm` is the centre, and the cone has dimension 2. The route has no
  indefinite Voronoi.
- The Lipschitz order over `Q(sqrt 2)` is totally definite. Its cone also
  has dimension 2.
- The symplectic and Hermitian quaternion branches, and the
  `_BinaryHermitianOrder` branch, cost one `V_Q` each.

**Witness transporters** (`witnesses/test_commutant_groups.sage`). Each test
makes two `equivariant_transporter` calls. The `3B` call stops at the
determinant or norm ratio.

- Levi tests (`r = 6`, `A / rad = M_2(Q)`): `levi_transporter` runs
  isometries of rank-2 Levi forms and lifts from the radical. The nonmaximal
  Levi order is the Eichler order of level 2, so it adds
  `[GL_2(Z) : Gamma_0(2)] = 3` coset lifts.
- Hermitian quaternion and isotropic skew-Hermitian tests: one `V_Q`, charged
  once. Each call adds a reduced-norm point search.
- Definite Hermitian test: the units come from a Fincke-Pohst enumeration.
- Golden test (`r = 16`): the commutant costs 3 s. Each call adds `H + H`
  models of both forms and a coset search over about 32 Eichler generators.
