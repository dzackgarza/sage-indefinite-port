# Complexity of the Clifford, Hermitian, quaternion, number-field genus, semisimple Voronoi and Witt-index-one Hilbert routes

These routines serve issues #51 and #54. The pinned trees under
`references/vendor/` contain no routine for them (no match for Clifford,
Hermitian or Voronoi in `Indefinite.jl@374a5ebb` or
`polyhedral_common@1592b246`). The authority for each route is the
literature cited in the module docstring. Costs below are counts derived from
the mathematics. Unit costs of backend primitives are in the last section.

## Even Clifford centre

Source: `src/sage_indefinite_port/indefinite/clifford_even_rank.py:19-98`.
Builds `CliffordAlgebra(Q)` of dimension `2^n` and its centre `Z(C_0)`, which
is `Q[z]` with `z^2 = (-1)^{m(2m-1)} prod q(e_i)` for `n = 2m`
(Knus, *Quadratic and Hermitian Forms over Rings*, IV.2; Lam, *Introduction to
Quadratic Forms over Fields*, V.2.4).

Branches, selected by the square class of `z^2`:

- **Split** (`z^2` a square): idempotents `(1 +- z/sqrt(z^2))/2` and the
  dimension of each component. Cost: `2 * 2^{n-1}` products of a
  `2^{n-1}`-term element with a basis monomial, so `O(4^n)` monomial
  products. For `n = 6`: about `6.5e4` monomial products.
- **Nonsplit** (`z^2` not a square): build the quadratic field
  `Q(sqrt(z^2))`. Cost: one PARI `nfinit` of degree 2.
- **Involution kind**: the reversal acts on `z` by `(-1)^{n(n-1)/2}`. Cost:
  one sign, independent of the form.

No branch loops. Termination is unconditional.

## Rank-five model and transfers from Sp_4 and SL_4

Sources: `clifford_rank_five.py:43-154`, `clifford_split_rank_six.py:34-151`.
Klein correspondence `SL_4 -> SO(Lambda^2)` and `Sp_4 -> SO_5` on the
pfaffian complement of `Omega` (Knus-Merkurjev-Rost-Tignol,
*The Book of Involutions*, 15.D).

- `second_compound`: 36 `2x2` minors of a `4x4` matrix.
- `pfaffian_action`: compound plus restriction to the 5-dimensional
  complement.
- Generators: 12 elementary matrices of `SL_4(Z)`, or the standard generators
  of `Sp_4(Z)`; 3 coset reflections of characters `(det, theta)`.
- Spinor character: one reflection factorization over `Q` (Cartan-Dieudonne,
  at most `n` reflections).

The verification of an image modulo 3 is a matrix-group order in GAP. The
groups are `O_6^+(3)` of order `24261120` and `O_5(3)` of order `103680`
(Atlas of Finite Groups, conventions for `O_n(q)`). Schreier-Sims on the
natural module of degree 6 over `GF(3)` has base length at most 6 and
orbit lengths at most `3^6 - 1 = 728`.

## Hermitian lattices of the hyperbolic plane

Sources: `binary_hyperbolic_unitary.py:15-87`,
`eisenstein_hyperbolic_unitary.py:31-228`.

- `word` reduces a unitary matrix by Euclidean division in `O_E`. This
  branch requires `O_E` norm-Euclidean: `E = Q(sqrt -d)`, `d = 1, 2, 3`,
  and `Z[omega]`. Each step decreases `|N(c)|` by a factor at most
  `max_x min_q |N(x - q)| < 1`, so the word length is
  `O(log |N(entries)|)`. For another `d` the step raises `ArithmeticError`.
- The public centralizer of `zeta` on the trace lattice is computed by the
  indefinite orthogonal-group route of `centralizers.md` (owner B). The
  trace lattice of `H` over `Q(zeta_3)` or `Q(zeta_4)` has rank 4,
  signature `(2, 2)` and Witt index 2.

## Quaternion norm lattices

Source: `quaternion_lattices.py:24-127`. `(-1, 3)_Q` is ramified at 2 and 3.

- Multiplication, conjugation and pair actions: constant work in dimension 4.
- `normalizer_classes_in_box(r)`: `(2r + 1)^4` elements, each with 4
  conjugations and at most `c` class comparisons, where `c` is bounded by
  `|N(Lambda) / Q^x Lambda^x| = prod_{p | disc} [local index] = 2 * 6 = 12`
  (Voight, *Quaternion Algebras*, 23.2 and 18.4). For `r = 3`: 2401 elements.

## Lattice isometry over a number field

Source: `genus_number_field.py:182-220`.

- **Definite branch** (totally definite): the trace forms over `Q` of rank
  `d n` are positive definite; PARI `qfisom` decides isometry
  (Plesken-Souvignier, J. Symb. Comp. 24 (1997)). Cost is governed by the
  number of vectors of the minimum and by `|Aut|`. For `I_9 tensor O_K`,
  `K = Q(sqrt 5)`: 18 vectors of norm 1, `|Aut| = 2^9 9!`.
- **Indefinite binary branch**: commutant of the trace forms
  (`_commutant_isometry`, l.153-180). Linear algebra of size `(d n)^2`.
- **Neighbour construction** (`Cusps`, l.225-380): `good_prime` iterates
  ideals of bounded norm until one is coprime to the data; `_residues`
  enumerates `prod D_ii` residues of a Smith form, equal to the index of the
  sublattice (for the test `p = (sqrt 5)`: index 5).

## Semisimple Voronoi

Source: `semisimple_voronoi.py:1-514`; Braun-Coulangeon-Nebe-Schoennenbeck,
J. Algebra 435 (2015), Thm. 3.4; Opgenorth, Exp. Math. 10 (2001).

- **Rational branch** (`splitting is None`, l.126-131): the algebra is
  transpose-stable; the cone is the real symmetric forms of the order.
- **Split branch**: user-supplied real places; herm dimension
  `sum_v n_v(n_v + 1)/2`.

Cost: (number of perfect classes) x (facets per class) x (transporter
candidates). For `M_2(Z)`: one perfect class (`A_2`), 3 facets, `|Min| = 6`
up to sign in each column choice; `_candidates` produces at most
`|Min| (2|Min|)^{k-1}` products with `k` the number of chosen minimal
vectors. Termination follows from finiteness of perfect classes (BCNS
Thm. 3.4). Each step calls PARI `qfminim` on a form of rank `dim O`.

## Witt-index-one Hilbert route

Source: `witt_one_hilbert.py:39-237`. `L = U + (O_E, Tr)` for `E = Q(sqrt D)`,
`O(L)` related to `SL_2(O_E)` (Vaserstein, Mat. Sb. 1972).

- **Recognition** (`recognized_model`, l.185-220): requires `|det L|` equal
  to `disc(E)`. Enumerates the box `[-3, 3]^4` (2400 vectors), pairs of null
  vectors with `b(e, f) = -1`, and at most `13^4 = 28561` changes of basis of
  the complement per pair. Failure returns `None`; it is no statement on `L`.
- **Class number one**: `<SL_2(O), -1, Stab(e)>` is `O(L)`; one cusp.
- **Class number > 1**: cusp links from `witt_one_fixed_norm_classes`; the
  number of cusps is `h(E)`.
- `transporter_to_e`: Euclidean reduction; requires `O_E` norm-Euclidean
  (`D = 2, 3, 5` are). Length `O(log` height`)`.
- `norm_one_unit`: one PARI `bnfinit` of degree 2.

## Spacelike vector classes at positive index one

Source: `marked_lattices.py:583-600`. When `sign * square < 0` the classes are
the primitive embeddings of `<square>` (Nikulin 1980, Prop. 1.15.1), through
the generic genus owner. Cost is
(number of complement genera) x (class number of each) x `|O(q_L)|`. For
`U + <-2>`: `A_L = Z/2`, `|O(q_L)| = 1`. See
`discriminant-groups-and-congruence-orbits.md` (owner F).

## Unit costs assumed

| Primitive | Cost |
| --- | --- |
| Sage exact matrix operation, size at most 8 | 50 microseconds |
| Number-field element operation, degree at most 2 | 50 microseconds |
| Clifford monomial product | 20 microseconds |
| PARI `nfinit` or `bnfinit`, degree 2 | 10 milliseconds |
| PARI `qfisom` or `qfauto`, rank at most 18, few minimal vectors | 5 seconds |
| GAP `Size` of a classical group over `GF(3)`, degree at most 6 | 5 seconds |
| Julia startup and Oscar load | 60 seconds |

## Semisimple Voronoi over a cubic division algebra

Source: `semisimple_voronoi.py:121-150` (split branch), `:251-370`; witness
`witnesses/test_semisimple_voronoi_cubic_division.sage`. The algebra is the
cyclic algebra `D = (L/Q, sigma, 2)` with `L` the cubic subfield of
`Q(zeta_7)`. The prime 2 is inert in `L`, so 2 is not a norm and `D` is a
division algebra of degree 3 (Pierce, *Associative Algebras*, 15.1 and 15.6).
`D` has no positive involution over `Q` (Albert's classification, Mumford,
*Abelian Varieties*, Sec. 21), so the split branch applies: one real place,
`h : D -> M_3(L)`, cone `Sym_3(L)` of dimension 6, order `O` of `Z`-rank 9
from `maximal_order`.

Counts that govern the cost:

- **Unit group.** `O^1 = {nrd = 1}` is a cocompact lattice in `SL_3(R)`
  (Borel-Harish-Chandra; Godement compactness criterion, since `D` is a
  division algebra). `O^x = +-O^1`: `nrd(O^x)` is a unit of `Z`, and `-1` has
  `nrd(-1) = -1` in odd degree. The left-regular determinant on `M_9(Q)` is
  `nrd^3`, which is the `+-1` test.
- **Perfect forms.** Finite up to `O^x` (BCNS Thm. 3.4; Opgenorth 2001).
  The closest tabulated case is BCNS (2015), Sec. 6.4: the degree-3 division
  algebra over `Q(zeta_9 + zeta_9^{-1})`, with invariants `1/3` at 2 and `2/3` at 3
  and a maximal order of discriminant `2^6 3^6`, has 431 classes of perfect forms
  (410 for the opposite algebra, Remark 6.3). The witness algebra is a different
  algebra. `L` is ramified only at 7, and 2 is not a cube modulo 7, so `D` ramifies
  at 2 and 7, not at 2 and 3. The class count `P` is estimated by covolume
  scaling: the covolume of `SL_1` of a maximal order is proportional to
  `prod_{p ram} (p - 1)(p^2 - 1)` (Prasad volume formula, Publ. IHES 69 (1989),
  for inner forms of `SL_3`). The ratio between ramification at 7 and at 3 is
  `(6 * 48) / (2 * 8) = 18`, which gives `P` of about `431 * 18`, roughly 7800.
  This is an estimate. Neither the witness file nor issue #54 records a count.
  `F`, the number of facets per class, is not tabulated in BCNS. It is assumed to
  be at most 30 in the cone of dimension 6.
  Sources searched: BCNS, J. Algebra 435 (2015) and arXiv:1407.6234, Sec. 6;
  Nebe, slides "Computing unit groups of orders" (2014); the witness file;
  the comments on issue #54.
- **Transporters.** For a division algebra one minimal vector generates
  `A` as a right module, so `_candidates` returns `|Min(s)|` candidates
  `m x_1^{-1}` per test (module docstring, Equivalence).
- **Neighbour step.** Schuermann Alg. 2: each probe is a Gram matrix
  `9 x 9` over `L` plus one `qfminim` of rank 9 on a rational minorant.
  The step doubles and then shrinks `rho`. The number of probes is
  logarithmic in the ratio of the final `rho` to the initial bound.

Cost: `P * F * (probes * Gram + |forms| * |Min| * act)`. The Gram matrix
needs 81 entries, each a sum of traces of `3 x 3` products over `L`. That is
about 2500 cubic-field operations. Polyhedron facet enumeration over `L`
(`backend="field"`) uses `|Min|` rays in dimension 6.
