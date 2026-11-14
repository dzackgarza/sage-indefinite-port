# Isometry testing: complexity

## Algorithm and pinned source

The public consumers are `L.is_isometric(M)` and `L.Isom(M).an_element()`.
The port's entry is `RECURSION.isometry`
(`src/sage_indefinite_port/indefinite/recursive.py:1204-1236`) and
`_indefinite_isometry` (`recursive.py:1238-1308`). Marked lattices dispatch
through `marked_lattices.isometry` (`marked_lattices.py:255-347`).

The pinned source is `INDEF_FORM_TestEquivalence`
(`references/vendor/polyhedral_common@1592b246/src_indefinite/CombinedAlgorithms.h:979-1238`,
`:1322`) and the covering search of `ApproximateModels.h:1305-1390`. The
Lorentzian kernel is `LORENTZ_TestEquivalence`
(`src_lorentzian/lorentzian_perfect.h:1176`) and
`LORENTZ_TestEquivalenceMatrices` (`:1887-2006`).

## Branches

| Branch | Condition | Decision | Witness |
| --- | --- | --- | --- |
| Degenerate | radical `R != 0` | ranks, then nondegenerate quotients | radical extension |
| Prefilter | always | `LatticePrefilter`, genus for even `L` | none |
| Definite | `m = 0` | PARI `qflllgram`, then `qfisom` | same |
| Low Witt | `m >= 2`, Witt index `< 2` | genus and spinor genus | `low_witt_isometry_witness` |
| Lorentzian | `m = 1` | spinor genus when `n >= 3` | perfect-cell transporter |
| Eichler | `m >= 2`, Witt index `>= 2` | spinor genus when `n >= 3` | covering search |

For an indefinite genus of rank `>= 3` each spinor genus is one class
(Eichler and Kneser; Cassels, *Rational Quadratic Forms*, Ch. 11, Thm. 1.3;
Conway and Sloane, *SPLAG*, Ch. 15, §9). The port uses this to decide a
negative verdict without the covering search (`recursive.py:1255-1285`,
issue #68). The positive verdict still needs a witness.

## Decision cost

The genus symbol needs the local Jordan decomposition at each prime of
`2 det L`: `O(n^3)` arithmetic per prime. The spinor-genus test
(`isometric_within_genus`, sage-extensions `lattice_engines.py:1586`) adds
the local spinor norms at the same primes. Budget:
`T_CALL n (1 + ell)` seconds per genus symbol: one PARI local call per prime
and Jordan block.

## Witness cost

- Lorentzian: traverse the complex of the source (`C(L)` cells), compute the
  initial cell of the target, and test it against each source cell. Cost
  `T_traverse(L) + T_cell + C(L) (T_CALL + n^4 U_ELT)`. The theory of `C(L)` and
  `T_cell` is in `orthogonal-group-generators.md`, "Lorentzian route", and in
  `perfect-domains-and-lorentzian-cells.md` (owner C).
- Eichler: build the target model, enumerate at most `K` covering vectors at
  the norm of the source splitting vector, and test pointed isometries until
  one exists. The spinor-genus decision guarantees that one exists. Worst
  case `K` pointed tests; each one recurses into `v^⊥` of rank `n - 1`.
- Binary: compare the reduction cycles of the two forms, `O(sqrt(D) log D)`.
- Definite over `ZZ`: `_witness_matrix` (sage-extensions
  `categories/lattice_morphisms.py:5271-5330`) LLL-reduces both Gram matrices
  with `qflllgram`. Equal reduced matrices give the witness at once; otherwise
  PARI `qfisom` runs. The cost is `T_definite` of
  `orthogonal-group-generators.md`, "Definite kernel". Oscar's
  `definite_isometry` (`:5327-5331`) runs only when the base ring is not `ZZ`.

## Pointed isometries

`RECURSION.pointed_isometry` (`recursive.py:2014`) decides whether some
isometry carries `v` to `w`. First it compares invariants: `q`, content,
`div` and the class `v/div(v)` in `A_L`. For isotropic `v` it then
transports through the isotropic-stabilizer construction
(`orthogonal-group-generators.md`, "Isotropic stabilizers"). For
anisotropic `v` it tests isometry of `v^⊥` with its gluing. A rejection by
an invariant costs `O(n^3)`. A positive answer costs one isometry test of
the complements and, when the rational extension of that isometry is not
integral, one integral correction (`MatrixIntegral_Equivalence_Bis`,
`CombinedAlgorithms.h:1813-1890`): the source computes `Aut(v^⊥)` for it on
every positive answer, the port only when the correction is needed.

Both hard steps are Plesken-Souvignier backtracks (Plesken and Souvignier,
J. Symbolic Comput. 24 (1997), §3): `qfisom` for the isometry, `qfauto` for
`Aut(v^⊥)`. Their cost is `T_definite(n, S, log2|Aut|)` of
`orthogonal-group-generators.md`, "Definite kernel", where `S` is the number
of vectors of norm at most the largest diagonal entry of the LLL-reduced
Gram matrix. For E10 and `v^2 = 2d > 0` the complement is `E8 ⊕ <2d>`
(negated). On the undecomposed Gram matrix the reduced diagonal contains
`2d`, so `S(d) = 1 + sum_{k <= d} 240 sigma_3(k)` (the theta series of E8):
`S(7) = 199,921`, and `S` grows as `d^4`.
Measured `qfauto` there (14 November 2026, local container): 26.5 s at
`d = 7`, 62.5 s at `d = 8`, 148 s at `d = 9`. On the decomposition into
indecomposables (`E8` and `<2d>`, unique by Eichler) `S = 240 + 2` and
`|Aut| = 2 |W(E8)|` for every `d`; sage-extensions computes `Aut` that way
(`_indecomposable_summands`), measured 0.02 s at every `d` from 7 to 10.
`qfisom` on the LLL-reduced pair measured 0.001-0.035 s for every E10
representative: the reduced pair of a reflected representative is
usually equal, so no backtrack runs.

## Pairing slices

`pairing_slices.py` enumerates orbits of tuples of vectors with prescribed
Gram matrix, slice by slice in the values of a fixed pairing. Each slice is a
finite set of lattice points, because the pairing bounds a positive
definite quadratic form on that slice. An infinite family is separated
when every slice is finite. The tests use rank `<= 4` and slice sizes below
`10^2`. Budget: 4 s per test, 6 s with a centralizer (one extra finite
image computation).

## Shared work across cases

The witness search memoizes the reduced source lattice and its complex
(`IsometryGroupoidCache`, `groupoid_cache.py`; `recursive.py:1238-1252`). The
145 cases of `test_lorentzian_pairs_are_isometric` use only 2 distinct source
lattices `mat1`: 142 with `|A| = 2^8` and 3 with `|A| = 2^10`, both rank 10 and
2-elementary. So `d = 4` or `d = 1`, and `C <= 8` or `C <= 5`. The source
traversal is charged once, to the first case of each source lattice in
collection order. Every case pays two genus symbols, the target's initial cell
and `C` cell-transporter tests.

## Known worst cases and non-termination

- Same-genus pairs of large 2-adic index without the spinor shortcut: the
  covering expansion of the source exceeds `10^6` cosets for
  `diag(1, -2, -64, 128)` (`recursive.py:1238-1252`). With the shortcut the
  negative verdict costs one spinor-genus test.
- Lorentzian pairs with many cells: the cost is quadratic in `C(L)`.

## Unit costs of backend primitives

The unit costs are those of `orthogonal-group-generators.md`, "Unit costs of
backend primitives".
