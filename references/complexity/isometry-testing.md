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
| Definite | `m = 0` | Plesken–Souvignier isometry | same |
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
`0.1 + 0.05 n (1 + ell)` seconds per genus symbol.

## Witness cost

- Lorentzian: traverse the complex of the source (`C(L)` cells), compute the
  initial cell of the target, and test it against each source cell. Cost
  `T_traverse(L) + T_cell + C(L) * 0.05 n`. The theory of `C(L)` and
  `T_cell` is in `orthogonal-group-generators.md`, "Lorentzian route", and in
  `perfect-domains-and-lorentzian-cells.md` (owner C).
- Eichler: build the target model, enumerate at most `K` covering vectors at
  the norm of the source splitting vector, and test pointed isometries until
  one exists. The spinor-genus decision guarantees that one exists. Worst
  case `K` pointed tests; each one recurses into `v^⊥` of rank `n - 1`.
- Binary: compare the reduction cycles of the two forms, `O(sqrt(D) log D)`.
- Definite: one Plesken–Souvignier isometry search.

## Known worst cases and non-termination

- Same-genus pairs of large 2-adic index without the spinor shortcut: the
  covering expansion of the source exceeds `10^6` cosets for
  `diag(1, -2, -64, 128)` (`recursive.py:1238-1252`). With the shortcut the
  negative verdict costs one spinor-genus test.
- Lorentzian pairs with many cells: the cost is quadratic in `C(L)`.

## Unit costs of backend primitives

The unit costs are those of `orthogonal-group-generators.md`, "Unit costs of
backend primitives".
