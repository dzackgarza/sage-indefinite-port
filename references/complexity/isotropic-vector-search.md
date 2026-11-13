# Isotropic vectors and orbits of totally isotropic sublattices

This file states the cost of deciding whether an indefinite lattice `L` of rank
`n` is isotropic, of finding a primitive isotropic vector, and of listing one
representative per `O(L)`-orbit of primitive totally isotropic sublattices of
rank `k` and of isotropic flags. The Lorentzian perfect tessellation is owned
by `perfect-domains-and-lorentzian-cells.md`; isometry tests by
`isometry-testing.md`; splitting of `O(L)`-orbits along a subgroup by
`discriminant-groups-and-congruence-orbits.md`.

No measured runtimes are used.

## Sources

| Step | Pinned source | Port |
| --- | --- | --- |
| Dispatcher by positive index | `polyhedral_common@1592b246/src_indefinite/CombinedAlgorithms.h:1702-1706, 2473-2492` (`INDEF_FORM_GetOrbit_IsotropicKplane`) | `src/sage_indefinite_port/indefinite/isotropic_induction.py:160-200` |
| Lorentzian null vertices | `src_lorentzian/lorentzian_perfect.h:1625-1644` | `isotropic_induction.py:186-189` |
| Eichler recursion, `k`-planes and flags | `CombinedAlgorithms.h:180-491, 822-1400` (`INDEF_FORM_Rec_IsotropicKplane`, `f_stab`, `f_equiv`) | `isotropic_induction.py:250-780`, `exact_isotropic.py` |
| Public entry | | `tests/acceptance/consumer.py:212-230` (`primitive_marking_classes`) |

## Unit costs

`OP`, `STAB(n)` and `BASE` are the units of `vinberg-and-edgewalk.md`. A local
Hilbert symbol at one prime costs `20 OP`.

## Branch selection

| Branch | Condition | Cost driver |
| --- | --- | --- |
| Definite | positive index 0 | none; empty result |
| Anisotropic indefinite | `n <= 4` and a Hasse–Minkowski local condition fails | local symbols at `p | 2 det L` |
| Binary hyperbolic | `n = 2`, signature `(1, 1)`, `-det L` a square | two isotropic lines up to sign; `O(1)` |
| Lorentzian | positive index 1, isotropic | perfect-cell orbits with null vertices |
| Eichler recursion | positive index `>= 2` | candidate classes from `A_L`, reductions `I^perp/I` |

Isotropy: an indefinite form of rank `n >= 5` is isotropic (Meyer 1884; Serre,
*A Course in Arithmetic*, IV.3.2, Cor. 2). For `n = 3, 4` it is isotropic iff
it is isotropic over every `Q_p` (Hasse–Minkowski, Serre IV.3.2, Thm. 8), and
only `p | 2 det L` can fail (Serre IV.2.2, Thm. 6). The decision costs
`BASE + omega(2 |A_L|) * 20 OP`. A vector is found by Simon's algorithm in time
polynomial in `log |det L|` once `det L` is factored (Simon, *Solving quadratic
equations using reduced unimodular quadratic forms*, Math. Comp. 74 (2005)).

## Lorentzian branch

`L` of signature `(n-1, 1)`. After indefinite reduction, the null vertices of
the perfect tessellation of the positive cone are enumerated up to `O(L)`.
Each cusp orbit of `O(L)` is a null vertex orbit. The work is the number `C`
of perfect-cell orbits visited times the work per cell, plus the transporter
tests to identify cells.

The port follows `MarkedCellOrbitAlgorithm.global_orbits`
(`src/sage_indefinite_port/indefinite/lorentzian_cells.py:3514-3560`): one
local stabilizer per quotient cell, one transport of marks per adjacency, and a
union of local orbits along adjacencies. The quotient cells and their
identification are owned by `perfect-domains-and-lorentzian-cells.md`.

For a reflective `L` (every fixture form in `isotropic_cases_8821.json` is a
reflective form of `ListReflect`), the estimate takes the number `C` of
quotient cells to be at most the number of vertices of the Coxeter polyhedron,
computed exactly from the recorded simple roots as in
`vinberg-and-edgewalk.md`, "Main walk"; without recorded roots, `C <= R`.
A cell has at most `n` adjacencies, and a new cell is identified against at
most two representatives of its invariant class.

```
T = BASE + omega(2|A_L|) 20 OP + 0.2 n + C (STAB(n) + 40 OP n) + C n (2 STAB(n) + 10 OP)
```

The rank-18 unimodular form has the largest bound, through `STAB(18)`. A non-reflective
Lorentzian lattice has no such bound from `R`; then `C` is bounded by the mass
of the genus, which grows as `|det L|^((n+1)/2)` up to local densities
(Siegel mass formula; Conway–Sloane, *Low-dimensional lattices IV*, Proc. R.
Soc. A 419 (1988)).

## Eichler branch

Positive index `>= 2`. For a primitive isotropic `v`, the divisor
`div(v)` with `b(v, L) = div(v) Z` and the class `v/div(v)` in `A_L` determine
the `O(L)`-orbit when `L` contains `U + U` (Eichler criterion; Gritsenko–Hulek–
Sankaran, *Abelianisation of orthogonal groups*, Prop. 3.3). So the orbits of
primitive isotropic vectors are in bijection with pairs `(d, x)`, `d | exp(A_L)`,
`x` an isotropic element of order `d` of `A_L` up to `O(q_L)`. Their number is
at most `sum_{d | exp} 2^omega(d)` for the fixtures, and is 2 for
`N = U + U(2) + E8(-2)` (`d = 1, 2`).

Planes and flags recurse: each line class `v` gives the reduction
`v^perp/v` of rank `n-2` and the lines of the reduction modulo the stabilizer
of `v`. The plane classes of `N` number 2 under `O(N)`; under the stable group
they split into `|O(q_N)| / |image of the stabilizer|` orbits per class
(24242 planes, 72199 flags, from Dutour Sikirić–Hulek, arXiv:2302.01679,
lines 1159-1166). That split is owned by
`discriminant-groups-and-congruence-orbits.md`.

For unimodular `L` (`I_{2,10}`) the Eichler criterion gives one orbit per
isomorphism class of `v^perp/v`, two for lines and two for planes
(Allcock fixture `allcock_i_2_10_orbits.json`). The cost is two reductions of
rank 10 and 8 and one isometry test per pair of candidates.

## Worst cases

- `|A_L|` with many primes and large `p`-ranks: the candidate count grows with
  the number of `O(q_L)`-orbits of isotropic elements, and each isometry test
  carries a lift through `O(q_L)`.
- Lattices without `U + U`: the Eichler criterion does not apply, and the
  recursion must test isometry between all candidates of equal invariant.
- The algorithm always terminates: each orbit set is finite (Borel–Harish-Chandra
  finiteness of cusps).
