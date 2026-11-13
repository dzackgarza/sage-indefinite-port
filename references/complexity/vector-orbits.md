# O(L)-orbits of primitive vectors: complexity

## Algorithm and pinned source

The public consumer is `category.primitive_marking_classes(<n>, L)` for the
coslice `<n>/Lattices`. The port's entry is `_pointed_classes`
(`src/sage_indefinite_port/indefinite/marked_lattices.py:560-614`). Two
classes are compared by `RECURSION.pointed_isometry`
(`recursive.py:2014`).

The pinned sources are `INDEF_FORM_GetOrbitRepresentative`
(`references/vendor/polyhedral_common@1592b246/src_indefinite/INDEF_FORM_GetOrbitRepresentative.cpp`)
and the Lorentzian enumeration of `src_lorentzian/lorentzian_perfect.h`.

## Branches

`marked_lattices.py:560-614` selects one branch:

| Branch | Condition | Route |
| --- | --- | --- |
| Degenerate | `L` degenerate, `s != 0` | `_degenerate_pointed_classes` |
| Anisotropic binary | `n = 2`, signature `(1,1)`, anisotropic | `binary_fixed_norm_representatives` |
| Split binary | `n = 2`, signature `(1,1)`, isotropic | `binary_split_primitive_vector_classes` |
| Isotropic | `m = 1`, `s = 0` | `INDUCTION.isotropic_vector_classes` |
| Embedding | `m = 1`, `sign * s < 0` | Nikulin primitive embeddings (Oscar bridge) |
| Lorentzian cone | `m = 1`, `sign * s > 0` | `MarkedCellOrbitAlgorithm` on the perfect complex |
| Eichler | `m >= 2` | `RECURSION.primitive_vector_classes` |

`s` is the square, `m` the positive index of the signed view and `sign` its
sign.

## Lorentzian cone

Every vector of positive square in the closed positive cone lies in some
perfect cell. The orbits of norm-`s` vectors are the orbits of the vectors in
one cell of each cell orbit under the cell stabilizers. The enumeration is a
Fincke–Pohst search of the vectors of norm `s` with bounded pairing against
the cell's interior vector (Cohen, GTM 138, Alg. 2.7.5). Its volume grows as
`s^((n-1)/2)`. Budget, with one cell for `E10 = U ⊕ E8(-1)`:

`T_traverse(L, 1) + 10 (s/2)^((n-1)/2) U_FP + k (0.05 n^2 + 0.2 T_definite(n-1))`,

`k` being the recorded number of orbits. Each orbit pays one stabilizer
computation through the definite complement `v^⊥`.

`test_published_representatives_are_one_per_orbit` adds `k^2` pointed
isometry tests. Each test is one definite isometry of the rank `n - 1`
complements.

## Embedding branch

Rank-one primitive embeddings `<s> -> L` are classified by the genus of the
complement and the double cosets of the glue (Nikulin, Izv. 14 (1980),
Prop. 1.15.1). The complement has rank `n - 1` and is definite when `L` is
Lorentzian. The cost is the enumeration of its genus (Kneser neighbours, `h`
classes, `h^2` definite isometry tests) plus the glue double cosets in
`O(q)`, which are small for `|A| = |s|`. The route uses the Oscar bridge.

## Eichler branch

When `L` contains `U ⊕ U`, the orbits of primitive vectors of square `s` are
in bijection with the elements `x ∈ A_L` with `q(x) = s/2 mod 2Z` modulo
`O(L)` (Eichler criterion; Gritsenko, Hulek and Sankaran, Compositio 143
(2007), Lemma 3.5). For a unimodular lattice there is exactly one orbit for
each represented `s`. The cost is the model construction, `0.01 n^3`, and one
pointed check per class. For an odd lattice the port uses the even envelope
(`EichlerEnvelope`, `eichler.py`).

## Isotropic branch

Primitive isotropic vectors of a Lorentzian lattice correspond to the cusps
of the perfect complex. The cost is the traversal, plus one check per cusp.

## Known worst cases

- Large `s` in the Lorentzian cone: the enumeration volume grows as
  `s^((n-1)/2)`; at `n = 10` this is `s^4.5`.
- Negative squares with a complement genus of large class number.

## Unit costs of backend primitives

The unit costs are those of `orthogonal-group-generators.md`, "Unit costs of
backend primitives".
