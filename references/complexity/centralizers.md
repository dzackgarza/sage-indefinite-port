# Centralizer images in O(q_L): complexity

## Algorithm and pinned source

The consumer is `f.centralizer_discriminant_image()` for an isometry `f` of
finite order. The provider calls the Oscar adapter
`DzackResearchOscarLatticeAdapter.centralizer_discriminant_image`
(sage-extensions `categories/lattice_engines.py:1393-1418`). The algorithm is
hermitian Miranda–Morrison theory (Brandhorst and Hofmann, "Finite subgroups
of automorphisms of K3 surfaces", §6) when the minimal polynomial of `f` is a
quadratic cyclotomic polynomial on an indefinite lattice. Otherwise it glues
the images of the centralizers of the cyclotomic components
(dzackgarza/research#397).

## Branches

| Branch | Condition | Route |
| --- | --- | --- |
| Hermitian | minimal polynomial `Φ_k`, `φ(k) = 2`, `L` indefinite | `U(L, f)` image by local determinants and spinor norms |
| Gluing | otherwise | images of the components' centralizers, then the glue stabilizer |

## Cost

- Hermitian branch: the image of the unitary group in `O(q_L)` is computed
  from the local unitary groups at the primes of `det L`, through the
  determinant and spinor-norm exact sequences (Miranda and Morrison,
  *Embeddings of integral quadratic forms*, Ch. VIII). The cost is
  polynomial in `n` and `log |det L|` per prime, plus the construction of
  the image as a permutation group on `|A|` points.
- Gluing branch: one centralizer image per cyclotomic component, then the
  stabilizer of the glue subgroup in the product of the images. The
  stabilizer is a permutation-group computation of degree `|A|`.

Budget: `2 + 0.5 n^2 + T_Oq(|A|, log2 |image|) + ord(f) n^3 U_ELT` seconds,
plus the Oscar bridge start. `T_Oq` is defined in
`discriminant-groups-and-congruence-orbits.md` (owner F).

## Known worst cases

Indefinite unitary groups have no generator algorithm in the literature. The
image is therefore not computed from generators of the centralizer.

## Unit costs of backend primitives

The unit costs are those of `orthogonal-group-generators.md`, "Unit costs of
backend primitives".
