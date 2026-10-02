# TODO

The work is the GitHub issue tree rooted at
[#4 Ledger: Indefinite lattice port](https://github.com/dzackgarza/sage-indefinite-port/issues/4).
Each milestone below is one phase; each line is one work unit, a leaf issue whose body holds its
obligation, acceptance, and proof. **Needs** lists the work units that must close first (GitHub
blocked-by links carry the same edges). Already complete: the capability inventory and oracle
corpus, the shared `LatticePrefilter`, and dependency wiring (T0).

Take the next unit with `uvx --from git+https://github.com/dzackgarza/itree itree next
dzackgarza/sage-indefinite-port`: the first open work unit in tree order whose Needs are closed.
T1, T2 and T3 are independent; T2 is listed first because it is in progress.

## T2: Integral structures and double cosets ([#25](https://github.com/dzackgarza/sage-indefinite-port/issues/25))

- [#10](https://github.com/dzackgarza/sage-indefinite-port/issues/10) `IntegralStructureAction`: in progress; the remaining work is the two oracle families. **Needs:** none.
- [#11](https://github.com/dzackgarza/sage-indefinite-port/issues/11) Construction-aware arithmetic subgroup carriers: in progress. **Needs:** none.

## T1: Subobjects, reductions, and exact lifts ([#24](https://github.com/dzackgarza/sage-indefinite-port/issues/24))

- [#5](https://github.com/dzackgarza/sage-indefinite-port/issues/5) `IsotropicReduction` K_I = I^⊥/I. **Needs:** none.
- [#6](https://github.com/dzackgarza/sage-indefinite-port/issues/6) Codimension-one extension, extension equation, `IsometryExtensionTorsor`. **Needs:** #5.
- [#7](https://github.com/dzackgarza/sage-indefinite-port/issues/7) Pointwise perpendicular kernel. **Needs:** #5, #6.
- [#8](https://github.com/dzackgarza/sage-indefinite-port/issues/8) `VectorOrthogonalSection` and rational lifts. **Needs:** #5, #6.
- [#9](https://github.com/dzackgarza/sage-indefinite-port/issues/9) `FlagType` and the preamble `IsotropicFlag`. **Needs:** none.

## T3: Lorentzian cell backend ([#26](https://github.com/dzackgarza/sage-indefinite-port/issues/26))

- [#12](https://github.com/dzackgarza/sage-indefinite-port/issues/12) Configuration isomorphism, cell stabilizers, bucket keys through Bliss. **Needs:** none.
- [#13](https://github.com/dzackgarza/sage-indefinite-port/issues/13) Facets, rays, facet orbits through Sage polyhedra. **Needs:** none.
- [#14](https://github.com/dzackgarza/sage-indefinite-port/issues/14) Definite leaf through the preamble. **Needs:** none.
- [#15](https://github.com/dzackgarza/sage-indefinite-port/issues/15) Lorentzian perfect-cell backend, complex traversal, marked-cell orbits. **Needs:** #12, #13, #14.

## T4: Eichler covers and the recursive full group ([#27](https://github.com/dzackgarza/sage-indefinite-port/issues/27))

- [#16](https://github.com/dzackgarza/sage-indefinite-port/issues/16) Eichler transvections, `OrbitCoverModel`, `EichlerEnvelope`, hyperbolic pairs, splitting vectors. **Needs:** T1 (#7, #8, #9), T2 (#10, #11), T3 (#15).
- [#17](https://github.com/dzackgarza/sage-indefinite-port/issues/17) Reduced presentations, bucket keys, isometry-groupoid cache. **Needs:** T1 (#7, #8, #9), T2 (#10, #11), T3 (#15).
- [#18](https://github.com/dzackgarza/sage-indefinite-port/issues/18) `IndefiniteOrthogonalAlgorithm`: O(L), isometry, vector stabilizers and transporters. **Needs:** #16, #17.
- [#19](https://github.com/dzackgarza/sage-indefinite-port/issues/19) Vector orbit decompositions and primitive isotropic orbits. **Needs:** #16, #18.

## T5: Exact parabolics and rank-two isotropic planes ([#28](https://github.com/dzackgarza/sage-indefinite-port/issues/28))

- [#20](https://github.com/dzackgarza/sage-indefinite-port/issues/20) `IntegralParabolicDatum` and exact parabolic stabilizers. **Needs:** #19.
- [#21](https://github.com/dzackgarza/sage-indefinite-port/issues/21) Inductive isotropic sublattice and flag orbits by double cosets. **Needs:** #20.
- [#3](https://github.com/dzackgarza/sage-indefinite-port/issues/3) Vector-stabilizer and isotropic-k-plane-equivalence kernels as `research` port realizations. **Needs:** #18, #21.

## T6: Finite-index arithmetic groups ([#29](https://github.com/dzackgarza/sage-indefinite-port/issues/29))

- [#22](https://github.com/dzackgarza/sage-indefinite-port/issues/22) Kernels, finite preimages, and orbit splitting. **Needs:** #21.

## T7: Centralizers and intersections ([#30](https://github.com/dzackgarza/sage-indefinite-port/issues/30))

- [#23](https://github.com/dzackgarza/sage-indefinite-port/issues/23) `EquivariantLattice` and centralizers. **Needs:** #22.
