# TODO

The work is the GitHub issue tree rooted at [#4 Ledger: Indefinite lattice port](https://github.com/dzackgarza/sage-indefinite-port/issues/4). Each line is one work unit, a leaf issue whose body holds its obligation and acceptance.
**Needs** lists the work units that must close first (GitHub blocked-by links carry the same edges).

Correctness is measured only against external facts: the reference implementation's recorded outputs, published results, and the oracle corpus derived from them.
There is one acceptance suite (#33), written up front; every capability unit's acceptance is a named set of its cases going green.
No unit adds tests.

Take the next unit with `uvx --from git+https://github.com/dzackgarza/itree itree next dzackgarza/sage-indefinite-port`: the first open work unit in tree order whose Needs are closed.

## Oracle corpus and acceptance suite

- [#32](https://github.com/dzackgarza/sage-indefinite-port/issues/32) Oracle corpus: re-derive every fixture from its primary source.
  **Needs:** none.

- [#33](https://github.com/dzackgarza/sage-indefinite-port/issues/33) Acceptance suite: every consumer operation against the corpus, written up front.
  **Needs:** #32.

## Closed on internal tests; accepted only when their #33 cases are green

T1 (#5–#9), T2 (#10, #11), T3 (#12–#15), and #31 (edgewalk).
A red #33 case for Lorentzian perfect domains (#15), edgewalk (#31), integral structures and double cosets (#10, #11), or the perpendicular kernel (#7) reopens that unit.

## T4: Eichler covers and the recursive full group ([#27](https://github.com/dzackgarza/sage-indefinite-port/issues/27))

- [#18](https://github.com/dzackgarza/sage-indefinite-port/issues/18) O(L), isometry, vector stabilizers and transporters.
  **Needs:** #33.

- [#19](https://github.com/dzackgarza/sage-indefinite-port/issues/19) Vector orbit decompositions and primitive isotropic orbits.
  **Needs:** #18, #33.

## T5: Exact parabolics and rank-two isotropic planes ([#28](https://github.com/dzackgarza/sage-indefinite-port/issues/28))

- [#20](https://github.com/dzackgarza/sage-indefinite-port/issues/20) IntegralParabolicDatum and exact parabolic stabilizers.
  **Needs:** #19, #33.

- [#21](https://github.com/dzackgarza/sage-indefinite-port/issues/21) Isotropic sublattice and flag orbits by double cosets.
  **Needs:** #20, #33.

- [#3](https://github.com/dzackgarza/sage-indefinite-port/issues/3) Preamble entry points for vector stabilizers and isotropic subspaces.
  **Needs:** #18, #21, #33.

## T6: Finite-index arithmetic groups ([#29](https://github.com/dzackgarza/sage-indefinite-port/issues/29))

- [#22](https://github.com/dzackgarza/sage-indefinite-port/issues/22) Kernels, finite preimages, and orbit splitting.
  **Needs:** #21, #33.

## T7: Centralizers and intersections ([#30](https://github.com/dzackgarza/sage-indefinite-port/issues/30))

- [#23](https://github.com/dzackgarza/sage-indefinite-port/issues/23) EquivariantLattice and centralizers.
  **Needs:** #22, #33.
