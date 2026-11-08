# Public complete isotropic flag adapter

Provider source: `012f397fa94b87e5a674cd8cbf33e2480d216ffb`, with the additional impossible-rank assertion at `391b2cd3d9337ae4e47b1e40d0497828b727b738`.
Port source: the published general flag implementation through `f5adaef`; subsequent port changes in this check are assertions and documentation.

`IsotropicFlag` retains every complete prefix, and the port already accepts general finite-ordinal flag markings. The provider's public orbit and stabilizer adapter still rejected every top rank except two. The committed red assertions exercise those literal refusals. The adapter now builds the full ordinal shape, routes it through the existing `primitive_marking_classes` operation, and reconstructs an adapted basis from each successive primitive rank-one quotient. This is the public `seq_dims_flag(rank)` domain of pinned `CombinedAlgorithms.h:2424–2479`.

Unimodularity and evenness give an integral hyperbolic completion to a primitive maximal null basis in nU. Its integral Levi group carries each complete primitive flag in that null lattice to the coordinate flag. Thus each displayed nU has one full-O complete flag orbit. The assertions also check the actual forward/inverse term transporters and form preservation. The stabilizer assertion checks every returned generator on every term in both directions and on the ambient form. It does not separately prove completeness of that generating family. Signature (2,2) makes a rank-three null flag impossible.

| Assertion in `tests/lattices/test_complete_isotropic_flag_orbits.sage` | Mathematical result | Call seconds | Existing 20-second gate |
| --- | --- | ---: | --- |
| `test_complete_null_flag_orbits_and_transporters_in_split_unimodular_lattices[rank_one]` | pass, U | 7.45 | within limit |
| `test_complete_null_flag_orbits_and_transporters_in_split_unimodular_lattices[rank_two]` | pass, 2U | 34.20 | failure |
| `test_complete_null_flag_orbits_and_transporters_in_split_unimodular_lattices[rank_three]` | pass, 3U | 165.83 | failure |
| `test_complete_rank_three_flag_stabilizer_preserves_every_term` | pass, 3U | 26.11 | failure |
| `test_two_hyperbolic_planes_have_no_rank_three_isotropic_flag` | pass, empty set | 3.48 | within limit |

The rank-one/rank-two/stabilizer session executes every assertion, reports three passes and exits 1 on its unchanged time gate. A separate rank-three diagnostic uses a 240-second observation bound, executes every assertion, reports one pass and exits 1 on the unchanged 20-second time gate. Its earlier default 60-second attempt supplied no rank-three result. The impossible-rank session exits 0. These owner performance defects remain open; this is mathematical functional evidence, not a green owner QC result.

The frozen suite is unchanged. Partial flags retain the canonical finite-ordinal marking interface. Full source parity, requested subgroup actions, all quotient/gluing branches and final-pair frozen evidence remain required. The separate SO(4U) moving-marking assertion still lacks a completed ambient generator computation.
