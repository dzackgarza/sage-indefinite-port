# Based Eichler envelope action

The pinned `polyhedral_common@1592b246/src_indefinite/ApproximateModels.h:1309-1349` constructs the stabilizer of the supplied embedding, then conjugates every stabilizer generator into that embedding's source frame. The complete integral subgroup includes the kernel of its finite action.

`EichlerEnvelope.stabilizer_of_lattice` retains the actual module morphism from L into the restriction of the rational envelope. Its basis images are precisely the rational isometry's images of the supplied basis of L. Its image, rational group and stable containing module are the same mathematical data as before. The frame now records L itself.

For a proper envelope, `OrbitCoverModel.approximate_generators` uses the existing complete Schreier factor realization in that frame, followed by the existing source Nielsen reduction, before constructing owned isometries. Each factor product is conjugated by the same basis comparison; the full group therefore agrees with the source. Nielsen products, inverse choices, identity removal and duplicate removal preserve the generated integral subgroup, including its finite-image kernel. The integral-identification and represented-split routes retain their source families.

| Mathematical assertion | Result |
| --- | --- |
| Every ordered native Schreier action on U(2)+U equals its owned rational action on each supplied basis vector | Pass |
| Every inverse commutes with the based inclusion | Pass |
| Every action preserves the Gram form and has determinant a unit | Pass |
| Existing integral kernel, primitive cyclic, permutation and large exact-integer assertions | Pass |
| Existing full family realization and noncommuting lattice-action assertions | Pass |

The first three assertions pass in 40.68 seconds. The existing normalization assertions pass in 6.56 seconds; the other complete-family assertions pass in 20.85 seconds. These establish the exercised source operations. The original O_plus_N and Dawes public consumers remain required at the published port/provider pair.

The prior actual rank-six Dawes model input retains its original Gram, rational embedding, envelope Gram and complete source generating matrices. Its optional reconstruction benchmark uses actual owned lattice and integral-structure objects; it does not fabricate an envelope model. Timing of that operation does not establish an entire public boundary computation.
