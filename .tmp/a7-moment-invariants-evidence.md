# Exact A7 moment-colour computation

Source: 604d76248237f710c4741ef9911692a7b4845932.
Provider: d0431b3f5c3ca0031dd74ea6fe3df649c269cfe9.

The A7 discriminant-image consumer remains incomplete. Its terminal 2deaefb/2dfd20b6 stack reaches the moment-colour norm sum in `_integral_configuration_invariants`. That stack does not establish the cost of the whole consumer.

A bounded diagnostic at port 975fc42 captures the first large production cell on entry to that operation, then stops. It retains the actual row configuration and actual Gram matrix together, before any invariant computation. This avoids reconstruction through randomized reduction. The optional capture is `.tmp/a7-paired-production-input.pickle`.

| Captured quantity | Value |
| --- | ---: |
| Configuration rows | 49,861 |
| Ambient rank | 10 |
| Gram signature | (1,9) |
| Moment entries compared | 498,610 |
| Moment contents compared | 49,861 |
| Moment norms compared | 49,861 |

The actual Gram is -A7 plus diag(-1,-1,1). For rows R and Gram G, the moment matrix is M=R^t R. A moment row is b=r G M. Its content is gcd(b), and its norm is b G b^t. The matrix identity `(B elementwise_product (B G)) * 1 = diag(B G B^t)` computes every norm over exact ZZ, without a full square product.

The previous moment block size used the number of configuration rows, although these matrices have only the ambient rank as their width. The revised computation uses the same rank-bounded block rule and native exact diagonal product already used by `_ConfigurationPairings.diagonal_blocks`. Each content, norm, class-sum colour and transporter certificate keeps its meaning.

| Same captured input, complete moment contents and norms | Seconds |
| --- | ---: |
| Original scalar sums, 16-row blocks | 1.59108 |
| Scalar sums, 1024-row blocks | 1.26707 |
| Native exact products, 1024-row blocks | 0.64731 |

Every result agrees entry by entry. The entire revised invariant computation also agrees with the saved original output: every full refined row spectrum, Smith factor, total row sum and row-moment entry. Its unprofiled call takes 2.77650 seconds. The original complete call was profiled, so its 6.01209-second duration is not a comparable whole-operation speed measurement.

The existing owner assertion `tests/test_configuration_pairings.py::test_large_coarse_hash_preserves_gcd_colours_under_integral_isometry` passes under the normal local checks. It exercises the moment branch and an integral shear larger than 2^80. The selected session passes in 3.99 seconds.

Required next evidence is the complete unchanged A7 discriminant-image assertion under its original bound, followed by its final-pair frozen result. The measured component improvement does not establish either result.

The unchanged A7 node at that pair is terminal under its original 1800-second bound, in immutable `sage-indefinite-proof-604d762`. Its optional terminal log is `.tmp/a7-604d762-d0431b3f.log`. The process is gone. The final stack reaches the recursive facet-orbit walk, `register -> libgap.RepresentativeAction(group, positions, target, OnSets)` at facet_orbits.py:443. The full image assertion remains incomplete. A preceding live sample reaches small-cell transporters, including exact configuration-preserving seed actions. Optional `.tmp/a7-604d762-active-transporter.pickle` retains an actual source/target pair and both Grams; `.tmp/a7-604d762-active-state.json` identifies that earlier snapshot. It is not the terminal facet-comparison input. Retain the actual group and compared incidences before selecting a repair to that operation.
