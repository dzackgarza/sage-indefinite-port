# Full Schreier action normalization

The observer of the original O_plus_N line/plane count computation at port 9691908/provider 438c23a5 retains a complete first rank-twelve family in `.tmp/o-plus-N-full-schreier-inputs/family-0.pickle`. Its native Gram and every emitted integral column matrix are retained together, with `complete=True`. This is the original computation's family, distinct from the earlier small published-plane input.

The pinned `MatrixGroupSimplification.h:650-751,1386-1487` tests four pair products and their inverses under entrywise L1 complexity. The production normalization body is measured between native column realization and owned isometry materialization. Its diagnostic shell supplies the actual columns and rank; it constructs no substitute lattice or integral-structure action.

| Complete native normalization | Wall seconds | CPU seconds |
| --- | ---: | ---: |
| Original individual products/reductions | 520.354280 | 508.868164 |
| Batched same products/reductions | 400.818500 | 394.104640 |

Every final ordered integral action key agrees exactly between the two executions. This is complete output equality, including the integral kernel actions. The candidate preserves the same pair order, L1 norms, inverse convention, tie ordering and Nielsen replacements.

Batch multiplication uses integer arrays. If both input directions have total L1 norm at most B=floor(sqrt(2^63-1)), every product's absolute partial sum and total L1 norm is at most B squared. Thus every machine-integer intermediate is exact. Larger directions use Python integer arrays, as in the existing source route. The change batches these same arithmetic operations; it does not bound or discard generators.

All seven existing normalization assertions pass in 9.75 seconds, including the complete unipotent kernel, primitive cyclic parameters, equal-norm permutation generators and large noncommuting integral actions. The original unchanged O_plus_N node remains required at the published repair/provider pair. The earlier observer is terminal at its original 2400-second bound and has no whole-case pass. Its retained family 1 is complete. Family 2 started, but its native input is not retained in the inspected capture files. The later 1a21c88/cfee1892 whole node also reaches its original bound, with its terminal stack in batched normalization. The first-family equality does not establish the later computation or whole node.

Optional timing/output inputs: `.tmp/o-plus-N-normalization-{baseline,batched}.{py,json,log}` and `.tmp/o-plus-N-normalization-{baseline,batched}-keys.pickle`. Cloud continuation uses the source and original consumer without these local captures.

## Necessary norm condition on the current complete input

The complete native family retained by the original 8e5436f/cfee1892 run in `.tmp/o-plus-N-8e5436f-native-stream/integral_columns-18.pickle` has an explicit completion footer. It is not identical to the older coupled Gram/family input above. The new measurement uses the current input throughout.

For a normalized pair product c, its stored complexity is the minimum of the direct and inverse L1 norms. It can affect the existing source step only by strictly lowering a pair's maximum norm or by equaling an already retained full integral action. The latter equality forces equal stored norm. Products failing both necessary conditions cannot change the step's result, so their full matrix keys need no construction. All source products and inverse products are still computed exactly. An exact multiplicity index follows every insertion and retirement in the retained action dictionary. This preserves its pair order, redundancy decisions, strict improvements and tie rules.

| Complete current native normalization | Wall seconds | CPU seconds |
| --- | ---: | ---: |
| Original source body | 394.795511 | 391.801627 |
| Necessary norm condition before key construction | 351.009311 | 348.096427 |

Every ordered final integral action key agrees. The source input, full integral kernel and each arithmetic branch are retained. The existing kernel, primitive cyclic, equal-norm permutation and large noncommuting integer assertions pass in 7.99 seconds at source 1e3ef33. Optional complete data: `.tmp/o-plus-N-norm-index-{baseline,candidate}.{log,json}`, their key files and `.tmp/o-plus-N-norm-index-equality.json`.

A distinct prior-action batching probe is slower on its scoped prefix and is rejected; its prefix is not whole-family evidence. The original public O_plus_N node remains required.

## Later public-consumer operation

The original 8e5436f/cfee1892 O_plus_N node is terminal at its unchanged bound. Its final stack is in `FiniteGSets.orbit_generator_permutations`, repeatedly ranking point-map images while realizing the later plane stabilizer's Schreier family. This is separate from the completed first native family above. The provider now admits its already retained integer transition table directly through the existing symmetric-group constructor; `results/orbit-generator-transitions.md` owns the defining action and point-label naturality checks. The final same-pair public result remains required.
