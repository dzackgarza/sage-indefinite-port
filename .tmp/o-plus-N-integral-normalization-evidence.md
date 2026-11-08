# Full Schreier action normalization

The observer of the original O_plus_N line/plane count computation at port 9691908/provider 438c23a5 retains a complete first rank-twelve family in `.tmp/o-plus-N-full-schreier-inputs/family-0.pickle`. Its native Gram and every emitted integral column matrix are retained together, with `complete=True`. This is the original computation's family, distinct from the earlier small published-plane input.

The pinned `MatrixGroupSimplification.h:650-751,1386-1487` tests four pair products and their inverses under entrywise L1 complexity. The production normalization body is measured between native column realization and owned isometry materialization. Its diagnostic shell supplies the actual columns and rank; it constructs no substitute lattice or integral-structure action.

| Complete native normalization | Wall seconds | CPU seconds |
| --- | ---: | ---: |
| Original individual products/reductions | 520.354280 | 508.868164 |
| Batched same products/reductions | 400.818500 | 394.104640 |

Every final ordered integral action key agrees exactly between the two executions. This is complete output equality, including the integral kernel actions. The candidate preserves the same pair order, L1 norms, inverse convention, tie ordering and Nielsen replacements.

Batch multiplication uses integer arrays. If both input directions have total L1 norm at most B=floor(sqrt(2^63-1)), every product's absolute partial sum and total L1 norm is at most B squared. Thus every machine-integer intermediate is exact. Larger directions use Python integer arrays, as in the existing source route. The change batches these same arithmetic operations; it does not bound or discard generators.

All seven existing normalization assertions pass in 9.75 seconds, including the complete unipotent kernel, primitive cyclic parameters, equal-norm permutation generators and large noncommuting integral actions. The original unchanged O_plus_N node remains required at the published repair/provider pair. The earlier observer is terminal at its original 2400-second bound and has no whole-case pass. Its retained family 1 is complete; family 2 is partial.

Optional timing/output inputs: `.tmp/o-plus-N-normalization-{baseline,batched}.{py,json,log}` and `.tmp/o-plus-N-normalization-{baseline,batched}-keys.pickle`. Cloud continuation uses the source and original consumer without these local captures.
