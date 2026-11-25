# Retained perfect-cell stabilizer action

Pinned `Lorentzian.g:1320-1322` retains `PermutationStabilizer` in the cell stabilizer record. The port already retained the integral generators but converted their action on the same ordered cell rows again for facet walks, cached-target transporter seeds and `_CellFaces`.

The original A7 cumulative diagnostic at port 9691908/provider 438c23a5 terminates inside this repeated conversion. Its captured first cell has 49861 rows in ambient rank ten, with the actual Gram, allowed squares and row order retained together. The isolated conversion uses that production input and provider aa75796e. It first computes the complete integral stabilizer, then measures the complete ordered action twice.

| Operation | Seconds |
| --- | ---: |
| First conversion before retention | 6.294495 |
| Repeated conversion before retention | 5.386605 |
| First conversion with retention | 6.465282 |
| Repeated retrieval with retention | 0.00001792 |

Both executions return exactly the same complete permutation tuple, including generator order and every ray image. Every action is independently checked as a permutation of the entire row-index set. The general operation for a caller-supplied ordered point family still evaluates that family directly. `cell_stabilizer_permutations` retains the action only for the immutable cell rows and the already retained generator family, using Sage's cached method.

The selected octahedral cached-target transporter, exact facet-orbit/transversal and owned-halfspace assertions pass together in 7.42 seconds. They establish the exercised transporter, facet and face consumers. The full original A7 discriminant-image computation remains required; these timings establish the cost of this operation, not completion of the quotient traversal.

Optional native input and timing records: `.tmp/a7-cumulative-phase-capture/facet-cell-0.pickle`, `.tmp/a7-cell-action-retained.pickle`, `.tmp/a7-cell-action-retention-before.json`, `.tmp/a7-cell-action-retention-after.json`. Source and public consumers suffice to resume without these local captures.
