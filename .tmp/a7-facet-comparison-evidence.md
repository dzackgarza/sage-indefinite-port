# A7 facet-comparison diagnosis

The unchanged frozen A7 node at port 604d76248237f710c4741ef9911692a7b4845932/provider d0431b3f5c3ca0031dd74ea6fe3df649c269cfe9 is terminal under its original bound. Its last frame is the recursive facet walk's RepresentativeAction comparison. That frame alone does not establish cumulative cost.

A preceding read-only snapshot retains an actual source/target pair and both Grams in `.tmp/a7-604d762-active-transporter.pickle`. The target has 1505 rows of rank ten. It is not verified as the terminal facet input.

The same source's facet walk on that retained target completes under a separate in-process bound. A local observer retains the actual permutation actions, rays and compared incidences of the recursive subproblems in `.tmp/a7-facet-capture-9691908/`. Its root result retains the exact chosen generator permutations and returned facet-orbit data. The following timings include diagnostic instrumentation in the whole walk; the native comparison timers exclude capture writes.

| Measured operation | Seconds |
| --- | ---: |
| Complete isolated facet walk | 30.22028 |
| All native set-transporter comparisons | 10.40423 |
| Comparisons on exactly repeated incidences | 0.73702 |

An exact duplicate cache would address only the last measured part. This result does not justify selecting it as the repair for the full original timeout. The individual slowest captured comparison in the first isolated execution takes 0.16078 seconds. A complete isolated walk also does not establish the original quotient traversal or its discriminant image.

The cumulative-phase diagnostic at port 9691908/provider 438c23a5 is terminal at its 300-second bound. Its completed facet phase includes stabilizer computations, so completed phase durations overlap. The final stack is in `point_permutations_from_rows`, converting the cached processed representative's full stabilizer action again. The first retained cell has 49,861 rows and its actual Gram is retained in `.tmp/a7-cumulative-phase-capture/facet-cell-0.pickle`. This identifies a repeated owning operation for separate measurement. `.tmp/a7-cell-action-retention-evidence.md` records that measurement and its source-faithful repair.

## Exact certificate pairing measurement

A later live A7 sample at bb69c92/aa75796e includes Python ray pairings in facet certification as well as native set-transporter comparisons. This is a window sample, not a cumulative attribution. A separate probe reconstructs the facet normal from the retained twenty-ray incidence of the earlier 49,861-row A7 cell, orients it against every ray and verifies its complete zero set and exact facet certificate. Its actual rows and incidence are retained in optional `.tmp/a7-bf43-current-cell.pickle` and `.tmp/a7-bf43-orbit-input.pickle`.

Replacing only the Python pairing expression with a freshly constructed exact Sage integer matrix/vector product preserves the certificate but increases the measured operation cost:

| Complete certificate | Seconds, repeated executions |
| --- | --- |
| Current Python pairings | 0.276526, 0.221446, 0.195215 |
| Fresh Sage matrix/vector product | 0.814945, 0.546577, 0.592287 |

This result does not justify that replacement. Production certificate code is unchanged. The full original A7 computation at bb69c92/aa75796e remains active under its original bound. Optional probe and timing: `.tmp/a7-facet-certificate-pairing-benchmark.py`, `.tmp/a7-facet-certificate-pairing-benchmark.json`.
