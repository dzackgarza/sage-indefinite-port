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
