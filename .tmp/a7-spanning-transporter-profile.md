# Spanning perfect-cell transporter

The unchanged `tests/test_perfect_cell_transporter.py::test_integral_transporter_of_spanning_perfect_cells` passes at port c2398eb172c46b3eb7438149a9f2d48cc5d0dfa8 / provider 259e62d61cedcd6c1a3dc4deda593b2ef568e6e5. Its explicit unimodular matrix independently certifies the two complete configurations. This specimen exercises the owned small-configuration transporter used in the terminal original A7 stack; it does not establish the complete A7 discriminant image.

The cProfile measurement includes profiler overhead. Cumulative times overlap and must not be added.

| Operation | Cumulative seconds |
| --- | ---: |
| Test call | 71.893 |
| Owned coloured-family transporter | 71.546 |
| Integral-structure transporter | 44.955 |
| Invariant module construction | 20.453 |
| Finite residue representation | 17.912 |
| Source and target containing bounds | 14.280 |
| Source and target spans | 12.569 |
| Lattice Gram contraction | 12.004 |

The selected lattice's finite pairing repeats the coordinate contraction in `Lattices.ParentMethods._evaluate_gram`. The tensor owner also implements this contraction and retains private matrix/vector adapters. Check exact tensor contraction and its lattice consumer before changing the configuration search. All framings, exact coefficients and owned scalar results must remain intact. The profile does not establish that this contraction is the only cost.

The unprofiled same-spanning-transporter comparison at providers 18cc2773 and 16bfac88 passes with call times 37.91 and 37.59 seconds. It gives no material whole-transporter improvement from scalar-action factory reuse.

Provider eb2ab6dd3b27ebed220ad391b87c347898bed516 delegates the finite lattice pairing to the tensor owner. Its exact ZZ/QQ contraction uses Sage matrix/vector products. The full 2 A19 root-difference assertion passes within the unchanged owner call limit; named framings, the zero lattice, polynomial coefficients and rectangular exact contractions pass too. The same spanning perfect-cell transporter passes with call time 38.65 seconds. This gives no material whole-transporter gain over the earlier unprofiled comparison. Complete A7 evidence remains required.

At provider 6b0a47415ae76eeb561d4b45e300c111e3b2c555 the same assertion passes with call time 28.93 seconds. The unchanged original A7 computation at port 1716980/provider 6b0a4741 remains separate evidence. Its live stack reaches both the facet walk's point-permutation construction and the neighbouring-cell transporter.

At port 57faa7359724ebd7f1a19304ad8a26d3912ab969/provider 9f578022486b6e6771ee7df320fe467b5bed9549 the same assertion passes with call time 36.50 seconds. This measurement runs beside the active original cases. It establishes the exact transporter assertion under the transition-table implementation; it supplies no measured whole-transporter gain over 6b0a4741. The additional owner index-two commensurability test did not reach its assertions because its session import exceeded the existing ceiling. It supplies no mathematical result.

A ten-second sample during that original computation places the sampled main thread in the neighbouring-cell transporter. Inclusive times overlap:

| Operation during the sampled interval | Seconds |
| --- | ---: |
| Coloured-family transporter | 9.95 |
| Integral-structure finite representation | 8.30 |
| Finite orbit sweep | 5.75 |
| Induced finite quotient automorphisms | 5.70 |
| Intermediate-lattice image keys | 3.80 |

These interval measurements do not give whole-computation proportions. The selected-generator transition-table realization published in provider 9f578022 addresses reuse after a completed sweep; the sweep and each distinct generator's quotient automorphism still require their exact computations. Preserve the full original assertion before claiming a gain on A7.

For live sampling, the [py-spy documentation](https://github.com/benfred/py-spy#how-can-i-avoid-pausing-the-python-program) states that `--nonblocking` permits read races and incomplete frames. Both nonblocking and default dumps can fail to locate this running interpreter; a standard ten-second recording with py-spy 0.4.3 produced the interval above without sampling errors. A failed dump supplies no stack evidence and does not show that the mathematical process ended.


Port 1078c61 retains the actual indexed family on each perfect cell. Both source and target roles use that family; the provider's identity-based configuration realization can then recognize repeated comparisons. The unchanged spanning-cell assertion passes with call time 39.65 seconds. A separate probe repeats the full comparison of the same source and target, checking the explicit form and vector-set certificate on every returned witness.

| Port / provider | First comparison, seconds | Repeated comparison, seconds |
| --- | ---: | ---: |
| 57faa73 / 9f578022 | 35.530857 | 32.717218 |
| 1078c61 / 9f578022 | 42.436372 | 25.293352 |

These sequential measurements run beside the active original cases. They establish the displayed exact comparisons and do not isolate all wall-time variation. The configuration owner's realization cache uses weak values; retaining the input family does not by itself keep every derived realization alive. Complete A7 evidence remains required.


At port afa1176/provider 25288942, a ten-second live sample reaches the facet-orbit walk. The sampled main thread spends 8.15 seconds in `facet_orbits._orbit`, including 4.15 seconds in its set-image generator. These are interval measurements, not whole-computation proportions.

The [GAP orb manual](https://gap-packages.github.io/orb/doc/chap3.html) supplies complete enumeration and `TraceSchreierTreeForward` in the original generators. A standalone exact comparison used a cyclic permutation of 1024 positions acting on a consecutive half-set. Both routes return exactly the same complete orbit and every word. The boundary conversion dominates the measured GAP route:

| Operation | Seconds |
| --- | ---: |
| Existing Python orbit and words | 0.2181 |
| GAP enumeration | 0.1370 |
| Bulk export of points and Schreier words | 3.7036 |
| Construction of Python incidences and words | 0.2686 |

This specimen does not support replacing the existing orbit realization with that export route. The production facet-orbit implementation is unchanged. Its full A7 result remains required.


A later ten-second sample at afa1176/25288942 reaches large-cell identification throughout. Pairing-column extraction accounts for 4.90 seconds, including 2.60 seconds in its Python conversion generator; inclusive times overlap. Direct export of the one-column integer matrix avoids the intermediate Sage vector and retains the same Python integer tuple.

An isolated measurement on an exact 49861-by-1 integer matrix compares ten complete column exports per repetition:

| Export | Repetition times, seconds |
| --- | --- |
| Vector extraction and generator conversion | 0.2804, 0.2601, 0.2864 |
| Flat entries and builtin integer conversion | 0.2409, 0.2144, 0.2095 |

The complete hyperbolic-quadrant transporter and its coarse-colour stabilizer assertions pass with the flat export. A separate direct production assertion checks every pairing among the six roots of the provider's negative A2 lattice scaled by `2^70+1`, against the owned bilinear form, with and without column reuse. Its simple-root square is exactly `-2*(2^70+1)^2`. These establish the exercised exact export and transporter contracts; the live original A7 computation predates this change.


At port 12b3d2da8986cd659fcdf45b7fbc402038b309f7/provider bc6bd572, the unchanged full spanning-cell transporter assertion passes with call time 28.92 seconds (33.24 seconds total). This pair includes direct exact pairing-column export and native Sage/PARI definite integral groups. The explicit unimodular witness and complete configuration equality remain the assertions. The original A7 discriminant-image assertion remains separate evidence.


Provider 0580f98d evaluates a module map at a selected generator by returning its admitted image when the framing lift is that single coefficient one. This is the defining equation of linear extension. The QQ, QQ[t] and M2(QQ) assertions verify the supplied images, zero, negation and left coefficients with reversed target labels. Module admission, tensor bifunctoriality and direct sums pass as well.

An isolated comparison of the same rational map and basis input, with 500 evaluations per repetition, gives:

| Provider source | Repetition times, seconds |
| --- | --- |
| bc6bd572 | 2.0117, 2.4550, 2.0147 |
| 0580f98d | 0.3551, 0.3888, 0.3263 |

These measure basis evaluation, not a whole transporter or the complete A7 computation. The active A7, I_1_1 and stable Dawes cases predate this provider change.

At port a170e2e91c0759649d48059e2cce9ea5c9893a56/provider 0580f98daaf1e698078429b99f93a0d2e0513ba8, the unchanged full spanning-cell transporter assertion passes with call time 35.72 seconds (44.28 seconds total). This confirms the complete witness under direct selected-generator evaluation. Compared with the 28.92-second call at 12b3d2d/bc6bd572, this measurement supplies no whole-transporter timing gain; the isolated basis-evaluation comparison remains a different measurement. Full A7 evidence remains required.
