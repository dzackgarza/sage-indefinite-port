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
