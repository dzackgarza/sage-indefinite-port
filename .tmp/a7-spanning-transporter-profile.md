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
