# Exact Eichler covering-class action

Pinned `polyhedral_common@1592b246/src_indefinite/ApproximateModels.h:313-405` reduces covering discriminant classes under GRPeasy. Lines 337-356 apply each easy generator to lifted class coordinates and then reduce modulo integers. The integral class order is preserved by every lattice isometry.

The port forms all scaled lift rows in the supplied complement basis and multiplies them by the actual complement inclusion. For each actual easy isometry, its owned tensor gives its matrix on the envelope basis. Multiplication of the complete row family by that matrix equals evaluation of that isometry on every included lift, by linearity. Entrywise reduction modulo the preserved class order therefore gives exactly the same class permutation. The existing finite group computes its orbits and selects the same representative positions.

The live original O_plus_N window at port 8e5436f/provider cfee1892 spends time in the old owned-map class action. A later snapshot finds no active lift context; it supplies no claim about retained live lifts. A separate snapshot retains the actual cached model's Gram, complement Gram, based inclusion and complete easy matrices. The measured reconstruction uses these exact mathematical data and regenerates the norm-zero discriminant lifts. Every ordered class image agrees between the owned and native computations, and each image family is a permutation.

| Complete action on the retained model | Value |
| --- | ---: |
| Admissible discriminant classes | 136 |
| Easy generators | 7 |
| Owned-map action, seconds | 12.148070 |
| Native action, seconds | 0.099791 |

These scoped action timings exclude cached-model capture and lift reconstruction. They are not whole-consumer timings. The old full computations continue at their immutable source/provider pair; their results remain required.

The independent mathematical assertion uses 2U plus the negative orthogonal root lattice. Signed coordinate permutations are transitive on each fixed-weight discriminant family. The null and norm-two covering vectors have the prescribed squares, are primitive, have the required divisors and represent the exact nonzero weight orbit. Odd squares have no representatives in this even lattice. A separate unimodular 2U assertion checks the trivial class with a zero complement, including its empty odd-square family. Both assertions pass in 21.85 seconds.

Optional evidence: `.tmp/o-plus-N-8e5436f-first-live.json`, `.tmp/o-plus-N-cached-covering-model.pickle`, `.tmp/o-plus-N-covering-class-action-benchmark.py`, `.log` and `.json`. The benchmark uses actual owned lattices, their based embedding and the retained source actions. It defines no substitute model object.
