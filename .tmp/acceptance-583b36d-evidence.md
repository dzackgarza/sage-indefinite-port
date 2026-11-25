# Original I_1_1 and spanning-cell assertions at 583b36d / 16c798db

Port: `583b36d625bcbadcf1ef397b4140230ce3fa300e`.
Provider: `16c798dbea898555a1373686bec1bb893edb973a`.
Both revisions were published; the committed dependency and lock select this provider.

| Unchanged node | Result | Original limit seconds | JUnit case seconds | Mathematical value |
| --- | --- | ---: | ---: | --- |
| `tests/acceptance/test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1]` | pass | 2400 | 2208.500 | Stabilizer discriminant-image index 1 |
| `tests/test_perfect_cell_transporter.py::test_integral_transporter_of_spanning_perfect_cells` | pass | unchanged | 33.137 | Integral unimodular witness carrying the complete configurations |

The full I_1_1 assertion takes 2218.490 seconds total. Its JUnit has no failures, errors or skips. It forms the actual discriminant image of the O+(N) stabilizer of Ze1 in N=U+U(2)+E8(-2), then compares its order with O(D(N)). The exact fixture and assertion are unchanged. This includes the native definite groups, retained coset transition tables, selected-generator evaluation and finite-coefficient span lifts.

A one-shot debugger snapshot reads full integral basis actions during that proof. Its measured overhead is 0.0169 seconds, included in the original limit. The mathematical state and source are unchanged. The input is only a prefix of the lazy Schreier family. The process and process group are terminal. No local I_1_1 repetition is needed before the final full frozen run.

Optional local I_1_1 log/JUnit are `.tmp/enriques-I_1_1-583b36d-16c798db.{log,xml}` in `sage-indefinite-proof-583b36d`. The spanning-cell log/JUnit are `.tmp/spanning-transporter-583b36d-16c798db.{log,xml}` in the active checkout. Their mathematical scope differs from the full A7 group assertion, which remains required. Final closure still requires the complete frozen suite at one published final port/provider pair.
