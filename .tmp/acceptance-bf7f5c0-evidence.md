# Selected assertions at bf7f5c0 / cab1d797

Port: `bf7f5c0e68738729c1653b8f61b60376acd15fc0`.
Provider: `cab1d79758606479ca271f67be49ff67f2c488a4`.
Both revisions were published on main. The committed dependency and uv.lock select this provider; `uv lock --check` passes in the shared Sage container.

| Unchanged node | Result | JUnit call seconds | Mathematical assertion |
| --- | --- | ---: | --- |
| `tests/acceptance/test_double_cosets.sage::test_double_cosets_partition_the_rational_group[DoubleCoset_n55_big5760_sma1440_vf47_idx0]` | pass | 2.001 | Actual group orders and complete double-coset partition |
| `tests/test_general_isotropic_markings.py::test_public_pointwise_maximal_stabilizer_in_three_u_has_all_three_skew_parameters` | pass | 12.183 | Every full integral skew parameter occurs in the pointwise maximal isotropic stabilizer of 3U |
| `tests/test_general_isotropic_markings.py::test_complete_flag_transporters_respect_the_requested_special_orthogonal_group[negative-determinant]` | pass | 102.931 | The O(3U) transporter exists and its intersection with SO(3U) is empty |
| `tests/test_general_isotropic_markings.py::test_complete_flag_transporters_respect_the_requested_special_orthogonal_group[positive-determinant]` | pass | 8.134 | The complete flag transporter belongs to SO(3U), with all naturality and inverse equations |

The frozen node finishes in 4.94 seconds total. The selected flag consumers finish in 128.69 seconds total. JUnit has no failures, errors or skips. Optional local logs/JUnit are `.tmp/double-coset-bf7f5c0-cab1d797.{log,xml}` and `.tmp/flag-consumers-bf7f5c0-cab1d797.{log,xml}`. These selected assertions do not establish the remaining frozen corpus or the original large computations.

Provider cab1d797 also passes the S3/H full stabilizer and double-coset assertions, both supplied free-cover orders, the C4→C2 kernel and infinite cyclic parity kernel under the owner limits. A separate attempt to collect `tests/groups/test_orbit_generator_permutations.sage` at that committed provider exits before collection because the session import takes 5.65 seconds against its unchanged 3-second ceiling. That attempt supplies no mathematical result; the assertions remain required under the owner gate.
