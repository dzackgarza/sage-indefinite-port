# Frozen consumer evidence

Port: `e66eb94b1743a3265503ef4b917b64098e3e86d0`.
Provider: `40642a1720ee295dab4596e66c324623fe38569e`.
Runtime: shared Sage container, Python 3.14.10. The test assertions and fixtures are unchanged.

| Exact node | Result | Mathematical assertion |
| --- | --- | --- |
| `tests/acceptance/test_orthogonal_group.sage::test_mertens_omega_and_minus_one_generate_o_l` | Pass; JUnit has no failure, error or skip | The computed orthogonal generators have the recorded discriminant image and finite reductions of the Mertens generators together with minus identity. |

The local JUnit is `.tmp/mertens-e66eb94-40642a1.xml` in the fixed port proof worktree. Its testcase time is 54.752 seconds. The test exits successfully. This result establishes the assertion of this specimen, not full integral group equality from finite images alone.

The unchanged I_1_1 consumer is still open at this pair. No complete stabilizer-family record or final assertion is available. All other required final-pair cases remain open until their exact reports are inspected. A7 at b56d77c/40642a1 and Dawes at ab180f8/4b74095 are separate fixed proof pairs.
