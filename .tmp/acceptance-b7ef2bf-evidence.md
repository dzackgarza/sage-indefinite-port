# Frozen consumer evidence

Port: `b7ef2bf74cc871d2c8a9957f9320dc33f8887f41`.
Provider: `16485a149e2ac30c660d2eb3083da2b01731fab2`.
Runtime: shared Sage container, Python 3.14.10. The test assertions and fixtures are unchanged.

| Exact node | Result | Mathematical assertion | Testcase seconds |
| --- | --- | --- | ---: |
| `tests/acceptance/test_orthogonal_group.sage::test_lorentzian_generators_match_polyhedral_common[LORENTZ_OrthogonalGroup2]` | Pass; JUnit has no failure, error or skip | The computed orthogonal generators have all recorded discriminant and prime-reduction images required by the frozen assertion. | 404.058 |

The optional local JUnit is `.tmp/lorentz2-b7ef2bf-16485a1.xml` in the fixed port proof worktree. The process exits successfully. This establishes the unchanged assertion of this specimen.

A7 remains open at this pair. Other required final-pair cases need their complete case reports.
