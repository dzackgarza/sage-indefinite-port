# Finitely supported coordinate evaluation

Provider: `01cfd2db26e5bbd4c794819781d80d3417575dfe`.

`_SparseFreeModuleElement.__call__` normalizes the basis label, returns a stored coefficient directly and constructs the declared ring zero for an absent coefficient. The previous dictionary default constructed that zero even when the coefficient was present. This operation occurs in owned tensor coordinate evaluation in the unresolved original computations.

The owner assertions verify present and absent coordinates, their scalar parents and every coordinate of the zero vector over ZZ and QQ. Together with the existing public free-module category and image/kernel assertions, the selected checks pass. The existing category test has deferred TL01 lint debt; the check uses its established `--no-test-lint` route. Mathematical call and hang bounds remain unchanged.

| Selected check | Result |
| --- | --- |
| Coordinate values over ZZ and QQ | pass |
| Public free-module category and image/kernel assertions | pass |
| Selected session total seconds | 1.79 |

An isolated measurement reads a present coordinate of the same owned QQ free-module element. Each repetition performs 3000 reads, with equality verified outside the timed interval. The conditional implementation retains the same normalization and label membership check.

| Repetition | Eager-default seconds | Conditional seconds |
| --- | ---: | ---: |
| 1 | 0.880855 | 0.487035 |
| 2 | 0.878019 | 0.517655 |
| 3 | 0.994128 | 0.577865 |

These measurements concern coordinate lookup only. The full original A7, Allcock, O_plus_N and Dawes assertions remain required at their recorded immutable pairs and at the final published pair.
