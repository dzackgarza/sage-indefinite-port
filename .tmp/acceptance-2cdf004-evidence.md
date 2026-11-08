# Original assertions at 2cdf004 / bc6bd572

Port: `2cdf004966bb1d16c79f0a19ddc79a0c0f7e63ae`.
Provider: `bc6bd57262e576194e3c8764c83f3bb1d081ba52`.
Both source revisions were published. The dedicated proof worktree and provider remain unchanged.

| Unchanged node | Result | Original limit | JUnit case seconds | Required mathematical value |
| --- | --- | ---: | ---: | --- |
| `tests/acceptance/test_isotropic_orbits.sage::test_stabilizer_images_have_dutour_sikiric_hulek_indices[I_1_1]` | pass | 2400 | 2272.321 | Index of the O+(N) line stabilizer image in O(D(N)) is 1 |
| `tests/acceptance/test_orthogonal_group.sage::test_discriminant_image_of_a7_plus_diagonal_is_all_of_o_q` | terminal at limit | 1800 | unavailable | Full discriminant-image assertion remains unproved |

For I_1_1, N=U+U(2)+E8(-2) and the marked line is the first basis line of U, as fixed by the unchanged fixture. The actual image is formed from discriminant morphisms of the returned full integral stabilizer generators. The assertion compares its actual order with O(D(N)); it is not a count of supplied generators. JUnit has no failure, error or skip. Total time is 2280.657 seconds. The process and its process group are terminal.

A7's final stack reaches small-configuration transport, integral restriction and selected-generator linear extension, including fraction-zero equality. Its final assertion is unreached. The pair predates the direct generator-image law in 0580f98d; the active bf7f5c0/cab1d797 proof includes that law.

Optional logs/JUnit are `.tmp/enriques-I_1_1-2cdf004-bc6bd57.{log,xml}` and `.tmp/a7-2cdf004-bc6bd57.{log,xml}` in `sage-indefinite-proof-2cdf004`. Final-pair evidence remains required. The full unchanged I_1_1 assertion is now active at 583b36d/16c798db under the same original limit, including the direct generator-image law, retained coset transitions and finite-coefficient span lifts. Preserve that live proof until terminal.
