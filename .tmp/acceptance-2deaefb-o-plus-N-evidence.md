# Original O_plus_N assertion: verified line count and unresolved plane count

Port: `2deaefb8cb80e2cdaa85b91c84a70feeebcff923`.
Provider: `2dfd20b685cdb57361f280553345f207864bd75f`.

The unchanged node `tests/acceptance/test_isotropic_orbits.sage::test_o_plus_n_has_two_line_and_two_plane_orbits` is terminal under its original 2400-second bound. The full case fails to complete. Its process group is empty.

The final main-thread stack reaches test line 59, inside `isotropic_classes(n, 2)`. The immediately preceding line 58 is the line-orbit assertion. That sequential assertion therefore completed successfully. `tests/fixtures/unpolarized_enriques.json:183–190` supplies its expected count for the component-preserving group. The plane assertion has not completed and no case pass or JUnit success follows from the line result.

| Required result | Verified status |
| --- | --- |
| O_plus_N primitive line orbit count | 2, assertion completed |
| O_plus_N primitive plane orbit count | required 2, assertion incomplete |
| Full unchanged frozen node | terminal at original bound |

The failed plane enumeration follows `subspace_extensions -> extension_stabilizer_images -> pointed_automorphism_generators -> isotropic_stabilizer_action -> _extend_perpendicular_isometry -> partial.witt_extension`. The provider's `witt_extension` reaches `points.inclusion()(points(extension))`; admission into the extension locus compares two actual module morphisms. Their composites are evaluated on the basis through owned rational-coordinate arithmetic, with the last recorded frame in fraction equality and scalar addition.

Inspect that precise extension-locus admission and its full endpoint framings before selecting a provider repair. Preserve the comparison itself, every partial-isometry equation, exact subgroup meaning, full integral kernels and all source extension witnesses. The native full-Schreier-family realization has its separate independent evidence in `.tmp/integral-schreier-columns-evidence.md`. The full plane classification and final-pair unchanged node remain required.
