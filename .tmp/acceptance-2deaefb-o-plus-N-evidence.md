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

## Extension-admission diagnosis

Provider d0431b3f5c3ca0031dd74ea6fe3df649c269cfe9 already caches `_framed_coordinate_morphism`. The extension solver lowers the partial map and both inclusions before solving the exact pairing equations. The terminal stack alone does not show that extension admission dominates the plane computation.

An independent diagnostic uses V=U+<-2>^(n-2), H=<e,a_1,...,a_(n-2)> and the partial isometry fixing e and negating every a_i. H is coisotropic, with H-perpendicular=<e>. Each admitted extension satisfies the whole prescribed square and the full inverse law. Times below include cProfile overhead for admission only.

| Ambient rank | Locus seconds | Native solve seconds | Raise extension seconds | Full admission seconds |
| --- | ---: | ---: | ---: | ---: |
| 4 | 0.05185 | 0.04373 | 0.06805 | 0.09641 |
| 12 | 0.03708 | 0.14028 | 0.31802 | 0.12034 |

These simple maps do not reproduce the original timeout. They give no basis for replacing the module-morphism comparison or weakening any extension equation.

The current diagnostic executes the original public O_plus_N plane computation at port 8a19c4440a5e223a302e59ac74e7c1cf26509daa/provider d0431b3f5c3ca0031dd74ea6fe3df649c269cfe9, with a 2400-second in-process bound. Python monitoring records each actual Witt-extension duration and retains the solver's native source/target Gram matrices, inclusion/prescribed-image rows and solved extension. The latest and slowest completed inputs are optional local captures in `.tmp/o-plus-N-witt-capture/`; the log is `.tmp/o-plus-N-witt-capture.log`. It skips the separate line computation and therefore cannot establish the whole frozen node. Inspect the actual captured equations and duration distribution before choosing an owner repair. No source repair is justified by the small examples alone.
