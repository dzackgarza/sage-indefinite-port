# Arithmetic cusp and incidence retention

Provider: `aa9b99f01f1107af240f88205272112f1c78c493`.
Cusp/incidence retention computation: `f1584783903c4f5de06a94385b0f1ab6dd124da8`.

The original stable Dawes graph at 2e36cff/7960e1ae is terminal under its original 3400-second bound. Its final stack is abbreviated: matrix coefficient extraction through `_matrix_column_coordinates`, `_matrix_index` and finite-index-set membership. A preceding complete live stack reaches `tits_building_incidence -> cusps(2) -> finite_ordered_set -> cusp equality -> plane stabilizer -> integral orbit-span construction`. The unchanged consumer already requested both cusp ranks before incidence. That second construction repeated classification and exact orbit comparisons.

`PredicateSubgroups.cusps` now retains its constructed quotient on the subgroup/rank. Incidence matching retains the exact subgroup transporter obtained for each successful membership decision and supplies that same witness to the incidence constructor. The original universe check, unique vertex decisions and complete flag family remain. The subgroup isotropic transporter also returns an ambient witness already in the subgroup before constructing a target parabolic; the existing adjustment operation previously returned that same witness after the eager parabolic calculation.

## Independent SO(2U) building

The determinant-one maximal parabolic splits the O orbit of maximal null planes into two SO orbits. A line parabolic contains the determinant-minus-one exchange of its perpendicular hyperbolic pair, so its O orbit remains one SO orbit. GL2(Z) in each maximal parabolic is transitive on its primitive lines, giving one flag per plane orbit.

The owner assertion `test_special_orthogonal_split_building_retains_exact_flag_incidence_and_transporters` verifies the complete counts, every incidence endpoint, subgroup membership of both term transporters, both directions of sublattice transport, inverse laws and every ambient basis pairing.

| Mathematical quantity | Verified value |
| --- | ---: |
| Line vertices | 1 |
| Plane vertices | 2 |
| Flag incidences | 2 |
| Cusp-retention source, complete diagnostic call seconds | 90.115 |
| With admitted-witness return, complete diagnostic call seconds | 77.521 |
| Separate diagnostic bound seconds | 180 |

Both bounded direct runs execute the entire unchanged owner assertion and exit successfully. These are functional mathematical results. The normal owner checks remain red at their unchanged 60-second hang bound and 20-second call gate; no claim of owner QC success follows from the diagnostic runs. An initial import also exceeds the unchanged three-second import gate. The diagnostic times are separate executions, not a whole-consumer speed certificate.

## Self-orbit subgroup witnesses

`test_maximal_null_marking_has_an_exact_self_transporter_in_the_special_group` verifies maximal coordinate null markings in U, 2U, 3U and 4U. Each witness has the exact endpoints, lies in SO, preserves the null sublattice in both directions, satisfies the inverse law and preserves the full ambient pairing.

| Owner assertion input | Call seconds | Existing gates |
| --- | ---: | --- |
| U | 1.30 | pass |
| 2U | 0.73 | pass |
| 3U | 1.36 | pass |
| 4U | 2.12 | pass |
| Selected session total | 7.54 | pass |

The generated initial parameter labels rank0 through rank3 denote the listed inputs; the committed test uses explicit U/2U/3U/4U labels. The assertions are unchanged.

The original stable Dawes graph still requires its complete frozen figure at the new published port/provider pair. Its graph vertices, integral stabilizer kernels, incidence multiplicities and exact subgroup lifts remain required. The Gamma2 fixture decision and all final-pair/full-source/arithmetic obligations remain open.
