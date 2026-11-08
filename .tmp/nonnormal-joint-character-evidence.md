# Nonnormal joint finite-character groups

Provider: `d0431b3f5c3ca0031dd74ea6fe3df649c269cfe9`.
Residual finite-form backend computation: `6ae02422`.

The standard root lattice A1 cubed has q(x)=2(x_0^2+x_1^2+x_2^2), so its minimal vectors force every integral isometry to be a signed permutation. In its discriminant form the classes e_i/2 are exactly the elements of square 1/2. Thus O(D) is S3, and independent coordinate signs give both determinants above every discriminant permutation. The joint image is the full product C2 times S3.

Let K fix e_0/2 in S3, and H impose both rho(g) in K and det(g)=1. The transposition subgroup K is nonnormal; the assertion verifies an explicit conjugate outside K. The cosets of H have coordinates (det(g),rho(g)e_0/2). The owner assertion checks their entire coordinate set, every live coset lift, every chosen generator's left action and inverse on all cosets, and the noncommuting joint-character laws.

The full group H restricts to signed permutations of the other two axes. The determinant condition determines the remaining axis sign, so this restriction is faithful and every signed permutation of the other axes lifts. It identifies H with D8. The owner generator assertion checks that each returned generator lies in the specified H, its action on those signed roots has full order eight and its discriminant image has order two. Faithfulness and the independent description of H prove completeness of the returned integral group, rather than only a quotient image. The kernel order follows by the first isomorphism theorem.

| Mathematical quantity | Verified value |
| --- | ---: |
| Discriminant group order | 8 |
| Discriminant orthogonal group order | 6 |
| Joint image order | 12 |
| Nonnormal target K order | 2 |
| Left cosets of H | 6 |
| Full generated H order | 8 |
| Discriminant image of H | 2 |
| Integral discriminant kernel in H | 4 |

## Maintained finite-form computation

Residual primary parts now use [Sage's full torsion-form orthogonal group](https://doc.sagemath.org/html/en/reference/modules/sage/modules/torsion_quadratic_module.html#sage.modules.torsion_quadratic_module.TorsionQuadraticModule.orthogonal_group). Omitting supplied generators requests the complete group. The operation receives the actual primary module, preserving both value moduli and its Smith framing. Native matrices record row images of those Smith generators and enter the existing primary-extension construction. The classical odd and integer-valued elementary-two structural routes retain their existing computations. Sage's general residual enumeration can be expensive; large residual cases still require actual evidence.

The nonnormal case completes under the unchanged owner gates. Existing exponent-four quadratic and elementary bilinear specimens also pass, including preservation of the quadratic values on every point of the exponent-four module.

| Selected owner check | Mathematical outcome | Call seconds | Gate |
| --- | --- | ---: | --- |
| Nonnormal joint cosets and lifts at published d0431b3f | complete | 9.55 | pass |
| Exponent-four hyperbolic discriminant | group order four, all quadratic values preserved | 5.98 | pass |
| Bilinear hyperbolic discriminant | full linear group order six | 2.11 | pass |
| Full nonnormal subgroup generator action | full group and kernel retained | 24.62 | call-time failure |
| Actual named point-set construction data | retained | 0.14 | pass |

The generator session reports both mathematical assertions passed, but exits with the unchanged 20-second call-time gate red. It is not owner QC success. The finite permutation computation supplies its ordered root presentation explicitly; an abstract unordered finite set does not supply a selected ranking. Encoding reads the symmetric-group owner's selected ranking while retaining its actual point set.

The original A7 proof at 2deaefb/2dfd20b6 is terminal in configuration moment-colour norms and still requires a measured repair and its full image assertion. O_plus_N and Dawes retain their live immutable pairs. General finite-image subgroup coverage across every marking and the full final-pair frozen suite remain open.
