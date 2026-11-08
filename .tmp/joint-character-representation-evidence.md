# Joint finite-character representations and cosets

Provider source: `a27ed65c`.
Provider with the signed-action assertions: `2dfd20b685cdb57361f280553345f207864bd75f`.

For chi:G->Q and H=chi^-1(K), the left G-set G/H is the orbit of K in Q/K through chi. `OrthogonalCharacterQuotient.character_representation` now supplies the owned joint homomorphism into the categorical group product of its discriminant and sign factors. The defining subgroup uses the retained discriminant targets and sign kernels. Complete finite-character cosets reuse `_preimage_coset_space`, preserving the actual image orbit, ambient lifts and Schreier section.

The product construction transports finite native GAP factors through their permutation isomorphisms, retaining the comparisons in every projection and universal mediating map. The character map is obtained by factoring its cone of component homomorphisms. Thus the product is formed in groups and its representation remains private.

The independent O(U(3)) assertion uses its four signed coordinate isometries. Reduction modulo three distinguishes them all. It checks the homomorphism on every pair, the exact kernel condition, all cosets and their live transporter lifts. The target product has order eight while the actual joint image has order four. The U assertions separately check determinant, real spinor and their joint kernel: the exact coset stabilizers, inverse action and all transporter lifts.

| Owner assertion | Mathematical result | Call seconds |
| --- | --- | ---: |
| O(U(3)), joint discriminant/determinant | Faithful image order four in target order eight, four cosets | 6.35 |
| O(U), determinant kernel | Index two with exact action and lifts | 0.59 |
| O(U), real spinor kernel | Index two with exact action and lifts | 1.25 |
| O(U), joint sign kernel | Index four with exact action and lifts | 1.29 |

The selected session passes under unchanged lint, call and hang gates; total 11.15 seconds. The existing product-of-two-swap-groups assertion also passes. The selected archived `Groups().product` placement assertion and mixed free/cyclic finite-presentation placement assertion fail identically at both the changed source and immutable provider 01cfd2db. They are existing product-presentation defects; logs retain the exact assertions. Their files also have existing TL13 debt. They do not provide passing product evidence.

The observed rank-four SO moving-maximal diagnostic originally reaches coset equality in full lattice-isometry composition at 2e36cff/7960e1ae. The new finite-character route still needs that public marking assertion. General subgroup coverage, including congruence/level conditions and nonnormal finite targets, remains required. These small examples do not discharge the full retained finite-image scope or the original large computations.
