# Stable Dawes line coverage

At published port 616a7feefac0b9abe91cf3d01aa347245b28504b and provider cfee18926582528f6f67b3c4b6fa5d8acd329e2b, the public `L.stable_orthogonal_group().intersection(L.O_plus()).cusps(1)` computation finishes for the unchanged `dawes_2U2_A2_stable` Gram. The standalone observer uses the original boundary bound and checks exact divided discriminant classes.

| Assertion | Result |
| --- | --- |
| Line classes | 9 |
| Equality with the independently required divided-class set | Pass |

The independent null vectors are (1,0,0,0,0,0), (0,1,0,0,0,0), (0,0,1,0,0,0), (0,0,0,1,0,0), (1,0,1,0,0,0), (1,0,0,1,0,0), (0,1,1,0,0,0), (0,1,0,1,0,0), and (1,1,1,-1,0,0). Their distinct divided classes give the independent lower bound; `references/dawes-2U2-A2-boundary.md` records the complete pinned-source classification.

Planes, incidence, all subgroup witnesses and the whole frozen boundary assertion remain open. The protected fixture correction still requires explicit approval. Nested diagnostic phase times overlap and are not a total runtime.

Optional reproduction: `.tmp/dawes-line-phases-616a7fe.py`, `.log` and its phase/input directory.

## Complete plane coverage at the next public pair

The independent public building computation at port 8e5436fafefea1ca652390a29879992ce2dee700/provider cfee18926582528f6f67b3c4b6fa5d8acd329e2b also establishes complete plane coverage. It checks exact equality of all returned plane divided-class sets with the independently constructed maximal isotropic planes of the 2-primary discriminant quadric. The full line divided-class set agrees too.

| Public arithmetic-subgroup assertion | Result |
| --- | --- |
| Line classes and complete independent divided-class set | 9; pass |
| Plane classes and complete independent discriminant-plane set | 6; pass |
| Full incidence and all subgroup transporters with inverses | Open |

The same process reaches its original bound while `tits_building_incidence` matches a plane cusp through `isotropic_subspace_isometry -> IntegralStructureAction.transporter -> FiniteCommensurabilityQuotient.restricted_automorphism`. The successful plane and line assertions appear explicitly in `.tmp/dawes-independent-building-8e5436f.log`. That terminal frame does not establish cumulative cost. A new full independent building computation uses published 3f1af64/b2b3478c and still requires every incidence and every integral subgroup witness. The protected fixture correction remains pending.
