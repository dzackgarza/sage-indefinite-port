# Stable Dawes line coverage

At published port 616a7feefac0b9abe91cf3d01aa347245b28504b and provider cfee18926582528f6f67b3c4b6fa5d8acd329e2b, the public `L.stable_orthogonal_group().intersection(L.O_plus()).cusps(1)` computation finishes for the unchanged `dawes_2U2_A2_stable` Gram. The standalone observer uses the original boundary bound and checks exact divided discriminant classes.

| Assertion | Result |
| --- | --- |
| Line classes | 9 |
| Equality with the independently required divided-class set | Pass |

The independent null vectors are (1,0,0,0,0,0), (0,1,0,0,0,0), (0,0,1,0,0,0), (0,0,0,1,0,0), (1,0,1,0,0,0), (1,0,0,1,0,0), (0,1,1,0,0,0), (0,1,0,1,0,0), and (1,1,1,-1,0,0). Their distinct divided classes give the independent lower bound; `references/dawes-2U2-A2-boundary.md` records the complete pinned-source classification.

Planes, incidence, all subgroup witnesses and the whole frozen boundary assertion remain open. The protected fixture correction still requires explicit approval. Nested diagnostic phase times overlap and are not a total runtime.

Optional reproduction: `.tmp/dawes-line-phases-616a7fe.py`, `.log` and its phase/input directory.
