# Exact finite-character transporter correction

Provider: 438c23a5bda7504ae2717d6bb4ab477f731438a1.
Committed red assertions: 662374dc for the stable-component line search and daea14d0 for the nonnormal false-negative decision.
Owner proof: `sage-extensions/results/stable-line-transporter.md` at the provider revision.

Let Gamma=chi^-1(S). The input supplies an ambient witness g:x->y and full generators of the target stabilizer P_y. A correction t exists exactly when chi(P_y) intersects S*chi(g)^-1 in the owned composition law. The native character model can reverse that law. In that case compute u in the native stabilizer image intersect S*a and use u^-1. The stabilizer's actual word section raises the selected element to an integral isometry; the exact character comparison and subgroup admission remain.

In the named negative A1 cubed, H imposes determinant one and fixing [e/2] in the discriminant. The signed exchange (-e,h,f) is an explicit H witness f->h. A supplied three-cycle and the full signed-permutation stabilizer of h previously returned None. The same target's e->h request is negative because H fixes [e/2]. Both corrected decisions pass, including full inverse and form equations.

The stable-component line assertion uses a supplied Eichler transvection in 2U(2)+A2 and in the A2.twist(-1) variant. The named A2 has Gram [[-2,1],[1,-2]]. Each public search must return an actual subgroup witness with both line inclusions, inverse law and every ambient basis pairing.

| Assertion | Mathematical result | Normal gate |
| --- | --- | --- |
| Nonnormal positive/negative correction | pass, call 3.94 seconds | pass |
| Stable-component line, signature (2,4), named A2 | every witness equation passes, call 51.91 seconds | 20-second call gate red |
| Signature (4,2), A2.twist(-1) | precursor full assertion passes mathematically, call 38.61 seconds; complete exact published function also passes in 17.22694 seconds | precursor call gate red; final explicit pytest node stops at the star-import gate before collection |
| Joint proper image, determinant/spinor/joint conditions and nonnormal left cosets/lifts | selected assertions pass, session 20.20 seconds | pass |

The named-A2/nonnormal selected session has both mathematical assertions pass in 58.89 seconds and exits red on the named-A2 call gate. The exact final A2.twist(-1) pytest node at immutable 438c23a5 stops before collection: its star import takes 5.11 seconds against the three-second gate. A separate in-process bounded diagnostic executes that exact published function unchanged, with the same parameter and every assertion; it passes in 17.22694 seconds. Optional logs are `.tmp/stable-positive-A2-line-transporter-438c23a5.log` and `.tmp/stable-positive-A2-line-transporter-438c23a5-direct.log`. The latter is functional mathematical evidence, not owner QC success.

The full original O_plus_N line-and-plane assertion, the failed original Dawes graph and the final frozen suite remain required. The stable Dawes subgroup is a normal kernel; the nonnormal false-negative specimen alone cannot explain its missing incidence endpoint. Its actual failed term and membership decisions still require the active old-pair capture.
